#
#    Copyright (C) 2021  Michelle Blom
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.


import argparse
import os
import numpy as np
import sys

from functools import partial

from multiprocessing import Pool
from multiprocessing.pool import ThreadPool

from utils import read_outcome, sample_size, ssm_sample_size,\
    tally_vs_tally_sample_size, read_ballots_txt, index_of, next_cand, \
    read_ballots_json, read_ballots_txt, read_ballots_blt, \
    read_ballots_stv, simulate_stv



# Determine if there is a candidate 'd' for which 'd' AG loser, based 
# on the AG relationships in neb_matrix, preferenced before 'loser' in
# the preference ranking 'prefs'. Return a boolean indicating whether
# there is such a candidate 'd', and the ASN of the cheapest AG 
# relationship (if there are multiple such 'd's).
#
#     ag_matrix[d][loser] will give ASN of 'd' AG loser or None if 
#         such an assertion does not exist.
#
# The function returns pair:
#
#     ag_present, ag_min_ss
#
# where neb_present is a boolean indicating whether we found a 
# 'd' such that 'd' AG loser, and 'd' is preferenced before loser
# in prefs.
#
# If loser is not present in the ranking 'prefs' the function will
# return:
#
#     False, np.inf
def rule_out_for_max(prefs, loser, ag_matrix, winners, candidates):
    ag_min_ss = np.inf
    ag_present = False

    ags_used = set()

    for p in prefs:
        if p == loser:
            return ag_present, ag_min_ss, ags_used

        if p in winners:
            # 'loser' could still get this vote as it could skip over 'p'
            continue

        ag = ag_matrix[p][loser]
        if ag != None and ag != np.inf:
            ag_present = True
            ags_used.add((p,"AG",loser, ag))
            ag_min_ss = min(ag_min_ss, ag)

    return False, np.inf, set()


        

# Merge a sequence of AG relationships that could be used to increase the
# assorter margin of an NL.
#
# Input list 'helpful_ags' will be a list of (asn, extra),
# where 'asn' is the cost of auditing the AG assertion, and 'extra' is 
# the increase to the assorter total if we incorporate those AGs. This 
# function takes a list of these AGs, and merges consecutive entries if: the
# ASN is the same. 
def merge_helpful_ags(helpful_ags, exp_merged_total):
    helpful_ags.sort()
    merged_helpful_ags = []

    merged_total = 0 
    if helpful_ags != []:
        cntr = 1
        curr_ag_asn = helpful_ags[0][0]                      
        curr_extra = helpful_ags[0][1]
        curr_desc = helpful_ags[0][2]

        lhelpfuls = len(helpful_ags)

        to_add = False
        while cntr < lhelpfuls:
            ag_asn, extra, desc = helpful_ags[cntr]
            if ag_asn != curr_ag_asn:
                merged_helpful_ags.append((curr_ag_asn, curr_extra, curr_desc))
                merged_total += curr_extra

                curr_ag_asn = ag_asn
                curr_extra = extra
                curr_desc = desc
            else:
                curr_extra += extra
                curr_desc.update(desc)

            cntr += 1                    

        merged_helpful_ags.append((curr_ag_asn, curr_extra, curr_desc))
        merged_total += curr_extra

    assert(merged_total >= exp_merged_total-0.00001 \
        and merged_total <= exp_merged_total+0.00001)

    return merged_helpful_ags


def compute_ag_matrix(candidates, ballots, ag_matrix, \
    INVALID, args, log=None):

    for cand_w in candidates:
        w = cand_w.num

        for cand_c in candidates:
            c = cand_c.num

            if c == w:
                continue

            min_w = 0
            max_c = 0

            # assertion: fpc(w) > maxc
            assorter = INVALID*0.5 # h(b) = ((b_w - b_c) + 1)/2
            for b in ballots:
                if b.prefs != []:
                    if b.prefs[0] == w:
                        # assorter value is 1 per vote
                        assorter += b.votes
                        min_w += b.votes
                        continue

                    if b.prefs[0] == c: 
                        # assorter value is 0 per vote
                        max_c += b.votes
                        continue

                # Default contribution of each instance of ballot type
                # to assorter.
                contrib = 0.5 * b.votes 

                for p in b.prefs:
                    if p == w:
                        break

                    if p == c:
                        contrib = 0
                        max_c += b.votes
                        break

                assorter += contrib

            # Compute assorter margin
            amean = assorter/args.voters

            if amean > 0.5:
                ss = sample_size(amean, args)

                ag_matrix[cand_w.num][c] = ss

                if log != None:
                    print("AG({},{}) = {}".format(cand_w.id, cand_c.id, \
                        ss), file=log)



def compute_ag_stars(winners, losers, candidates, ballots, ag_matrix, \
    first_winners, min_tvs, aud_tvs, INVALID, args, desc):

    for w in winners:
        cand_w = candidates[w]

        for c in losers:
            if c == w:
                continue

            cand_c = candidates[c]

            min_w = 0
            max_c = 0

            # assertion: fpc(w) > maxc
            assorter = INVALID*0.5 # h(b) = ((b_w - b_c) + 1)/2
            for b in ballots:
                if b.prefs != []:
                    if b.prefs[0] == w:
                        # assorter value is 1 per vote
                        assorter += b.votes
                        min_w += b.votes
                        continue

                    if b.prefs[0] == c: 
                        # assorter value is 0 per vote
                        max_c += b.votes
                        continue

                # Default contribution of each instance of ballot type
                # to assorter.
                contrib = 0.5 * b.votes 

                counted = False
                if w in b.prefs:
                    # Check if all preceding candidates are first winners
                    prior = b.prefs[:b.prefs.index(w)]

                    if prior != [] and all([p in first_winners for p in prior]):
                        # Candidate c gets votes at smallest of transfer
                        # values across winneres in prior
                        min_tv = min([min_tvs[p] for p in prior])
                        min_w += min_tv * b.votes
                        contrib = b.votes * ((1 + min_tv)/2.0) 
                        counted = True

                if not counted:
                    weight = 1

                    if c in b.prefs:
                        prior = b.prefs[:b.prefs.index(c)]
                        atvs = [aud_tvs[p] for p in prior if p in first_winners]
                        if atvs != []:
                            weight = max(atvs)

                    for p in b.prefs:
                        if p == w:
                            break

                        if p == c:
                            contrib = b.votes * ((1 - weight)/2.0)
                            max_c += weight*b.votes
                            break

                assorter += contrib

            # Compute assorter margin
            amean = assorter/args.voters

            if amean > 0.5:
                ss = sample_size(amean, args)

                ag_matrix[w][c] = ss

                desc += "AG*({},{}) = {}\n".format(cand_w.id, \
                    cand_c.id, ss)


def get_max_ballot_weight(prefs, first_winners, aud_tv):
    weights = []
    fw_prefix = []

    for p in prefs:
        if p in first_winners:
            weights.append(aud_tv[p])
            fw_prefix.append(p)    
        else:
            break

    if weights == []:
        return 1, []
    else:
        return max(weights), fw_prefix


def inner_loop(winners_on_fp, args, candidates, cands, valid_ballots,\
    INVALID,max_sample_size,ballots,min_tvs,ows,fws,losers,winners, \
    straight_iqx_verified, aud_tvs):

    np.seterr(all='ignore')

    ag_matrix = [[None for c in candidates] for o in candidates]

    desc = "---------------------------------------------\n"
    desc += "START INNER LOOP\n"

    inner_loop_assertions = set()
    max_this_loop = 0
    max_ss_mt = 0

    winners_verified = []

    ut_asns = {}

    # Check that TVs of the first winners i are at most aud_TV[i]
    for f in winners_on_fp:
        fw_i = candidates[f]
        aud_tv_i = aud_tvs[f]

        a = 1 / (1 - aud_tv_i)
        thresholdT = a * valid_ballots / (args.seats + 1)
        thresholdProp = thresholdT / valid_ballots
        threshold = 1 - thresholdProp
        
        tally_others = valid_ballots - fw_i.fp_votes

        ss_i = ssm_sample_size(threshold,tally_others,INVALID,args)
        ut_asns[f] = (aud_tv_i, ss_i)

        if ss_i != np.inf:
            winners_verified.append(f)

        desc += "AUD TV for {}, {}, ss {}\n".format(fw_i.id, \
            aud_tv_i, ss_i)

        inner_loop_assertions.add((fw_i.num, "MT", None, \
            (aud_tv_i, ss_i)))

        max_ss_mt = max(max_ss_mt, ss_i)

    partial_ss = max_this_loop

    # For remaining winners, is it the case that they cannot be 
    # eliminated prior to all reported losers? We will
    # look at this from TWO perspectives. (1) compute AG*'s
    # between the remaining winner r and opponents using our upper
    # and lower bounds on the initial winner's transfer value. If 
    # we can form enough AG*'s, we can verify the winner.
    # (2) Forming NL assertions, potentially making use of our
    #  set of AG* relationships. We will use the 
    # cheaper of the two options.
    compute_ag_stars([ow_i.num for ow_i in ows], losers, \
        candidates, ballots, ag_matrix, winners_on_fp, min_tvs, \
        aud_tvs, INVALID, args, desc)

    compute_ag_stars(losers, losers,  candidates, ballots, \
        ag_matrix, winners_on_fp, min_tvs, aud_tvs, INVALID, \
        args, desc)

    qthresh = 1.0/(args.seats + 1);

    for ow_i in ows:
        ags = {c : ag_matrix[ow_i.num][c] for c in losers \
            if ag_matrix[ow_i.num][c] != None}
               
        # Can we show that ow_i's first preferences plus all the 
        # ballots that would flow from candidates l for which
        # AG*(ow_i, l) is greater than a quota? 
        totvotes = 0
        ss_ag_iqx = 0
        ss_iqx = np.inf

        iqx_assertions = set()
        helpful_ags = []
        pot_margin_inc = 0
        for b in ballots:
            if b.prefs[0] == ow_i.num:
                totvotes += b.votes
                continue
            
            if not ow_i.num in b.prefs:
                continue

            value = 1

            if b.prefs[0] in winners_on_fp:
                value = min_tvs[b.prefs[0]]

            prefs = [p for p in b.prefs if not p in winners_on_fp]

            if prefs[0] == ow_i.num:
                totvotes += value*b.votes
                continue

            descs = set()
            max_ags_here = 0
                                           
            o_idx = prefs.index(ow_i.num)                 
            for d,dval in ags.items():
                if d in prefs:
                    idx_d = prefs.index(d)
                    if idx_d < o_idx:
                        prefs.remove(d)
                        o_idx -= 1
                        rag = (ow_i.num,"AG*",d,dval)
                        descs.add(rag)
                        max_ags_here=max(max_ags_here,dval)

            if prefs[0] == ow_i.num:
                helpful_ags.append((max_ags_here, value*b.votes, descs))
                pot_margin_inc += value*b.votes

        ss_ag_iqx = 0
        merged_helpful_ags = merge_helpful_ags(helpful_ags, \
            pot_margin_inc)

        while totvotes < args.quota and \
            merged_helpful_ags != []:
            
            ag_asn, extra_contrib, descs = merged_helpful_ags.pop(0)

            totvotes += extra_contrib
            ss_ag_iqx = max(ss_ag_iqx, ag_asn)

            iqx_assertions.update(descs)


        if totvotes > args.quota:
            ss = ssm_sample_size(qthresh, totvotes, INVALID, args)

            for ag_asn, extra_contrib, descs in merged_helpful_ags:
                if ss < ag_asn:
                    break

                totvotes += extra_contrib
                ss = ssm_sample_size(qthresh, totvotes, INVALID, args)

                ss_ag_iqx = max(ss_ag_iqx, ag_asn)
                iqx_assertions.update(descs)


            iqx_assertions.add((ow_i.num, "IQX", None, ss))
            desc += "Can form IQX({}) with sample size {}, AG*s\n".format(\
                ow_i.id, ss, ss_ag_iqx) 
            ss_iqx = max(ss, ss_ag_iqx)
            

        nl_assertion_set = set()
        max_with_nls_i = 0

        # Determine NL's between original losers and winner ow_i
        for c in cands:
            if c in winners:
                continue

            cand_c = candidates[c]

            curr_nl = (ow_i.num, "NL", c)

            min_ow_i = 0 # Min tally for ow_i.
            max_c = 0 # Max tally for candidate c.

            # Keep running tally of total votes we can increase the
            # margin of assertion 'ow_i' NL 'c' with 'AG*' 
            # relationships
            pot_margin_inc = 0

            helpful_ags = []

            # Assertion: min_ow_i > maxc
            assorter = INVALID*0.5
            for b in ballots:
                if b.prefs != []:
                    if b.prefs[0] == ow_i.num:
                        assorter += b.votes
                        min_ow_i += b.votes
                        continue

                if b.prefs[0] == c:
                    max_c += b.votes
                    continue

                # is this a ballot for 'c' over 'ow_i'?
                c_in = c in b.prefs
                ow_in = ow_i.num in b.prefs

                weight, fwprefix = get_max_ballot_weight(b.prefs,\
                    winners_on_fp, aud_tvs)

                c_idx=b.prefs.index(c) if c_in else np.inf
                o_idx=b.prefs.index(ow_i.num) if ow_in else np.inf

                if c_in and c_idx < o_idx:
                    # These ballots could not go to 'ow' but may go
                    # to 'c'.

                    # Could we exclude ballots by using an AG*?
                    is_ag, ag_asn, descs = rule_out_for_max(\
                        b.prefs, c, ag_matrix, winners, candidates)
                            
                    contrib = b.votes*((1-weight)/2.0)

                    if is_ag:
                        alt_contrib = b.votes*0.5
                        helpful_ags.append((ag_asn, alt_contrib -\
                            contrib, descs))
                        pot_margin_inc += alt_contrib-contrib

                    assorter += contrib
                    max_c += b.votes

                elif ow_in:
                    # Check if all candidates before ow_i are
                    # first winners.
                    if b.prefs[len(fwprefix)] == ow_i.num:
                        # Give votes to ow_i at minimum of the 
                        # transfer values of first winners.
                        minval = min([min_tvs[f] for f in fwprefix])
                        min_ow_i += minval * b.votes
                        assorter += b.votes * ((1 + minval)/2.0)
                    else:
                        # If we remove all cand 'd' for which ow AG*
                        # d that appear before ow in the ballot,
                        # will 'ow' be the first ranked cand? If so,
                        # we could add these ballots to the min
                        # tally of 'ow'.
                        prefs = b.prefs[:]

                        descs = set()
                        max_ags_here = 0
                                           
                        for d,dval in ags.items():
                            if d in prefs:
                                idx_d = prefs.index(d)
                                if idx_d < o_idx:
                                    prefs.remove(d)
                                    o_idx -= 1
                                    rag = (ow_i.num,"AG*",d,dval)
                                    descs.add(rag)
                                    max_ags_here=max(max_ags_here,\
                                        dval)

                        assorter += 0.5*b.votes

                        if prefs != [] and prefs[0] == ow_i.num:
                            # These ballots would have originally
                            # had a contribution of 0.5 to the
                            # assorter. By utilising these AG*'s we
                            # can increase the contribution of these
                            # ballots to 1 each.
                            helpful_ags.append((max_ags_here,  \
                                0.5*b.votes, descs))

                            pot_margin_inc += b.votes*0.5

                        else:
                            fwprefix = []
                            minws = []
                            for p in prefs:
                                if p in winners_on_fp:
                                    fwprefix.append(p)
                                    minws.append(min_tvs[p])
                                else:
                                    break

                            if fwprefix != [] and \
                                prefs[len(fwprefix)] == ow_i.num:
                                minval = min(minws) if b.prefs[0]\
                                    in winners_on_fp else 1

                                base_contrib = 0.5*b.votes
                                alt_contrib = b.votes * ((1 + \
                                    minval)/2.0)
                                dconfig = alt_contrib-base_contrib

                                if dconfig > 0:
                                    helpful_ags.append(\
                                        (max_ags_here, \
                                        dconfig, descs))

                                    pot_margin_inc += dconfig
                else:
                    assorter += 0.5*b.votes

            # Max ASN of any AG*'s used to increase assorter 
            # margins when forming NLs.
            max_ags_used = 0  
            merged_helpful_ags = merge_helpful_ags(helpful_ags, \
                pot_margin_inc)

            # Incorporate use of  AG*'s that either make the
            # assertion possible, or whose ASN is already
            # within/equal to current lower bound on audit
            # difficulty.
            while assorter/args.voters <= 0.5 and \
                merged_helpful_ags != []:
            
                ag_asn, extra_contrib, descs = \
                    merged_helpful_ags.pop(0)

                assorter += extra_contrib
                max_ags_used = max(max_ags_used, ag_asn)

                nl_assertion_set.update(descs)


            # Can we reduce the sample size required for the
            # assertion by including more AG*'s?
            amean = assorter/args.voters
            if amean > 0.5:
                # Current sample size for assertion
                ss = sample_size(amean, args) 

                # Go through remaining AG*'s 
                for ag_asn, extra_contrib, descs in \
                    merged_helpful_ags:
                    if ss < ag_asn:
                        break

                    # would reducing margin by votes reduce sample
                    # size by more than increase caused by including
                    # AG*?
                    assorter += extra_contrib
                    amean = assorter/args.voters

                    ss  = sample_size(amean, args)
                    max_ags_used = max(max_ags_used, ag_asn)
            
                    nl_assertion_set.update(descs)

                max_with_nls_i=max(max_with_nls_i,max(max_ags_used,ss))
                desc += "NL({},{}) = {}, AG*'s used {}\n".format(\
                    ow_i.id, cand_c.id, ss, max_ags_used)

                nl_assertion_set.add((ow_i.num, "NL", c, \
                    max(ss,max_ags_used)))
            else:
                max_with_nls_i = np.inf
                desc += "NL({},{}) NOT POSSIBLE\n".format(\
                    ow_i.id, cand_c.id)

        straight_iqx_asn = straight_iqx_verified[ow_i.num][0] if \
            ow_i.num in straight_iqx_verified else np.inf

        if ss_iqx < args.voters or max_with_nls_i < args.voters or \
            straight_iqx_asn < args.voters:
            winners_verified.append(ow_i.num)
            partial_ss = max(partial_ss, min([ss_iqx, max_with_nls_i,\
                straight_iqx_asn])) 

        if straight_iqx_asn < ss_iqx and straight_iqx_asn < max_with_nls_i:
            desc += "CHOOSE STRAIGHT IQX\n"
            max_this_loop = max(max_this_loop, straight_iqx_asn)
            inner_loop_assertions.update(straight_iqx_verified[ow_i.num][1])
        
        elif ss_iqx < max_with_nls_i:
            desc += "CHOOSE IQX over NLs\n"
            max_this_loop = max(max_this_loop, ss_iqx)
            inner_loop_assertions.update(iqx_assertions)
        else:
            max_this_loop = max(max_this_loop, max_with_nls_i)
            inner_loop_assertions.update(nl_assertion_set)
            
       
    max_in_loop = max(max_this_loop, max(max_ss_mt, max_sample_size)) 
    partial_ss = max(partial_ss, max(max_ss_mt, max_sample_size))

    desc += "Max in loop {}, {}, partial {}\n".format(max_in_loop, {f.id : \
        aud_tvs[f.num] for f in fws}, partial_ss)

    return max_in_loop, max_ss_mt, max_this_loop, desc, aud_tvs, \
        inner_loop_assertions, winners_verified, partial_ss, ut_asns


def create_upper_neighbours(curr_aud_tvs, winners_on_fp, deltat, maxtv):
    nboors = []

    for fw in winners_on_fp:
        new_aud_tvs = curr_aud_tvs.copy()
        if new_aud_tvs[fw] + deltat < maxtv:
            new_aud_tvs[fw] += deltat
            nboors.append(new_aud_tvs)

    return nboors


def outer_loop(args, ows, fws, losers, winners, winners_on_fp,  \
    candidates, cands, valid_ballots, INVALID, ballots, max_ss, deltat, \
    deltam, maxtv, act_tvs, starting_uts, straight_iqx_verified, min_tvs):
   
    np.seterr(all='ignore')
 
    outerdesc = "------------------------------------------------\n"
    outerdesc += "START OUTER LOOP, min_tvs {}\n".format({f.id : \
        min_tvs[f.num] for f in fws})
    mintv_ss = 0
   
    lts = set()

    for fw_i in fws:
        min_tv_i = min_tvs[fw_i.num]

        if min_tv_i > 0:
            mintally = args.quota / (1 - min_tv_i)
            thresh = mintally/valid_ballots

            mintv_ss_i = ssm_sample_size(thresh, fw_i.fp_votes, \
                INVALID, args)
            outerdesc += "Sample size to show min tv of {} is {}\n".format(\
                min_tv_i, mintv_ss_i)

            lts.add((fw_i.num, "LT", None, (min_tv_i, mintv_ss_i)))
            mintv_ss = max(mintv_ss, mintv_ss_i)

    max_in_loop = np.inf

    best_inner_assertions = None

    # TODO: Keep track of sample sizes of TV points for winners_on_fp. 
    # Use to start each inner loop at a setting that is near the best
    # found partial audit ASN. 

    if starting_uts == None:
        curr_aud_tvs = {fw : act_tvs[fw] + deltat for fw in winners_on_fp}
    else:
        curr_aud_tvs = starting_uts.copy()

    tv_ub_nboors = [curr_aud_tvs] + create_upper_neighbours(curr_aud_tvs, \
        winners_on_fp, deltat, maxtv)

    curr_winners_verified = []
    best_partial_ss = np.inf

    improved = True
    partial_improved = True

    all_ut_asns = {f : set() for f in winners_on_fp}

    while improved and tv_ub_nboors != []:

        # Run inner loop for each setting of tv upper bounds
        # we are interested in. Establish what leads to the 
        # best audit.
        with ThreadPool(processes=args.cpus) as pool:
            func = partial(inner_loop, winners_on_fp, args, \
                candidates, cands, valid_ballots, INVALID, max_sample_size,\
                ballots,min_tvs,ows,fws,losers,winners,straight_iqx_verified)

            results = pool.map(func, tv_ub_nboors)

        improved = False
        partial_improved = False
        min_this_loop = np.inf

        for max_this_loop, uts_ss, asn_exc_lts, desc, nb_aud_tvs, \
            inner_loop_assertions, winners_verified, partial_ss, \
            ut_asns in results:

            for f,utfs in ut_asns.items():
                all_ut_asns[f].add(utfs)

            min_this_loop = min(max_this_loop, min_this_loop)

            outerdesc += desc
            if max_this_loop < max_in_loop or (max_this_loop == \
                max_in_loop and partial_ss <= best_partial_ss):
                best_inner_assertions = inner_loop_assertions
                max_in_loop = max_this_loop
                curr_aud_tvs = nb_aud_tvs
                curr_winners_verified = winners_verified
            
                if partial_ss < best_partial_ss:
                    partial_improved = True
                
                best_partial_ss = partial_ss
                improved = True

        if mintv_ss > best_partial_ss:
            break

        tv_ub_nboors = create_upper_neighbours(curr_aud_tvs,\
            winners_on_fp, deltat, maxtv)

    max_in_loop = max(max_in_loop, mintv_ss)

    outerdesc += "Output of outer loop is: {}, partial exc lts {}\n".format(\
        max_in_loop, best_partial_ss)

    return max_in_loop, min_tvs, outerdesc, lts, mintv_ss,\
        best_inner_assertions, curr_winners_verified,  max(best_partial_ss, \
        mintv_ss), best_partial_ss, all_ut_asns


def create_lower_neighbours(curr_min_tvs, winners_on_fp, \
    deltam, act_tvs):

    nboors = []

    for fw in winners_on_fp:
        new_min_tvs = curr_min_tvs.copy()
        if new_min_tvs[fw] > 0:
            new_min_tvs[fw] = max(0, new_min_tvs[fw] - deltam)
            nboors.append(new_min_tvs)
            
    return nboors
        


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    # Input: stv data file
    parser.add_argument('-d', dest='data')

    # Input: anticipated error rate (default value is 0)
    parser.add_argument('-e1', dest='erate1', default=0.002, type=float)
    parser.add_argument('-e2', dest='erate2', default=0, type=float)

    # Input: risk limit (default is 5%)
    parser.add_argument('-r', dest='rlimit', default=0.05, type=float)
    
    # Input: number of repetitions to perform in simulation to determine
    # an initial sample size estimation -- the quantile of the sample
    # size is computed (default is 1 repetition -- no error rate)
    parser.add_argument('-reps', dest='reps', default=20)

    # Input: seed (default is 9368663)
    parser.add_argument('-s', dest='seed', default=9368663, type=int)
    
    # Input: number of cpus to use if parallelising tasks.
    parser.add_argument('-cpus', dest='cpus', default=9, type=int)

    # Input: number of seats in election (default is 2)
    parser.add_argument('-seats', dest='seats', default=2, type=int)

    # Input: election outcome file
    parser.add_argument('-outcome', dest='outcome')
    
    # Input: Quota
    parser.add_argument('-quota', dest='quota', type=int)

    # Input: Total voters (used to express total valid+informal ballots)
    parser.add_argument('-voters', dest='voters', type=int)

    # Output: Log file 
    parser.add_argument('-log', dest='log', type=str)
    
    # Input: increment to use when exploring candidate lower and upper 
    # bounds on first winner's transfer value (for 1-quota method).
    parser.add_argument('-deltat', dest='deltat', default=0.005, type=float)
    parser.add_argument('-deltam', dest='deltam', default=0.005, type=float)
    
    args = parser.parse_args()

    log = open(args.log, "w")

    maxtv = args.seats/(args.seats + 1)

    # Read STV data file
    candidates, ballots, id2group, cid2num, valid_ballots = \
        None, None, None, None, None

    # Check for given input data type
    if args.data.endswith(".json"):
        candidates, ballots, id2group, cid2num, valid_ballots = \
            read_ballots_json(args.data)

    elif args.data.endswith(".blt"):
        candidates, ballots, id2group, cid2num, valid_ballots = \
            read_ballots_blt(args.data)

    elif args.data.endswith(".txt") or args.data.endswith(".us"):
        candidates, ballots, id2group, cid2num, valid_ballots = \
            read_ballots_txt(args.data)

    elif args.data.endswith(".stv"):
        candidates, ballots, id2group, cid2num, valid_ballots = \
            read_ballots_stv(args.data)

    else:
        print("Unsupported data file type.")
        exit()


    np.seterr(all='ignore')

    # Read STV outcome file
    outcome = read_outcome(args.outcome, cid2num)

    # Count of informal votes.
    INVALID = args.voters - valid_ballots

    ncand = len(outcome.cand)
    cands = []
    winners = [] # Reported winners
    winners_on_fp = [] # Identifiers for candidates that win seat in 1st round
    losers = [] # Reported losers

    for i in range(ncand):
        c = outcome.cand[i]
        cands.append(c)
        cd = candidates[c]

        if outcome.action[i]:
            winners.append(c)
            if cd.fp_votes > args.quota:
                winners_on_fp.append(c)
        else:
            losers.append(c)

    max_sample_size = 0

    # List of text descriptions, and associated sample sizes, of all
    # assertions needed in an audit.
    assertions_used = set()

    if outcome.action[0] != 1:
        print("Election does not satisfy first winner criterion.")
        exit()

    else: 
        # For storing computing AG*'s
        ag_matrix = [[None for c in candidates] for o in candidates]


        # First option -- can we just show that IQX(w) for each reported
        # winner?
        compute_ag_matrix(candidates, ballots, ag_matrix, INVALID, args, \
            log = log)

        straight_iqx_asn = 0
        straight_iqx_verified = {}

        qthresh = 1.0/(args.seats + 1);

        iqx_assertions = set()
        for w in winners:
            ags = {c : ag_matrix[w][c] for c in losers \
                if ag_matrix[w][c] != None}
               
            totvotes = 0
            ss_iqx = np.inf

            assertions = set()

            helpful_ags = []
            pot_margin_inc = 0
            for b in ballots:
                if b.prefs[0] == w:
                    totvotes += b.votes
                    continue

                prefs = b.prefs[:]
                if not w in prefs:
                    continue

                descs = set()
                max_ags_here = 0
                          
                o_idx = prefs.index(w)                 
                for d,dval in ags.items():
                    if d in prefs:
                        idx_d = prefs.index(d)
                        if idx_d < o_idx:
                            prefs.remove(d)
                            o_idx -= 1
                            rag = (w,"AG",d,dval)
                            descs.add(rag)
                            max_ags_here=max(max_ags_here,dval)

                if prefs[0] == w:
                    helpful_ags.append((max_ags_here, b.votes, descs))
                    pot_margin_inc += b.votes

            ss_ag_iqx = 0
            merged_helpful_ags = merge_helpful_ags(helpful_ags, \
                pot_margin_inc)

            # Incorporate use of  AG*'s that either make the
            # assertion possible, or whose ASN is already
            # within/equal to current lower bound on audit
            # difficulty.
            while totvotes < args.quota and \
                merged_helpful_ags != []:
            
                ag_asn, extra_contrib, descs = \
                    merged_helpful_ags.pop(0)

                totvotes += extra_contrib
                ss_ag_iqx = max(ss_ag_iqx, ag_asn)

                assertions.update(descs)

            
            if totvotes > args.quota:
                ss = ssm_sample_size(qthresh, totvotes, INVALID, args)

                for ag_asn, extra_contrib, descs in \
                    merged_helpful_ags:
                    if ss < ag_asn:
                        break

                    totvotes += extra_contrib
                    ss = ssm_sample_size(qthresh, totvotes, INVALID, args)

                    ss_ag_iqx = max(ss_ag_iqx, ag_asn)
                    assertions.update(descs)

                assertions.add((w, "IQX", None, ss))
                ss_iqx = max(ss, ss_ag_iqx)
                print("Can form IQX({}) with sample size {} AGs {}\n".format(\
                    candidates[w].id, ss, ss_ag_iqx), file=log) 

                if ss_iqx < args.voters:
                    straight_iqx_verified[w] = (ss_iqx, assertions)

                iqx_assertions.update(assertions)


            straight_iqx_asn = max(straight_iqx_asn, ss_iqx);


        # Check that all winners that won on first preferences have a FP
        # greater than a quota.
        fws = [candidates[w] for w in winners_on_fp]
        ows = [candidates[w] for w in winners if not w in winners_on_fp]
       
        print("Winners on first preferences: {}".format([fw.id for fw in fws]),\
            file=log)
        print("Other winners: {}".format([ow.id for ow in ows]), file=log)
 
        thresh = 1.0/(args.seats + 1);
    
        qts = set()
        qts_verified = {}

        act_tvs = {}
        min_tvs = {}

        remove_from_fws = []
        for fw in fws:
            ss = ssm_sample_size(thresh, fw.fp_votes, INVALID, args)
            print("{} ballot checks required to assert that winner {} "\
                " has a quota's worth of votes".format(ss, fw.id), file=log)

            if ss >= args.voters:
               remove_from_fws.append(fw)
            else: 
                qts.add((fw.num, "QT", None, ss))
                max_sample_size = max(max_sample_size, ss)

                qts_verified[fw.num] = (ss, (fw.num, "QT", None, ss))

                act_tvs[fw.num] = (fw.fp_votes - args.quota)/fw.fp_votes
                min_tvs[fw.num] = 0

        for fw in remove_from_fws:
            fws.remove(fw)
            winners_on_fp.remove(fw.num)
            ows.append(fw)

        qts_sample_size = max_sample_size

        assertions_used.update(qts)

        if max_sample_size >= args.voters:
            max_sample_size = np.inf 

        if max_sample_size == np.inf or fws == []:
            if straight_iqx_asn != np.inf:
                print("--------------------------------------------",file=log)
                print("Assertion set:", file=log)
                for asstn in iqx_assertions:
                    print(asstn, file=log)
 
                print("1Q,{},{},{},{},{},{},{},{},{}".format(args.data,\
                    ncand, valid_ballots, args.quota, len(qts),\
                        qts_sample_size, straight_iqx_asn, len(winners),\
                        straight_iqx_asn))
            else:
                print("No audit possible.", file=log)

                print("1Q,{},{},{},{},{},{},{},{},{}".format(args.data, \
                    ncand,valid_ballots,args.quota,\
                    len(qts), qts_sample_size, max_sample_size, len(fws),\
                    max_sample_size))
            exit()

        if len(fws) == args.seats:

            print("Sample size required for audit is {} ballots".format(\
                max_sample_size), file=log)

            print("----------------------------------------------",file=log)
            for asstn in assertions_used:
                print(asstn, file=log)
            print("----------------------------------------------",file=log)

            print("1Q,{},{},{},{},{},{},{},{},{}".format(args.data, \
                ncand,valid_ballots,args.quota,\
                len(qts), qts_sample_size, max_sample_size, len(fws),\
                max_sample_size))
            exit()

        # Increment to use when stepping through candidate upper/lower
        # bounds on first winners' transfer values.
        deltat = args.deltat 
        deltam = args.deltam 

        max_in_outer_loop = np.inf

        best_outer_assertions = set()
        best_partial_ss = np.inf

        nfws = len(fws)

        curr_min_tvs = {fw : max(0, act_tvs[fw]-deltam) for fw in winners_on_fp}
        tv_lb_nboors = [curr_min_tvs] + create_lower_neighbours(curr_min_tvs, \
            winners_on_fp, deltam, act_tvs)

        curr_winners_verified = []
 
        outer_improved = True
        partial_improved = True

        starting_uts = None 

        while outer_improved and tv_lb_nboors != []:

            with ThreadPool(processes=args.cpus) as pool:
                func = partial(outer_loop, args, ows, fws, losers, winners,\
                    winners_on_fp, candidates, cands, valid_ballots, \
                    INVALID, ballots, max_sample_size, deltat, deltam, \
                    maxtv, act_tvs, starting_uts, straight_iqx_verified)

                results = pool.map(func, tv_lb_nboors)

            outer_improved = False
            partial_improved = False
            min_lts = np.inf
            min_partial_exc_lts = np.inf

            all_ut_asns = {f : set() for f in winners_on_fp}

            for max_in_loop, nb_min_tvs, outerdescs, lts, lts_ss, \
                best_inner_assertions, winners_verified, partial_ss,\
                partial_ss_exc_lts, ut_asns in results:

                for f,asns in ut_asns.items():
                    all_ut_asns[f].update(asns)

                min_lts = min(min_lts, lts_ss)
                min_partial_exc_lts = min(min_partial_exc_lts, \
                    partial_ss_exc_lts)

                print(outerdescs, file=log)
                if max_in_loop < max_in_outer_loop or (max_in_loop == \
                    max_in_outer_loop and partial_ss <= best_partial_ss):
                    max_in_outer_loop = max_in_loop
                    best_outer_assertions = lts
                    best_outer_assertions.update(best_inner_assertions)
                    curr_winners_verified = winners_verified
                    curr_min_tvs = nb_min_tvs

                    if partial_ss < best_partial_ss:
                        partial_improved = True

                    best_partial_ss = partial_ss
                    outer_improved = True

            if min_lts <= min_partial_exc_lts:
                break

            tv_lb_nboors = create_lower_neighbours(curr_min_tvs,\
                winners_on_fp, deltam, act_tvs)

            starting_uts = {f : None for f in winners_on_fp}

            for f in winners_on_fp:
                # Find the smallest ut whose asn is <= best partial asn
                filtered = [(utf,utf_asn) for (utf,utf_asn) in \
                    all_ut_asns[f] if utf_asn <= best_partial_ss]
                filtered.sort()

                starting_uts[f] = filtered[0][0]

        
        assertions_used.update(best_outer_assertions)

        max_sample_size = max(max_sample_size, max_in_outer_loop)

        if max_sample_size >= args.voters:
            max_sample_size = np.inf 

        if straight_iqx_asn >= args.voters:
            straight_iqx_asn = np.inf

        if straight_iqx_asn < max_sample_size:
            print("Straight IQX assertions best.", file=log)
            assertions_used = iqx_assertions
            max_sample_size = straight_iqx_asn


        if len(straight_iqx_verified) > len(curr_winners_verified):
            print("WEIRD SITUATION: more in siv than cwv")

        if [c for c in straight_iqx_verified if not 
            c in curr_winners_verified] != []:
            print("WEIRD SITUATION: c in siv not in cwv");


        s1 = set([k for k in straight_iqx_verified])
        s2 = set(curr_winners_verified)
        s3 = set([k for k in qts_verified])
        
        qts_verified_asn = max([qts_verified[w][0] for w in \
                qts_verified])

        if s2 == s3 and best_partial_ss > qts_verified_asn:
            print("Revert to straight IQ partial audit.", file=log)
            best_partial_ss = qts_verified_asn

        partial_iqx = max([straight_iqx_verified[w][0] for w in \
                straight_iqx_verified])

        if s1 == s2:
            if partial_iqx < best_partial_ss:
                print("Revert to straight IQX partial audit.", file=log)
                best_partial_ss = partial_iqx


        assertions_text = [(candidates[w].id, n, candidates[l].id \
            if l != None else None) for w,n,l,_ in assertions_used]

        # Filter out redundant assertions (eg. an A NL B when we have 
        # A AG* B already in the audit).
        
        final_assertions = []
        for assrt in assertions_used:
            if assrt[1] != "NL":
                final_assertions.append(assrt)
                continue
                   
            if (assrt[0], "AG*", assrt[2]) in assertions_text:
                continue

            final_assertions.append(assrt)

        print("Sample size required for audit is {} ballots".format(\
            max_sample_size), file=log)

        if max_sample_size != np.inf:
            print("------------------------------------------------",file=log)
            print("Assertion set BEFORE filtering:", file=log)
            for asstn in assertions_used:
                print(asstn, file=log)
            
            print("------------------------------------------------",file=log)
            print("Final set of assertions AFTER filtering:", file=log)
            for asstn in final_assertions:
                print(asstn, file=log)

 
        print("1Q,{},{},{},{},{},{},{},{},{}".format(args.data,\
            ncand, valid_ballots, args.quota, len(qts),\
            qts_sample_size, max_sample_size, len(curr_winners_verified),\
            best_partial_ss))


    log.close() 
