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


def compute_ag_stars(winners, losers, candidates, ballots, ag_matrix, \
    first_winners, min_tvs, aud_tvs, INVALID, args, log=None):

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

                if log != None:
                    print("AG*({},{}) = {}".format(cand_w.id, cand_c.id, \
                        ss), file=log)


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
    parser.add_argument('-cpus', dest='cpus', default=6, type=int)

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
    parser.add_argument('-deltat', dest='deltat', default=0.05, type=float)
    
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

        # Check that all winners that won on first preferences have a FP
        # greater than a quota.
        fws = [candidates[w] for w in winners_on_fp]
        ows = [candidates[w] for w in winners if not w in winners_on_fp]
       
        print("Winners on first preferences: {}".format([fw.id for fw in fws]),\
            file=log)
        print("Other winners: {}".format([ow.id for ow in ows]), file=log)
 
        thresh = 1.0/(args.seats + 1);
    
        act_tvs = {}
        min_tvs = {}
        for fw in fws:
            ss = ssm_sample_size(thresh, fw.fp_votes, INVALID, args)
            print("{} ballot checks required to assert that winner {} "\
                " has a quota's worth of votes".format(ss, fw.id), file=log)

            assertions_used.add((fw.num, "QT", None, ss))
            max_sample_size = max(max_sample_size, ss)

            act_tvs[fw.num] = (fw.fp_votes - args.quota)/fw.fp_votes
            min_tvs[fw.num] = 0

        if max_sample_size >= args.voters:
            max_sample_size = np.inf 

        if max_sample_size == np.inf:
            print("No audit possible.", file=log)
            exit()

        if len(fws) == args.seats:

            print("Sample size required for audit is {} ballots".format(\
                max_sample_size), file=log)

            print("----------------------------------------------",file=log)
            for asstn in assertions_used:
                print(asstn, file=log)
            print("----------------------------------------------",file=log)

            print("1Q,{},{},{},{},{}".format(args.data, ncand, valid_ballots, \
                args.quota, max_sample_size))
            exit()

        # Increment to use when stepping through candidate upper/lower
        # bounds on first winners' transfer values.
        deltat = args.deltat 

        max_in_outer_loop = np.inf

        best_outer_assertions = set()
        best_outer_exp = set()

        nfws = len(fws)

        while any([min_tvs[f] < act_tvs[f] for f in winners_on_fp]):
            print("------------------------------------------------",file=log)
            print("START OUTER LOOP, min_tvs {}".format({f.id : min_tvs[f.num]\
                for f in fws}), file=log)
            mintv_ss = 0
           
            lts = set()

            breakout = False
            for fw_i in fws:
                min_tv_i = min_tvs[fw_i.num]

                if min_tv_i > 0:
                    mintally = args.quota / (1 - min_tv_i)
                    thresh = mintally/valid_ballots

                    mintv_ss_i = ssm_sample_size(thresh, fw_i.fp_votes, \
                        INVALID, args)
                    print("Sample size to show min tv of {} is {}".format(\
                        min_tv_i, mintv_ss), file=log)

                    if mintv_ss_i >= max_in_outer_loop:
                        breakout = True
                        break

                    lts.add((fw_i.num, "LT", None, (min_tv_i, mintv_ss_i)))
                    mintv_ss = max(mintv_ss, mintv_ss_i)

            if breakout:
                break

            aud_tvs = {fw : act_tvs[fw] + deltat for fw in winners_on_fp} 

            max_in_loop = np.inf

            best_inner_assertions = None
            best_inner_explanations = None

            # MAIN LOOP OF ONE QUOTA METHOD
            while any([aud_tvs[fw] < maxtv for fw in winners_on_fp]):
                print("---------------------------------------------",file=log)
                print("START INNER LOOP", file=log)

                inner_loop_assertions = set()
                max_this_loop = 0
                max_ss_mt = 0

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

                    print("AUD TV for {}, {}, ss {}".format(fw_i.id, \
                        aud_tv_i, ss_i), file=log)

                    inner_loop_assertions.add((fw_i.num, "MT", None, \
                        (aud_tv_i, ss_i)))

                    max_ss_mt = max(max_ss_mt, ss_i)

                max_this_loop = max(max_ss_mt, max(mintv_ss, max_sample_size))

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
                    aud_tvs, INVALID, args, log=log)

                compute_ag_stars(losers, losers,  candidates, ballots, \
                    ag_matrix, winners_on_fp, min_tvs, aud_tvs, INVALID, \
                    args, log=log)

                max_with_nls = 0
                # explanations: map between AG* and the NL's that they
                # are supporting. NOTE: Explanations currently not working.
                explans = {}
                for ow_i in ows:
                    ags = {c : ag_matrix[ow_i.num][c] for c in losers \
                        if ag_matrix[ow_i.num][c] != None}
                
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
                                                if rag in explans:
                                                    explans[rag].add(curr_nl)
                                                else:
                                                    explans[rag]=set([curr_nl])
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

                            inner_loop_assertions.update(descs)

                            for ag in descs:
                                if ag in explans:
                                    explans[ag].add(curr_nl)
                                else:
                                    explans[ag] = set([curr_nl])


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
                        
                                inner_loop_assertions.update(descs)

                                for ag in descs:
                                    if ag in explans:
                                        explans[ag].add(curr_nl)
                                    else:
                                        explans[ag] = set([curr_nl])
                                
                            max_with_nls=max(max_with_nls,max(max_ags_used,ss))
                            print("NL({},{}) = {}, AG*'s used {}".format(\
                                ow_i.id, cand_c.id, ss, max_ags_used), file=log)
    
                            inner_loop_assertions.add((ow_i.num, "NL", c, \
                                max(ss,max_ags_used)))
                        else:
                            max_with_nls = np.inf
                            print("NL({},{}) NOT POSSIBLE".format(\
                                ow_i.id, cand_c.id),file=log)

                max_this_loop=max(max_this_loop,max(max_ss_mt,max_with_nls))

                print("Max in loop {}, {}".format(max_this_loop, {f.id : \
                    aud_tvs[f.num] for f in fws}), file=log)

                if max_this_loop <= max_in_loop:
                    best_inner_assertions = inner_loop_assertions
                    best_inner_explanations = explans
                    max_in_loop = max_this_loop

                    if max_with_nls >= max_ss_mt:
                        # Contuining in loop will not improve overall
                        # audit.
                        print("Breaking out of inner loop {} > {}".format(\
                            max_this_loop, max_in_loop), file=log)
                        break


                    # Add deltat to the smallest of aud_tvs
                    fmin, fminv = None, None
                    for f in winners_on_fp:
                        if fmin == None:
                            fmin = f
                            fminv = aud_tvs[f]
                        else:
                            v = aud_tvs[f]
                            if v < fminv:
                                fmin = f
                                fminv = v

                    aud_tvs[fmin] += deltat
                else:
                    print("Breaking out of inner loop {} > {}".format(\
                        max_this_loop, max_in_loop), file=log)
                    break

            
            if max_in_loop <= max_in_outer_loop:
                max_in_outer_loop = max_in_loop
                best_outer_assertions = lts
                best_outer_assertions.update(best_inner_assertions)
                best_outer_exp = best_inner_explanations
            else:
                print("Breaking out of outer loop", file=log)
                break    

            # If there is a min_tv for a candidate that is 0, set it to
            # act_tv/2 for that candidate
            zeros = [f for f in winners_on_fp if min_tvs[f] == 0]
            if zeros != []:
                f = zeros[0]
                min_tvs[f] = act_tvs[f]/2.0
            else:
                fmin, fminv = None, None
                for f in winners_on_fp:
                    if fmin == None:
                        fmin = f
                        fminv = min_tvs[f]
                    else:
                        v = min_tvs[f]
                        if v < fminv:
                            fmin = f
                            fminv = v

                min_tvs[fmin] += deltat
                print("New mintv, {} = {}".format(candidates[fmin].id, \
                    min_tvs[fmin]),file=log)

        assertions_used.update(best_outer_assertions)
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

        for asstn in best_outer_exp:
            filtered = []
            for exp in best_outer_exp[asstn]:
                if exp in final_assertions:
                    filtered.append(exp)

            best_outer_exp[asstn] = filtered


        max_sample_size = max(max_sample_size, max_in_outer_loop)

        if max_sample_size >= args.voters:
            max_sample_size = np.inf 

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

            print("------------------------------------------------",file=log)
            print("Explanations for AG* inclusion:", file=log)
            for asstn in final_assertions:
                if asstn in best_outer_exp:
                    print("   Explanation for {}: {}".format(asstn, \
                        best_outer_exp[asstn]), file=log)
            print("------------------------------------------------",file=log)

        print("1Q,{},{},{},{},{}".format(args.data, ncand, valid_ballots, \
            args.quota, max_sample_size))


    log.close() 
