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



# Determine if there is a reported loser 'd' for which 'd' AG loser, based 
# on the AG relationships in ag_matrix, preferenced before 'loser' in
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
# where ag_present is a boolean indicating whether we found a 
# 'd' such that 'd' AG loser, and 'd' is preferenced before loser
# in prefs.
#
# If loser is not present in the ranking 'prefs' the function will
# return:
#
#     False, np.inf
#
# Note that the set winners contains the set of reported winners for the 
# contest.
def rule_out_for_max(prefs, loser, ag_matrix, nl_matrix, winners, candidates):
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
        nl,nl_supp = nl_matrix[p][loser] if nl_matrix != None else (None,set())
        if ag != None and ag != np.inf:
            ag_present = True
            ags_used.add((p,"AG",loser, ag))
            ag_min_ss = min(ag_min_ss, ag)

        elif nl != None and nl != np.inf:
            ag_present = True
            ags_used.update(nl_supp)
            ag_min_ss = min(ag_min_ss, nl)

    return False, np.inf, set()
        

# Merge a sequence of AG relationships that could be used to increase the
# assorter margin of an NL.
#
# Input list 'helpful_ags' will be a list of (asn, extra, desc),
# where 'asn' is the cost of auditing the AG assertion, 'extra' is 
# the increase to the assorter total if we incorporate those AGs, and 'desc'
# is a set of textual descriptions of the assertions. This 
# function takes a list of these AGs, and merges consecutive entries if the
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


def get_max_ballot_weight(prefs, first_winners, aud_tv, elected, maxtv,\
    eliminated, standing):

    weights = []
    fw_prefix = []

    if prefs[0] in elected:
        return maxtv, [] # Reducing this will be helpful! 

    for p in prefs:
        if p in eliminated:
            continue

        if p in standing:
            return 0, []

        if p in first_winners:
            weights.append(aud_tv[p])
            fw_prefix.append(p) 
        else:
            break

    if weights == []:
        return 1, []
    else:
        return weights[0], fw_prefix


def establish_max_tvalue_nonfw(w, candidates, ballots, winners_on_fp, \
    min_tvs, aud_tvs, ag_matrix, maxtv, INVALID, args, curr_asn):
    
    info = " curr_asn = {}\n".format(curr_asn)

    # Verify that 'w' did not receive a quota on first preference NIQ(w)
    # We need to do this as this approach of bounding w's transfer value
    # relies on them receiving a quota by the distribution of votes from
    # another candidate (either via a surplus or elimination).
    cand_w = candidates[w]

    if cand_w.fp_votes >= args.quota:
        info = "Cannot form NIQ({}) as they have a quota on FPs\n".format(\
            cand_w.id)
        return maxtv, 0, set(), info

    # Show that total of everyone ELSE's fp votes is more than seats-1 quotas.
    thresh = 1 - (1.0/(args.seats + 1));
    
    valid_ballots = args.voters - INVALID
    others_total = valid_ballots - cand_w.fp_votes

    asn = ssm_sample_size(thresh, others_total, INVALID, args)

    if asn >= args.voters:
        info = "Cannot form NIQ({}); too expensive.\n".format(cand_w.id)
        return maxtv, 0, set(), info
        
    # Find how much each candidate would transfer to w on their elimination
    # (if they are not a candidate in winners_on_fp, or election otherwise).
    max_surpluses = {}
    max_overall = 0

    tot_vote = 0

    # Assuming a context where each ballot may have a value < 1 
    weighted_ballots = []
    for b in ballots:
        if b.prefs == []:
            continue

        p = b.prefs[0]
        weight = aud_tvs[p] if p in winners_on_fp else 1

        tot_vote += weight*b.votes

        weighted_ballots.append((b, weight))

    for cand in candidates:
        c = cand.num

        if c == w:
            continue
            
        # Use aud_tvs as maximum transfer value for any ballot that 
        # has a candidate in winners_on_fp as first preference and 1 
        # otherwise.
        transfer = 0
        tot_not_t = 0

        # Record contribution of ballot type to "Not transfer" category, to
        # total vote, and number of ballots of that type.
        details = []

        for (b,v) in weighted_ballots:
            if b.prefs == []:
                continue

            win = b.prefs.index(w) if w in b.prefs else np.inf
            cin = b.prefs.index(c) if c in b.prefs else np.inf

            if w in b.prefs and c in b.prefs:
                widx = b.prefs.index(w) 
                cidx = b.prefs.index(c) 

                if cidx < widx:
                    details.append((0, v, b.votes)) 
                    transfer += v * b.votes
                else:
                    details.append((v, v, b.votes))
                    tot_not_t += v * b.votes

            else:
                details.append((v, v, b.votes))
                tot_not_t += v * b.votes


        max_surpluses[c] = (transfer, tot_not_t, details)
        max_overall = max(max_overall, transfer)
           
    # Try and establish a max_tr_tv that would cost less than or equal
    # to 'curr_asn' to verify in an audit.
    delta = 50
    max_tr_tv_ss = np.inf
    bar = max_overall
    aset = set()
    desc = ""

    while max_tr_tv_ss > max(asn, curr_asn) and bar + delta < tot_vote:

        bar += delta

        # should we use total original vote or total vote in the context
        # we are interested in? (ie. with first winners seated and their
        # transfer values set to their maximum value).
        thresh = (tot_vote - bar)/tot_vote

        # We want to show that the bucket of votes in the category "will
        # transfer to w in the event of election/elimination" is less than
        # bar.

        max_tr_tv_ss = 0
        aset = set()

        desc += "MAX SURPLUS ASSERTIONS:\n"

        for c,(s,tnt,details) in max_surpluses.items():
            ss = ssm_sample_size(thresh, tnt, INVALID, args)

            # We want to show that p_s < thresh
            # so p_tnt > t where t = 1 - thresh
            #amean = 0.5*INVALID

            #for cont_tnt, cont_v, num in details:
            #    amean += num*(((cont_tnt - thresh*cont_v) + thresh)/(2*thresh))

            #amean /= args.voters
            
            #ss = sample_size(amean, args)
  
            max_tr_tv_ss = max(ss, max_tr_tv_ss)

            aset.add((c, "SP[{}]".format(bar), None, ss))
            desc += "SP[{},{}]({}) with ASN {}\n".format(bar, s, \
                candidates[c].id, ss)

    max_tr_tv_ss = max(asn, max_tr_tv_ss)
        
    max_tr_tv = maxtv
    if max_tr_tv_ss <= curr_asn:
        max_tr_tv = (bar/(bar+args.quota))
        info += desc
        info += "Reduce maxtv for {} to {} ASN {}\n".format(cand_w.id, \
            max_tr_tv, max_tr_tv_ss)

    if max_tr_tv < maxtv: 
        return max_tr_tv, max_tr_tv_ss, aset, info
 
    # Return established maximum on transfer value, asn to verify this,
    # set of assertions required to establish this, logging info
    return maxtv, 0, set(), info

def form_NL(candidates, c, ow_i, ags, nls, ballots, INVALID, winners_on_fp, \
    min_tvs, aud_tvs, ag_matrix, nl_matrix, winners, maxtv, elected, \
    eliminated,standing):
    
    cand_c = candidates[c]

    min_ow_i = 0 # Min tally for ow_i.
    max_c = 0 # Max tally for candidate c.

    # Keep running tally of total votes we can increase the
    # margin of assertion 'ow_i' NL 'c' with 'AG*' 
    # relationships
    pot_margin_inc = 0

    helpful_ags = []
    nl_assertion_set = set()

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

        # Reducing upper bound on 'elected's transfer values down
        # from maxtv will likely be helpful.
        weight, fwprefix = get_max_ballot_weight(b.prefs,\
            winners_on_fp, aud_tvs, elected, maxtv, \
            eliminated, standing)

        c_idx=b.prefs.index(c) if c_in else np.inf
        o_idx=b.prefs.index(ow_i.num) if ow_in else np.inf

        if c_in and c_idx < o_idx:
            # These ballots could not go to 'ow' but may go to 'c'.
            if weight == 0:
                # There is someone still assumed to be standing before 'c'
                assorter += 0.5*b.votes
            else:
                # Could we exclude ballots by using an AG*?
                is_ag, ag_asn, descs = rule_out_for_max(\
                    b.prefs, c, ag_matrix, nl_matrix, winners, candidates)
                    
                contrib = b.votes*((1-weight)/2.0)

                if is_ag:
                    alt_contrib = b.votes*0.5
                    helpful_ags.append((ag_asn, alt_contrib - contrib, descs))
                    pot_margin_inc += alt_contrib-contrib

                assorter += contrib
                max_c += b.votes

        elif ow_in:
            # Check if all candidates before ow_i are
            # first winners.
            exprefs = [p for p in b.prefs if not p in fwprefix and \
                not p in eliminated]

            if exprefs[0] == ow_i.num and fwprefix == []:
                # Only candidates before ow_i are assumed to be eliminated
                min_ow_i += b.votes
                assorter += b.votes
             
            elif exprefs[0] == ow_i.num and b.prefs[0] in fwprefix:
                minval = min_tvs[fwprefix[0]]
                min_ow_i += minval * b.votes
                assorter += b.votes * ((1 + minval)/2.0)
            else:
                # If we remove all cand 'd' for which ow AG*
                # d that appear before ow in the ballot,
                # will 'ow' be the first ranked cand? If so,
                # we could add these ballots to the min
                # tally of 'ow'.
                prefs = [p for p in b.prefs if not p in eliminated\
                    and not p in winners_on_fp]

                descs = set()
                max_ags_here = 0

                minval = 1
                if b.prefs != [] and b.prefs[0] in winners_on_fp:
                    minval = min_tvs[b.prefs[0]]

                other_winners = [w for w in winners if not w in winners_on_fp]
                skip_winner = False
                                   
                for d,dval in ags.items():
                    if d in prefs:
                        idx_d = prefs.index(d)
                        if idx_d < o_idx:
                            prefs.remove(d)
                            o_idx -= 1
                            rag = (ow_i.num,"AG*",d,dval)
                            descs.add(rag)
                            max_ags_here=max(max_ags_here,dval)
                            if d in other_winners:
                                skip_winner = True
                                break


                for d,(dval,dset) in nls.items():
                    if d in prefs:
                        idx_d = prefs.index(d)
                        if idx_d < o_idx:
                            prefs.remove(d)
                            o_idx -= 1
                            descs.update(dset)
                            max_ags_here=max(max_ags_here,dval)
                            if d in other_winners:
                                skip_winner = True
                                break

                assorter += 0.5*b.votes

                if not skip_winner:
                    if prefs != [] and prefs[0] == ow_i.num:
                        base_contrib = 0.5*b.votes
                        alt_contrib = b.votes * ((1 + minval)/2.0)
                        dconfig = alt_contrib-base_contrib

                        if dconfig > 0:
                            helpful_ags.append((max_ags_here, dconfig, descs))

                            pot_margin_inc += dconfig
                    #else:
                    #    fwprefix = []
                    #    minws = {}
                    #    for p in prefs:
                    #        if p in winners_on_fp:
                    #            fwprefix.append(p)
                    #            minws[p] = min_tvs[p]
                    #        else:
                    #            break

                    #    if fwprefix != [] and prefs[len(fwprefix)] == ow_i.num:
                    #        minval = minws[b.prefs[0]] if b.prefs[0]\
                    #            in winners_on_fp else 1

                    #        base_contrib = 0.5*b.votes
                    #        alt_contrib = b.votes * ((1 + minval)/2.0)
                    #        dconfig = alt_contrib-base_contrib

                    #        if dconfig > 0:
                    #            helpful_ags.append((max_ags_here, dconfig, descs))

                    #            pot_margin_inc += dconfig
        else:
            assorter += 0.5*b.votes

    # Max ASN of any AG*'s used to increase assorter 
    # margins when forming NLs.
    max_ags_used = 0  
    merged_helpful_ags = merge_helpful_ags(helpful_ags, pot_margin_inc)

    # Incorporate use of  AG*'s that either make the
    # assertion possible, or whose ASN is already
    # within/equal to current lower bound on audit
    # difficulty.
    while assorter/args.voters <= 0.5 and merged_helpful_ags != []:
        ag_asn, extra_contrib, descs = merged_helpful_ags.pop(0)

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


        asn = max(ss, max_ags_used)

        if asn >= args.voters:
            return False, np.inf, None, "NL would require full hand count"

        desc = "NL({},{}) = {}, AG*'s used {}\n".format(\
            ow_i.id, cand_c.id, ss, max_ags_used)

        nl_assertion_set.add((ow_i.num, "NL", c, \
            max(ss,max_ags_used)))

        return True, max(max_ags_used,ss), nl_assertion_set, desc
    else:
        desc = "NL({},{}) NOT POSSIBLE, AMEAN {}\n".format(\
            ow_i.id, cand_c.id, amean)

        return False, np.inf, None, desc


def compute_iqx(candidates, ow_i, ag_matrix, ballots, INVALID, args, min_tvs,
    winners_on_fp, winners):

    desc = ""
    qthresh = 1.0/(args.seats + 1);

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

    winners_not_fp = [w for w in winners if not w in winners_on_fp]

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

                    if d in winners_not_fp:
                        # value of ballot uncertain
                        value = 0
                        break

        if value > 0 and prefs[0] == ow_i.num:
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
        desc = "Can form IQX({}) with sample size {}, AG*s {}\n".format(\
            ow_i.id, ss, ss_ag_iqx) 
        ss_iqx = max(ss, ss_ag_iqx)

        return True, ss_iqx, iqx_assertions, desc

    return False, ss_iqx, set(), desc
    

def form_SB(nv, ow, candidates, ballots, winners_on_fp, winners, args, INVALID):
    max_v_ow = 0

    desc = ""

    # These are winners that while they may have won on first preferences,
    # we have not demonstrated that they have done so. Excluding nv and ow.
    winners_not_on_fp = [w for w in winners if not w in winners_on_fp \
        and w != nv and w != ow]

    for b in ballots:
        if b.prefs == []:
            continue

        ow_in = ow in b.prefs

        if not ow_in:
            continue

        nv_in = nv in b.prefs

        value = 1

        if b.prefs[0] in winners_on_fp:
            value = min_tvs[b.prefs[0]]

        if ow_in and not nv_in:
            max_v_ow += value*b.votes

        else:
            ow_idx = b.prefs.index(ow)
            nv_idx = b.prefs.index(nv)
        
            if ow_in < nv_in:
                max_v_ow += value*b.votes

            else:
                # could ballot have skipped over 'nv'?
                for wnfp in winners_not_on_fp:
                    if wnfp in b.prefs:
                        idx = b.prefs.index(wnfp)
                        if idx < ow_idx:
                            max_v_ow += value*b.votes
                            break
                            

#    if max_v_ow < args.quota:
        #todo
#    else:
#        desc += "Cannot show that {} SB {}.\n".format(candidates[nv].id,\
#            candidates[ow].id)

#        return None, set(), desc        



def inner_loop(winners_on_fp, args, candidates, cands, valid_ballots,\
    INVALID,max_sample_size,mintv_ss,ballots,min_tvs,ows,fws,losers,winners, \
    straight_iqx_verified, runner_up, aud_tvs):

    np.seterr(all='ignore')

    ag_matrix = [[None for c in candidates] for o in candidates]

    desc = "---------------------------------------------\n"
    desc += "START INNER LOOP\n"

    inner_loop_assertions = set()
    max_this_loop = 0
    max_ss_mt = 0

    winners_verified = winners_on_fp[:]

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

        desc += "AUD TV for {}, {}, ss {}\n".format(fw_i.id, \
            aud_tv_i, ss_i)

        inner_loop_assertions.add((fw_i.num, "MT", None, \
            (aud_tv_i, ss_i)))

        max_ss_mt = max(max_ss_mt, ss_i)

    partial_ss = 0

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

    newly_verified = []
    could_not_verify = []

    nl_matrix = [[(None,set()) for c in candidates] for o in candidates]

    for ow_i in ows:
        ags = {c : ag_matrix[ow_i.num][c] for c in losers \
            if ag_matrix[ow_i.num][c] != None}
               
        success, ss_iqx, iqx_assertions, info = compute_iqx(candidates, \
            ow_i, ag_matrix, ballots, INVALID, args, min_tvs, winners_on_fp,\
            winners)

        desc += "{}-{}\n".format(success, ss_iqx)
        desc += info
            

        nl_assertion_set = set()
        max_with_nls_i = 0

        # Determine NL's between original losers and winner ow_i
        for c in cands:
            if c in winners:
                continue

            successNL, asn, aset, info = form_NL(candidates, c, ow_i, ags, \
                {}, ballots,INVALID,winners_on_fp, min_tvs,aud_tvs,ag_matrix,\
                None, winners, maxtv, [], [], [])

            desc += info

            max_with_nls_i = max(max_with_nls_i, asn)

            if successNL:
                nl_assertion_set.update(aset)
                nl_matrix[ow_i.num][c] = (asn, aset)


        straight_iqx_asn = straight_iqx_verified[ow_i.num][0] if \
            ow_i.num in straight_iqx_verified else np.inf

        if ss_iqx < args.voters or max_with_nls_i < args.voters or \
            straight_iqx_asn < args.voters:
            winners_verified.append(ow_i.num)
            newly_verified.append(ow_i.num)
            partial_ss = max(partial_ss, min([ss_iqx, max_with_nls_i,\
                straight_iqx_asn])) 
            desc += "{}-{}-{}\n".format(ss_iqx,max_with_nls_i,straight_iqx_asn)

            if straight_iqx_asn < ss_iqx and straight_iqx_asn < max_with_nls_i:
                desc += "CHOOSE STRAIGHT IQX\n"
                max_this_loop = max(max_this_loop, straight_iqx_asn)
                partial_ss = max(partial_ss, straight_iqx_asn)
                inner_loop_assertions.update(straight_iqx_verified[ow_i.num][1])
        
            elif ss_iqx < max_with_nls_i:
                desc += "CHOOSE IQX over NLs\n"
                max_this_loop = max(max_this_loop, ss_iqx)
                partial_ss = max(partial_ss, ss_iqx)
                inner_loop_assertions.update(iqx_assertions)

            else:
                max_this_loop = max(max_this_loop, max_with_nls_i)
                inner_loop_assertions.update(nl_assertion_set)
           
        else:
            could_not_verify.append(ow_i.num)

    nv_maxtvs = {}

    ub = max(partial_ss, max_sample_size)
    outer_asn = max(mintv_ss, max_ss_mt)

    if outer_asn < args.voters:
        ub = max(ub, outer_asn)

    for nv in newly_verified:
        max_tr_tv, max_tr_tv_asn, max_tr_tv_aset, max_tr_tv_info = \
            establish_max_tvalue_nonfw(nv, candidates, ballots, winners_on_fp,\
            min_tvs, aud_tvs, ag_matrix, maxtv, INVALID, args, ub)

        desc += max_tr_tv_info

        nv_maxtvs[nv] = (max_tr_tv, max_tr_tv_asn, max_tr_tv_aset)

    if could_not_verify != [] and newly_verified != []:
        max_this_loop = partial_ss

        # Can we show that any of the candidates in could_not_verify
        # cannot achieve a quota before one of the candidates in
        # newly_verified is seated?
        #for ow in could_not_verify:
        #    for nv in newly_verified:

                # Compute max vote of ow in context where 
                # nv is still standing.
        #        form_SB(nv, ow, candidates, ballots, )

        improvement = True
        while improvement and could_not_verify != []:

            improvement = False
            new_could_not_verify = []

            for ow in could_not_verify:
                successNL = True
                asn = 0
                aset = set()

                nls = {c : nl_matrix[ow][c] for c in losers \
                    if nl_matrix[ow][c][0] != None}

                for c in losers:
                    successNLc = False

                    if nl_matrix[ow][c][0] != None:
                        continue
                
                    for nv in newly_verified:
                        # Can we generate NL under two assumptions about the
                        # candidate nv? ie. that they have already been elected, 
                        # or are still standing?

                        mtv, mtv_asn, mtv_aset = nv_maxtvs[nv]

                        # For the first case, it will likely help if we can reduce
                        # maxtv for 'nv'
                        s2, asn2, aset2, info2 = form_NL(candidates, c, ow_i, ags, \
                            nls, ballots, INVALID, winners_on_fp, min_tvs, aud_tvs,\
                            ag_matrix, nl_matrix, winners, mtv, [nv], [], []) # elected
                        s3, asn3, aset3, info3 = form_NL(candidates, c, ow_i, ags, \
                            nls, ballots, INVALID, winners_on_fp, min_tvs, aud_tvs, \
                            ag_matrix, nl_matrix, winners, mtv, [], [], [nv]) # standing

                        if s2 and s3:
                            desc += "Can show {} NL-R[{}] {}\n".format(ow_i.id, \
                                candidates[nv].id, candidates[c].id)

                            successNLc = True
                            asn = max(mtv_asn, max(asn2, asn3))
                            aset = aset2
                            aset.update(aset3) 
                            aset.update(mtv_aset)

                            desc += info2
                            desc += info3

                            nl_matrix[ow][c] = (asn, set(aset))
                            improvement = True


                        else:
                            desc += "CANNOT show {} NL-R[{}] {}\n".format(\
                                ow_i.id, candidates[nv].id, candidates[c].id)
                            desc += info2
                            desc += info3

                        if successNLc:
                            break

                    if not successNLc:
                        successNL = False
                        break


                if successNL:
                    winners_verified.append(ow)
                    newly_verified.append(ow)
                    max_this_loop = max(max_this_loop, asn)
                    partial_ss = max(partial_ss, asn)
                    inner_loop_assertions.update(aset)

                else:
                    new_could_not_verify.append(ow)        
       
            could_not_verify = new_could_not_verify[:]

    if could_not_verify != []:
        max_this_loop = np.inf

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
    deltam, maxtv, act_tvs, starting_uts, straight_iqx_verified, \
    runner_up, min_tvs):
   
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
                mintv_ss, ballots, min_tvs, ows, fws, losers, winners,\
                straight_iqx_verified, runner_up)

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

    runner_up = None

    seats = 0
    for i in range(ncand):
        c = outcome.cand[i]
        cands.append(c)
        cd = candidates[c]

        if outcome.action[i]:
            winners.append(c)
            if cd.fp_votes > args.quota:
                winners_on_fp.append(c)
            seats += 1
        else:
            if seats < args.seats:
                runner_up = c
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
                    maxtv, act_tvs, starting_uts, straight_iqx_verified,\
                    runner_up)

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

            if len(s3) == args.seats:
                max_sample_size = best_partial_ss

        partial_iqx = max([straight_iqx_verified[w][0] for w in \
                straight_iqx_verified])

        if s1 == s2:
            if partial_iqx < best_partial_ss:
                print("Revert to IQX partial audit.", file=log)
                best_partial_ss = partial_iqx
           
                if len(s1) == args.seats:
                    max_sample_size = best_partial_ss


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
