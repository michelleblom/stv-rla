#
#    Copyright (C) 2024  Michelle Blom
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
import itertools

from functools import partial

from multiprocessing import Pool

from utils import read_outcome, sample_size, ssm_sample_size,\
    tally_vs_tally_sample_size, read_ballots_txt, index_of, next_cand, \
    read_ballots_json, read_ballots_txt, read_ballots_blt, \
    read_ballots_stv, simulate_stv


# Returns 'c1' if in the ballot ranking 'prefs', 'c1' appears earlier 
# than 'c2', or if 'c1' appears and 'c2' does not. Returns 'c2' if 
# that candidates appears before 'c1', or 'c2' appears and 'c1' does not.
# Returns 'None' if neither candidates is mentioned in the ranking. 
def vote_for_cand_over_cand(c1, c2, prefs):
    for p in prefs:
        if p == c1:
            return c1

        if p == c2:
            return c2

    return None


# This function is used to determine whether a ballot should contribute
# to a lower bound on the tally of a candidate 'c' when it is acting as the
# 'winner' in a NL assertion.
#
# Determines if candidate 'c' is the first ranked candidate in the ranking
# 'prefs' once we remove all candidates 'd' for which 'c' AG 'd'. These
# AGs are stored in the map 'c_ag'. 
#
#     c_ag[d] gives the ASN of the AG ('c' AG 'd')
#
# If we use any AGs in the determination of whether the ballot should be
# awarded to 'c', then we return the maximum ASN of those AGs. If the 
# provided threshold is not None, we do not make use of any AG assertions that
# have an ASN of greater than the threshold.
#
# If the ballot should be awarded to 'c', we return True,Max ASN of AGs used.
#
# If the ballot should not be awarded to 'c', we return False,None.
def vote_for_cand_ags1(c, prefs, c_ag, candidates):
    ags_used = 0

    descs = set()

    for p in prefs:
        if p == c:
            return True,ags_used,descs

        if p != c:
            if p in c_ag:
                ags_used = max(ags_used, c_ag[p])
                descs.add((c,"AG",p,c_ag[p]))
            else:
                return False, None, set()

    return False, None, set()


# This function is used to determine whether a ballot should contribute
# to the upper bound on the tally of a candidate 'c1' when it is acting as the
# 'loser' in a NL assertion with winner 'c2'.
#
# The ballot in question has the ranking 'prefs'.
#
# We assume that the candidates in 'winners' are seated at some point.
#
# Candidate c1 will not get this ballot if:
# -- It does not appear (returns False,None)
# -- Candidate 'c2' appears before 'c1'
# -- A candidate 'd', for which 'd' AG 'c1' appears before 'c1'. If we 
#    use such an AG to determine that this ballot does not go to 'c1', 
#    we keep track of the maximum AG assertion ASN used. 
#
# Note, we do not use an AG assertion with an ASN above a certain threshold
# when determining if we should count ballots with ranking 'prefs' in the
# maximum tally of 'c1'.
#
# The input 'ag_matrix' contains AG relationships between candidates.
#
#      ag_matrix[d][c1] gives the ASN of the AG assertion 'd' AG 'c1'
#          if it exists (ag_matrix[d][c1] will be None otherwise).
#
# Returns False,None or False,Max ASN of AGs used if the ballot should
# not be counted toward the maximum tally of 'c1' in the NL assertion
# 'c2' NL 'c1'.
#
# Return True,_ if the ballot should be counted toward 'c1's max tally.
def vote_for_cand_ags2(c1, c2, prefs, ag_matrix, winners, candidates):
    ag_present = False
    ag_min_ss = np.inf

    descs = set()
    for p in prefs:
        if p == c1:
            return True, ag_present, ag_min_ss, descs

        if p == c2:
            return False, False, None, set()

        if p in winners:
            # c1 could still get this vote as it could skip over 'p'
            continue

        ag = ag_matrix[p][c1]
        if ag != None and ag != np.inf:
            ag_present = True
            ag_min_ss = min(ag_min_ss, ag)

            descs.add((p, "AG", c1, ag))

    return False, False, None, set()



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
    first_winner, min_tv, aud_tv, INVALID, args, log=None):

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
                    if b.prefs[0] == cand_w.num:
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

                if min_tv > 0 and len(b.prefs) > 1 and b.prefs[:2] == \
                    [first_winner.num, cand_w.num]:
                    min_w += min_tv * b.votes
                    contrib = b.votes * ((1 + min_tv)/2.0)
                else:
                    weight = 1

                    if b.prefs != [] and b.prefs[0] == first_winner.num:
                        weight = aud_tv

                    for p in b.prefs:
                        if p == cand_w.num:
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


def rule_out_n(candidates, cands, ruled_out, args, INVALID, ag_matrix, \
    ballots, tvmax, nset):

    np.seterr(all='ignore')

    cset = [candidates[c] for c in nset]
    idcset = [cand.id for cand in cset]
    info = "ALTERNATE OUTCOME SET: {}\n".format(idcset)

    # ==============================================================
    # No assumption around order in which candidates
    # are elected required.
    # --------------------------------------------------------------
    # For each subset of cset of size n-1, subcset, if we assume 
    # they get a seat at some stage, does there exist a candidate o 
    # \notin cset that cannot be eliminated before the one 
    # candidate 'c' that is in cset but not in subcset? If so, 
    # this outcome -- cset at the set of winners -- cannot be true.

    failure = True

    N = len(cset)
    others = [o for o in cands if not (o in nset or o in ruled_out)]

    best_asn = np.inf
    best_assertions = set()

    for c in nset:
        cand_c = candidates[c]
        subcset = [can for can in nset if can != c]

        # Assuming subcset is seated.
        for o in others:
            o_assertions = set()
            o_asn = np.inf
            max_ags_used = 0

            # Find set of candididates 'd' for which o AG d.
            # Keep track of cost of each of those AGs
            cand_o = candidates[o]
            o_ag = {}
            for d in others:
                if d != o and ag_matrix[o][d] != None and \
                ag_matrix[o][d] != np.inf:
                    o_ag[d] = ag_matrix[o][d]

            # Consider max vote of 'c' given subcset seated at some stage and 
            # 'o' is still standing.
            assorter = INVALID*0.5

            helpful_ags = []
            pot_margin_inc = 0

            for b in ballots:
                awarded, used_ss, descs = vote_for_cand_ags1(o, b.prefs,\
                    o_ag, candidates)

                if awarded:
                    if used_ss != None and used_ss > 0:
                        dconfig = b.votes*0.5
                        helpful_ags.append((used_ss, dconfig, descs))
                        pot_margin_inc += dconfig
                        assorter += b.votes*0.5
                        continue 
                        
                    assorter += b.votes
                    continue    
                        
                if not c in b.prefs: 
                    assorter += b.votes*0.5
                    continue

                weight = 1.0

            if b.prefs != [] and b.prefs[0] in subcset:
                weight = tvmax

            # In this analysis, we are assuming that no one other
            # than a candidate in cset has won -- so if AG(d, c) and 'd'
            # is not in cset, and 'd' is preferenced before c, then
            # c does not get these ballots. 
            awarded, ag_present, used_ss, descs = vote_for_cand_ags2(\
                c, o, b.prefs, ag_matrix, nset, candidates)

            if awarded:
                contrib = b.votes*((1 - weight)/2.0)
                if ag_present:
                    altcontrib = b.votes*0.5
                    dconfig = altcontrib - contrib
                    helpful_ags.append((used_ss, dconfig, descs))
                    pot_margin_inc += dconfig

                contrib = b.votes*((1 - weight)/2.0)
                assorter += contrib
            else:
                assorter += b.votes*0.5

            amean = assorter/args.voters         
            merged_helpful_ags=merge_helpful_ags(helpful_ags,pot_margin_inc)

            G = set()
            O = set()

            # Incorporate use of  AG's that either make the assertion
            # possible, or whose ASN is already within/equal to current
            # lower bound on audit difficulty.
            while amean <= 0.5 and merged_helpful_ags != []:
                ag_asn, extra_contrib, descs = merged_helpful_ags.pop(0)

                assorter += extra_contrib
                max_ags_used = max(max_ags_used, ag_asn)
                o_assertions.update(descs)

                for w,_,l,_ in descs:
                    if w == o:
                        G.add(l)

                    elif l == c:
                        O.add(w)

                amean = assorter/args.voters

            
            info += "   (1) can we show that {} NL {}?\n".format(cand_o.id,\
                cand_c.id)
            info += "        a. margin {}\n".format(2*amean - 1)

            # Is the minimum tally for 'o' larger than the maximum
            # possible tally for 'c'? This means 'o' cannot be 
            # eliminated before 'c'
            if amean > 0.5:
                ss = sample_size(amean, args)

                # Go through remaining AG's 
                for ag_asn, extra_contrib, descs in merged_helpful_ags:
                    if ss < ag_asn:
                        break

                    # would reducing margin by votes reduce sample size
                    # by more than increase caused by including AG*?
                    assorter += extra_contrib
                    amean = assorter/args.voters

                    ss  = sample_size(amean, args)
                    max_ags_used = max(max_ags_used, ag_asn)
                        
                    o_assertions.update(descs)

                    for w,_,l,_ in descs:
                        if w == o:
                            G.add(l)

                        elif l == c:
                            O.add(w)


                failure = False
                o_asn = max(ss, max_ags_used)
                o_assertions.add(\
                    "{} NL {} = {} [rules out {}, G = {}, O = {}]\n".format(\
                    cand_o.id, cand_c.id, o_asn, idcset, G, O))

                if ss != np.inf:
                    info += "      yes, {}, ags used {} = {}\n".format(ss, \
                         max_ags_used, o_asn)
                else:
                    info += "      no\n"
            else:
                info += "      no\n"

            
            if o_asn < best_asn: 
                best_assertions = o_assertions
                best_asn = o_asn

        if not failure:
            info += "      Best assertion(s) require sample: {}\n".format(\
                best_asn)

    return failure, nset, best_asn, best_assertions, info


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    # Input: stv data file
    parser.add_argument('-d', dest='data')

    # Input: anticipated error rate (default value is 0)
    parser.add_argument('-e1', dest='erate1', default=0.002, type=float)
    parser.add_argument('-e2', dest='erate2', default=0, type=float)

    # Input: risk limit (default is 10%)
    parser.add_argument('-r', dest='rlimit', default=0.1, type=float)
    
    # Input: number of repetitions to perform in simulation to determine
    # an initial sample size estimation -- the quantile of the sample
    # size is computed (default is 1 repetition -- no error rate)
    parser.add_argument('-reps', dest='reps', default=20)

    # Input: seed (default is 9368663)
    parser.add_argument('-s', dest='seed', default=9368663, type=int)
    
    # Input: number of cpus to use if parallelising tasks.
    parser.add_argument('-cpus', dest='cpus', default=6, type=int)

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
    losers = [] # Reported losers

    for i in range(ncand):
        c = outcome.cand[i]
        cands.append(c)
        cd = candidates[c]

        if outcome.action[i]:
            winners.append(c)
        else:
            losers.append(c)

    winnerset = set(winners)

    nseats = len(winners)
    maxtv = nseats/(nseats + 1)

    print("VALID BALLOTS {}".format(valid_ballots), file=log)

    # This matrix will be used to store AG relationships between candidates.
    ag_matrix = [[None for c in cands] for o in cands]

    # ag_matrix[c][o] will give you the sample size required to show
    # that c AG o. It will be None if no such AG relationship exists.
    rule_out = []
    max_rule_out_ss = 0
    ruled_out = []

    # possibilities[c] will give you all the sample sizes of 
    # AG relationships where 'c' is the loser.
    possibilities = [[] for c in cands]

    print("Computing AG relationships", file=log)
    # Compute "Always greater" AG relationships.
    for c in cands:
        # Cand 'c' is always greater than 'o' if c's first preference
        # tally is greater than the maximum tally for 'o'
        cand_c = candidates[c]
        fpc = cand_c.fp_votes

        for o in cands:
            if c == o: continue

            if ag_matrix[c][o] != None:
                continue

            # The maximum tally for 'o' is equal to the number of
            # ballots that preference 'o' before 'c' or preference
            # 'o' and not 'c'. 
            max_vote_o = 0
            cand_o = candidates[o]

            for b in ballots:
                awarded = vote_for_cand_over_cand(c, o, b.prefs)

                if awarded == o:
                    max_vote_o += b.votes

            if fpc > max_vote_o:
                ss  = tally_vs_tally_sample_size(fpc, max_vote_o,\
                    valid_ballots, args) 

                if ss != np.inf and ss < args.voters:
                    ag_matrix[c][o] = ss
                    possibilities[o].append((ss,(c,"AG", o, ss)))

                    print("{} AG {} = {}".format(candidates[c].id, \
                        candidates[o].id, ss), file=log)

    rule_out_assertions = {}

    # Now, we are looking for candidates c where there are at least
    # nseat other candidates o such that o AG c. This is the case where
    # possibilities[c] contains at least nseat sample sizes.
    for c in cands:
        if len(possibilities[c]) > nseats-1:
            sorted_poss = sorted(possibilities[c])[:nseats]
                
            sample = max([ss for ss,_ in sorted_poss])
            rule_out_assertions[c] = set([d for _,d in sorted_poss])

            # We can rule candidate 'c' out with sample size 'sample'
            rule_out.append((c, sample))
            ruled_out.append(c)
            max_rule_out_ss = max(max_rule_out_ss, sample)

    # 'max_rule_out_ss' gives us the sample size required to rule
    # out the candidates in 'rule_out'.
    print("Ruling out candidates: {}".format([(candidates[c].id, ss) for\
        c,ss in rule_out]), file=log)

    print("STAGE 1 ASN (Ruling out candidates): {}".format(\
        max_rule_out_ss), file=log)

    asn_overall = 0
    cset_stage_asn = 0
    assertions_used = set()

    posswinners = [c for c in cands if not c in ruled_out]

    csets = [cs for cs in list(itertools.combinations(posswinners, nseats)) \
        if not(set(cs) == winnerset)]

    cannot_rule_out = []

    with Pool(args.cpus) as pool:
        cset_proc = partial(rule_out_n, candidates, cands, ruled_out,\
            args, INVALID, ag_matrix, ballots, maxtv)

        results = pool.map(cset_proc, csets)

        for failure, cset, best_asn, best_assertions, info in results:
            print(info, file=log)
            if failure:
                asn_overall = np.inf
                cannot_rule_out.append(cset)
            else:
                cset_stage_asn = max(cset_stage_asn, best_asn)
                asn_overall = max(asn_overall, best_asn)
                assertions_used.update(best_assertions)

        print("RULE OUT CSET STAGE ASN {}".format(cset_stage_asn), file=log)
                 
    # Can we reduce sample size of audit?
    print("------------------------------------------------",file=log)
    print("EXAMINING WAYS TO REDUCE COST OF AUDIT", file=log)
    rule_out.sort(key=lambda x: x[1], reverse=True)

    while max_rule_out_ss > cset_stage_asn and rule_out != []:
        r, ss = rule_out[0]

        # New csets that can be made with r as a member.
        not_ruled_out = [c for c in cands if not c in ruled_out]

        new_csets = [(r,)+ccset for ccset in \
            list(itertools.combinations(not_ruled_out, nseats-1))]

        new_cset_asn = 0
        new_assertions = set()

        break_out = False

        with Pool(args.cpus) as pool:
            cset_proc = partial(rule_out_n, candidates, cands, \
                ruled_out, args, INVALID, ag_matrix, ballots, maxtv)
            results = pool.map(cset_proc, new_csets)

            for failure, cset, best_asn, assertions, info in results:
                print(info, file=log)
                if failure:
                    break_out = True
                    break
                else:
                    new_cset_asn = max(new_cset_asn, best_asn)
                    new_assertions.update(assertions)

        if break_out:
            break

        if new_cset_asn < ss:
            print("{} ruled out by csets, sample size {}".format(\
                r, new_cset_asn), file=log)
            max_rule_out_ss = 0 if ruled_out == [r] else \
                max([ros for _,ros in rule_out[1:]])
            rule_out = [] if rule_out == [r] else rule_out[1:]
            ruled_out.remove(r)
            del rule_out_assertions[r]

            print("STAGE 1 ASN NOW {}".format(max_rule_out_ss), file=log)
            
            assertions_used.update(new_assertions)
            cset_stage_asn = max(cset_stage_asn, new_cset_asn)

            print("NEW RULE OUT PAIR STAGE ASN {}".format(\
                cset_stage_asn),file=log)
        else:
            print("BREAK OUT OF STAGE 4", file=log)
            break

    for _,alist in rule_out_assertions.items():
        assertions_used.update(alist)

    asn_overall = max(asn_overall, max(max_rule_out_ss, cset_stage_asn))
 
    print("------------------------------------------------",file=log)
    for cset in cannot_rule_out:
        cset_cands = [candidates[c].id for c in cset]
        print("Alternate outcome cannot be ruled out: {}".format(\
            cset_cands), file=log)

    partial_asn = asn_overall
    if asn_overall >= args.voters:
        asn_overall = np.inf
        print("Full audit not possible", file=log)

        partial_asn = max(max_rule_out_ss,cset_stage_asn)

        print("PARTIAL audit sample size: {}".format(partial_asn), file=log)
    else:        
        print("FULL audit sample size: {}".format(asn_overall), file=log)
        

    print("------------------------------------------------",file=log)
    print("Final set of assertions generated:", file=log)
    for asstn in set(assertions_used):
        print(asstn, file=log)

    if cannot_rule_out == []:
        print("All winners are verified in audit.", file=log)
    else:
        known_winners = set()

        for w in winners:
            in_all = True
            for cset in cannot_rule_out:
                if not w in cset:
                    in_all = False
                    break
            if in_all:
                known_winners.add(w)

        if known_winners != []:
            print("{} are verified winners with these assertions.".format(\
                [candidates[w].id for w in known_winners]), file=log)

        candidate_winners = set(winners)
        for cset in cannot_rule_out:
            candidate_winners.update(set(cset))
            
        candidate_winners = set([w for w in candidate_winners if \
            not w in known_winners])

        if candidate_winners != []:
            print("{} remain as potential winners in this audit.".format(\
                [candidates[w].id for w in candidate_winners]), file=log)

        def_losers = [c for c in cands if not(c in candidate_winners or \
            c in known_winners)]
        print("Definite losers: {}".format([candidates[c].id for c in \
            def_losers]), file=log)
    
    print("------------------------------------------------",file=log)
 

    print("NGEN,{},{},{},{},{},{}".format(args.data, ncand, \
            valid_ballots, args.quota, asn_overall, partial_asn))

    log.close() 
