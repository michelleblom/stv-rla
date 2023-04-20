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


from utils import read_outcome, ssm_sample_size,\
    c_vs_c_sample_size_simple, read_ballots_txt, index_of, next_cand, \
    read_ballots_json, simulate_stv, c_vs_c_sample_size_amargin

# Make sure shangrla is in your PYTHONPATH
from shangrla.NonnegMean import NonnegMean


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
def vote_for_cand_ags1(c, prefs, c_ag):
    ags_used = 0

    for p in prefs:
        if p == c:
            return True,ags_used

        if p != c:
            if p in c_ag:
                ags_used = max(ags_used, c_ag[p])
            else:
                return False, None

    return False,None


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
def vote_for_cand_ags2(c1, c2, prefs, ag_matrix, winners):
    ag_present = False
    ag_min_ss = np.inf

    for p in prefs:
        if p == c1:
            return True, ag_present, ag_min_ss

        if p == c2:
            return False, False, None

        if p in winners:
            # c1 could still get this vote as it could skip over 'p'
            continue

        ag = ag_matrix[p][c1]
        if ag != None and ag != np.inf:
            ag_present = True
            ag_min_ss = min(ag_min_ss, ag)

    return False, False, None



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
def rule_out_for_max(prefs, loser, ag_matrix, winners):
    ag_min_ss = np.inf
    ag_present = False

    for p in prefs:
        if p == loser:
            return ag_present, ag_min_ss

        if p in winners:
            # 'loser' could still get this vote as it could skip over 'p'
            continue

        ag = ag_matrix[p][loser]
        if ag != None and ag != np.inf:
            ag_present = True
            ag_min_ss = min(ag_min_ss, ag)

    return False, np.inf


# Compute max possible tally of z upon its elimination and the portion of
# that tally that could move to a given candidate w.
#
# Compute number of ballots on which 'z' is preferenced before any of
# the candidates in 'standing', or where 'z' appears and standing == [].
#
# Exclude ballots where a 'c' is preferenced before 'z' and 'c' AG 'z' if
# the cost of that AG is less than current asn of audit. 
def max_tally_z_exc(z, w, ballots, ag_matrix, standing, asn_overall):
    count_max = 0
    count_z_to_w = 0
    for blt in ballots:
        for p in blt.prefs:
            if p in standing:
                break

            if p == z:
                count_max += blt.votes
                if w in blt.prefs:
                    count_z_to_w += blt.votes
                break

            ag = ag_matrix[p][z] 
            if ag != None: # and ag <= asn_overall:
                break

    return count_max, count_z_to_w


# Compute an upper bound on the possible transfer value for a candidate c1
# assuming they are elected with a quota. Return upper bound and ASN required
# to show that the transfer value for c1, if they are at any point elected
# with a quota, is less than that upper bound. Idea: compute maximum surplus
# for c1 given that someone has to be eliminated to elect them. Look at
# the maximum tally that someone could have when [c1, c2, o] are still
# standing, and the amount of that tally that could go to c1. 
def max_tv_by_z(c1, c2, o, cands, candidates, ballots, delta_ut, args, log,
    valid_ballots, INVALID, asn_overall, ag_matrix, test):

    # NOTE: C1 is the assumed winner for which we are computing TV.

    zeds_maxtv = []
    max_maxtv_v = -1

    # 'z' is the candidate who will give 'c1' its quota 
    for z in cands:
        if z == c1 or z == c2 or z == o: continue

        if ag_matrix[z][c1] != None and ag_matrix[z][c1] <= asn_overall:
            continue
        if ag_matrix[z][c2] != None and ag_matrix[z][c2] <= asn_overall:
            continue
        if ag_matrix[z][o] != None and ag_matrix[z][o] <= asn_overall:
            continue

        # Compute maximum tally of z when [c1, c2, o] are still standing.
        # Also return the portion of those votes that could move to c1. 
        max_t_z, max_t_z_c = max_tally_z_exc(z, c1, ballots, \
            ag_matrix, [c1, c2, o], asn_overall)

        if max_t_z_c > args.quota:
            max_t_z = args.quota
            max_t_z_c = args.quota


        # If cand z gives c1 its quota, maximum transfer 
        # value will be: max_t_z_c / (max_t_z_c + args.quota)
        max_tv_poss = max_t_z_c / (max_t_z_c + args.quota)

        zeds_maxtv.append((z, max_t_z, max_t_z_c))
        max_maxtv_v = max(max_maxtv_v, max_tv_poss)

    if max_maxtv_v >= args.maxtv or max_maxtv_v == -1:
        return args.maxtv, 0

    cmax = max_maxtv_v + delta_ut
    final_ss = -1
                
    while cmax < args.maxtv:
        # Estimate difficulty of ensuring that for all z in
        # zeds_maxtv, their maximum "to c1" tally Vz,c <= cmax * Q / (1-cmax) 

        # Upper boud on max z-c tally for the current cmax
        Vzc = (cmax * args.quota) / (1 - cmax)

        prop_other = 1 - (Vzc / valid_ballots) 

        # Continue to increase cmax until we can form all 
        # |zeds_maxtv| with less than certain cost
        max_ss = 0
        for z,mtz,mtzc in zeds_maxtv:
            # Want a mtzc < Vzc assertion

            # Imagine all ballots that we do not categorise
            # as belonging in the maximum tally of z->c belong
            # to an imaginary candidate "other". We want
            # to show that "other" has more than prop_other
            # proportion of the ballots.
            tally_other = valid_ballots - mtzc

            ss, _, _ = ssm_sample_size(prop_other, \
                tally_other, INVALID, args, test)

            #print("With z = {}, ASN = {}".format(candidates[z].id,\
            #    ss), file=log)
            max_ss = max(max_ss, ss)

        if max_ss <= asn_overall:
            final_ss = max_ss
            break

        cmax += delta_ut

    if cmax >= args.maxtv or final_ss == np.inf:
        return args.maxtv, 0

    print("    REV MAX TV for C1,C2 = {}, {}, o = {}, is {}, asn {}".format(\
        candidates[c1].id, candidates[c2].id, candidates[o].id, \
        cmax, final_ss), file=log)

    return cmax, final_ss



        

# Merge a sequence of AG relationships that could be used to increase the
# assorter margin of an NL.
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

        lhelpfuls = len(helpful_ags)

        to_add = False
        while cntr < lhelpfuls:
            ag_asn, extra = helpful_ags[cntr]
            if ag_asn != curr_ag_asn:
                merged_helpful_ags.append((curr_ag_asn, curr_extra))
                merged_total += curr_extra

                curr_ag_asn = ag_asn
                curr_extra = extra
            else:
                curr_extra += extra

            cntr += 1                    

        merged_helpful_ags.append((curr_ag_asn, curr_extra))
        merged_total += curr_extra

    assert(merged_total >= exp_merged_total-0.00001 \
        and merged_total <= exp_merged_total+0.00001)

    return merged_helpful_ags


def reduce_margin(amargin, assorter_total, helpful_ags, voters, w, l, \
    asn_overall):

    # Incorporate use of all AGs that either make the assertion
    # possible, or whose ASN is already within/equal to current
    # lower bound on audit difficulty.
    curr_amargin = amargin
    curr_assorter_total = assorter_total

    max_ags_used = 0

    add_min, add_min_papers, del_max, del_max_papers = 0, 0, 0, 0
    while amargin <= 0 and helpful_ags != []:                
        ag_asn, extra, c, val, papers = helpful_ags.pop(0)

        curr_assorter_total += extra

        if c == w:
            add_min += val
            add_min_papers += papers
        else:
            del_max += val
            del_max_papers += papers

        curr_amargin = 2*(curr_assorter_total/voters) - 1

        max_ags_used = max(max_ags_used, ag_asn)

    for ag_asn, extra, c, val, papers in helpful_ags:
        if ag_asn > max(max_ags_used, asn_overall):
            break

        curr_assorter_total += extra

        if c == w:
            add_min += val
            add_min_papers += papers
        else:
            del_max += val
            del_max_papers += papers

        curr_amargin = 2*(curr_assorter_total/voters) - 1

        max_ags_used = max(max_ags_used, ag_asn)

    return curr_amargin, max_ags_used, add_min, add_min_papers, \
        del_max, del_max_papers



if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    # Input: stv data file
    parser.add_argument('-d', dest='data')

    # Input: anticipated error rate (default value is 0)
    parser.add_argument('-e1', dest='erate1', default=0.002, type=float)
    parser.add_argument('-e2', dest='erate2', default=0, type=float)

    # Input: risk limit (default is 5%)
    parser.add_argument('-r', dest='rlimit', default=0.10, type=float)
    
    # Input: parameter for some risk functions
    parser.add_argument('-g', dest='g', default=0.1, type=float)

    # Input: parameter for some risk functions
    parser.add_argument('-t', dest='t', default=1/2, type=float)
    
    # Risk function to use for estimating sample size
    parser.add_argument('-rf', dest='rfunc', default="kaplan_kolmogorov")

    # Input: number of repetitions to perform in simulation to determine
    # an initial sample size estimation -- the quantile of the sample
    # size is computed (default is 1 repetition -- no error rate)
    parser.add_argument('-reps', dest='reps', default=100)

    # Input: seed (default is 9368663)
    parser.add_argument('-s', dest='seed', default=9368663, type=int)

    # Input: number of seats in election (default is 2)
    parser.add_argument('-seats', dest='seats', default=2, type=int)

    # Input: election outcome file
    parser.add_argument('-outcome', dest='outcome')
    
    # Input: Quota
    parser.add_argument('-quota', dest='quota', type=int)

    # Input: Total voters (used to express total valid+informal ballots)
    parser.add_argument('-voters', dest='voters', type=int)

    # Input: maximum assumed value of transferred ballots after seating
    parser.add_argument('-maxtv', dest='maxtv', default=2/3, type=float)

    # Output: Log file 
    parser.add_argument('-log', dest='log', type=str)
    
    # Flags
    # Just simulate STV election
    parser.add_argument('-justsim', dest='justsim', default=False,\
        action='store_true')
    
    # Run two quota method, instead of 1 quota method, if the instance is 
    # suitable
    parser.add_argument('-twoq',dest='twoq',default=False,action='store_true')

    # Use old 1 quota method when 1 quota method is applicable.
    parser.add_argument('-old1q',dest='old1q',default=False,action='store_true')

    # Run general method for instance irrespective of outcome structure
    parser.add_argument('-gen', dest='gen', default=False,action='store_true')

    args = parser.parse_args()

    log = open(args.log, "w")

    assert(args.seats == 2) # Methods only work for 2 seats currently.

    # Read STV data file
    candidates, ballots, id2group, cid2num, valid_ballots = \
        None, None, None, None, None

    # Check for given input data type
    if args.data.endswith(".json"):
        candidates, ballots, id2group, cid2num, valid_ballots = \
            read_ballots_json(args.data)

    elif args.data.endswith(".txt") or args.data.endswith(".us"):
        candidates, ballots, id2group, cid2num, valid_ballots = \
            read_ballots_txt(args.data)
    else:
        print("Unsupported data file type.")
        exit()

    if args.justsim:
        order_c = []
        order_a = []

        simulate_stv(ballots, candidates, args.seats, order_c, order_a,log=log)   

        print("{}".format(candidates[order_c[0]].id), end='')

        for i in range(1, len(candidates)):
            print(",{}".format(candidates[order_c[i]].id), end='')

        print("")

        print("{}".format(order_a[0]), end='')
        for i in range(1, len(candidates)):
            print(",{}".format(order_a[i]), end='')
        exit()


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

    print("VALID BALLOTS {}".format(valid_ballots), file=log)


    # Create test for estimating sample sizes (use default settings)
    nnm = NonnegMean(test=NonnegMean.kaplan_kolmogorov, \
        estim=NonnegMean.optimal_comparison, N = args.voters, t=0.5, g=0.1)
    np.seterr(all="ignore")

    if args.twoq and outcome.action[0] == 1 and outcome.action[1] == 1:
        # TWO QUOTA METHOD: run if applicable and args.twoq flag set
        max_sample_size = 0

        # Check that first winner's FP is greater than a quota
        first_winner = candidates[winners[0]]
        thresh = 1.0/(args.seats + 1);
    
        ss, _, _ = ssm_sample_size(thresh,first_winner.fp_votes,\
            INVALID, args, nnm)
        print("{} ballot checks required to assert that first winner"\
            " has a quota's worth of votes".format(ss), file=log)

        max_sample_size = max(max_sample_size, ss)

        # Check that second winner's FP is greater than a quota
        second_winner = candidates[winners[1]]
    
        ss, _, _ = ssm_sample_size(thresh,second_winner.fp_votes,\
            INVALID, args, nnm)
        print("{} ballot checks required to assert that second winner"\
            " has a quota's worth of votes".format(ss), file=log)

        max_sample_size = max(max_sample_size, ss)

        print("2Q,{},{},{},{},{}".format(args.data, ncand, valid_ballots, \
            args.quota, max_sample_size))


    elif (not args.gen) and outcome.action[0] == 1:
        # NEW 1 QUOTA METHOD
        # CASE: 2 seats and first winner has more than a quota on
        # first preferences.

        max_sample_size = 0

        # Check that first winner's FP is greater than a quota
        first_winner = candidates[winners[0]]
        thresh = 1.0/(args.seats + 1);
    
        ss, _, _ = ssm_sample_size(thresh,first_winner.fp_votes,\
            INVALID, args, nnm)
        print("{} ballot checks required to assert that first winner"\
            " has a quota's worth of votes".format(ss), file=log)

        max_sample_size = max(max_sample_size, ss)

        # For candidate upper bounds on the transfer value of the ballots
        # leaving the first winner's tally pile -- denoted aud_tv -- find
        # the ASN of assertions that show that the second winner cannot be
        # eliminated ('never loses' NL's) before/to any reported loser. 
        
        # We start with aud_tv equalling the reported transfer value plus
        # a small delta (0.01).

        # For each choice of upper bound, we use it in the construction of
        # NL's between the second winner and reported loser. Having a lower
        # upper bound on the transfer value makes it easier to create these
        # NL's. 

        # We consider the overall ASN of an audit configuration for a given
        # choice of aud_tv. This includes:
        # -- ASN for determining that the first winner has quota in round 1
        # -- ASN for checking that the first winners' transfer value is 
        #    below aud_tv
        # -- The ASN for 'd' NL 'second winner' for each reported loser 'd'

        # We keep increasing aud_tv while the overall ASN of the resulting
        # audit is decreasing. Once it starts increasing, we stop. 

        act_tv = (first_winner.fp_votes - args.quota)/first_winner.fp_votes

        # Compute basic AG relationships between candidates, without
        # any assumptions around who is seated.
        ag_matrix = [[None for c in cands] for o in cands]

        for c in cands:
            fpc = candidates[c].fp_votes

            for o in cands:
                if c == o: continue

                max_vote_o = 0

                for b in ballots:
                    # Function will return which of 'c' and 'o' should
                    # be awarded this ballot (None if neither).
                    awarded = vote_for_cand_over_cand(c, o, b.prefs)

                    if awarded == o:
                        max_vote_o += b.votes

                if fpc > max_vote_o:
                    ss, m, _ = c_vs_c_sample_size_simple(fpc,max_vote_o,\
                        valid_ballots, args, nnm) 

                    if ss != np.inf:
                        ag_matrix[c][o] = ss
                        print("AG({},{}) = {}".format(candidates[c].id, \
                            candidates[o].id, ss), file=log), 

        min_tv = 0

        max_in_outer_loop = np.inf

        while min_tv < act_tv - 0.01:
            mintv_ss = 0
           
            # First imagine that min tv for first winner is 0, so we are       
            # not generating any LT assertions. 
            if min_tv > 0:
                mintally = args.quota / (1 - min_tv)
                thresh = mintally/valid_ballots

                mintv_ss, _, _ = ssm_sample_size(thresh,first_winner.fp_votes,\
                    INVALID, args, nnm)
                print("Sample size to show min tv of {} is {}".format(min_tv,\
                    mintv_ss), file=log)

                if mintv_ss >= max_in_outer_loop:
                    break

            aud_tv = act_tv + 0.01

            max_in_loop = np.inf

            # MAIN LOOP OF ONE QUOTA METHOD
            while aud_tv < args.maxtv:

                # Check that TV of first winner is at most aud_TV
                a = 1 / (1 - aud_tv)
                thresholdT = a * valid_ballots / (args.seats + 1)
                thresholdProp = thresholdT / valid_ballots
                threshold = 1 - thresholdProp
        
                tally_others = valid_ballots - first_winner.fp_votes

                ss, _, _ = ssm_sample_size(threshold, tally_others, \
                    INVALID, args, nnm)

                max_this_loop = max(ss, max(mintv_ss, max_sample_size))

                # For second winner, is it the case that they cannot be 
                # eliminated prior to all remaining candidates?
                sw = candidates[winners[1]]

                ags = {}

                max_with_ags = ss
                max_with_nls = ss

                print("AUD TV {}, ss {}".format(aud_tv, ss), file=log)
                # Compute AGs between original losers and second winner
                # Note: we do this for each different upper bound on the 
                # transfer value as we can form more AGs this way -- utilising
                # the fact that ballots leaving the first winner will be reduced
                # in value to aud_tv.
                for c in cands:
                    if c in winners:
                        continue

                    cand_c = candidates[c]

                    min_sw = 0
                    max_c = 0

                    # assertion: fpc(sw) > maxc
                    assorter = INVALID*0.5 # h(b) = ((b_sw - b_c) + 1)/2
                    for b in ballots:
                        if b.prefs[0] == sw.num:
                            # assorter value is 1 per instance of ballot type
                            assorter += b.votes
                            min_sw += b.votes
                            continue

                        if b.prefs[0] == c: 
                            # assorter value is 0 per instance of ballot type
                            max_c += b.votes
                            continue

                        # Default contribution of each instance of ballot type
                        # to assorter.
                        contrib = 0.5 * b.votes 

                        if min_tv > 0 and len(b.prefs) > 1 and b.prefs[:2] == \
                            [first_winner.num, sw.num]:
                            min_sw += min_tv * b.votes
                            contrib = b.votes * ((1 + min_tv)/2.0)
                        else:
                            weight = 1

                            if b.prefs[0] == first_winner.num:
                                weight = aud_tv

                            for p in b.prefs:
                                if p == sw.num:
                                    break

                                if p == c:
                                    contrib = b.votes * ((1 - weight)/2.0)
                                    max_c += weight*b.votes
                                    break

                        assorter += contrib

                    # Compute assorter margin
                    amargin = 2 * (assorter/args.voters) - 1

                    if amargin > 0:
                        ss, _, _ = c_vs_c_sample_size_amargin(amargin,args,nnm) 

                        ags[c] = ss

                        max_with_ags = max(max_with_ags, ss)
                        print("AG({},{}) = {}".format(sw.id, cand_c.id, ss),\
                            file=log)
                    else:
                        max_with_ags = np.inf
                        print("AG({},{}) NOT POSSIBLE".format(sw.id, \
                            cand_c.id),file=log)
                        print("min {}, max {}, amargin {}".format(min_sw,\
                            max_c, amargin), file=log)


                # Determine NL's between original losers and second winner
                for c in cands:
                    if c in winners:
                        continue

                    cand_c = candidates[c]

                    min_sw = 0 # Min tally for second winner.
                    max_c = 0 # Max tally for candidate c.

                    # Keep running tally of total votes we can increase the
                    # margin of assertion 'sw' NL 'c' with 'AG' relationships
                    pot_margin_inc = 0

                    helpful_ags = []

                    # Assertion: min_sw > maxc
                    assorter = INVALID*0.5
                    for b in ballots:
                        if b.prefs[0] == sw.num:
                            assorter += b.votes
                            min_sw += b.votes
                            continue

                        if b.prefs[0] == c:
                            max_c += b.votes
                            continue

                        # is this a ballot for 'c' over 'sw'?
                        c_in = c in b.prefs
                        s_in = sw.num in b.prefs

                        weight = 1
                        if b.prefs[0] == first_winner.num:
                            weight = aud_tv

                        c_idx = b.prefs.index(c) if c_in else np.inf
                        s_idx = b.prefs.index(sw.num) if s_in else np.inf

                        if c_in and c_idx < s_idx:
                            # These ballots could not go to 'sw' but may go
                            # to 'c'.

                            # Could we exclude these ballots by using an AG?
                            is_ag, ag_asn = rule_out_for_max(b.prefs, c, \
                                ag_matrix, winners)
                            
                            contrib = b.votes*((1-weight)/2.0)

                            if is_ag:
                                alt_contrib = b.votes*0.5
                                helpful_ags.append((ag_asn,alt_contrib-contrib))
                                pot_margin_inc += alt_contrib-contrib

                            assorter += contrib
                            max_c += b.votes

                        elif s_in:
                            if min_tv > 0 and len(b.prefs) > 1 and \
                                b.prefs[:2] == [first_winner.num, sw.num]:
                                min_sw += min_tv * b.votes
                                assorter += b.votes * ((1 + min_tv)/2.0)
                            else:
                                # If we remove all cand 'd' for which sw AG d
                                # that appear before sw in the ballot, will 'sw'
                                # be the first ranked cand? If so, we could add
                                # these ballots to the min tally of 'sw'.
                                prefs = b.prefs[:]
           
                                max_ags_here = 0
                                for d,dval in ags.items():
                                    if d in prefs:
                                        idx_d = prefs.index(d)
                                        if idx_d < s_idx:
                                            prefs.remove(d)
                                            s_idx -= 1
                                            max_ags_here=max(max_ags_here,dval)

                                assorter += 0.5*b.votes

                                if prefs[0] == sw.num:
                                    # These ballots would have originally had a 
                                    # contribution of 0.5 to the assorter. By
                                    # utilising these AGs we can increase the 
                                    # contribution of these ballots to 1 each.
                                    helpful_ags.append((max_ags_here,  \
                                        0.5*b.votes))

                                    pot_margin_inc += b.votes*0.5
                                elif min_tv > 0 and  len(prefs) > 1 and \
                                    prefs[:2] ==  [first_winner.num, sw.num]:
                                    base_contrib = 0.5*b.votes
                                    alt_contrib = b.votes * ((1 + min_tv)/2.0)
                                    dconfig = alt_contrib - base_contrib

                                    if dconfig > 0:
                                        helpful_ags.append((max_ags_here, \
                                            dconfig))

                                        pot_margin_inc += dconfig

                        else:
                            assorter += 0.5*b.votes

                    # Max ASN of any AGs used to increase assorter margins 
                    # when forming NLs.
                    max_ags_used = 0  
                    merged_helpful_ags = merge_helpful_ags(helpful_ags, \
                        pot_margin_inc)
      
                    # Incorporate use of all AGs that either make the assertion
                    # possible, or whose ASN is already within/equal to current
                    # lower bound on audit difficulty.
                    while assorter/args.voters <= 0.5 and \
                        merged_helpful_ags != []:
                        
                        ag_asn, extra_contrib = merged_helpful_ags.pop(0)

                        assorter += extra_contrib
                        max_ags_used = max(max_ags_used, ag_asn)

                    cntr = 0
                    for ag_asn, extra_contrib in merged_helpful_ags:
                        if ag_asn > max(max_ags_used, max_with_nls):
                            break

                        assorter += extra_contrib
                        cntr += 1

                    if cntr > 0 and merged_helpful_ags != []:
                        merged_helpful_ags = merged_helpful_ags[cntr:]
            
                    # Can we reduce the sample size required for the assertion
                    # by including more AGs?
                    amargin = 2*(assorter/args.voters) - 1
                    if amargin > 0:
                        # Current sample size for assertion
                        ss, _, _ = c_vs_c_sample_size_amargin(amargin,args,nnm) 

                        # Go through remaining AGs 
                        for ag_asn, extra_contrib in merged_helpful_ags:
                            if ss < ag_asn:
                                break

                            # would reducing margin by votes reduce sample size
                            # by more than increase caused by including AG?
                            assorter += extra_contrib
                            amargin = 2*(assorter/args.voters)-1

                            ss,_, _=c_vs_c_sample_size_amargin(amargin,args,nnm)
                            max_ags_used = max(max_ags_used, ag_asn)
                                
                        max_with_nls = max(max_with_nls, max(max_ags_used,ss))
                        print("NL({},{}) = {}, AGs used {}".format(sw.id, \
                            cand_c.id, max(max_ags_used,ss), max_ags_used),\
                            file=log)
                    else:
                        max_with_nls = np.inf
                        print("NL({},{}) NOT POSSIBLE".format(\
                            sw.id, cand_c.id),file=log)


                max_this_loop=max(max_this_loop,min(max_with_ags,max_with_nls))

                print("Max in loop {}, {}".format(max_this_loop, aud_tv), \
                    file=log)

                if max_this_loop <= max_in_loop:
                    max_in_loop = max_this_loop
                    aud_tv += 0.01
                else:
                    break
            
            if max_in_loop <= max_in_outer_loop:
                max_in_outer_loop = max_in_loop
            else:
                break    

            if min_tv == 0:
                min_tv = act_tv/2.0
            else:
                min_tv += 0.01

        max_sample_size = max(max_sample_size, max_in_outer_loop)

        if max_sample_size >= args.voters:
            max_sample_size = np.inf 

        print("Sample size required for audit is {} ballots".format(\
            max_sample_size), file=log)

        print("1Q,{},{},{},{},{}".format(args.data, ncand, valid_ballots, \
            args.quota, max_sample_size))

    else: 
        # CASE: 2 seats and no candidate has a quota on first preferences
        # according to reported results.  
        overall_asn = 0

        remcase = outcome.action[ncand-2:] == [1,1]

        for w in winners:
            print("Winner {}".format(candidates[w].id), file=log)

        # First look at whether we can batch eliminate some candidates at
        # the start of counting. If we can show that E are eliminated first,
        # before any candidate can get a quota, then for all other candidates
        # 'c' (not in E) we can give 'c' votes from any candidate in E in
        # their minimum tally as part of any NL calculation. 
        batch_asn = np.inf
        first_elim = [None for c in candidates]

        if ncand > 4 and outcome.action[:3] == [0,0,0]:
            idx = -1

            # Can we batch eliminate the first few candidates?
            for i in range(2,ncand):
                if outcome.action[i] == 1:
                    break

                ilosers = outcome.cand[:i]
                iwinners = outcome.cand[i:]
       
                # compute total fp's of the losers
                tally_losers = 0
                for c in ilosers:
                    tally_losers += candidates[c].fp_votes

                # minimum fp of the winners
                tally_winner = candidates[iwinners[0]].fp_votes

                # Can we audit the winner vs losers comparison?
                ss, m, _ = c_vs_c_sample_size_simple(tally_winner, \
                    tally_losers, valid_ballots, args, nnm) 

                if ss < batch_asn and ss < args.voters/10:
                    batch_asn = ss
                    idx = i

            if idx != -1:
                # we can perform a batch elimination
                ilosers = outcome.cand[:idx]
                iwinners = outcome.cand[idx:]

                loser_ids = [candidates[c].id for c in ilosers]

                print("Batch elimination of {} with ASN {}".format(loser_ids,\
                    batch_asn), file=log)

                # check that no candidate has a quota upon the elimination
                # of ilosers
                tallies_on_be = [0 for c in candidates]
        
                for blt in ballots:
                    # Who is the first ranked candidate excluding ilosers
                    next_c = next_cand(blt.prefs, ilosers)

                    if next_c != None:
                        tallies_on_be[next_c] += blt.votes

                for c in iwinners:
                    tally = tallies_on_be[c]

                    # Can we show that 'c' does not have a quota? We need 
                    # to show that the super candidate C \ {c} has a prop
                    # of the votes that is more than 1 - 1/(seats + 1)
                    thresh = 1 -(1.0/(args.seats + 1));

                    ss, _, _ = ssm_sample_size(thresh, valid_ballots-tally, \
                        INVALID, args, nnm)

                    print("Show that {} does not have a quota after "\
                        "batch elim of {}, ASN = {}".format(candidates[c].id,\
                        loser_ids, ss), file=log)

                    batch_asn = max(batch_asn, ss)
            
                    if batch_asn is np.inf:
                        break

                
                # If we can show that no one has a quota after the elimination
                # of ilosers, then we know that there must be at least one
                # more IRV elimination! Create an assertion to show that 
                # the next eliminated candidate in the reported outcome is 
                # the correct one.
                if batch_asn != np.inf:
                    for c in ilosers:
                        first_elim[c] = batch_asn

                    next_c = outcome.cand[idx]

                    asn_next_elim = 0
                
                    # Assert that next_c has less votes at this point than
                    # all other remaining candidates
                    tally_next_c = tallies_on_be[next_c]
                
                    for i in range(idx+1, ncand):
                        next_w = outcome.cand[i]

                        ss,_,_ =  c_vs_c_sample_size_simple(\
                            tallies_on_be[next_w],tally_next_c, \
                            valid_ballots, args, nnm) 

                        print("{} vs {}, ASN {}".format(candidates[next_w].id,\
                            candidates[next_c].id, ss), file=log)

                        asn_next_elim = max(asn_next_elim, ss)
                        if ss is np.inf:
                            break

                    if asn_next_elim < np.inf:
                        comb_asn = max(batch_asn, asn_next_elim)
                        first_elim[next_c] = comb_asn

                        print("Candidate {}".format(candidates[next_c].id)\
                            + " is the next eliminated, ASN {}".format(\
                            comb_asn), file=log)

        # ag_matrix[c][o] will give you the sample size required to show
        # that c AG o. It will be None if no such AG relationship exists.
        ag_matrix = [[None for c in cands] for o in cands]

        rule_out = []
        max_rule_out_ss = 0
        ruled_out = []

        delta_ut = 0.01

        # possibilities[c] will give you all the sample sizes of 
        # AG relationships where 'c' is the loser.
        possibilities = [[] for c in cands]

        print("Computing AG relationships", file=log)
        # Compute "Always greater" AG relationships.
        for c in cands:
            # Cand 'c' is always greater than 'o' if c's first preference
            # tally is greater than the maximum tally for 'o'
            fpc = candidates[c].fp_votes

            for o in cands:
                if c == o: continue

                # The maximum tally for 'o' is equal to the number of
                # ballots that preference 'o' before 'c' or preference
                # 'o' and not 'c'. 
                max_vote_o = 0

                for b in ballots:
                    awarded = vote_for_cand_over_cand(c, o, b.prefs)

                    if awarded == o:
                        max_vote_o += b.votes

                if fpc > max_vote_o:
                    ss, _, _ = c_vs_c_sample_size_simple(fpc,max_vote_o,\
                        valid_ballots, args, nnm) 

                    if ss != np.inf and ss < args.voters:
                        ag_matrix[c][o] = ss
                        possibilities[o].append(ss)

                        print("{} AG {} = {}".format(candidates[c].id, \
                            candidates[o].id, ss), file=log)

        # Supplement AG matrix with information obtained from the batch
        # elimination matrix. 
        for i in range(ncand):
            iv = first_elim[i]

            if iv is None: continue

            for j in range(ncand):
                if j == i: continue

                jv = first_elim[j]

                if jv != None: continue

                # Here 'i' is known to be first eliminated, and 'j' not.
                agv = ag_matrix[j][i]

                if agv is None or iv < agv:
                    ag_matrix[j][i] = iv
                    possibilities[i].append(iv)

        # Can we establish a winner based on the fact that they are AG
        # at least NCAND - 1 candidates? If there is one candidate that 
        # we cannot establish an AG relationship with that is OK. Then
        # we can rule out all alternate outcomes that do not involve that
        # known winner.
        known_winners = []

        for w in winners:
            # for how many candidates c can we show that w AG c
            cost = 0
            num_ag = 0
            for c in cands:
                if c == w:
                    continue

                if ag_matrix[w][c] != None:
                    num_ag += 1
                    cost = max(cost, ag_matrix[w][c])

            if num_ag >= ncand - 2:
                print("Candidate {} is a known winner, cost {}".format(\
                    candidates[w].id, cost), file = log)
                known_winners.append((w,cost))

        if len(known_winners) == args.seats:
            print("All winners can be established with AG relationships.",\
                file=log)
            print("ASN: {}".format(max([cost for _,cost in known_winners])),\
                file=log)

            
            log.close()
            exit(0)   

        # Now, we are looking for candidates c where there are at least
        # 2 other candidates o such that o AG c. This is the case where
        # possibilities[c] contains at least two sample sizes.
        for c in cands:
            if len(possibilities[c]) > 1:
                sorted_poss = sorted(possibilities[c])
                
                sample = max(sorted_poss[0],sorted_poss[1])

                # We can rule candidate 'c' out with sample size 'sample'
                rule_out.append((candidates[c].id, sample))
                ruled_out.append(c)
                max_rule_out_ss = max(max_rule_out_ss, sample)

        # 'max_rule_out_ss' gives us the sample size required to rule
        # out the candidates in 'rule_out'.
        print("Ruling out candidates: {}".format(rule_out), file=log)
        print("Sample size for ruling out cands: {}".format(max_rule_out_ss),\
            file=log)

        # Initially, our ASN consists of the ballot samples required to
        # establish the set of alternative winner pairs, ruling out candidates
        # c where there are at least 2 candidates o such that o AG c.
        asn_overall = max_rule_out_ss

        # Form pairs of alternate winner outcomes
        pairs = []
        assertions = []

        must_win = None
        if known_winners != []:
            must_win = known_winners[0][0]
            asn_overall = max(asn_overall, known_winners[0][1]) 

        for i in range(ncand):
            c = cands[i]
            if c in ruled_out: continue

            for j in range(i+1, ncand):
                o = cands[j]

                if o in ruled_out: continue
                if c in winners and o in winners: continue

                if must_win != None and c != must_win and o != must_win:
                    continue

                pairs.append((c,o))

        assertions = []
        for c1,c2 in pairs:
            print("ALTERNATE OUTCOME PAIR: {},{}".format(candidates[c1].id,\
                candidates[c2].id), file=log)

            # ==============================================================
            # No assumption around order in which candidates
            # are elected required.
            # -------------------------------------------------------------
            # If we assume 'c1' gets a seat at some stage, does there
            # exist a candidate 'o' \= 'c1' that cannot be eliminated before 
            # 'c2'? If so, 'c2' cannot get a seat.

            # We can use as a maximum tally for c1, for the purposes
            # of the o NL c2 test, the number of ballots on which
            # c1 is mentioned before c2 or c1 appears and c2 does not.

            best_asn1 = np.inf # assuming 'c1' gets a seat
            best_asn2 = np.inf # assuming 'c2' gets a seat

            not_applicable = True

            for o in cands:
                if o == c1 or o == c2 or o in ruled_out: continue
           
                ctvmax1, ctvmax_ss1 = max_tv_by_z(c1, \
                    c2, o, cands, candidates, ballots, 0.005, args, log, \
                    valid_ballots, INVALID, asn_overall, ag_matrix, nnm)

                ctvmax2, ctvmax_ss2 = max_tv_by_z(c2, \
                    c1, o, cands, candidates, ballots, 0.005, args, log, \
                    valid_ballots, INVALID, asn_overall, ag_matrix, nnm)

                # Find set of candididates 'c' for which o AG c.
                # Keep track of cost of each of those AGs
                cand_o = candidates[o]
                o_ag = {}
                for c in cands:
                    if c != c1 and c != c2 and ag_matrix[o][c] != None \
                        and ag_matrix[o][c] != np.inf:
                        o_ag[c] = ag_matrix[o][c]

                #=============================================================
                # ASSUMING C1 IS SEATED

                # Max ASN of any AGs used to increase/decrease the minimum
                # /maximum tallies of second winner/candidate c when forming NL.
                max_ags_used1 = ctvmax_ss1 

                # Consider max vote of 'c2' given 'c1' seated at some stage and 
                # 'o' is still standing.
                assorter_o_v_c2 = INVALID*0.5
                assorter_vals = [INVALID*0.5]

                mint, minpapers, maxt, maxpapers = 0, 0, 0, 0
                for b in ballots:
                    awarded, used_ss = vote_for_cand_ags1(o, b.prefs, o_ag)

                    if awarded:
                        if used_ss != None and used_ss > 0:
                            max_ags_used1 = max(max_ags_used1, used_ss)
                            
                        assorter_o_v_c2 += b.votes
                        mint += b.votes
                        minpapers += b.votes
                        
                        assorter_vals.append(b.votes)
                        continue    
                        
                    if not c2 in b.prefs: 
                        assorter_o_v_c2 += b.votes*0.5
                        assorter_vals.append(b.votes*0.5)
                        continue

                    weight = 1.0

                    if b.prefs[0] == c1:
                        weight = ctvmax1

                    # In this analysis, we are assuming that no one other
                    # than c1 and c2 has won -- so if AG(d, c2) and 'd'
                    # is not c1, and 'd' is preferenced before c2, then
                    # c2 does not get these ballots. 
                    awarded, ag_present, used_ss = vote_for_cand_ags2(c2, o, \
                        b.prefs, ag_matrix, [c1])

                    if awarded:
                        if ag_present:
                            assorter_o_v_c2 += b.votes*0.5
                            assorter_vals.append(b.votes*0.5)
                            max_ags_used1 = max(max_ags_used1, used_ss)
                        else:
                            contrib = b.votes*((1 - weight)/2.0)
                            assorter_o_v_c2 += contrib
                            assorter_vals.append(contrib)
                            maxt += weight*b.votes
                            maxpapers += b.votes
                    else:
                        assorter_o_v_c2 += b.votes*0.5
                        assorter_vals.append(b.votes*0.5)

                amargin = 2*(sum(assorter_vals)/args.voters) - 1
                amargin1 = 2*(assorter_o_v_c2/args.voters) - 1

                if  amargin != amargin1:
                    print("{} vs {}".format(amargin, amargin1))

                print("   (1) can we show that {} NL {}? ".format(cand_o.id,\
                    candidates[c2].id), file=log)
                print("        a. margin {}".format(amargin), file=log)

                # Is the minimum tally for 'o' larger than the maximum
                # possible tally for 'c2'? This means 'o' cannot be 
                # eliminated before 'c2'
                if amargin > 0:
                    ss, _, _ = c_vs_c_sample_size_amargin(amargin,args,nnm)

                    if ss != np.inf:
                        best_asn1 = min(max(ss, max_ags_used1), best_asn1)
                        not_applicable = False
                        print("      yes, {}, ags used {} = {}".format(ss, \
                             max_ags_used1, max(ss, max_ags_used1)), file=log)
                    else:
                        print("      no", file=log)
                else:
                    print("      no, min {:.2f},{} vs max {:.2f},{}".format(\
                        100*(mint)/args.voters, minpapers,\
                        100*(maxt)/args.voters, maxpapers), file=log)

                #=============================================================
                # ASSUMING C2 IS SEATED

                # Consider max vote of 'c1' given 'c2' seated at some stage and 
                # 'o' is still standing.
                assorter_o_v_c1 = INVALID*0.5
                mint, minpapers, maxt, maxpapers = 0, 0, 0, 0

                assorter_vals = [INVALID*0.5]

                # Max ASN of any AGs used to increase/decrease the minimum/
                # maximum tallies of second winner/candidate c when forming NL.
                max_ags_used2 = ctvmax_ss2 

                for b in ballots:
                    awarded, used_ss = vote_for_cand_ags1(o, b.prefs, o_ag)

                    if awarded:
                        if used_ss != None and used_ss > 0:
                            max_ags_used2 = max(max_ags_used2, used_ss)
                            
                        assorter_o_v_c1 += b.votes
                        minpapers += b.votes
                        mint += b.votes
                        assorter_vals.append(b.votes)
                        continue    
                        
                    if not c1 in b.prefs: 
                        assorter_o_v_c1 += b.votes*0.5
                        assorter_vals.append(b.votes*0.5)
                        continue

                    weight = 1.0

                    if b.prefs[0] == c2:
                        weight = ctvmax2

                    awarded, ag_present, used_ss = vote_for_cand_ags2(c1, o, \
                        b.prefs, ag_matrix, [c2])

                    if awarded:
                        if ag_present:
                            assorter_o_v_c1 += b.votes*0.5
                            assorter_vals.append(b.votes*0.5)
                            max_ags_used2 = max(max_ags_used2, used_ss)
                        else:
                            contrib = b.votes*((1 - weight)/2.0)
                            assorter_o_v_c1 += contrib
                            assorter_vals.append(contrib)
                            maxt += weight*b.votes
                            maxpapers += b.votes
                    else:
                        assorter_o_v_c1 += b.votes*0.5
                        assorter_vals.append(b.votes*0.5)

                amargin = 2*((sum(assorter_vals)/args.voters)) - 1
 
                amargin1 = 2*(assorter_o_v_c1/args.voters) - 1

                if  amargin != amargin1:
                    print("{} vs {}".format(amargin, amargin1))
 
                print("   (1) can we show that {} NL {}? ".format(cand_o.id,\
                    candidates[c1].id), file=log)
                print("        a. margin {}".format(amargin), file=log)


                # Is the minimum tally for 'o' larger than the maximum
                # possible tally for 'c1'? This means 'o' cannot be 
                # eliminated before 'c1'
                if amargin > 0:
                    ss, _, _ = c_vs_c_sample_size_amargin(amargin,args,nnm)

                    if ss != np.inf:
                        best_asn2 = min(max(ss, max_ags_used2), best_asn2)
                        not_applicable = False
                        print("      yes, {}, ags used {} = {}".format(ss, \
                             max_ags_used2, max(ss, max_ags_used2)), file=log)
                    else:
                        print("      no", file=log)
                else:
                    print("      no, min {:.2f},{} vs max {:.2f},{}".format(\
                        100*(mint)/args.voters, minpapers, \
                        100*(maxt)/args.voters, maxpapers), file=log)

            best_asn = min(best_asn1, best_asn2)
            assertions.append((best_asn, (c1,c2), not_applicable)) 

            asn_overall = max(asn_overall, best_asn)

            if not not_applicable:
                print("Best assertion(s) require sample: {}".format(\
                    best_asn), file=log)
                   
        for _,(c1,c2),na in assertions:
            if na:
                c1o = candidates[c1]
                c2o = candidates[c2]
                print("Alternate outcome cannot be ruled out: {},{}".format(\
                    c1o.id, c2o.id), file=log)

        if asn_overall >= args.voters:
            asn_overall = np.inf
        

        print("Best sample size: {}".format(asn_overall), file=log)
        if remcase:
            print("CASEC,{},{},{},{},{}".format(args.data, ncand, \
                valid_ballots, args.quota, asn_overall))
        else:
            print("CASEB,{},{},{},{},{}".format(args.data, ncand, \
                valid_ballots, args.quota, asn_overall))

    log.close() 
