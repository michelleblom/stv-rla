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
                descs.add((candidates[c].id,"AG",candidates[p].id,c_ag[p]))
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

            descs.add((candidates[p].id, "AG", candidates[c1].id, ag))

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
            ags_used.add((candidates[p].id,"AG",candidates[loser].id, ag))
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
    
    # Input: increment to use when exploring candidate lower and upper 
    # bounds on first winner's transfer value (for 1-quota method).
    parser.add_argument('-deltat', dest='deltat', default=0.05, type=float)
    
    # Flags
    # Just simulate STV election
    parser.add_argument('-justsim', dest='justsim', default=False,\
        action='store_true')
    
    # Run two quota method, instead of 1 quota method, if the instance is 
    # suitable
    parser.add_argument('-twoq',dest='twoq',default=False,action='store_true')

    # Run general method for instance irrespective of outcome structure
    parser.add_argument('-gen', dest='gen', default=False,action='store_true')

    # If present, and running general method -- just check validity of 
    # first winner.
    parser.add_argument('-cfw', dest='cfw', default=False,action='store_true')

    # If present, first N candidates in the recorded outcome are eliminated 
    # as part of a group elimination as it is mathematically impossible for 
    # them to win.
    parser.add_argument('-gelim', dest='gelim', type=int, default=0)
    
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

    if args.justsim:
        order_c = []
        order_a = []
        order_q = {}
        sim_winners = []

        simulate_stv(ballots, candidates, args.seats, order_c, order_a,\
            order_q, sim_winners, log=log)   

        if log != None:
            print("{}".format(candidates[order_c[0]].id), end='', file=log)

            for i in range(1, len(candidates)):
                print(",{}".format(candidates[order_c[i]].id), end='',file=log)

            print("", file=log)

            print("{}".format(order_a[0]), end='', file=log)
            for i in range(1, len(candidates)):
                print(",{}".format(order_a[i]), end='', file=log)
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

    print("VALID BALLOTS {}".format(valid_ballots), file=log)

    # This matrix will be used to store AG relationships between candidates.
    ag_matrix = [[None for c in cands] for o in cands]

    max_sample_size = 0

    # List of text descriptions, and associated sample sizes, of all
    # assertions needed in an audit.
    assertions_used = set()

    # List of candidates that are eliminated in an initial group elimination
    # phase.
    geliminated = []

    # If a group elimination has taken place, generate assertions to verify
    # that the group elimination was valid. 
    if args.gelim > 0:
        geliminated = outcome.cand[:args.gelim]
        
        # A candidate is mathematically not viable for winning according 
        # to the Minneapolis rules if: their total mentions is not enough
        # to surpass the candidate with the next higher current vote total; or
        # in a multiple seat context, their total mentions is less than the
        # current tally of args.seats other candidates.    
        failed_to_verify = []

        options = {g : [] for g in geliminated}

        for c in outcome.cand[args.gelim:]:
            fpc = candidates[c].fp_votes

            for g in geliminated:
                max_vote_g = 0

                for b in ballots:
                    # Function will return which of 'c' and 'g' should
                    # be awarded this ballot (None if neither).
                    awarded = vote_for_cand_over_cand(c, g, b.prefs)

                    if awarded == g:
                        max_vote_g += b.votes

                if fpc > max_vote_g:
                    ss = tally_vs_tally_sample_size(fpc, max_vote_g, \
                        valid_ballots, args) 

                    if ss != np.inf:
                        ag_matrix[c][g] = ss
                        print("GELIM AG({},{}) = {}".format(candidates[c].id, \
                            candidates[g].id, ss), file=log)

                        options[g].append((ss,(candidates[c].id, "AG",\
                            candidates[g].id, ss)))

        for g in geliminated:
            ags = options[g]
            if len(ags) < args.seats:
                failed_to_verify.append(g)
                max_sample_size = np.inf
            else:
                ags.sort()
                for ss,assrtn in ags[:args.seats]:
                    assertions_used.add(assrtn)
                    max_sample_size = max(max_sample_size, ss)
                

        print("Group elimination sample size {}".format(max_sample_size), \
            file=log)

        if failed_to_verify != []:
            print("    Failed to verify elimination of {}".format(\
                [candidates[i].id for i in failed_to_verify]))
             
        outcome.cand = outcome.cand[args.gelim:]
        outcome.action = outcome.action[args.gelim:]

        if len(outcome.cand) == args.seats:
            print("------------------------------------------------",file=log)
            print("Final set of assertions generated:", file=log)
            for asstn in set(assertions_used):
                print(asstn, file=log)
            print("------------------------------------------------",file=log)


            print("GE,{},{},{},{},{}".format(args.data, ncand, valid_ballots, \
                args.quota, max_sample_size))

            exit(0)

        # Rework election profile to remove these eliminated candidates so
        # that the remainder of the analysis is based on next phase of
        # the election. Even if we could not verify that the group 
        # elimination was correct, looking at what other assertions
        # could be formed if we assume it was correct gives insights.
        for cand in candidates:
            cand.reset()

        for blt in ballots:
            newprefs = []

            for p in blt.prefs:
                if p in geliminated:
                    continue
                newprefs.append(p)

            blt.prefs = newprefs
            
            if blt.prefs != []:
                fcand = candidates[blt.prefs[0]]
                fcand.fp_votes += blt.votes
                fcand.ballots.append(blt.num)

                for p in blt.prefs:
                    candidates[p].mentions.append(blt.num)

        # NOTE: remaining logic has to deal with cases where a ballot
        # may have an empty preference list, but still represents a valid
        # set of votes.
        cands = outcome.cand[:]
        ncand = len(cands)
        

    if args.twoq and outcome.action[0] == 1 and outcome.action[1] == 1:
        # TWO QUOTA METHOD: run if applicable and args.twoq flag set

        # Check that first winner's FP is greater than a quota
        first_winner = candidates[winners[0]]
        thresh = 1.0/(args.seats + 1);
    
        ss = ssm_sample_size(thresh,first_winner.fp_votes, INVALID, args)
        print("{} ballot checks required to assert that first winner"\
            " has a quota's worth of votes".format(ss), file=log)

        assertions_used.add((first_winner.id, "QT", ss))

        max_sample_size = max(max_sample_size, ss)

        # Check that second winner's FP is greater than a quota
        second_winner = candidates[winners[1]]
    
        ss = ssm_sample_size(thresh,second_winner.fp_votes,INVALID, args)
        print("{} ballot checks required to assert that second winner"\
            " has a quota's worth of votes".format(ss), file=log)

        max_sample_size = max(max_sample_size, ss)
        
        assertions_used.add((second_winner.id, "QT", ss))

        print("------------------------------------------------",file=log)
        print("Final set of assertions generated:", file=log)
        for asstn in set(assertions_used):
            print(asstn, file=log)
        print("------------------------------------------------",file=log)


        print("2Q,{},{},{},{},{}".format(args.data, ncand, valid_ballots, \
            args.quota, max_sample_size))


    elif (not args.gen) and outcome.action[0] == 1:
        # 1 QUOTA METHOD
        # CASE: 2 seats and first winner has more than a quota on
        # first preferences.

        # Initial check: can we show that both winners are always 
        # greater than all other remaining candidates?

        # NOTE: some of below appears in the general case as well -
        # at some point refactor to move common code into functions.


        # Check that first winner's FP is greater than a quota
        first_winner = candidates[winners[0]]
        second_winner = candidates[winners[1]]
        thresh = 1.0/(args.seats + 1);
    
        ss = ssm_sample_size(thresh, first_winner.fp_votes, INVALID, args)
        print("{} ballot checks required to assert that first winner"\
            " has a quota's worth of votes".format(ss), file=log)

        assertions_used.add((first_winner.id, "QT", None, ss))

        max_sample_size = max(max_sample_size, ss)

        # Increment to use when stepping through candidate upper/lower
        # bounds on first winner transfer value.
        deltat = args.deltat 

        act_tv = (first_winner.fp_votes - args.quota)/first_winner.fp_votes

        min_tv = 0

        max_in_outer_loop = np.inf

        best_outer_assertions = set()
        best_outer_exp = set()

        while min_tv < act_tv - deltat:
            mintv_ss = 0
           
            lts = set()

            # First imagine that min tv for first winner is 0, so we are       
            # not generating any LT assertions. 
            if min_tv > 0:
                mintally = args.quota / (1 - min_tv)
                thresh = mintally/valid_ballots

                mintv_ss = ssm_sample_size(thresh,first_winner.fp_votes,\
                    INVALID, args)
                print("Sample size to show min tv of {} is {}".format(min_tv,\
                    mintv_ss), file=log)

                if mintv_ss >= max_in_outer_loop:
                    break

                lts.add((first_winner.id, "LT", min_tv, mintv_ss))

            aud_tv = act_tv + deltat

            max_in_loop = np.inf

            best_inner_assertions = None
            best_inner_explanations = None

            # MAIN LOOP OF ONE QUOTA METHOD
            while aud_tv < args.maxtv:

                inner_loop_assertions = set()

                # Check that TV of first winner is at most aud_TV
                a = 1 / (1 - aud_tv)
                thresholdT = a * valid_ballots / (args.seats + 1)
                thresholdProp = thresholdT / valid_ballots
                threshold = 1 - thresholdProp
        
                tally_others = valid_ballots - first_winner.fp_votes

                ss = ssm_sample_size(threshold, tally_others, INVALID, args)

                inner_loop_assertions.add((first_winner.id, "MT", \
                    aud_tv, ss))

                max_this_loop = max(ss, max(mintv_ss, max_sample_size))

                # For second winner, is it the case that they cannot be 
                # eliminated prior to all remaining candidates? We will
                # look at this from TWO perspectives. (1) compute AG*'s
                # between the second winner and opponents using our upper
                # and lower bounds on the first winner's transfer value. If 
                # we can form enough AG*'s, we can verify the second winner.
                # (2) Forming NL assertions, potentially making use of our
                #  set of AG* relationships. We will use the 
                # cheaper of the two options.
                sw = candidates[winners[1]]

                ags = {}

                max_with_nls = ss

                print("AUD TV {}, ss {}".format(aud_tv, ss), file=log)
                # Compute AG*'s between original losers and second winner
                # Note: we do this for each different upper/lower bound on the 
                # transfer value as we can form more AG*'s this way -- utilising
                # the fact that ballots leaving the first winner will be reduced
                # in value to aud_tv, and will have a minimum value of min_tv
                for c in cands:
                    if c in winners:
                        continue

                    cand_c = candidates[c]

                    min_sw = 0
                    max_c = 0

                    # assertion: fpc(sw) > maxc
                    assorter = INVALID*0.5 # h(b) = ((b_sw - b_c) + 1)/2
                    for b in ballots:
                        if b.prefs != []:
                            if b.prefs[0] == sw.num:
                                # assorter value is 1 per vote
                                assorter += b.votes
                                min_sw += b.votes
                                continue

                            if b.prefs[0] == c: 
                                # assorter value is 0 per vote
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

                            if b.prefs != [] and b.prefs[0] == first_winner.num:
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
                    amean = assorter/args.voters

                    if amean > 0.5:
                        ss = sample_size(amean, args)

                        ags[c] = ss
                        ag_matrix[sw.num][c] = ss

                        print("AG*({},{}) = {}".format(sw.id, cand_c.id, ss),\
                            file=log)

                
                # Determine NL's between original losers and second winner
                # explanations: map between AG* and the NL's that they
                # are supporting. 
                explanations = {}
                for c in cands:
                    if c in winners:
                        continue

                    cand_c = candidates[c]

                    curr_nl = (sw.id, "NL", cand_c.id)

                    min_sw = 0 # Min tally for second winner.
                    max_c = 0 # Max tally for candidate c.

                    # Keep running tally of total votes we can increase the
                    # margin of assertion 'sw' NL 'c' with 'AG*' relationships
                    pot_margin_inc = 0

                    helpful_ags = []

                    # Assertion: min_sw > maxc
                    assorter = INVALID*0.5
                    for b in ballots:
                        if b.prefs != []:
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
                        if b.prefs != [] and b.prefs[0] == first_winner.num:
                            weight = aud_tv

                        c_idx = b.prefs.index(c) if c_in else np.inf
                        s_idx = b.prefs.index(sw.num) if s_in else np.inf

                        if c_in and c_idx < s_idx:
                            # These ballots could not go to 'sw' but may go
                            # to 'c'.

                            # Could we exclude these ballots by using an AG?
                            is_ag, ag_asn, descs = rule_out_for_max(b.prefs,c, \
                                ag_matrix, winners, candidates)
                            
                            contrib = b.votes*((1-weight)/2.0)

                            if is_ag:
                                alt_contrib = b.votes*0.5
                                helpful_ags.append((ag_asn,alt_contrib-contrib,\
                                    descs))
                                pot_margin_inc += alt_contrib-contrib

                            assorter += contrib
                            max_c += b.votes

                        elif s_in:
                            if min_tv > 0 and len(b.prefs) > 1 and \
                                b.prefs[:2] == [first_winner.num, sw.num]:
                                min_sw += min_tv * b.votes
                                assorter += b.votes * ((1 + min_tv)/2.0)
                            else:
                                # If we remove all cand 'd' for which sw AG* d
                                # that appear before sw in the ballot, will 'sw'
                                # be the first ranked cand? If so, we could add
                                # these ballots to the min tally of 'sw'.
                                prefs = b.prefs[:]

                                descs = set()
                                max_ags_here = 0
                                            
                                for d,dval in ags.items():
                                    if d in prefs:
                                        idx_d = prefs.index(d)
                                        if idx_d < s_idx:
                                            prefs.remove(d)
                                            s_idx -= 1
                                            rag = (second_winner.id, "AG*", \
                                                candidates[d].id, dval)
                                            descs.add(rag)
                                            if rag in explanations:
                                                explanations[rag].add(curr_nl)
                                            else:
                                                explanations[rag] = \
                                                    set([curr_nl])
                                            max_ags_here=max(max_ags_here,dval)

                                assorter += 0.5*b.votes

                                if prefs != [] and prefs[0] == sw.num:
                                    # These ballots would have originally had a 
                                    # contribution of 0.5 to the assorter. By
                                    # utilising these AG*'s we can increase the 
                                    # contribution of these ballots to 1 each.
                                    helpful_ags.append((max_ags_here,  \
                                        0.5*b.votes, descs))

                                    pot_margin_inc += b.votes*0.5

                                elif min_tv > 0 and  len(prefs) > 1 and \
                                    prefs[:2] ==  [first_winner.num, sw.num]:
                                    base_contrib = 0.5*b.votes
                                    alt_contrib = b.votes * ((1 + min_tv)/2.0)
                                    dconfig = alt_contrib - base_contrib

                                    if dconfig > 0:
                                        helpful_ags.append((max_ags_here, \
                                            dconfig, descs))

                                        pot_margin_inc += dconfig

                        else:
                            assorter += 0.5*b.votes

                    # Max ASN of any AG*'s used to increase assorter margins 
                    # when forming NLs.
                    max_ags_used = 0  
                    merged_helpful_ags = merge_helpful_ags(helpful_ags, \
                        pot_margin_inc)

                    # Incorporate use of  AG*'s that either make the assertion
                    # possible, or whose ASN is already within/equal to current
                    # lower bound on audit difficulty.
                    while assorter/args.voters <= 0.5 and \
                        merged_helpful_ags != []:
                        
                        ag_asn, extra_contrib, descs = merged_helpful_ags.pop(0)

                        assorter += extra_contrib
                        max_ags_used = max(max_ags_used, ag_asn)

                        inner_loop_assertions.update(descs)

                        for ag in descs:
                            if ag in explanations:
                                explanations[ag].add(curr_nl)
                            else:
                                explanations[rag] = set([curr_nl])


                    # Can we reduce the sample size required for the assertion
                    # by including more AG*'s?
                    amean = assorter/args.voters
                    if amean > 0.5:
                        # Current sample size for assertion
                        ss = sample_size(amean, args) 

                        # Go through remaining AG*'s 
                        for ag_asn, extra_contrib, descs in merged_helpful_ags:
                            if ss < ag_asn:
                                break

                            # would reducing margin by votes reduce sample size
                            # by more than increase caused by including AG*?
                            assorter += extra_contrib
                            amean = assorter/args.voters

                            ss  = sample_size(amean, args)
                            max_ags_used = max(max_ags_used, ag_asn)
                        
                            inner_loop_assertions.update(descs)

                            for ag in descs:
                                if ag in explanations:
                                    explanations[ag].add(curr_nl)
                                else:
                                    explanations[ag] = set([curr_nl])
                                
                        max_with_nls = max(max_with_nls, max(max_ags_used,ss))
                        print("NL({},{}) = {}, AG*'s used {}".format(sw.id, \
                            cand_c.id, ss, max_ags_used), file=log)
    
                        inner_loop_assertions.add((sw.id, "NL", cand_c.id, ss))
                    else:
                        max_with_nls = np.inf
                        print("NL({},{}) NOT POSSIBLE".format(\
                            sw.id, cand_c.id),file=log)

                max_this_loop=max(max_this_loop,max_with_nls)

                print("Max in loop {}, {}".format(max_this_loop, aud_tv), \
                    file=log)

                if max_this_loop <= max_in_loop:
                    best_inner_assertions = inner_loop_assertions
                    best_inner_explanations = explanations
                    max_in_loop = max_this_loop
                    aud_tv += deltat
                else:
                    break

            
            if max_in_loop <= max_in_outer_loop:
                max_in_outer_loop = max_in_loop
                best_outer_assertions = lts
                best_outer_assertions.update(best_inner_assertions)
                best_outer_exp = best_inner_explanations
            else:
                break    

            if min_tv == 0:
                min_tv = act_tv/2.0
            else:
                min_tv += deltat

        assertions_used.update(best_outer_assertions)
        assertions_text = [(w,n,l) for w,n,l,_ in assertions_used]

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

    else: 
        # CASE: 2 seats and no candidate has a quota on first preferences
        # according to reported results.  

        remcase = outcome.action[ncand-2:] == [1,1]

        for w in winners:
            print("Winner {}".format(candidates[w].id), file=log)

        # First look at whether we can batch eliminate some candidates at
        # the start of counting. If we can show that E are eliminated first,
        # before any candidate can get a quota, then for all other candidates
        # 'c' (not in E) we can give 'c' votes from any candidate in E in
        # their minimum tally as part of any NL calculation. 
        batch_asn = np.inf

        fw = winners[0]

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
                        possibilities[o].append((ss,(cand_c.id,"AG",\
                            cand_o.id,ss)))

                        print("{} AG {} = {}".format(candidates[c].id, \
                            candidates[o].id, ss), file=log)


        # Can we establish a winner based on the fact that they are AG
        # at least NCAND - 2 candidates? If there is one candidate that 
        # we cannot establish an AG relationship with that is OK. Then
        # we can rule out all alternate outcomes that do not involve that
        # known winner.
        known_winners = []

        for w in winners:
            # for how many candidates c can we show that w AG c
            cost = 0
            num_ag = 0
        
            descs = []

            w_cand = candidates[w]
            for c in cands:
                if c == w:
                    continue

                if ag_matrix[w][c] != None:
                    num_ag += 1
                    cost = max(cost, ag_matrix[w][c])
                    descs.append((w_cand.id,"AG",candidates[c].id,\
                        ag_matrix[w][c]))

            if num_ag >= ncand - 2:
                print("Candidate {} is a known winner, cost {}".format(\
                    candidates[w].id, cost), file = log)
                known_winners.append((w,cost))

                assertions_used.update(descs)

        if len(known_winners) == args.seats:
            print("All winners can be established with AG relationships.",\
                file=log)
            print("ASN: {}".format(max([cost for _,cost in known_winners])),\
                file=log)

            print("Final set of assertions required:", file=log)
            for asstn in assertions_used:
                print(asstn, file=log)
            
            log.close()
            exit(0)   

        # Now, we are looking for candidates c where there are at least
        # 2 other candidates o such that o AG c. This is the case where
        # possibilities[c] contains at least two sample sizes.
        for c in cands:
            if len(possibilities[c]) > 1:
                sorted_poss = sorted(possibilities[c])
                
                ss1, d1 = sorted_poss[0]
                ss2, d2 = sorted_poss[1]
    
                sample = max(ss1, ss2)

                assertions_used.add(d1)
                assertions_used.add(d2)

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
        asn_overall = max(max_sample_size, max_rule_out_ss)

        # Form pairs of alternate winner outcomes
        pairs = []
        assertions = set()

        must_win = None
        if known_winners != []:
            must_win = known_winners[0][0]
            asn_overall = max(asn_overall, known_winners[0][1]) 

        for i in range(ncand):
            c = cands[i]
            if c in ruled_out: continue

            if args.cfw and c == fw: continue

            for j in range(i+1, ncand):
                o = cands[j]

                if o in ruled_out: continue
                if c in winners and o in winners: continue

                if args.cfw and o == fw: continue

                if must_win != None and c != must_win and o != must_win:
                    continue

                pairs.append((c,o))

        cannot_rule_out = []
        for c1,c2 in pairs:
            cand_c1 = candidates[c1]
            cand_c2 = candidates[c2]

            print("ALTERNATE OUTCOME PAIR: {},{}".format(cand_c1.id, \
                cand_c2.id), file=log)

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
            best_asn = np.inf 

            failure = True

            best_assertions = set()

            for o in cands:
                if o == c1 or o == c2 or o in ruled_out: continue
           
                ctvmax1, ctvmax_ss1 = args.maxtv, 0 #max_tv_by_z(c1, \
                    #c2, o, cands, candidates, ballots, 0.005, args, log, \
                    #valid_ballots, INVALID, asn_overall, ag_matrix)

                ctvmax2, ctvmax_ss2 = args.maxtv, 0 #max_tv_by_z(c2, \
                    #c1, o, cands, candidates, ballots, 0.005, args, log, \
                    #valid_ballots, INVALID, asn_overall, ag_matrix)

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

                c1_assertions = set()
                asn_c1 = np.inf

                # Max ASN of any AGs used to increase/decrease the minimum
                # /maximum tallies of second winner/candidate c when forming NL.
                max_ags_used1 = ctvmax_ss1 

                # Consider max vote of 'c2' given 'c1' seated at some stage and 
                # 'o' is still standing.
                assorter_vals = [INVALID*0.5]

                mint, minpapers, maxt, maxpapers = 0, 0, 0, 0
                for b in ballots:
                    awarded, used_ss, descs = vote_for_cand_ags1(o, b.prefs,\
                        o_ag, candidates)

                    if awarded:
                        if used_ss != None and used_ss > 0:
                            max_ags_used1 = max(max_ags_used1, used_ss)
                            c1_assertions.update(descs)
                            
                        mint += b.votes
                        minpapers += b.votes
                        
                        assorter_vals.append(b.votes)
                        continue    
                        
                    if not c2 in b.prefs: 
                        assorter_vals.append(b.votes*0.5)
                        continue

                    weight = 1.0

                    if b.prefs != [] and b.prefs[0] == c1:
                        weight = ctvmax1

                    # In this analysis, we are assuming that no one other
                    # than c1 and c2 has won -- so if AG(d, c2) and 'd'
                    # is not c1, and 'd' is preferenced before c2, then
                    # c2 does not get these ballots. 
                    awarded, ag_present, used_ss, descs = vote_for_cand_ags2(\
                        c2, o, b.prefs, ag_matrix, [c1], candidates)

                    if awarded:
                        if ag_present:
                            assorter_vals.append(b.votes*0.5)
                            max_ags_used1 = max(max_ags_used1, used_ss)
                            c1_assertions.update(descs)
                        else:
                            contrib = b.votes*((1 - weight)/2.0)
                            assorter_vals.append(contrib)
                            maxt += weight*b.votes
                            maxpapers += b.votes
                    else:
                        assorter_vals.append(b.votes*0.5)

                amean = sum(assorter_vals)/args.voters


                print("   (1) can we show that {} NL {}? ".format(cand_o.id,\
                    cand_c2.id), file=log)
                print("        a. margin {}".format(2*amean - 1), file=log)

                # Is the minimum tally for 'o' larger than the maximum
                # possible tally for 'c2'? This means 'o' cannot be 
                # eliminated before 'c2'
                if amean > 0.5:
                    ss = sample_size(amean, args)

                    c1_assertions.add("{} NL {} = {} [rules out {}]".format(\
                        cand_o.id, cand_c2.id, ss, (cand_c1.id, cand_c2.id)))

                    if ss != np.inf:
                        asn_c1 = max(ss, max_ags_used1)

                        failure = False
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

                c2_assertions = set()
                asn_c2 = np.inf

                # Consider max vote of 'c1' given 'c2' seated at some stage and 
                # 'o' is still standing.
                mint, minpapers, maxt, maxpapers = 0, 0, 0, 0

                assorter_vals = [INVALID*0.5]

                # Max ASN of any AGs used to increase/decrease the minimum/
                # maximum tallies of second winner/candidate c when forming NL.
                max_ags_used2 = ctvmax_ss2 

                for b in ballots:
                    awarded, used_ss, descs = vote_for_cand_ags1(o, b.prefs, \
                        o_ag, candidates)

                    if awarded:
                        if used_ss != None and used_ss > 0:
                            max_ags_used2 = max(max_ags_used2, used_ss)
                            c2_assertions.update(descs)
                            
                        minpapers += b.votes
                        mint += b.votes
                        assorter_vals.append(b.votes)
                        continue    
                        
                    if not c1 in b.prefs: 
                        assorter_vals.append(b.votes*0.5)
                        continue

                    weight = 1.0

                    if b.prefs != [] and b.prefs[0] == c2:
                        weight = ctvmax2

                    awarded, ag_present, used_ss, descs = vote_for_cand_ags2(\
                        c1, o, b.prefs, ag_matrix, [c2], candidates)

                    if awarded:
                        if ag_present:
                            assorter_vals.append(b.votes*0.5)
                            max_ags_used2 = max(max_ags_used2, used_ss)
                            c2_assertions.update(descs)

                        else:
                            contrib = b.votes*((1 - weight)/2.0)
                            assorter_vals.append(contrib)
                            maxt += weight*b.votes
                            maxpapers += b.votes
                    else:
                        assorter_vals.append(b.votes*0.5)

                amean = sum(assorter_vals)/args.voters
 
                print("   (1) can we show that {} NL {}? ".format(cand_o.id,\
                    candidates[c1].id), file=log)
                print("        a. margin {}".format(2*amean - 1), file=log)


                # Is the minimum tally for 'o' larger than the maximum
                # possible tally for 'c1'? This means 'o' cannot be 
                # eliminated before 'c1'
                if amean > 0.5:
                    ss = sample_size(amean, args)
                    
                    c2_assertions.add("{} NL {} = {} [rules out {}]".format(\
                        cand_o.id, cand_c1.id, ss, (cand_c1.id, cand_c2.id)))

                    if ss != np.inf:
                        asn_c2 = max(ss, max_ags_used2)
                        failure = False
                        print("      yes, {}, ags used {} = {}".format(ss, \
                             max_ags_used2, max(ss, max_ags_used2)), file=log)
                    else:
                        print("      no", file=log)
                else:
                    print("      no, min {:.2f},{} vs max {:.2f},{}".format(\
                        100*(mint)/args.voters, minpapers, \
                        100*(maxt)/args.voters, maxpapers), file=log)

                if asn_c1 < best_asn: 
                    best_assertions = c1_assertions
                    best_asn = asn_c1

                if asn_c2 < best_asn:
                    best_assertions = c2_assertions
                    best_asn = asn_c2
                    
            asn_overall = max(asn_overall, best_asn)

            if failure:
                cannot_rule_out.append((c1,c2))
            else:
                print("Best assertion(s) require sample: {}".format(\
                    best_asn), file=log)

                assertions_used.update(best_assertions)
                 
          
        for (c1,c2) in cannot_rule_out:
            c1o = candidates[c1]
            c2o = candidates[c2]
            print("Alternate outcome cannot be ruled out: {},{}".format(\
                c1o.id, c2o.id), file=log)

        if asn_overall >= args.voters:
            asn_overall = np.inf
        
        print("Best sample size: {}".format(asn_overall), file=log)


        print("------------------------------------------------",file=log)
        print("Final set of assertions generated:", file=log)
        for asstn in set(assertions_used):
            print(asstn, file=log)

        if cannot_rule_out == []:
            print("All winners are verified in audit.", file=log)
        else:
            known_winners = set([w for w,_ in known_winners])

            for w in winners:
                in_all = True
                for (c1,c2) in cannot_rule_out:
                    if w != c1 and w != c2:
                        in_all = False
                        break
                if in_all:
                    known_winners.add(w)

            if known_winners != []:
                print("{} are verified winners with these assertions.".format(\
                    [candidates[w].id for w in known_winners]), file=log)

            candidate_winners = set(winners)
            for (c1,c2) in cannot_rule_out:
                candidate_winners.add(c1)
                candidate_winners.add(c2)
            
            candidate_winners = set([w for w in candidate_winners if \
                not w in known_winners])

            if candidate_winners != []:
                print("{} remain as potential winners in this audit.".format(\
                    [candidates[w].id for w in candidate_winners]), file=log)
    
        print("------------------------------------------------",file=log)
 

        if remcase:
            print("CASEC,{},{},{},{},{}".format(args.data, ncand, \
                valid_ballots, args.quota, asn_overall))
        else:
            print("CASEB,{},{},{},{},{}".format(args.data, ncand, \
                valid_ballots, args.quota, asn_overall))

    log.close() 
