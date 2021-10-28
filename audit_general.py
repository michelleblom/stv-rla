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


from utils import read_ballots_stv, read_outcome, subsupermajority_sample_size,\
    cand_vs_cand_sample_size, read_ballots_txt, index_of


def vote_for_cand_over_cand(c1, c2, prefs):
    for p in prefs:
        if p == c1:
            return c1

        if p == c2:
            return c2

    return None



# 'c_neb' denotes the candidates 'd' for which candidate 'c' cannot be 
# eliminated before (ie. c NEB d). Is 'c' the first preferenced candidate 
# in the ranking 'prefs' if we exclude those in 'c_neb'. Keep track of the
# audit difficulty of any NEBs relied on.
def vote_for_cand_nebs1(c, prefs, c_neb, threshold):
    nebs_used = 0

    for p in prefs:
        if p == c:
            return True,nebs_used

        if p != c:
            if p in c_neb and (threshold != None and c_neb[p] <= threshold):
                nebs_used = max(nebs_used, c_neb[p])
            else:
                return False, None

    return False,None



# Does c1 get this vote over c2 assuming winning candidates are in
# list 'winners'
def vote_for_cand_nebs2(c1, c2, prefs, neb_matrix, winners, threshold):
    neb_present = False
    neb_min_ss = np.inf

    for p in prefs:
        if p == c1:
            if neb_present:
                return False,neb_min_ss
            else:
                return True,0

        if p == c2:
            return False,None

        if p in winners:
            # c1 could still get this vote as it could skip over 'p'
            continue

        neb = neb_matrix[p][c1]
        if neb != None and neb <= threshold: #np.inf:
            neb_present = True
            neb_min_ss = min(neb_min_ss, neb)

    return False,None



if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    # Input: stv data file
    parser.add_argument('-d', dest='data')

    # Input: anticipated error rate (default value is 0)
    parser.add_argument('-e', dest='erate', default=0.02, type=float)

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

    args = parser.parse_args()

    log = open(args.log, "w")

    assert(args.seats == 2)

    # Read STV data file
    candidates, ballots, id2group, cid2num, valid_ballots = \
        None, None, None, None, None

    if args.data.endswith(".stv"):
        candidates, ballots, id2group, cid2num, valid_ballots = \
            read_ballots_stv(args.data)
    else:
        candidates, ballots, id2group, cid2num, valid_ballots = \
            read_ballots_txt(args.data)

    # Read STV outcome file
    outcome = read_outcome(args.outcome, cid2num)

    INVALID = args.voters - valid_ballots

    ncand = len(outcome.cand)
    cands = []
    winners = []
    winners_on_fp = []
    losers = []

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

    # TODO: add parameter to control whether we want to run 2 quota 
    # method and replace 'False' below with that.
    if False and args.seats == 2 and outcome.action[0] == 1 and \
        outcome.action[1] == 1:

        max_sample_size = 0

        # Check that first winner's FP is greater than a quota
        first_winner = candidates[winners[0]]
        thresh = args.quota/valid_ballots
    
        ss, _, _ = subsupermajority_sample_size(thresh,first_winner.fp_votes,\
            INVALID, args)
        print("{} ballot checks required to assert that first winner"\
            " has a quota's worth of votes".format(ss), file=log)

        max_sample_size = max(max_sample_size, ss)

        # Check that second winner's FP is greater than a quota
        second_winner = candidates[winners[1]]
    
        ss, _, _ = subsupermajority_sample_size(thresh,second_winner.fp_votes,\
            INVALID, args)
        print("{} ballot checks required to assert that second winner"\
            " has a quota's worth of votes".format(ss), file=log)

        max_sample_size = max(max_sample_size, ss)

        print("2Q,{},{},{},{},{}".format(args.data, ncand, valid_ballots, \
            args.quota, max_sample_size))

    elif args.seats == 2 and outcome.action[0] == 1:
        # CASE: 2 seats and first winner has more than a quota on
        # first preferences.

        max_sample_size = 0

        # Check that first winner's FP is greater than a quota
        first_winner = candidates[winners[0]]
        thresh = args.quota/valid_ballots
    
        ss, _, _ = subsupermajority_sample_size(thresh,first_winner.fp_votes,\
            INVALID, args)
        print("{} ballot checks required to assert that first winner"\
            " has a quota's worth of votes".format(ss), file=log)

        max_sample_size = max(max_sample_size, ss)

        act_tv = (first_winner.fp_votes - args.quota)/first_winner.fp_votes
        aud_tv = act_tv + 0.01

        max_in_loop = np.inf

        neb_matrix = [[None for c in cands] for o in cands]

        for c in cands:
            fpc = candidates[c].fp_votes

            for o in cands:
                if c == o: continue

                max_vote_o = 0

                for b in ballots:
                    awarded = vote_for_cand_over_cand(c, o, b.prefs)

                    if awarded == o:
                        max_vote_o += b.votes

                if fpc > max_vote_o:
                    ss, m, _ = cand_vs_cand_sample_size(fpc, max_vote_o, \
                        valid_ballots, args) 

                    if ss != np.inf:
                        neb_matrix[c][o] = ss

        while aud_tv < 2/3:
            # Check that TV of first winner is at most aud_TV
            T = args.quota/(1 - aud_tv)

            # Check that first winner's FP is less than T
            threshold = 1 - (T - 1)/valid_ballots
            tally_others = valid_ballots - first_winner.fp_votes

            ss, _, _ = subsupermajority_sample_size(threshold, tally_others, \
                INVALID, args)

            max_this_loop = ss

            # For second winner, is it the case that they cannot be 
            # eliminated prior to all remaining candidates?
            sw = candidates[winners[1]]
            min_sw = sw.fp_votes

            nlts = {}

            max_with_nlts = ss
            max_with_nebs = ss

            # Compute NLTs between original losers and second winner
            for c in cands:
                if c in winners:
                    continue

                cand_c = candidates[c]
                max_c = cand_c.fp_votes

                for b in ballots:
                    if b.prefs[0] == sw.num or b.prefs[0] == c: 
                        continue

                    weight = 1
                    for p in b.prefs:
                        if p == sw.num:
                            break

                        if p == c:
                            max_c += b.votes * weight
                            break

                        if p in winners_on_fp:
                            weight = aud_tv

                if max_c < min_sw:
                    ss, m, _ = cand_vs_cand_sample_size(min_sw, max_c, \
                        valid_ballots, args) 

                    nlts[c] = ss

                    max_with_nlts = max(max_with_nlts, ss)
                else:
                    max_with_nlts = np.inf


            # Determine NEBs between original losers and second winner
            for c in cands:
                if c in winners:
                    continue

                cand_c = candidates[c]

                min_sw_w_extra = 0

                max_c = cand_c.fp_votes

                max_nlts_used = 0

                for b in ballots:
                    if b.prefs[0] == c or not c in b.prefs:
                        continue

                    weight = 1

                    if b.prefs[0] == first_winner.num:
                        weight = aud_tv

                    awarded, used_ss = vote_for_cand_nebs2(c, sw.num, b.prefs,\
                        neb_matrix, [first_winner.num], 500)

                    if awarded:
                        max_c += weight * b.votes

                        if used_ss != None:
                            max_nlts_used = max(max_nlts_used, used_ss)

                for b in ballots:
                    if not sw.num in b.prefs:
                        continue

                    idx_sw = b.prefs.index(sw.num)
                    if first_winner.num in b.prefs and \
                        b.prefs.index(first_winner.num) < idx_sw:
                        # this might be a bit conservative to not give
                        # the second winner votes from the first winner. 
                        # The only way these votes skip the second winner
                        # is if the second winner gets a quota. 
                        continue

                    prefs = b.prefs[:]

                    # remove candidates 'd' for which sw NLT 'd' where
                    # d != any of {'c', sw.num}
                    for d,dval in nlts.items():
                        if d != c and d in prefs and dval < 500:
                            idx_d = prefs.index(d)
                            if idx_d < idx_sw:
                                prefs.remove(d)
                                idx_sw -= 1
                                max_nlts_used = max(max_nlts_used, dval)

                    if prefs[0] == sw.num:
                        min_sw_w_extra += b.votes      
                        

                if max_c < min_sw_w_extra:
                    ss, m, _ = cand_vs_cand_sample_size(min_sw_w_extra, max_c,\
                        valid_ballots, args) 

                    max_with_nebs = max(max_with_nebs, max(max_nlts_used,ss))

                else:
                    max_with_nebs = np.inf

            max_this_loop=max(max_this_loop,min(max_with_nlts,max_with_nebs))

            print("Max in loop {}, {}".format(max_this_loop, aud_tv), file=log)
            if max_this_loop < max_in_loop:
                max_in_loop = max_this_loop
                act_tv = aud_tv
                aud_tv += 0.01
            else:
                break

        max_sample_size = max(max_sample_size, max_in_loop)
    
        if max_in_loop is np.inf:
            print("Not possible to show that second winner could not be "\
                "eliminated prior to all remaining candidates.", file=log)
        else:
            print("Demonstrating that first winner's TV is less"\
                " than {}".format(act_tv), file=log)

        print("Sample size required for audit is {} ballots".format(\
            max_sample_size), file=log)

        print("1Q,{},{},{},{},{}".format(args.data, ncand, valid_ballots, \
            args.quota, max_sample_size))

    else: # args.seats == 2 and outcome.action[0] == 0:
        # CASE: 2 seats and no candidate has a quota on first preferences
        # according to reported results.
        neb_matrix = [[None for c in cands] for o in cands]

        rule_out = []
        max_rule_out_ss = 0
        ruled_out = []

        possibilities = [[] for c in cands]
        for c in cands:
            fpc = candidates[c].fp_votes

            for o in cands:
                if c == o: continue

                max_vote_o = 0

                for b in ballots:
                    awarded = vote_for_cand_over_cand(c, o, b.prefs)

                    if awarded == o:
                        max_vote_o += b.votes

                if fpc > max_vote_o:
                    ss, m, _ = cand_vs_cand_sample_size(fpc, max_vote_o, \
                        valid_ballots, args) 

                    if ss != np.inf:
                        neb_matrix[c][o] = ss
                        possibilities[o].append(ss)

                        print("{} NEB {} = {}".format(candidates[c].id, \
                            candidates[o].id, ss), file=log)

        for c in cands:
            if len(possibilities[c]) > 1:
                sorted_poss = sorted(possibilities[c])
                
                sample = max(sorted_poss[0],sorted_poss[1])

                # We can rule candidate 'c' out with sample size 'sample'
                rule_out.append((candidates[c].id, sample))
                ruled_out.append(c)
                max_rule_out_ss = max(max_rule_out_ss, sample)

        print("Ruling out candidates: {}".format(rule_out), file=log)
        print("Sample size for ruling out cands: {}".format(max_rule_out_ss),\
            file=log)

        # Form pairs of alternate winner outcomes
        pairs = []
        assertions = []

        for i in range(ncand):
            c = cands[i]
        
            if c in ruled_out: continue

            for j in range(i+1, ncand):
                o = cands[j]

                if o in ruled_out: continue

                if c in winners and o in winners: continue

                pairs.append((c,o))

        assertions = []
        for c1,c2 in pairs:
            print("ALTERNATE OUTCOME PAIR: {},{}".format(candidates[c1].id,\
                candidates[c2].id), file=log)

            # If we assume 'c1' gets a seat at some stage, does there
            # exist a candidate 'o' \= 'c1' that cannot be eliminated before 
            # 'c2'? If so, 'c2' cannot get a seat.

            best_asn1 = np.inf

            best_asn2 = np.inf

            not_applicable = True

            for o in cands:
                if o == c1 or o == c2 or o in ruled_out: continue
               
                # Find set of candididates 'c' for which o NEB c.
                # Keep track of cost of each of those NEBs
                cand_o = candidates[o]
                o_neb = {}
                for c in cands:
                    if c != c1 and c != c2 and neb_matrix[o][c] != None \
                        and neb_matrix[o][c] != np.inf:
                        o_neb[c] = neb_matrix[o][c]

                min_tally_o = 0
                used_nebs_mt_o = 0
 
                used_nebs1 = 0

                # Max vote of 'c2' given 'c1' seated at some stage and 
                # 'o' is still standing.
                max_vote_c2 = 0

                for b in ballots:
                    awarded, used_ss = vote_for_cand_nebs1(o, b.prefs, o_neb,
                        500)

                    if awarded:
                        min_tally_o += b.votes
                        
                        if used_ss != None:
                            used_nebs_mt_o = max(used_nebs_mt_o, used_ss)
                        continue

                    if not c2 in b.prefs: 
                        continue

                    weight = 1

                    if b.prefs[0] == c1:
                        weight = args.maxtv
                   
                    # In this analysis, we are assuming that no one other
                    # than c1 and c2 has won -- so if NEB(d, c2) and 'd'
                    # is not c1, and 'd' is preferenced before c2, then
                    # c2 does not get these ballots. 
                    awarded, used_ss = vote_for_cand_nebs2(c2, o, b.prefs,\
                        neb_matrix, [c1], 500)

                    if awarded:
                        max_vote_c2 += weight * b.votes

                        if used_ss != None:
                            used_nebs1 = max(used_nebs1, used_ss)

                used_nebs1 = max(used_nebs1, used_nebs_mt_o)

                print("   (1) can we show that {} NEB* {}? ".format(cand_o.id,\
                    candidates[c2].id), file=log)
                print("      min tally {} is {}".format(cand_o.id, \
                    min_tally_o), file=log)
                print("      max tally {} is {}".format(candidates[c2].id, \
                    max_vote_c2), file=log)


                # Is the minimum tally for 'o' larger than the maximum
                # possible tally for 'c2'? This means 'o' cannot be 
                # eliminated before 'c2'
                if min_tally_o > max_vote_c2:
                    ss, m, _ = cand_vs_cand_sample_size(min_tally_o, \
                        max_vote_c2, valid_ballots, args) 

                    best_asn1 = min(max(ss,used_nebs1), best_asn1)
                    not_applicable = False
                    print("      yes, {}, nebs used {} = {}".format(ss, \
                        used_nebs1, max(ss,used_nebs1)), file=log)

                else:
                    print("      no", file=log)

                # Max vote of 'c1' given 'c2' seated at some stage and 
                # 'o' is still standing.
                max_vote_c1 = 0

                used_nebs2 = used_nebs_mt_o

                for b in ballots:
                    if not c1 in b.prefs: continue

                    weight = 1

                    if b.prefs[0] == c2:
                        weight = args.maxtv
                    
                    awarded, used_ss = vote_for_cand_nebs2(c1, o, b.prefs,\
                        neb_matrix, [c2], 500)

                    if awarded:
                        max_vote_c1 += weight * b.votes

                        if used_ss != None:
                            used_nebs2 = max(used_nebs2, used_ss)

                print("   (2) can we show that {} NEB* {}? ".format(cand_o.id,\
                    candidates[c1].id), file=log)
                print("      min tally {} is {}".format(cand_o.id, \
                    min_tally_o), file=log)
                print("      max tally {} is {}".format(candidates[c1].id, \
                    max_vote_c1), file=log)


                # Is the minimum tally for 'o' larger than the maximum
                # possible tally for 'c1'? This means 'o' cannot be 
                # eliminated before 'c1'
                if min_tally_o > max_vote_c1:
                    ss, m, _ = cand_vs_cand_sample_size(min_tally_o, \
                        max_vote_c1, valid_ballots, args) 

                    best_asn2 = min(max(ss,used_nebs2), best_asn2)
                    not_applicable = False
                    print("      yes, {}, nebs used {} = {}".format(ss, \
                        used_nebs2, max(ss,used_nebs2)), file=log)

                else:
                    print("      no", file=log)
            
            best_asn = min(best_asn1, best_asn2)
            assertions.append((best_asn, (c1,c2), not_applicable)) 

            if not not_applicable:
                print("Best assertion(s) require sample: {}".format(\
                    best_asn), file=log)
                   
        best_asn = max_rule_out_ss
        for asn,(c1,c2),na in assertions:
            if na:
                c1o = candidates[c1]
                c2o = candidates[c2]
                print("Alternate outcome cannot be ruled out: {},{}".format(\
                    c1o.id, c2o.id), file=log)

            best_asn = max(best_asn, asn)

        print("Best sample size: {}".format(best_asn), file=log)
        #print("{}, Sample size: {}".format(args.data, best_asn))
        print("CASEB,{},{},{},{},{}".format(args.data, ncand, valid_ballots,\
            args.quota, best_asn))

    log.close() 
