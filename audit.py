import argparse
import os


from utils import read_ballots_stv, read_outcome, subsupermajority_sample_size,\
    cand_vs_cand_sample_size, read_ballots_txt


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    # Input: stv data file
    parser.add_argument('-d', dest='data')

    # Input: anticipated error rate (default value is 0)
    parser.add_argument('-e', dest='erate', default=0, type=float)

    # Input: risk limit (default is 5%)
    parser.add_argument('-r', dest='rlimit', default=0.05, type=float)
    
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

    # Input: Allowed shift in transfer value 
    parser.add_argument('-shift', dest='shift', type=float)

    # Input: Modelling of denominator in transfer value equation
    parser.add_argument('-dchoice', dest='dchoice', type=int, default=0)
    
    # Output: Log file 
    parser.add_argument('-log', dest='log', type=str)

    args = parser.parse_args()

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

    assert(args.seats == 2)
    assert(outcome.action[0] == 1)

    INVALID = args.voters - valid_ballots

    assertions_na = False

    # Form assertion (and estimate sample size required to check) that 
    # the first winner acheived a quota on first preferences.
    threshold = args.quota / valid_ballots

    first_winner = candidates[outcome.cand[0]]
    second_winner = None

    for i in range(1, len(outcome.cand)):
        if outcome.action[i] == 1:
            second_winner = candidates[outcome.cand[i]]
            break

    assert(second_winner != None)

    max_sample = 0
    sample_size, m, mub = subsupermajority_sample_size(threshold, \
        first_winner.fp_votes, INVALID, args)

    max_sample = max(max_sample, sample_size)

    with open(args.log, "w") as log:
        print("{},{},{},{}".format(sample_size, m, mub,\
            "Check first winner > Quota"), file=log)

        if outcome.action[1] == 1 and second_winner.fp_votes > args.quota:
            # Our second candidate also achieves a quota in the first round
            sample_size, m, mub = subsupermajority_sample_size(threshold, \
                second_winner.fp_votes, INVALID, args)
            print("{},{},{},{}".format(sample_size, m, mub, \
                "Check second winner > Quota"), file=log)
            max_sample = max(max_sample, sample_size)
            print("{},{},{},{},NA,NA,{}".format(args.data, args.erate, \
                args.rlimit,args.shift,max_sample))

            exit()
    
        # Look at difficulty of assertions that state that all other candidates
        # do not achieve a quota in the first around. We do not need to do
        # this for the second winner.
        for c in candidates:
            if c.num == first_winner.num or c.num == second_winner.num:
                continue

            # For assertion (and check sample size) that this candidate did 
            # not achieve a quota in the first round. 
            threshold = 1 - ((args.quota - 1) / valid_ballots)

            tally_others = valid_ballots - c.fp_votes

            sample_size, m, mub = subsupermajority_sample_size(threshold, \
                tally_others, INVALID, args)

            print("{},{},{},{}".format(sample_size, m, mub, \
                "Check {} does not have a quota in first round".format(c.num)),\
                file=log)
            max_sample = max(max_sample, sample_size)

        # Now look at the difficulty of showing that the second winner could
        # not be eliminated before anyone else remaining -- with the first
        # winner already elected. 

        # Assume case where second winner gets transferred votes from first
        # winner
    
        # First, establish a lower and upper bound on the transfer value of 
        # votes from winner 1.

        # Reported surplus
        rep_surplus = first_winner.fp_votes - args.quota

        # Denominator of transfer value equation.
        rep_tfable = 0

        for bnum in first_winner.ballots:
            prefs = ballots[bnum].prefs

            if args.dchoice == 0 and prefs != [first_winner.num]:
                rep_tfable += ballots[bnum].votes
            elif args.dchoice >= 1:
                rep_tfable += ballots[bnum].votes

        rep_tv = min(1.0, rep_surplus/rep_tfable)

        rep_tv_lb = min(1.0, max(0, rep_tv - args.shift))
        rep_tv_ub = min(1.0, max(0, rep_tv + args.shift))

        print("Reported transfer value for winner 1: {}".format(rep_tv),\
            file=log)
        print("Reported bounds on transfer value for winner 1: {},{}".format(\
            rep_tv_lb, rep_tv_ub), file=log)

        # Compute difficultly of auditing the bounds on transfer values
        # We need to show that tally_first_winner >= rep_surplus_lb+quota,
        # tally_first_winner <= rep_surplus_ub+quota.

        shiftv = (rep_surplus - rep_tv_lb * rep_tfable)/(1 + rep_tv_lb);
        rep_surplus_lb = max(0, rep_surplus - shiftv)
        rep_surplus_ub = rep_surplus + shiftv

        # Check first winner surplus is above a lower bound
        threshold = (rep_surplus_lb + args.quota)/valid_ballots

        sample_size, m, mub = subsupermajority_sample_size(threshold, \
            first_winner.fp_votes, INVALID, args)
        
        print("{},{},{},{},{}".format(sample_size, m, mub, threshold,\
            "Check first winner surplus >= lower bound"), file=log)
        max_sample = max(max_sample, sample_size)

        if rep_tv < 1:
            # Check first winner surplus is below an upper bound
            # We want to check that tally_first_winner <= quota+rep_surplus_lb
            #
            # This would mean that all other candidates should get at least
            # valid_ballots - (quota+rep_surplus+shift) + 1 votes.
            rep_tally_others = valid_ballots - first_winner.fp_votes

            threshold = 1 - ((args.quota + rep_surplus_ub)/valid_ballots) 

            sample_size, m, mub = subsupermajority_sample_size(threshold, \
                rep_tally_others, INVALID, args)
            max_sample = max(max_sample, sample_size)
        
            print("{},{},{},{},{}".format(sample_size, m, mub, threshold,\
                "Check first winner surplus <= upper bound"), file=log)

        # Check that number of transferrable ballots in first winner's tally 
        # is greater than a lower bound.
        if rep_tv < 1:
            threshold = max(0.0, (rep_tfable - shiftv))/valid_ballots
   
            if threshold > 0: 
                sample_size, m, mub = subsupermajority_sample_size(threshold, \
                    rep_tfable, INVALID, args)
                max_sample = max(max_sample, sample_size)
        
                print("{},{},{},{},{}".format(sample_size, m, mub, threshold,\
                    "Check first winner tf ballots >= lower bound"),file=log)
        
        # What would transferrable ballots <= upper bound look like?
        # Imagine we have two "virtual candidates". The first gets awarded
        # all the ballots that have the first winner as the first
        # preference, and that can be transferred to someone else. The
        # second gets all other ballots. For the first category of ballots
        # to be less than or equal to the upper bound rep_tfable +
        # args.shift, we want the second category to contain at least
        # valid_ballots - (rep_tfable+args.shift) + 1 ballots.
        threshold = 1 - ((rep_tfable + shiftv)/valid_ballots)

        sample_size, m, mub = subsupermajority_sample_size(threshold, \
            valid_ballots - rep_tfable, INVALID, args)
        
        print("{},{},{},{},{}".format(sample_size, m, mub, threshold,\
            "Check first winner tf ballots <= upper bound"),file=log)
        max_sample = max(max_sample, sample_size)

        # Now look at assertions to show that no remaining candidate could
        # have been eliminated before the second winner

        # Find min tally of second_winner: their first preferences + the minimum
        # transfer from first_winner
        tally_sw_c1 = second_winner.fp_votes

        for b in first_winner.ballots:
            # Is second_winner the next most preferred candidate?
            ballot = ballots[b]

            if len(ballot.prefs) > 1 and ballot.prefs[1] == second_winner.num:
                tally_sw_c1 += rep_tv_lb * ballot.votes

        for cand in candidates:
            if cand.num == first_winner.num or cand.num == second_winner.num:
                continue  

            # Assertions to show that second_winner could not have been
            # eliminated before cand.

            # Find maximum tally of cand -- include their first preference votes
            # and any votes that preference cand before second_winner
            tally_cand_c1 = cand.fp_votes

            for blt in ballots:
                if blt.prefs[0] == cand.num:
                    continue

                fw_before_c = False
                sw_before_c = False
                c_present = False

                for p in blt.prefs:
                    if p == cand.num:
                        c_present = True
                        break

                    if not c_present:
                        if p == first_winner.num:
                            fw_before_c = True

                        elif p == second_winner.num:
                            sw_before_c = True

                if c_present is False:
                    continue

                # if sw_before_c is True, don't give cand the votes.
                if sw_before_c is False:
                    if fw_before_c:
                        tally_cand_c1 += blt.votes * rep_tv_ub
                    else:
                        tally_cand_c1 += blt.votes

            ss1, m1, mub1 = cand_vs_cand_sample_size(tally_sw_c1, \
                tally_cand_c1, valid_ballots, args) 
            print("Second winner vs {}, {}, {}, {}, {},{}".format(cand.num,\
                tally_sw_c1, tally_cand_c1, m1, mub1,ss1), file=log)
            max_sample = max(max_sample, ss1) 

            if m1 < 0:
                assertions_na = True

        print("{},{},{},{},{},{},{},NEB N/A: {}".format(args.data, args.erate,\
            args.rlimit, args.shift, rep_tv_lb, rep_tv_ub, max_sample,\
            assertions_na))
              
