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


# Importing code from Philip's SHANGRLA repository:
# https://github.com/pbstark/SHANGRLA/tree/main/Code
from assertion_audit_utils import TestNonnegMean

import os
import numpy as np
import statistics
import math
import re
import json

class Ballot:
    def __init__(self, num, votes, prefs, atl=False):
        self.num = num
        self.votes = votes
        self.prefs = prefs[:]
        self.atl = atl

        self.papers = votes

class Candidate:
    def __init__(self, num, idn):
        self.num = num
        self.id = idn
        self.name = None
        self.group_id = None
        self.position = None

        self.ballots = []
        self.fp_votes = 0

        self.mentions = []

        # For simulation purposes
        self.sim_votes = 0
        self.max_votes = 0
        self.bweights = []
        self.standing = 1
        self.seat = -1
        self.surplus = 0

class Outcome:
    def __init__(self):
        self.cand = []
        self.action = []

class Group:
    def __init__(self, idstr):
        self.id = idstr
        self.cands = []


def read_outcome(path, cid2num):
    outcome = Outcome()

    with open(path, "r") as otc:
        lines = otc.readlines()

        # First line contains candidate numbers
        # in order of seating/elimination. After
        # all seats have been filled, the remaining
        # candidates will be listed in no particular
        # order.
        toks = lines[0].strip().split(',')

        outcome.cand = [cid2num[int(c)] for c in toks]

        # The second line contains a series of 0s or
        # 1s, separated by commas, indicating whether
        # a seating or elimination happened in that
        # position.
        toks = lines[1].strip().split(',')

        outcome.action = [int(a) for a in toks]

    return outcome

def read_ballots_txt(path):
    ballots = []
    candidates = []
    cid2num = {}

    total_votes = 0
    us_ver = True if path.endswith(".us") else False

    with open(path, "r") as cvr:
        lines = cvr.readlines()

        clist = [int(c.strip()) for c in lines[0].strip().split(',')]
        ncands = len(clist)

        for i in range(len(clist)):
            cand = Candidate(i, clist[i])
            cand.name = ""
            cand.group_id = -1
            cand.position = -1

            candidates.append(cand)
            cid2num[clist[i]] = i


        bcntr = 0
        nextline = 2 if us_ver else 5

        for i in range(nextline, len(lines)):
            line = lines[i].strip()
            toks = [t.strip() for t in line.split(':')]

            pre_prefs = toks[0][1:-1].split(',') 

            prefs = [int(p.strip()) for p in pre_prefs if p != '']

            if prefs == []:
                continue

            votes = int(toks[1])

            cprefs = [cid2num[p] for p in prefs]
            ballot = Ballot(bcntr, votes, cprefs, atl=False)
            ballots.append(ballot)

            fpcand = candidates[cprefs[0]]
            fpcand.ballots.append(bcntr)
            fpcand.fp_votes += votes

            total_votes += votes

            for p in cprefs:
                candidates[p].mentions.append(bcntr)

            bcntr += 1

    return candidates,ballots,{},cid2num,total_votes

def read_ballots_json(path):
    ballots = []
    candidates = []
    id2group = {}
    cid2num = {}

    total_votes = 0

    with open(path, "r") as cvr:
        data = json.load(cvr)

        # get candidates
        cid = 0
        for cand in data["metadata"]["candidates"]:
            name = cand["name"]
            party = int(cand["party"])
            pos = int(cand["position"])

            cobj = Candidate(cid, cid)
            cobj.name = name
            cobj.group_id = party
            cobj.position = pos

            candidates.append(cobj)

            # Candidate id to number mapping not needed for this
            # file type, but required for consistency across file
            # types
            cid2num[cid] = cid
            cid += 1


        # get party info
        pid = 0
        for party in data["metadata"]["parties"]:
            group = Group(pid)

            for num in party["candidates"]:
                group.cands.append(num)

            id2group[pid] = group

            pid += 1

        # process atl ballots
        bcntr = 0
        for atl in data["atl"]:
            n = int(atl["n"])

            groups = [id2group[int(g)] for g in atl["parties"]]

            prefs = []
            for g in groups:
                prefs.extend(g.cands)
                

            blt = Ballot(bcntr, n, prefs, atl=True)

            fcand = candidates[prefs[0]]
            fcand.ballots.append(bcntr)
            fcand.fp_votes += n

            total_votes += n

            ballots.append(blt)
            bcntr += 1

        # process btl ballots
        for btl in data["btl"]:
            n = int(btl["n"])

            # no need to use cid2num since candidate ids range from 0 to 
            # num_cands -1 already.
            prefs = [int(c) for c in btl["candidates"]]

            blt = Ballot(bcntr, n, prefs, atl=False)

            fcand = candidates[prefs[0]]
            fcand.ballots.append(bcntr)
            fcand.fp_votes += n

            total_votes += n
            ballots.append(blt)

            bcntr += 1
    
    return candidates,ballots,id2group,cid2num,total_votes
            

def read_ballots_stv(path):
    ballots = []
    candidates = []
    id2group = {}
    cid2num = {}

    total_votes = 0

    with open(path, "r") as cvr:
        lines = cvr.readlines()

        # Skip the first 3 lines, the fourth line
        # indicates the number of candidates
        ncands = int(lines[3].strip())

        # The next 'ncands' lines represent candidate 
        # details
        cntr = 0
        for i in range(4, 4+ncands):
            toks = lines[i].strip().split('\t')

            # toks = [Name, Group, Position in Group]
            cand = Candidate(cntr,cntr)
            cand.name = toks[0]
            cand.group_id = toks[1]
            cand.position = int(toks[2])

            cid2num[cntr] = cntr
            candidates.append(cand)
            cntr += 1

        # Get group info
        ngroups = int(lines[5+ncands].strip())

        for i in range(6+ncands, 6+ncands+ngroups):
            toks = lines[i].strip().split('\t')

            # toks = [Group ID, Group name]
            group = Group(toks[0])
            group.name = "" if len(toks) < 2 else toks[1]
    
            id2group[group.id] = group

        # Add candidates to their groups
        for cand in candidates:
            id2group[cand.group_id].cands.append(cand.num)

        # Continue until we get to RATLS (above the line entries)
        lcntr = 6+ncands+ngroups
        numratls = 0
        for i in range(6+ncands+ngroups, len(lines)):
            line = lines[i].strip()

            if line.startswith("RATLs"):
                # Next line details the number of RATLs
                numratls = int(lines[i+1].strip())
                
                lcntr = i+2
                break

        bcntr = 0

        assert(lcntr > 0)

        # Read above the line votes
        for i in range(lcntr, lcntr + numratls):
            toks = lines[i].strip().split()
               
            # Last element of toks is the number of votes with
            # the given ranking of groups. We translate the above
            # the line vote into the sequence of candidates that the
            # vote would move between.
            votes = int(toks[-1])
            prefs = []

            for gid in toks[:-1]:
                group = id2group[gid]
                for c in group.cands:
                    prefs.append(c)

            ballot = Ballot(bcntr, votes, prefs, atl=True)
            ballots.append(ballot)

            fpcand = candidates[prefs[0]]
            fpcand.ballots.append(bcntr)
            fpcand.fp_votes += votes

            total_votes += votes

            bcntr += 1

        # lcntr+numratls+1 is the line detailing the nubmer of BTL entries
        numbtls = int(lines[lcntr+numratls+1].strip())

        for i in range(lcntr+numratls+2,lcntr+numratls+2+numbtls):
            toks = lines[i].strip().split()

            votes = int(toks[-1])
            prefs = [int(c) for c in toks[0].split(',')]

            ballot = Ballot(bcntr, votes, prefs)
            ballots.append(ballot)
            
            fpcand = candidates[prefs[0]]
            fpcand.ballots.append(bcntr)
            fpcand.fp_votes += votes

            total_votes += votes

            bcntr += 1

    return candidates,ballots,id2group,cid2num,total_votes


# This function extracts code from audit_assertion_utils.py in the 
# SHANGRLA repository.
def sample_size_kaplan_kolgoromov(margin, prng, N, error_rate, rlimit, t=1/2, \
    g=0.1, upper_bound=1, quantile=0.5, reps=20):

    clean = 1.0/(2 - margin/upper_bound)
    one_vote_over = (1-0.5)/(2-margin/upper_bound) 

    samples = [0]*reps

    for i in range(reps):
        pop = clean*np.ones(N)
        inx = (prng.random(size=N) <= error_rate)  # randomly allocate errors
        pop[inx] = one_vote_over

        sample_total = 0
        mart = (pop[0]+g)/(t+g) if t > 0 else  1
        p = min(1.0/mart,1.0)
        j = 1

        while p > rlimit and j < N:
            mart *= (pop[j]+g)*(1-j/N)/(t+g - (1/N)*sample_total)
    
            if mart < 0:
                break
            else:
                sample_total += pop[j] + g

            
            p = min(1.0/mart,1.0) if mart > 1.0 else 1.0

            j += 1;

        if p <= rlimit:
            samples[i] = j
        else:
            return np.inf 

    return np.quantile(samples, quantile)


def subsupermajority_sample_size(threshold, tally, invalid, args):
    # supermajority assertion: party p1 achieved 
    # more than 'threshold' of the vote.
    share = 1.0/(2*threshold)
    
    amean = (tally*share + 0.5*invalid)/args.voters

    m = 2*amean - 1

    sample_size = estimate(args.seed, m, args.voters, args.erate, \
        args.rlimit, args.t, args.g, share, args.rfunc, args.reps)
    return sample_size, m, share


def cand_vs_cand_sample_size(tally1, tally2, valid_votes, args):
    # assorter for a ballot b yields 1 for a vote for cand 1,
    # 0 for a vote for cand 2, and 0.5 for all other votes
    other = args.voters - (tally1 + tally2)
    amean = (tally1 + 0.5*other)/args.voters

    m = 2*amean - 1

    # Estimate sample size via simulation
    sample_size = estimate(args.seed, m, args.voters, args.erate, \
        args.rlimit, args.t, args.g, 1, args.rfunc, args.reps)
    return sample_size, m, 1


def estimate(seed, m, total_voters, erate, rlimit, t, g, ub, rfunc, REPS):
    sample_size = np.inf
    if rfunc == "kaplan_kolmogorov":
        prng = np.random.RandomState(seed) 
        sample_size =  sample_size_kaplan_kolgoromov(m, prng, total_voters, \
            erate, rlimit, t=t, g=g, upper_bound=ub, quantile=0.5, reps=REPS)
    else:
        # Use kaplan martingale
        risk_fn=lambda x: TestNonnegMean.kaplan_martingale(x,N=total_voters)[0]
                
        sample_size =  TestNonnegMean.initial_sample_size(risk_fn, \
            total_voters, m, erate, alpha=rlimit, t=t, upper_bound=ub,\
            reps=REPS, bias_up=True, quantile=0.5, seed=seed)

    if sample_size is np.inf:
        return np.inf

    return math.ceil(sample_size)


def index_of(item, values):
    idx = 0
    for i in values:
        if i == item:
            return idx
        idx += 1

    return None 

def next_cand(prefs, excluded):
    for p in prefs:
        if p in excluded:
            continue

        return p

    return None


def simulate_stv(ballots, candidates, nseats, order_c, order_a, log):
    totvotes = 0
    print("First preference tallies: ", file=log)

    for cand in candidates:
        cand.sim_votes = 0
        cand.max_votes = 0
        cand.bweights = []
        cand.standing = 1
        cand.seat = -1
        cand.surplus = -1

        for bid in cand.ballots:
            cand.bweights.append((bid, 1))
            cand.sim_votes += ballots[bid].votes

        cand.max_votes = cand.sim_votes
        totvotes += cand.sim_votes

        print(f"    Candidate {cand.id} {cand.sim_votes}", file=log)

    # Step 1: Determine quota
    quota = (int)(1.0 + (totvotes/(nseats+1.0))) 

    print(f"The quota for election is {quota}", file=log)

    surpluses = []      

    currseat = 0

    while currseat < nseats:
        standing = 0

        # if a candidate has a quota, add them to the list of 
        # candidates with a surplus
        for cand in candidates:
            if cand.standing:
                standing += 1

            if cand.surplus != -1 or not cand.standing:
                continue

            if cand.standing and cand.sim_votes >= quota:
                cand.surplus = max(0, cand.sim_votes - quota)
                surpluses.append(cand)

        if standing == nseats - currseat:
            print("Number of candidates left standing equals number of "\
                "remaining seats", file=log)

            slist = []
            for cand in candidates:
                if cand.standing:
                    inserted = False
                    for i in range(len(slist)):
                        if cand.sim_votes > candidates[slist[i]].sim_votes:
                            slist.insert(i)
                            inserted = True
                            break

                    if not inserted:
                        standing_list.append(cand.num)

            for cnum in slist:
                cand = candidates[cnum]

                print("Candidate {} elected (votes {})".format(\
                    cand.name, cand.sim_votes), file=log)

                cand.seat = currseat
                currseat += 1

                cand.standing = 0
                order_c.append(cnum)
                order_a.append(1)

        if surpluses == []:
            # Eliminated candidate with fewest votes.
            # Distribute votes at their current value.
            leastvotes = -1
            toeliminate = -1

            for cand in candidates:
                if cand.standing:
                    print("Candidate {} has {} votes".format(cand.name,\
                        cand.sim_votes), file=log)
                    
                    if leastvotes == -1 or cand.sim_votes < leastvotes:
                        leastvotes = cand.sim_votes
                        toeliminate = cand

            
            order_c.append(toeliminate.num)
            order_a.append(0)

            print("Candidate {} eliminated on {} votes".format(\
                toeliminate.name, toeliminate.sim_votes), file=log)

            eliminate_candidate(toeliminate, candidates, ballots, log)

        else:
            # Start with candidate with the largest surplus
            elect = surpluses.pop(0)

            elect.seat = currseat
            currseat += 1

            order_c.append(elect.num)
            order_a.append(1)

            print("Candidate {} elected (votes {})".format(elect.name,\
                elect.sim_votes), file=log)

            if currseat < nseats:
                # Distribute surplus
                distribute_surplus(elect, candidates, ballots, log)


        if currseat == nseats:
            # All seats filled.
            if len(order_c) != len(candidates):
                for cand in candidates:
                    if cand.standing:
                        order_c.append(cand.num)
                        order_a.append(0)

            break


def next_candidate(prefs, cnum, candidates):
    idx = prefs.index(cnum)

    for p in prefs[idx+1:]:
        cand = candidates[p]
        if not cand.standing or cand.surplus != -1:
            continue

        return p

    return -1 


def distribute_surplus(elect, candidates, ballots, log):
    elect.standing = 0

    if elect.surplus < 0.001: return


    # Compute total number of papers in candidates tally
    totalpapers = sum([ballots[bid].votes for bid,_ in elect.bweights])

    tvalue = elect.surplus/totalpapers

    print("Transfer value is {}".format(tvalue), file=log)

    # Each ballot in elect's tally now has value of 'tvalue'
    totransfer = [[] for c in candidates]

    for bid,_ in elect.bweights:
        blt = ballots[bid]

        nextc = next_candidate(blt.prefs, elect.num, candidates)

        if nextc != -1:
            totransfer[nextc].append((bid, tvalue))

    for cand in candidates:
        tlist = totransfer[cand.num]

        total = 0
        for bid,weight in tlist:
            blt = ballots[bid]

            cand.ballots.append(bid)
            cand.sim_votes += weight*blt.votes

            total += weight*blt.votes

            cand.bweights.append((bid,weight))

        print("{} votes distributed from {} to {}".format(\
            total, elect.name, cand.name), file=log)

    elect.sim_votes -= elect.surplus
    elect.surplus = -1
        

def eliminate_candidate(toelim, candidates, ballots, log):
    toelim.standing = 0

    totransfer = [[] for c in candidates]

    # Distribute all ballots (at their current value) to rem candidates
    for bid,weight in toelim.bweights:
        nextc = next_candidate(ballots[bid].prefs, toelim.num, candidates)

        if nextc != -1:
            totransfer[nextc].append((bid, weight))

    toelim.sim_votes = 0

    for cand in candidates:
        tlist = totransfer[cand.num]

        total = 0
        for bid,weight in tlist:
            blt = ballots[bid]

            cand.ballots.append(bid)
            cand.sim_votes += weight*blt.votes

            total += weight*blt.votes

            cand.bweights.append((bid,weight))

        if total > 0:
            print("{} votes distributed from {} to {}".format(\
                total, toelim.name, cand.name), file=log)

