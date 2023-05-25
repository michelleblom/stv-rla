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


import os
import numpy as np
import statistics
import math
import re
import json

# Make sure shangrla is in your PYTHONPATH
from shangrla.NonnegMean import NonnegMean


class Ballot:
    def __init__(self, num, votes, prefs):
        self.num = num
        self.votes = votes
        self.prefs = prefs[:]

        self.papers = votes

class Candidate:
    def __init__(self, num, idn):
        self.num = num
        self.id = idn
        self.name = None
        self.group_id = None
        self.position = None

        self.reset()

    def reset(self):
        self.mentions = []
        self.ballots = []
        self.fp_votes = 0

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


def read_ballots_blt(path):
    ballots = []
    candidates = []
    cid2num = {}

    total_votes = 0

    with open(path, "r") as cvr:
        lines = [l.strip() for l in cvr.readlines() if l.strip() != ""]

        cands,seats = [int(tok) for tok in lines[0].split()]

        split_idx = lines.index("0")

        ballot_strs = lines[:split_idx]
        cand_strs = lines[split_idx+1:-1]

        for i in range(len(cand_strs)):
            cand = Candidate(i, cand_strs[i])
            cand.name = str(clist[i])
            cand.group_id = -1
            cand.position = -1

            candidates.append(cand)
            cid2num[clist[i]] = i


        bcntr = 0

        for bline in ballot_strs:
            toks = [int(t) for t in bline.split()]

            n = toks[0]
            prefs = toks[1:-1]

            ballot = Ballot(bcntr, n, prefs)
            ballots.append(ballot)

            fpcand = candidates[cprefs[0]]
            fpcand.ballots.append(bcntr)
            fpcand.fp_votes += votes

            total_votes += votes

            for p in prefs:
                candidates[p].mentions.append(bcntr)

            bcntr += 1

    return candidates,ballots,{},cid2num,total_votes


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

            ballot = Ballot(bcntr, votes, prefs)
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
            cand.name = str(clist[i])
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
            ballot = Ballot(bcntr, votes, cprefs)
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

    prefs2ballot = {}

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
                
            blt = None
            fcand = candidates[prefs[0]]
            tprefs = tuple(prefs)
            if tprefs in prefs2ballot:
                blt = prefs2ballot[tprefs]
                blt.votes += n
            else:
                blt = Ballot(bcntr, n, prefs)
                prefs2ballot[tprefs] = blt
                ballots.append(blt)
                fcand.ballots.append(blt.num)
                bcntr += 1

            fcand.fp_votes += n

            total_votes += n

        # process btl ballots
        for btl in data["btl"]:
            n = int(btl["n"])

            # no need to use cid2num since candidate ids range from 0 to 
            # num_cands -1 already.
            prefs = [int(c) for c in btl["candidates"]]

            blt = None
            tprefs = tuple(prefs)
            fcand = candidates[prefs[0]]
            if tprefs in prefs2ballot:
                blt = prefs2ballot[tprefs]
                blt.votes += n
            else:
                blt = Ballot(bcntr, n, prefs)
                prefs2ballot[tprefs] = blt
                ballots.append(blt)
                fcand.ballots.append(blt.num)
                bcntr += 1

            fcand.fp_votes += n

            total_votes += n

    return candidates,ballots,id2group,cid2num,total_votes
            


def sample_size(mean, args, upper_bound=1): # just for comparison audits
    N = args.voters
    margin = 2*mean - 1
    u = 2/(2-(margin/upper_bound))

    test = NonnegMean(test=NonnegMean.alpha_mart, \
        estim=NonnegMean.optimal_comparison, N=N, u=u, eta=mean)

    # over: (1 - o/u)/(2 - v/u)

    # where o is the overstatement, u is the upper bound on the value
    # assorter assigns to any ballot, v is the assorter margin.
    big =  (1 / (2 - margin/upper_bound)) # o=0
    small = (0.5 / (2 - margin/upper_bound)) # o=0.5

    r1 = args.erate1
    r2 = args.erate2

    x = big*np.ones(N)

    rate_1_i = np.arange(0, N, step=int(1/r1), dtype=int) if r1 else []
    rate_2_i = np.arange(0, N, step=int(1/r2), dtype=int) if r2 else []

    x[rate_1_i] = small
    x[rate_2_i] = 0

    return test.sample_size(x, alpha=args.rlimit, reps=args.reps, \
        seed=args.seed, random_order=True)


def ssm_sample_size(threshold, tally, invalid, args):
    # supermajority assertion: party p1 achieved 
    # more than 'threshold' of the vote. 
    share = 1.0/(2*threshold)
    
    amean = (tally*share + 0.5*invalid)/args.voters

    sam_size = sample_size(amean, args, upper_bound = share)
    return sam_size


def tally_vs_tally_sample_size(tally1, tally2, valid_votes, args):
    # assorter for a ballot b yields 1 for a vote for cand 1,
    # 0 for a vote for cand 2, and 0.5 for all other votes
    other = args.voters - (tally1 + tally2)

    amean = (tally1 + 0.5*other)/args.voters

    # Estimate sample size via simulation
    sam_size = sample_size(amean, args)
    return sam_size 



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


def simulate_stv(ballots, candidates, nseats, order_c, order_a, order_q, \
    winners, log=None):
    
    totvotes = 0
    if log != None:
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

        if log != None:
            print(f"    Candidate {cand.id} {cand.sim_votes}", file=log)

    # Step 1: Determine quota
    quota = (int)(1.0 + (totvotes/(nseats+1.0))) 

    if log != None:
        print(f"The quota for election is {quota}", file=log)

    surpluses = []      

    currseat = 0

    r = 0

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
                insert_surplus(surpluses, cand)

                order_q[cand.num] = r

        if standing == 0:
            break

        if standing == nseats - currseat:
            if log != None:
                print("Number of candidates left standing equals number of "\
                    "remaining seats", file=log)

            slist = []
            for cand in candidates:
                if cand.standing:
                    inserted = False
                    for i in range(len(slist)):
                        if cand.sim_votes > candidates[slist[i]].sim_votes:
                            slist.insert(i, cand.num)
                            inserted = True
                            break

                    if not inserted:
                        slist.append(cand.num)

            for cnum in slist:
                cand = candidates[cnum]

                if log != None:
                    print("Candidate {} elected (votes {})".format(\
                        cand.name, cand.sim_votes), file=log)

                cand.seat = currseat
                currseat += 1

                cand.standing = 0
                order_c.append(cnum)
                order_a.append(1)

                winners.append(cand.num)

        elif surpluses == []:
            # Eliminated candidate with fewest votes.
            # Distribute votes at their current value.
            leastvotes = -1
            toeliminate = -1

            for cand in candidates:
                if cand.standing:
                    if log != None:
                        print("Candidate {} has {} votes".format(cand.name,\
                            cand.sim_votes), file=log)
                    
                    if leastvotes == -1 or cand.sim_votes < leastvotes:
                        leastvotes = cand.sim_votes
                        toeliminate = cand

           
            if toeliminate != -1: 
                order_c.append(toeliminate.num)
                order_a.append(0)

                if log != None:
                    print("Candidate {} eliminated on {} votes".format(\
                        toeliminate.name, toeliminate.sim_votes), file=log)

                eliminate_candidate(toeliminate, candidates, ballots, log)

            r += 1

        else:
            new_surpluses = []

            while surpluses != []:
                # Start with candidate with the largest surplus
                elect = surpluses.pop(0)

                elect.seat = currseat
                currseat += 1

                order_c.append(elect.num)
                order_a.append(1)

                winners.append(elect.num)

                if log != None:
                    print("Candidate {} elected (votes {})".format(elect.name,\
                        elect.sim_votes), file=log)

                if currseat < nseats:
                    # Distribute surplus
                    distribute_surplus(elect, candidates, ballots, log)

                next_surpluses = []
                for cand in candidates:
                    if cand.surplus != -1 or not cand.standing:
                        continue

                    if cand.sim_votes >= quota:
                        cand.surplus = max(0, cand.sim_votes - quota)
                        insert_surplus(next_surpluses, cand)
                        order_q[cand.num] = r

                new_surpluses.extend(next_surpluses)

                r += 1

            surpluses = new_surpluses

        if currseat == nseats:
            # All seats filled.
            if len(order_c) != len(candidates):
                for cand in candidates:
                    if cand.standing:
                        order_c.append(cand.num)
                        order_a.append(0)

            break

    return quota

def insert_surplus(surpluses, cand):
    for i in range(len(surpluses)):
        if cand.surplus >= surpluses[i].surplus:
            surpluses.insert(i, cand)
            return

    surpluses.append(cand)


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

    if log != None:
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

        if log != None:
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
            if log != None:
                print("{} votes distributed from {} to {}".format(\
                    total, toelim.name, cand.name), file=log)


def print_summary(candidates,id2group, seats, quota, order_c, order_a,\
    order_q, winners):
    print(f"Candidates,{len(candidates)},Groups,{len(id2group)}")
    print(f"Seats,{seats}")
    print(f"Quota,{quota}")
    
    order_c_ids = [str(candidates[c].id) for c in order_c]

    order_c_str = "Outcome-ns"
    order_a_str = "Outcome-ns"

    for i in range(len(candidates)):
        order_c_str += "," + order_c_ids[i]
        order_a_str += "," + str(order_a[i])

    
    print(order_c_str)
    print(order_a_str)

    for w in winners:
        print("Quota,{},{}".format(w, order_q[w]))

    for i in range(len(id2group)):
        gstr = "Group,{},Candidates".format(i)

        group = id2group[i] 
        for cnum in group.cands:
            gstr += "," + str(candidates[cnum].id)

        print(gstr)

    for cand in candidates:
        islast = 1 if id2group[cand.group_id].cands[-1] == cand.num else 0

        print("Candidates,{},{},{},{},{},{},({})".format(cand.id,cand.num_atls,\
            cand.num_btls,cand.group_id,cand.position-1,islast,cand.name))
    
