import argparse
import os
import numpy as np
import sys

from utils import read_ballots_json, read_ballots_blt, simulate_stv, \
    print_summary, get_stats_blt, get_stats_json

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    # Input: stv data file
    parser.add_argument('-d', dest='data')

    parser.add_argument('-stats', dest='stats', default=False, \
        action='store_true')

    args = parser.parse_args()

    candidates, ballots, id2group, cid2num, valid_ballots = \
        None, None, None, None, None

    prefix = None

    actual_seats = 0
    if args.data.endswith(".json"):
        candidates, ballots, id2group, cid2num, valid_ballots = \
            read_ballots_json(args.data)
        prefix = args.data[:-5]
        _, actual_seats = get_stats_json(args.data)

    elif args.data.endswith(".blt"):
        candidates, ballots, id2group, cid2num, valid_ballots = \
            read_ballots_blt(args.data)
        prefix = args.data[:-4]
        _, actual_seats = get_stats_blt(args.data)

    else:
        print("Unrecognised ballot file type.")
        exit()

    order_c = []
    order_a = []

    order_q = [-1]*len(candidates)
    winners = []

    quota = simulate_stv(ballots, candidates, actual_seats, order_c, order_a,\
        order_q, winners)

    if args.stats:
        # Candidates,Actual Seats,Case
        print("{},{},{}".format(len(candidates), actual_seats , order_a[0]))

    else:
        # Note: total number of voters != valid_ballots but for lack of data,
        # treat valid_ballots as representing the total number of voters for
        # the purposes of experiments.
        print_summary(candidates, id2group, actual_seats, quota, order_c, \
            order_a, order_q, winners, valid_ballots, prefix)
