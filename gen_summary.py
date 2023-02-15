import argparse
import os
import numpy as np
import sys

from utils import read_ballots_json, simulate_stv, print_summary

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    # Input: stv data file
    parser.add_argument('-d', dest='data')

    # Input: number of seats in election (default is 2)
    parser.add_argument('-seats', dest='seats', default=2, type=int)

    args = parser.parse_args()

    candidates, ballots, id2group, cid2num, valid_ballots = \
        None, None, None, None, None

    candidates, ballots, id2group, cid2num, valid_ballots = \
        read_ballots_json(args.data)

    order_c = []
    order_a = []

    order_q = [-1]*len(candidates)
    winners = []

    quota = simulate_stv(ballots, candidates, args.seats, order_c, order_a,\
        order_q, winners) 

    print_summary(candidates, id2group, args.seats, quota, order_c, order_a,\
        order_q, winners)
