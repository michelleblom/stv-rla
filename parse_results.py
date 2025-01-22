import os
import sys
import numpy as np
import statistics

if __name__ == '__main__':

    max_samples = int(sys.argv[2])

    # map of seats : (num instances, num fully verified, partial map)
    # where the partial map is {num verified : (num instances, [ASNs])}
    # max sample size for full verified given by 'max_samples'
    summary = {}

    with open(sys.argv[1], 'r') as data:
        for line in data.readlines():
            toks = line.strip().split(',')

            seats = int(toks[1])
            verified = int(toks[9])
            samples = int(toks[10]) if toks[10] != "inf" else np.inf

            if not seats in summary:
                if samples > max_samples or verified == 0:
                    summary[seats] = (1, 0, {}, [])
                elif verified < seats:
                    summary[seats] = (1, 0, {verified : (1, [samples])})
                else:
                    summary[seats] = (1, 1, {}, [samples])

            else:
                insts, numfull, pmap, samples_full = summary[seats]

                if samples > max_samples or verified == 0:
                    summary[seats] = (insts+1, numfull, pmap, samples_full)

                elif verified < seats:
                    if not verified in pmap:
                        pmap[verified] = (1, [samples])
                    else:
                        num, spls = pmap[verified]
                        pmap[verified] = (num + 1, spls + [samples])

                    summary[seats] = (insts+1, numfull, pmap, samples_full)
                else:
                    summary[seats] = (insts+1, numfull+1, pmap, samples_full + [samples])


    for seats in summary:
        insts, numfull, pmap, samples_full = summary[seats]
        
        total = numfull + sum([n for (n,_) in [vm for (v, vm) in pmap.items()]])
 
        avg_full = statistics.mean(samples_full)
        min_full = min(samples_full)
        max_full = max(samples_full)

        print("{} instances with {} seats".format(insts, seats))
        print("    {}/{} fully verified within {} samples (avg {}, min {}, max {})".format(
            numfull, insts, max_samples, avg_full, min_full, max_full))

        for ps,vm in pmap.items():
            avg_p = statistics.mean(vm[1])
            min_p = min(vm[1])
            max_p = max(vm[1])
            print("    {} winners verified in {}/{} within {} samples  (avg {}, min {}, max {})".format(\
                ps, vm[0], insts, max_samples, avg_p, min_p, max_p))

        print("    {}/{} instances, no audit".format(insts-total, insts))
