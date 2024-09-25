import os
import sys
import numpy as np

if __name__ == '__main__':

    max_samples = int(sys.argv[2])

    # map of seats : (num instances, num fully verified, partial map)
    # where the partial map is {num verified : num instances}
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
                    summary[seats] = (1, 0, {})
                elif verified < seats:
                    summary[seats] = (1, 0, {verified : 1})
                else:
                    summary[seats] = (1, 1, {})

            else:
                insts, numfull, pmap = summary[seats]

                if samples > max_samples or verified == 0:
                    summary[seats] = (insts+1, numfull, pmap)

                elif verified < seats:
                    if not verified in pmap:
                        pmap[verified] = 1
                    else:
                        pmap[verified] = pmap[verified] + 1

                    summary[seats] = (insts+1, numfull, pmap)
                else:
                    summary[seats] = (insts+1, numfull+1, pmap)


    for seats in summary:
        insts, numfull, pmap = summary[seats]
        
        total = numfull + sum([n for (v,n) in pmap.items()])
 
        print("{} instances with {} seats".format(insts, seats))
        print("    {}/{} fully verified within {} samples".format(numfull, \
            insts, max_samples))

        for ps,n in pmap.items():
            print("    {} winners verified in {}/{} within {} samples".format(\
                ps, n, insts, max_samples))

        print("    {}/{} instances, no audit".format(insts-total, insts))
