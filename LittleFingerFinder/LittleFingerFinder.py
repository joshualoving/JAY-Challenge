from RegionDef import *

from BackgroundFinder import background_finder                               
from Bio import SeqIO
lfs = []

def filterer(stuff):
    s = stuff[0].score()
    if background[r.name][round(s)]/float(background[stuff[1].name]['total']) > 0.05:
        return False
    return True
background = background_finder()
for record in SeqIO.parse(open("superfam/ecoli.txt", "r"), "fasta"):
    candidates = {}
    for r in regiondefs:
        #print(str(record.seq))
        candidates[r] = region_finder(r, str(record.seq))
        candidates[r].process()
        print(len(candidates[r].prospects))
        candidates[r].prospects = [e[0] for e in filter(filterer, [(c, r) for c in candidates[r].prospects])]
        print(len(candidates[r].prospects))
    for b4 in candidates[beta4].prospects:
        for a2 in candidates[alpha2].prospects:
            if a2.delim[1] > b4.delim[0]:# or (b4.delim[0] - a2.delim[1]) > 30:
                break
            for b3 in candidates[beta3].prospects:
                if b3.delim[1] > a2.delim[0]:# or (a2.delim[0] - b3.delim[1]) > 30:
                    break
                for b2 in candidates[beta2].prospects:
                    if b2.delim[1] > b3.delim[0]:# or (b3.delim[0] - b2.delim[1]) > 30:
                        break
                    for a1 in candidates[alpha1].prospects:
                        if a1.delim[1] > b2.delim[0]:# or (b2.delim[0] - a1.delim[1]) > 30:
                            break
                        for b1 in candidates[beta1].prospects:
                            if b1.delim[1] > a1.delim[0]:# or (a1.delim[0] - b1.delim[1]) > 30:
                               break
                            lfs.append(lf(candidates[r].seq, b1, a1, b2, b3, a2, b4, b1.score() + a1.score() + b2.score() + b3.score() + a2.score() + b4.score()))

print(len(lfs))
cool = sorted([(l.score, l) for l in lfs])
cool.reverse()

print(str(cool[0][1]))
print(str(cool[1][1]))