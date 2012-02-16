import pylab as pl
from RegionDef import *

def background_finder():
    from Bio import SeqIO
    grouped = {'alpha1': {}, 'alpha2': {}, 'beta1': {}, 'beta2': {}, 'beta3': {}, 'beta4': {}}

    for record in SeqIO.parse(open("superfam/aggregate_filtered.fasta", "r"), "fasta"):
        candidates = {}
        for r in regiondefs:
            #print(str(record.seq))
            candidates[r] = region_finder(r, str(record.seq))
            candidates[r].process()
        for b4 in candidates[beta4].prospects:
            try:
                grouped['beta4'][b4.score()] += 1
            except:
                grouped['beta4'][b4.score()] = 1
        for a2 in candidates[alpha2].prospects:
            try:
                grouped['alpha2'][a2.score()] += 1
            except:
                grouped['alpha2'][a2.score()] = 1
        for b3 in candidates[beta3].prospects:
            try:
                grouped['beta3'][b3.score()] += 1
            except:
                grouped['beta3'][b3.score()] = 1
        for b2 in candidates[beta2].prospects:
            try:
                grouped['beta2'][b2.score()] += 1
            except:
                grouped['beta2'][b2.score()] = 1
        for a1 in candidates[alpha1].prospects:
            try:
                grouped['alpha1'][a1.score()] += 1
            except:
                grouped['alpha1'][a1.score()] = 1
        for b1 in candidates[beta1].prospects:
            try:
                grouped['beta1'][b1.score()] += 1
            except:
                grouped['beta1'][b1.score()] = 1
    # the histogram of the data with histtype='step'
    n, bins, patches = pl.hist(grouped['alpha1'], 50, normed=1, histtype='stepfilled')
    pl.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
    pl.figure()
    pl.plot()

    for g in grouped.keys():
        grouped[g]['total'] = sum([grouped[g][k] for k in grouped[g].keys()])
            
    return grouped


