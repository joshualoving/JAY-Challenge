from RegionDef import *
import pylab as pl
from BackgroundFinder import background_finder                               
from Bio import SeqIO
lfs = []

def filterer(stuff):
    s = stuff[0].score()
    if background[r.name][round(s)]/float(background[stuff[1].name]['total']) > 0.10:
        return False
    return True
number = input("Please input the number of little fingers to find")
background = background_finder()
# the histogram of the data with histtype='step'
outregs = open("top20.csv", 'w')
for record in SeqIO.parse(open("../../../../Dropbox/JAY/superfam/ecoli.txt", "r"), "fasta"):
    candidates = {}
    for r in regiondefs:
        #print(str(record.seq))
        candidates[r] = region_finder(r, str(record.seq))
        candidates[r].process()
        outregs.write(r.name)
        outregs.write(", \t")
        for c in candidates[r].top20():
            outregs.write(str(c))
            outregs.write(",\t")
        outregs.write("\n")
        #print(len(candidates[r].prospects))
        candidates[r].prospects = [e[0] for e in filter(filterer, [(c, r) for c in candidates[r].prospects])]
        #print(len(candidates[r].prospects))
        
    ordering = [beta4, alpha2, beta3, beta2, alpha1, beta1]
    traceback = [[] for i in range(number)]
    for path in traceback:
        path.extend([[c] for c in candidates[betaclamp].sample(len(str(record.seq)))])
    
    for reg in ordering:
        for path in traceback:
            for p in path:
                #print str(p[-1].delim[0])
                
                p.append(max([(bc.score() - 0.05*(p[-1].delim[0] - bc.delim[1]), bc) for bc in candidates[reg].sample(p[-1].delim[0]-2)])[1])
    for path in traceback:
        score, select = max([(sum([r.score() for r in p]), p) for p in path])
        print select[6].seq, select[5].seq, select[4].seq, select[3].seq, select[2].seq, select[1].seq, select[0].seq
        lfs.append(lf(str(record.seq), select[6], select[5], select[4], select[3], select[2], select[1], select[0], score))
outhtml = open("lf.html", "w")
outhtml.write("<html><body>")
for l in lfs:
    print str(l)
    outhtml.write(l.html())
    outhtml.write("<p>")
outhtml.write("</html></body>")
