#Module of region stuff
class regiondef:
    hydro = {'A': 1, 'C': 1, 'F': 1, 'I': 1, 'L': 1, 'M': 1, 'P': 1, 'V': 1, 'W': 1, 'Y': .75, 'D': 0, 'E': 0, 'K': 0, 'N': 0, 'Q': 0, 'R': 0, 'S': 0, 'T': .25, 'G': 0.5, 'H': 0.5}
    red = {'A': 0, 'C': 0, 'F': 0, 'I': 0, 'L': 0, 'M': 0, 'P': 0, 'V': 0, 'W': 0, 'Y': 0, 'D': 0, 'E': 0, 'K': 1, 'N': 0.5, 'Q': 1, 'R': 1, 'S': 0, 'T': 0, 'G': 0, 'H': 0}
    pink = {'A': 0, 'C': 0, 'F': 0, 'I': 0, 'L': 0, 'M': 0, 'P': 0, 'V': 0, 'W': 0, 'Y': 0, 'D': 0, 'E': 0, 'Q': 0.5, 'K': 0.5, 'N': 1, 'R': 0.5, 'S': 1, 'T': 1, 'G': 0, 'H': 0} 
    def __init__(self, name, length, phobics, red = None, pink = None):
        self.name = name
        self.length = length
        self.phobics = phobics
        self.red = red
        self.pink = pink
        
class candidate:
    def __init__(self, seq, region, delims):
        self.seq = seq
        self.delim = delims
        self.region = region
    def score(self):
        #print (self.seq)
        phscore = sum([regiondef.hydro[self.seq[p]] for p in self.region.phobics])
        rscore = sum([regiondef.red[self.seq[r]] for r in self.region.red])
        pscore = sum([regiondef.pink[self.seq[p]] for p in self.region.pink])
        return phscore + rscore + pscore

class region_finder:
    def __init__(self, region, seq):
        self.region = region
        self.seq = seq
        self.prospects = []
        
    def process(self):
        for i in range(200, len(self.seq) - self.region.length):
            #print(self.seq[i:i + self.region.length])
            self.prospects.append(candidate(self.seq[i:i + self.region.length], self.region, (i,i + self.region.length)))

class lf:
    def __init__(self, seq, beta1, alpha1, beta2, beta3, alpha2, beta4, betaclamp, score):
        self.seq = seq
        self.beta1 = beta1
        self.alpha1 = alpha1
        self.beta2 = beta2
        self.beta3 = beta3
        self.alpha2 = alpha2
        self.beta4 = beta4
        self.betaclamp = betaclamp
        self.score = score
    def __str__(self):
        aln = len(self.seq)*'-'
        for reg in [self.beta1, self.beta2, self.beta3, self.beta4, self.alpha1, self.alpha2, self.betaclamp]:
            aln = aln[:reg.delim[0]] + reg.seq + aln[reg.delim[1]:]
        return "Score: " + str(self.score) + "\n" + self.seq + "\n" + aln
 
beta1 = regiondef('beta1', 10, [3],[0, 1, 5], [2, 8])
alpha1 = regiondef('alpha1', 19, [0, 4, 7, 10, 11, 14, 18], [17], [])
beta2 = regiondef('beta2', 10, [2, 4, 6], [],[1])
beta3 = regiondef('beta3', 10, [],[5], [4,8])
alpha2 = regiondef('alpha2', 15,[4, 5, 7, 8, 11, 12], [], [14])
beta4 = regiondef('beta4', 11, [2, 4, 6], [0, 1, 5, 8], [])
betaclamp = regiondef('betaclamp', 5,[1,2,3,4], red=[0], pink=[]) 
regiondefs = [beta1, alpha1, beta2, beta3, alpha2, beta4, betaclamp]
