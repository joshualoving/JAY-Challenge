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
    def __str__(self):
        return str(self.score()) + self.seq + str(self.delim[0]) + "-" + str(self.delim[1])
class region_finder:
    def __init__(self, region, seq):
        self.region = region
        self.seq = seq
        self.prospects = []
        
    def process(self):
        for i in range(0, len(self.seq) - self.region.length):
            #print(self.seq[i:i + self.region.length])
            self.prospects.append(candidate(self.seq[i:i + self.region.length], self.region, (i,i + self.region.length)))
    def sample(self, start, stop=0):
        if stop == 0:
            stop = len(self.seq)
        from random import randint
        for i in range(len(self.prospects)):
            if self.prospects[i].delim[1] > start:
                i = i - 1
                break
        print len(self.prospects)
        if i == 0 or i == -1:
            i  = len(self.prospects) - 1
        print i
        return [self.prospects[randint(0, i)]for j in self.prospects]
    def top20(self):
        sorts = sorted([(p.score(), p) for p in self.prospects])
        sorts.reverse()
        sorts = [s[1] for s in sorts[:20]]
        return sorts
        
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
    def html(self):
        aln = self.seq
        counter = 0
        for reg in [self.beta1, self.alpha1, self.beta2, self.beta3, self.alpha2, self.beta4]:
            aln = aln[:reg.delim[0] + counter]+"<SPAN style=\"BACKGROUND-COLOR: #ffff00\">" + reg.seq  + "</SPAN>" + aln[reg.delim[1]+counter:]
            counter += 46
        return """Sequence: <br/>
                  MRKIIHVDMDCFFAAVEMRDNPALRDIPIAIGGSRERRGVISTANYPARKFGVRSAMPTGMALKLCPHLTLLPGRFDAYKEASNHIREIFSRYTSRIEPLSLDEAYLDVTDSVHCHGSATLIAQEIRQTIFNELQLTASAGVAPVKFLAKIASDMNKPNGQFVITPAEVPAFLQTLPLAKIPGVGKVSAAKLEAMGLRTCGDVQKCDLVMLLKRFGKFGRILWERSQGIDERDVNSERL<SPAN style="BACKGROUND-COLOR: #ffff00">RKSVGVERTM</SPAN>AEDIH<SPAN style="BACKGROUND-COLOR: #ffff00">HWSECEAIIERLYPELERRL</SPAN>AKVKPDLLI<SPAN style="BACKGROUND-COLOR: #ffff00">ARQGVKLKF</SPAN>DD<SPAN style="BACKGROUND-COLOR: #ffff00">FQQTTQEHVW</SPAN>PRL<SPAN style="BACKGROUND-COLOR: #ffff00">NKADLIATARKTWDE</SPAN>RRGGRG<SPAN style="BACKGROUND-COLOR: #ffff00">VRLVGLHVTL</SPAN>LDPQMERQLVLGL <p>
                  Aligned: <br/>
                  """ + aln 
beta1 = regiondef('beta1', 10, [3],[0, 1, 5], [2, 8])
alpha1 = regiondef('alpha1', 19, [0, 4, 7, 10, 11, 14, 18], [17], [])
beta2 = regiondef('beta2', 10, [2, 4, 6], [],[1])
beta3 = regiondef('beta3', 10, [],[5], [4,8])
alpha2 = regiondef('alpha2', 15,[4, 5, 7, 8, 11, 12], [], [14])
beta4 = regiondef('beta4', 11, [2, 4, 6], [0, 1, 5, 8], [])
betaclamp = regiondef('betaclamp', 5,[1,2,3,4], red=[0], pink=[]) 
regiondefs = [beta1, alpha1, beta2, beta3, alpha2, beta4, betaclamp]
