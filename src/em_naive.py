import sys
import csv
import math
from math import log as log

class NaiveEM:
    #This is the latent variable - genotype for every x_i
    gts = []
    #List of x_i s
    qts = []

    def __init__(self, qt_file):
        with open(qt_file) as qtfh:
            qtfh_dict = csv.DictReader(qtfh, delimiter = ",")
            for line in qtfh_dict:
                self.qts.append(float(line['QTp']))

    def e_step(self):
        "In the E step, figure out the latent variable, i.e genotype for every x"
        self.sigmasq = self.sigma ** 2
        print("sigmasq", self.sigmasq)
        print("q", self.q)
        for i, x in enumerate(self.qts):
            #Get log lik of AA
            if self.q == 1.0: #No A alleles
                aalik = 0
            else:
                aalik = math.exp(2 * log(1.0 - self.q) - (x - self.mu1)**2/(2 * self.sigmasq) - 0.5 * log(2 * math.pi * self.sigmasq))
            #Get log lik of AB
            if self.q == 1.0 or self.q == 0: #No A allele [or] No B allele
                ablik = 0
            else:
                ablik = math.exp(log(2) + log(1.0 - self.q) + log(self.q) - (x - self.mu2)**2/(2 * self.sigmasq) - 0.5 * log(2 * math.pi * self.sigmasq))
            #Get log lik of BB
            if self.q == 0.0: #No B allele
                bblik = 0
            else:
                bblik = math.exp(2 * log(self.q) - (x - self.mu3)**2/(2 * self.sigmasq) - 0.5 * log(2 * math.pi * self.sigmasq))
            print(aalik, ablik, bblik)
            #Normalize
            aalik_norm = aalik/(aalik + ablik + bblik)
            ablik_norm = ablik/(aalik + ablik + bblik)
            bblik_norm = bblik/(aalik + ablik + bblik)
            print("normalized", aalik_norm, ablik_norm, bblik_norm)
            self.gts.insert(i, 0) #Set to AA by default
            #Get max lik genotype
            if aalik_norm >= ablik_norm and aalik_norm >= bblik_norm:
                self.gts[i] = 0
            if ablik_norm >= aalik_norm and ablik_norm >= bblik_norm:
                self.gts[i] = 1
            elif bblik_norm >= aalik_norm and bblik_norm >= ablik_norm:
                self.gts[i] = 2
            print("gt", self.gts[i])

    def m_step(self):
        n1 = 0 #Number of AA genotypes
        n2 = 0 #Number of AB genotypes
        n3 = 0 #Number of BB
        #Total sum for each genotype
        aasum = 0
        absum = 0
        bbsum = 0
        totalsum = 0
        totalmean = 0
        n = len(self.qts)
        for i, x in enumerate(self.qts):
            if self.gts[i] == 0:
                aasum += x
                totalsum +=x
                n1 += 1
            elif self.gts[i] == 1:
                absum += x
                totalsum +=x
                n2 += 1
            elif self.gts[i] == 2:
                bbsum += x
                totalsum +=x
                n3 += 1
            else:
                raise RuntimeError("Invalid genotype!")
        print("n, n1, n2, n3", n, n1, n2, n3)
        totalmean = totalsum/n
        self.sigmanew = 0
        for i, x in enumerate(self.qts):
            self.sigmanew += (totalmean - x) ** 2
        self.sigmanew = self.sigmanew/(n-1)
        if n1 == 0:
            self.mu1 = 0
        else:
            self.mu1 = aasum/n1
        if n2 == 0:
            self.mu2 = 0
        else:
            self.mu2 = absum/n2
        if n3 == 0:
            self.mu3 = 0
        else:
            self.mu3 = bbsum/n3
        self.sigma = self.sigmanew
        self.q = (n2 + 2 * n3) / (2 * n)

    def run_em(self):
        #Initial condition
        self.mu1 = -0.5
        self.mu2 =  0.0
        self.mu3 =  0.5
        self.q = 0.5
        self.sigma = 1
        for i in range(1000):
            self.e_step()
            self.m_step()
            print([self.mu1, self.mu2, self.mu3, self.sigma, self.q])

def main():
    qt_file = sys.argv[1]
    em1 = NaiveEM(qt_file)
    em1.run_em()

main()
