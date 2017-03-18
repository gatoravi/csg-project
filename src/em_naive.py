import sys
import csv
import math
import scipy.stats
from math import log as log

class NaiveEM:
    #This is the latent variable - genotype for every x_i
    gts = []
    gtsaa = []
    gtsab = []
    gtsbb = []
    #List of x_i s
    qts = []

    def __init__(self, qt_file):
        with open(qt_file) as qtfh:
            qtfh_dict = csv.DictReader(qtfh, delimiter = ",")
            for line in qtfh_dict:
                self.qts.append(float(line['QTp']))
        for i in range(len(self.qts)):
            self.gts.append(-1) #Initialize with null GT
            self.gtsaa.append(-1) #Initialize with null GT
            self.gtsab.append(-1) #Initialize with null GT
            self.gtsbb.append(-1) #Initialize with null GT

    def e_step(self):
        "In the E step, figure out the latent variable, i.e genotype for every x"
        self.sigmasq = self.sigma ** 2
        if self.q == 0:
            self.q = 0.01
        if self.q == 1.0:
            self.q = 0.99
        print("q", self.q)
        for i, x in enumerate(self.qts):
            print("x", x)
            print("sigmasq", self.sigmasq)
            print("q", self.q)
            print("mu1", self.mu1)
            print("mu2", self.mu2)
            print("mu3", self.mu3)
            #Hom Ref
            aalik = log((1.0 - self.q) * (1.0 - self.q)) + scipy.stats.norm.logpdf(x, self.mu1, self.sigma)
            #Het
            ablik = log(2 * self.q * (1.0 - self.q)) + scipy.stats.norm.logpdf(x, self.mu2, self.sigma)
            #Hom Alt
            bblik = log(self.q) + log(self.q) + scipy.stats.norm.logpdf(x, self.mu3, self.sigma)
            print(aalik, ablik, bblik)
            aalik_norm = math.exp(aalik)
            ablik_norm = math.exp(ablik)
            bblik_norm = math.exp(bblik)
            print(aalik_norm, ablik_norm, bblik_norm)
            aalik_norm = aalik_norm/(aalik_norm + ablik_norm + bblik_norm + sys.float_info.min)
            ablik_norm = ablik_norm/(aalik_norm + ablik_norm + bblik_norm + sys.float_info.min)
            bblik_norm = bblik_norm/(aalik_norm + ablik_norm + bblik_norm + sys.float_info.min)
            print(aalik_norm, ablik_norm, bblik_norm)
            print(aalik_norm + ablik_norm + bblik_norm)
            print("un-normalized", aalik, ablik, bblik)
            #Get max lik genotype
            if aalik >= ablik and aalik >= bblik:
                self.gts[i] = 0
            elif ablik >= aalik and ablik >= bblik:
                self.gts[i] = 1
            elif bblik >= aalik and bblik >= ablik:
                self.gts[i] = 2
            self.gtsaa[i] = aalik_norm
            self.gtsab[i] = ablik_norm
            self.gtsbb[i] = bblik_norm
            print("gt[i]", self.gts[i], self.gtsaa[i], self.gtsab[i], self.gtsbb[i])
        print("gt", self.gts, len(self.gts))

    def m_step(self):
        n = len(self.qts)
        n1 = 0 #Number of AA genotypes
        n2 = 0 #Number of AB genotypes
        n3 = 0 #Number of BB
        #Total sum for each genotype
        aasum = 0
        absum = 0
        bbsum = 0
        for i, x in enumerate(self.qts):
            if self.gts[i] == 0:
                aasum += x
                n1 += 1
            elif self.gts[i] == 1:
                absum += x
                n2 += 1
            elif self.gts[i] == 2:
                bbsum += x
                n3 += 1
            else:
                raise RuntimeError("Invalid genotype!")
        print("n, n1, n2, n3", n, n1, n2, n3)
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
        self.sigmanew = 0
        for i, x in enumerate(self.qts):
            if self.gts[i] == 0:
                mean = self.mu1
            elif self.gts[i] == 1:
                mean = self.mu2
            elif self.gts[i] == 2:
                mean = self.mu3
            self.sigmanew += (mean - x) ** 2
        self.sigma = math.sqrt(self.sigmanew/(n-1))
        self.q = (n2 + 2 * n3) / (2 * n)
        print("qnew", self.q)

    def run_em(self):
        #Initial condition
        self.mu1 = 1
        self.mu2 = 2
        self.mu3 = 3
        self.q = 0.1
        self.sigma = 1
        for i in range(10000):
            print("Running E step")
            self.e_step()
            print("Running M step")
            self.m_step()
            print([self.mu1, self.mu2, self.mu3, self.sigma, self.q])

def main():
    qt_file = sys.argv[1]
    em1 = NaiveEM(qt_file)
    em1.run_em()

main()
