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
        self.n1 = 0
        self.n2 = 0
        self.n3 = 0
        for i, x in enumerate(self.qts):
            print("x", x)
            print("sigmasq", self.sigmasq)
            print("q", self.q)
            print("mu1", self.mu1)
            print("mu2", self.mu2)
            print("mu3", self.mu3)
            #Hom Ref
            aalik = log(1.0 - self.q) + log(1.0 - self.q) + scipy.stats.norm.logpdf(x, self.mu1, self.sigma)
            #Het
            ablik = log(2 * self.q * (1.0 - self.q)) + scipy.stats.norm.logpdf(x, self.mu2, self.sigma)
            #Hom Alt
            bblik = log(self.q) + log(self.q) + scipy.stats.norm.logpdf(x, self.mu3, self.sigma)
            print(aalik, ablik, bblik)
            aalik_norm = math.exp(aalik)
            ablik_norm = math.exp(ablik)
            bblik_norm = math.exp(bblik)
            print(aalik_norm, ablik_norm, bblik_norm)
            sumlik_norm = self.mu1 * aalik_norm + self.mu2 * ablik_norm + self.mu3 * bblik_norm
            print("sumlik_norm", aalik_norm + ablik_norm + bblik_norm)
            if sumlik_norm == 0:
                aalik_norm = 1.0/3
                ablik_norm = 1.0/3
                bblik_norm = 1.0/3
            else:
                aalik_norm = (self.mu1 * aalik_norm)/sumlik_norm
                ablik_norm = (self.mu2 * ablik_norm)/sumlik_norm
                bblik_norm = (self.mu3 * bblik_norm)/sumlik_norm
            #Temporary
            aalik_norm = aalik
            ablik_norm = ablik
            bblik_norm = bblik

            print(aalik_norm, ablik_norm, bblik_norm)
            print(aalik_norm + ablik_norm + bblik_norm)

            self.gtsaa[i] = 0
            self.gtsab[i] = 0
            self.gtsbb[i] = 0
            if aalik_norm > ablik_norm and aalik_norm > bblik_norm:
                self.gtsaa[i] = 1.0
                self.n1 += 1
            if ablik_norm > aalik_norm and ablik_norm > bblik_norm:
                self.gtsab[i] = 1.0
                self.n2 += 1
            if bblik_norm > aalik_norm and bblik_norm > ablik_norm:
                self.gtsbb[i] = 1.0
                self.n3 += 1
            print("gt[i]", self.gtsaa[i], self.gtsab[i], self.gtsbb[i])
        print("gt", self.gtsaa, self.gtsab, self.gtsbb, len(self.gtsaa), len(self.gtsab))

    def m_step(self):
        n = len(self.qts)
        #Total sum for each genotype
        aasum = 0
        absum = 0
        bbsum = 0
        self.q = 0
        for i, x in enumerate(self.qts):
            aasum += x * self.gtsaa[i]
            absum += x * self.gtsab[i]
            bbsum += x * self.gtsbb[i]
            self.q += 2 * self.gtsbb[i] + self.gtsab[i]
        self.q = self.q/(2 * n) + sys.float_info.min
        if self.n1 == 0:
            self.mu1 = 0
        else:
            self.mu1 = aasum/(self.n1)
        if self.n2 == 0:
            self.mu2 = 0
        else:
            self.mu2 = absum/(self.n2)
        if self.n3 == 0:
            self.mu3 = 0
        else:
            self.mu3 = bbsum/(self.n3)
        self.sigmanew = 0
        for i, x in enumerate(self.qts):
            self.sigmanew += ((x * self.gtsaa[i]) ** 2 - (x * self.gtsaa[i]))
            self.sigmanew += ((x * self.gtsab[i]) ** 2 - (x * self.gtsab[i]))
            self.sigmanew += ((x * self.gtsbb[i]) ** 2 - (x * self.gtsbb[i]))
        print("sigmanew", self.sigmanew)
        self.sigma = math.sqrt(self.sigmanew/(n-1))
        self.sigma = 0.1
        self.q = 0.5
        print("sigma", self.sigma)
        print("qnew", self.q)
        print("gts", sum(self.gtsaa), sum(self.gtsab), sum(self.gtsbb))

    def run_em(self):
        #Initial condition
        self.mu1 = 10
        self.mu2 = 10
        self.mu3 = 10
        self.q = 0.5
        self.sigma = 0.1
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
