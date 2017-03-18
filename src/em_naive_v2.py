import sys
import csv
import math
import random
import scipy.stats
import statistics
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
        self.total_like_new = 0
        self.total_like_old = 0
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
        print("q", self.q)
        self.n1 = 0
        self.n2 = 0
        self.n3 = 0
        self.total_like_old = self.total_like_new
        self.total_like_new = 0
        for i, x in enumerate(self.qts):
            """
            print("x", x)
            print("sigmasq", self.sigmasq)
            print("q", self.q)
            print("1-q", 1 - self.q)
            print("mu1", self.mu1)
            print("mu2", self.mu2)
            print("mu3", self.mu3)
            print("fractions: ", log(1.0 - self.q) + log(1.0 - self.q), log(2) + log(self.q) + log(1.0 - self.q), log(self.q) + log(self.q))
            """
            #Hom Ref
            aalik = log(1.0 - self.q) + log(1.0 - self.q) + scipy.stats.norm.logpdf(x, self.mu1, self.sigma)
            #Het
            ablik = log(2) + log(self.q) + log(1.0 - self.q) + scipy.stats.norm.logpdf(x, self.mu2, self.sigma)
            #Hom Alt
            bblik = log(self.q) + log(self.q) + scipy.stats.norm.logpdf(x, self.mu3, self.sigma)
            aalik_new = aalik - max(aalik, ablik, bblik)
            ablik_new = ablik - max(aalik, ablik, bblik)
            bblik_new = bblik - max(aalik, ablik, bblik)
            aalik = aalik_new
            ablik = ablik_new
            bblik = bblik_new
            if aalik < -100:
                aalik = -100
            if ablik < -100:
                ablik = -100
            if bblik < -100:
                bblik = -100
            #print(scipy.stats.norm.logpdf(x, self.mu1, self.sigma))
            aalik_norm = math.exp(aalik)
            ablik_norm = math.exp(ablik)
            bblik_norm = math.exp(bblik)
            #print("after exp", aalik_norm, ablik_norm, bblik_norm)
            sumlik_norm = aalik_norm + ablik_norm + bblik_norm
            #print("sumlik_norm", aalik_norm + ablik_norm + bblik_norm)
            self.gtsaa[i] = math.exp(math.log(aalik_norm) - math.log(sumlik_norm))
            self.gtsab[i] = math.exp(math.log(ablik_norm) - math.log(sumlik_norm))
            self.gtsbb[i] = math.exp(math.log(bblik_norm) - math.log(sumlik_norm))
            self.n1 += self.gtsaa[i]
            self.n2 += self.gtsab[i]
            self.n3 += self.gtsbb[i]
            self.total_like_new += math.log(math.exp(aalik) + math.exp(ablik) + math.exp(bblik))
            """
            Ignore for now
            if sumlik_norm == 0:
                print("Inside norm")
                aalik_norm = 1.0/3
                ablik_norm = 1.0/3
                bblik_norm = 1.0/3
            else:
                aalik_norm = math.log(aalik_norm) - math.log(sumlik_norm)
                ablik_norm = math.log(ablik_norm) - math.log(sumlik_norm)
                bblik_norm = math.log(bblik_norm) - math.log(sumlik_norm)
            self.gtsaa[i] = 0
            self.gtsab[i] = 0
            self.gtsbb[i] = 0
            if aalik_norm >= ablik_norm and aalik_norm >= bblik_norm:
                self.gtsaa[i] = 1.0
                self.n1 += 1
            elif ablik_norm >= aalik_norm and ablik_norm >= bblik_norm:
                self.gtsab[i] = 1.0
                self.n2 += 1
            elif bblik_norm >= aalik_norm and bblik_norm >= ablik_norm:
                self.gtsbb[i] = 1.0
                self.n3 += 1
            else:
                raise("Invalid genotypes " + str(aalik_norm) + "," + str(ablik_norm) + "," + str(bblik_norm))
            """
            #print("gt[i]", self.gtsaa[i], self.gtsab[i], self.gtsbb[i])
        #print("gt", self.gtsaa, self.gtsab, self.gtsbb)
        #print("gtsum", [x + y + z for x, y, z in zip(self.gtsaa, self.gtsab, self.gtsbb)])
        print("liks", self.total_like_old, self.total_like_new, self.total_like_new - self.total_like_old)

    def m_step(self):
        n = len(self.qts)
        n = (self.n1 + self.n2 + self.n3)
        #Total sum for each genotype
        aasum = 0
        absum = 0
        bbsum = 0
        self.q = 0
        for i, x in enumerate(self.qts):
            #print("i", i, "x", x, "gts", self.gtsaa[i], self.gtsab[i], self.gtsbb[i])
            aasum += x * self.gtsaa[i]
            absum += x * self.gtsab[i]
            bbsum += x * self.gtsbb[i]
            gtsum = self.gtsaa[i] + self.gtsab[i] + self.gtsbb[i]
            if gtsum < 0.95 or gtsum > 1.05:
                print("gtsum", gtsum)
            self.q += 2 * self.gtsbb[i] + self.gtsab[i]
        print("aasum", "absub", "bbsum", aasum, absum, bbsum, "n1, n2, n3", self.n1, self.n2, self.n3)
        self.q = self.q/(2 * n)
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
            self.sigmanew += ((x - self.mu1)) ** 2 * self.gtsaa[i]
            self.sigmanew += ((x - self.mu2)) ** 2 * self.gtsab[i]
            self.sigmanew += ((x - self.mu3)) ** 2 * self.gtsbb[i]
        self.sigma = math.sqrt(self.sigmanew/n)
        print("sigma", self.sigma)
        print("qnew", self.q)
        print("gts", sum(self.gtsaa), sum(self.gtsab), sum(self.gtsbb))

    def run_em(self):
        #Initial condition
        self.mu1 = self.qts[random.randint(0, len(self.qts) - 1)]
        self.mu2 = self.qts[random.randint(0, len(self.qts) - 1)]
        self.mu3 = self.qts[random.randint(0, len(self.qts) - 1)]
        self.sigma = 10
        self.q = 0.01
        for i in range(1000):
            print("params mus, sigma, q", [self.mu1, self.mu2, self.mu3, self.sigma, self.q])
            print("Iteration: ", i+1)
            print("Running E step")
            self.e_step()
            print("Running M step")
            self.m_step()

def main():
    qt_file = sys.argv[1]
    em1 = NaiveEM(qt_file)
    em1.run_em()

main()
