import sys
import csv
import math
from math import log as log

def read_qt(qt_file):
    qts = []
    with open(qt_file) as qtfh:
        qtfh_dict = csv.DictReader(qtfh, delimiter = ",")
        for line in qtfh_dict:
            qts.append(float(line['QTp']))
    return(qts)

def ml_grid_search(qts):
    max_lik = float("-inf")
    for q in [0.6, 0.4]:
    #for q in [x / 100 for x in range(10, 100, 10)]:
        for mu1 in [14, 10, -5]:
        #for mu1 in [x / 10 for x in range(-10, 11, 1)]:
            #for mu2 in [x / 10 for x in range(-10, 11, 1)]:
            for mu2 in [14, -5, 10]:
                #for mu3 in [x / 10 for x in range(-10, 11, 1)]:
                for mu3 in [10, -5, 14]:
                    #if mu1 == mu2 or mu2 == mu3 or mu1 == mu3:
                        #continue
                    for sigma in [0.1]:
                        sigmasq = sigma ** 2
                        #print(q, mu1, mu2, mu3, sigmasq)
                        params = "\t".join([str(q), str(mu1), str(mu2), str(mu3), str(sigmasq)])
                        lik = 0
                        for x in qts:
                            lik1 = 0
                            #print(mu3, log(1.0 - q), ((x - mu3)**2)/(2 * sigmasq))
                            #Hom Ref
                            lik += log(1.0/6.0)- ((x - mu1)**2)/(2 * sigmasq) - 0.5 * log(2 * math.pi * sigmasq)
                            #print("lik1", 2 * log(1.0 - q) - ((x - mu1)**2)/(2 * sigmasq) - 0.5 * log(2 * math.pi * sigmasq))
                            #Het
                            lik += log(2.0/6.0) - ((x - mu2)**2)/(2 * sigmasq) - 0.5 * log(2 * math.pi * sigmasq)
                            #print("lik2", log(2) + log(1.0 - q) + log(q) - ((x - mu2)**2)/(2 * sigmasq) - 0.5 * log(2 * math.pi * sigmasq))
                            #Hom Alt
                            lik += log(3.0/6.0) - ((x - mu3)**2)/(2 * sigmasq) - 0.5 * log(2 * math.pi * sigmasq)
                            #print("lik3", 2 * log(q) - ((x - mu3)**2)/(2 * sigmasq) - 0.5 * log(2 * math.pi * sigmasq))
                        if lik > max_lik:
                            max_lik = lik
                            maxlik_params = params
                        print(lik, params)
    return maxlik_params

def main():
    qt_file = sys.argv[1]
    qts = read_qt(qt_file)
    maxlik_params = ml_grid_search(qts)

main()
