"""
    Copyright 2015 Andreas Dahl, Christian Muf

    This file is part of MufDahlClust.

    MufDahlClust is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MufDahlClust is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MufDahlCLust.  If not, see <http://www.gnu.org/licenses/>.
"""

__author__ = 'Andreas Dahl'

import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.optimize import curve_fit
from collections import Counter


def load_data(file_path, delimiter=' '):
    """
    Load a data file into a numpy array.
    :param file_path:  File path of the data file to load.
    :return:  Numpy array of loaded data.
    """
    with open(file_path) as f:
        reader = csv.reader(f, delimiter=delimiter)
        data = []
        for row in reader:
            data.append([value for value in row])
        return np.array(data)


def plot_hits(data, name='hits'):
    n = 10
    plt.figure()
    plt.title(name)
    plt.xlabel("Cluster Index")
    plt.ylabel("Hit Count")
    plt.hist(data, 50, color='g')
    common = Counter(data).most_common(n)
    median = np.median(data)
    per95 = np.percentile(data, 95)
    plt.plot([median, median], plt.ylim(), c='b')
    plt.plot([per95, per95], plt.ylim(), c='r')

    print "Stats for '%s'" % name
    print "median: %.2f" % median
    print "95th percentile: %.2f" % per95
    print "Most common clusters:", common
    print "%.2f%% of values hit the same %d cluster" % (float(sum([v for _, v in common])) / len(data) * 100, n)
    print


def plot_lookups(data, name='lookups'):
    plt.figure()
    plt.title(name)
    plt.xlabel("Lookups to hit")
    plt.hist(data, 50)
    median = np.median(data)
    per95 = np.percentile(data, 95)
    print "%s 95th percentile %.2f" % (name, per95)
    plt.plot([median, median], plt.ylim(), c='g')
    plt.plot([per95, per95], plt.ylim(), c='r')


def analyse_data():
    # load data
    full = load_data("res/length_muf_full.csv")
    lkmer = load_data("res/lkmer_muf_full.csv")
    kmer = load_data("res/kmer_muf_full.csv")
    greed = load_data("res/length_muf_greedy.csv")
    lru = load_data("res/length_muf_lru.csv")

    full.astype(int)
    lkmer.astype(int)
    kmer.astype(int)
    greed.astype(int)
    lru.astype(int)

    # Remove new clusters from backtrack
    nb = [x for x in full[:, 1] if x != 0]
    gb = [x for x in greed[:, 1] if x != 0]
    lb = [x for x in lru[:, 1] if x != 0]

    plot_hits(full[:, 0], "Best Hits")
    plot_hits(kmer[:, 0], "Kmer Hits")
    plot_hits(lkmer[:, 0], "Length-Kmer Hits")
    plot_hits(greed[:, 0], "Greedy Hits")

    plot_lookups(nb, "Full lookups")
    plot_lookups(gb, "Greedy lookups")
    plot_lookups(lb, "LRU lookups")


def analyse_comparisons():
    data = load_data("res/compare.csv")
    data.astype(float)
    plt.axis([0, 1, 0, 1])
    plt.xlabel("k-mers")
    plt.ylabel("levenshtein")

    plt.scatter(data[:, 0], data[:, 1], marker='x')
    plt.plot(plt.xlim(), plt.ylim(), c='r')


def analyse_comparisons2():
    data = load_data("res/compare3.csv", ';')
    data = data.astype(float)
    plt.figure()
    plt.axis([0.8, 1, 0.8, 1])
    plt.plot([0.95, 0.95], plt.ylim(), c='black')
    plt.plot(plt.xlim(), [0.95, 0.95], c='black')
    plt.xlabel("levenshtein")
    plt.ylabel("k-mer-dac")

    # Distance comparison
    plt.scatter(data[:, 2], data[:, 0], marker='o')
    plt.plot(plt.xlim(), plt.ylim(), c='r')


    # Box plots
    plt.figure()
    plt.boxplot(data[:,[1,3]], whis=[5,95])

    print "MufDifference Speeds\nMedian: %f\nAvg: %f" % (np.median(data[:,1]), np.average(data[:,1]))
    print "Levenshtein Speeds\nMedian: %f\nAvg: %f" % (np.median(data[:,3]), np.average(data[:,3]))

    print "Levenshtein takes", np.median(data[:, 3]) / np.median(data[:, 1]), "times the time of MufDifference"

    print [(muf, lev) for [muf, _, lev, _] in data]


def pow_fit(x, a, b, c):
    return a + b * (x ** c)


def hyp_fit(x, a, b):
    return a / x + b


def r_squared(expected, actual):
    ssreg = np.sum((expected - actual) ** 2)
    sstot = np.sum((np.mean(actual) - actual) ** 2)
    return 1.0 - ssreg / sstot

def cache_analysis():
    data = load_data("../stats.csv", ";")
    data = data.astype(float)

    cache_total = data[:,0] + data[:,1]
    time_data = data[:,2] / 1000000

    c = []
    for i in range(len(data[:,0])):
        l = data[i,0]
        r = data[i,1]
        if (l > r):
            c.append((r - l) / max(data[:,1]))
        else:
            c.append(-(l - r) / max(data[:,0]))

    popt, pcov = curve_fit(pow_fit, cache_total, time_data)

    x = np.arange(0, 300, 1)


    r2 = r_squared(pow_fit(cache_total, *popt), time_data)
    print "R squared:", r2
    plt.figure("Duration")
    plt.xlabel("Total Cache Size")
    plt.ylabel("Clustering Duration in seconds")
    plt.axis([0, 300, 0, 300])
    plt.scatter(cache_total, data[:,2] / 1000000, marker='o', c=c)
    fit = plt.plot(x, pow_fit(x, *popt), c="grey",
                   label="%f + %fx^%f, R^2: %f" % (popt[0], popt[1], popt[2], r2))
    plt.legend(handles=fit, loc=4)
    plt.gray()

    popt, pcov = curve_fit(hyp_fit, cache_total, data[:,3])
    r2 = r_squared(hyp_fit(cache_total, *popt), data[:,3])
    print "R squared:", r2
    plt.figure("Cluster Count")
    plt.xlabel("Total Cache Size")
    plt.ylabel("Resulting Number of Clusters")
    plt.axis([0, 300, 0, 150000])
    plt.scatter(cache_total, data[:,3], marker='o')
    fit = plt.plot(x, hyp_fit(x, *popt), c="grey",
                   label="%.02f/x + %.02f,  R^2: %.02f" % (popt[0], popt[1], r2))
    plt.legend(handles=fit, loc=1)


def cluster_quality():
    data = load_data("res/caecum.uc", "\t")
    clusters = data[:,1].astype(int)

    labels, values = zip(*Counter(clusters).items())

    plt.hist(values, 100)



if __name__ == "__main__":
    # analyse_comparisons2()
    cache_analysis()
    # cluster_quality()
    plt.show()