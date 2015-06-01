__author__ = 'Andreas Dahl'

import csv
import numpy as np
import matplotlib.pyplot as plt
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
    data = load_data("res/muf_lev_silva_unfinal_200.csv", ';')
    data = data.astype(float)
    plt.figure()
    plt.axis([0, 1, 0, 1])
    plt.plot([0.95, 0.95], plt.ylim(), c='r')
    plt.plot(plt.xlim(), [0.95, 0.95])
    plt.xlabel("levenshtein")
    plt.ylabel("k-mer")

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




if __name__ == "__main__":
    analyse_comparisons2()

    plt.show()