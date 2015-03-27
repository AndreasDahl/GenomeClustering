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
            data.append([int(value) for value in row])
        return np.array(data)


def plot_hits(data, name='hits'):
    n = 10
    plt.figure()
    plt.title(name)
    print len(data)
    plt.hist(data, 50, color='g')
    common = Counter(data).most_common(n)
    median = np.median(data)
    per95 = np.percentile(data, 95)
    print "%s 95 percentile %.2f" % (name, per95)
    plt.plot([median, median], plt.ylim(), c='b')
    plt.plot([per95, per95], plt.ylim(), c='r')
    print name, "most common clusters:", common
    print "%.2f%% of values hit the same %d cluster" % (float(sum([v for _, v in common])) / len(data) * 100, n)


def plot_lookups(data, name='lookups'):
    plt.figure()
    plt.title(name)
    plt.hist(data, 50)
    median = np.median(data)
    per95 = np.percentile(data, 95)
    print "%s 95 percentile %.2f" % (name, per95)
    plt.plot([median, median], plt.ylim(), c='g')
    plt.plot([per95, per95], plt.ylim(), c='r')


if __name__ == "__main__":
    # load data
    full = load_data("res/length_muf_full.csv")
    lkmer = load_data("res/dkmer_muf_full.csv")
    kmer = load_data("res/kmer_muf_full.csv")
    greed = load_data("res/length_muf_greedy.csv")
    lru = load_data("res/length_muf_lru.csv")

    # Remove new clusters from backtrack
    nb = [x for x in full[:,1] if x != 0]
    gb = [x for x in greed[:,1] if x != 0]
    lb = [x for x in lru[:,1] if x != 0]

    plot_hits(full[:,0], "Best Hits")
    plot_hits(kmer[:,0], "Kmer Hits")
    plot_hits(lkmer[:,0], "Length-Kmer Hits")
    plot_hits(greed[:,0], "Greedy Hits")

    plot_lookups(nb, "Full lookups")
    plot_lookups(gb, "Greedy lookups")
    plot_lookups(lb, "LRU lookups")

    plt.show()