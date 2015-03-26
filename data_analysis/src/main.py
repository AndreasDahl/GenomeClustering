__author__ = 'Andreas Dahl'

import csv
import numpy as np
import matplotlib.pyplot as plt

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


if __name__ == "__main__":
    # load data
    greed = load_data("res/muf_greed.csv")
    nogreed = load_data("res/muf_nogreed.csv")
    lru = load_data("res/muf_lru.csv")

    # Remove new clusters from backtrack
    nb = [x for x in nogreed[:,1] if x != 0]
    gb = [x for x in greed[:,1] if x != 0]
    lb = [x for x in lru[:,1] if x != 0]

    print "Full lookups median", np.median(nb)
    print "Greedy lookups median", np.median(gb)
    print "LRU lookups median", np.median(lb)

    print "LRU lookups 95% percentile", np.percentile(lb, 95)

    plt.figure()
    plt.title("Best Hits")
    plt.hist(nogreed[:,0], 50, color='g')

    plt.figure()
    plt.title("Greedy Hits")
    plt.hist(greed[:,0], 50, color='g')

    plt.figure()
    plt.title("Full backtrack")
    plt.hist(nb, 50)

    plt.figure()
    plt.title("Greedy backtrack")
    plt.hist(gb, 50)

    plt.figure()
    plt.title("LRU backtrack")
    plt.hist(lb, 50)

    plt.show()