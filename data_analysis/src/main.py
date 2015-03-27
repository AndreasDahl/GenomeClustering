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
    greed = load_data("res/muf_greed.csv")
    nogreed = load_data("res/muf_nogreed.csv")
    lru = load_data("res/muf_lru.csv")

    # Remove new clusters from backtrack
    nb = [x for x in nogreed[:,1] if x != 0]
    gb = [x for x in greed[:,1] if x != 0]
    lb = [x for x in lru[:,1] if x != 0]

    plt.figure()
    plt.title("Best Hits")
    plt.hist(nogreed[:,0], 50, color='g')

    plt.figure()
    plt.title("Greedy Hits")
    plt.hist(greed[:,0], 50, color='g')

    plot_lookups(nb, "Full lookups")
    plot_lookups(gb, "Greedy lookups")
    plot_lookups(lb, "LRU lookups")

    plt.show()