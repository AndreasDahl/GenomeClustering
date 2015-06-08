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

from __future__ import division
import random
import numpy as np
from matplotlib import pyplot

__author__ = 'Andreas Dahl'


l = 10
k = 2
n = 5


class Tree(object):
    value = 0
    parent = None
    coverage = []
    branches = []

    def __init__(self, parent, value, coverage):
        self.parent = parent
        self.value = value
        self.coverage = coverage
        self.branches = []

    def __str__(self):
        return "Tree:\n%s\n%s\n%s" % (self.value, self.coverage, self.branches)

    def add_branch(self, branch):
        self.branches.append(branch)

    def total_value(self):
        if (self.parent != None):
            return self.value + self.parent.total_value()
        else:
            return self.value


def find_hits(spot, coverage, k):
    new_coverage = [t for t in coverage]
    for target in range(spot, spot + k):
        new_coverage[target] = True
    return new_coverage


def calculate(l, k, n):
    spots = l - (k - 1)

    root = Tree(None, 0, [False] * l)
    branches = [root]

    print "String Length: %d" % l

    for i in range(n):
        new_branches = []
        for branch in branches:
            for spot in range(spots):
                new_coverage = find_hits(spot, branch.coverage, k)
                new_hits = np.sum(new_coverage) - np.sum(branch.coverage)
                new_branch = Tree(branch, new_hits, new_coverage)
                branch.add_branch(new_branch)
                new_branches.append(new_branch)
        branches = new_branches

    # for i in range(n):
    print "leaves:", len(branches)
    total = 0
    full_total = 0
    for branch in branches:
        total += branch.value
        full_total += branch.total_value()
    print "level: %.2f" % (total / len(branches))
    print "full : %.2f" % (full_total / len(branches))
    print "-----------"
    return (full_total / len(branches))


def test_k():
    pyplot.figure("K")
    x = range(1, 11)
    y = [calculate(l, cur, n) / l for cur in x]
    pyplot.plot(x, y)


def test_l():
    pyplot.figure("L")
    x = range(k, 15)
    y = [calculate(l, k, n) / l for l in x]
    pyplot.plot(x, y)


def test_n():
    pyplot.figure("N")
    x = range(1, 6)
    y = [calculate(l, k, n) / l for n in x]
    pyplot.plot(x, y)


if __name__ == "__main__":
    test_k()
    test_l()
    test_n()
    pyplot.show()
