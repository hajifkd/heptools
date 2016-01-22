#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#

from sympy import Matrix, sqrt, zeros, diag
from itertools import izip

class Root(object):

    def __init__(self, q, weight, cartan_matrix):
        self.q = q
        self.weight = weight
        self.cartan_matrix = cartan_matrix

        self.row = sum((w * cartan_matrix.row(j)
                        for j, w in enumerate(self.weight)),
                       zeros(1, len(self.q)))

        self.p = [qi - ri for qi, ri in izip(q, self.row)]

    def ancestors(self):
        ''' get all the possible roots '''
        for i, cand in enumerate(self.p):
            if not cand:
                continue

            # if cand != 0, there must be a root
            nweight = tuple(w + 1 if i == j else w
                            for j, w in enumerate(self.weight))

            yield (nweight, i, self.q[i] + 1)


class Algebra(object):
    pass


class SimpleLieAlgebra(Algebra):

    def __init__(self, cartan_matrix):
        self.cartan_matrix = cartan_matrix
        self.rank = cartan_matrix.rows

    def construct_roots(self):
        # first, we look for all the positive roots
        current_roots = [Root((0,) * i + (2,) + (0,) * (self.rank - 1 - i), 
                              (0,) * i + (1,) + (0,) * (self.rank - 1 - i),
                              self.cartan_matrix) for i in xrange(self.rank)]
        roots = [(0,) * self.rank] + [x.weight for x in current_roots]

        while True:
            ancestors = {}
            for root in current_roots:
                for w, i, qi in root.ancestors():
                    if w in ancestors:
                        ancestors[w][i] = qi
                    else:
                        ancestors[w] = [0] * i + [qi] + [0] * (self.rank - 1 - i)

            if not ancestors:
                break

            current_roots = [None] * len(ancestors)
            for i, w in enumerate(ancestors):
                roots.append(w)
                current_roots[i] = Root(ancestors[w], w, self.cartan_matrix)

        return roots
