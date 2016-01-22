#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#

from sympy import Matrix, sqrt, zeros, diag
from itertools import izip

class Weight(object):

    def __init__(self, p, row, cartan_matrix):
        self.p = p
        self.row = row
        self.cartan_matrix = cartan_matrix
        self.q = [pi + ri for pi, ri in izip(p, row)]

    def descendants(self):
        for i, cand in enumerate(self.q):
            if not cand:
                continue

            row = tuple(qi - pi - ci
                        for qi, pi, ci in izip(self.q, self.p,
                                               self.cartan_matrix.row(i)))

            yield (row, i, self.p[i] + 1)

class Root(object):

    def __init__(self, q, coeff, cartan_matrix):
        self.q = q
        self.coeff = Matrix([coeff])
        self.cartan_matrix = cartan_matrix

        self.row = self.coeff * cartan_matrix

        self.p = [qi - ri for qi, ri in izip(q, self.row)]

    def ancestors(self):
        ''' get all the possible root coefficients '''
        for i, cand in enumerate(self.p):
            if not cand:
                continue

            # if cand != 0, there must be a root
            coeff = tuple(w + 1 if i == j else w
                            for j, w in enumerate(self.coeff))

            yield (coeff, i, self.q[i] + 1)


class Algebra(object):
    pass


class SimpleLieAlgebra(Algebra):

    def __init__(self, cartan_matrix):
        self.cartan_matrix = cartan_matrix
        self.rank = cartan_matrix.rows
        self.simple_roots = Matrix(self.construct_simple_roots())

        pcoeff = self.positive_root_coeffcients()

        coeff = ([[0] * self.rank] * self.rank + pcoeff +
                 [[-c for c in e] for e in pcoeff])

        self.root_coefficients = Matrix(coeff)
        self.roots = self.root_coefficients * self.simple_roots
        self.fund_weights = self.cartan_matrix.inv() * self.simple_roots

    def weights(self, highest):
        return Matrix(self.weight_coefficients(highest)) * self.fund_weights

    def weight_coefficients(self, highest):
        current_weights = [Weight((0, 0), highest, self.cartan_matrix)]
        weight_coeffs = [highest]
        while True:
            descendants = {}
            for weights in current_weights:
                for r, i, pi in weights.descendants():
                    if r in descendants:
                        descendants[r][i] = pi
                    else:
                        descendants[r] = ([0] * i + [pi] +
                                             [0] * (self.rank - 1 - i))
            if not descendants:
                break

            current_weights = [None] * len(descendants)
            for i, r in enumerate(descendants):
                weight_coeffs.append(r)
                current_weights[i] = Weight(descendants[r], r, self.cartan_matrix)

        return weight_coeffs

    def positive_root_coeffcients(self):
        # first, we look for all the positive roots
        current_roots = [Root((0,) * i + (2,) + (0,) * (self.rank - 1 - i),
                              (0,) * i + (1,) + (0,) * (self.rank - 1 - i),
                              self.cartan_matrix) for i in xrange(self.rank)]
        positive_roots = [x.coeff for x in current_roots]

        while True:
            ancestors = {}
            for root in current_roots:
                for w, i, qi in root.ancestors():
                    if w in ancestors:
                        ancestors[w][i] = qi
                    else:
                        ancestors[w] = ([0] * i + [qi] +
                                        [0] * (self.rank - 1 - i))

            if not ancestors:
                break

            current_roots = [None] * len(ancestors)
            for i, w in enumerate(ancestors):
                positive_roots.append(w)
                current_roots[i] = Root(ancestors[w], w, self.cartan_matrix)

        return positive_roots

    def construct_simple_roots(self):
        simple_roots = [None] * self.rank
        isused = [False] * self.rank
        prev_new_dim = -1

        def findroot(i=0, prev_i=-1, prev_new_dim=-1, prev_len=1):
            if i == 0:
                simple_roots[i] = [1] + [0] * (self.rank - 1)
                isused[i] = True
                new_dim = i
                length = 1
            else:
                angle = sqrt(self.cartan_matrix[i, prev_i] *
                             self.cartan_matrix[prev_i, i]) / -2
                lengthsq = (self.cartan_matrix[i, prev_i] /
                            self.cartan_matrix[prev_i, i])
                length = sqrt(lengthsq)
                inner_prod = angle * prev_len * length
                root = [0] * self.rank
                root[prev_new_dim] = (inner_prod /
                                      simple_roots[prev_i][prev_new_dim])

                for j in xrange(self.rank):
                    if not isused[j]:
                        new_dim = j
                        break

                root[new_dim] = sqrt(lengthsq - root[prev_new_dim]**2)
                simple_roots[i] = root

            for j, e in enumerate(self.cartan_matrix.row(i)):
                if e < 0 and j > i:
                    findroot(j, i, new_dim, length)

        findroot()
        return simple_roots

# for debug use
su3 = SimpleLieAlgebra(Matrix([[2, -1], [-1, 2]]))
g2 = SimpleLieAlgebra(Matrix([[2, -1], [-3, 2]]))
c3 = SimpleLieAlgebra(Matrix([[2, -1, 0], [-1, 2, -1], [0, -2, 2]]))
