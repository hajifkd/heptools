#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#

from sympy import Matrix, ImmutableMatrix, sqrt, zeros, diag, eye
from itertools import izip, combinations


class Weight(object):

    def __init__(self, p, row, cartan_matrix):
        self.p = p
        self.row = row
        self.cartan_matrix = cartan_matrix
        self.q = p + row

    def descendants(self):
        for i, cand in enumerate(self.q):
            if not cand:
                continue

            row = self.q - self.p - self.cartan_matrix.row(i)

            yield (row.as_immutable(), i, self.p[i] + 1)


class Root(object):

    def __init__(self, q, coeff, cartan_matrix):
        self.q = q
        self.coeff = Matrix([coeff])
        self.cartan_matrix = cartan_matrix

        self.row = self.coeff * cartan_matrix

        self.p = q - self.row

    def ancestors(self):
        ''' get all the possible root coefficients '''
        for i, cand in enumerate(self.p):
            if not cand:
                continue

            # if cand != 0, there must be a root
            coeff = self.coeff.copy()
            coeff[i] += 1

            yield (coeff.as_immutable(), i, self.q[i] + 1)


class Algebra(object):
    pass


class SimpleLieAlgebra(Algebra):

    def __init__(self, cartan_matrix):
        self.cartan_matrix = cartan_matrix
        self.rank = cartan_matrix.rows
        self.simple_roots = Matrix(self.construct_simple_roots())

        pcoeff, self.levels = self.positive_root_coeffcients()

        coeff = (pcoeff + [[0] * self.rank] * self.rank +
                 [[-c for c in e] for e in pcoeff])

        self.root_coefficients = Matrix(coeff)
        self.roots = self.root_coefficients * self.simple_roots
        self.fund_weights = self.cartan_matrix.inv() * self.simple_roots
        self.proots = Matrix(self.roots.tolist()[:len(pcoeff)])
        self.delta = sum((self.proots.row(i) for i in xrange(len(pcoeff))),
                         zeros(1, self.rank)) / 2

    def index(self, highest):
        return (self.quad_casimir(highest) * self.dimension(highest) /
                self.roots.rows)

    def dimension(self, highest):
        d = 1
        highest_weight = Matrix([highest]) * self.fund_weights
        for i in xrange(self.proots.rows):
            alpha = self.proots.row(i)
            d *= (1 + alpha.dot(highest_weight) / alpha.dot(self.delta))

        return d

    def quad_casimir(self, highest):
        highest_weight = Matrix([highest]) * self.fund_weights
        return highest_weight.norm()**2 + 2 * highest_weight.dot(self.delta)

    def weights(self, highest):
        current_weights = [Weight(ImmutableMatrix([(0,) * len(highest)]),
                                  ImmutableMatrix([highest]),
                                  self.cartan_matrix)]
        weights = [Matrix([highest]) * self.fund_weights]
        highest_weight = weights[0]
        hd = highest_weight + 2 * self.delta
        d = {highest_weight.as_immutable(): 1}
        depth = 0

        while True:
            depth += 1
            descendants = {}
            degtocalc = []

            for weight in current_weights:
                for r, i, pi in weight.descendants():
                    if r in descendants:
                        descendants[r][i] = pi
                        if r not in degtocalc:
                            degtocalc.append(r)
                    else:
                        descendants[r] = ([0] * i + [pi] +
                                          [0] * (self.rank - 1 - i))

            if not descendants:
                break

            current_weights = [None] * len(descendants)
            for i, r in enumerate(descendants):
                weight = (Matrix([r]) * self.fund_weights).as_immutable()

                if r in degtocalc:
                    # Freudenthal recursion formula
                    deg = 0
                    denom = (hd + weight).dot(highest_weight - weight)

                    for j in xrange(self.proots.rows):
                        alpha = self.proots.row(j)
                        level = self.levels[j]
                        tw = weight
                        # TODO: any root can be a simple root so that no two roots are parallel.
                        # that means once we find k for a difference, no other root must we find.
                        # moreover, we do NOT compare all the element. first, we just see some 
                        # non-zero element and devide it by the one of alpha. if k < depth / level,
                        # then we compare all the element. this reduces #(summing)
                        for k in xrange(1, depth / level + 1):
                            tw += alpha
                            if tw in d:
                                deg += d[tw] * tw.dot(alpha)

                    deg = 2 * deg / denom
                else:
                    deg = 1

                d[weight] = deg

                for j in xrange(deg):
                    weights.append(weight)

                current_weights[i] = Weight(ImmutableMatrix([descendants[r]]),
                                            r, self.cartan_matrix)

        return ImmutableMatrix(weights)

    def positive_root_coeffcients(self):
        # first, we look for all the positive roots
        current_roots = [None] * self.rank
        for i in xrange(self.rank):
            cm = ImmutableMatrix([[0] * i + [1] + [0] * (self.rank - 1 - i)])
            cm2 = (2 * cm).as_immutable()
            current_roots[i] = Root(cm2, cm, self.cartan_matrix)

        positive_roots = [x.coeff for x in current_roots]
        levels = [1] * self.rank
        level = 1

        while True:
            level += 1
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
                levels.append(level)
                current_roots[i] = Root(ImmutableMatrix([ancestors[w]]),
                                        w, self.cartan_matrix)

        return (positive_roots, levels)

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
                        isused[j] = True
                        break

                root[new_dim] = sqrt(lengthsq - root[prev_new_dim]**2)
                simple_roots[i] = root

            for j, e in enumerate(self.cartan_matrix.row(i)):
                if e < 0 and j > i:
                    findroot(j, i, new_dim, length)

        findroot()
        return simple_roots


class A(SimpleLieAlgebra):
    @staticmethod
    def cartan_matrix(n):
        base = eye(n) * 2
        for i in xrange(1, n):
            base[i, i - 1] = -1
            base[i - 1, i] = -1

        return base

    def __init__(self, n):
        super(A, self).__init__(A.cartan_matrix(n))


def su(n):
    return A(n - 1)

# for debug use
su3 = SimpleLieAlgebra(Matrix([[2, -1], [-1, 2]]))
g2 = SimpleLieAlgebra(Matrix([[2, -1], [-3, 2]]))
c3 = SimpleLieAlgebra(Matrix([[2, -1, 0], [-1, 2, -1], [0, -2, 2]]))
