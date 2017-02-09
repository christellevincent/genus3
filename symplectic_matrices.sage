"""
Recip -- Sage package for using Shimura's reciprocity law

See the file README.txt for version information and instructions.

#*****************************************************************************
# Copyright (C) 2010,2011,2012,2013,2016 Marco Streng <marco.streng@gmail.com>
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#*****************************************************************************

This file implements the arithmetic of matrices in GSp_2g(ZZ), see the class

 * GSp_element

"""

from sage.matrix.matrix import is_Matrix

#from sage.matrix.all import (is_Matrix, Matrix, identity_matrix)
#from sage.rings.all import ZZ
#from sage.structure.factorization import Factorization


def diag(a):
    """
    Given a square matrix a, returns the list of diagonal entries of a.

    EXAMPLES::

        sage: load("recip.sage")
        sage: g = Matrix(Zmod(8),[(5,0,7,1),(6,6,0,7),(1,2,2,3),(2,3,5,4)])
        sage: diag(g)
        [5, 6, 2, 4]
    """
    n = a.ncols()
    assert n == a.nrows()
    return [a[k,k] for k in range(n)]


def matrix_from_blocks(A, B, C, D):
    """
    Given matrices A,B,C,D, returns a single matrix
    ( A B )
    ( C D )
    """
    m1 = A.nrows()
    m2 = C.nrows()
    if not m1 == B.nrows():
        raise ValueError, "A and B do not have the same number of rows in " \
                          "matrix_from_blocks A = %s, B = %s." % (A, B)
    if not m2 == D.nrows():
        raise ValueError, "C and D do not have the same number of rows in " \
                          "matrix_from_blocks C = %s, D = %s" % (C, D)
    if not A.ncols() == C.ncols():
        raise ValueError, \
              "A and C do not have the same number of columns in " \
              "matrix_from_blocks A = %s, C = %s" % (A, C)
    if not B.ncols() == D.ncols():
        raise ValueError, "B and D do not have the same number of columns " \
                          "in matrix_from_blocks B = %s, D = %s" % (B, D)
    return Matrix([Sequence(A[i]) + Sequence(B[i]) for i in range(m1)] +\
                  [Sequence(C[i]) + Sequence(D[i]) for i in range(m2)])


def zero_matrix(m, n = None, base_ring = ZZ):
    """
    Returns the m x n matrix with all zero coefficients over base_ring.

    If n is unspecified, then n=m,
    if base_ring is unspecified, then base_ring = ZZ.
    """
    if n == None:
        n = m
    return Matrix([[base_ring(0) for j in range(n)] for i in range(m)])

    
def Omega(g):
    """
    Returns the 2g x 2g matrix given in terms of gxg blocks as
    ((0, -1) (1, 0))
    """
    one = identity_matrix(g)
    zero = zero_matrix(g, g)
    return matrix_from_blocks(zero, -one, one, zero)

    
def nu(m):
    """
    Given a 2g x 2g matrix m, returns a scalar nu such that
    m.transpose() * Omega(g) * m = nu * Omega(g) if it exists.
    Otherwise returns None.
    """
    g = ZZ(m.ncols()/2)
    if 2*g != m.nrows():
        raise ValueError, "Non-square matrix in nu: %s" % m
    a = m.transpose()*Omega(g)*m
    n = a[g,0]
    if a != n*Omega(g):
        return None
    return n


def ABCD(M, m = None, n = None):
    """
    Returns matrices A, B, C, D such that M equals
    ( A B )
    ( C D )
    and A is mxn.
    If m (resp. n) is not specified, then it is half the number
    of rows (resp. columns) of M.
    """
    if m == None:
        m = M.nrows()
        if m % 2 == 1:
            raise ValueError, "m not specified in ABCD while M (=Matrix(%s))" \
                              " has an odd number of rows" % Sequence(M)
        m = ZZ(m/2)
    if n == None:
        n = M.ncols()
        if n % 2 == 1:
            raise ValueError, "n not specified in ABCD while M (=Matrix(%s))" \
                              " has an odd number of columns" % Sequence(M)
        n = ZZ(n/2)
    M.subdivide(m,n)
    A = M.subdivision(0,0)
    B = M.subdivision(0,1)
    C = M.subdivision(1,0)
    D = M.subdivision(1,1)
    return A,B,C,D

    
def Zmod_to_SL(a):
    """
    Given a in (ZZ/mu*ZZ)^*, outputs A in SL_2(ZZ) with
    A = ((a 0) (0 a^-1)) mod mu
    """
    mu = a.parent().order()
    if not a.parent() == Zmod(mu):
        raise ValueError, "a (=%s) not in ZZ/muZZ for mu = a.parent()." \
                          "order() = %s" % (a, mu)
    p = lift_small(a)
    q = lift_small(a**-1)
    e = 1 - p*q
    A = Matrix([[p, -e], [e, q*(1+e)]])
    assert mat_convert(A, a.parent()) == Matrix([[a, 0],[0, a**-1]])
    assert A.determinant() == 1
    assert A.base_ring() == ZZ
    return A


class Sp_group:
    """
    The group Sp(n, R) for R = Zmod(m).
    
    EXAMPLES::
    
        sage: l = Sp_group(4, 2).list() # long time, 2 seconds
        sage: len(l) # depends on previous line with long time
        720

    """
    _n = None
    _m = None
    _gens = None
    _elements_known = None
    _elements_not_exhausted = None
    
    def __init__(self, n, m):
        if not (n in ZZ and n > 0 and m in ZZ and m > 0 and ZZ(n) % 2 == 0):
            raise ValueError, "n (=%s) and m (=%s) must be positive integers" \
                              " with n even" % (n, m)
        self._n = ZZ(n)
        self._m = ZZ(m)
        R = Zmod(m)
        self._gens = [g.matrix(R) for g in symplectic_generators(n/2)] + \
                    [g.matrix(R).inverse() for g in symplectic_generators(n/2)]
        id = identity_matrix(ring=Zmod(m), n=n)
        self._elements_not_exhausted = [id]
        self._elements_known = [id]
        
    def list(self):
        """
        Returns a list of all elements of self. Slow!
        
        This uses that Sp_2n(ZZ) maps surjectively to Sp_2n(ZZ/mZZ),
        hence the generators of Sp_2n(ZZ) generate Sp_2n(ZZ/mZZ)

        EXAMPLES::
        
            sage: Sp_group(2, 2).list()
            [
            [1 0]  [0 1]  [1 1]  [1 1]  [1 0]  [0 1]
            [0 1], [1 0], [0 1], [1 0], [1 1], [1 1]
            ]

        """
        elements_not_exhausted = self._elements_not_exhausted
        gens = self._gens
        elements_known = self._elements_known
        while elements_not_exhausted != []:
            e = elements_not_exhausted.pop()
            for g in gens:
                h = e*g
                if not h in elements_known:
                    elements_known.append(h)
                    elements_not_exhausted.append(h)
            if get_verbose():
                print "found: %s, to check: %s" % (len(elements_known),
                                                   len(elements_not_exhausted))
        return elements_known
        
    def order(self):
        """
        Returns the order of self. Way too slow!!
        
        EXAMPLES::
        
            sage: Sp_group(2, 5).order()
            120
        """
        return len(self.list())
            

class GSp_group:
    """
    The group GSp(n, R) for R = Zmod(m).
    
    EXAMPLES::
    
        sage: l = GSp_group(4, 2).list() # long time, 2 seconds
        sage: len(l) # long time
        720

    """
    _n = None
    _m = None
    _S = None
    
    def __init__(self, n, m):
        self._n = ZZ(n)
        self._m = ZZ(m)
        self._S = Sp_group(n, m)
        
    def list(self):
        """
        Returns a list of all elements of self. Slow!
        
        EXAMPLES::
        
            sage: GSp_group(2, 2).list()
            [
            [1 0]  [0 1]  [1 1]  [1 1]  [1 0]  [0 1]
            [0 1], [1 0], [0 1], [1 0], [1 1], [1 1]
            ]

        """
        ret = []
        m = self._m
        for a in Zmod(m):
            if gcd(a, m) == 1:
                for M in self._S.list():
                    ret.append(GSp_element(M, a))
        return ret
        
    def order(self):
        """
        Returns the order of self. Way too slow!!
        
        EXAMPLES::
        
            sage: GSp_group(2, 5).order() # long time, 1 second
            480
        """
        return len(self.list())


class GSp_element:
    """
    A matrix
        M = iota(nu^-1) * S
    where iota(t) = diag(1,1,...,1,t^-1,t^-1,...,t^-1) and S is in Sp_{2g}
    It acts from the right on theta_ring, where iota(t) acts as t on
    coefficients, and S acts via the usual action.
    """
    _matrix = None
    _Sp_part = None
    _nu = None
    _g = None

    def __init__(self, arg1, arg2 = None, ring = None):
        """
        INPUT:
        
         - `ring` a ring or None (to derive the ring from `arg1`
        
         - `arg1` a GSp-matrix, or
         
         - `arg1` an Sp-matrix and `arg2` an element of the unit group of `ring`
         
        """
        if not is_Matrix(arg1):
            arg1 = Matrix(arg1, ring = ring)
        elif ring != None:
            arg1 = mat_convert(arg1, ring)
        g = arg1.ncols()/2
        if not g in ZZ:
            raise ValueError, "Number of columns of arg1 (=%s) must be even in GSp_element" % arg1
        g = ZZ(g)
        if not arg1.nrows() == 2*g:
            raise ValueError, "arg1 (=%s) must be a square matrix in GSp_element" % arg1
        self._g = g
        assert 2*g == arg1.nrows()
        if arg2 == None:
            self._matrix = arg1
            self._nu = nu(arg1)
            if self._nu == None:
                raise ValueError, "Input matrix arg1 (= %s) does not satisfy arg1.transpose()*Omega*arg1 = Omega*nu for any nu" % arg1
            self._Sp_part = diagonal_matrix([1 for i in range(g)] + [self._nu**-1 for i in range(g)]) * arg1
        else:
            if ring != None:
                arg2 = ring(arg2)
            else:
                assert sage.structure.element.is_RingElement(arg2)
            assert nu(arg1) == 1
            self._nu = arg2
            self._Sp_part = arg1
            self._matrix = diagonal_matrix([1 for i in range(g)] + [self._nu for i in range(g)]) * arg1
          
          
    def g(self):
        return self._g

    def _matrix_(self, base_ring = None):
        if base_ring == None:
            return self._matrix
        else:
            return self._matrix._matrix_(base_ring)
        
    matrix = _matrix_

    def nu(self):
        return self._nu

    def Sp_part(self):
        return self._Sp_part

    def __str__(self):
        return self.matrix().__str__()
        
    def __repr__(self):
        return self.matrix().__repr__()

    def __call__(self, tau):
        if is_Matrix(tau):
            return Sp_action(self.matrix(), tau)
#        if tau in theta_ring:
#            f = cycl_galois_action_on_theta(self.nu()**-1, tau)
#            return Sp_action_on_forms(self.Sp_part(), f)

    @cached_method
    def action_on_theta_generators(self, den):
        """
        Returns a sequence L such that L[i] is the image of the i-th generator
        of the theta ring under self.
        
        EXAMPLES::
        
            sage: load("recip.sage")
            sage: M = GSp_element([[7,0,2,6],[0,7,6,0],[0,6,5,0],[6,6,0,5]], ring=Zmod(8))
            sage: M.Sp_part()
            [7 0 2 6]
            [0 7 6 0]
            [0 2 3 0]
            [2 2 0 3]
            sage: M.action_on_theta_generators(2)
            [t0, t1, (zeta8^2)*t2, (-zeta8^2)*t3, (-zeta8^2)*t4, (-zeta8^2)*t5, t6, -t7, t8, t9, (zeta8^2)*t10, (-zeta8^2)*t11, (zeta8^2)*t12, (zeta8^2)*t13, -t14, t15]
        """
        g = self._g
        M = self.matrix()
        B = M.base_ring()
        if B == QQ:
            M = mat_convert(M, ZZ)
        elif B != ZZ:
            mu = B.order()
            if not (B == Zmod(mu) and mu % 2*den**2 ==0 and mu % 8 == 0):
                raise ValueError
            M = mat_convert(M, lift_small)
        nu_inv = lift_small(self.nu()**-1)
        ds = [theta_action_without_kappa(M, nu_inv, num_to_c(i, g, den)) \
              for i in range(den**(2*g))]
        P = theta_ring(g, den)[0]
        zeta = P.base_ring().gen()
        k = zeta.multiplicative_order()
        ims = [zeta**(k*a[1]) * P.gens()[c_to_num(a[0], den)] for a in ds]
        return ims

    def parent(self):
        return self.matrix().parent()

    def __pow__(self, n):
        """
        Return self to the power n
        """
        return GSp_element(self.matrix().__pow__(n))

    def __eq__(self, rhs):
        """
        Returns True if self equals rhs as a matrix, and False otherwise.
        """
        return self._matrix_() == rhs._matrix_(rhs.base_ring())

    def __mul__(self, rhs):
        """
        Returns the product self * rhs.
        """
        return GSp_element(self._matrix_() * rhs._matrix_(rhs.base_ring()))

    def base_ring(self):
        """
        Returns the base ring of self.
        """
        return self._matrix_().base_ring()
        
    def transpose(self):
        """
        Returns the transpose of self.
        
        EXAMPLE::
        
            sage: load("recip.sage")
            sage: M = GSp_element([[0,-1,1,1],[5,0,4,-5],[0,-2,2,1],[1,0,1,-1]])
            sage: M.transpose()
            [ 0  5  0  1]
            [-1  0 -2  0]
            [ 1  4  2  1]
            [ 1 -5  1 -1]

        """
        return GSp_element(self._matrix_().transpose())


def symplectic_generators(g, subgroup=None, level=None):
    """
    Returns a set of generators of Sp_2g(ZZ)
    
    WARNING: I doubt the correctness of this function, I don't get all the
    invariants that I expected to get.
    """
    if g == 1:
        return [GSp_element(Matrix([[0,1],[-1,0]])),
                GSp_element(Matrix([[1,1],[0,1]]))]
    if g == 2:
        ret = []
        Sbig = GSp_element(Matrix([[0,0,1,0],[0,0,0,1],[-1,0,0,0],[0,-1,0,0]]))
        T1 = GSp_element(Matrix([[1,0,1,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]))
        T2 = GSp_element(Matrix([[1,0,0,1],[0,1,1,0],[0,0,1,0],[0,0,0,1]]))
        T3 = GSp_element(Matrix([[1,0,0,0],[0,1,0,1],[0,0,1,0],[0,0,0,1]]))
        Tul = GSp_element(Matrix([[1,1,0,0],[0,1,0,0],[0,0,1,0],[0,0,-1,1]]))
        neg = GSp_element(Matrix([[1,0,0,0],[0,-1,0,0],[0,0,1,0],[0,0,0,-1]]))
        Sul = GSp_element(Matrix([[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]]))
                
        if subgroup is None:
            assert level is None
            return [Sbig, T1, T2, T3, Tul, neg, Sul]
        n = level
        assert not level is None
        T1n = GSp_element(Matrix([[1,0,n,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]))
        T2n = GSp_element(Matrix([[1,0,0,n],[0,1,n,0],[0,0,1,0],[0,0,0,1]]))
        T3n = GSp_element(Matrix([[1,0,0,0],[0,1,0,n],[0,0,1,0],[0,0,0,1]]))
        if subgroup == '_0':
            return [T1, T2, T3, Tul, neg, Sul, T1n.transpose(), T2n.transpose(), T3n.transpose()]
        if subgroup == '^0':
            return [T1n, T2n, T3n, Tul, neg, Sul, T1.transpose(), T2.transpose(), T3.transpose()]
        Tuln = GSp_element(Matrix([[1,n,0,0],[0,1,0,0],[0,0,1,0],[0,0,-n,1]]))
        if subgroup == '_1' or subgroup == '^1':
            ret = [T1, T2, T3, Tuln, Tuln.transpose(), T1n.transpose(), T2n.transpose(), T3n.transpose()] + [neg for i in range(level<=2)] 
            if subgroup == '_1':
                return ret
            return [g.transpose() for g in ret]
        if subgroup == '':
            return [T1n, T2n, T3n, Tuln, Tuln.transpose(), T1n.transpose(), T2n.transpose(), T3n.transpose()] + [neg for i in range(level<=2)]
            
        raise ValueError
                
    raise NotImplementedError, "symplectic_generators only implemented for g=2, not for g=%s" % g


def group_generators_to_list(gens, G = None):
    """
    Given a finite list `gens` of elements of a group `G`,
    assuming `gens` generates a finite subgroup `H` of `G`,
    returns a list of all elements of `H`.
    """
    if G != None:
        unit_element = G.one()
    else:
        unit_element = gens[0].parent().one()
    H = [unit_element]
    Hpreviousround = [unit_element]
    while len(Hpreviousround) > 0:
        Hthisround = []
        for h in Hpreviousround:
            for g in gens:
                for e in [1, -1]:
                    n = (g**e)*h
                    if not n in H:
                        Hthisround += [n]
                        H += [n]
        if get_verbose() == 2:
            print "Found %s new period matrices in H" % len(Hthisround)
        Hpreviousround = Hthisround
    return H


def minimal_generating_subset(gens, G=None, order=None, k=0):
    """
    Given a finite list `gens` of elements of a group `G`,
    assuming `gens` generates a finite subgroup `H` of `G`,
    returns a minimal subset `s` of `gens` that generates `H`.

    Here minimal is with respect to inclusion, so it may
    not be the smallest such subset.

    `order` is the order of `H`, and `s` must gens[:k]
    """
    if G == None:
        G = gens[0].parent()
    if order == None:
        order = len(group_generators_to_list(gens, G))
    for l in range(k, len(gens)):
        s = gens[:l] + gens[l+1:]
        if len(group_generators_to_list(s, G)) == order:
            return minimal_generating_subset(s, G, order, l)
    return gens
    

def random_symplectic_matrix(g, n, subgroup=None, level=None):
    """
    Returns the product of n random elements of symplectic_generators(g)
    and inverses of such.
    
    See symplectic_generators for which subgroups are possible.
    
    EXAMPLE::
    
        sage: load("recip.sage")
        sage: M = random_symplectic_matrix(2, 20)
        sage: M # random output
        [ 0 -1  1  1]
        [ 5  0  4 -5]
        [ 0 -2  2  1]
        [ 1  0  1 -1]
        sage: is_symplectic(M)
        True

    """
    gens = symplectic_generators(g, subgroup, level)
    gens = gens + [a**-1 for a in gens]
    ret = GSp_element(identity_matrix(2*g))
    for i in range(n):
        ret = ret * gens[randrange(len(gens))]
    return ret





