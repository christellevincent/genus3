
#*****************************************************************************
# Copyright (C) 2016 Marco Streng <marco.streng@gmail.com> and 
#                   Pinar Kilicer <pinarkilicer@gmail.com>
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

# source url: https://bitbucket.org/pkilicer/period-matrices-for-genus-3-cm-curves/src/master/period_matrices_genus3.sage


def period_matrix(bas, embs):
    M = Matrix([[phi(b) for b in bas] for phi in embs])
    Omega1 = M[:,:3]
    Omega2 = M[:,3:]
    return Omega2.inverse()*Omega1
    
def split(bas):
    g = ZZ(bas.degree()/2)
    return vector(bas[:g]), vector(bas[g:])
    
def join(bas1,bas2):
    return vector(list(bas1)+list(bas2))
    
def B_to_M(B):
    g = B.nrows()
    assert B.ncols() == g
    return ABCD_to_M(identity_matrix(g), B, zero_matrix(g, g), identity_matrix(g))
    
def ABCD_to_M(A,B,C,D):
    g = A.ncols()
    for E in [A,B,C,D]:
        assert E.nrows() == E.ncols() == g
    return Matrix([join(A[i],B[i]) for i in range(A.nrows())] +
                  [join(C[i],D[i]) for i in range(A.nrows())])

def reduce_Siegel_epsilon(bas, embs, verbose=False):
    """
    Given a symplectic basis and the embeddings from a CM-type, reduces to the domain F_S^epsilon of Kilicer et al:
    
      * Re(tau) has coeffs in [-1/2,1/2]
      * Im(tau) is LLL-reduced
      * Im(tau[0,0]) >= sqrt(3/4) - 0.01

    returns a new basis and a basis transformation
    """
    input_bas = bas
    g = ZZ(bas.degree()/2)
    M = identity_matrix(2*g)
    
    while True:
        # the imaginary part:
        bas, N = reduce_imag_part_LLL(bas, embs)
        M = N * M
        assert act_by_M_on_basis(M, input_bas) == bas
        if verbose:
            print "imag reduction:"
            print N
        # the real part:
        bas, N = reduce_real_part(bas, embs)
        M = N * M
        assert act_by_M_on_basis(M, input_bas) == bas
        if verbose:
            print "real reduction:"
            print N
        # the upper left imaginary entry
        tau = period_matrix(bas, embs)
        if abs(tau[0,0]) < 0.99:
            N = M_for_reduce_using_z11(g)
            M = N * M
            bas = act_by_M_on_basis(N, bas)
            assert act_by_M_on_basis(M, input_bas) == bas
            if verbose:
                print "y11 reduction:"
                print N

        else:
            assert act_by_M_on_basis(M, input_bas) == bas

            return bas, M
            
def reduce_real_part(bas, embs):
    """
    Given a symplectic basis and the embeddings from a CM-type,
    change the basis such that
    
      * Re(tau) has coeffs in [-1/2,1/2]
    
    returns a new basis and a basis transformation
    """
    tau = period_matrix(bas, embs)
    g = tau.ncols()
    B = Matrix(ZZ, [[0 for i in range(g)] for j in range(g)])
    for i in range(g):
        for j in range(g):
            c = - ZZ(tau[i,j].real().round())
            B[i,j] = c
            B[j,i] = c
            # TODO: maybe test whether tau is "sufficiently symmetric"
            # (tau should be symmetric, but rounding may make it non-symmetric)
    M = B_to_M(B)
    bas1, bas2 = split(bas)
    bas1 = bas1 + bas2*B
    bas_new = join(bas1, bas2)
    assert bas_new == bas * M.transpose()
    return bas_new, M

def reduce_imag_part_LLL(bas, embs, prec=None):
    """
    Given a symplectic basis and the embeddings from a CM-type,
    change the basis such that
    
      * Im(tau) is LLL-reduced

    returns a new basis and a basis transformation
    """
    g = ZZ(bas.degree() / 2)
    if prec is None:
        prec = embs[0].codomain().precision()
    tau = period_matrix(bas, embs)
    Y = Matrix([[ZZ((2^prec*t.imag()).round()) for t in u] for u in tau])
    U = Y.LLL_gram()
    # U^T * Y * U is reduced
    # So want: M = (U^T  0 )
    #              ( 0  U^-1 )
    M = ABCD_to_M(U.transpose(), zero_matrix(g, g), zero_matrix(g, g), U.inverse())
    return bas*M.transpose(), M
    
def M_for_reduce_using_z11(g):
    M = zero_matrix(2*g,2*g)
    M[0,g] = -1
    M[g,0] = 1
    for i in range(1,g):
        M[i,i] = 1
        M[g+i,g+i] = 1
    return M
    
def act_by_M_on_basis(M, bas):
    return bas*M.transpose()
