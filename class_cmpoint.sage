class CMPoint:

    def __init__(self, K, CM_type, ideal, xi, prec=None):
        """
        K: a CMField (see class CMFieldfromPoly)
        CM_type: can be generated from the CMField, see documentation 
        ideal: can be generated from the CMField, see documentation
        xi: can be generated from the CMField, see documentation
        prec : bits of precision for complex-valued calculations, if None, inherited from CMField
        """
        self._K = K
        self._CM_type = CM_type
        self._ideal = ideal
        self._xi = xi
        if prec == None:
            self._prec = self._K._prec
        else:
            self._prec = prec
            
    def __repr__(self):
        return "CM point for CM field %s \n given by the ideal %s \n and the CM-type %s,\n with polarization given by %s"%(self._K, self._ideal, self._CM_type, self._xi)

    def cubic_fld(self):
        """
        Gives CMpoint its cubic field K0, which is a subfield of the optimal sextic
        """
        self._K0 = self._K._K0
        return self._K0

    def sextic_fld(self):
        """
        Gives CMpoint its sextic field
        """
        return self._K

    def period_matrix(self, test = False):
        """
        Set test to True to see the period matrix and the eigenvalues; can be used to check if it looks symmetric and positive definite
        """
        K = self._K
        phi0, phi1, phi2 = self._CM_type
        ideal = self._ideal
        xi = self._xi

        basis = ideal.basis()
        riemann_form = Matrix(ZZ,[[(conjugate(x)*xi*y).trace() for y in basis] for x in basis])
        symplectic_basis = Sequence(riemann_form.symplectic_form()[1]*vector(basis))
        big_period_matrix = Matrix([[phi(b) for b in symplectic_basis] for phi in self._CM_type])
        big_period_matrix.subdivide(3,3)
        Omega1 = big_period_matrix.subdivision(0,0)
        Omega2 = big_period_matrix.subdivision(0,1)
        Z = Omega2.adjoint()*Omega1/Omega2.det()
        if test:
            print Z
            
        if test:
            CMvalimag=matrix(RR,3,3)
            for i in range(3):
                 for j in range(3):
                        CMvalimag[i,j]=Z[i,j].imag()
            print CMvalimag.eigenvalues()

        A= Z[0][0]
        B= Z[0][1]
        C= Z[0][2]
        D= Z[1][1]
        E= Z[1][2]
        F= Z[2][2]

        M =[A, B, C, D, E, F]
        self._period_matrix = M
        return self._period_matrix
