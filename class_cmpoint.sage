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
        
        
    def all_thetas(self, start_bound = 20, prec = None, bound = True):
        """
        if bound is set to true, uses theta_with_bound for each theta characteristic. if bound is set to false, uses theta_without_bound.
        """
        try:
            period_matrix = self._period_matrix
        except:
            period_matrix = self.period_matrix()
        if prec == None:
            prec = self._prec

        all_evens = [[[0,0,0],[0,0,0]],[[1,0,0],[0,0,0]],[[0,1,0],[0,0,0]],[[0,0,1],[0,0,0]],[[0,0,1],[1,0,0]],[[0,0,1],[0,1,0]],[[1,1,0],[0,0,0]],[[1,0,1],[0,0,0]],[[0,1,1],[0,0,0]],[[0,0,0],[1,0,1]],[[0,0,0],[0,1,1]],[[0,0,0],[1,1,0]],[[1,1,1],[0,0,0]],[[0,0,0],[1,1,1]],[[0,1,1],[1,0,0]],[[1,0,1],[0,1,0]],[[1,1,0],[0,0,1]],[[0,1,0],[1,0,1]],[[0,0,1],[1,1,0]],[[1,0,0],[0,1,1]],[[1,0,1],[1,0,1]],[[1,1,0],[1,1,0]],[[0,1,1],[0,1,1]],[[1,0,1],[1,1,1]],[[1,1,0],[1,1,1]],[[1,1,1],[0,1,1]],[[1,1,1],[1,0,1]],[[1,1,1],[1,1,0]],[[0,1,1],[1,1,1]],[[0,0,0],[1,0,0]],[[0,0,0],[0,1,0]],[[0,0,0],[0,0,1]],[[1,0,0],[0,1,0]],[[1,0,0],[0,0,1]],[[0,1,0],[0,0,1]],[[0,1,0],[1,0,0]]]

        all_values = []

        for even in all_evens:
            if bound == True:
                all_values.append([even, theta_with_bound(period_matrix,even[0],even[1],start_bound,prec)])
            elif bound == False:
                all_values.append([even, theta_without_bound(period_matrix,even[0],even[1],start_bound,False,prec)])
                
        self._all_thetas = all_values
        return self._all_thetas
    
    def counting(self, epsilon = 10.^(-2),bound = True):
        """
        this is meant to be used on the output of cmpoint.all_thetas() to count how many theta values are zero
        epsilon can be changed depending on what value one wants to consider to be "zero"
        """
        count = 0
        try:
            all_values = self._all_thetas
        except:
            all_values = self.all_thetas(bound = bound)
        for value in all_values:
            if value[1].abs() < epsilon:
                count += 1
                return count
        
    def vanishing_char(self, bound = True, epsilon = 10.^(-2)):
        """
        inputs:
        bound is passed to all_thetas
        epsilon can be changed depending on what value one wants to consider to be "zero"
        outputs:
        if the period matrix is plane quartic, returns None
        if the period matrix is hyperelliptic, returns the theta characteristic "delta" such that theta[delta](Z) = 0
        if there are more than one vanishing characteristics, raises an error
        """
        count = 0
        try:
            all_values = self._all_thetas
        except:
            all_values = self.all_thetas(bound = bound)
        
        for value in all_values:
            if value[1].abs() < epsilon:
                count += 1
        if count == 0:
            return None
        elif count > 1:
            raise TypeError('The entries of this period matrix are too large, the theta functions don\'t converge well')
        elif count == 1:
            for value in all_values:
                if value[1].abs() < epsilon:
                    self._vanishing_char = value[0]
                    return self._vanishing_char
                    
    def eta_dict(self, bound = True):
        """
        bound is passed to all_thetas (ultimately) True is theta_with_bound and False is theta_without_bound
        returns a dictionary giving values eta_1, eta_2, ... eta_7 for an eta-map associated to the period matrix computed for this cm point
        """
        try:
            vanishing_char = self._vanishing_char
        except:
            vanishing_char = self.vanishing_char(bound = bound)
        if vanishing_char == None:
            raise TypeError('This is a plane quartic Jacobian')
        else:
            delta = matrix(GF(2),[[vanishing_char[0][0]],[vanishing_char[0][1]],[vanishing_char[0][2]],[vanishing_char[1][0]],[vanishing_char[1][1]],[vanishing_char[1][2]]])
                    
        if delta == matrix(GF(2)):
            self._eta_dict = eta_bar
            return self._eta_dict
        else:
            delta.set_immutable()
            M = pairs[delta]
            self._eta_dict = {1: M*mumford_eta[1], 2: M*mumford_eta[2], 3: M*mumford_eta[3], 4: M*mumford_eta[4], 5:M*mumford_eta[5], 6: M*mumford_eta[6], 7:M*mumford_eta[7]}
            return self._eta_dict
            
    def U_set(self, bound = True):
        """
        returns U = {2, 4, 6} if delta is non zero and U = {1, 2, 3, 4, 5, 6, 7} if delta is zero (infinity is implicitly in both sets)
        """
        try:
            vanishing_char = self._vanishing_char
        except:
            vanishing_char = self.vanishing_char(bound = bound)
        if vanishing_char == None:
            raise TypeError('This is a plane quartic Jacobian')
        else:
            delta = matrix(GF(2),[[vanishing_char[0][0]],[vanishing_char[0][1]],[vanishing_char[0][2]],[vanishing_char[1][0]],[vanishing_char[1][1]],[vanishing_char[1][2]]])
            
        if delta == matrix(GF(2)):
            self._U_set = Set([1,2,3,4,5,6,7])
            return self._U_set
        else:
            self._U_set = Set([2,4,6])
            return self._U_set


    def one_rosenhain_coeff(self, j, prec = None, start_bound = 20, bound = True):
        """
        This function computes the Rosenhain coefficient j, with 3 <= j <= 7. We assume a1 = 0, a2 = 1. Three values are computed, for the three ways in which you can split the set in Takase's paper. This serves as an additional check that the period matrix is truly hyperelliptic.
        """
        try:
            eta_dict = self._eta_dict
        except:
            eta_dict = self.eta_dict(bound = bound)
        try:
            U = self._U_set
        except:
            U = self.U_set(bound = bound)
        if prec == None:
            prec = self._prec
        try:
            all_values = self._all_thetas
        except:
            all_values = self.all_thetas(start_bound, prec, bound)

            
        T = Set([1,2,3,4,5,6,7])
        S = Set([1,2,j])
        Y = T.difference(S)
        X = Y.subsets(2)

        ajvec = []
        # we introduce an auxiliary variable a, this only makes it so we compute the value for the decomposition given by V and W once (rather than twice, for V and W interchanged)
        if (j >= 3) and (j <= 6):
            a = j + 1
        elif j == 7:
            a = 3
        else:
            raise ValueError('j is not between 3 and 7 inclusively')

        for V in X:
            if a in V:
                W = Y.difference(V)
                setA = U.symmetric_difference(V.union(Set([1,2])))
                setB = U.symmetric_difference(W.union(Set([1,2])))
                setC = U.symmetric_difference(V.union(Set([1,j])))
                setD = U.symmetric_difference(W.union(Set([1,j])))
                A = theta_from_char_and_list(all_values, compute_characteristic_sum_from_set_and_etas(setA,eta_dict))
                B = theta_from_char_and_list(all_values, compute_characteristic_sum_from_set_and_etas(setB,eta_dict))
                C = theta_from_char_and_list(all_values, compute_characteristic_sum_from_set_and_etas(setC,eta_dict))
                D = theta_from_char_and_list(all_values, compute_characteristic_sum_from_set_and_etas(setD,eta_dict))
                aj = ((A*B)^2)/((C*D)^2)
                ajvec.append(aj)
        return ajvec
        
    def all_rosenhain_coeffs(self, prec = None, start_bound = 20, bound = False):
        """
        returns all of the vectors of Rosenhain coefficients
        """
        if prec == None:
            prec = self._prec
            
        all_coeffs = []
        for j in range(3,8):
            all_coeffs.append(self.one_rosenhain_coeff(j,prec, start_bound, bound))
        return all_coeffs

