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
        
    def symplectic_basis(self):
        K = self._K
        ideal = self._ideal
        xi = self._xi
        prec = self._prec 
        
        basis = ideal.basis()
        riemann_form = Matrix(ZZ,[[(x.conjugate()*xi*y).trace() for y in basis] for x in basis])
        symplectic_basis = Sequence(riemann_form.symplectic_form()[1]*vector(basis))
        return symplectic_basis       
        
    def acc_period_matrix(self, test = False):
        """
        This computes the period matrix several times, refining the precision of the CM type until two consecutive period matrices have entries that agree up to the precision of the CM point
        """
        K = self._K
        ideal = self._ideal
        xi = self._xi
        prec = self._prec
        from sage.rings.number_field.number_field import refine_embedding
        
        basis = ideal.basis()
        riemann_form = Matrix(ZZ,[[(x.conjugate()*xi*y).trace() for y in basis] for x in basis])
        symplectic_basis = Sequence(riemann_form.symplectic_form()[1]*vector(basis))
        phis = self._CM_type
        big_period_matrix = Matrix([[phi(b) for b in symplectic_basis] for phi in phis])
        big_period_matrix.subdivide(3,3)
        Omega1 = big_period_matrix.subdivision(0,0)
        Omega2 = big_period_matrix.subdivision(0,1)
        Zs = [Omega2.adjoint()*Omega1/Omega2.det()]
        
        equality = False
        iterates = 0        
        
        while equality == False:
            iterates += 1
            phis = [refine_embedding(phi,prec + iterates*20) for phi in phis]
            big_period_matrix = Matrix([[phi(b) for b in symplectic_basis] for phi in phis])
            big_period_matrix.subdivide(3,3)
            Omega1 = big_period_matrix.subdivision(0,0)
            Omega2 = big_period_matrix.subdivision(0,1)
            Zs.append(Omega2.adjoint()*Omega1/Omega2.det())
        
            if test:
                print "computing iteration number {0}".format(iterates) 
                print Zs[iterates]
                CMvalimag=matrix(RR,3,3)
                for i in range(3):
                    for j in range(3):
                        CMvalimag[i,j] = Zs[iterates][i,j].imag()
                print CMvalimag.eigenvalues()
            
            if all([compare(Zs[iterates][i,j],Zs[iterates-1][i,j],prec+10) for i in range(3) for j in range(3)]):
                equality = True 
        
        Z = Zs[iterates]
        
        CC = ComplexField(prec)
        A= CC(Z[0][0])
        B= CC(Z[0][1])
        C= CC(Z[0][2])
        D= CC(Z[1][1])
        E= CC(Z[1][2])
        F= CC(Z[2][2])

        M =[A, B, C, D, E, F]
        self._period_matrix = M
        return self._period_matrix

    def reciprocity_map_image(self, level):
        """
        Given a CM period matrix Z, returns Z.epsilon(b) mod level
        for a set of generators A of the group (H(1) cap I(modulus)) / H(modulus).

        Here I(N) is the group of ideals of the reflex field coprime to N,
        and H(N) is the subgroup of ideals A such that the reflex type norm of A
        is principal and generated by mu with mu*mubar in QQ.
    
        The b in Z.epsilon(b) is b=mu with mu as in the definition of H(1).
    
        If modulus is None, then take modulus=level
    
        EXAMPLE::
    
            sage: load("recip.sage")
            sage: k = CM_Field((x^2+5)^2-4*5)
            sage: Z = k.one_period_matrix(); Z
            Period Matrix
            [0.6909830056251? + 0.9510565162952?*I 1.1909830056251? + 0.5877852522925?*I]
            [1.1909830056251? + 0.5877852522925?*I          0.?e... + 1.1755705045850?*I]
            sage: reciprocity_map_image(Z, 6) # not tested, is this answer correct?
            [[1 1 4 1]
            [2 4 1 5]
            [2 2 0 3]
            [1 2 4 2], [3 2 2 0]
            [2 3 0 2]
            [4 0 3 0]
            [2 4 0 5]]
            sage: gammas = reciprocity_map_image(Z, 8)
            sage: len(gammas)
            4
            sage: gammas[0] # not tested, is this answer correct?
            [0 2 1 6]
            [0 1 6 7]
            [7 1 2 1]
            [2 7 7 2]
        """
        Phi = self._CM_type
        K = self._K
        Psi = K.reflex(Phi)
        
        elts = K.principal_type_norms(Psi, level)
        
        gammas = [GSp_element(mat_convert(self.epsilon(b), Zmod(level))) for b in elts]
        gammas = [g for g in gammas if not g == identity_matrix(6)]
        return gammas
        
    def Bt(self):
        return Matrix([b.vector() for b in self.symplectic_basis()])
        
    def complex_conjugation_symplectic_matrix(self,level,mu,A):
        raise NotImplementedError
    
    def epsilon(self, x):
        """
        The map epsilon of page 57 of Shimura's "on certain reciprocity laws..."
        Returns the transpose of the matrix of multiplication by x wrt the basis self.basis()
        """
        Mt = Matrix([(x*b).vector() for b in self.symplectic_basis()])
        return Mt*(self.Bt().inverse())
        
    def all_thetas(self, start_bound = 20, prec = None, bound = False):
        """
        if bound is set to true, uses theta_with_bound for each theta characteristic. if bound is set to false, uses theta_without_bound.
        """
        try:
            period_matrix = self._period_matrix
        except:
            period_matrix = self.acc_period_matrix()
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
    
    def counting(self, epsilon = 10.^(-2),bound = False):
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
        
    def vanishing_char(self, bound = False, epsilon = 10.^(-2)):
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
                    
    def eta_dict(self, bound = False, epsilon = 10.^(-2)):
        """
        bound is passed to all_thetas (ultimately) True is theta_with_bound and False is theta_without_bound
        returns a dictionary giving values eta_1, eta_2, ... eta_7 for an eta-map associated to the period matrix computed for this cm point
        """
        try:
            vanishing_char = self._vanishing_char
        except:
            vanishing_char = self.vanishing_char(bound = bound, epsilon = epsilon)
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
            
    def U_set(self, bound = False, epsilon = 10.^(-2)):
        """
        returns U = {2, 4, 6} if delta is non zero and U = {1, 2, 3, 4, 5, 6, 7} if delta is zero (infinity is implicitly in both sets)
        """
        try:
            vanishing_char = self._vanishing_char
        except:
            vanishing_char = self.vanishing_char(bound = bound, epsilon = epsilon)
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


    def one_rosenhain_coeff(self, j, prec = None, start_bound = 20, bound = False, epsilon = 10.^(-2)):
        """
        This function computes the Rosenhain coefficient j, with 3 <= j <= 7. We assume a1 = 0, a2 = 1. Three values are computed, for the three ways in which you can split the set in Takase's paper. This serves as an additional check that the period matrix is truly hyperelliptic.
        """
        try:
            eta_dict = self._eta_dict
        except:
            eta_dict = self.eta_dict(bound = bound, epsilon = epsilon)
        try:
            U = self._U_set
        except:
            U = self.U_set(bound = bound, epsilon = epsilon)
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
        
    def all_rosenhain_coeffs(self, prec = None, start_bound = 20, bound = False, epsilon = 10.^(-2)):
        """
        returns all of the vectors of Rosenhain coefficients
        """
        if prec == None:
            prec = self._prec
            
        all_coeffs = []
        for j in range(3,8):
            all_coeffs.append(self.one_rosenhain_coeff(j,prec, start_bound, bound, epsilon))
        return all_coeffs
        
    def all_four_chars(self, j, prec = None, start_bound = 20, bound = False, epsilon = 10.^(-2)):
        try:
            all_values = self._all_thetas
        except:
            all_values = self.all_thetas(start_bound, prec, bound)
        try:
            eta_dict = self._eta_dict
        except:
            eta_dict = self.eta_dict(bound = bound, epsilon = epsilon)
        try:
            U = self._U_set
        except:
            U = self.U_set(bound = bound, epsilon = epsilon)
        if prec == None:
            prec = self._prec
            
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
                charA = compute_characteristic_sum_from_set_and_etas(setA,eta_dict)
                charB = compute_characteristic_sum_from_set_and_etas(setB,eta_dict)
                charC = compute_characteristic_sum_from_set_and_etas(setC,eta_dict)
                charD = compute_characteristic_sum_from_set_and_etas(setD,eta_dict)
                ajvec.append([charA,charB,charC,charD])
    
        return ajvec
    
    def four_chars(self, j, prec = None, start_bound = 20, bound = False, epsilon = 10.^(-2)):
        try:
            all_values = self._all_thetas
        except:
            all_values = self.all_thetas(start_bound, prec, bound)
        try:
            eta_dict = self._eta_dict
        except:
            eta_dict = self.eta_dict(bound = bound, epsilon = epsilon)
        try:
            U = self._U_set
        except:
            U = self.U_set(bound = bound, epsilon = epsilon)
        if prec == None:
            prec = self._prec
            
        T = Set([1,2,3,4,5,6,7])
        S = Set([1,2,j])
        Y = T.difference(S)
        X = Y.subsets(2)
    
        V = X[0]
        W = Y.difference(V)
        setA = U.symmetric_difference(V.union(Set([1,2])))
        setB = U.symmetric_difference(W.union(Set([1,2])))
        setC = U.symmetric_difference(V.union(Set([1,j])))
        setD = U.symmetric_difference(W.union(Set([1,j])))
        charA = compute_characteristic_sum_from_set_and_etas(setA,eta_dict)
        charB = compute_characteristic_sum_from_set_and_etas(setB,eta_dict)
        charC = compute_characteristic_sum_from_set_and_etas(setC,eta_dict)
        charD = compute_characteristic_sum_from_set_and_etas(setD,eta_dict)
    
        return [charA,charB,charC,charD]    
                
    def conjugate_chars(self, good_H, ros_char):
        conjugate_chars = []
        for pair in good_H:
            action = theta_action_on_rosenhain(good_ros_chars,pair)
            conjugate_chars = add_to_list(conjugate_chars,action)
        return conjugate_chars
        
    def lambda_conj(self, chars_factor, prec = None, start_bound = 20, bound = False):
        try:
            all_values = self._all_thetas
        except:
            all_values = self.all_thetas(start_bound, prec, bound)
            
        A = theta_from_char_and_list(all_values, c_to_char(chars_factor[0]))
        B = theta_from_char_and_list(all_values, c_to_char(chars_factor[1]))
        C = theta_from_char_and_list(all_values, c_to_char(chars_factor[2]))
        D = theta_from_char_and_list(all_values, c_to_char(chars_factor[3]))
        
        if A == None or B == None or C == None or D == None:
            return None
            
        aj = chars_factor[4]*((A*B)^2)/((C*D)^2)
        return aj    
        
    def all_lambda_conjs(self, good_H, ros_chars, prec = None, start_bound = 20, bound = False):
        try:
            all_values = self._all_thetas
        except:
            all_values = self.all_thetas(start_bound, prec, bound)
            
        good_ros_chars = [char_to_c(char) for char in ros_chars]
        
        conjugate_chars = []
        for pair in good_H:
            action = theta_action_on_rosenhain(good_ros_chars,pair)
            conjugate_chars = add_to_list(conjugate_chars,action)
            
        conjugate_lambda = []
        for c in conjugate_chars:
            conjugate  = lambda_conj(all_values,c)
            if conjugate == None:
                pass
            else:
                conjugate_lambda.append(conjugate)
                
        return conjugate_lambda


def mat_convert(M, ring_or_map):
    """
    Applies `ring_or_map` to the coefficients of M,
    i.e. given M = (m_ij)_ij, returns
    (ring_or_map(m_ij))_ij
    
    EXAMPLES::
    
        sage: load("recip.sage")
        sage: M = Matrix([[1/2, 1/5], [7, 9/2]])
        sage: mat_convert(M, floor)
        [0 0]
        [7 4]
        sage: mat_convert(M, GF(3))
        [2 2]
        [1 0]
    """
    return Matrix([[ring_or_map(M[i,j]) for j in range(M.ncols())] \
                    for i in range(M.nrows())])
                    
def char_to_c(char):
    c = [QQ(i)*1/2 for i in char[0]] + [QQ(i)*1/2 for i in char[1]]
    return c
    
def theta_action_on_rosenhain(good_chars,pair):
    almost_ros = [theta_action_without_kappa(pair[0],pair[1],good_chars[i]) for i in range(4)]
    factor = prod([exp(2*pi*I*ros[1]) for ros in almost_ros])
    return [ros[0] for ros in almost_ros] + [factor]
    
def add_to_list(conjugates,new_elt):
    print 'new element' #takeout
    print new_elt #takeout
    if new_elt in conjugates:
        print 'already in list' #takeout
        return conjugates
    else:
        print 'not in list' #takeout
        conjugates.append(new_elt)
        return conjugates
        
def c_to_char(c):
    char = [[ZZ(c[i]*2) for i in range(3)],[ZZ(c[i]*2) for i in range(3,6)]]
    return char
    
