"""

This is the latest version of the code computing the Rosenhain invariants, with the sign corrected.

"""


load('reduced_period_matrix.sage')

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

    def old_period_matrix(self, test = False):
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
        riemann_form = Matrix(ZZ,[[(conjugate(x)*xi*y).trace() for y in basis] for x in basis])
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
        
        #print Z
	
        return self._period_matrix, Omega2


    def reduced_period_matrix(self, prec=None):
        """
        This method was written by Kilicer and Streng (see reduced_period_matrix.sage for full copyright statement)
        """

        K = self._K   
        Phi=self._CM_type
        ideal = self._ideal
        xi = self._xi
        prec = self._prec        
        
        basis = ideal.basis()
        riemann_form = Matrix(ZZ,[[(conjugate(x)*xi*y).trace() for y in basis] for x in basis])
        symplectic_basis = Sequence(riemann_form.symplectic_form()[1]*vector(basis))
        bas=vector(symplectic_basis)
        (redbas, M) = reduce_Siegel_epsilon(bas, Phi, False)        
        bigmatrix = Matrix([[phi(b) for b in redbas] for phi in Phi])
        Omega1 = bigmatrix[:,:3]
        Omega2 = bigmatrix[:,3:]
        Z = Omega2.inverse()*Omega1
        CC = ComplexField(prec)
        A= CC(Z[0][0])
        B= CC(Z[0][1])
        C= CC(Z[0][2])
        D= CC(Z[1][1])
        E= CC(Z[1][2])
        F= CC(Z[2][2])

        M =[A, B, C, D, E, F]
        self._period_matrix = M
        return self._period_matrix, Omega2
       
        
    def power(a,n):
       '''calcule la n-ieme puissance de a mod q par square and multiply, en recursif'''
       if n==0: return 1
       else:
         k=n/2 
         b=power(a, k)
       if (n%2==0):
          return b*b 
       else:
          return ((b*b)*a)

    def sigma_140(self,prec=None, start_bound=20, bound=False, epsilon=10.^(-2)):

        try:
            all_values = self._all_thetas
        except:
            all_values = self.all_thetas(bound = bound)
	all_powered_values=[]    
        for value in all_values: 
           all_powered_values.append(power(value[1],8))
        sigma=0
        for i in range(36):
              prod=1
              for j in range(36):
                   if (i!=j):
                       prod=prod*all_powered_values[j]
              sigma=sigma+prod
	return sigma


    def h4(self,prec=None, start_bound=20, bound=False, epsilon=10.^(-2)):
        try:
            all_values = self._all_thetas
        except:
            all_values = self.all_thetas(bound = bound)
	    sigma=0    
        for value in all_values:
            sigma=sigma+power(value[1],8)
        
        return sigma
    
    def all_thetas(self, start_bound = 20, prec = None, bound = False):
        """
        if bound is set to true, uses theta_with_bound for each theta characteristic. if bound is set to false, uses theta_without_bound.
        """
        try:
            period_matrix = self._period_matrix
        except:
            #period_matrix = self.acc_period_matrix()[0]
            period_matrix = self.reduced_period_matrix()[0]
        if prec == None:
            prec = self._prec

        all_evens = [[[0,0,0],[0,0,0]],[[1,0,0],[0,0,0]],[[0,1,0],[0,0,0]],[[0,0,1],[0,0,0]],[[0,0,1],[1,0,0]],[[0,0,1],[0,1,0]],[[1,1,0],[0,0,0]],[[1,0,1],[0,0,0]],[[0,1,1],[0,0,0]],[[0,0,0],[1,0,1]],[[0,0,0],[0,1,1]],[[0,0,0],[1,1,0]],[[1,1,1],[0,0,0]],[[0,0,0],[1,1,1]],[[0,1,1],[1,0,0]],[[1,0,1],[0,1,0]],[[1,1,0],[0,0,1]],[[0,1,0],[1,0,1]],[[0,0,1],[1,1,0]],[[1,0,0],[0,1,1]],[[1,0,1],[1,0,1]],[[1,1,0],[1,1,0]],[[0,1,1],[0,1,1]],[[1,0,1],[1,1,1]],[[1,1,0],[1,1,1]],[[1,1,1],[0,1,1]],[[1,1,1],[1,0,1]],[[1,1,1],[1,1,0]],[[0,1,1],[1,1,1]],[[0,0,0],[1,0,0]],[[0,0,0],[0,1,0]],[[0,0,0],[0,0,1]],[[1,0,0],[0,1,0]],[[1,0,0],[0,0,1]],[[0,1,0],[0,0,1]],[[0,1,0],[1,0,0]]]
        
        all_values = [[even,theta_function(period_matrix,even[0],even[1],prec)] for even in all_evens]
        #all_values = [[[val[0][0][1],val[0][0][2]],val[1]] for val in sorted(list(theta_function([(period_matrix,even[0],even[1],prec) for even in all_evens])))]
                
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

    def vareps(self, j, eta_dict):
        v1 = eta_dict[j]-eta_dict[2]
        v2 = eta_dict[1]
        J = matrix(GF(2),[[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]])
        value = v1.transpose()*J*v2
        if value == 0:
            return 1
        elif value == 1:
            return -1


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
        #X = Y.subsets(2)
        jk = Y.subsets(2)[0]
        D1 = jk.union(Set([j]))
        D2 = jk.union(Set([2]))
        
        #ajvec = []
        # we introduce an auxiliary variable a, this only makes it so we compute the value for the decomposition given by V and W once (rather than twice, for V and W interchanged)
        #if (j >= 3) and (j <= 6):
        #    a = j + 1
        #elif j == 7:
        #    a = 3
        #else:
        #    raise ValueError('j is not between 3 and 7 inclusively')

        Cec=ComplexField(prec)
        #for V in X:
        #    if a in V:
        #        W = Y.difference(V)
        #        setA = U.symmetric_difference(V.union(Set([1,j]))) #was 1,2
        #        setB = U.symmetric_difference(W.union(Set([1,j]))) #was 1,2
        #        setC = U.symmetric_difference(V.union(Set([1,2]))) #was 1,j
        #        setD = U.symmetric_difference(W.union(Set([1,2]))) #was 1,j
        setA = U.symmetric_difference(D1.union(Set([1])))
        setB = U.symmetric_difference(D2)
        setC = U.symmetric_difference(D1)
        setD = U.symmetric_difference(D2.union(Set([1])))
        A = Cec(theta_from_char_and_list(all_values, compute_characteristic_sum_from_set_and_etas(setA,eta_dict)))
        B = Cec(theta_from_char_and_list(all_values, compute_characteristic_sum_from_set_and_etas(setB,eta_dict)))
        C = Cec(theta_from_char_and_list(all_values, compute_characteristic_sum_from_set_and_etas(setC,eta_dict)))
        D = Cec(theta_from_char_and_list(all_values, compute_characteristic_sum_from_set_and_etas(setD,eta_dict)))
        aj = Cec(self.vareps(j,eta_dict)*((A*B)^2)/((C*D)^2))
                #ajvec.append(aj)
        #return ajvec
        return aj
        
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




    

    
    
