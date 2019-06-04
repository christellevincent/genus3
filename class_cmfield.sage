from sage.rings.number_field.number_field import NumberField_absolute
from sage.numerical.mip import MIPSolverException

class CMFieldfromPoly(NumberField_absolute):
    
    def __init__(self, polynomial, names, prec = 664):
        NumberField_absolute.__init__(self, polynomial, names)
        self._prec = prec
        self._K0 = self.subfields(3)[0][0]
        self._embedK0toK = self.subfields(3)[0][1]
        self._CMPoints = [] #here we can cache the CM points of K
        
    def K0(self):
        return self._K0
    
    def embedK0toK(self):
        return self._embedK0toK
        
    def one_CMtype(self, prec=None):
        """
        Given a CMField K, returns a PRIMITIVE CM type for K, as a list of embeddings of K into the complex numbers.
        """
        if prec == None:
            prec = self._prec
        embeddings = []
        for embed in self.complex_embeddings(prec):
            embeddings.append(embed)
        phi0 = embeddings.pop(0)
        phi0bar = find_conj(embeddings,phi0)
        embeddings.pop(embeddings.index(phi0bar))
        phi1= embeddings.pop(0)
        phi1bar = find_conj(embeddings,phi1)
        embeddings.pop(embeddings.index(phi1bar))
        phi2 = embeddings.pop(0)
        candidate = [phi0,phi1,phi2]
        if is_primitive_CMtype(candidate, prec):
            return candidate
        else:
            phi2bar = find_conj(embeddings,phi2)
            candidate = [phi0,phi1,phi2bar]
            assert is_primitive_CMtype(candidate, prec)
            return candidate

    def all_CMtypes(self, prec = None):
        """
        Given a CM sextic field K, returns a list of all PRIMITIVE CM types for K
        """
        if prec == None:
            prec = self._prec
        primitives = []
        embeddings = set(self.complex_embeddings(prec))
        from itertools import combinations
        all_triples = list(map(list, combinations(embeddings, 3)))
        for triple in all_triples:
            try:
               if is_primitive_CMtype(triple, prec):
                 primitives.append(triple)
            except TypeError:
                continue
        self._allCMtypes = primitives
        return self._allCMtypes
    
    def generators_U_plus(self, case = False):
        """
        Returns a list of free generators of the group U^+, which is the group of totally positive units in K0, the cubic subfield of K
        if case = True, returns the case (1, 2a, 2b, 3 or 4, see Koike-Weng) for use in the units_epsilon algorithm
        """
        K0 = self._K0
        embedK0toK = self._embedK0toK
        u1,u2 = UnitGroup(K0).fundamental_units()
        # we make sure that each unit is positive under the first embedding (to ensure it's not totally negative)
        if K0.complex_embeddings()[0](u1) < 0:
            u1 = -u1
        if K0.complex_embeddings()[0](u2) < 0:
            u2 = -u2
        #no precision is given here for the embeddings since we only care if they are positive or negative numbers
        u1values = [embed(u1) > 0 for embed in K0.complex_embeddings()]
        u2values = [embed(u2) > 0 for embed in K0.complex_embeddings()]
        generators = []
        case_number = ''
        if (all(u1values) and all(u2values)): #both fundamental units are totally positive
            generators = [u1, u2]
            case_number = '1'
        elif all(u1values): #the first fundamental unit is totally positive
            generators = [u1, u2^2]
            case_number = '2a'
        elif all(u2values): #the second fundamental unit is totally positive
            generators = [u1^2,u2]
            case_number = '2b'
        elif all([u1values[i] == u2values[i] for i in range(3)]): #both fundamental units have the same embedding type
            generators = [u1*u2,u2^2]
            case_number = '3'
        else:
            generators = [u1^2,u2^2]
            case_number = '4'
        self._generators_U_plus = generators
        if case == False:
            return self._generators_U_plus
        elif case == True:
            return [self._generators_U_plus,case_number]

    def unit_maker(self, tietze_list):
        """
        In units_for_xi we use a free group to compute a quotient of unit groups. This converts an element of the free group in Tietze form to an element of the unit group
        """
        unit = self.one()
        for el in tietze_list:
            if el == 1:
                unit *= u0
            elif el == -1:
                unit *= u0**(-1)
            elif el == 2:
                unit *= u1
            elif el == -2:
                unit *= u1**(-1)
            elif el == 3:
                unit *= u2
            elif el == -3:
                unit *= u2**(-1)
        return unit

    def units_for_xi(self):
        """
        Returns a list of representatives of the group U_K/U^+, where U_K are all units in K, and U^+ are all totally positive units in K0.
        """
        K0 = self._K0
        embedK0toK = self._embedK0toK
        UK = self.unit_group()
        try:
            generators_U_plus = self._generators_U_plus
        except:
            generators_U_plus = self.generators_U_plus()
        logs = [UK.log(embedK0toK(u)) for u in generators_U_plus] #write generators of U^+ in terms of powers of generators of U_K
        #we now use a free group to compute the quotient group
        F.<a,b,c> = FreeGroup()
        commutators = [F([1, 2, -1, -2]),F([1,3,-1,-3]),F([2,3,-2,-3])] #unit group is commutative
        roots_of_unity = [a^UK.zeta_order()] #first generator is a root of unity
        pos_units = [prod((u**e for u,e in zip([a,b,c],logs[0])), F(1) ),prod((u**e for u,e in zip([a,b,c],logs[1])), F(1) )] #write generators of U^+ into the free group
        relations = commutators + pos_units + roots_of_unity
        G = F/relations #quotient group
        elements = [el.Tietze() for el in G.list()] #output elements of quotient group in readable form
        UK.inject_variables(verbose=False)
        unit_list = [self.unit_maker(element) for element in elements]
        self._units_for_xi = unit_list
        return self._units_for_xi
    
    def in_U_1(self,u):
        """
        u is an element of self.unit_group()
        returns True if u is in U_1, the elements of the form epsilon * bar(epsilon) for epsilon a unit in K
        returns False otherwise
        """
        UK = self.unit_group()
        v0, v1, v2 = UK.gens()
        conj = self.complex_conjugation()
        e10,e11,e12 = UK.log(conj(v1))
        e20,e21,e22 = UK.log(conj(v2))
        a0,a1,a2 = UK.log(u)
        try:
            p = MixedIntegerLinearProgram(maximization=False, solver = "GLPK")
            w = p.new_variable(integer=True, nonnegative=True)
            p.add_constraint((1+e11)*w[0] + e21*w[1] == a1)
            p.add_constraint(e12*w[0] + (1+e22)*w[1] == a2)
            p.set_min(w[0], None)
            p.set_min(w[1], None)
            p.set_objective(None)
            p.solve()
            n1,n2 = [int(v[1]) for v in p.get_values(w).iteritems()]
            assert n1*e10 + n2*e20 == a0
            return True
        except MIPSolverException:
            return False

    def units_epsilon(self):
        """
        returns a list of representatives of the group U^+/U_1, where U^+ are all totally positive units in K and U_1 are all elements of the form epsilon * bar(epsilon) for epsilon a unit in K
        Based on the algorithm in Koike-Weng
        """
        generators, case = self.generators_U_plus(case = True)
        embedK0toK = self._embedK0toK
        if case == '4':
            return [self.one()]
        elif case == '2a':
            if self.in_U_1(embedK0toK(generators[0])):
                return [self.one()]
            else:
                return [self.one(), embedK0toK(generators[0])]
        elif case == '2b':
            if self.in_U_1(embedK0toK(generators[1])):
                return [self.one()]
            else:
                return [self.one(), embedK0toK(generators[1])]
        elif case == '3':
            if self.in_U_1(embedK0toK(generators[0])):
                return [self.one()]
            else:
                return [self.one(), embedK0toK(generators[0])]
        elif case == '1':
            bool1 = self.in_U_1(embedK0toK(generators[0]))
            bool2 = self.in_U_1(embedK0toK(generators[1]))
            if bool1 and bool2:
                return [self.one()]
            elif bool1:
                return [self.one(), embedK0toK(generators[1])]
            elif bool2:
                return [self.one(), embedK0toK(generators[0])]
            else:
                if self.in_U_1(embedK0toK(generators[0]*generators[1])):
                    return [self.one(), embedK0toK(generators[0])]
                else:
                    return [self.one(), embedK0toK(generators[0]), embedK0toK(generators[1]),embedK0toK(generators[0]*generators[1])]

    def suitable_ideals(self):
        """
        Given a CMField K, gives a list of representatives of elements of the ideal class group such that different_K * ideal * bar(ideal) is principal, where bar is the complex conjugation on K.
        For each of these ideal classes, there could be a primitive CM type Phi such that CC^3/Phi(ideal) is a principally polarized abelian variety
        """
        ideals = [J.ideal() for J in list(self.class_group())]
        conj = self.complex_conjugation()
        suitables = []
        for ideal in ideals:
            princ_ideal = self.different() * ideal * conj(ideal)
            if princ_ideal.is_principal():
                suitables.append(ideal)
        self._suitable_ideals = suitables
        return self._suitable_ideals
    
    def good_generator(self,CM_type,princ_ideal):
        """
        Given a CMField K, a primitive CM type and PRINCIPAL ideal, returns b, a generator of the ideal that is totally imaginary and that has imaginary part negative under each embedding in the CM type
        IMPORTANT: This function might fail to produce a good generator even if one exists. The generator delta might be totally imaginary, but embed into CC with a very small real part because of rounding errors. If one of the units is very large, this could make the real part of u*delta too large and fail to return this generator. One way to catch this mistake is to run the computation for all CM types in a given equivalence class.
        """
        from sage.rings.number_field.number_field import refine_embedding
        prec = self._prec
        
        delta = princ_ideal.gens_reduced()[0]
        
        try:
            units_for_xi = self._units_for_xi
        except:
            units_for_xi = self.units_for_xi()
        for u in units_for_xi:
            #we refine the embedding to get more precision and hopefully keep the real part small
            new_CM_type = [refine_embedding(phi,2*prec) for phi in CM_type]
            if all([new_CM_type[i](u*delta).imag()<0 for i in range(3)]) and all([compare(new_CM_type[i](u*delta).real(),0.0,prec) for i in range(3)]):
                return u*delta
        return None
    
    def princ_polarized(self,CM_type):
        """
        Given a CMField K and a CM-type Phi, returns a list of all pairs (ideal, xi) such that CC^3/Phi(ideal) is a principally polarized abelian variety with p.p. given by xi.
        IMPORTANT: This might miss some pairs because of a rounding error in good_generator (see note there for more details)
        """
        all_pairs = []
        conj = self.complex_conjugation()
        try:
            suitable_ideals = self._suitable_ideals
        except:
            suitable_ideals = self.suitable_ideals()
        for ideal in suitable_ideals:
            princ_ideal = self.different() * ideal * conj(ideal)
            b = self.good_generator(CM_type,princ_ideal)
            if not(b == None):
                xi = b**(-1)
                for u in self.units_epsilon():
                    all_pairs.append([ideal,u*xi])
        return all_pairs
    
def CMField(cubic,triple,prec=664):
    """
    creates a CMField from a cubic and a triple
    """
    return CMFieldfromPoly(polynomial_from_cubic_triple(cubic,triple),'a',prec)

def polynomial_from_cubic_triple(cubic,triple):
    """
    creates the optimal polynomial for a CM sextic field from a cubic and a triple
    """
    k0.<d1> = NumberField(cubic)
    a0 = k0.gen()
    S.<x> = k0['x']
    A,B,C = triple
    fpol = x**2+A+B*a0+C*a0**2
    K.<d2> = NumberField(fpol)
    R_abs = K.absolute_polynomial() 
    s = 'newf = polredabs(%s);' %R_abs
    gp(s)
    newpoly = gp.eval('newf')
    S.<x> = PolynomialRing(QQ,'x')
    return S(newpoly)

def compare(num1, num2, prec = None):
    """
    Because we embedd everything into CC, there are losses of precision. We need to check if two complex numbers are "close enough"
    """
    if prec == None:
        prec1 = num1.prec()
        prec2 = num2.prec()
        prec = min(prec1, prec2)/2 #this is ad hoc
    B=2^(-prec)
    if (num1.real() - num2.real()).abs() < B and (num1.imag() - num2.imag()).abs() < B:
        return True
    else:
        return False
    
def find_conj(embeddings,phi):
    """
    This finds the conjugate of a complex embedding within a list of complex embeddings
    """
    for embed in embeddings:
        if embed.im_gens()[0].conjugate() == phi.im_gens()[0]:
            return embed

def is_CMtype(phis):
    """
    Given a list of three complex embeddings of a CM sextic field, checks if they form a CM type
    """
    if (phis[0].im_gens()[0].conjugate()==phis[1].im_gens()[0] or phis[1].im_gens()[0].conjugate()==phis[2].im_gens()[0] or phis[2].im_gens()[0].conjugate()==phis[0].im_gens()[0]):
        return False
    else:
        return True
        
def is_primitive_galois(phis):
    """
    Given a list of three complex embeddings of a CM sextic field with Galois group ZZ/6, checks if it is a primitive CM type
    """
    if is_CMtype(phis) == False:
        raise TypeError('This is not even a CM type')
    K = phis[0].domain()
    c = K.gen()
    phi0,phi1,phi2 = phis
    G = K.galois_group()
    #H is a group of order 3. if the CM-type "is" this group H, then it is not primitive. in turn, we set each element of the CM-type to be the identity and see if the other two elements can be the other two elements, non-trivial, elements of H. See Weng for proof
    H = []
    for g in G:
        if g.order() == 3:
            H.append(g)
    if compare(phi0(c), phi1(H[0](c))) and compare(phi0(c), phi2(H[1](c))):
        return False
    elif compare(phi0(c), phi2(H[0](c))) and compare(phi0(c), phi1(H[1](c))):
        return False
    elif compare(phi1(c), phi0(H[0](c))) and compare(phi1(c), phi2(H[1](c))):
        return False
    elif compare(phi1(c), phi2(H[0](c))) and compare(phi1(c), phi0(H[1](c))):
        return False
    elif compare(phi2(c), phi0(H[0](c))) and compare(phi2(c), phi1(H[1](c))):
        return False
    elif compare(phi2(c), phi1(H[0](c))) and compare(phi2(c), phi0(H[1](c))):
        return False
    else:
        return True

def is_H(psi0,psi1,psi2,rho,c,H):
    """
    Given an ordered triple psi0, psi1, psi2, a generator c of the field K, an embedding rho of the Galois closure of K into CC, and a subgroup H of the Galois group of K, checks if psi0 = rho H[0], psi1 = rho H[1] and psi2 = rho H[2]
    """
    if compare(rho(H[0](c)),psi0(c)) and compare(rho(H[1](c)),psi1(c)) and compare(rho(H[2](c)),psi2(c)):
        return True
    elif compare(rho(H[0](c)).conjugate(),psi0(c)) and compare(rho(H[1](c)).conjugate(),psi1(c)) and compare(rho(H[2](c)).conjugate(),psi2(c)):
        return True
    else:
        return False 

def is_primitive_nongalois(phis, prec):
    """
    Given a list of three complex embeddings of a CM sextic field with Galois group ZZ/2 * S3, checks if it is a primitive CM type
    """
    if is_CMtype(phis) == False:
        raise TypeError('This is not even a CM type')
    K = phis[0].domain()
    c = K.gen()
    phi0,phi1,phi2 = phis
    L.<t> = K.galois_closure()
    G = K.galois_group(names = 't')
    #H is a group of order 3. if the CM-type "is" this group H, then it is not primitive. in turn, we set each element of the CM-type to be the identity and see if the other two elements can be the other two elements, non-trivial, elements of H. See Weng for proof
    H = []
    for g in G:
        if g.order() == 3 or g.order() == 1:
            H.append(g)
    for rho in L.complex_embeddings(prec):
        if is_H(phi0,phi1,phi2,rho,c,H) or is_H(phi0,phi2,phi1,rho,c,H) or is_H(phi1,phi0,phi2,rho,c,H) or is_H(phi1,phi2,phi0,rho,c,H) or is_H(phi2,phi0,phi1,rho,c,H) or is_H(phi2,phi1,phi0,rho,c,H):
            return False
    return True
       
def is_primitive_CMtype(phis, prec):
    """
    Given a list of three complex embeddings of a CM sextic field, checks if it is a primitive CM type
    """
    
    if is_CMtype(phis) == False:
        raise TypeError('This is not even a CM type')
    K = phis[0].domain()
    G = K.galois_group(type = 'pari')
    group = str(G.group().__pari__()[3])
    #the case where G = ZZ/6
    if group =='C(6) = 6 = 3[x]2':
        return is_primitive_galois(phis)
    #the case where G = ZZ/2 * S3
    elif group =='D(6) = S(3)[x]2':   
        return is_primitive_nongalois(phis, prec)
    #all others all CM types are primitive
    else:
        return True

