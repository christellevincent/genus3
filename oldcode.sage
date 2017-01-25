def get_tau3_from_cubic_and_abc(cubic, abc_list):
    """
    Inputs a cubic (where the cubic number field is assuemd to have class number 1), creates sextic as a relative extension, then defines it as an absolute extension; returns the absolute extension and the ring of integers

    Here "abc_list" is sticking with Christelle's convention -- e.g., [2,2,2], which corresponds to -1* each entry in the list we generated earlier.

    EXAMPLES::
        sage: S.<y> = QQ['y']
        sage: g = y^3+y^2-2*y-1
        sage: get_tau3_from_cubic_and_abc(g, [2,2,2])
    """
    a, b, c = abc_list
    S.<y> = NumberField(cubic)
    if S.class_number() != 1:
        raise TypeError

    #creates the cubic in pari
    s = 'k = bnfinit(%s);'%cubic

    # creates the sextic as a relative extension
    s += 'a = %s;'%a
    s += 'b = %s;'%b
    s += 'c = %s;'%c
    s += 'cubic = %s;'%cubic
    s += 'R = rnfinit(k, x^2+(Mod(a + b*y + c*y^2, cubic)));'

    # defines the polynomial of the sextic
    s += 'f = R[11][1];'

    # defines the sextic as an absolute extension
    s += 'K = bnfinit(f);'

    #K.gen is the list of generators of the ideal class group of K. K.gen[1] is the first generator. Multiplying by K.zk is, as far as I can tell, only a technicality. This is thinking of this ideal from K as an ideal in R, i.e. a module over O_k, the ring of integers of k.
    s += 'ideal = K.zk * K.gen[1];'
    s += 'RC = rnfidealabstorel(R, ideal);'

    #the first part of RC is a matrix A = [[1,6],[0,1]], and each entry is in k. It is Pari trying to give us two numbers, z1 and z2, that are in R, and the entries of the matrix give us the zi's in term of the relative basis from k to R.
    s += 'M = R.zk[1]*RC[1];'
    s += 'z1= M[1]; z2= M[2];'

    # the second part of RC are two ideals in k. Since k has class number 1, they must be principal ideals. This returns the generator of each ideal, but in terms of the integer basis for k.
    s += 'F1 = bnfisprincipal(k,RC[2][1])[2];'
    s += 'F2 = bnfisprincipal(k,RC[2][2])[2];'

    # here we get the "actual" generators (as numbers, not in terms of the basis)
    s += 'f1 = k.zk*F1;  f2 = k.zk*F2;'

    # the ideal is z1*(f1)+z2*(f2), but we want one tau (i.e. we want the ideal in the form O_k + tau O_k). This gets the tau.
    s += 'tau2 = nfeltdiv(K, z1*f1,z2*f2);'

    # we are interested in the number tau3, where x is a root of x^6+14x^4+56x^2+56 (the f defining K in the very first cell) We want a tau that is suitable, i.e. that has negative imaginary part, and positive imaginary part under the two automorphisms of K that we have chosen to be the CM type. We will do this in sage, in good_tau.sagews
    s += 'tau3 = K.zk*tau2;'
    gp(s)
    tau3 = gp.eval('tau3')
    R.<x> =QQ['x']
    return R(tau3)

def get_period_matrix_entries_from_tau3(cubic, abc_list, tau3):
    """
    computes the "right" tau from tau3, computes the symplectic basis, finds the period matrix, and returns the list [A, B, C, D, E, F] (the distinct entries in the period matrix, as elements of CC) in Christelle's notation

    WARNING: we hard-coded the value of "a" here. fix this!
    """
    R.<y> = QQ['y']
    K0.<a> = NumberField(cubic)
    S.<x> =  K0['x']
    A,B,C = abc_list
    fpol=x^2+A+B*a+C*a^2
    R.<b> = NumberField(fpol)
    R_abs= R.absolute_polynomial()
    K.<c> = NumberField(R_abs.change_variable_name('y'))

    tau = K(tau3.change_variable_name('c'))
    L = []
    ri = K.automorphisms()
    for i in range(len(ri)):
        if (ri[i]^3 != K.automorphisms()[0]) and (ri[i]^2 != K.automorphisms()[0]):
            L.append(i)
    r0 = K.automorphisms()[0]
    rgen = K.automorphisms()[L[0]]

    Gal_conj=[tau, rgen(tau), (rgen^2)(tau), (rgen^3)(tau), (rgen^4)(tau),(rgen^5)(tau), tau, (rgen)(tau)]

    CCGal_conj=[]
    for i in range(len(Gal_conj)):
        CCGal_conj.append(Gal_conj[i].complex_embeddings()[0])
    suitables = []
    for i in range(6):
        if CCGal_conj[i].imag() > 0:
            if CCGal_conj[i+1].imag() > 0:
                if CCGal_conj[i+2].imag() > 0:
                    suitables.append([i,Gal_conj[i],CCGal_conj[i]])

    tau = suitables[0][1]

    #to write a in terms of c
    Z.<w> = K['w']
    cubic_over_Z = cubic.change_variable_name('w')
    #print 'this is cubic over Z:' ,cubic_over_Z
    #print 'this is with the coercion: ', Z(cubic_over_Z)
    #print 'roots without the coercion: ', cubic_over_Z.roots()
    #print 'roots with the coercion: ', Z(cubic_over_Z).roots()
    a =  Z(cubic_over_Z).roots()[0][0]
    # a = (-1 + sqrt(-3 - 2*c^2))/2
    #print cubic_over_Z(a)
    ac=a
    print ac
    relative=w^2+A+B*a+C*a^2

    b=relative.roots()[0][0]
    #the symplectic basis
    a1=tau*a^2
    a2=tau*a
    a3=tau
    a4=1
    a5=-rgen(a)-(rgen^2)(a)
    a6=rgen(a)*(rgen^2)(a)

    A1=matrix([[a1,a2,a3],[rgen(a1),rgen(a2),rgen(a3)],[(rgen^2)(a1),(rgen^2)(a2),(rgen^2)(a3)]])
    A2=matrix([[a4,a5,a6],[rgen(a4),rgen(a5),rgen(a6)],[(rgen^2)(a4),(rgen^2)(a5),(rgen^2)(a6)]])
    CM=A1.inverse()*A2
    print "CM is positive definite: ",  CM.is_positive_definite()
    print "CM is symmetric: ", CM.is_symmetric()

    #the entries we need for chi_18 are:
    A= CM[0][0]
    B= CM[0][1]
    C= CM[0][2]
    D= CM[1][1]
    E= CM[1][2]
    F= CM[2][2]
    diff=(a-rgen(a))*(a-(rgen^2)(a))

    print diff.complex_embeddings()
    print (rgen)(diff).complex_embeddings()
    print (rgen^2)(diff).complex_embeddings()

    return [A.complex_embeddings()[0], B.complex_embeddings()[0], C.complex_embeddings()[0], D.complex_embeddings()[0], E.complex_embeddings()[0],    F.complex_embeddings()[0]]
