"""

 This file has the functions to compute values of theta functions

"""

def theta_with_bound(period_matrix_entries, vec1, vec2, bound, prec = 664):
    """
    INPUT:
    period matrix entries (as a tuple of six entries)
    vec1 and vec2 give the theta characteristic
    bound gives the size of the rectangular box over which we will sum in the theta function
    prec gives the precision IN BITS (which is then converted to digits for pari/gp)

    OUTPUT:
    value of theta[vec1,vec2](period_matrix)
    """
    a, b, c, d, e, f= period_matrix_entries
    dig = floor(RR(prec * log(2,10)))
    
    S =  'A = %s;'%a
    S += 'B = %s;'%b
    S += 'C = %s;'%c
    S += 'D = %s;'%d
    S += 'E = %s;'%e
    S += 'F = %s;'%f
    S += 'bound = %s;' %bound
    S += 'dig = %s;' %dig
    S +='default(realprecision,dig); \
    thetaval(del,eps,a,b,c,d,e,f,bd)= \
    {s=0; t=0; \
    cutoff = -dig*log(10.0)-2*log(bd); \
    for(i=-bd,bd, \
        for(j=-bd,bd, \
            for(k=-bd,bd, \
                t = Pi*I*(i^2*a+i*j*b+i*k*c+1/2*i*(del[1]*a+del[2]*b+del[3]*c)+1/2*del[1]*(i*a+j*b+k*c)+ \
                          1/4*del[1]*(del[1]*a+del[2]*b+del[3]*c)+i*j*b+j^2*d+j*k*e+1/2*j*(del[1]*b+del[2]*d+del[3]*e)+\
                          1/2*del[2]*(i*b+j*d+k*e)+1/4*del[2]*(del[1]*b+del[2]*d+del[3]*e)+i*k*c+j*k*e+k^2*f+ \
                          1/2*k*(del[1]*c+del[2]*e+del[3]*f)+1/2*del[3]*(i*c+j*e+k*f)+ \
                          1/4*del[3]*(del[1]*c+del[2]*e+del[3]*f))+ \
            2*Pi*I*(1/2*i*eps[1]+1/4*del[1]*eps[1]+1/2*j*eps[2]+1/4*del[2]*eps[2]+1/2*k*eps[3]+1/4*del[3]*eps[3]);\
    if(real(t)>cutoff,s = s + exp(t)))));  return(s);}'

    gp(S)
    V =  'v1 = %s;'%vec1[0]
    V += 'v2 = %s;'%vec1[1]
    V += 'v3 = %s;'%vec1[2]
    V += 'v4 = %s;'%vec2[0]
    V += 'v5 = %s;'%vec2[1]
    V += 'v6= %s;'%vec2[2]
    V+= 'Vec1=[v1,v2,v3];'
    V+= 'Vec2=[v4,v5,v6];'
    V+= 'thetan=thetaval(Vec1,Vec2,A,B,C,D,E,F,bound);'
    gp(V)
    theta=gp.eval('thetan')
    Cec=ComplexField(prec)
    return Cec(theta)

def theta_without_bound(period_matrix_entries, vec1, vec2, start_bound = 20, check = False, prec = 664):
    """
    input:
    period matrix entries (as a tuple of six entries)
    vec1 and vec2 give the theta characteristic
    start_bound gives the bound for the first rectangular box over which we will sum in the theta function
    if check is set to True, will return all of the values obtained, not only the last, stable one
    prec gives the precision IN BITS (which is then converted to digits for pari/gp)

    This algorithm increases the size of the rectangular box over which we will sum in the theta function to see if there is convergence

    output:
    either the value of theta[vec1,vec2](period_matrix), to enough precision that increasing the box did not change the value (if check = False)
    or as many values of theta[vec1,vec2](period_matrix) as were necessary for the value to stabilizes
    """

    values = []
    values.append(theta_with_bound(period_matrix_entries, vec1, vec2, start_bound, prec))
    iterates = 0
    equality = False
    while equality == False:
        iterates += 1
        values.append(theta_with_bound(period_matrix_entries, vec1, vec2, start_bound+20*iterates, prec))
        if check:
            print "computing iteration number {0}".format(iterates)
            print values[iterates]
        if compare(values[iterates],values[iterates-1],prec+10):
            equality = True
            
    return values[iterates]


def compute_characteristic_sum_from_set_and_etas(S,eta_dict):
    """
    Given a dictionary of values eta_1, eta_2, ... eta_7 (giving a map eta), computes eta_S = sum_{i in S} eta_i
    Returns a list [[a,b,c],[d,e,f]], a,b,c,d,e,f are in QQ
    """
    sum = [[0,0,0],[0,0,0]]
    for i in S:
        sum[0][0] += QQ(ZZ(eta_dict[i][0][0])/2)
        sum[0][1] += QQ(ZZ(eta_dict[i][1][0])/2)
        sum[0][2] += QQ(ZZ(eta_dict[i][2][0])/2)
        sum[1][0] += QQ(ZZ(eta_dict[i][3][0])/2)
        sum[1][1] += QQ(ZZ(eta_dict[i][4][0])/2)
        sum[1][2] += QQ(ZZ(eta_dict[i][5][0])/2)
    return sum


def theta_from_char_and_list(all_values, eta_dict, characteristic):
    """
    inputs:
    the list of all theta values computed already for a given period matrix (outputted by all_thetas)
    eta_dict
    a vector [[a,b,c],[d,e,f]] obtained via compute_characteristic_sum_from_set_and_etas
    output:
    returns the value of theta[[a,b,c],[d,e,f]](Z)
    """
    double_char = [[i*2 for i in characteristic[0]],[i*2 for i in characteristic[1]]]
    reduced_char = [[i % 2 for i in double_char[0]],[i % 2 for i in double_char[1]]]
    int_vec = [[0,0,0],[0,0,0]]
    for i in range(2):
        for j in range(3):
            int_vec[i][j] = (double_char[i][j] - reduced_char[i][j])/2
    eta1 = [QQ(ZZ(list(i)[0]))/2 for i in eta_dict[1]]
    int_list = list(int_vec[0]) + list(int_vec[1])
    v1 = matrix(QQ,eta1)
    v2 = matrix(QQ,int_list)
    J = matrix(ZZ,[[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]])
    num_mat = v1*J*v2.transpose()
    num = num_mat[0][0]
    sign = exp(2*pi*I*num_mat[0][0])
    for pair in all_values:
        if pair[0] == reduced_char:
            return sign*pair[1]
