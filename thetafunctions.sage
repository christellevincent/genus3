"""

 This file has the functions to compute values of theta functions

"""
import warnings
warnings.simplefilter("ignore", UserWarning) #In order to remove warnings when computing eigenvalues


def theta_function(period_matrix, vec1, vec2, prec = 664, as_tuple = True):
    if as_tuple:
        a, b, c, d, e, f= period_matrix
        Z = Matrix(a.parent(),[[a,b,c],[b,d,e],[c,e,f]])
        
    elif not as_tuple:
        Z = period_matrix
        a= Z[0][0]
        b= Z[0][1]
        c= Z[0][2]
        d= Z[1][1]
        e= Z[1][2]
        f= Z[2][2]
    
    eig = min((Z.apply_map(lambda aux: aux.imag_part())).eigenvalues())
    
    dig = floor(RR(prec * log(2,10)))
    
    bound = ceil(RR(1./2 + sqrt(1./4 + dig*ln(10)/(pi.n(prec)*eig))))  
    
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


def theta_from_char_and_list(all_values, characteristic):
    """
    inputs:
    the list of all theta values computed already for a given period matrix (outputted by all_thetas)
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
    int_list = list(int_vec[0]) + list(int_vec[1])
    reduced_list = list(reduced_char[0]) + list(reduced_char[1])
    v1 = matrix(QQ,reduced_list)
    v2 = matrix(QQ,int_list)
    J = matrix(ZZ,[[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]])
    num_mat = v1*J*v2.transpose()
    num = num_mat[0][0]
    if num % 2 == 1:
        sign = -1
    elif num % 2 == 0:
        sign = 1
    for pair in all_values:
        if pair[0] == reduced_char:
            return sign*pair[1]
