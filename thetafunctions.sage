"""

 This file has the functions to compute values of theta functions

"""

print 'heyoyo'

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

    t1=cputime()
    S +='default(realprecision,dig); \
    thetaval(del,eps,a,b,c,d,e,f, B)= \
    {s=0; t=0; \
    cutoff = -dig*log(10.0)-2*log(B); \
    for(i=-B,B, \
        for(j=-B,B, \
            for(k=-B,B, \
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
    t2 = cputime()
    #print 'sage time: ', t2-t1
    #print "gp time? : ", gp.eval('##')
    theta=gp.eval('thetan')
    Cec=ComplexField(prec)
    #print Cec(theta)
    return Cec(theta)

def theta_without_bound(period_matrix_entries, vec1, vec2, start_bound = 20, check = False, prec = 53):
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
    a=theta_with_bound(period_matrix_entries, vec1, vec2, start_bound, prec)
    values.append(a)
    iterates = 0
    equality = False
    while equality == False:
        iterates += 1
        a=theta_with_bound(period_matrix_entries, vec1, vec2, start_bound+20*iterates, prec)
        values.append(a)
        if compare(values[iterates],values[iterates-1]):
            equality = True
    return values[iterates]


def all_thetas(period_matrix,start_bound =20, prec = 664):

    all_evens = [[[0,0,0],[0,0,0]],[[1,0,0],[0,0,0]],[[0,1,0],[0,0,0]],[[0,0,1],[0,0,0]],[[0,0,1],[1,0,0]],[[0,0,1],[0,1,0]],[[1,1,0],[0,0,0]],[[1,0,1],[0,0,0]],[[0,1,1],[0,0,0]],[[0,0,0],[1,0,1]],[[0,0,0],[0,1,1]],[[0,0,0],[1,1,0]],[[1,1,1],[0,0,0]],[[0,0,0],[1,1,1]],[[0,1,1],[1,0,0]],[[1,0,1],[0,1,0]],[[1,1,0],[0,0,1]],[[0,1,0],[1,0,1]],[[0,0,1],[1,1,0]],[[1,0,0],[0,1,1]],[[1,0,1],[1,0,1]],[[1,1,0],[1,1,0]],[[0,1,1],[0,1,1]],[[1,0,1],[1,1,1]],[[1,1,0],[1,1,1]],[[1,1,1],[0,1,1]],[[1,1,1],[1,0,1]],[[1,1,1],[1,1,0]],[[0,1,1],[1,1,1]],[[0,0,0],[1,0,0]],[[0,0,0],[0,1,0]],[[0,0,0],[0,0,1]],[[1,0,0],[0,1,0]],[[1,0,0],[0,0,1]],[[0,1,0],[0,0,1]],[[0,1,0],[1,0,0]]]

    all_values = []

    for even in all_evens:
        #all_values.append([even, theta_with_bound(period_matrix,even[0],even[1],start_bound,prec)])
        all_values.append([even, theta_without_bound(period_matrix,even[0],even[1],start_bound,False,prec)])

    return all_values

def counting(values_list):
    count = 0
    for value in values_list:
        if value[1].abs() < 10.^(-2):
            count += 1
    return count



def compute_characteristic_sum_for_set_and_lists(S,eta):
  #same function as before but we work with lists
  sum=[[GF(2)(0),GF(2)(0),GF(2)(0)],[GF(2)(0),GF(2)(0),GF(2)(0)]]
  for i in S:
    for j in range(3):

      sum[0][j]=sum[0][j]+eta[i][0][j]
      sum[1][j]=sum[1][j]+eta[i][1][j]
  return sum


def thetaval(period_matrix_entries, S, eta):
  sum=compute_characteristic_sum_for_set_and_lists(S,eta)
  print sum
  return theta_with_bound(period_matrix_entries,sum[0],sum[1], bound)

def compute_rosenhain_coeffs(period_matrix_entries, U, T, eta, j):
    """
    This function computes the Rosenhain coefficient j, with j\geq
    2. We assume a0=0, a1=1. Three values are computed, for
    the three ways in which you can split the set in Takase's paper.
    """
    S=Set([0,1,j])
    Y=T.difference(S)
    X=Y.subsets()
    ajvec=[]
    signe=0
    aj=0
    for x in X:
        if (j<6):
	        a = j+1
        else:
	        a = 2
        if ((x.cardinality()==2) and (a in x)):
            M=x
            N=Y.difference(M)
            A=thetaval(period_matrix_entries, U.symmetric_difference(M.union(Set([0,1]))),eta)
            B=thetaval(period_matrix_entries, U.symmetric_difference(N.union(Set([0,1]))),eta)
            if ((IsZero(A)==0) and (IsZero(B)==0)):
                C=thetaval(period_matrix_entries, U.symmetric_difference(M.union(Set([0,j]))),eta)
                D=thetaval(period_matrix_entries, U.symmetric_difference(N.union(Set([0,j]))),eta)
                aj=((D*C)^2)/((A*B)^2)
                ajvec.append(aj)
    return ajvec


def IsZero(A,prec=664):
    n=A.imag()
    if (n<0):
        n=-n
    m=A.real()
    if (m<0):
        m=-m
    if (m<10^(-190) and n<10^(-190)):
        a=1
    else:
        a=0
    return a


def rosenhain_coeffs(period_matrix, M, prec=664):
    """M is the matrix we use to correct the system of characteristics"""


    lam=[]
    etaal=[]
    etaal.append(vector([GF(2)(1),GF(2)(0),GF(2)(0),GF(2)(0),GF(2)(0),GF(2)(0)]))
    etaal.append(vector([GF(2)(1),GF(2)(0),GF(2)(0),GF(2)(1),GF(2)(0),GF(2)(0)]))
    etaal.append(vector([GF(2)(0),GF(2)(1),GF(2)(0),GF(2)(1),GF(2)(0),GF(2)(0)]))
    etaal.append(vector([GF(2)(0),GF(2)(1),GF(2)(0),GF(2)(1),GF(2)(1),GF(2)(0)]))
    etaal.append(vector([GF(2)(0),GF(2)(0),GF(2)(1),GF(2)(1),GF(2)(1),GF(2)(0)]))
    etaal.append(vector([GF(2)(0),GF(2)(0),GF(2)(1),GF(2)(1),GF(2)(1),GF(2)(1)]))
    etaal.append(vector([GF(2)(0),GF(2)(0),GF(2)(0),GF(2)(1),GF(2)(1),GF(2)(1)]))
    #etaal.append(vector([GF(2)(0),GF(2)(0),GF(2)(0),GF(2)(0),GF(2)(0),GF(2)(0)]))



    U=Set([0,2,4,6])
    T1=Set([0,1,2,3,4,5,6])
    new_eta=[]
    eta_correction=[]
    for i in range(7):
           eta_correction.append(M*etaal[i])
           new_eta.append(matrix([[eta_correction[i][j] for j in range(3)],[eta_correction[i][j+3] for j in range(3)]]))

    for i in range(5):

           sys.stdout.flush()
           a=compute_rosenhain_coeffs(period_matrix,U,T1, new_eta,i+2)
           lam.append(a)


    return lam


#M=identity_matrix(GF(2),6,6)
bound=50
