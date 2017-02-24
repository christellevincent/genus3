load('class_cmfield.sage')
load('class_cmpoint.sage')
load('thetafunctions.sage')
load('precomputepairs.sage')
load('etamaps.sage')

S.<x> = PolynomialRing(QQ)
#K = CMFieldfromPoly(x^6+x^5+x^4+x^3+x^2+x+1,'a', prec = 25) # low precision so this is fast

K = CMField(the_fields[0][0],the_fields[0][1], prec = 600)
CMtype = K.all_CMtypes()[1] # pick a CM type among primitive CM types
polarization = K.princ_polarized(CMtype)[0] # pick a lattice with a principal polarization

CMpoint = CMPoint(K,CMtype,polarization[0],polarization[1]) # we construct the CM point

#print CMpoint.period_matrix()
#print CMpoint.counting() # counts how many even theta constants are 0
#print CMpoint.all_rosenhain_coeffs()

big_ros = CMpoint.all_rosenhain_coeffs()
little_ros = [ros[0] for ros in big_ros]
C = little_ros[0].parent()
P.<x>=C[]
f = x*(x-1) * prod([x-l for l in little_ros])
S = magma(f).ShiodaInvariants(nvals=2)
s = S[0].sage()
#weights = [2,3,4,5,6,7,8,9,10]
#J2 = s[0]
#J3 = s[1]
#l = J2/J3

#t = [s[i] / J2**(weights[i]/2) if weights[i] % 2 == 0 else 0 for i in range(9)] # if i in the field
t = s
pols = [algdep(c,1) for c in t]
T = [-p[0]/p[1] for p in pols]
#C_magma = magma(T).HyperellipticCurveFromShiodaInvariants()
#C = C_magma.sage()
t = [s[i].real_part() if weights[i] % 2 == 0 else 0 for i in range(9)]
