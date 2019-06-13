load('class_cmfield.sage')
load('class_cmpoint.sage')
load('thetafunctions.sage')
load('precomputepairs.sage')
load('etamaps.sage')

S.<x> = PolynomialRing(QQ)
K = CMFieldfromPoly(x^6 + 29*x^4 + 180*x^2 + 64,'a', prec = 200) # medium precision so this is fast

CMtype = K.all_CMtypes()[1] # pick a CM type among primitive CM types
polarization = K.princ_polarized(CMtype)[0] # pick a lattice with a principal polarization

# we construct the CM point
CMpoint = CMPoint(K,CMtype,polarization[0],polarization[1])

#print CMpoint.reduced_period_matrix()
print CMpoint.counting() # counts how many even theta constants are 0
print CMpoint.all_rosenhain_coeffs()

# x^6-14*x^3+63*x^2+168*x+161 pinar's field
# x^6 - x^5-11*x^4+15*x^3+73*x^2-127*x+211
