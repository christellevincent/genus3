load('class_cmfield.sage')
load('class_cmpoint.sage')
load('thetafunctions.sage')
load('precomputepairs.sage')
load('etamaps.sage')

S.<x> = PolynomialRing(QQ)
K = CMFieldfromPoly(x^6 - 12*x^5 + 39*x^4 - 62*x^3 + 744*x^2 + 228*x + 3473,'a', prec = 200) # medium precision so this is fast

CMtype = K.all_CMtypes()[1] # pick a CM type among primitive CM types
polarization = K.princ_polarized(CMtype)[0] # pick a lattice with a principal polarization

# we construct the CM point
CMpoint = CMPoint(K,CMtype,polarization[0],polarization[1])

#print CMpoint.reduced_period_matrix()
print CMpoint.counting() # counts how many even theta constants are 0
print CMpoint.all_rosenhain_coeffs()
