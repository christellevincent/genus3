load('class_cmfield.sage')
load('class_cmpoint.sage')
load('thetafunctions.sage')
load('precomputepairs.sage')
load('etamaps.sage')

S.<x> = PolynomialRing(QQ)
K = CMFieldfromPoly(x^6+x^5+x^4+x^3+x^2+x+1,'a', prec = 25) # low precision so this is fast

CMtype = K.all_CMtypes()[1] # pick a CM type among primitive CM types
polarization = K.princ_polarized(CMtype)[0] # pick a lattice with a principal polarization

CMpoint = CMPoint(K,CMtype,polarization[0],polarization[1]) # we construct the CM point

#print CMpoint.period_matrix()
#print CMpoint.counting() # counts how many even theta constants are 0
#print CMpoint.all_rosenhain_coeffs()
