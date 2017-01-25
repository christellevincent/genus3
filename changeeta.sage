def list_of_Sp(N):
    A = matrix(GF(2),[[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,0,0,1],[0,0,0,0,1,0],[0,0,0,1,0,0]]) 
    elements_of_Sp = []
    i = 0
    for s in Sp(6,GF(2)):
        i += 1
        elements_of_Sp.append(A.inverse()*s*A)
        if i > N:
            break
    return elements_of_Sp

def change_eta_matrices(CMPoint,mats):
    try:
        eta_dict = CMPoint._eta_dict
    except:
        eta_dict = CMPoint.eta_dict()
    try:
        vanishing_char = CMPoint._vanishing_char
    except:
        vanishing_char = CMPoint.vanishing_char()
        
    delta = matrix(GF(2),[[vanishing_char[0][0]],[vanishing_char[0][1]],[vanishing_char[0][2]],[vanishing_char[1][0]],[vanishing_char[1][1]],[vanishing_char[1][2]]])
    
    good_matrices = []
    
    for mat in mats:
        eta_U = matrix(GF(2),6,1)
        U_set = []
        new_dict = {1: mat*eta_dict[1], 2: mat*eta_dict[2], 3: mat*eta_dict[3], 4: mat*eta_dict[4], 5:mat*eta_dict[5], 6: mat*eta_dict[6], 7: mat*eta_dict[7]}
        for j in range(1,8):
            if is_even(new_dict[j]):
                eta_U += new_dict[j]
                U_set.append(j)        
        if delta == eta_U:
            good_matrices.append(mat)
    return good_matrices
    
def get_the_stuff(CMPoint,mat,old_dict):
    eta_dict = old_dict
    try:
        vanishing_char = CMPoint._vanishing_char
    except:
        vanishing_char = CMPoint.vanishing_char()
    
    delta = matrix(GF(2),[[vanishing_char[0][0]],[vanishing_char[0][1]],[vanishing_char[0][2]],[vanishing_char[1][0]],[vanishing_char[1][1]],[vanishing_char[1][2]]])    
                
    eta_U = matrix(GF(2),6,1)
    U_set = []
    new_dict = {1: mat*eta_dict[1], 2: mat*eta_dict[2], 3: mat*eta_dict[3], 4: mat*eta_dict[4], 5:mat*eta_dict[5], 6: mat*eta_dict[6], 7: mat*eta_dict[7]}
    for j in range(1,8):
        if is_even(new_dict[j]):
            eta_U += new_dict[j]
            U_set.append(j)
    assert delta == eta_U
    return new_dict, Set(U_set)
    
def all_in_one(CMPoint,mat,old_dict):
    CMPoint._eta_dict, CMPoint._U_set = get_the_stuff(CMPoint,mat,old_dict)
    print CMPoint._eta_dict
    print CMPoint._U_set
    bigros = CMPoint.all_rosenhain_coeffs()
    littleros = [bigros[i][0].real() for i in range(5)]
    polys = [algdep(littleros[i],20,CMPoint._prec) for i in range(5)]
    for ros in littleros:
        print ros
    for poly in polys:
        print poly
    return 'all done'

    
        
    
