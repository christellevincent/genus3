def check_splitting(K,p):
    P = K.ideal(p)
    fact = list(P.factor())
    if len(fact) == 6:
        return p
    else:
        return None
        
def find_p(K,N):
    e = 0
    Pr = Primes()
    found = None
    while found == None:
        p = P.next(N+e*5)
        found = check_splitting(K,p)
    return found
 
def find_conj(K,L,id):
    conj = K.complex_conjugation()
    for l in L:
        if conj(l) == id:
            return True
    return False
          
def good_ideals(K,p):
    P = K.ideal(p)
    ideals = [id[0] for id in P.factor()]
    assert len(ideals) == 6
    one_of_each = []
    for id in ideals:
        if not(find_conj(K,one_of_each,id)):
            one_of_each.append(id)
    assert len(one_of_each) == 3
    return one_of_each
    
def list_of_units(K):
    G = UnitGroup(K)
    v,u1,u2 = G.gens_values()
    exps = [[j,k] for j in range(-99,100) for k in range(-99,100)]
    units = []
    for exp in exps:
        units.append(prod((u^e for u,e in zip([u1,u2],exp)), K(1)))
    return units
  
      
    
def find_norms(K,p,id):
    if not(id.is_principal()):
        return 'the ideal is not principal'
    else:
        norms = []
        pi = id.gens_reduced()[0]
        G = UnitGroup(K)
        v,u1,u2 = G.gens_values()
        fin = v.multiplicative_order()
        for u in list_of_units(K):
            newpi = u*pi
            if all([phi(newpi).abs()==sqrt(p) for phi in K.complex_embeddings()]):                    
                    norms.append([(1-v^e*pi).norm() for e in range(fin)])
        return norms
            
        
    
    
