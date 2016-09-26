def is_even(M):
    M1 = vector([M[0][0],M[1][0],M[2][0]])
    M2 = vector([M[3][0],M[4][0],M[5][0]])
    num = M1.dot_product(M2)
    if num == 0 % 2:
        return True
    if num == 1 % 2:
        return False

A = matrix(GF(2),[[1,0,0,0,0,0],[0,0,0,1,0,0],[0,1,0,0,0,0],[0,0,0,0,1,0],[0,0,1,0,0,0],[0,0,0,0,0,1]])
mumford_delta = matrix(GF(2),[[1],[1],[1],[1],[0],[1]])

        
def look_for_pairs(N):
    elements = []
    i = 0
    for s in SO(6, GF(2), e = 1):
        i += 1
        elements.append(A.inverse()*s*A)
        if i > N:
            break
    pairs ={}
    found_deltas = []
    for g in elements:
        delta = g*mumford_delta
        if is_even(delta):
            found = False
            for h in pairs.keys():
                if delta in found_deltas:
                    found = True
                    break
            if found == False:
                found_deltas.append(delta)
                delta.set_immutable()
                pairs[delta] = g
        if len(pairs) == 35:
            break
    if len(pairs) == 35:
        return pairs
    elif len(pairs) < 35:
        return []
            
N = 500
try_pairs = look_for_pairs(N)
while len(try_pairs) == 0:
    N *= 10
    try_pairs = look_for_pairs(N)
    print N
pairs = try_pairs
