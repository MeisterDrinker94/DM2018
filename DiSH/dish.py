from random import shuffle

"""
def SDist(p, q, wp, wq):
"""
    #Subspace distance between p and q, with wp = w(p), wq = w(q)
"""
    d1 = 0 # = lambda(p,q) + Delta(p,q)
    d2 = 0 # = Dist(p,q)  # distance in subspace defined by inverse of w(p,q)
    return (d1, d2)

def ReachDist(p, q, r, wp, wq, wr):
"""
    #subspace reachability. wp = w(p), wq = w(q), wr = w(r).
    #r is the miu-Nearest Neighbor of p
"""
    s1 = SDist(p, r, wp, wr)
    s2 = SDist(p, q, wp, wq)
    # s1 <= s2
    if s1[0] < s2[0] or (s1[0] == s2[0] and s1[1] <= s2[1]):
        return s2
    # s1 > s2
    else:
        return s1
"""

def epsilonNeighborhood(featureData,epsi,miu):
    """
    Find the epsi-neighborhoods for each point along the dimension given
    by the vector featureData. Neighborhoods with less than miu points
    are excluded.
    """
    #Array with index and Datavalue
    features = [(x,y) for x,y in enumerate(featureData)]
    #Sort the array by Datavalue
    features.sort(key=lambda tup: tup[1])
    
    #create Dictonary with objects for this candidate Subspace
    neighbourInfo = dict()
    
    for o in range(0,len(features)):
        i = o
        while(i >= 0 and abs(features[o][1]-features[i][1]) <= epsi):
            #adds indices to dictonary
            if features[o][0] in neighbourInfo.keys():
                neighbourInfo[features[o][0]].add(features[i][0])
            else:
                neighbourInfo[features[o][0]] = {features[i][0]}
                
            if features[i][0] in neighbourInfo.keys():
                neighbourInfo[features[i][0]].add(features[o][0])
            else:
                neighbourInfo[features[i][0]] = {features[o][0]}
            i -= 1;
    
    #Trim dict s.t. only neighbourhoods with more than miu items remain
    #THEY SHALL PASS!!!!        
    for key,val in neighbourInfo.items():
        if len(val) < miu:
            del neighbourInfo[key]
                    
    return neighbourInfo

def computeAllNeighborhoods(data, epsi, miu):
    """
    Compute neighborhoods for all features. Returns a list of dicts.
    """
    numFeatures = np.shape(data[0,:])
    neighborhood = []
    for f in range(0, numFeatures):
        neighborhood.append(epsilonNeighborhood(data[:,f], epsi, miu))
    return neighborhood;


def bestSubspaceForDataPoint(neighborList, index, epsi, miu):
    """
    neighborList: result from computeAllNeighborhoods()
    index:        index of the data point
    epsi, miu:    neighborhood parameters
    """
    subspace = set()
    # determine candidate attributes of the data point (1. step)
    candidates = set()
    for i in range(0, len(neighborList)):
        if index in neighborlist[i]:
            candidates.add(i)
    # find attribute with greatest neighborhood, delete from candidates and add to subspace (2. step)
    maxSize = -1
    maxAttrib = -1
    for ai in candidates:
        if len(neighborList[ai][index]) > maxSize:
            maxAttrib = ai
            maxSize = len(neighborList[ai][index])
    subspace.add(maxAttrib)
    candidates = candidates - maxAttrib
    # set current intersection (3.step)
    I = neighborList[maxAttrib][index]
    # 4./5. step
    while not len(candidates):
        # find attribute with greatest neighborhood
        maxSize = -1
        maxAttrib = -1
        maxIntersect = set()
        for ai in candidates:
            currentIntersect = I & neighborList[ai][index]
            if len(currentIntersect) > maxSize:
                maxAttrib = ai
                maxSize = len(currentIntersect)
                maxIntersect = currentIntersect
        # Add maxAttrib if intersection is large enough (4.a)
        if len(maxIntersect) >= miu:
            subspace.add(maxAttrib)
            I = maxIntersect
            candidates = candidates - maxAttrib
        # 4.b
        else:
            break
    return subspace


def subspacePreference(subspace, numFeatures):
    w = [1 if i in subspace else 0 for i in range(0, numFeatures)]


def main():
    test = [0,0.1,0.2,0.5,1.1,1.15,1.3,1.9,2,2.1,5000,5000.1,5000.2]
    shuffle(test)
    print(test)
    epsi = .5
    miu = 3
    
    d = epsilonNeighborhood(test,epsi,miu)
    
    print(d)

if __name__ == '__main__':
    main()

"""
algorithm DiSH(D, m, e):
    co <- cluster order     # initially empty
    pq <- priority queue ordered by ReachDist
    for p in D:
        compute w(p)        # Frequent itemset mining (algorithm on page 5)
        p.ReachDist <- infinity
        insert p into pq
    while pq not empty:
        o <- pq.pop()
        r <- m-NN of o w.r.t. SDist
        for p in pq:
            new_sr <- ReachDist(o,p,r)
            pq.update(p, new_sr)
        append o to co
    return co
"""
