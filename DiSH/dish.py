from heapq import heappush, heappop
import math
import numpy as np

def SDist(p, q, wp, wq,epsi):
    """
    #Subspace distance between p and q, with wp = w(p), wq = w(q)
    """
    #Builds the common subspace preference vector by joining gigedy the two lists
    wpq = [int(wpi and wqi) for wpi,wqi in zip(wp,wq)]
    
    #Lambda = the amounts of zeros in wpq
    Lambda = len(wpq) - sum(wpq)

    d1 = Lambda + Delta(p,q,wp,wq,wpq,epsi)
    
    #creates the inverse of wpq
    invwpq = [(1-x) for x in wpq]

    d2 = distSubspace(p,q,invwpq)  # distance in subspace defined by inverse of w(p,q)
    
    return (d1, d2)

def Delta(p,q,wp,wq,wpq,epsi):
    #return value
    delta_p_q = 0

    if (wpq == wp or wpq == wq) and distSubspace(p,q,wpq) > 2*epsi:
        delta_p_q = 1

    return delta_p_q

def distSubspace(p,q,wpq,pnorm=2):
    """
    Calculates a distance for two points in their subspace preference
    """
    dist = 0
    
    for pi,qi,wpqi in zip(p,q,wpq):
        if wpqi == 1:
            dist += abs(pi-qi)**pnorm


    return dist**(1/pnorm)

def ReachDist(p, q, r, wp, wq, wr, epsi):
    """
    #subspace reachability. wp = w(p), wq = w(q), wr = w(r).
    #r is the miu-Nearest Neighbor of p
    """
    s1 = SDist(p, r, wp, wr, epsi)
    s2 = SDist(p, q, wp, wq, epsi)
    # s1 <= s2
    if s1[0] < s2[0] or (s1[0] == s2[0] and s1[1] <= s2[1]):
        return s2
    # s1 > s2
    else:
        return s1

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
    #for key,val in neighbourInfo.items():
    #    if len(val) < miu:
    #        del neighbourInfo[key]
    neighbourInfo = { k:v for k,v in neighbourInfo.items() if len(v) >= miu }
                    
    return neighbourInfo

def computeAllNeighborhoods(data, epsi, miu):
    """
    Compute neighborhoods for all features. Returns a list of dicts.
    """
    numFeatures = data.shape[1]
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
        if index in neighborList[i]:
            candidates.add(i)
    # find attribute with greatest neighborhood, delete from candidates and add to subspace (2. step)
    maxSize = -1
    maxAttrib = -1
    for ai in candidates:
        if len(neighborList[ai][index]) > maxSize:
            maxAttrib = ai
            maxSize = len(neighborList[ai][index])
    subspace.add(maxAttrib)
    candidates.remove(maxAttrib)
    # set current intersection (3.step)
    I = neighborList[maxAttrib][index]
    # 4./5. step
    while len(candidates):
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
            candidates.remove(maxAttrib)
        # 4.b
        else:
            break
    return subspace

def miuNearestNeighbor(indexP, data, preferences, epsi, miu):
    """
    Calculate the miu neareset neighboor for point p in respect to SDist
    """
    neighboorsDistance = []

    for i in range(data.shape[0]):
            q = data[i,:]
            wq = preferences[i]
            p = data[indexP,:]
            wp = prefrences[indexP]
            d = SDist(p,q,wp,wq,epsi)
            neighboorsDistance.append((d,i))
    
    #Sort all neigboors to size
    neighboorsDistance.sort()

    return neighboorsDistance[miu][1]
    

def subspacePreference(subspace, numFeatures):
    return [1 if i in subspace else 0 for i in range(0, numFeatures)]

def dish(data, epsi, miu):
    clusterOrder = []
    # precompute all epsi-neighborhoods with more than miu points
    numFeatures = data.shape[1]
    neighborList = computeAllNeighborhoods(data, epsi, miu)
    # find subspace preference vectors for each point
    preferences = []
    pq = [] # empty priority queue
    for o in range(0, data.shape[0]):
        subspace = bestSubspaceForDataPoint(neighborList, o, epsi, miu)
        print("subspace:", subspace)
        wo = subspacePreference(subspace, numFeatures)
        preferences.append(wo)
        heappush(pq, (math.inf, o))
    while pq:
        o = heappop(pq)
        r = miuNeearestNeighbor # TODO: nearest neighbor
        for idx, p in enumerate(pq):
            newSr = ReachDist(data[o[1],:], data[p[1],:], data[r,:], preferences[o[1]], preferences[p[1]], preferences[r], epsi)
            pq[idx] = (newSr, p[1])
        clusterOrder.append(o[1])
    return clusterOrder

def testDish():
    data = np.array([[1.0,  3.0],
                     [1.5,  0.0],
                     [0.0,  3.5],
                     [1.3, 10.0]])
    epsi = 0.5
    miu = 2
    print(dish(data, epsi, miu))

def main():
    testDish()

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
