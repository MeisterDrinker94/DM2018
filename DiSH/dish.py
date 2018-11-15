import csv
from heapq import heappush, heappop, heapify
import math
import numpy as np
from SyntheticData import createSynthetic
import sys
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def subspaceDimensionality(w):
    """
    The subspace dimensionality is the number of zeros in the
    subspace preference vector w.
    """
    return len(w) - sum(w)

def combinedSubspace(wp, wq):
    """
    Given preference vector w(p) and w(q), compute the combined
    subspace w(p, q) by attribute-wise logical and.
    """
    return [int(wpi and wqi) for wpi, wqi in zip(wp, wq)]

def Delta(p,q,wp,wq,wpq,epsi):
    """
    Compute Delta for use in SDist(p, q, wp, wq, epsi)
    """
    if (wpq == wp or wpq == wq) and distSubspace(p,q,wpq) > 2*epsi:
        return 1
    else:
        return 0

def SDist(p, q, wp, wq,epsi):
    """
    Compute the subspace distance between p and q, with
    preference vectors wp and wq.
    """
    wpq = combinedSubspace(wp, wq)
    Lambda = subspaceDimensionality(wpq)
    d1 = Lambda + Delta(p,q,wp,wq,wpq,epsi)
    # distance in subspace defined by inverse of w(p,q) 
    invwpq = [0 if x == 1 else 1 for x in wpq]
    d2 = distSubspace(p,q,invwpq)  
    return (d1, d2)

def distSubspace(p,q,wpq,pnorm=2):
    """
    Calculates a distance for two points in their preferred subspace 
    """
    dist = 0
    for pi,qi,wpqi in zip(p,q,wpq):
        if wpqi == 1:
            dist = dist + abs(pi-qi)**pnorm
    return dist**(1/pnorm)

def ReachDist(p, q, r, wp, wq, wr, epsi):
    """
    Calculate the subspace reachability of q from p, where r
    is the miu-nearest neighbor of p. wp, wq and wr are the
    subspace preference vectors of p, q and r, respectively.
    Computes max(SDist(p,r),SDist(p,q))
    """
    s1 = SDist(p, r, wp, wr, epsi)
    s2 = SDist(p, q, wp, wq, epsi)
    return max(s1,s2)

def epsilonNeighborhood(featureData,epsi,miu):
    """
    Find the epsi-neighborhoods for each point along the dimension given
    by the vector featureData. Neighborhoods with less than miu points
    are excluded.
    """
    # Array with index and Datavalue
    features = [(x,y) for x,y in enumerate(featureData)]
    # Sort the array by Datavalue
    features.sort(key=lambda tup: tup[1])
    
    # create Dictionary with objects for this candidate Subspace
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
            i = i - 1;
    
    #Trim dict s.t. only neighbourhoods with more than miu items remain
    #THEY SHALL PASS!!!!        
    neighbourInfo = { k:v for k,v in neighbourInfo.items() if len(v) >= miu }
                    
    return neighbourInfo

def computeAllNeighborhoods(data, epsi, miu):
    """
    Compute neighborhoods for all features. Returns a list of dicts,
    such that the neighborhoods associated with attribute i are stored
    in the i-th element of the list.
    """
    numFeatures = data.shape[1]
    neighborhood = []
    for f in range(0, numFeatures):
        neighborhood.append(epsilonNeighborhood(data[:,f], epsi, miu))
    return neighborhood

def bestSubspaceForDataPoint(neighborList, index, epsi, miu):
    """
    Calculate, for some data point, the set of attributes that
    form the best subspace.

    neighborList: result from computeAllNeighborhoods()
    index:        index of the data point
    epsi, miu:    neighborhood parameters
    """
    subspace = set()
    # determine candidate attributes of the data point (1. step),
    # i.e., attributes for which there is a large enough neighborhood
    candidates = set()
    for i, N in enumerate(neighborList):
        if index in N:
            candidates.add(i)
    # if there are no candidate attributes, there is no best subspace
    if not candidates:
        return set()
    # find the attribute with the greatest neighborhood,
    # delete it from the candidates and add it to the subspace (2. step)
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
        # find attribute for which the intersection with I is largest
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

def miuNearestNeighbor(pIdx, data, preferences, epsi, miu):
    """
    Calculate the miu neareset neighbor for point p with respect to SDist
    """
    neighborsDistance = []
    p = data[pIdx,:]
    wp = preferences[pIdx]
    for i in range(data.shape[0]):
            q = data[i,:]
            wq = preferences[i]
            d = SDist(p,q,wp,wq,epsi)
            neighborsDistance.append((d,i))
    # find the miu-th nearest neighbor 
    neighborsDistance.sort()
    return neighborsDistance[miu-1][1]
    
def subspacePreference(subspace, numFeatures):
    """
    Build the subspace preference vector based on the best subspace
    as found in bestSubspaceForDataPoint

    subspace:       best subspace of some data point
    numFeatures:    dimensionality of the data set
    """
    return [1 if i in subspace else 0 for i in range(0, numFeatures)]

def dish(data, epsi, miu):
    """
    Perform the DiSH algorithm.

    data:   the data to cluster
    epsi:   the neighborhood size
    miu:    the minimum number of points in neighborhoods
    """
    clusterOrder = []
    numFeatures = data.shape[1]
    # precompute all epsi-neighborhoods with more than miu points
    neighborList = computeAllNeighborhoods(data, epsi, miu)
    # subspace preference vectors of all points
    preferences = []
    # empty priority queue to store subspace reachabilities of points
    # as tuples of (ReachDist, point index)
    pq = []
    print("Compute subspace preferences")
    sys.stdout.flush()
    # find subspace preference vectors and initialize pq
    for o in range(0, data.shape[0]):
        subspace = bestSubspaceForDataPoint(neighborList, o, epsi, miu)
        wo = subspacePreference(subspace, numFeatures)
        preferences.append(wo)
        heappush(pq, ((math.inf,math.inf), o))
    print("Compute cluster order")
    sys.stdout.flush()
    # compute the cluster order
    while pq:
        o = heappop(pq)
        r = miuNearestNeighbor(o[1], data, preferences, epsi, miu)
        # update subspace reachability for each point in the queue
        for idx, p in enumerate(pq):
            newSr = ReachDist(data[o[1],:], data[p[1],:], data[r,:],
                              preferences[o[1]], preferences[p[1]], preferences[r], epsi)
            pq[idx] = (newSr, p[1])
            heapify(pq)
        clusterOrder.append(o[1])
    return clusterOrder, preferences

def extractCluster(clusterOrder,preferences, data, epsi):
    """
        Finds the clusters given the cluster order and subspace preferences

        clusterOrder:   cluster order computed in dish()
        preferences:    preference vectors of all points computed in dish()
        data:           data to be clustered
        epsi:           neighborhood size

        Returns a list of lists, where
            clusters[i][0] are the points in cluster i
            clusters[i][1] is the preference vector
            clusters[i][2] is the centroid
    """
    print("extracting Clusters")
    sys.stdout.flush()
    # Clusters with points, preference Vector, center.
    # The first point has no predecessor, so we create a new cluster.
    # To simplify updating the cluster center, we only store the sum of the points
    # in the cluster for now. Division by the number of points yields the cluster center.
    clusters = []
    clusters.append([[clusterOrder[0]],preferences[0],data[clusterOrder[0],:]])

    for i in range(1,len(clusterOrder)):
        oIdx = clusterOrder[i]
        pIdx = clusterOrder[i-1]
        foundCluster = False
        for c in clusters:
            wc = c[1]
            wop = combinedSubspace(preferences[oIdx], preferences[pIdx])
            if wc == wop and distSubspace(data[oIdx,:], c[2]/len(c[0]), wop) <= 2*epsi:
                c[0].append(oIdx)
                c[2] = (c[2] + data[oIdx,:])
                foundCluster = True
        # No appropriate cluster found, create a new one
        if not foundCluster:
            clusters.append([[oIdx],preferences[oIdx],data[oIdx,:]])

    # compute the cluster centers from the sum of the points in a cluster
    for c in clusters:
        c[2] = c[2]/len(c[0])
    return clusters

def buildHierarchy(clusters, epsi):
    """
    Build the subspace cluster hierarchy from the clusters computed in extractCluster()
    """
    # dimensionality of the data
    dim = len(clusters[0][2])

    # dictionary of clusters and their parents
    parentDict = {}
     
    for ciIdx, ci in enumerate(clusters):
        lambdaci = subspaceDimensionality(ci[1]) 
        for cjIdx, cj in enumerate(clusters):
            lambdacj = subspaceDimensionality(cj[1])
            if lambdacj > lambdaci:
                wij = combinedSubspace(ci[1], cj[1])
                dist = distSubspace(ci[2],cj[2],wij)
                
                # Check if there is a cluster that is a parent of ci
                # and has lower dimensionality lambda
                ciHasParents = False
                if ciIdx in parentDict.keys():
                    for c in parentDict[ciIdx]:
                        lambdac = subspaceDimensionality(clusters[c][1])
                        if lambdac < lambdacj:
                            ciHasParents = True
                            break

                if lambdacj == dim or (dist <= 2*epsi and not ciHasParents):
                    if cjIdx in parentDict.keys():
                        parentDict[cjIdx].append(ciIdx)
                    else:
                        parentDict[cjIdx] = [ciIdx]
    return parentDict


def testDish():
#    data = np.array([[1.0,  3.0],
#                     [1.5,  0.0],
#                     [0.0,  3.5],
#                     [1.3,  3.0],
#                     [1.2,  3.1]])
#    data = createSynthetic(noisePoints=0)
    data = np.array([[0,0,0],
                    [0.25,0,0],
                    [0.5,0,0],
                    [0,1,1],
                    [.5,1,1],
                    [1,1,1]])
#    data = [];
#    with open('breast-cancer-wisconsin.data', 'r') as csvfile:
#        reader = csv.reader(csvfile, delimiter=',')
#        for row in reader:
#            data.append(row[1:])
#    data = np.array(data, dtype=np.float64)
    epsi = 0.3
    miu = 3 
    order, prefs = dish(data, epsi, miu)
    print("Cluster Order", order)
    clusters = extractCluster(order,prefs,data,epsi)
    hierarchy = buildHierarchy(clusters, epsi)
    #print("parentDict:", hierarchy)
    print("len(parentDict):", len(hierarchy))
    print("len(clusters):", len(clusters))
    #print("clusters:", clusters)
    plotData(clusters,data)

def plotData(clusters, Data):
    fig = plt.figure()
    
    ax = fig.add_subplot(111, projection='3d')
    
    markerStyle = ['g','r','b','y','m']
    
    col = -1
    for c in clusters:
        col = col + 1
        if col == 5:
            break
        for dindex in c[0]:
            ax.scatter(Data[dindex,0],Data[dindex,1],Data[dindex,2], color=markerStyle[col])
    plt.xlabel("x-label")
    plt.ylabel("y-label")
    plt.show()

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
