def SDist(p, q, wp, wq):
    """
    Subspace distance between p and q, with wp = w(p), wq = w(q)
    """
    d1 = 0 # = lambda(p,q) + Delta(p,q)
    d2 = 0 # = Dist(p,q)  # distance in subspace defined by inverse of w(p,q)
    return (d1, d2)

def ReachDist(p, q, r, wp, wq, wr):
    """
    subspace reachability. wp = w(p), wq = w(q), wr = w(r).
    r is the µ-Nearest Neighbor of p
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
