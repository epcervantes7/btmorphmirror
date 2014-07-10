# -*- coding: utf-8 -*-

class WeightedQuickUnionUF:
    """
    Direct translation of java code from here:
    http://algs4.cs.princeton.edu/15uf/
    |  The *WeightedQuickUnionUF* class represents a union-find data structure.
    |  It supports the *union* and *find* operations, along with
    methods for determinig whether two objects are in the same component
    and the total number of components.
    |  This implementation uses weighted quick union by size (without path compression).
    |  Initializing a data structure with *N* objects takes linear time.
    |  Afterwards, *union*, *find*, and *connected* take
    logarithmic time (in the worst case) and *count* takes constant
    time.
    |  For additional documentation, see <a href="http://algs4.cs.princeton.edu/15uf">Section 1.5</a> of
    Algorithms, 4th Edition by Robert Sedgewick and Kevin Wayne.
    """
    def __init__(self, N):
        """
        Initializes an empty union-find data structure with N isolated components 0 through N-1
        complexity: N
        
        Paremeters:
        --------
        N : the number of objects
        """
        self.id = [0]*N # id[i] = parent of i
        self.sz = [0]*N # sz[i] = number of objects in subtree rooted at i
        self.count = N # number of components
        for i in range(0, N):
            self.id[i] = i
            self.sz[i] = 1

    def find(self, p):
        """
        Returns the component identifier for the component containing site *p*
        complexity: ln N
        
        Parameters:
        -------
        p : the integer representing one site
        
        Returns:
        -------
        the component identifier for the component containing site *p*
        """
        while p != self.id[p]:
            p = self.id[p]
        return p
        
    def connected(self, p, q):
        """
        Are the two sites *p* and *q* in the same component?
        complexity: 2 * ln N
        
        Parameters:
        -------
        p : the integer representing one site
        q : the integer representing the other site
        
        Returns:
        -------
        **true** if the two sites *p* and *q* are in the same component, 
        and **false** otherwise
        """
        return self.find(p) == self.find(q)
        
    def union(self, p, q):
        """
        Merges the component containing site *p* with the component containing site *q*
        complexity: ln N
        
        Parameters:
        -------
        p : the integer representing one site
        q : the integer representing the other site        
        """
        rootP = self.find(p)
        rootQ = self.find(q)
        if rootP == rootQ:
            return
        if self.sz[rootP] < self.sz[rootQ]:
            self.id[rootP] = rootQ
            self.sz[rootQ] += self.sz[rootP]
        else:
            self.id[rootQ] = rootP
            self.sz[rootP] += self.sz[rootQ]
        self.count = self.count - 1
        
    def __str__(self):
        print "The number of components is: {}".format(self.count)



