# -*- coding: utf-8 -*-
from btmorph import WeightedQuickUnionUF
import itertools
import matplotlib.pyplot as plt
class Percolation:
    """
    Author: Irina Reshodko
    Direct translation from Java code
    
    We model a percolation system using an N-by-M grid of sites. 
    Each site is either open or blocked. A full site is an open site that can be 
    connected to an open site in the top row via a chain of neighboring 
    (left, right, up, down) open sites. We say the system percolates if 
    there is a full site in the bottom row. In other words, a system percolates 
    if we fill all open sites connected to the top row and that 
    process fills some open site on the bottom row. 
    (For the insulating/metallic materials example, 
    the open sites correspond to metallic materials, 
    so that a system that percolates has a metallic path from top to bottom, 
    with full sites conducting. For the porous substance example, 
    the open sites correspond to empty space 
    through which water might flow, so that a system 
    that percolates lets water fill open sites, 
    flowing from top to bottom.)
    """
    
    def __init__(self, sz):
        """
        Initialize
        Create a grid with given dimension, with all sites blocked
        
        Parameters:
        ------
        sz : Dimensions of a grid
        """
        # NxN grid to maintain open/blocked information
        self.sites = None
        sizeOfUnion = 1
        for i in range(0, len(sz)):
            sizeOfUnion *= sz[i]
            if sz[i] <= 0:
                raise IndexError(str(sz) + " : grid size must me greater than zero")
        # Create union object
        sizeOfUnion += 2 # + Two virtual sites
        # Union object for maintaining connections
        self.sitesUnion = WeightedQuickUnionUF(sizeOfUnion)
        # Union object to keep track of connection to the bottom
        self.sitesUnion_btm = WeightedQuickUnionUF(sizeOfUnion)
        # Virtual top site index
        self.top = 0;
        # Virtual bottom site index
        self.bottom = sizeOfUnion - 1
        # Create open/blocked info array; Set to false by default
        self.sites = [False]*(sizeOfUnion - 2)
        # Set grid size
        self.size = sz
        
    def to_1D(self, c):
        """
        3D or 2D coordinates to 1D
        
        Parameters:
        -------
        c : iterable
        2D or 3D coordinates
        
        Returns:
        -------
        index: 1D coordinate
        """
        if len(c) < 2 or len(c) > 3:
            raise TypeError(str(c) + " : should be 2d or 3d")
        index = -1
        #i + sz[0]*(j + sz[1]*k)
        if len(c) == 2: # 2 dimensions => 1 dimension
            index = c[0] + self.size[0]*c[1]
        elif len(c) == 3: #3 dimensions => 1 dimension
            index = c[0] + self.size[0]*(c[1] + self.size[1]*c[2])
        return index
    
    @staticmethod    
    def get_neighbours(coord):
        """
        Defines what sites are considered as neighbouring
        Should be implemented according to the application
        
        Parameters:
        coord : iterable
        Site coordinates
        
        Returns:
        neighbours : list of coordinates of neighbouring sites
        """
        # TODO : define what is neighbour according to your application !!!!
        sets = []
        for i in range(0, len(coord)):
            sets.append([coord[i] - 1, coord[i], coord[i] + 1])
        neighbours = itertools.product(*sets)
        return neighbours
        
    def open(self, coord):
        """
        open site with given coordinates *coords* if it is not open already
        
        Parameters:
        ------
        coord : coordinates of the site
        """
        #coord = list(coord)
        if len(coord) < 2 or len(coord) > 3:
            raise TypeError(str(coord) + " : should be 2d or 3d")
        for i in range(0, len(coord)):
            if coord[i] < 1 or coord[i] > self.size[i]:
                raise IndexError(str(coord) + ": Index is out of range (1..N) where N is" + str(self.size))
        c = list(coord)
        for i in range(0, len(c)):
            c[i] -= 1 # [1..N] => [0..N-1]
        index = self.to_1D(c)
        if self.sites[index] ==  True:
            return
        self.sites[index] = True
        # Connect to the virtual top
        index += 1
        if c[0] == 0:
            self.sitesUnion.union(index, self.top)
            self.sitesUnion_btm.union(index, self.top)
        # Connect to the virtual bottom    
        if c[0] == self.size[0] - 1:#map(lambda x: x - 1, self.size[:-1]):
            self.sitesUnion_btm.union(index, self.bottom)
        # if has open neighbours then connect sites
        neighbours = Percolation.get_neighbours(coord)
        for n in neighbours:
            n = list(n)
            # Don't check oneself
            if n == coord:
                continue
            ok = True
            for i in range(0, len(n)):
                if not (n[i] >= 1 and n[i] <= self.size[i]):
                    ok = False
                    break
            if not ok:
                continue            
            if not self.is_open(n):
                continue
            for m in range(0, len(n)):
                n[m] -= 1
            neigh_site = self.to_1D(n) + 1
            self.sitesUnion.union(index, neigh_site)
            self.sitesUnion_btm.union(index, neigh_site)
                        
    def is_open(self, coord):
        """
        Is site with the given coordinates *coord* open?
        Each site is either open or blocked.
        
        Parameters:
        ------
        coord : coordinates of the site
        
        Returns:
        ------
        True if open and False otherwise
        """
        coord = list(coord)
        if len(coord) < 2 or len(coord) > 3:
            raise TypeError(str(coord) + " : should be 2d or 3d")
        for i in range(0, len(coord)):
            if coord[i] < 1 or coord[i] > self.size[i]:
                raise IndexError(str(coord) + ": Index is out of range (1..N)")
        c = coord
        for i in range(0, len(c)):
            c[i] -= 1 # [1..N] => [0..N-1]
        index = self.to_1D(c)
        return self.sites[index]
        
    def is_full(self, coord):
        """
        is site with the given coordinates *coord* full?
        A full site is an open site that can be 
        connected to an open site in the top row via a chain of 
        neighboring open sites
        
        Parameters:
        ------
        coord : coordinates of the site
        
        Returns:
        ------
        True if full and False otherwise
        """
        if len(coord) < 2 or len(coord) > 3:
            raise TypeError(str(coord) + " : should be 2d or 3d")
        for i in range(0, len(coord)):
            if coord[i] < 1 or coord[i] > self.size[i]:
                raise IndexError(str(coord) + ": Index is out of range (1..N)")
        c = list(coord)
        for i in range(0, len(c)):
            c[i] -= 1 # [1..N] => [0..N-1]
        index = self.to_1D(c) + 1
        if not self.is_open(coord):
            return False
        return self.sitesUnion.connected(self.top, index)
    
    def percolates(self):
        """
        does the system percolate?
        We say the system percolates if 
        there is a full site in the bottom row. 
        In other words, a system percolates 
        if we fill all open sites connected to the top row and that 
        process fills some open site on the bottom row.
        """
        return self.sitesUnion_btm.connected(self.top, self.bottom)
    
    def plot(self):
        """
        Print out the percolation grid
        """
        for i in range(1, self.size[0]+1):
            for j in range(1, self.size[1]+1):
                if len(self.size) == 3:
                    if self.is_full((i, j, 1)):
                        print '*',
                    elif self.is_open((i,j, 1)):
                        print '0',
                    else:
                        print '-',
                else:
                    if self.is_full((i, j)):
                        print '*',
                    elif self.is_open((i,j)):
                        print '0',
                    else:
                        print '-',
            print('\n')
            
        
