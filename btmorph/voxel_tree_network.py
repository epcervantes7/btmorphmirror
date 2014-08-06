"""
Created on Fri Jul 11 17:32:44 2014
@author: Irina Reshodko

Defines data structure to store a network of voxelized trees
"""
import math
import btmorph
from btmorph import VoxelGrid
from btmorph import Percolation
from numpy import random as rnd
from numpy import mean
from matplotlib import pyplot as plt

class VoxelizedTreeNetwork:
    """
    Defines data structure to store a network of voxelized trees
    """
    def __init__(self, dim, voxelSz, trees = None, percDirection = 0):
        """
        Initialize
        
        Parameters:
        -------
        dim : iterable
        Global dimensions of the network in x, y and z direction in micrometers.
        |  If dim[z] == 0 then a 2d voxel grid will be generated
        voxelSz : number
        Dimensions of one voxel in x, y and z direction in micrometers. Must be greater than 0 
        |  This parameter determines resolution of the voxel grid.
        |  Will be adjusted according to the requirements of the grid size (power of two)
        trees : iterable of :class:`btmorph.btstructs2.STree2`
        Trees to be added to the network. Can be added later one be one
        percDirection : int
        Integer >= 0 and <=2. Specifies direction of percolation. Can only be set once.
        """
        if not len(dim) == 3:
            raise TypeError('Dimensions: Iterable of length 3 expected, got ' + type(dim) + "of length " + len(dim))
        if voxelSz <= 0:
            raise TypeError('Voxel size must be greater than zero. Got ' + str(voxelSz))
        for el in dim:
            if el < 0:
                raise TypeError('Dimensions must be greater than or equal to zero. Got ' + str(dim))
        self.dim = dim
        self.__voxelSz__ = voxelSz
        dx,dy,dz = self.dim
        res = [int(2**round(math.log(dx/voxelSz, 2))), int(2**round(math.log(dy/voxelSz, 2))), 1]
        if dz > 0.0:
            res[2] = int(2**round(math.log(dz/voxelSz, 2)))
        self.vg = VoxelGrid(self.dim, res)
        self.dim = self.vg.dim
        self.__percDirection__ = percDirection
        # Put percolation index in first place
        percRes = list(self.vg.res)
        t = percRes[self.__percDirection__]
        percRes[self.__percDirection__] = percRes[0]
        percRes[0] = t
        self.perc = Percolation(percRes)
        if trees == None:
            return
        for tree in trees:
            self.add_tree(tree)
        
    def add_tree(self, tree, offset = None):
        """
        Add a tree to the network
        
        Parameters:
        -------
        tree : :class:`btmorph.btstructs2.STree2`
        A tree to be voxelized and added to the network
        offset : iterable
        Real coordinates (micrometers in 3 directions)
        Tree position in the grid
        """
        if tree == None:
            return
        newKeys = self.vg.add_tree(tree, offset)
        for key in newKeys:
            # Put percolation index in first place
            key = map(lambda x: x+1, list(key))
            t = key[self.__percDirection__]
            key[self.__percDirection__] = key[0]
            key[0] = t
            self.perc.open(key)
    
    def calc_percolation_value(self, tree, averagingN, maxIterNum = 10**5):
        """
        Calculate percolation value of specified set of trees.
        |  The algorithm is as follows:
        |   we add each tree in the percolation space with random position 
        |   we check if the space percolates
        |   If it does, we stop and record the number of iterations used
        |   If it does not, we continue adding
        
        Parameters:
        -----
        tree : :class:`btmorph.btstructs2.STree2`
        A tree to  be added
        averagingN : int
        Number of obtained percolation numbers to average
        maxIterNum : int
        Maximum number of iterations in one experiment
        details : list
        Peroclation values for every experiment will be stored in this list if not None
        
        Returns:
        -----
        percNumber : number
        Average of percolation numbers of *averagingN* experiments
        """
        if tree == None:
            return None
        percNums = []
        stats = btmorph.BTStats(tree)
        tree_dims = stats.total_dimension()
        percMargin = 1
        for ex in range(0, averagingN):
            self.__init__(self.dim, self.__voxelSz__, None, self.__percDirection__)
            print "Experiment number " + str(ex)
            for i in range(0, maxIterNum):
                r = [0, 0, 0]
                # Random offset
                for j in range(0, len(self.dim)):
                    r[j] = rnd.rand()*self.dim[j]
                r[self.__percDirection__] = 0
                r[self.__percDirection__] = - percMargin*tree_dims[self.__percDirection__] + rnd.rand()*(self.dim[self.__percDirection__] + percMargin*tree_dims[self.__percDirection__])
                self.add_tree(tree, r)
                if self.perc.percolates():
                    percNums.append(i)
                    break
        percNumber = mean(percNums)
        return percNumber, percNums
    
    def plot(self):
        if self.vg == None:
            return
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        keys = self.vg.grid.keys();
        xs = map(lambda (x,y,z): x, keys)
        ys = map(lambda (x,y,z): y, keys)
        zs = map(lambda (x,y,z): z, keys)
        print 'xs', min(xs), max(xs), 'ys', min(ys), max(ys), 'zs', min(zs), max(zs)
        ax.scatter(xs, ys, zs, zdir="z")
        line = [[0]*self.vg.res[0],[0]*self.vg.res[1], [0]*self.vg.res[2]]
        # Top
        for j in range(0, len(line)-1):
            for i in range(0, len(line[j])):
                line[j][i] = i
        line[self.__percDirection__] = [0]*len(line[self.__percDirection__])
        ax.scatter(line[0], line[1], line[2], c='r')
        #Bottom
        line[self.__percDirection__] = [self.vg.res[self.__percDirection__] - 1]*len(line[self.__percDirection__])
        ax.scatter(line[0], line[1], line[2], c='g')        
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        
    
