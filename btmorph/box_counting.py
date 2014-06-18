"""
3D box counting method implementation
Contains
    BoxCounter
    
by Irina Reshodko
"""

from btstructs2 import VoxelGrid
import math

class BoxCounter:
    """
    3D box counting functionality
    """
    def __init__(self, vg):
        """
        Initialize
        
        Parameters
        -------
        vg : :class:'btmorph.btstructs2.VoxelGrid'
            Ready to use voxel grid
        """
        # Resolution is power of two
        self.vg = vg
        [dx, dy, dz] = vg.res
        nSizes = int(math.log(min(vg.res), 2)) + 1
        self.countVals = [[]]*nSizes
        self.coverageVals = [0]*nSizes
        for i in range(0, len(self.countVals)):
            self.countVals[i] = []
            
    def gridCoverage(self, startDim, coords = None):
        """
        Box counting for coverage (standard grid method)
        Works faster than gridCount
        
        Parameters
        ----------
        startDim : int
        Dimension of the box in voxels
        coords : tuple of ints
        Grid position
        """
        if startDim < 1:
            return -1
        # Should be a power of two
        if not VoxelGrid.isPowerOfTwo(startDim):
            return -1
        # startDim should not be greater than smallest dimension
        if startDim > min(self.vg.res):
            return -1
        if coords != None:
            if 1 == startDim:
                if True == self.vg[coords]:
                    return True
                else:
                    return False
            else:
                m = int(math.log(startDim, 2))
                newDim = startDim/2
                new_c = [coords[0]*2, coords[1]*2, coords[2]*2]
                tmp = False
                for i in [new_c[0], new_c[0] + 1]:
                    for j in [new_c[1], new_c[1] + 1]:
                        for k in [new_c[2], new_c[2] + 1]:
                            tmp = self.gridCoverage(newDim, (i,j,k)) or tmp
                self.coverageVals[m] += tmp
                return tmp
        else:
            dx = self.vg.res[0]/startDim
            dy = self.vg.res[1]/startDim
            dz = self.vg.res[2]/startDim
            for i in range(0, dx):
                for j in range(0, dy):
                    for k in range(0, dz):                        
                        self.gridCoverage(startDim, (i, j, k))
    
    def gridCount(self, startDim, coords = None):
        """
        Box counting for grid method (standard)
        
        Parameters
        ----------
        startDim : int
        Dimension of the box in voxels
        coords : tuple of ints
        Grid position
        """
        if startDim < 1:
            return -1
        # Should be a power of two
        if not VoxelGrid.isPowerOfTwo(startDim):
            return -1
        # startDim should not be greater than smallest dimension
        if startDim > min(self.vg.res):
            return -1
        if coords != None:
            if 1 == startDim:
                if True == self.vg[coords]:
                    return 1
                else:
                    return 0
            else:
                m = int(math.log(startDim, 2))
                newDim = startDim/2
                new_c = [coords[0]*2, coords[1]*2, coords[2]*2]
                s = 0
                for i in [new_c[0], new_c[0] + 1]:
                    for j in [new_c[1], new_c[1] + 1]:
                        for k in [new_c[2], new_c[2] + 1]:
                            s += self.gridCount(newDim, (i,j,k))
                self.countVals[m].append(s)
                return s
        else:
            dx = self.vg.res[0]/startDim
            dy = self.vg.res[1]/startDim
            dz = self.vg.res[2]/startDim
            for i in range(0, dx):
                for j in range(0, dy):
                    for k in range(0, dz):                        
                        self.gridCount(startDim, (i, j, k))
        
