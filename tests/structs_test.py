'''
Test routines for the btstructs2.py file
'''

import numpy
import random as rndm
from  nose.tools import raises
import math
from nose.tools import with_setup
#sys.path.append('..')

#from btstructs import STree, SNode,P3D
#import btstructs2

import btmorph
from btmorph import STree2
from btmorph import VoxelGrid
from btmorph import BoxCounter

def test_load_swc():
    '''
    Test whether SWC files are correctly loaded
    '''
    swc_tree = btmorph.STree2()
    swc_tree.read_SWC_tree_from_file('tests/v_e_moto1.CNG.swc')
    all_nodes1 = swc_tree.get_nodes()

    print '\nlen(swc_tree1)', len(all_nodes1) 

    assert(len(all_nodes1) == 562)

def test_load_and_write_swc():
    '''
    Test whether SWC treesare correctly written to file
    '''
    swc_tree = STree2()
    swc_tree.read_SWC_tree_from_file('tests/v_e_moto1.CNG.swc')
    swc_tree.write_SWC_tree_to_file('tests/moto_1_outputted.swc')
    swc_tree2 = STree2()
    swc_tree2.read_SWC_tree_from_file('tests/moto_1_outputted.swc')
    print 'len(swc_tree2)', len(swc_tree2.get_nodes()) 
    
    assert(len(swc_tree2.get_nodes()) == 562)
    

if __name__ == '__main__' :
    test_load_swc()

@raises(TypeError)
def VoxelGrid_adjustDim_junk(d, r):
    """
    Test if adjustDimensions raises TypeError if fed with junk
    """
    VoxelGrid.adjustDimensions(d, r)

@raises(IndexError, Exception)
def VoxelGrid_adjustDim_negative(d, r):
    """
    Test if adjustDimensions raises IndexError if fed with negative value
    """
    VoxelGrid.adjustDimensions(d, r)

def test_voxelGrid_adjust_dimensions():
    """
    Test whether x:y:z = rx:ry:rz, where (x,y,z) - output dimensions, (rx,ry,rz) - resolution.
        And none of the original dimensions can be bigger than the output (expansion only)
    Special cases: One dimension and corresponding resolution are zero => no modification needed
                   More than one dimension is zero => return nothing
                   One or more dimension/resolution is zero while corresponding resolution/dimension is not => return nothing
    """   
    epsilon = 0.001
    # Should not accept junk and values less than zero
    junk = [None, [], [4, 't', 4], "str"]
    dimensions = numpy.array([rndm.randint(1000, 2000)/10.0, rndm.randint(1000, 2000)/10.0, rndm.randint(1000, 2000)/10.0]) #[x, y, z]
    resolution = [rndm.randint(1000, 2000), rndm.randint(1000, 2000), rndm.randint(1000, 2000)]    
    for el in junk:
        yield VoxelGrid_adjustDim_junk, dimensions, junk
        yield VoxelGrid_adjustDim_junk, junk, resolution
    for i in range(0,3):
        dimensions = numpy.array([rndm.randint(1000, 2000)/10.0, rndm.randint(1000, 2000)/10.0, rndm.randint(1000, 2000)/10.0]) #[x, y, z]
        resolution = [rndm.randint(1000, 2000), rndm.randint(1000, 2000), rndm.randint(1000, 2000)]    
        resolution[i] = -100
        yield VoxelGrid_adjustDim_negative, dimensions, resolution
        resolution = [rndm.randint(1000, 2000), rndm.randint(1000, 2000), rndm.randint(1000, 2000)]    
        dimensions[i] = -100
        yield VoxelGrid_adjustDim_negative, dimensions, resolution
    # More than one dimension is zero => return nothing
    dimensions = numpy.array([0, 0, 0])
    resolution = [0, 0, 0]
    assert(VoxelGrid.adjustDimensions(dimensions, resolution) == None)
    for i in range(0, 3):
        dimensions = numpy.array([0, 0, 0])
        resolution = [0, 0, 0]
        #One of the dimensions is not zero
        dimensions[i] = rndm.randint(1000, 2000)/10.0
        resolution[i] = rndm.randint(1000, 2000)
        assert(VoxelGrid.adjustDimensions(dimensions, resolution) == None)
    # One or more dimension/resolution is zero while corresponding resolution/dimension is not => return nothing
    # All
    dimensions = numpy.array([0, 0, 0])
    resolution = [rndm.randint(1000, 2000), rndm.randint(1000, 2000), rndm.randint(1000, 2000)]
    assert(VoxelGrid.adjustDimensions(dimensions, resolution) == None)
    dimensions = numpy.array([rndm.randint(1000, 2000)/10.0, rndm.randint(1000, 2000)/10.0, rndm.randint(1000, 2000)/10.0]) #[x, y, z]
    resolution = [0, 0, 0]
    assert(VoxelGrid.adjustDimensions(dimensions, resolution) == None)
    # 1 (2 are taken care of in earlier)
    for i in range(0,3):
        dimensions = numpy.array([rndm.randint(1000, 2000)/10.0, rndm.randint(1000, 2000)/10.0, rndm.randint(1000, 2000)/10.0]) #[x, y, z]
        resolution = [rndm.randint(1000, 2000), rndm.randint(1000, 2000), rndm.randint(1000, 2000)] #[rx, ry, rz]        
        dimensions[i] = 0
        assert(VoxelGrid.adjustDimensions(dimensions, resolution) == None)
        dimensions = numpy.array([rndm.randint(1000, 2000)/10.0, rndm.randint(1000, 2000)/10.0, rndm.randint(1000, 2000)/10.0]) #[x, y, z]
        resolution = [rndm.randint(1000, 2000), rndm.randint(1000, 2000), rndm.randint(1000, 2000)] #[rx, ry, rz]        
        resolution[i] = 0
        assert(VoxelGrid.adjustDimensions(dimensions, resolution) == None)
    # Everything is not zero
    for i in range(0,10):
        # Random floating point from 100 to 200 with 1 digit after comma
        dimensions = numpy.array([rndm.randint(1000, 2000)/10.0, rndm.randint(1000, 2000)/10.0, rndm.randint(1000, 2000)/10.0]) #[x, y, z]
        resolution = [rndm.randint(1000, 2000), rndm.randint(1000, 2000), rndm.randint(1000, 2000)] #[rx, ry, rz]
        [x,y,z] = dimensions
        [rx,ry,rz]=resolution
        [nx,ny,nz] = VoxelGrid.adjustDimensions(dimensions, resolution)
        # None of the original dimensions can be bigger than the output (expansion only)
        assert(x <= nx) #x
        assert(y <= ny) #y
        assert(z <= nz) #z
        # x:y:z =?= rx:ry:rz
        assert(abs(nx/ny - float(rx)/float(ry)) < epsilon)
        assert(abs(ny/nz - float(ry)/float(rz)) < epsilon)    
    
def test_isPowerOfTwo() :
    """
    Tests VoxelGrid.isPowerOfTwo function
    Must return True for 2^m and False otherwise
    """
    assert(VoxelGrid.isPowerOfTwo(-100) == False)
    assert(VoxelGrid.isPowerOfTwo(120.56) == False)
    for i in range(0,32):
        res_true = VoxelGrid.isPowerOfTwo(pow(2, i))
        res_false = VoxelGrid.isPowerOfTwo(pow(2, i)*3)
        assert(res_true == True)
        assert(res_false == False)
        
@raises(TypeError)     
def voxelGrid_check_key_bad_type(key):
    """
    The key must be a tuple of tree integers greater or equal than zero and 
    less than the cossesponding resolution
    Here testing: bad type
    """
    VoxelGrid.checkKey([100,100, 100], key)
    
@raises(IndexError)
def voxelGrid_check_key_out_of_boundary(key):
    """
    The key must be a tuple of tree integers greater or equal than zero and 
    less than the cossesponding resolution
    Here testing: out of boundary
    """
    VoxelGrid.checkKey([200,200,200], key)

def test_gen_voxelgrid_checkKey():
    """
    The key must be a tuple of tree integers greater or equal than zero and 
    less than the cossesponding resolution
    """    
    bad_type = [100.4, ['r', 'r'], "str", None, (10.8, 30, 190), (40, 15.7, 90), (78, 15, 87.9)]
    bad_num = [(100, -1, 100), (-1, 100, 100), (100, 100, -1)]
    bad_index = [(100, 1000, 100), (1000, 100, 100), (100, 100, 1000)]
    good = [(100, 120, 40), (11, 33, 90)]
    #bad type 
    for el in bad_type:
        yield voxelGrid_check_key_bad_type, el
    #bad number
    for el in bad_num:        
        yield voxelGrid_check_key_out_of_boundary, el
    #out of boundary
    for el in bad_index:
        yield voxelGrid_check_key_out_of_boundary, el
    #good
    for el in good:
        assert(VoxelGrid.checkKey([200,200,200], el) == True)

@raises(IndexError)
def VoxelGrid_init_not_pow2(dims, res):
    """
    Test if VoxelGrid initialization raises an exception if not power of two
    """
    VoxelGrid(dims, res)

@raises(TypeError)
def VoxelGrid_init_wrongType(dims, res):
    """
    Test if VoxelGrid initialization raises an exception if fed with junk
    """
    VoxelGrid(dims, res)

def test_VoxelGrid_init():
    """
    Test if VoxelGrid initialization was successfull
    """
    epsilon = 0.001
    #Should not accept junk
    junk = [None, [], [4, 't', 4], "str"]
    for el in junk:
        yield VoxelGrid_init_wrongType, junk, junk
    # Random floating point from 100 to 200 with 1 digit after comma
    dimensions = numpy.array([rndm.randint(100, 200)/10.0, rndm.randint(100, 200)/10.0, rndm.randint(100, 200)/10.0]) #[x, y, z]
    # Should accept only power of two as resolution
    for i in range(0,3):
        resolution = [1024, 1024, 1024] #[rx, ry, rz]
        resolution[i] = 5*1024
        yield VoxelGrid_init_not_pow2, dimensions, resolution        
    resolution = [64, 128, 64] #[rx, ry, rz]
    vg = VoxelGrid(dimensions, resolution)
    good_dim = VoxelGrid.adjustDimensions(dimensions, resolution)
    #check dimensions
    for i in range(0, 3):
        assert(abs(good_dim[i] - vg.dim[i]) < epsilon)
    #check grid - should be all blank
    for i in range(0, resolution[0]):
        for j in range(0, resolution[1]):
            for k in range(0, resolution[2]):
                assert(vg[(i,j,k)] == False)

@raises(TypeError)
def VoxelGrid_setitem_badtype(vg, tpl):
    """
    VoxelGrid.setitem method should raise an error if the type of the key is invalid
    """
    vg[tpl] = True

@raises(TypeError)
def VoxelGrid_setitem_badassign(vg, tpl, val):
    """
    VoxelGrid.setitem method should raise an error if the type of the value to be assigned is not boolean
    """
    vg[tpl] = val

@raises(IndexError)
def VoxelGrid_setitem_badindex(vg, tpl):
    """
    VoxelGrid.setitem method should raise an error if index is out of range
    """
    vg[tpl] = True
    
@raises(TypeError)
def VoxelGrid_getitem_badtype(vg, tpl):
    """
    VoxelGrid.getitem method should raise an error if the type of the key is invalid
    """
    vg[tpl]

@raises(IndexError)
def VoxelGrid_getitem_badindex(vg, tpl):
    """
    VoxelGrid.getitem method should raise an error if index is out of range
    """
    vg[tpl]

def test_VoxelGrid_setitem():
    """
    Test if VoxelGrid setitem method works properly
    """
    # Random floating point from 100 to 200 with 1 digit after comma
    dimensions = numpy.array([rndm.randint(100, 200)/10.0, rndm.randint(100, 200)/10.0, rndm.randint(100, 200)/10.0]) #[x, y, z]
    resolution = [64, 128, 64]
    vg = VoxelGrid(dimensions, resolution)    
    # Bad type
    bad_type = [1, 20.5, (1,2), ["s"], ("1", 2, 3), (1, "2", 3), (1, 2, "3"), "str", None]
    for el in bad_type:
        yield VoxelGrid_setitem_badtype, vg, el
    # bad index
    for i in range(0, 3):
        key = [30, 100, 60]
        key[i] = -40
        yield VoxelGrid_setitem_badindex, vg, tuple(key)
    key = (30, 100, 60)
    # Value to be assigned is not boolean
    for el in bad_type:
        yield VoxelGrid_setitem_badassign, vg, key, el
     # Good example
    vg[key] = True            
    assert(key in vg.grid and vg.grid[key] == True)
    for i in range(0, resolution[0]):
        for j in range(0, resolution[1]):
            for k in range(0, resolution[2]):
                if i == key[0] and j == key[1] and k == key[2]:
                    continue
                assert(not (i,j,k) in vg.grid)
    vg[key] = False
    assert(len(vg.grid) == 0)
    vg[key] = False
    assert(len(vg.grid) == 0)
    
def test_VoxelGrid_getitem():
    """
    Test if VoxelGrid getitem method works properly
    """
    # Random floating point from 100 to 200 with 1 digit after comma
    dimensions = numpy.array([rndm.randint(100, 200)/10.0, rndm.randint(100, 200)/10.0, rndm.randint(100, 200)/10.0]) #[x, y, z]
    resolution = [64, 128, 64]
    vg = VoxelGrid(dimensions, resolution)    
    # Bad type
    bad_type = [1, 20.5, (1,2), ["s"], ("1", 2, 3), (1, "2", 3), (1, 2, "3"), "str", None]
    for el in bad_type:
        yield VoxelGrid_getitem_badtype, vg, el
    # bad index
    for i in range(0, 3):
        key = [30, 100, 60]
        key[i] = -40
        yield VoxelGrid_getitem_badindex, vg, tuple(key)
    # Good example
    key = (30, 100, 60)
    vg.grid[key] = True
    r = vg[key]
    assert(r == True)
    for i in range(0, resolution[0]):
        for j in range(0, resolution[1]):
            for k in range(0, resolution[2]):
                if i == key[0] and j == key[1] and k == key[2]:
                    continue
                assert(vg[(i,j,k)] == False)

def test_VoxelGrid_getsetitem():
    """
    Test if setitem/getitem combination is correct for VoxelGrid
    """
    # Random floating point from 100 to 200 with 1 digit after comma
    dimensions = numpy.array([rndm.randint(100, 200)/10.0, rndm.randint(100, 200)/10.0, rndm.randint(100, 200)/10.0]) #[x, y, z]
    resolution = [64, 128, 64]
    vg = VoxelGrid(dimensions, resolution)    
    key = (30, 100, 60)
    vg[key] = True
    assert(vg[key] == True)

def test_VoxelGrid_plot():
    """
    Plotting the Voxel Grid
    """
    # Random floating point from 100 to 200 with 1 digit after comma
    dimensions = numpy.array([rndm.randint(100, 200)/10.0, rndm.randint(100, 200)/10.0, rndm.randint(100, 200)/10.0]) #[x, y, z]
    resolution = [64, 128, 64**2]
    vg = VoxelGrid(dimensions, resolution)
    for i in range(0, 64):
        for j in range(0, 128):
            vg[(i,j,i**2)] = True
    vg.plot()
    assert(True)

def test_VoxelGrid_calcEncBox_sphere():
    """
    Test if the encompassing box for a sphere is calculated properly
    """
    epsilon = 1.1
    dimensions = numpy.array([64, 128, 64])
    resolution = [2*64, 2*128, 2*64]
    vg = VoxelGrid(dimensions, resolution)
    # If radius is less than zero => return None
    center = (dimensions[0]/2, dimensions[1]/2, dimensions[2]/2)
    radius = -10
    assert(vg.calcEncompassingBox_sphere(center, radius) == None)
    # If radius is zero => one point
    radius = 0
    res = vg.calcEncompassingBox_sphere(center, radius)
    assert(res != None)
    [(x1,x2), (y1,y2), (z1,z2)] = res
    assert(x1 == x2 == resolution[0]/2)
    assert(y1 == y2 == resolution[1]/2)
    assert(z1 == z2 == resolution[2]/2)
    # If the sphere is completely out of the grid, should return None
    center = (1000, 1000, 1000)
    radius = 100
    assert(vg.calcEncompassingBox_sphere(center, radius) == None)
    # If the sphere is completely inside the grid, should return:
    # x: [(x0-r)*rx/dx, (x0+r)*rx/dx], analog. for y and z
    center = (dimensions[0]/2, dimensions[1]/2, dimensions[2]/2) # <- the very center
    radius = dimensions[0]/4 # <- well inside the grid
    res = vg.calcEncompassingBox_sphere(center, radius)
    assert(res != None)
    [(x1,x2), (y1,y2), (z1,z2)] = res
    assert(abs(x1 - 2*(center[0] - radius)) < epsilon) # Resolution is twice the corr. dimension
    assert(abs(x2 - 2*(center[0] + radius)) < epsilon)
    assert(abs(y1 - 2*(center[1] - radius)) < epsilon)
    assert(abs(y2 - 2*(center[1] + radius)) < epsilon)
    assert(abs(z1 - 2*(center[2] - radius)) < epsilon)
    assert(abs(z2 - 2*(center[2] + radius)) < epsilon)
    # Partial intersection cases
    center1 = (0, 0, 0)
    center2 = (dimensions[0], dimensions[1], dimensions[2])
    res = vg.calcEncompassingBox_sphere(center1, radius)
    assert(res != None)    
    [(x1,x2), (y1,y2), (z1,z2)] = res
    assert(abs(x1 - 0) < epsilon)
    assert(abs(x2 - 2*(center1[0] + radius)) < epsilon)
    assert(abs(y1 - 0) < epsilon)
    assert(abs(y2 - 2*(center1[1] + radius)) < epsilon)
    assert(abs(z1 - 0) < epsilon)
    assert(abs(z2 - 2*(center1[2] + radius)) < epsilon)
    res = vg.calcEncompassingBox_sphere(center2, radius)
    assert(res != None)    
    [(x1,x2), (y1,y2), (z1,z2)] = res
    assert(abs(x1 - 2*(center2[0] - radius)) < epsilon)
    assert(abs(x2 - resolution[0]) < epsilon)
    assert(abs(y1 - 2*(center2[1] - radius)) < epsilon)
    assert(abs(y2 - resolution[1]) < epsilon)
    assert(abs(z1 - 2*(center2[2] - radius)) < epsilon)
    assert(abs(z2 - resolution[2]) < epsilon)

def test_VoxelGrid_fallsIntoSphere():
    """
    Test VoxelGrid.fallsIntoSphere
    """
    dimensions = numpy.array([64, 128, 64])
    resolution = [2*64, 2*128, 2*64]
    vg = VoxelGrid(dimensions, resolution)
    center = (dimensions[0]/2, dimensions[1]/2, dimensions[2]/2) # <- the very center
    radius = dimensions[0]/4 # <- well inside the grid
    not_inside = [(0,0,0), (-100,0,0), (2*resolution[0], 0, 0)]
    voxel_center = (resolution[0]/2, resolution[1]/2, resolution[2]/2)
    inside = [voxel_center, (resolution[0]/2 + resolution[0]/4 -3, voxel_center[1], voxel_center[2])]
    # If radius is < 0 => False
    assert(vg.fallsIntoSphere((1,1,1), center, -100) == False)
    # If radius is == 0 => only the center
    assert(vg.fallsIntoSphere(voxel_center, center, 0) == True)
    assert(vg.fallsIntoSphere((voxel_center[0], voxel_center[1]+1, voxel_center[2]), center, 0) == False)
    # Radius > 0
    for el in inside:
        assert(vg.fallsIntoSphere(el, center, radius) == True)
    for el in not_inside:
        assert(vg.fallsIntoSphere(el, center, radius) == False)


def test_VoxelGrid_fallsIntoFrustum():
    """
    Test VoxelGrid.fallsIntoFrustum
    """
    dimensions = numpy.array([64, 128, 64])
    resolution = [2*64, 2*128, 2*64]
    vg = VoxelGrid(dimensions, resolution)
    # A frustum completely inside the grid
    c1 = (dimensions[0]/4, dimensions[1]/4, dimensions[2]/4)
    c2 = (dimensions[0]/3, dimensions[1]/3, dimensions[2]/3)
    c1_v = vg.dimensionToVoxel(((c1[0]+c2[0])/2.0, (c1[1]+c2[1])/2.0, (c1[2]+c2[2])/2.0))
    #c2_v = (resolution[0]/2, resolution[1]/2, resolution[2]/2)
    r1 = 5
    r2 = 10
    not_inside = [(0,0,0), (-100,0,0), (2*resolution[0], 0, 0)]
    inside = [c1_v, vg.dimensionToVoxel(c1), vg.dimensionToVoxel(c2)]
    # Normal case
    for el in inside:
        assert(vg.fallsIntoFrustum(el, c1, r1, c2, r2) == True)
    for el in not_inside:
        assert(vg.fallsIntoFrustum(el, c1, r1, c2, r2) == False)
    #c1 == c2
    not_inside = [(0,0,0), (-100,0,0), (2*resolution[0], 0, 0), vg.dimensionToVoxel(c2), c1_v]
    inside = [vg.dimensionToVoxel(c1)]
    for el in inside:
        assert(vg.fallsIntoFrustum(el, c1, r1, c1, r2) == True)
    for el in not_inside:
        assert(vg.fallsIntoFrustum(el, c1, r1, c1, r2) == False)

def test_VoxelGrid_addFrustum():
    """
    Test if a frustum is added properly
    """
    dimensions = numpy.array([64, 128, 64])    
    resolution = [2*64, 2*128, 2*64]
    vg = VoxelGrid(dimensions, resolution)    
    # A frustum completely inside the grid
    c1 = (dimensions[0]/4, dimensions[1]/4, dimensions[2]/4)
    c2 = (dimensions[0]/3, dimensions[1]/3, dimensions[2]/3)
    h = math.sqrt((c2[0] - c1[0])**2 + (c2[1] - c1[1])**2 + (c2[2] - c1[2])**2)
    r1 = 10.0
    r2 = 5.0
    # If one of the radii < zero => nothing is added
    vg.addFrustum(c1, -r1, c2, r2)
    assert(len(vg.grid) == 0)
    vg.addFrustum(c1, r1, c2, -r2)
    assert(len(vg.grid) == 0)
    # if c2 = c1 and both of the radii == 0 => only one point is added
    # Reset
    vg = VoxelGrid(dimensions, resolution)
    vg.addFrustum(c1, 0, c1, 0)
    assert(len(vg.grid) == 1)
    assert(vg[vg.dimensionToVoxel(c1)] == True)
    # If frustum is too far => nothing is added
    vg = VoxelGrid(dimensions, resolution)
    vg.addFrustum((-1000, -1000, -1000), r1, (-500, -900, -100), r2)
    assert(len(vg.grid) == 0)
    # if both of the radii == 0 => only central axis of height h is added
    vg.addFrustum(c1, 0 , c2, 0)
    assert(vg[vg.dimensionToVoxel(c1)] == True)
    assert(vg[vg.dimensionToVoxel(c2)] == True)
    assert(len(vg.grid) >= abs(vg.dimensionToVoxel(c2)[2]-vg.dimensionToVoxel(c1)[2]) and len(vg.grid) < 3*abs(vg.dimensionToVoxel(c2)[2]-vg.dimensionToVoxel(c1)[2]))
    # A frustum completely inside the grid
    epsilon = 0.02
    vg = VoxelGrid(dimensions, resolution)
    V_f = math.pi*h*(r1**2 + r1*r2 + r2**2)/3.0
    for i in range(1, 6):
        resolution = [(2**i)*64, (2**i)*128, (2**i)*64]
        voxel_vol = (dimensions[0]/float(resolution[0])) * (dimensions[1]/float(resolution[1])) * (dimensions[2]/float(resolution[2]))
        vg = VoxelGrid(dimensions, resolution)
        vg.addFrustum(c1, r1, c2, r2)
       # print ("diff",len(vg.grid), voxel_vol , V_f, len(vg.grid)*voxel_vol, epsilon*V_f)
        if abs(len(vg.grid)*voxel_vol - V_f) < epsilon*V_f:
            assert(True)
            return
    assert(False)
    
def test_VoxelGrid_addSphere():
    """
    Test if a sphere is added properly
    """
    dimensions = numpy.array([64, 128, 64])    
    center = (dimensions[0]/2, dimensions[1]/2, dimensions[2]/2)
    resolution = [2*64, 2*128, 2*64]
    vg = VoxelGrid(dimensions, resolution)
    voxel_center = (resolution[0]/2, resolution[1]/2, resolution[2]/2)
    # If radius < zero => nothing is added
    vg.addSphere(center, -100)
    assert(len(vg.grid) == 0)
    # If radius == 0 => only one point is added
    vg = VoxelGrid(dimensions, resolution)
    vg.addSphere(center, 0)
    assert(len(vg.grid) == 1)
    assert(vg[voxel_center] == True)
    # If sphere is too far => nothing is added
    vg = VoxelGrid(dimensions, resolution)
    center = (1000, 1000, 1000)
    radius = 100
    vg.addSphere(center, radius)
    assert(len(vg.grid) == 0)
    # Sphere is completely inside the grid
    epsilon = 0.01
    center = (dimensions[0]/2, dimensions[1]/2, dimensions[2]/2) # <- the very center
    radius = 20.0 # <- well inside the grid    
    V_r = (4*numpy.pi*radius**3)/3.0
    for i in range(1, 6):
        resolution = [(2**i)*64, (2**i)*128, (2**i)*64]
        voxel_vol = (dimensions[0]/float(resolution[0])) * (dimensions[1]/float(resolution[1])) * (dimensions[2]/float(resolution[2]))
        vg = VoxelGrid(dimensions, resolution)
        vg.addSphere(center, radius)
        #print ("diff",len(vg.grid), voxel_vol , V_r, len(vg.grid)*voxel_vol)
        if abs(len(vg.grid)*voxel_vol - V_r) < epsilon*V_r:
            assert(True)
            return
    assert(False)
        

def test_VoxelGrid_calcEncBox_frustum():
    """
    Test if encompassing box of a frustum is calculated properly
    """
    dimensions = numpy.array([64, 128, 64])
    resolution = [64, 128, 64]
    vg = VoxelGrid(dimensions, resolution)
    # If one of the parameters is None then must return None
    c1 = (2,1,1)
    c2 = (1, 1,1)
    r1 = 1
    r2 = 2
    assert(vg.calcEncompassingBox_frustum(None, r1, c2, r2) == None)
    assert(vg.calcEncompassingBox_frustum(c1, None, c2, r2) == None)
    assert(vg.calcEncompassingBox_frustum(c1, r1, None, r2) == None)
    assert(vg.calcEncompassingBox_frustum(c1, r1, c2, None) == None)
    # if r1 or r2 is less than zero => None
    r1 = -1
    r2 = 2
    assert(vg.calcEncompassingBox_frustum(c1, r1, c2, r2) == None)
    assert(vg.calcEncompassingBox_frustum(c1, r2, c2, r1) == None) 
    # A frustum completely inside the grid
    c1 = (dimensions[0]/4, dimensions[1]/4, dimensions[2]/4)
    c2 = (dimensions[0]/3, dimensions[1]/3, dimensions[2]/3)
    r1 = 5
    r2 = 10
    res = vg.calcEncompassingBox_frustum(c1, r1, c2, r2)
    assert(res != None)
    [(x1,x2), (y1,y2), (z1,z2)] = res
    left = [x1, y1, z1]
    right = [x2, y2, z2]
    for i in range(0,3):        
        if c1[i] - r1 < 0 or c2[i] - r2 < 0:
            assert(left[i] == 0)
        else:
            assert(left[i] <= round((c1[i] - r1)*vg.res[i]/vg.dim[i]) and left[i] <= round((c2[i] - r2)*vg.res[i]/vg.dim[i]))
    for i in range(0,3):
        if c1[i] + r1 > vg.dim[i] or c2[i] + r2 > vg.dim[i]:
            assert(right[i] == vg.res[i])
        else:
            assert(right[i] >= round((c1[i]+r1)*vg.res[i]/vg.dim[i]) and right[i] >= round((c2[i]+r2)*vg.res[i]/vg.dim[i]))
    #A  frustum partially intersecting the grid
    # A frustum partially intersecting the grid will be properly dealt with only after rotation
    c1 = (-dimensions[0]/4, dimensions[1]/4, dimensions[2]/4)
    c2 = (dimensions[0]/2, dimensions[1]/2, dimensions[2]/2)
    r1 = 5
    r2 = 10
    res = vg.calcEncompassingBox_frustum(c1, r1, c2, r2)
    assert(res != None)
    [(x1,x2), (y1,y2), (z1,z2)] = res  
    left = [x1, y1, z1]
    right = [x2, y2, z2]
    for i in range(0,3):        
        if c1[i] - r1 < 0 or c2[i] - r2 < 0:
            assert(left[i] == 0)
        else:
            assert(left[i] <= round((c1[i] - r1)*vg.res[i]/vg.dim[i]) and left[i] <= round((c2[i] - r2)*vg.res[i]/vg.dim[i]))
    for i in range(0,3):
        if c1[i] + r1 > vg.dim[i] or c2[i] + r2 > vg.dim[i]:
            assert(right[i] == vg.res[i])
        else:
            assert(right[i] >= round((c1[i]+r1)*vg.res[i]/vg.dim[i]) and right[i] >= round((c2[i]+r2)*vg.res[i]/vg.dim[i]))
    
""" New fucntions by Irina - test"""
test_trees = []
test_stats = []
    
def setup_func_small_tree():
    """
    Setup function for tree initialization and loading
    """
    global test_trees
    global test_stats
    #0 - Only soma tree
    #test_trees.append(btmorph.STree2().read_SWC_tree_from_file("tests/soma_only.swc")) 
    #1 - Wiki test tree moto_1_outputted
    #test_trees.append(btmorph.STree2().read_SWC_tree_from_file("tests/horton-strahler_test_wiki_3pointsoma.swc"))    
    test_trees.append(btmorph.STree2().read_SWC_tree_from_file("tests/moto_1_outputted.swc"))        
    test_stats = [btmorph.BTStats(test_trees[0])]

def teardown_func_small_tree():
    """
    Teardown function for tree initialization and loading
    """
    global test_trees
    global test_stats
    test_trees = []
    test_stats = []

@with_setup(setup_func_small_tree, teardown_func_small_tree) 
def test_VoxelizeTree():
    """
    Crudely test if a tree is voxelized properly
    """
    global test_trees
    global test_stats
    res = [256,128,256]
    for i in range(0, len(test_trees)):
        dx,dy,dz = test_stats[i].total_dimension()
        dims = [dx,dy,dz]
        print(dims, res)
        vg = VoxelGrid(dims, res)
        vg.addTree(test_trees[i])
        print vg
        print "stats volume:" + str(test_stats[i].total_volume()[0])
        vg.plot()
        assert(len(vg.grid) > 0)
        
@with_setup(setup_func_small_tree, teardown_func_small_tree) 
def test_BoxCounter_init():
    """
    Test if BoxCounter is initialized properly
    """
    global test_trees
    global test_stats
    res = [256,128,256]
    dx,dy,dz = test_stats[0].total_dimension()
    dims = [dx,dy,dz]
    vg = VoxelGrid(dims, res, test_trees[0])
    bc = BoxCounter(vg)
    assert(len(bc.countVals) == 8)

@with_setup(setup_func_small_tree, teardown_func_small_tree) 
def test_gridCount():
    """
    Test if standard box counting works properly
    """
    # if startDim < 0 => return -1
    startDim = -20
    global test_trees
    global test_stats
    res = [64,32,64]
    dx,dy,dz = test_stats[0].total_dimension()
    dims = [dx,dy,dz]
    vg = VoxelGrid(dims, res, test_trees[0])
    bc = BoxCounter(vg)
    assert(bc.gridCount(startDim) == -1)
    # if not power of two => -1
    startDim = 14
    assert(bc.gridCount(startDim) == -1)
    # if > smallest dim => -1
    assert(bc.gridCount(res[0]) == -1)
    # Normal case
    # All sums must be the same
    bc.gridCount(res[1])
    s_t = len(bc.vg.grid)
    print(bc.vg)
    for i in range(1,len(bc.countVals)):
       s = sum(bc.countVals[i])
       print(s_t, s)
       assert(s == s_t)

def test_coverage_and_count():
    """
    Test coverage count and box count together
    """
    res = [32, 32, 32]
    dim = [32.0, 32.0, 32.0]
    vg = VoxelGrid(dim, res)
    bc = BoxCounter(vg)
    # Generate solid box inside the grid
    a = 16
    for i in range(8,8+a):
        for j in range(8, 8+ a):
            for k in range(8, 8 + a):
                bc.vg[(i,j,k)] = True
    startDim = res[0]/2
    bc.gridCoverage(startDim)
    bc.gridCount(startDim)
    for i in range(1,len(bc.countVals)):
        s = sum(1 for e in bc.countVals[i] if e)
        assert(s == bc.coverageVals[i])
       
def test_coverageCount():
    """
    Test if gridCoverage method in BoxCounter works properly
    """     
    res = [32, 32, 32]
    dim = [32.0, 32.0, 32.0]
    vg = VoxelGrid(dim, res)
    bc = BoxCounter(vg)
    # if startDim < 0 => return -1
    startDim = -20
    assert(bc.gridCount(startDim) == -1)
    # if not power of two => -1
    startDim = 14
    assert(bc.gridCount(startDim) == -1)
    # if > smallest dim => -1
    assert(bc.gridCount(res[0]*2) == -1)
    # Generate solid box inside the grid
    a = 16
    for i in range(8,8+a):
        for j in range(8, 8+ a):
            for k in range(8, 8 + a):
                bc.vg[(i,j,k)] = True
    startDim = res[0]/2
    bc.gridCoverage(startDim)
    print(bc.coverageVals)
    assert(bc.coverageVals[1] == 8**3)
    assert(bc.coverageVals[2] == 4**3)
    assert(bc.coverageVals[3] == 2**3)
    assert(bc.coverageVals[4] == 2**3)