"""
Tests for percolation feature
Author: Irina Reshodko
"""
from btmorph import Percolation
from  nose.tools import raises

def test_open_normal():
    """
    Test if open procedure works fine
    |  after open(coord) isOpen(coord) must be true;
    |  No other sites should change
    |  Should throw IndexError if coordinates are out of range
    """    
    N = 10
    # Constructor should be tested in other test
    p = Percolation((N,N,N))
    # Normal 
    # after open(coord) isOpen(coord) must be true;
    # No other sites should change
    coord = (5,5,5)
    p.open(coord)
    assert(p.is_open(coord))
    for i in range(1, N+1):
        for j in range(1, N+1):
            for k in range(1, N+1):
                if (i, j, k) != coord:
                    assert(not p.is_open((i, j ,k)))
    # Erroneous
    errIndex = [-1, 0, N+1]
    for i in errIndex:
        for j in errIndex:
            for k in errIndex:
                yield bad_index_open_err, p, (i, j, k)
    # Boundary cases
    # N = 1
    p = Percolation((1,1,1)) 
    p.open((1,1,1))
    assert(p.is_open((1,1,1)))
    # N = 2
    N = 2
    p = Percolation((N,N,N))
    coord = (1,1,1)
    p.open(coord)
    for i in range(1, N+1):
        for j in range(1, N+1):
            for k in range(1, N+1):
                if (i, j, k) != coord:
                    assert(not p.is_open((i, j ,k)))
 
@raises(IndexError)                   
def bad_index_open_err(p, coord):
    """
    Test if open procedure works fine
    |  Should throw IndexError if coordinates are out of range
    """
    p.is_open(coord)

def test_open_connections():
    """
    Test if opened sites are connected properly
    |  Test via isFull
    """
    N = 10
    sz = (N, N)
    p = Percolation(sz)
    # 11*222*0***
    # 1***2*000**
    # **3****0***
    # *333*******
    coord = (1, 1) # 1
    p.open(coord)
    p.open((coord[0], coord[1] + 1))
    p.open((coord[0] + 1, coord[1]))
    assert(p.is_full(coord))
    assert(p.is_full((coord[0], coord[1] + 1)))
    assert(p.is_full((coord[0] + 1, coord[1])))
    p.open((coord[0] + 1, coord[1] + 2))
    assert(not p.is_full((coord[0] + 1, coord[1] + 3)))
    coord = (1, 6) # 2
    p.open(coord)
    p.open((coord[0], coord[1] + 1))
    p.open((coord[0], coord[1] - 1))
    p.open((coord[0]+1, coord[1]))
    assert(p.is_full(coord))
    assert(p.is_full((coord[0], coord[1] + 1)))
    assert(p.is_full((coord[0], coord[1] - 1)))
    assert(p.is_full((coord[0]+1, coord[1])))
    coord = (2, 9) # 0
    p.open(coord)
    p.open((coord[0], coord[1] + 1))
    p.open((coord[0], coord[1] - 1))
    p.open((coord[0]+1, coord[1]))
    p.open((coord[0]-1, coord[1]))
    assert(p.is_full(coord))
    assert(p.is_full((coord[0], coord[1] + 1)))
    assert(p.is_full((coord[0], coord[1] - 1)))
    assert(p.is_full((coord[0]+1, coord[1])))
    assert(p.is_full((coord[0]-1, coord[1])))
    # Bounding cases
    # N = 1 3D
    N = 1
    sz = (N, N, N)
    coord = (1, 1, 1)
    p = Percolation(sz)
    p.open(coord)
    assert(p.is_full(coord))
    # N = 2 2D
    N = 3
    coord = (1,1)
    sz = (N, N)
    p = Percolation(sz)
    p.open(coord)
    assert(p.is_full(coord))
    p.open((3,3))
    assert(not p.is_full((3,3)))
    p.open((2,2))
    assert(p.is_full((2,2)))
    assert(not p.is_full((2,1)))
    assert(p.is_full((3,3)))
    # N = 2 3D
    sz = (N, N, N)
    p = Percolation(sz)
    p.open((1,1,1))
    p.is_full((1,1,1))
    p.open((3,3,3))
    assert(not p.is_full((3,3,3)))
    p.open((1,2,1))
    p.open((1,2,2))
    p.open((2,2,2))
    assert(p.is_full((2,2,2)))
    assert(p.is_full((3,3,3)))

@raises(IndexError)
def wrong_size_init(sz):
    """
    Test if constructor works fine
    |  Should throw IndexError if N <= 0
    """
    p = Percolation(sz)
    
def test_init():
    """
    Test if constructor works fine
    |  is_open should be false for all sites
    |  Should throw IndexError if N <= 0
    """
    N = 10
    # Normal
    p = Percolation((N, N, N))
    for i in range(1, N+1):
        for j in range(1, N+1):
            for k in range(1, N+1):
                assert(not p.is_open((i,j,k)))
    bad_sz = [(-1, 1, 1), (1, -1, 1), (1, 1, -1), (0, 1, 1), (1, 0, 1), (1, 1, 0)]
    for sz in bad_sz:
        yield wrong_size_init, sz
        
def test_ND_to_1D():
    """
    Test if index conversion works properly
    """
    N = 10
    p = Percolation((N, N, N))
    a = [0]*N**3
    for i in range(0,N):
        for j in range(0,N):
            for k in range(0,N):
                ind = p.to_1D((i,j,k))
                a[ind] = 1
    assert(sum(a) == len(a))

def test_percolates():
    """
    Test if percolates function works correctly
    |  Should return true if there are open sites on column j for all i; false otherwise
    """   
    for N in range(1, 8):
        for j in range(1, N+1):
            for k in range(1, N+1):
                # Percolates
                p = Percolation((N, N, N))
                for i in range(1, N+1):
                    p.open((i, j, k))
                assert(p.percolates())
                p = Percolation((N, N))
                for i in range(1, N+1):
                    p.open((i, j))
                    assert(p.is_full((i,j)))
                assert(p.percolates())
                # Does not percolate
                p = Percolation((N, N, N))
                for i in range(1, N):
                    p.open((i, j, k))
                assert(not p.percolates())
                p = Percolation((N, N, N))
                for i in range(2, N+1):
                    p.open((i, j, k))
                assert(not p.percolates())

def test_is_full():
    """
    Test if is_full works properly
    """
    for N in range(1, 8):
        for j in range(1, N+1):
            for k in range(1, N+1):
                # Full
                p = Percolation((N, N, N))
                for i in range(1, N+1):
                    p.open((i, j, k))
                    assert(p.is_full((i,j,k)))
                p = Percolation((N, N))
                for i in range(1, N+1):
                    p.open((i, j))
                    assert(p.is_full((i,j)))
                # Not full
                p = Percolation((N, N, N))
                for i in range(2, N):
                    p.open((i, j, k))
                    assert(not p.is_full((i,j,k)))
                p.open((1, j, k))
                for i in range(1, N):
                    assert(p.is_full((i,j,k)))
                p = Percolation((N, N, N))
                for i in range(2, N+1):
                    p.open((i, j, k))