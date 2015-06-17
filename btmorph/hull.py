#import time
import numpy as np
#import matplotlib.pyplot as plt
import scipy
from scipy.spatial import ConvexHull as CH

def get_convex_hull_from_tree(tree):
    """
    Compute the convex hull around some 2D collection of points

    Parameters
    ----------
    tree : :class:`STree2`
        Neuronal tree for which to compute the convex hull.

        .. warning:: Works correctly when the tree is projected along
           the PCA axes. Use :py:func:`btmorph.btviz.pca_translate_tree`

    Returns
    -------
    hull_points : array-like
        List of points on the convex hull
    
    """
    nodes = tree.get_nodes()
    N = len(nodes)
    coords = map(lambda n: n.content['p3d'].xyz, nodes)
    points = np.transpose(coords)
    
    xz_pts  = np.array([points[0,:],points[2,:]]).T
    yz_pts  = np.array([points[1,:],points[2,:]]).T
    
    xz_hull = CH(xz_pts)
    yz_hull = CH(yz_pts)

    xz_h_pts = np.zeros((len(xz_hull.vertices)+1,2))
    for i,point in zip(range(len(xz_hull.vertices)),xz_hull.vertices):
        xz_h_pts[i,0] = xz_pts[point][0]
        xz_h_pts[i,1] = xz_pts[point][1]
    xz_h_pts[-1,0] = xz_h_pts[0,0] # first point again to close the polygon
    xz_h_pts[-1,1] = xz_h_pts[0,1]

    yz_h_pts = np.zeros((len(yz_hull.vertices)+1,2))
    for i,point in zip(range(len(yz_hull.vertices)),yz_hull.vertices):
        yz_h_pts[i,0] = yz_pts[point][0]
        yz_h_pts[i,1] = yz_pts[point][1]
    yz_h_pts[-1,0] = yz_h_pts[0,0] # first point again to close the polygon
    yz_h_pts[-1,1] = yz_h_pts[0,1]        
    
    return xz_h_pts,yz_h_pts,xz_pts,yz_pts

def polygon_area(corners):
    """
    Compute the area of convex polygons according to the Shoelace algorthm
    Code from: http://stackoverflow.com/questions/24467972/calculate-area-of-polygon-given-x-y-coordinates

    Parameters
    ----------
    corners: array-like:
        list of 2D points of the vertices (array should be shaped: Nx2).
        .. warning:: Corners must be ordered counter-clockwise

    Returns
    -------
    area : total area in the same units as the polygon's vertices
    """
    x = corners[:,0]
    y=corners[:,1]
    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))
