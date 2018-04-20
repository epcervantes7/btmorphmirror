import sys
import sys
import btmorph
import btmorph.btstructs2
import copy

sys.setrecursionlimit(10000)
import numpy as np
import math
import scipy
import scipy.optimize
import os
from shutil import copyfile
import matplotlib.pyplot as plt
import pylab as p, time
from grids import VoxelGrid
from box_counting import BoxCounter
# from numpy import mean,cov,double,cumsum,dot,linalg,array,rank
from numpy import mean, cov, cumsum, dot, linalg, array, rank
import copy
import math
# import btmorph
from pylab import plot, subplot, axis, stem, show, figure
import Queue


class BTStats(object):
    '''
    Compute morphometric features and statistics of a single morphology

    Assume the "3 point" soma of the curated NeuroMorpho format. (`website <http://neuromorpho.org/neuroMorpho/SomaFormat.html>`_)

    B. Torben-Nielsen (legacy code)
    '''

    def __init__(self, tree):
        """
        Constructor.

        Parameters
        -----------
        tree : :class:`STree2`
            Neuronal tree for which to compute morphometrics
        """
        self._tree = tree
        self._all_nodes = self._tree.get_nodes()

        # compute some of the most used stats +
        self._soma_points, self._bif_points, self._end_points = \
            self.get_points_of_interest()

    def get_end_branches(self):
        """
        Get lists containting the "terminal branches", i.e., from end point to first bifurcation \

        Returns
        -------
        end_branches : list

         """
        end_branches = []

        for node in self._end_points:
            n = self._tree.path_to_bifurcation(node)
            # if not node.index in (1,2,3) :
            end_branches.append(n)
        return end_branches

    def get_sub_tree(self):
        """
        Find and return the nodes of all branches


        Returns
        -------
        path : list of :class:`STree`
             .
        """
        sub_tree = []
        sub_tree.append(self._tree)
        for node in self._bif_points:
            for child in node.children:
                sub_tree.append(self._tree.get_sub_tree(child))
        return sub_tree

    def get_main_trees(self):
        trees = []
        for node in self._tree.root.children:
            # if len(node.children) > 1  :
            #   if not node.parent is None :
            # branches.append(node) # the root is not a bifurcation
            trees.append(self._tree.get_sub_tree(node))

        return trees

    def get_brain_branches(self):
        children = []
        if self._tree.root != None:
            for node in self._tree.root.children:
                tmp = self._tree

        return children

    def write_main_branches(self, path, neuronName, class_number):

        trees = []

        if self._tree.root != None:
            for node in self._tree.root.children:
                if len(node.children) > 0 & node.index != 2:
                    tmp_tree = btmorph.STree2()
                    tmp_tree.root = copy.copy(self._tree.root)
                    tmp_tree.root.children = [node]
                    trees.append(tmp_tree)
            file_num = 1

            for sub_tree in trees:
                # print (len(sub_tree.root.children))

                # sub_tree_name="./"+path+"/"+neuronName+"_branch"
                # if not os.path.exists(sub_tree_name):
                #     os.makedirs(sub_tree_name)
                # sub_tree_name="./"+path+"/"+neuronName+"_branch"
                store_path = "./main_branches/" + path
                if not os.path.exists(store_path):
                    os.makedirs(store_path)
                # n=sub_tree.other_write_SWC_tree_to_file(str(store_path) + '/' + str(class_number) + '__' + neuronName + '_' + str(file_num) + '.swc')

                # if len(sub_tree.get_nodes()) > 1:
                #     # print(str(sub_tree_name)+str(class_number)+'__'+neuronName+'_'+
                #     #                                       new_file_name+'.swc')
                n = sub_tree.other_write_SWC_tree_to_file(
                    str(store_path) + '/' + str(class_number) + '__' + neuronName + '_' + str(file_num) + '.swc')
                # else:
                #     sub_tree = None
                file_num += 1
                # copyfile(path+"/"+neuronName , sub_tree_name+'/'+neuronName)
        return trees

    def write_subtree_by_level(self, level, path, neuronName):
        sub_tree_name = "./" + path + "/" + str(level) + "/"
        # sub_tree_name="./"+path+"/"
        if not os.path.exists(sub_tree_name):
            os.makedirs(sub_tree_name)
        sub_tree = self.get_tree_by_level(level)
        sub_tree.other_write_SWC_tree_to_file(str(sub_tree_name) + neuronName + "_" + str(level) + '.swc')
        return sub_tree

    def write_subtree_by_range(self, from_level, to_level, path, neuronName, class_number):
        """
        This method gets slices from an specified level from_level to the terminals of the neuron.
        :rtype: subtree
        """
        # new_file_name = str(from_level)
        new_file_name = str(to_level)
        # sub_tree_name="./"+path+"/inv"+new_file_name+"/"
        sub_tree_name = "./" + path + "/acu" + new_file_name + "/"
        # sub_tree_name="./"+path+"/viz/"
        if not os.path.exists(sub_tree_name):
            os.makedirs(sub_tree_name)

        if from_level == 0:
            sub_tree = self._tree
        else:
            sub_tree = self.get_tree_by_range(from_level, to_level)
        if (len(sub_tree.get_nodes()) > 1):
            sub_tree.other_write_SWC_tree_to_file(
                str(sub_tree_name) + str(class_number) + '__' + neuronName + '_' + new_file_name + '.swc')
        return sub_tree

    def write_subtree_by_slice_range(self, slice_begin, slice_end, path, neuronName, class_number):
        """
        This method gets slices from in a range specified by slice_begin and slice_end.
        :rtype: subtree
        """
        # if(slice_begin<slice_end and 0<=slice_begin and slice_end <= self._tree.max_height):
        if (slice_begin <= slice_end and 0 <= slice_begin):
            new_file_name = str(slice_begin) + "_" + str(slice_end)
            sub_tree_name = "./" + path + "/slice_" + new_file_name + "/"  # <-------
            # print(sub_tree_name)
            # sub_tree_name="./"+path+"/"
            if not os.path.exists(sub_tree_name):
                os.makedirs(sub_tree_name)

            sub_tree = self.get_tree_by_range(slice_begin, slice_end)
            if len(sub_tree.get_nodes()) > 1:
                sub_tree.other_write_SWC_tree_to_file(str(sub_tree_name) + str(class_number) + '__' + neuronName + '_' +
                                                      new_file_name + '.swc')
            else:
                sub_tree = None
            return sub_tree
        else:
            return None

    def re_write_tree(self, path, name_type, neuronName):
        sub_tree_name = "./" + path + "/" + name_type + "/"
        # sub_tree_name="./"+path+"/"
        if not os.path.exists(sub_tree_name):
            os.makedirs(sub_tree_name)
        sub_tree = self._tree
        sub_tree.re_write_SWC_tree_to_file(str(sub_tree_name) + neuronName + '_apical.swc')
        return sub_tree

    def set_level_to_tree(self):
        self._tree.set_level_to_node()

    def divide_neuron_into_apical_basal(self):
        self._tree.set_level_to_node()

    def search_limit(self, node, level):
        if node.level <= level:
            for child in node.children:
                self.search_limit(child, level)
        else:
            node.children = []

    def get_apical_dendrites(self, level):
        # node=self._tree.get_root()

        for node in self._tree:
            for child in node.children:
                if child.level > level:  # and node != self._tree.root:
                    node.children = []
                    # print str(node.content['p3d'])+" "+ str(node.level
        # print "quedan ---------------------"
        # for node in self._tree:
        #        print str(node.content['p3d'])+" - "+str(node.level)
        '''
        print "------------"
        print "root "+str(node.content['p3d'])
        print "filhos " + str(len(node.children))
        if node.level <= level:
            for child in node.children:
                    self.search_limit(child,level)
                else:
                    node.children = []
        else:
            node = []
        '''
        return self._tree

    def get_tree_by_level(self, level):

        for node in self._tree:
            for child in node.children:
                if child.level > level:  # and node != self._tree.root:
                    node.children = []

        return self._tree

    def get_tree_by_range(self, from_level, to_level):
        new_top_level = []
        root_radius = self._tree.root.content["p3d"].radius
        [rx, ry, rz] = self._tree.root.content["p3d"].xyz
        for node in self._tree:
            if node != self._tree.root:
                if node.level == from_level and node.get_parent().level == from_level - 1:  # the initial level
                    # ---------OK-cone----------
                    # node.set_parent(self._tree.root)
                    # new_top_level.append(node)
                    # ---------OK-cone----------
                    # #---------OK-move-the-point-near-to-soma------ to improve visualization
                    [nx, ny, nz] = node.content["p3d"].xyz
                    x = nx - rx
                    y = ny - ry
                    z = nz - rz
                    m = math.sqrt(x ** 2 + y ** 2 + z ** 2)
                    node.content["p3d"].xyz = array(
                        [round((x / m) * root_radius, 2) + rx, round((y / m) * root_radius, 2) + ry,
                         round((z / m) * root_radius, 2) + rz])
                    node.set_parent(self._tree.root)
                    new_top_level.append(node)
                    # #---------OK-move-the-point-near-to-soma------



                elif node.level > to_level:  # after the final level are removed
                    node.get_parent().children = []
            else:
                if node.level == from_level:  # case from_level 0 consider all current children of root
                    new_top_level = self._tree.root.children

        if (from_level == 0):  # case from_level==0 do not add anything just consider the current children on node
            soma_children = []
        else:
            soma_children = self.get_soma_children()
        new_children = soma_children + new_top_level  # join old soma children and new soma children
        self._tree.root.children = []
        self._tree.root.children = new_children

        return self._tree

    def get_tree_by_range_new(self, from_level, to_level):
        new_top_level = []
        root_radius = self._tree.root.content["p3d"].radius
        for node in self._tree:
            if node != self._tree.root:
                if node.level == from_level and node.get_parent().level == from_level - 1:  # the initial level
                    [x, y, z] = node.content["p3d"].xyz
                    tmp_node = btmorph.SNode2(node.index)
                    tmp_node = copy.copy(node)
                    m = math.sqrt(x ** 2 + y ** 2 + z ** 2)
                    tmp_node.content["p3d"].xyz = [round((x / m) * root_radius, 2), round((y / m) * root_radius, 2),
                                                   round((z / m) * root_radius, 2)]
                    tmp_node.set_parent(self._tree.root)
                    new_top_level.append(tmp_node)
                    node.set_parent(tmp_node)

                elif node.level > to_level:  # after the final level are removed
                    node.get_parent().children = []
            else:
                if node.level == from_level:
                    new_top_level = self._tree.root.children

        soma_children = self.get_soma_children()
        new_children = soma_children + new_top_level  # join old soma children and new soma children
        self._tree.root.children = []
        self._tree.root.children = new_children

        return self._tree

    def get_soma_children(self):
        children = self._tree.root.children
        soma_children = []
        for child in children:
            swc_type, t_node, parent_index = self._tree.total_all_nodes[child.index]
            if (swc_type == 1):  # if the children is soma copy
                soma_children.append(child)
        return soma_children

    def get_main_branches(self):

        trees = []

        for node in self._tree.root.children:
            if len(node.children) > 0:
                trees.append(self._tree.get_sub_tree(node))

        file_num = 1
        for sub_tree in trees:
            print (sub_tree.get_root())
            a = sub_tree.path_to_bifurcation(sub_tree.get_root())
            print (len(a))
            n = sub_tree.write_SWC_tree_to_file('sub_trees/sub_tree_' + str(file_num) + '.swc')
            file_num += 1

        return trees

    def get_points_of_interest(self):
        """
        Get lists containting the "points of interest", i.e., soma points, \
        bifurcation points and end/terminal points.

        Returns
        -------
        soma_points : list
        bif_points : list
        end_points : list

        """
        soma_points = []
        bif_points = []
        end_points = []

        # upated 2014-01-21 for compatibility with new btstructs2
        for node in self._all_nodes:
            if len(node.children) > 1:
                if not node.parent is None:
                    bif_points.append(node)  # the root is not a bifurcation
            if len(node.children) == 0:
                if node.parent.index != 1:  # "3 point soma", avoid the two side branches
                    end_points.append(node)
            if node.parent is None:
                soma_points = node

        return soma_points, bif_points, end_points

    """
    Global measures (1 for each tree)
    """

    def approx_soma(self):
        """
        *Scalar, global morphometric*

        By NeuroMorpho.org convention: soma surface ~ 4*pi*r^2, \
        where r is the abs(y_value) of point 2 and 3 in the SWC file


        Returns
        -------
        surface : float
             soma surface in micron squared

        """

        r = self._tree.get_node_with_index(1).content['p3d'].radius
        return 4.0 * np.pi * r * r

    def no_bifurcations(self):
        """
        *Scalar, global morphometric*

        Count the number of bifurcations points in a complete moprhology

        Returns
        -------
        no_bifurcations : int
             number of bifurcation
        """
        return len(self._bif_points)

    def no_terminals(self):
        """
        *Scalar, global morphometric*

        Count the number of temrinal points in a complete moprhology

        Returns
        -------
        no_terminals : int
            number of terminals
        """
        return len(self._end_points)

    def no_stems(self):
        """
        *Scalar, global morphometric*

        Count the number of stems in a complete moprhology (except the three \
        point soma from the Neuromoprho.org standard)


        Returns
        -------
        no_stems : int
            number of stems

        """
        return len(self._tree.root.children) - 2

    def total_length(self):
        """
        *Scalar, global morphometric*

        Calculate the total length of a complete morphology


        Returns
        -------
        total_length : float
            total length in micron

        """
        L = 0
        # upated 2014-01-21 for compatibility with new btstructs2
        for node in self._all_nodes:
            n = node.content['p3d']
            if not node.index in (1, 2, 3):
                p = node.parent.content['p3d']
                d = np.sqrt(np.sum((n.xyz - p.xyz) ** 2))
                L += d

        return L

    def total_surface(self):
        """
        *Scalar, global morphometric*

        Total neurite surface (at least, surface of all neurites excluding
        the soma. In accordance to the NeuroMorpho / L-Measure standard)

        Returns
        -------
        total_surface : float
            total surface in micron squared

        """
        total_surf = 0
        all_surfs = []
        # upated 2014-01-21 for compatibility with new btstructs2
        for node in self._all_nodes:
            n = node.content['p3d']
            if not node.index in (1, 2, 3):
                p = node.parent.content['p3d']
                H = np.sqrt(np.sum((n.xyz - p.xyz) ** 2))
                surf = 2 * np.pi * n.radius * H
                all_surfs.append(surf)
                total_surf = total_surf + surf
        return total_surf, all_surfs

    def total_volume(self):
        """
        *Scalar, global morphometric*

        Total neurite volume (at least, surface of all neurites excluding
        the soma. In accordance to the NeuroMorpho / L-Measure standard)

        Returns
        -------
        total_volume : float
            total volume in micron cubed

        """
        total_vol = 0
        all_vols = []
        # upated 2014-01-21 for compatibility with new btstructs2
        for node in self._all_nodes:
            n = node.content['p3d']
            if not node.index in (1, 2, 3):
                p = node.parent.content['p3d']
                H = np.sqrt(np.sum((n.xyz - p.xyz) ** 2))
                vol = np.pi * n.radius * n.radius * H
                all_vols.append(vol)
                total_vol = total_vol + vol
        return total_vol, all_vols

    def total_dimension(self):
        """
        *Scalar, global morphometric* Overall dimension of the morphology

        Returns
        -------
        dx : float
            x-dimension
        dy : float
            y-dimension
        dz : float
            z-dimension

        """
        dx, dy, dz, values = self.total_dimensions_verbose()
        return dx, dy, dz

    def total_dimensions_verbose(self):
        """
        *Scalar, global morphometric*

        Overall dimension of the whole moprhology. (No translation of the \
        moprhology according to arbitrary axes.)


        Returns
        -------
        dx : float
            x-dimension
        dy : float
            y-dimension
        dz : float
            z-dimension
        data : list
            minX,maxX,minY,maxY,minZ,maxZ

        """
        minX = sys.maxint
        maxX = -1 * sys.maxint
        minY = sys.maxint
        maxY = -1 * sys.maxint
        minZ = sys.maxint
        maxZ = -1 * sys.maxint
        for node in self._all_nodes:
            n = node.content['p3d']
            nx = n.xyz[0]
            ny = n.xyz[1]
            nz = n.xyz[2]
            minX = nx if nx < minX else minX
            maxX = nx if nx > maxX else maxX

            minY = ny if ny < minY else minY
            maxY = ny if ny > maxY else maxY

            minZ = nz if nz < minZ else minZ
            maxZ = nz if nz > maxZ else maxZ
        dx = np.sqrt((maxX - minX) * (maxX - minX))
        dy = np.sqrt((maxY - minY) * (maxY - minY))
        dz = np.sqrt((maxZ - minZ) * (maxZ - minZ))
        return dx, dy, dz, [minX, maxX, minY, maxY, minZ, maxZ]

    def global_horton_strahler(self):
        """
        Calculate Horton-Strahler number at the root
        See :func:`local_horton_strahler`

        Parameters
        ---------

        Returns
        ---------
        Horton-Strahler number at the root
        """
        return self.local_horton_strahler(self._tree.root)

    """
    Local measures
    """

    def get_diameters(self):
        """
        *Vector, local morphometric*

        Get the diameters of all points in the morphology
        """
        diams = []
        for node in self._all_nodes:
            if not node.index in (1, 2, 3):
                diams.append(node.content['p3d'].radius * 2.0)
        return diams

    def get_segment_pathlength(self, to_node):
        """
        *Vector, local morphometric*.

        Length of the incoming segment. Between this node and the soma or \
        another branching point. A path is defined as a stretch between \
        the soma and a bifurcation point, between bifurcation points, \
        or in between of a bifurcation point and a terminal point

        Parameters
        ----------
        to_node : :class:`btmorph.btstructs2.SNode2`
           Node *to* which the measurement is taken

        Returns
        -------
        length : float
            length of the incoming path in micron

        """
        # upated 2014-01-21 for compatibility with new btstructs2
        L = 0
        if self._tree.is_leaf(to_node):
            path = self._tree.path_to_root(to_node)
            L = 0
        else:
            path = self._tree.path_to_root(to_node)[1:]
            p = to_node.parent.content['p3d']
            n = to_node.content['p3d']
            d = np.sqrt(np.sum((n.xyz - p.xyz) ** 2))
            L = L + d

        for node in path:
            # print 'going along the path'
            n = node.content['p3d']
            if len(node.children) >= 2:  # I arrive at either the soma or a branchpoint close to the soma
                return L
            else:
                p = node.parent.content['p3d']
                d = np.sqrt(np.sum((n.xyz - p.xyz) ** 2))
                L = L + d

    def get_pathlength_to_root(self, from_node):
        """
        Length of the path between from_node to the root.
        another branching point

        Parameters
        ----------
        from_node : :class:`btmorph.btstructs2.SNode2`

        Returns
        -------
        length : float
            length of the path between the soma and the provided node

        """
        L = 0
        if self._tree.is_leaf(from_node):
            path = self._tree.path_to_root(from_node)
            L = 0
        else:
            path = self._tree.path_to_root(from_node)[1:]
            p = from_node.parent.content['p3d']
            n = from_node.content['p3d']
            d = np.sqrt(np.sum((n.xyz - p.xyz) ** 2))
            L = L + d

        for node in path[:-1]:
            # print 'going along the path'
            n = node.content['p3d']
            p = node.parent.content['p3d']
            d = np.sqrt(np.sum((n.xyz - p.xyz) ** 2))  # np.sqrt(np.sum((n.xyz-p.xyz)**2))
            L = L + d
        return L

    def get_segment_Euclidean_length(self, to_node):
        """
        Euclidean length to the incoming segment. Between this node and the soma or \
        another branching point

        Parameters
        ----------
        from_node : :class:`btmorph.btstructs2.SNode2`

        Returns
        -------
        length : float
            Euclidean distance *to* provided node (from soma or first branch point with lower order)

        """
        L = 0
        if self._tree.is_leaf(to_node):
            path = self._tree.path_to_root(to_node)
        else:
            path = self._tree.path_to_root(to_node)[1:]

        n = to_node.content['p3d']
        for node in path:
            if len(node.children) >= 2:
                return L
            else:
                p = node.parent.content['p3d']
                d = np.sqrt(np.sum((n.xyz - p.xyz) ** 2))
                L = d

    def get_Euclidean_length_to_root(self, from_node):
        """
        euclidean length between the from_node and the root

        Parameters
        ----------
        from_node : :class:`btmorph.btstructs2.SNode2`

        Returns
        -------
        length : float
            length of the path between the soma and the provided node

        """
        n = from_node.content['p3d']
        p = self._tree.root.content['p3d']
        d = np.sqrt(np.sum((n.xyz - p.xyz) ** 2))
        return d

    def degree_of_node(self, node):
        """
        Degree of a node. (The number of leaf node in the subtree mounted at \
        the provided node)

        Parameters
        ----------
        node : :class:`btmorph.btstructs2.SNode2`

        Returns
        -------
        degree : float
            degree of the subtree rooted at node

        """
        return self._tree.degree_of_node(node)

    def order_of_node(self, node):
        """
        Order of a node. (Going centrifugally away from the soma, the order \
        increases with 1 each time a bifurcation point is passed)

        Parameters
        ----------
        node : :class:`btmorph.btstructs2.SNode2`

        Returns
        -------
        order : float
            order of the subtree rooted at node

        """
        return self._tree.order_of_node(node)

    def partition_asymmetry(self, node):
        """
        *Vector, local morphometric*

        Compute the partition asymmetry for a given node.

        Parameters
        ----------
        node : :class:`btmorph.btstructs2.SNode2`

        Returns
        -------
        partition_asymmetry : float
            partition asymmetry of the subtree rooted at node (according to vanpelt and schierwagen 199x)

        """
        if node.children is None or len(node.children) == 1:
            return None
        d1 = self._tree.degree_of_node(node.children[0])
        d2 = self._tree.degree_of_node(node.children[1])
        if (d1 == 1 and d2 == 1):
            return 0  # by definition
        else:
            return np.abs(d1 - d2) / (d1 + d2 - 2.0)

    def bifurcation_angle_vec(self, node, where='local'):
        """
        *Vector, local morphometric*

        Only to be computed at branch points (_bif_points). Computes the angle
        between the two daughter branches in the plane defined by the \
        parent and the two daughters.

        cos alpha = :math:`(a \dot b) / (|a||b|)`

        Parameters
        -----------
        node : :class:`btmorph.btstructs2.SNode2`
        where : string
            either "local" or "remote". "Local" uses the immediate daughter \
            points while "remote" uses the point just before the next bifurcation or terminal point.

        Returns
        -------
        angle : float
            Angle in degrees
        """
        child_node1, child_node2 = self._get_child_nodes(node, where=where)
        scaled_1 = child_node1.content['p3d'].xyz - node.content['p3d'].xyz
        scaled_2 = child_node2.content['p3d'].xyz - node.content['p3d'].xyz
        amp = lambda a: np.sqrt(np.sum((a) ** 2))
        return np.arccos(np.dot(scaled_1, scaled_2) / (amp(scaled_1) * amp(scaled_2))) / (2 * np.pi / 360)

    def bifurcation_sibling_ratio(self, node, where='local'):
        """
        *Vector, local morphometric*

        Ratio between the diameters of two siblings.

        Parameters
        ----------
        node : :class:`btmorph.btstructs2.SNode2`
        where : string
            Toggle 'local' or 'remote'

        Returns
        -------
        result : float
            Ratio between the diameter of two siblings

        """
        child1, child2 = self._get_child_nodes(node, where=where)
        radius1 = child1.content['p3d'].radius
        radius2 = child2.content['p3d'].radius
        if radius1 > radius2:
            return radius1 / radius2
        else:
            return radius2 / radius1

    def _get_child_nodes(self, node, where):
        if where == 'local':
            return node.children[0], node.children[1]
        else:
            grandchildren = []
            for child in node.children:
                t_child = self._find_remote_child(child)
                grandchildren.append(t_child)
        return grandchildren[0], grandchildren[1]

    def _find_remote_child(self, node):
        t_node = node
        while len(t_node.children) < 2:
            if len(t_node.children) == 0:
                # print t_node, '-> found a leaf'
                return t_node
            t_node = t_node.children[0]
        # print t_node,' -> found a bif'
        return t_node

    def bifurcation_ralls_power_fmin(self, node, where='local'):
        """
        *Vector, local morphometric*

        Approximation of Rall's ratio using scipy.optimize.fmin.
        The error function is :math:`F={D_{d1}}^n+{D_{d2}}^n-{D_p}^n`

        Parameters
        ----------
        node : :class:`btmorph.btstructs2.SNode2`
        where : string
            either "local" or "remote". "Local" uses the immediate daughter \
            points while "remote" uses the point just before the next bifurcation or terminal point.

        Returns
        -------
        rr : float
            Appriximation of Rall's ratio
        """
        p_diam = node.content['p3d'].radius * 2
        child1, child2 = self._get_child_nodes(node, where=where)
        d1_diam = child1.content['p3d'].radius * 2
        d2_diam = child2.content['p3d'].radius * 2
        # print 'pd=%f,d1=%f,d2=%f' % (p_diam,d1_diam,d2_diam)

        if d1_diam >= p_diam or d2_diam >= p_diam:
            return np.nan

        import scipy.optimize
        mismatch = lambda n: np.abs(np.power(d1_diam, n) + np.power(d2_diam, n) - np.power(p_diam, n))

        p_lower = 0.0
        p_upper = 5.0  # THE associated mismatch MUST BE NEGATIVE

        best_n = scipy.optimize.fmin(mismatch, (p_upper - p_lower) / 2.0, disp=False)
        if 0.0 < best_n < 5.0:
            return best_n
        else:
            return np.nan

    def bifurcation_rall_ratio_classic(self, node, where='local'):
        """
        *Vector, local morphometric*

        The ratio :math:`\\frac{ {d_1}^p + {d_2}^p  }{D^p}` computed with :math:`p=1.5`

        Parameters
        -----------
        node : :class:`btmorph.btstructs2.SNode2`
        where : string
            either 'local or 'remote'. 'Local' uses the immediate daughter \
            points while "remote" uses the point just before the next bifurcation or terminal point.

        Returns
        -------
        rr : float
            Approximation of Rall's ratio

        """
        p_diam = node.content['p3d'].radius * 2
        child1, child2 = self._get_child_nodes(node, where=where)
        d1_diam = child1.content['p3d'].radius * 2
        d2_diam = child2.content['p3d'].radius * 2

        return (np.power(d1_diam, 1.5) + np.power(d2_diam, 1.5)) / np.power(p_diam, 1.5)

    def bifurcation_ralls_power_brute(self, node, where='local', min_v=0, max_v=5, steps=1000):
        """
        *Vector, local morphometric*

        Approximation of Rall's ratio.
        :math:`D^p = {d_1}^p + {d_2}^p`, p is approximated by brute-force checking the \
        interval [0,5] in 1000 steps (by default, but the exact search \
        dimensions can be specified by keyworded arguments.

        Parameters
        -----------
        node : :class:`btmorph.btstructs2.SNode2`
        where : string
            either 'local or 'remote'. 'Local' uses the immediate daughter \
            points while "remote" uses the point just before the next bifurcation or terminal point.

        Returns
        -------
        rr : float
            Approximation of Rall's power, p

        """
        p_diam = node.content['p3d'].radius * 2
        child1, child2 = self._get_child_nodes(node, where=where)
        d1_diam = child1.content['p3d'].radius * 2
        d2_diam = child2.content['p3d'].radius * 2
        # print 'pd=%f,d1=%f,d2=%f' % (p_diam,d1_diam,d2_diam)

        if d1_diam >= p_diam or d2_diam >= p_diam:
            return None

        test_v = np.linspace(min_v, max_v, steps)
        min_mismatch = 100000000000.0
        best_n = -1
        for n in test_v:
            mismatch = (np.power(d1_diam, n) + np.power(d2_diam, n)) - np.power(p_diam, n)
            # print "n=%f -> mismatch: %f" % (n,mismatch)
            if np.abs(mismatch) < min_mismatch:
                best_n = n
                min_mismatch = np.abs(mismatch)
        return best_n

    def _get_ampl_angle(self, node):
        """
        Compute the angle of this node on the XY plane and against the origin
        """
        pos_angle = lambda x: x if x > 0 else 180 + (180 + x)
        a = np.rad2deg(np.arctan2(node.content['p3d'].y, node.content['p3d'].x))
        return pos_angle(a)

    def local_horton_strahler(self, node):
        """
        We assign Horton-Strahler number to all nodes of a tree, in bottom-up order, as follows:

        If the node is a leaf (has no children), its Strahler number is one.
        If the node has one child with Strahler number i, and all other children have Strahler numbers less than i, then the Strahler number of the node is i again.
        If the node has two or more children with Strahler number i, and no children with greater number, then the Strahler number of the node is i + 1.
        *If the node has only one child, the Strahler number of the node equals to the Strahler number of the child
        The Strahler number of a tree is the number of its root node.

        See wikipedia for more information: http://en.wikipedia.org/wiki/Strahler_number

        Parameters
        ---------
        node : :class:`btmorph.btstructs2.SNode2`
            Node of interest
        Returns
        ---------
        hs : int
            The Horton-Strahler number (Strahler number) of the node
        """
        # Empy tree
        if node is None:
            return -1
        # Leaf => HS=1
        if len(node.children) == 0:
            return 1
        # Not leaf
        childrenHS = map(self.local_horton_strahler, node.children)
        return max(childrenHS + [(min(childrenHS) + 1)])

    def fractal_dimension_box_counting_core(self, vg):
        """
        Calculates fractal dimension of the given voxel grid by this formula:
        D = lim e -> 0 of (log(Ne)/log(e))
        http://rsbweb.nih.gov/ij/plugins/fraclac/FLHelp/Glossary.htm#db
        """
        # Box counting
        bc = BoxCounter(vg)
        if vg.res[2] == 0 or vg.res[2] == 1:
            startDim = min(vg.res[:-1]) / 2
        else:
            startDim = min(vg.res) / 2
        bc.grid_coverage(startDim)
        szs = map(lambda x: 2 ** x / float(max(vg.res)), range(1, int(math.log(startDim, 2)) + 1))
        cover = bc.coverageVals[1:-1]
        slope, intercept = np.polyfit(np.log(szs), np.log(cover), 1)
        return -slope

    def lacunarity_box_counting_core(self, vg):
        """
        Calculate lacunarity based on standard fixed grid box counting method with coef. of variation
        See wikipedia for more information: http://en.wikipedia.org/wiki/Lacunarity#equation_1
        Note: here we ignore orientations (all boxes start from (0,0,0)) and box sizes are always power of two

        Parameters
        ----------
        vg : :class:`btmorph.btstructs2.VoxelGrid`
            Ready to use voxel grid

        Returns
        -------
        lacunarity : float


        """
        bc = BoxCounter(vg)
        if vg.res[2] == 0 or vg.res[2] == 1:
            startDim = min(vg.res[:-1]) / 2
        else:
            startDim = min(vg.res) / 2
        bc.grid_count(startDim)
        lambdas = []
        for el in bc.countVals[1:-1]:
            lambdas.append((np.std(el) / np.mean(el)) ** 2)
        lc = np.mean(lambdas)
        return lc

    def fractal_dimension_lacunarity(self, voxelSize):
        """
        Calculate both lacunarity and fractal dimension of a tree.
        Faster than calling fractal_dim_box_counting and lacunarity_standard separately

        Parameters
        ----------
        voxelSize : number
            Desired voxel size, affects resolution. Both measures use voxelization of the 3D tree for calculations

        Returns
        ----------
        (lacunarity, fractal_dimension)
        """
        # Best resolution
        dx, dy, dz = self.total_dimension()
        res = [int(2 ** round(math.log(dx / voxelSize, 2))), int(2 ** round(math.log(dy / voxelSize, 2))), 1]
        if dz > 0.0:
            res[2] = int(2 ** round(math.log(dz / voxelSize, 2)))
        dim = [dx, dy, dz]
        self.vg = VoxelGrid(dim, res, self._tree)
        return self.frac_dim_lac(self.vg)

    def frac_dim_lac(self, vg=None):
        """
        Compute both lacunarity and fractal dimension
        Calculates lacunarity based on standard fixed grid box counting method with coef. of variation
        See wikipedia for more information: http://en.wikipedia.org/wiki/Lacunarity#equation_1
        Note: here we ignore orientations (all boxes start from (0,0,0)) and box sizes are always power of two
        Calculates fractal dimension of the given voxel grid by this formula:
        D = lim e -> 0 of (log(Ne)/log(e))
        http://rsbweb.nih.gov/ij/plugins/fraclac/FLHelp/Glossary.htm#db

        Parameters
        ----------
        vg : :class:`btmorph.btstructs2.VoxelGrid`
            Ready to use voxel grid

        Returns
        --------
        lacunarity, fractal_dimension : tuple
        """
        bc = BoxCounter(vg)
        if vg.res[2] == 0 or vg.res[2] == 1:
            startDim = min(vg.res[:-1]) / 2
        else:
            startDim = min(vg.res) / 2
        bc.grid_count(startDim)
        lambdas = []
        for el in bc.countVals[1:-1]:
            lambdas.append((np.std(el) / np.mean(el)) ** 2)
        lc = np.mean(lambdas)
        for i in range(1, len(bc.countVals)):
            s = sum(1 for e in bc.countVals[i] if e)
            bc.coverageVals[i] = s
        szs = map(lambda x: 2 ** (3 * x), range(1, int(math.log(startDim, 2)) + 1))
        szs2 = map(lambda x: 2 ** x / float(max(vg.res)), range(1, int(math.log(startDim, 2)) + 1))
        cover = bc.coverageVals[1:-1]
        slope, intercept = np.polyfit(np.log(szs2), np.log(cover), 1)
        # ----
        lc_slope, interc_lac = np.polyfit(np.log(szs), np.log(lambdas), 1)
        return (lc, -slope)

    def pca(self, A):
        """ performs principal components analysis
         (PCA) on the n-by-p data matrix A
         Rows of A correspond to observations, columns to variables.

         Returns :
          coeff :
        is a p-by-p matrix, each column containing coefficients
        for one principal component.
          score :
        the principal component scores; that is, the representation
        of A in the principal component space. Rows of SCORE
        correspond to observations, columns to components.
          latent :
        a vector containing the eigenvalues
        of the covariance matrix of A.
        source: http://glowingpython.blogspot.jp/2011/07/principal-component-analysis-with-numpy.html
        """
        # computing eigenvalues and eigenvectors of covariance matrix
        M = (A - mean(A.T, axis=1)).T  # subtract the mean (along columns)
        [latent, coeff] = linalg.eig(cov(M))  # attention:not always sorted
        score = dot(coeff.T, M)  # projection of the data in the new space
        return coeff, score, latent
