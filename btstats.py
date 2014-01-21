import sys
import numpy as np
import scipy
import scipy.optimize
import matplotlib.pyplot as plt
import pylab as p, time

class BTStats :
    '''
    Compute morphometric features and statistics of a single morphology

    Assume the "3 point" soma of the curated NeuroMorpho format. (see website)
    
    B. Torben-Nielsen, 2014-01 @ OIST/CNS (from BTN legacy code)
    '''
    
    def __init__(self,tree) :
        self._tree = tree
        self._all_nodes = self._tree.get_nodes()
        
        # compute some of the most used stats + 
        self._soma_points, self._bif_points, self._end_points = \
          self.get_points_of_interest()
    
    def get_points_of_interest(self) :
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
        for node in self._all_nodes :
            if len(node.get_child_nodes()) > 1  :
                if node._parent_node != None :
                    bif_points.append(node) # the root is not a bifurcation
            if len(node.get_child_nodes()) == 0  :
                if node.get_parent_node()._index != 1  :
                    end_points.append(node)
            if  node._parent_node == None :
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
        
        r = abs(self._tree.get_node_with_index(2).get_content()['p3d'].xyz[1]) # abs(y_value)
        print 'r=',r
        return 4.0*np.pi*r*r

    def no_bifurcations(self) :
        """
        *Scalar, global morphometric*
        
        Count the number of bifurcations points in a complete moprhology
        
        Returns
        -------
        no_bifurcations : int
             number of bifurcation
        """
        return len(self._bif_points)
        
    def no_terminals(self) :
        """
        *Scalar, global morphometric*
        
        Count the number of temrinal points in a complete moprhology
        
        Returns
        -------
        no_terminals : int
            number of terminals
        """        
        return len(self._end_points)
        
    def no_stems(self) :
        """
        *Scalar, global morphometric*
        
        Count the number of stems in a complete moprhology (except the three \
        point soma from the Neuromoprho.org standard)


        Returns
        -------
        no_stems : int
            number of stems
        
        """
        return len( self._tree.get_root().get_child_nodes() ) -2 
        
    def total_length(self) :
        """
        *Scalar, global morphometric*
        
        Calculate the total length of a complete morphology


        Returns
        -------
        total_length : float
            total length in micron
        
        """        
        L = 0
        # for node in self._all_nodes :
        #     n = node.get_content()['p3d']
        #     if(n.index not in [1,2,3]) :           
        #         p = node.get_parent_node().get_content()['p3d']
        #         d = np.sqrt( (n.x-p.x)*(n.x-p.x) + (n.y-p.y)*(n.y-p.y) + \
        #                      (n.z-p.z)*(n.z-p.z) )
        #         L += d
        
        # upated 2014-01-21 for compatibility with new btstructs2
        for node in self._all_nodes :
            n = node.get_content()['p3d']
            if(node._index not in [1,2,3]) :           
                p = node.get_parent_node().get_content()['p3d']
                d = np.sqrt(np.sum((n.xyz-p.xyz)**2))
                L += d
        
        return L

    def total_surface(self) :
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
        for node in self._all_nodes :
            n = node.get_content()['p3d']
            if node._index not in [1,2,3] :
                p = node.get_parent_node().get_content()['p3d']
                H = np.sqrt(np.sum((n.xyz-p.xyz)**2))
                surf = 2*np.pi*n.radius*H
                all_surfs.append(surf)
                total_surf = total_surf + surf
        return total_surf, all_surfs

    def total_volume(self) :
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
        for node in self._all_nodes :
            n = node.get_content()['p3d']
            if node._index not in [1,2,3] :
                p = node.get_parent_node().get_content()['p3d']
                H = np.sqrt(np.sum((n.xyz-p.xyz)**2))
                vol = np.pi*n.radius*n.radius*H
                all_vols.append(vol)
                total_vol = total_vol + vol
        return total_vol, all_vols

    def total_dimension(self) :
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
        dx,dy,dz,values = self.total_dimensions_verbose()
        return dx,dy,dz    


    def total_dimensions_verbose(self) :
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
        for node in self._all_nodes :
            n = node.get_content()['p3d']
            nx = n.xyz[0]
            ny = n.xyz[1]
            nz = n.xyz[2]            
            minX = nx if nx < minX else minX
            maxX = nx if nx > maxX else maxX

            minY = ny if ny < minY else minY
            maxY = ny if ny > maxY else maxY

            minZ = nz if nz < minZ else minZ
            maxZ = nz if nz > maxZ else maxZ
        dx = np.sqrt((maxX-minX)*(maxX-minX))
        dy = np.sqrt((maxY-minY)*(maxY-minY))
        dz = np.sqrt((maxZ-minZ)*(maxZ-minZ))
        return dx,dy,dz, [minX,maxX,minY,maxY,minZ,maxZ]

    """
    Local measures
    """
    def get_segment_pathlength(self,to_node) :
        """
        *Vector, local morphometric*. 

        Length of the incoming segment. Between this node and the soma or \
        another branching point. A path is defined as a stretch between \
        the soma and a bifurcation point, between bifurcation points, \
        or in between of a bifurcation point and a terminal point
        
        Parameters
        ----------
        to_node : SNode
           Node *to* which the measurement is taken

        Returns
        -------
        length : float
            length of the incoming path in micron
        
        """
        # upated 2014-01-21 for compatibility with new btstructs2
        L = 0
        if self._tree.is_leaf(to_node) :
            path = self._tree.path_to_root(to_node)
            L = 0
        else :
            path = self._tree.path_to_root(to_node)[1:]
            p = to_node.get_parent_node().get_content()['p3d']
            n = to_node.get_content()['p3d']
            # d = np.sqrt( (nx-px)*(nx-px) + (ny-py)*(ny-py) + (nz-pz)*(nz-pz) )
            d = np.sqrt(np.sum((n.xyz-p.xyz)**2))
            L = L + d
        
        for node in path :
            # print 'going along the path'
            n = node.get_content()['p3d']
            if len(node.get_child_nodes()) >= 2 :
                return L
            else :
                p = node.get_parent_node().get_content()['p3d']
                # d = np.sqrt( (nx-px)*(nx-px) + (ny-py)*(ny-py) + (nz-pz)*(nz-pz) )
                d = np.sqrt(np.sum((n.xyz-p.xyz)**2))
                L = L + d

    def get_pathlength_to_root(self,from_node) :
        """
        Length of the path between from_node to the root. 
        another branching point

        Parameters
        ----------
        from_node : SNode
        
        Returns
        -------
        length : float
            length of the path between the soma and the provided node
        
        """
        L = 0
        if self._tree.is_leaf(from_node) :
            path = self._tree.path_to_root(from_node)
            L = 0
        else :
            path = self._tree.path_to_root(from_node)[1:]
            p = from_node.get_parent_node().get_content()['p3d']
            n = from_node.get_content()['p3d']
            # px = p.xyz[0]
            # py = p.xyz[1]
            # pz = p.xyz[2]            
            # nx = n.xyz[0]
            # ny = n.xyz[1]
            # nz = n.xyz[2]            
            # d = np.sqrt( (nx-px)*(nx-px) + (ny-py)*(ny-py) + (nz-pz)*(nz-pz) ) 
            d = np.sqrt(np.sum((n.xyz-p.xyz)**2))
            L = L + d
        
        for node in path[:-1] :
            # print 'going along the path'
            n = node.get_content()['p3d']
            p = node.get_parent_node().get_content()['p3d']
            # d = np.sqrt( (nx-px)*(nx-px) + (ny-py)*(ny-py) + (nz-pz)*(nz-pz) )
            d = np.sqrt(np.sum((n.xyz-p.xyz)**2))
            L = L + d
        return L
                
    def get_segment_Euclidean_length(self,to_node) :
        """
        Euclidean length to the incoming segment. Between this node and the soma or \
        another branching point

        Parameters
        ----------
        from_node : SNode
        
        Returns
        -------
        length : float
            Euclidean distance *to* provided node (from soma or first branch point with lower order)
        
        """
        L = 0
        if self._tree.is_leaf(to_node) :
            path = self._tree.path_to_root(to_node)
        else :
            path = self._tree.path_to_root(to_node)[1:]

        n = to_node.get_content()['p3d']
        # nx = n.xyz[0]
        # ny = n.xyz[1]
        # nz = n.xyz[2]
        for node in path :
            if len(node.get_child_nodes()) >= 2 :
                return L
            else :
                p = node.get_parent_node().get_content()['p3d']
                # px = p.xyz[0]
                # py = p.xyz[1]
                # pz = p.xyz[2]            
                # d = np.sqrt( (nx-px)*(nx-px) + (ny-py)*(ny-py) + (nz-pz)*(nz-pz) )
                d = np.sqrt(np.sum((n.xyz-p.xyz)**2))
                L = d                

    def get_Euclidean_length_to_root(self,from_node) :
        """
        euclidean length between the from_node and the root

        Parameters
        ----------
        from_node : SNode
        
        Returns
        -------
        length : float
            length of the path between the soma and the provided node

        """
        n = from_node.get_content()['p3d']
        # nx = n.xyz[0]
        # ny = n.xyz[1]
        # nz = n.xyz[2]        
        p = self._tree.get_root().get_content()['p3d']
        # px = p.xyz[0]
        # py = p.xyz[1]
        # pz = p.xyz[2]                    
        # d = np.sqrt( (n.x-p.x)*(n.x-p.x) + (n.y-p.y)*(n.y-p.y) + (n.z-p.z)*(n.z-p.z) )
        d = np.sqrt(np.sum((n.xyz-p.xyz)**2))
        return d

    def degree_of_node(self,node) :
        """
        Degree of a node. (The number of leaf node in the subtree mounted at \
        the provided node)

        Parameters
        ----------
        node : SNode

        Returns
        -------
        degree : float
            degree of the subtree rooted at node
                
        """
        return self._tree.degree_of_node(node) 
        
    def order_of_node(self,node):
        """
        Order of a node. (Going centrifugally away from the soma, the order \
        increases with 1 each time a bifurcation point is passed)

        Parameters
        ----------
        node : SNode

        Returns
        -------
        order : float
            order of the subtree rooted at node        
        
        """        
        return self._tree.order_of_node(node)

    def partition_asymmetry(self,node) :
        """
        Compute the partition asymmetry for a given node. 

        Parameters
        ----------
        node : SNode

        Returns
        -------
        partition_asymmetry : float
            partition asymmetry of the subtree rooted at node (according to vanpelt and schierwagen 199x)
            
        """
        children =  node.get_child_nodes()
        if children == None or len(children) == 1 :
            return None 
        d1 = self._tree.degree_of_node(children[0])
        d2 = self._tree.degree_of_node(children[1])
        if(d1 == 1 and d2 == 1) :
            return 0 # by definition
        else :
            return np.abs(d1-d2)/(d1+d2-2.0)

    def bifurcation_angle(self,node,where='local') :
        """
        Compute the angles related to the parent - daughter constellation:
        amplitude, tilt and torque. 
        Can be computed locally, i.e., taking the parent segment and the first
        segment of the daughters, or, remote between the parent segment and
        next bifurcation or terminal point
        Follows: http://cng.gmu.edu:8080/Lm/help/index.htm
        
        Parameters
        ----------
        node : SNode
        where : string
            Toggle 'local' or 'remote' (see L-Measure website)

        Returns
        -------
        data : list of float
            amplitude, -1, -1 WHY?
        """
        child_node1,child_node2 = self._get_child_nodes(node,where=where)
        ampl_alpha = self._get_ampl_angle(child_node1)
        ampl_beta = self._get_ampl_angle(child_node2)
        if ampl_alpha > 360 or ampl_beta > 360 :
            print 'alpha=%f, beta=%f' % (ampl_alpha,ampl_beta)
            raw_input('ENTER')
        # ampl_gamma= ampl_alpha - ampl_beta if ampl_beta < ampl_alpha else \
        #  ampl_beta - ampl_alpha
        ampl_gamma = np.sqrt( (ampl_alpha-ampl_beta)**2  )  
        ampl_gamma = ampl_gamma if ampl_gamma <=180 else 360 -ampl_gamma
        # print "alpha=%f, beta=%f, gamma=%f" % (ampl_alpha,ampl_beta,ampl_gamma)

        return [ampl_gamma,-1,-1]

    def bifurcation_sibling_ratio(self,node,where='local') :
        """
        *Vector, local morphometric*

        Ratio between the diameters of two siblings. 

        Parameters
        ----------
        node : SNode
        where : string
            Toggle 'local' or 'remote'

        Returns
        -------
        result : float
            Ratio between the diameter of two siblings
        
        """
        child1,child2 = self._get_child_nodes(node,where=where)
        #print 'child1=',child1,', child2=',child2, ' @ ',where
        radius1 = child1.get_content()['p3d'].radius
        radius2 = child2.get_content()['p3d'].radius
        if radius1 > radius2 :
            return radius1 / radius2
        else :
            return radius2 / radius1
            
    def _get_child_nodes(self,node,where) :
        children = node.get_child_nodes()
        if where == 'local' : 
            return children[0],children[1]
        else :
            grandchildren = []
            for child in children :
                t_child = self._find_remote_child(child)
                grandchildren.append(t_child)
        return grandchildren[0],grandchildren[1]

    def _find_remote_child(self,node) :
        children = node.get_child_nodes()
        t_node = node
        while len(children) < 2 :
            if len(children) == 0 :
                # print t_node, '-> found a leaf'
                return t_node
            t_node = children[0]
            children = t_node.get_child_nodes()
        # print t_node,' -> found a bif'
        return t_node
                
    def bifurcation_ralls_ratio(self,node,precision=0.2,where='local') :
        """
        *Vector, local morphometric*

        Approximation of Rall's ratio. First implementation
        D^p = d1^p + d2^p, p being the approximated value of Rall's ratio
        """
        p_diam = node.get_content()['p3d'].radius*2
        child1,child2 = self._get_child_nodes(node,where=where)
        d1_diam = child1.get_content()['p3d'].radius*2
        d2_diam = child2.get_content()['p3d'].radius*2
        print 'pd=%f,d1=%f,d2=%f' % (p_diam,d1_diam,d2_diam)

        p_lower = 1.0
        p_upper = 5.0 # THE associated mismatch MUST BE NEGATIVE

        mismatch=100000000
        count_outer = 0
        while mismatch > precision and count_outer < 20:
            lower_mismatch = (np.power(d1_diam,p_lower) + np.power(d2_diam,p_lower))-np.power(p_diam,p_lower)
            upper_mismatch = (np.power(d1_diam,p_upper) + np.power(d2_diam,p_upper))-np.power(p_diam,p_upper)

            p_mid = (p_lower + p_upper)/2.0
            mid_mismatch = (np.power(d1_diam,p_mid) + np.power(d2_diam,p_mid))-np.power(p_diam,p_mid)

            # print 'p_lower=%f, p_mid=%f' % (p_lower,p_mid)   
            #print 'looking for the MID'
            count_inner = 0
            while mid_mismatch > 0 and count_inner < 20:
                p_lower = p_mid
                lower_mismatch = (np.power(d1_diam,p_lower) + np.power(d2_diam,p_lower))-np.power(p_diam,p_lower)
                p_mid = (p_lower + p_upper)/2.0
                mid_mismatch = (np.power(d1_diam,p_mid) + np.power(d2_diam,p_mid))-np.power(p_diam,p_mid)
                # print 'p_lower=%f, p_mid=%f, p_upper=%f' % (p_lower,p_mid,p_upper)
                # print 'lower=%f, mid=%f, upper=%f' % (lower_mismatch,mid_mismatch,upper_mismatch)
                count_inner = count_inner + 1
            #print 'found the MID'
            # print '\nlow_m=%f,mid_m=%f,up_m=%f' % (lower_mismatch,mid_mismatch,upper_mismatch)
            # print 'low_p=%f,mid_p=%f,up_p=%f' % (p_lower,p_mid,p_upper)

            if upper_mismatch < mid_mismatch :
                p_upper = p_mid
            else :
                p_lower = p_mid
            mismatch = np.abs(mid_mismatch)
            count_outer = count_outer + 1
            #print 'mismatch=%f, p_lower=%f,p_upper=%f' % (mismatch,p_lower,p_upper)
            # raw_input('Enter')
        if count_outer >=19 :
            return np.nan
        else :
            return p_mid

    def bifurcation_ralls_ratio2(self,node,precision=0.2,where='local') :
        """
        *Vector, local morphometric*

        Approximation of Rall's ratio. Second implementation
        D^p = d1^p + d2^p, p being the approximated value of Rall's ratio
        
        """
        p_diam = node.get_content()['p3d'].radius*2
        child1,child2 = self._get_child_nodes(node,where=where)
        d1_diam = child1.get_content()['p3d'].radius*2
        d2_diam = child2.get_content()['p3d'].radius*2
        print 'pd=%f,d1=%f,d2=%f' % (p_diam,d1_diam,d2_diam)

        if d1_diam >= p_diam or d2_diam >= p_diam :
            return 1

        p_lower = 1.0
        p_upper = 5.0 # THE associated mismatch MUST BE NEGATIVE

        mismatch=100000000
        count_outer = 0
        while mismatch > precision and count_outer < 100:
            lower_mismatch = (np.power(d1_diam,p_lower) + np.power(d2_diam,p_lower))-np.power(p_diam,p_lower)
            upper_mismatch = (np.power(d1_diam,p_upper) + np.power(d2_diam,p_upper))-np.power(p_diam,p_upper)

            p_mid = (p_lower + p_upper)/2.0
            mid_mismatch = (np.power(d1_diam,p_mid) + np.power(d2_diam,p_mid))-np.power(p_diam,p_mid)

            # print 'p_lower=%f, p_mid=%f' % (p_lower,p_mid)   
            #print 'looking for the MID'

            if upper_mismatch < mid_mismatch :
                p_upper = p_mid
            else :
                p_lower = p_mid
            mismatch = np.abs(mid_mismatch)
            count_outer = count_outer + 1
            #print 'mismatch=%f, p_lower=%f,p_upper=%f' % (mismatch,p_lower,p_upper)
            # raw_input('Enter')
        if count_outer >=19 :
            return np.nan
        else :
            return p_mid
        

    def _get_ampl_angle(self,node) :
        """
        Compute the angle of this node on the XY plane and against the origin
        """
        pos_angle = lambda x: x if x > 0 else 180 + (180+x)
        a = np.rad2deg(np.arctan2(node.get_content()['p3d'].y,node.get_content()['p3d'].x))
        return pos_angle(a)
