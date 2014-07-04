from btmorph import P3D2, SNode2

class STree2 :
    '''
    Simple tree for use with a simple Node (:class:`SNode2`).

    While the class is designed to contain binary trees (for neuronal morphologies) the number of children is not limited.
    As such, this is a generic implementation of a tree structure as a linked list.
    '''
    
    def __init__(self) :
        self.root = property(self.get_root, self.set_root, None, None)
        self.root = None
        
    def set_root(self,node) :
        """
        Set the root node of the tree

        Parameters
        -----------
        node : :class:`SNode2`
            to-be-root node
        """
        self._root = node
        if not self.root is None:
            self.root.parent = None
        
    def get_root(self) :
        """
        Obtain the root node

        Returns
        -------
        root : :class:`SNode2`
        """
        return self._root
        
    def is_root(self,node) :
        """
        Check whether a node is the root node

        Returns
        --------
        is_root : boolean
            True is the queried node is the root, False otherwise
        """
        if node.parent is None:
            return True
        else :
            return False
            
    def is_leaf(self,node) :
        """
        Check whether a node is a leaf node, i.e., a node without children

        Returns
        --------
        is_leaf : boolean
            True is the queried node is a leaf, False otherwise
        """
        if len(node.children) == 0  :
            return True
        else :
            return False
            
    def add_node_with_parent(self,node,parent) :
        """
        Add a node to the tree under a specific parent node

        Parameters
        -----------
        node : :class:`SNode2`
            node to be added
        parent : :class:`SNode2`
            parent node of the newly added node
        """
        node.parent = parent
        if not parent is None:
            parent.add_child(node)
        
    def remove_node(self,node) :
        """
        Remove a node from the tree

        Parameters
        -----------
        node : :class:`SNode2`
            node to be removed
        """
        node.parent.remove_child(node)
        self._deep_remove(node)
        
            
    def _deep_remove(self,node) :
        children = node.children
        node.make_empty()
        for child in children :
            self._deep_remove(child)        

    def get_nodes(self) :
        """
        Obtain a list of all nodes int the tree

        Returns
        -------
        all_nodes : list of :class:`SNode2`
        """
        n = []
        self._gather_nodes(self.root,n) 
        return n

    def get_sub_tree(self,fake_root) :
        """
        Obtain the subtree starting from the given node

        Parameters
        -----------
        fake_root : :class:`SNode2`
            Node which becomes the new root of the subtree

        Returns
        -------
        sub_tree :  STree2
            New tree with the node from the first argument as root node
        """
        ret = STree2()
        cp = fake_root.__copy__()
        cp.parent = None
        ret.root = cp
        return ret

    def _gather_nodes(self,node,node_list) :
        if not node is None:
            node_list.append(node)
            for child in node.children :
                self._gather_nodes(child,node_list)
    
    def get_node_with_index(self, index) :
        """
        Get a node with a specific name. The name is always an integer

        Parameters
        ----------
        index : int
            Name of the node to be found

        Returns
        -------
        node : :class:`SNode2`
            Node with the specific index
        """
        return self._find_node(self.root,index)
        
    def get_node_in_subtree(self,index,fake_root) :
        """
        Get a node with a specific name in a the subtree rooted at fake_root. The name is always an integer

        Parameters
        ----------
        index : int
            Name of the node to be found
        fake_root: :class:`SNode2`
            Root node of the subtree in which the node with a given index is searched for

        Returns
        -------
        node : :class:`SNode2`
            Node with the specific index
        """
        return self._find_node(fake_root,index)
        
    def _find_node(self,node,index) :
        """
        Sweet breadth-first/stack iteration to replace the recursive call. 
        Traverses the tree until it finds the node you are looking for.

        Parameters
        -----------

        
        Returns
        -------
        node : :class:`SNode2`
            when found and None when not found
        """
        stack = []; 
        stack.append(node)
        while(len(stack) != 0) :
            for child in stack :
                if child.index == index  :
                    return child
                else :
                    stack.remove(child)
                    for cchild in child.children :
                        stack.append(cchild)
        return None # Not found!
        
    def degree_of_node(self,node) :
        """
        Get the degree of a given node. The degree is defined as the number of leaf nodes in the subtree rooted at this node.

        Parameters
        ----------
        node : :class:`SNode2`
            Node of which the degree is to be computed.

        Returns
        -------
        degree : int
        """
        sub_tree = self.get_sub_tree(node)
        st_nodes = sub_tree.get_nodes()
        leafs = 0
        for n in st_nodes :
            if sub_tree.is_leaf(n) :
                leafs = leafs +1
        return leafs
        
    def order_of_node(self,node) :
        """
        Get the order of a given node. The order or centrifugal order is defined as 0 for the root and increased with any bifurcation.
        Hence, a node with 2 branch points on the shortest path between that node and the root has order 2.

        Parameters
        ----------
        node : :class:`SNode2`
            Node of which the order is to be computed.

        Returns
        -------
        order : int
        """
        ptr = self.path_to_root(node)
        order = 0
        for n in ptr :
            if len(n.children) > 1  :
                order = order+1
        # order is on [0,max_order] thus subtract 1 from this calculation
        return order-1 
                
    def path_to_root(self,node) :
        """
        Find and return the path between a node and the root.

        Parameters
        ----------
        node : :class:`SNode2`
            Node at which the path starts

        Returns
        -------
        path : list of :class:`SNode2`
            list of :class:`SNode2` with the provided node and the root as first and last entry, respectively.
        """
        n = []
        self._go_up_from(node,n)            
        return n
        
    def _go_up_from(self,node,n):
        n.append(node)
        if not node.parent is None :
            self._go_up_from(node.parent,n)

    def path_between_nodes(self,from_node,to_node) :
        """
        Find the path between two nodes. The from_node needs to be of higher \
        order than the to_node. In case there is no path between the nodes, \
        the path from the from_node to the soma is given.

        Parameters
        -----------
        from_node : :class:`SNode2`
        to_node : :class:`SNode2`
        """
        n = []
        self._go_up_from_until(from_node,to_node,n)
        return n

    def _go_up_from_until(self,from_node,to_node,n) :
        n.append(from_node)
        if from_node == to_node :
            return
        if not node.parent is None :
            self._go_up_from_until(node.parent,to_node,n)

    def write_SWC_tree_to_file(self,file_n) :
        """
        Non-specific for a tree. Used to write an SWC file from a morphology stored in this tree.

        Parameters
        ----------
        file_n : str
            name of the file to open
        """
        writer = open(file_n,'w')
        nodes = self.get_nodes()
        nodes.sort()

        # 3 point soma representation (See Neuromoprho.org FAQ)
        s1p = nodes[0].content["p3d"]
        s1_xyz = s1p.xyz
        s2p = nodes[1].content["p3d"]
        s2_xyz = s2p.xyz
        s3p = nodes[2].content["p3d"]
        s3_xyz = s3p.xyz
        soma_str = "1 1 " +str(s1_xyz[0]) + " " + str(s1_xyz[1]) + \
          " " + str(s1_xyz[2]) + " " + str(s1p.radius) + " -1\n" + \
          "2 1 " +str(s2_xyz[0]) + " " + str(s2_xyz[1]) + \
          " " + str(s2_xyz[2]) + " " + str(s2p.radius) + " 1\n" + \
          "3 1 " +str(s3_xyz[0]) + " " + str(s3_xyz[1]) + \
          " " + str(s3_xyz[2]) + " " + str(s3p.radius) + " 1\n"
        writer.write(soma_str)
        writer.flush()
        
        # add the soma compartment, then enter the loop
        for node in nodes[3:] :
            p3d = node.content['p3d'] # update 2013-03-08
            xyz = p3d.xyz
            radius = p3d.radius
            tt = p3d.type
            p3d_string = str(node.index)+' '+str(tt) + ' ' + \
                str(xyz[0]) + ' ' + str(xyz[1])+ ' ' + str(xyz[2]) + \
                ' ' + str(radius) + ' ' \
                + str(node.parent.index)
            # print 'p3d_string: ', p3d_string
            writer.write( p3d_string + '\n' )
            writer.flush()
        writer.close()          
        #print 'STree::writeSWCTreeToFile -> finished. Tree in >',fileN,'<'
        
    def read_SWC_tree_from_file(self,file_n,types=range(1,10)) :
        """
        Non-specific for a tree.
        Read and load a morphology from an SWC file. and parse it into an STree2 object. 
        On the NeuroMorpho.org website, 5 types of somadescriptions are 
        considered. (http://neuromorpho.org/neuroMorpho/SomaFormat.html)
        
        Parameters
        -----------
        file_n : str
        name of the file to open
        soma_type : int
            [1-5], see specs
        """
        import numpy as np
        file = open(file_n,'r')
        all_nodes = dict()
        for line in file :
            if not line.startswith('#') :
                split = line.split()
                index = int(split[0].rstrip())
                type = int(split[1].rstrip())
                x = float(split[2].rstrip())
                y = float(split[3].rstrip())
                z = float(split[4].rstrip())
                radius = float(split[5].rstrip())
                parent_index = int(split[6].rstrip())

                if type in types:
                    tP3D = P3D2(np.array([x,y,z]),radius,type)
                    t_node = SNode2(index)
                    t_node.content = {'p3d': tP3D}
                    all_nodes[index] = (t_node,parent_index)
                
        for index,(node,parent_index) in all_nodes.items() :
            if index == 1:
                print('Setting root node: %s.' % str(node))
                self.root = node
            elif index in (2,3):
                # the 3-point soma representation (http://neuromorpho.org/neuroMorpho/SomaFormat.html)
                self.add_node_with_parent(node,self.root)
            else:
                parent_node = all_nodes[parent_index][0]
                self.add_node_with_parent(node,parent_node)

        return self

    def __str__(self) :
        return "STree2 ("+str(len(self.get_nodes()))+" nodes)"

