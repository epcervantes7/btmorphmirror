"""
Module contains:
    P3D2
    SNode2   
    STree2

B. Torben-Nielsen, 2013-10 @ OIST,
2013-01 @ BBP
from BTN legacy code
"""


class P3D2 :
    def __init__(self,xyz,radius,type=7) :
        self.xyz = xyz
        self.radius = radius
        self.type = type

    def __str__(self) :
        return "P3DD " + self.xyz + ", R=",self.radius, ", T=",self.type

class SNode2 :
    """
    Simple Node for use with a simple Tree (STree)
    
    By design, the "_content" should be a dictionary. (2013-03-08)
    """
    
    def __init__(self,index) :
        self._parent_node = None
        self._index = index
        self._child_nodes = []
        
    def get_index(self) :
        """
        Return the index of this node

        Returns
        -------
        index : int
        """
        return self._index
        
    def get_parent_node(self) :
        return self._parent_node
        
    def get_child_nodes(self) :
        return self._child_nodes
        
    def get_content(self) :
        return self._content
    
    def set_index(index) :
        self._index = index
    
    def set_parent_node(self,parent_node) :
        self._parent_node = parent_node
        
    def set_content(self,content) :
        if isinstance(content,dict) :
            self._content = content 
        else :
            raise Exception("SNode2.set_content must receive a dict")    
        
    def add_child(self,child_node) :
        self._child_nodes.append(child_node)
            
    def make_empty(self):
        self._parent_node = None
        self._content = None
        self._child_nodes = []
            
    def remove_child(self, child_node) :
        self._child_nodes.remove(child_node)
        
    def __str__(self) :
        return 'SNode2 (ID: '+str(self._index)+')'

    def __lt__(self,other):
        if self._index < other._index :
            return True
    def __le__(self,other):
        if self._index <= other._index :
            return True
    def __gt__(self,other):
        if self._index > other._index :
            return True
    def __ge__(self,other):
        if self._index >= other._index :
            return True
    
    def __copy__(self) : # customization of copy.copy
        ret = SNode2(self._index)
        for child in self._child_nodes :
            ret.add_child(child)
        ret.content = self._content
        ret.set_parent_node(self._parent_node)
        return ret
        
class STree2 :
    """
    Simple tree for use with a simple Node (SNode2)
    """
    
    def __init__(self) :
        """ Initialize an empty tree by default
        """
        _root = None
        
    def set_root(self,node) :
        self._root = node
        self._root.set_parent_node(None)
        
    def get_root(self) :
        return self._root
        
    def is_root(self,node) :
        if node.get_parent_node() != None :
            return False
        else :
            return True
            
    def is_leaf(self,node) :
        if len(node.get_child_nodes()) == 0  :
            return True
        else :
            return False
            
    def add_node_with_parent(self,node,parent) :
        node.set_parent_node(parent)
        parent.add_child(node)
        
    def remove_node(self,node) :
        node.get_parent_node().remove_child(node)
        self._deep_remove(node)
        
            
    def _deep_remove(self,node) :
        children = node.get_child_nodes()
        node.make_empty()
        for child in children :
            self._deep_remove(child)        

    def get_nodes(self) :
        n = []
        self._gather_nodes(self._root,n) 
        return n 

    # def get_indexed_nodes(self) :
    #     n = {}
    #     self._gather_nodes(self._root,n) 
    #     return n 
            
    def get_sub_tree(self,fake_root) :
        ret = STree2()
        cp = fake_root.__copy__()
        cp.set_parent_node(None)
        ret.set_root(cp)
        return ret

    def _gather_nodes(self,node,node_list) :
        node_list.append(node)
        for child in node.get_child_nodes() :
            self._gather_nodes(child,node_list)
    
    def get_node_with_index(self, index) :
        return self._find_node(self._root,index)
        
    def get_node_in_subtree(self,index,fake_root) :
        return self._find_node(fake_root,index)
        
    def _find_node(self,node,index) :
        """
        Sweet breadth-first/stack iteration to replace the recursive call. 
        Traverses the tree until it finds the node you are looking for.     
        Returns SNode2 when found and None when not found
        """
        stack = []; 
        stack.append(node)
        while(len(stack) != 0) :
            for child in stack :
                if child.get_index() == index  :
                    return child
                else :
                    stack.remove(child)
                    for cchild in child.get_child_nodes() :
                        stack.append(cchild)
        return None # Not found!
        
    def degree_of_node(self,node) :
        sub_tree = self.get_sub_tree(node)
        st_nodes = sub_tree.get_nodes()
        leafs = 0
        for n in st_nodes :
            if sub_tree.is_leaf(n) :
                leafs = leafs +1
        return leafs
        
    def order_of_node(self,node) :
        # if node == self._root :
        #     return 0
        ptr =self.path_to_root(node)
        order = 0
        for n in ptr :
            if len(n.get_child_nodes()) > 1  :
                order = order +1
        """ order is on [0,max_order] thus subtract 1 from this calculation """
        return order -1 
                
    def path_to_root(self,node) :
        n = []
        self._go_up_from(node,n)            
        return n
        
    def _go_up_from(self,node,n):
        n.append(node)
        p_node = node.get_parent_node()
        if p_node != None :
            self._go_up_from(p_node,n)

    def path_between_nodes(self,from_node,to_node) :
        """
        Find the path between two nodes. The from_node needs to be of higher \
        order than the to_node. In case there is no path between the nodes, \
        the path from the from_node to the soma is given.
        """
        n = []
        self._go_up_from_until(from_node,to_node,n)
        return n

    def _go_up_from_until(self,from_node,to_node,n) :
        n.append(from_node)
        if from_node == to_node :
            return
        p_node = from_node.get_parent_node()
        if p_node != None :
            self._go_up_from_until(p_node,to_node,n)

    def write_SWC_tree_to_file(self,file_n) :
        """
        Save a tree to an SWC file
        TODO
        """
        # raise Exception("Not yet implemented")
        writer = open(file_n,'w')
        nodes = self.get_nodes()
        nodes.sort()

        # 3 point soma representation (See Neuromoprho.org FAQ)
        s1p = nodes[0].get_content()["p3d"]
        s1_xyz = s1p.xyz
        s2p = nodes[1].get_content()["p3d"]
        s2_xyz = s2p.xyz
        s3p = nodes[2].get_content()["p3d"]
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
            p3d = node.get_content()['p3d'] # update 2013-03-08
            xyz = p3d.xyz
            radius = p3d.radius
            tt = p3d.type
            p3d_string = str(node.get_index())+' '+str(tt) + ' ' + \
                str(xyz[0]) + ' ' + str(xyz[1])+ ' ' + str(xyz[2]) + \
                ' ' + str(radius) + ' ' \
                + str(node.get_parent_node().get_index())
            # print 'p3d_string: ', p3d_string
            writer.write( p3d_string + '\n' )
            writer.flush()
        writer.close()          
        #print 'STree::writeSWCTreeToFile -> finished. Tree in >',fileN,'<'
        
    def read_SWC_tree_from_file(self,file_n) :
        """
        Read and load a morphology from an SWC file. 
        On the NeuroMorpho.org website, 5 types of somadescriptions are 
        considered. (http://neuromorpho.org/neuroMorpho/SomaFormat.html)
        Arguments:
        - file_n : file name
        - soma_type : [1-5], see specs
        """
        import numpy as np
        file = open(file_n,'r')
        all_nodes = dict()
        for line in file :
            if not line.startswith('#') :
                split= line.split()
                index = int(split[0].rstrip())
                type = int(split[1].rstrip())
                x = float(split[2].rstrip())
                y = float(split[3].rstrip())
                z = float(split[4].rstrip())
                radius = float(split[5].rstrip())
                parent_index = int(split[6].rstrip())
                
                tP3D = P3D2(np.array([x,y,z]),radius,type)
                t_node = SNode2(index)
                t_node.set_content({'p3d':tP3D})
                all_nodes[index] = (t_node,parent_index)
                
        for index, (node,parent_index) in all_nodes.items() :
            if index == 1:
                self.set_root(node)
            elif index == 2 or index==3:
                # the 3-point soma representation (http://neuromorpho.org/neuroMorpho/SomaFormat.html)
                self.add_node_with_parent(node,self._root)
            else:
                parent_node = all_nodes[parent_index][0]
                self.add_node_with_parent(node,parent_node)
            
        return self

    def __str__(self) :
        return "STree2 ("+str(len(self.get_nodes()))+" nodes)"
