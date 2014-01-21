"""
Basic visualization of neurite morphologies. Color coding for individual \
sections is supported. Also synapse can be drawn.

B. Torben-Nielsen @ OIST/CNS (updated BTN legacy code)
"""
import sys,time
sys.setrecursionlimit(10000)

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D

import btmorph

""" internal constants required for the dendrogram generation """
H_SPACE = 20
V_SPACE = 0
C = 'k'

max_width = 0
max_height = 0

def plot_2D_SWC(file_name='P20-DEV139.CNG.swc',cs=None,synapses=None,syn_cs=None,outN=None,offset=None,show_axis=False) :
    '''
    Colors can be
    None: uniform/default matplotlib color
    Any color code: uniform but specified color
    array of values:
    colormap?
    '''
    # read the SWC into a dictionary: key=index, value=(x,y,z,d,parent)
    x = open(file_name,'r')
    SWC = {}
    for line in x :
        if(not line.startswith('#')) :
            splits = line.split()
            index = int(splits[0])
            n_type = int(splits[1])
            if offset == None :
                x = float(splits[2])
                y = float(splits[3])
                z = float(splits[4])
            else :
                x = offset[0] + float(splits[2])
                y = offset[1] + float(splits[3])
                z = float(splits[4])
            r = float(splits[5])
            parent = int(splits[-1])
            SWC[index] = (x,y,z,r,parent,n_type)

    my_color_list = ['r','g','b','c','m','y','r--','b--','g--']
            
    if cs == None :
        pass#plt.clf()
    else :
        norm = colors.normalize(np.min(cs),np.max(cs))
        Z = [[0,0],[0,0]]
        levels=range(int(np.min(cs)),int(np.max(cs)),100)
        levels = np.linspace(np.min(cs),np.max(cs),100)
        CS3 = plt.contourf(Z,levels,cmap=cm.jet)
        plt.clf()

    min_y = 100000.0

    syn_index = -1
    for index in SWC.keys() : # not ordered but that has little importance here
        # draw a line segment from parent to current point
        current_SWC = SWC[index]
        #print 'index: ', index, ' -> ', current_SWC
        c_x = current_SWC[0]
        c_y = current_SWC[1]
        c_z = current_SWC[2]
        c_r = current_SWC[3]
        parent_index = current_SWC[4]

        if(c_y < min_y) :
            min_y = c_y
                
        if(index <= 3) :
            print 'do not draw the soma and its CNG, 2 point descriptions'
        else :
            parent_SWC = SWC[parent_index]
            p_x = parent_SWC[0]
            p_y = parent_SWC[1]
            p_z = parent_SWC[2]
            p_r = parent_SWC[3]
            if(p_y < min_y) :
                min_y= p_y
            # print 'index:', index, ', len(cs)=', len(cs)
            if(cs == None) :
                pl = plt.plot([p_x,c_x],[p_y,c_y],my_color_list[current_SWC[5]-1],linewidth=c_r)
                #pl = plt.plot([p_x,c_x],[p_y,c_y],'r',linewidth=c_r)
            else :
                try :
                    pl = plt.plot([p_x,c_x],[p_y,c_y],c=cm.jet(norm(cs[index])),linewidth=c_r)
                except Exception :
                    pass# it's ok: it's the list size...

        # add the synapses
        if(synapses != None) :
            if index in synapses :
                syn_index = syn_index + 1
                if syn_cs :
                    plt.plot(c_x,c_y,syn_cs)
                else :
                    plt.plot(c_x,c_y,'ro')


    frame1 = plt.gca()

    if show_axis :
        frame1.axes.get_xaxis().set_ticks([])
        frame1.axes.get_yaxis().set_ticks([])
        frame1.axes.get_xaxis().set_visible(False)
        frame1.axes.get_yaxis().set_visible(False)
        
    # plt.xlabel('X')
    # plt.ylabel('Y')

    # draw a scale bar
    if offset == None :
        plt.plot([0,100],[min_y*1.1,min_y*1.1],'k',linewidth=5) # 250 for MN, 100 for granule
    else :
        pass # when offset is used, the figure is to scale and no scalebar needed
    
    if(cs != None) :
        cb = plt.colorbar(CS3) # bit of a workaround, but it seems to work
        ticks_f = np.linspace(np.min(cs)-1,np.max(cs)+1,5)
        ticks_i = map(int,ticks_f)
        cb.set_ticks(ticks_i)

    if(outN != None) :
        plt.savefig(outN)


def plot_3D_SWC(file_name='P20-DEV139.CNG.swc',cs=None,synapses=None,syn_cs=None,outN=None) :
    """
    Matplotlib rendering of a SWC described morphology in 3D
    Colors can be
    None: uniform/default matplotlib color
    Any color code: uniform but specified color
    array of values:
    colormap?
    """
    my_color_list = ['r','g','b','c','m','y','r--','b--','g--']
    
    if(cs == None) :
        pass
    else :
        norm = colors.normalize(np.min(cs),np.max(cs))
        Z = [[0,0],[0,0]]
        levels=range(int(np.min(cs)),int(np.max(cs)),100)
        levels = np.linspace(np.min(cs),np.max(cs),100)
        CS3 = plt.contourf(Z,levels,cmap=cm.jet)
        plt.clf()
    
    # read the SWC into a dictionary: key=index, value=(x,y,z,d,parent)
    x = open(file_name,'r')
    SWC = {}
    for line in x :
        if(not line.startswith('#')) :
            splits = line.split()
            index = int(splits[0])
            n_type = int(splits[1])
            x = float(splits[2])
            y = float(splits[3])
            z = float(splits[4])
            r = float(splits[5])
            parent = int(splits[-1])
            SWC[index] = (x,y,z,r,parent,n_type)

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    for index in SWC.keys() : # not ordered but that has little importance here
        # draw a line segment from parent to current point
        current_SWC = SWC[index]
        #print 'index: ', index, ' -> ', current_SWC
        c_x = current_SWC[0]
        c_y = current_SWC[1]
        c_z = current_SWC[2]
        c_r = current_SWC[3]
        parent_index = current_SWC[4]
                
        if(index <= 3) :
            print 'do not draw the soma and its CNG, !!! 2 !!! point descriptions'
        else :
            parent_SWC = SWC[parent_index]
            p_x = parent_SWC[0]
            p_y = parent_SWC[1]
            p_z = parent_SWC[2]
            p_r = parent_SWC[3]
            # print 'index:', index, ', len(cs)=', len(cs)
            if cs == None :
                pl = plt.plot([p_x,c_x],[p_y,c_y],[p_z,c_z],my_color_list[current_SWC[5]-1],linewidth=c_r)
            else :
                try :
                    pl = plt.plot([p_x,c_x],[p_y,c_y], \
                                  c=cm.jet(norm(cs[index])),linewidth=c_r)
                except Exception :
                    print 'something going wrong here'
                    # pass# it's ok: it's the list size...

        # add the synapses
        if(synapses != None) :
            if(index in synapses) :
                if syn_cs :
                    plt.plot(c_x,c_y,syn_cs)
                else :
                    plt.plot(c_x,c_y,'ro')

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

    if(cs != None) :
        cb = plt.colorbar(CS3) # bit of a workaround, but it seems to work
        ticks_f = np.linspace(np.min(cs)-1,np.max(cs)+1,5)
        ticks_i = map(int,ticks_f)
        cb.set_ticks(ticks_i)

    if(outN != None) :
        plt.savefig(outN)

    return ax

def plot_dendrogram(file_name,transform='plain',shift=0,c='k',radius=True,rm=20000.0,ra=200,outN=None) :
    global C, RM, RA, max_width, max_height # n.a.s.t.y.
    '''
    Generate a dendrogram from an SWC file. The SWC has to be formatted with a "three point soma"
    '''
    swc_tree = btmorph.STree2()
    swc_tree = swc_tree.read_SWC_tree_from_file(file_name)
    RM = rm
    RA = ra
    C = c
    max_height = 0
    max_width = 0
    plt.clf()
    print 'Going to build the dendrogram. This might take a while...'
    ttt = time.time()
    _expand_dendrogram(swc_tree.get_root(),swc_tree,shift,0,radius=radius,transform=transform)
    if(transform == 'plain') :
        plt.ylabel('L (micron)')
    elif(transform == 'lambda') :
        plt.ylabel('L (lambda)')
    print (time.time() - ttt), ' later the dendrogram was finished. '

    print 'max_widht=%f, max_height=%f' % (max_width,max_height)
    x_bound = (max_width / 2.0) + (0.1*max_width)
    max_y_bound = max_height + 0.1*max_height
    plt.axis([-1.0*x_bound,x_bound,-0.1*max_height,max_y_bound])

    plt.plot([x_bound,x_bound],[0,100],'k', linewidth=5) # 250 for MN, 100 for granule

    frame1 = plt.gca()
    frame1.axes.get_xaxis().set_visible(False)
    frame1.axes.get_yaxis().set_visible(False)
    
    if(outN != None) :
        plt.savefig(outN)

def _expand_dendrogram(cNode,swc_tree,off_x,off_y,radius,transform='plain') :
    global max_width,max_height # middle name d.i.r.t.y.
    '''
    Gold old fashioned recursion... sys.setrecursionlimit()!
    '''
    place_holder_h = H_SPACE
    max_degree = swc_tree.degree_of_node(cNode)
    required_h_space = max_degree * place_holder_h
    start_x = off_x-(required_h_space/2.0)
    if(required_h_space > max_width) :
        max_width = required_h_space

    
    children = cNode.get_child_nodes()

    if swc_tree.is_root(cNode) :
        print 'i am expanding the root'
        children.remove(swc_tree.get_node_with_index(2))
        children.remove(swc_tree.get_node_with_index(3))
    
    for cChild in  children :
        l = _path_between(swc_tree,cChild,cNode,transform=transform)
        r = cChild.get_content()['p3d'].radius

        cChild_degree = swc_tree.degree_of_node(cChild)
        new_off_x = start_x + ( (cChild_degree/2.0)*place_holder_h )
        new_off_y = off_y+(V_SPACE*2)+l
        r = r if radius  else 1
        plt.vlines(new_off_x,off_y+V_SPACE,new_off_y,linewidth=r,colors=C)
        if((off_y+(V_SPACE*2)+l) > max_height) :
            max_height = off_y+(V_SPACE*2)+l

        _expand_dendrogram(cChild,swc_tree,new_off_x,new_off_y,radius=radius,transform=transform)

        start_x = start_x + (cChild_degree*place_holder_h)
        plt.hlines(off_y+V_SPACE,off_x,new_off_x,colors=C)

def _path_between(swc_tree,deep,high,transform='plain') :
    path = swc_tree.path_to_root(deep)
    pl = 0
    pNode = deep
    for node in path[1:] :
        pPos = pNode.get_content()['p3d'].xyz
        cPos = node.get_content()['p3d'].xyz
        pl = pl + np.sqrt(np.sum((cPos-pPos)**2))
        #pl += np.sqrt( (pPos.x - cPos.x)**2 + (pPos.y - cPos.y)**2 + (pPos.z - cPos.z)**2 )
        pNode = node
        if(node == high) : break
        
    if(transform == 'plain'):
        return pl
    elif(transform == 'lambda') :
        DIAM = (deep.get_content()['p3d'].radius*2.0 + high.get_content()['p3d'].radius*2.0) /2.0 # naive...
        c_lambda = np.sqrt(1e+4*(DIAM/4.0)*(RM/RA))
        return pl / c_lambda
