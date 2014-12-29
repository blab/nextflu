import dendropy, datetime
import numpy as np



def color_BioTree_by_attribute(T,attribute, vmin=None, vmax = None, missing_val='min', transform = lambda x:x, cmap=None):
    '''
    simple function that assigns a color to each node in a biopython tree
    the color can be determined by any attribute of the nodes. missing attributes will be
    determined from the children, all children are assumed to have the attribute
    in addition, the attribute can be transformed for example by taking the log
    parameters:
    T               -- BioPython tree
    attribute       -- name of the attribute that is to be used to color the tree. 
    vmin            -- lower offset that is subtracted
    vmax            -- values are scaled as (val-vmin)/(vmax-vmin)
    missing val     -- if the attribute does not exist is a particular node, 
                       the min, max, or mean of the children is used
    transform       -- function mapping float to float, e.g. log
    cmap            -- colormap to be used
    '''
    # make a list of tranformed data
    vals = [transform(t.__getattribute__(attribute)) for t in 
            T.get_terminals()+T.get_nonterminals() if attribute in t.__dict__]
    if vmin is None:  # if vmin or vmax is not provided, use min or max of data
        vmin = min(vals)
        print "Set vmin to",vmin
    if vmax is None:
        vmax = max(vals)
        print "Set vmax to",vmax
    if cmap is None:
        from matplotlib.cm import jet
        cmap=jet

    # assign function used to determine missing values from children
    if missing_val=='min':
        missing_val_func = min
    elif missing_val=='mean':
        missing_val_func = mean
    elif missing_val=='max':
        missing_val_func = max
    else:
        missing_val_func = min

    # loop over all nodes, catch missing values and assign
    for node in T.get_nonterminals(order='postorder'):
        if attribute not in node.__dict__:
            node.__setattr__(attribute, missing_val_func([c.__getattribute__(attribute) for c in node.clades]))
            print "node", node,"has no",attribute,"Setting to min:", node.__getattribute__(attribute)

    # map value to color for each node
    for node in T.get_terminals()+T.get_nonterminals():
        node.color = map(int, np.array(cmap((transform(node.__getattribute__(attribute))-vmin)/(vmax-vmin))[:-1])*255)



def date_to_day(date):
    if isinstance(date, datetime.date):
        return date.toordinal()
    elif isinstance(date,basestring):
        return datetime.datetime.strptime(date, '%Y-%m-%d').toordinal()
    elif isinstance(date, int):
        return date
    else:
        print "unknown date format", date
        return np.nan

def from_json(json):
    '''
    read a json dictionary and make a dendropy tree from it. 
    '''
    tree = dendropy.Tree()
    tree.get_from_string(';', 'newick')
    root = tree.leaf_nodes()[0]
    from_sub_json(json, root)
    root.edge_length=0.01
    return tree

def from_sub_json(json, node):
    '''
    recursively calls itself for all children of node and 
    builds up the tree. entries in json are added as node attributes
    '''
    for attr,val in json.iteritems():
        if attr=='children':
            for sub_json in val:
                child_node = dendropy.Node()
                from_sub_json(sub_json, child_node)
                node.add_child(child_node, edge_length = child_node.xvalue - node.xvalue)
        elif attr=='date':
            node.date = date_to_day(val)
        else:
            try:
                node.__setattr__(attr, float(val))
            except:
                node.__setattr__(attr, val)
    if len(node.child_nodes())==0:
        node.taxon= True
    else:
        node.taxon = False

def get_average_T2(tree, dt=365):
    '''
    calculates the average pairwise distance between nodes from a time interval dt
    repeats this for a bunch of dt intervals, returns the average
    tree  -- dendropy tree
    dt    -- time interval
    '''

    leaves = tree.leaf_nodes()
    leaves.sort(key=lambda x: date_to_day(x.date), reverse = True)    
    upper_cutoff = date_to_day(leaves[0].date)
    print len(leaves), upper_cutoff
    ndt = 4
    T2 = 0
    ntrees = 0
    for ii in xrange(ndt):
        leaf_set = set([leaf for leaf in leaves 
                if date_to_day(leaf.date)<upper_cutoff-ii*dt 
                and date_to_day(leaf.date)>=upper_cutoff-(ii+1)*dt])
        if len(leaf_set)>10:
            tmp_T2= avg_T2_of_leaf_set(tree, leaf_set) 
            T2 += tmp_T2
            ntrees += 1
    return T2/ntrees

def avg_T2_of_leaf_set(tree, leaf_set):
    '''
    walk through the nodes in postorder and for each edge increment the T2
    count by the edge length multiplied by the number of node pairs separated
    by this split
    tree     -- dendropy tree containing the distances
    leaf_set -- subset of the leaves of tree whose pairwise distance is returned
    '''
    T2 = 0
    nleaves = len(leaf_set)
    for node in tree.postorder_node_iter():
        # get number of leaves in leaf_set that are descendants of node
        ndescendants = len(leaf_set.intersection(node.leaf_nodes()))
        # multiply edge_length by the number of pairs across the split
        T2 += node.edge_length*ndescendants*(nleaves-ndescendants)
    return 2*T2/nleaves/(nleaves-1)

def calc_LBI(tree, tau, window = None):   
    '''
    traverses the tree in postorder and preorder to calculate the
    up and downstream tree length exponentially weighted by distance.
    then adds them as LBI
    tree -- dendropy tree for whose node the LBI is being computed
    tau  -- the memory time scale of tree length along branches of the tree
    '''
    import numpy as np
    # traverse the tree in postorder (children first) to mark alive nodes
    if window is None:
        for node in tree.postorder_node_iter():
            node.alive=True
    else:
        for node in tree.postorder_node_iter():
            if len(node.child_nodes())==0 and (node.date>=window[0] and node.date<window[1]):
                node.alive=True
            else:
                node.alive = any([c.alive for c in node.child_nodes()])

    # traverse the tree in postorder (children first) to calculate msg to parents
    for node in tree.postorder_node_iter():
        node.down_polarizer = 0
        node.up_polarizer = 0
        if node.alive:
            for child in node.child_nodes():
                if child.alive:
                    node.up_polarizer += child.up_polarizer
            bl =  node.edge_length/tau
            node.up_polarizer *= np.exp(-bl)
            node.up_polarizer += tau*(1-np.exp(-bl))

    # traverse the tree in preorder (parents first) to calculate msg to children
    for node in tree.preorder_internal_node_iter():
        for child1 in node.child_nodes():
            child1.down_polarizer = node.down_polarizer
            for child2 in node.child_nodes():
                if child1!=child2:
                    child1.down_polarizer += child2.up_polarizer
            
            bl =  child1.edge_length/tau
            child1.down_polarizer *= np.exp(-bl)
            child1.down_polarizer += tau*(1-np.exp(-bl))

    # go over all nodes and calculate the LBI (can be done in any order)
    for node in tree.postorder_node_iter():
        node.LBI = node.down_polarizer
        for child in node.child_nodes():
            node.LBI += child.up_polarizer
    #done

def calc_delta_LBI(tree, tau, cutoff1, cutoff2=None, log_delta = False):
    '''
    function that calculates the LBI for nodes up to cutoff1 and cutoff2
    and then takes the difference betweem the two values. 
    parameters:
    tree      --  dendropy tree
    cutoff1   --  the earlier cutoff as a datetime object
    cutoff2   --  the later cutoff. optional, by default 01/01/2050
    '''
    calc_LBI(tree, tau = tau, window = [datetime.datetime(1900,1,1).toordinal(),cutoff1.toordinal()])
    for node in tree.postorder_node_iter():
        node.delta_LBI = -node.LBI
    if cutoff2 is None:
        cutoff2 =  datetime.datetime(2050,1,1)
    calc_LBI(tree, tau = tau, window = [datetime.datetime(1900,1,1).toordinal(),cutoff2.toordinal()])
    for node in tree.postorder_node_iter():
        if log_delta:
            node.delta_LBI = -(node.LBI)/(node.delta_LBI+1e-10)
        else:
            node.delta_LBI += node.LBI


def add_dates_to_internal(tree):
    for node in tree.postorder_node_iter():
        if len(node.child_nodes()):
            node.date = min([x.date for x in node.child_nodes()])

def to_Biopython(tree):
    from Bio import Phylo
    from StringIO import StringIO
    from itertools import izip
    bT  = Phylo.read(StringIO(tree.as_newick_string()), 'newick')

    for new_leaf, old_leaf in izip(bT.get_terminals(), tree.leaf_nodes()):
        for attr,val in old_leaf.__dict__.iteritems():
            try:
                new_leaf.__setattr__(attr, float(val))
            except:
                new_leaf.__setattr__(attr, val)

    for new_leaf, old_leaf in izip(bT.get_nonterminals(order='postorder'), tree.postorder_internal_node_iter()):
        for attr,val in old_leaf.__dict__.iteritems():
            try:
                new_leaf.__setattr__(attr, float(val))
            except:
                new_leaf.__setattr__(attr, val)
    return bT

def test():
    from io_util import read_json 
    tree = from_json(read_json('data/20141228_tree_auspice.json'))
    print "calculate local branching index"
    T2 = get_average_T2(tree, 365)
    tau =  T2*2**-4
    print "avg pairwise distance:", T2
    print "memory time scale:", tau
    calc_delta_LBI(tree, tau, datetime.datetime(2014,1,1))
    bioTree = to_Biopython(tree)
    color_BioTree_by_attribute(bioTree, 'date')
    return bioTree


if __name__=='__main__':
    tree = test()
    from Bio import Phylo
    Phylo.draw(tree)

