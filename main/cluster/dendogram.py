'''
Class for drawing dendrograms based on dissimilarity matrix
@author: Kadri Umbleja
'''

from cluster.classicalHCC import HCC
import scipy as sp
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform
from main.Common import Histogram,SymbolicObject,apply_permutation
import numpy as np
from scipy.cluster.hierarchy import dendrogram, to_tree,leaves_list


#class that draws dendograms based on HCC
#mode defines metric 0-histogram, 4-interval
def drawDendo(aList,nre,nrf,objects,titles,type,mode=4,till=1,l=[]):
    

    diss,extra=HCC(aList,nrf,nre,till,mode=mode,Fmin=[],Fmax=[],l=l)

    drawDendowithDiss(diss,objects)
    
def flatten(l):
    return [item for sublist in l for item in sublist]

def getOptimaLeave(diss):  
    diss=sp.array(diss)
    linkage_matrix = linkage(squareform(diss), "complete",optimal_ordering=True)
    return leaves_list(linkage_matrix)

def drawDendowithDiss(diss,objects,notes=[],filename="",tt="",nr_of_ann=0):
    

    diss=sp.array(diss)
    linkage_matrix = linkage(squareform(diss), "complete",optimal_ordering=True)
    dend =dendrogram(linkage_matrix, labels=objects,leaf_rotation=90.,leaf_font_size=8.,color_threshold=0, above_threshold_color='black')
    o2=objects=apply_permutation(objects,leaves_list(linkage_matrix)) #get objects from optimal order
    plt.subplots_adjust(bottom=0.3)
    plt.title(tt)
    
    
    #write notes on the node points
    if(len(notes)>0):
        X = flatten(dend['icoord'])
        Y = flatten(dend['dcoord'])
        leave_coords = [(x,y) for x,y in zip(X,Y) if y==0]
        order = np.argsort([x for x,y in leave_coords])
        id_to_coord = dict(zip(dend['leaves'], [leave_coords[idx] for idx in order])) # <- main data structure
        
        children_to_parent_coords = dict()
        for i, d in zip(dend['icoord'], dend['dcoord']):
            x = (i[1] + i[2]) / 2
            y = d[1] # or d[2]
            parent_coord = (x, y)
            left_coord = (i[0], d[0])
            right_coord = (i[-1], d[-1])
            children_to_parent_coords[(left_coord, right_coord)] = parent_coord
        
        # traverse tree from leaves upwards and populate mapping ID -> (x,y)
        root_node, node_list = to_tree(linkage_matrix , rd=True)
        ids_left = range(len(dend['leaves']), len(node_list))
        
        while len(ids_left) > 0:
        
            for ii, node_id in enumerate(ids_left):
                node = node_list[node_id]
                if (node.left.id in id_to_coord) and (node.right.id in id_to_coord):
                    left_coord = id_to_coord[node.left.id]
                    right_coord = id_to_coord[node.right.id]
                    id_to_coord[node_id] = children_to_parent_coords[(left_coord, right_coord)]
        
            ids_left = [node_id for node_id in range(len(node_list)) if not node_id in id_to_coord]
        
        kokku=len(dend['icoord'])
        ax = plt.gca()
        for node_id, (x, y) in id_to_coord.items():
            if not node_list[node_id].is_leaf():
                if(nr_of_ann==0):
                    if(node_id<len(objects)+kokku-1):
                        ax.plot(x, y, 'ro')
                        ax.annotate(notes[node_id-len(objects)],  (x, y), xytext=(0,8), textcoords='offset points',  va='bottom', ha='right', rotation=90) 
                    else:
                        ax.annotate(notes[node_id-len(objects)], (x, y), xytext=(0,12), textcoords='offset points',  va='top', ha='center') 
                elif(node_id>=kokku+len(objects)-nr_of_ann):
                    if(node_id<len(objects)+kokku-1):
                        ax.plot(x, y, 'ro')
                        ax.annotate(notes[node_id-len(objects)],  (x-15, y), xytext=(0,8), textcoords='offset points',  va='bottom', ha='left', size=8, rotation=90) 
              
                    
                
        dend['node_id_to_coord'] = id_to_coord
    
    
    
    if(len(filename)>0):
        plt.savefig('graphs/'+filename+'.png', format='png', dpi=300)
    else:
        plt.savefig('graphs/dendo.png', format='png', dpi=300)
        plt.show()
    plt.clf()