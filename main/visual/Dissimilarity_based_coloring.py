
import math
import re
from main.cluster.classicalHCC import HCC
from main.cluster.MHCC import MHCC

#SCROLL DOWN TO MAIN PART



#FUNCTIONS
#The main method
def dissimilarityBasedColoring(str,DIS,objects):
    #builds tree structured list from string representation of dendrogram
    tree = make_tree(str)

    
    #step 1 - initial shade values are found
    l=traverse_tree(tree,DIS,objects)
    #step 2 - initial values are normalized

    #find the maximum value of initial shades
    max=get_max(l)
    #normalized list is returned
    n=normalize(l,max)

    #step 3- match objects with shades
    i=0
    
    #match objects and shades
    for each_item in tree:
        match_object_and_shade(each_item,n[i])
        i=i+1
        
    shades=[]
    #produce list matching the positions of original objects in list 'object'

    for object in objects:
        shades.append(dict[object])
#         print(object+"-")
#         print(dict[object])

    return shades

#traverse the tree using post-order
def traverse_tree(the_list,DIS,objects):
    #if is node with children
    if isinstance(the_list, list):
        #list used to store information for current node
        r=[]
        nr=[]
        #call borth children
        for each_item in the_list:
            #get sublist [[v_1, v_2],[ v_3, v_4]] 
            r.append(traverse_tree(each_item,DIS,objects))
            #find two objects, one from each sublist
            nr.append(get_item(each_item,objects))
        #find dissimilarity from dissimilarity matrix using two object found from sublists
        #due to hierarchical clustering, the dissimilarity of current node equals 
        #the distance between two elements from different subtrees
        dis=DIS[int(nr[0])][int(nr[1])]
        #if both children were leaves
        if(r[0]==0 and r[1]==0):
            return [0,dis]
        #if at least one child was sublist
        else:
            #find minimum and maximum values from sublists
            min1=get_min(r[0])
            min2=get_min(r[1])
            max1=get_max(r[0])
            max2=get_max(r[1])
            #find midpoints for both subtrees
            m1=(min1+max1)/2
            m2=(min2+max2)/2
            #add difference of dissimilarity between midpoints of both subtree
            diff=m1+dis-m2
#             diff=max1+dis
            
            #ALTERINTIVE FOR MORE SEPARATION BETWEEN CLUSTERS
#             max1=get_max(r[0])
#             diff=max1+dis
            
            #prepare results            
            r2=[]
            r2.append(r[0])
            #recalculate second subtree
            r2.append(reCalculate(r[1],diff))
            #return repositioned list reflecting dissimilarity between two clusters at current node
            return r2
    #is leaf
    else:
        return 0  
    
#find minimum value of tree structured list    
def get_min(the_list):
    #set large positive value for initial value
    min=9999999
    #recursively call children
    if isinstance(the_list, list):
        for each_item in the_list:
            m=get_min(each_item)
            #if m is smaller than current min, update it
            if(m<min):
                min=m
        #return smallest value in current tree
        return min 
    #is leaf
    else:
        return the_list 
    
#find maximum value of tree structured list
def get_max(the_list):
    #set large negative value for initial value
    max=-9999999
     #has children
    if isinstance(the_list, list):
        #recursively call children
        for each_item in the_list:
            m=get_max(each_item)
            #if m is larger than current max, update it
            if(m>max):
                max=m
        #return largest value in current tree
        return max 
    #is leaf
    else:
        return the_list 
    
#returns first found leaf of a subtree for finding dissimilarity 
def get_item(the_list,objects):
    if isinstance(the_list, list):
        return get_item(the_list[0],objects)
    else:
        #returns index of object
        return objects.index(the_list)
    
#normalizes the feature space [0,max] to [0,1]    
def normalize(the_list,max):
    #has children
    if isinstance(the_list, list):
        r=[]
        #recursively call children
        for each_item in the_list:
            r.append(normalize(each_item,max))
        return r 
    #is leaf
    else:
        #actual normalization on single object shade
        return the_list/max 
    
#recursively traverse the sub.tree and add diff to every element
def reCalculate(the_list,diff):
    l=[]
    #has children
    if isinstance(the_list, list):
        for each_item in the_list:
            #recursively call children
            l.append(reCalculate(each_item,diff))
        return l
    #is leaf
    else:
        return the_list+diff 

# turns string representation of dendogram stucture to tree structured list
def make_tree(data):
    items = re.findall(r"\(|\)|\w+", data)

    def req(index):
        result = []
        item = items[index]
        while item != ")":
            if item == "(":
                subtree, index = req(index + 1)
                result.append(subtree)
            else:
                result.append(item)
            index += 1
            item = items[index]
        return result, index

    return req(1)[0]

#assumes that all objects' names/indicators  are stored in list 'objects'
#matches the one dimensional list of objects with tree structured list with names/indicators and calculated shades
def match_object_and_shade(the_list,shades):
    i=0
    if isinstance(the_list, list):
        for each_item in the_list:
            match_object_and_shade(each_item,shades[i])
            i=i+1
    else:
        #save name/indicator with corresponding shade value
        dict[the_list]=shades
        

def travel_tree(the_list,parentC,level,b):
    global i
    if(i%2==0):
        C=parentC-(math.pow(b,level+1)+1/math.pow(2,level+1)) 
    else:
        C=parentC+(math.pow(b,level+1)+1/math.pow(2,level+1)) 
 
    if isinstance(the_list, list):
        for each_item in the_list:
            i=i+1
            travel_tree(each_item,C,level+1,b)
    else:
        dict[the_list]=C
        
def translateClustersB(str,objects,b=0):
    tree = make_tree(str)
    global i
    i=0

    for each_item in tree:
       travel_tree(each_item,0.5,1,b)
       i=i+1
    hues=[]
    for object in objects:
        hues.append(dict[object])
    return hues

#mode: 0-histogram HCC, 1-interval HCC
#method 0-my dissimilarity based method, 1- previous method
def getHues(aList,nre,nrf,objects,titles,type,mode=1,method=0,quantiles=[0,0.1,0.25,0.5,0.75,0.9,1]):
    global dict

    till=1
    
    if(mode==0): #quantile microcsopic
        DIS,cluster_str=MHCC(aList,nre,nrf,objects,quantiles,till,r=0)
    elif(mode==1): #interval macroscopic
        DIS,cluster_str=HCC(aList,nrf,nre,till,mode=4)
    cluster_str=cluster_str[0]
    
    #call main function, find shade values
    print(cluster_str)
    dict= {}
    if(method==0):
        shades=dissimilarityBasedColoring(cluster_str,DIS,objects)
    elif(method==1):
        shades=translateClustersB(cluster_str,objects)

    return shades


    
    
