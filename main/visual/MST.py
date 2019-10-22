
import networkx as nx
import matplotlib.pyplot as plt  
import main.PCA.quantilePCA as sp
  
#Class to represent a graph 
class Graph: 
  
    def __init__(self,vertices): 
        self.V= vertices #No. of vertices 
        self.graph = [] # default dictionary  
                                # to store graph 
          
   
    # function to add an edge to graph 
    def addEdge(self,u,v,w): 
        self.graph.append([u,v,w]) 
  
    # A utility function to find set of an element i 
    # (uses path compression technique) 
    def find(self, parent, i): 
        if parent[i] == i: 
            return i 
        return self.find(parent, parent[i]) 
  
    # A function that does union of two sets of x and y 
    # (uses union by rank) 
    def union(self, parent, rank, x, y): 
        xroot = self.find(parent, x) 
        yroot = self.find(parent, y) 
  
        # Attach smaller rank tree under root of  
        # high rank tree (Union by Rank) 
        if rank[xroot] < rank[yroot]: 
            parent[xroot] = yroot 
        elif rank[xroot] > rank[yroot]: 
            parent[yroot] = xroot 
  
        # If ranks are same, then make one as root  
        # and increment its rank by one 
        else : 
            parent[yroot] = xroot 
            rank[xroot] += 1
  
    # The main function to construct MST using Kruskal's  
        # algorithm 
    def KruskalMST(self,objects,G,aList2): 
        nrf=aList2[0].nrf
        
  
        result =[] #This will store the resultant MST 
  
        i = 0 # An index variable, used for sorted edges 
        e = 0 # An index variable, used for result[] 
  
            # Step 1:  Sort all the edges in non-decreasing  
                # order of their 
                # weight.  If we are not allowed to change the  
                # given graph, we can create a copy of graph 
        self.graph =  sorted(self.graph,key=lambda item: item[2]) 
        
#         for i in range(len(self.graph)):
#             print(objects[self.graph[i][0]]+" "+objects[self.graph[i][1]]+" "+str(self.graph[i][2]))
  
        parent = [] ; rank = [] 
  
        # Create V subsets with single elements 
        for node in range(self.V): 
            parent.append(node) 
            rank.append(0) 
        i=0
        # Number of edges to be taken is equal to V-1 
        while e < self.V -1 : 
  
            # Step 2: Pick the smallest edge and increment  
                    # the index for next iteration 
            u,v,w =  self.graph[i] 
            i = i + 1
#             print(objects[u]+" "+objects[v]+" "+str(w))
            x = self.find(parent, u) 
            y = self.find(parent ,v) 
#             print(str(x)+" "+str(y))
  
            # If including this edge does't cause cycle,  
                        # include it in result and increment the index 
                        # of result for next edge 
            if x != y: 
                e = e + 1     
                result.append([u,v,w]) 
                self.union(parent, rank, x, y)             
            # Else discard the edge 
  
        # print the contents of result[] to display the built MST 
#         print("Following are the edges in the constructed MST")
        
        D=[]
        obj=[]
        for u,v,weight  in result: 
            G.add_edge(objects[u], objects[v], length = round(weight,3)) 
            
            #print str(u) + " -- " + str(v) + " == " + str(weight) 
#             print ("%s -- %s == %f" % (objects[u],objects[v],weight)) 
#             i=int(objects[u])-1
#             j=int(objects[v])-1
# #             print(aList2[i].name+"-"+aList2[j].name, end=" ")
#             new=[]
#             obj.append(aList2[i].name+"-"+aList2[j].name)
#             for f in range(nrf):
#                 new.append(aList2[i].combineR(aList2[j]).getHistogram(f).getSpan())
# #                 print(aList2[i].combineR(aList2[j]).getHistogram(f).getSpan(), end=" ")
# #             print()
#             D.append(new)
        
#         sp.mainWithMatrix(D,len(D),nrf,2,obj)
        
        pos = nx.spring_layout(G)
        nx.draw(G, pos, with_labels = True)  #with_labels=true is to show the node number in the output graph
        edge_labels = nx.get_edge_attributes(G,'length')
        nx.draw_networkx_edge_labels(G, pos, edge_labels = edge_labels, font_size = 11) #prints weight on all the edges
        
        for u,v,weight  in result:    
            nx.draw_networkx_edges(G, pos, edgelist = [(objects[u], objects[v])], width = 2.5, alpha = 0.6, edge_color = 'r')
        plt.savefig("graphs/MST.png")
        return G
 

def MST(aList2,nre,nrf,objects,quantiles,till=1): 
    # Driver code 
    gg = Graph(nre) 
    nrq=len(quantiles)
    
    G= nx.Graph()

#     objects=[x for x in range(nre)]
    X=[[0 for x in range(nre)] for y in range(nre)]
    for i in range(nre):
        for j in range(nre):
            if(i>j):
                diss=0
                for f in range(nrf):
                    for n in range(nrq):
                        diss+=abs(aList2[i].getHistogram(f).get_quantile(quantiles[n])-aList2[j].getHistogram(f).get_quantile(quantiles[n]))
                gg.addEdge(i, j,diss/(nrf*nrq))
                
#                 X[i][j]=diss/(nrf*nrq)
#                 X[j][i]=diss/(nrf*nrq)

    

    gg.KruskalMST(objects,G,aList2) 

    
    plt.show()
    




    
