'''
Instead of training the whole dataset by all vertex points (fractorial!!!), only min/max points are used
The coverage with 2 or 3 PCs is usually better with this approach than when using all vertex points
'''

from sklearn.preprocessing import scale
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import numpy as np
from main.PCA.PCA_common import itterateAll,itterateMin,itterateMax
from main.PCA.drawPCA import plot_resultsMAIN,plot_results


 
def prepare(aList,nre,nrf,objects,titles,type):
        
    
    resAll=[]
    resMin=[]
    resMax=[]
    for i in range(nre):
        itterateAll(aList[i],[],0,nrf,resAll)
        itterateMin(aList[i],[],0,nrf,resMin)
        itterateMax(aList[i],[],0,nrf,resMax)
   
    fitRes=[]
    for i in range(nre):
        fitRes.append(resMin[i])
        fitRes.append(resMax[i])

#     print(str(len(resAll))+" "+str(len(fitRes)))

    
    pca = PCA(n_components=nrf)
    fitRes=np.array(fitRes)
    sd=np.std(fitRes, axis=0)
    mu = np.mean(fitRes, axis=0)
    fitRes = scale(fitRes)
  
    pca.fit(fitRes)
    resAll=np.array(resAll)
    

    
    resAll-=mu
    resAll/=sd
    newX=pca.transform(resAll)
    varexp=pca.explained_variance_ratio_

    backRes=np.dot(pca.transform(resAll)[:,:2], pca.components_[:2,:])
#     backRes*=sd
#     backRes+= mu
    
#     newX=np.array(newX)
# X,newX,Y
    var_exp=varexp

    plot_results(resAll,pca,newX,nre,nrf,objects,titles,backRes=newX,Y=backRes,config=[1,1,1])
    
    #1 plot
    plot_resultsMAIN(pca,newX,objects,nrf)

def optimized_vertex(aList,nre,nrf,objects,titles,type):
    prepare(aList,nre,nrf,objects,titles,type)
