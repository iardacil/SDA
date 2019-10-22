from sklearn.preprocessing import scale
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import numpy as np
from main.PCA.PCA_common import itterateAll,itterateMin,itterateMax
from main.PCA.drawPCA import plot_resultsMAIN,plot_results


 
def prepare(aList,nre,nrf,objects,titles,type):
        
    
    resAll=[]

    for i in range(nre):
        itterateAll(aList[i],[],0,nrf,resAll)


    pca = PCA(n_components=nrf)
    resAll=np.array(resAll)
    resAll = scale(resAll)
    pca.fit(resAll)


    newX=pca.transform(resAll)
    varexp=pca.explained_variance_ratio_

    backRes=np.dot(pca.transform(resAll)[:,:2], pca.components_[:2,:])

    var_exp=varexp

    plot_results(resAll,pca,newX,nre,nrf,objects,titles,backRes=newX,Y=backRes,config=[1,1,1])
    
    #1 plot
    plot_resultsMAIN(pca,newX,objects,nrf)

def vertex(aList,nre,nrf,objects,titles,type):
    prepare(aList,nre,nrf,objects,titles,type)
