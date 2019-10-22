'''
classical PCA
'''


from sklearn.decomposition import PCA as sklearnPCA

import numpy as np
from sklearn.preprocessing import scale
from main.PCA.drawPCA import plot_resultsMAIN,plot_results
from main.PCA.PCA_common import itterateAll


#to apply method with already formed points
def getPCAPoints(D,nre,nrf,objects):
    return mainMethod(nre,nrf,objects,D=D)

#mode=0: draw graph
#mode=1: return center points
def mainMethod(nre,nrf,objects,titles,mode=0,aList=[],D=[]):

    n_samples = nre
    n_features = nrf
    n=objects

    
    #get center points
    X=[]
    if(len(D)==0):
        for i in range(nre):
            abi=[]
            for j in range(nrf):
                val=aList[i].getHistogram(j).get_median()
                abi.append(val)
            X.append(abi)
    else:X=D
        
    #will need these later
    sd=np.std(X, axis=0)
    mu2 = np.mean(X, axis=0)
    #scale the centers
    X=scale(X)

    #apply classical PCA
    pca = sklearnPCA(n_components=n_features)
    newX = pca.fit_transform(X)


    #eigenvector
#     print(pca.components_)
    
    #find all points
    backRes=[]
    backRes2=[]
    if(len(aList)>0):
        resAll=[]
        for i in range(n_samples):
            itterateAll(aList[i],[],0,nrf,resAll)

        resAll=np.array(resAll)
    
        #scale all vertex points accordingly
        resAll-=mu2
        resAll/=sd

        #apply transformation
        backRes=(pca.transform(resAll))
        
        #another option: use 2 PCs to get back to original cordinates
        backRes2=np.dot(pca.transform(resAll)[:,:2], pca.components_[:2,:])
    

    if(mode==0):
        #4 plots
        plot_results(X,pca,newX,n_samples,n_features,n,titles,backRes,Y=backRes2,config=[0,1,1])
        
        #1 plot
        plot_resultsMAIN(pca,newX,n,n_features,npoint=1)
        
    elif(mode==1):
        res=[]
        for j in range(2):
            abi=[]
            for i in range(nre):
                abi.append(newX[i][j])
            res.append(abi)
        return res
