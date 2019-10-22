
import numpy as np
from sklearn.preprocessing import scale
import seaborn; seaborn.set() # set plot style
from wpca import  WPCA
from main.PCA.PCA_common  import itterateAll,itterateMin,itterateMax
from main.PCA.drawPCA import plot_resultsMAIN,plot_results



#actual PCA happens here
def do_PCA(fitRes, allRes,n_features,n_samples,objects,titles,weights=None, ncomp=2):
    # Compute the standard/weighted PCA
    if weights is None:
        kwds = {}
    else:
        kwds = {'weights': weights} #uses compactness for weight
    
    #weighted PCA fitting
    pca = WPCA(n_components=n_features).fit(fitRes, **kwds)
    
    
    newX=pca.transform(fitRes,**kwds) #with weights
    W=pca.transform(fitRes) #without weights
    V2=pca.transform(allRes)

    varexp=pca.explained_variance_ratio_
    
#for back fitting
    pca2=WPCA(n_components=ncomp).fit(fitRes, **kwds)

    Y=pca2.inverse_transform(newX[:,:ncomp])
    W2=pca2.inverse_transform(W[:,:ncomp])
    V3=pca2.inverse_transform(V2[:,:ncomp])
    

    plot_results(allRes, pca,newX,n_samples,n_features,objects,titles,backRes=V2,Y=V3,config=[1,1,1],pperelem=2,imp1=0,imp2=2)
#     plot_resultsMAIN(pca,newX,objects,n_features,npoint=2)
    


#prepare data
def prepare(aList,nre,nrf,objects,titles,type):
    n_samples = nre
    n_features = nrf
    nrf=n_features
    nre=n_samples

    
    #make symbolic objects, compactness are used for weights - twice needed for non vertwex method
    W=[]
    for i in range(nre):
        abi=[]
        for j in range(nrf):
            abi.append(aList[i].getHistogram(j).compactness())
        W.append(abi)
        W.append(abi)
    
    #datasets for fitting find, min and max sets and all
    resAll=[]
    resMin=[]
    resMax=[]
    for i in range(nre):
        itterateAll(aList[i],[],0,nrf,resAll)
        itterateMin(aList[i],[],0,nrf,resMin)
        itterateMax(aList[i],[],0,nrf,resMax)
   
    #set for fitting use min and max
    fitRes=[]
    for i in range(nre):
        fitRes.append(resMin[i])
        fitRes.append(resMax[i])
        
        
    fitRes=np.array(fitRes)
    sd=np.std(fitRes, axis=0)
    mu = np.mean(fitRes, axis=0)
    fitRes = scale(fitRes)
    
    #all vertexes
    resAll=np.array(resAll)
    resAll-=mu
    resAll/=sd
        
#     P=fitRes
#     V=resAll
#     
#     X=P
#     
#     sd=np.std(fitRes, axis=0)
#     mu = np.mean(fitRes, axis=0)
#        
#     resAll-=mu
#     resAll/=sd
#     
#     X_scaled=scale(P) 
    
    do_PCA(fitRes,resAll,n_features,n_samples,objects,titles,W)

def WPCA_my(aList,nre,nrf,objects,titles,type):
    prepare(aList,nre,nrf,objects,titles,type)