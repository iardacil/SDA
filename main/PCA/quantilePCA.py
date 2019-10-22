'''
Using Spearmann
'''



import numpy as np
from scipy import stats
from main.PCA.drawPCA import plot_resultsQuantile
from operator import itemgetter 


def mainWithMatrix(D,nre,nrf,nrq,objects):

    aList2=D
    aList3=[[0 for x in range(nrf)] for y in range(nrf)] #nrf x nrf
    for i in range(nrf):
        for j in range(nrf):
            aList3[i][j]=stats.spearmanr(aList2[i], aList2[j])[0]
    

    eig_vals, eig_vecs = np.linalg.eig(aList3)
     
#     print('Eigenvectors \n%s' %eig_vecs)
#     print('\nEigenvalues \n%s' %eig_vals)
 
    # Make a list of (eigenvalue, eigenvector) tuples
    eig_pairs = [(np.abs(eig_vals[i]), eig_vecs[:,i]) for i in range(len(eig_vals))]
     
    eig_pairs=sorted(eig_pairs,key=itemgetter(0), reverse=True)
 


    tot = sum(eig_vals)
    var_exp = [(i / tot)*100 for i in sorted(eig_vals, reverse=True)]
    
    matrix_w = np.hstack((eig_pairs[0][1].reshape(nrf,1), eig_pairs[1][1].reshape(nrf,1)))
    V=np.array(aList2)
    Y_sklearn = V.transpose().dot(matrix_w)



    plot_resultsQuantile(Y_sklearn,nrq,objects,var_exp)
    


def mainMethod(aList,nre,nrf,objects,quantiles=[0,0.1,0.25,0.5,0.75,0.9,1]):

    nreO=nre
    nre=nre*len(quantiles) #every quantile for every object gets its own "object"
    aList2=[[0 for x in range(nre)] for y in range(nrf)]  #nre x nrf
    c=0
    for i in range(nreO):
        for j in range(len(quantiles)):
            for n in range(nrf):
                val=aList[i].getHistogram(n).get_quantile(quantiles[j])
                aList2[n][c]=val
            c=c+1
    
    mainWithMatrix(aList2,nre,nrf,len(quantiles),objects)
    