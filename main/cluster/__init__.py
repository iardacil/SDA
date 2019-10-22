'''
Main Class
Two datasets Oils and Hardwood are included
Two methods - classical HCC for symbolic data and MHCC can be called
Dendograms can be drawn according to dissimilarity matrixes produced by clustering

@author: Kadri Umbleja
'''
import numpy as np
from main.Common import Histogram,SymbolicObject,normalize
from cluster.classicalHCC import HCC
from cluster.MHCC import MHCC
from cluster.dendogram import drawDendowithDiss


def callMethod(X,objects,titles,type,macro_mode):
    nre= len(X)
    nrf = len(titles)
    
    quantiles=[0,0.1,0.25,0.5,0.75,0.9,1]
    quantiles=[0,1]

    nrq=len(quantiles)
    
    #transform data to symbolic object, normalize
    aList2,Fmin,Fmax=normalize(X,objects,titles,type)
    

    
    diss=MHCC(aList2,nre,nrf,objects,quantiles,till=1)
    drawDendowithDiss(diss,objects)
    
    diss=HCC(aList2,nrf,nre,1,mode=macro_mode)
    drawDendowithDiss(diss,objects)