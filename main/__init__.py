'''
Main Class

'''
import numpy as np
from main.Common import Histogram,SymbolicObject,normalize
from main.cluster.classicalHCC import HCC
from main.cluster.MHCC import MHCC
from main.cluster.dendogram import drawDendowithDiss
import main.PCA.classicalPCA as cPCA
import main.PCA.vertices_optimized as voPCA
import main.PCA.vertexPCA as vPCA
# import main.PCA.WPCA as wPCA #only import when needed. Used WPCA library messes with graph drawing :/
import main.PCA.quantilePCA as qPCA
import main.visual.ACGs as ACGs
import main.visual.Zoomstar as star
import main.visual.ShapeEncoding as se
from main.visual.Dissimilarity_based_coloring import getHues
import main.visual.MST as MST
import main.regression.predictKNN as KNN


def callMethod(X,objects,titles,type,macro_mode,Fmin=None,Fmax=None,nMode=0,cat=[],bins=[]):
    nre= len(X)
    nrf = len(titles)
    
    quantiles=[0,0.1,0.25,0.5,0.75,0.9,1]
#     quantiles=[0,1]

    nrq=len(quantiles)
    
    #transform data to symbolic object, normalize
    aList2,Fmin,Fmax=normalize(X,objects,titles,type,mode=nMode,Fmin=Fmin,Fmax=Fmax)
    
    #print all data
#     for i in range(nre):
#         print(aList2[i].getHistRounded(3))

    #clustering - microscopic and macroscopic 
#     diss,list=MHCC(aList2,nre,nrf,objects.copy(),quantiles,till=1,r=0)
#     drawDendowithDiss(diss,objects)
    
#     diss,st=HCC(aList2,nrf,nre,1,mode=macro_mode)
#     drawDendowithDiss(diss,objects)

    #PCA
    
    #centers by classical PCA
#     cPCA.mainMethod(nre,nrf,objects,titles,aList=aList2)
    
    #vertex method
#     vPCA.vertex(aList2,nre,nrf,objects,titles,type)
    
    #optimized vertex
#     voPCA.optimized_vertex(aList2,nre,nrf,objects,titles,type)

    #weighted PCA
#     wPCA.WPCA_my(aList2,nre,nrf,objects,titles,type)
    
    #quantile PCA
#     qPCA.mainMethod(aList2,nre,nrf,objects,quantiles)


    #VISUALS
    #Acumulated concept graphs
#     line_width=0.2
#     dot_size=1
#     ACGs.drawSimple_for_quantiles(aList2,objects,line_width,dot_size,quantiles)
#     ACGs.drawMinMax_for_quantiles(aList2,objects,line_width,dot_size,quantiles)
#     ACGs.drawUniliteral_for_quantiles(aList2,objects,line_width,dot_size,quantiles)
#     ACGs.drawBiliteral_for_quantiles(aList2,objects,line_width,dot_size,quantiles)
#     ACGs.drawQVACG_for_quantiles(aList2,objects,line_width,dot_size,quantiles)
#     ACGs.drawSQVACG_for_quantiles(aList2,objects,line_width,dot_size,quantiles)
#     ACGs.drawTACG_for_quantiles(aList2,objects,line_width,dot_size,quantiles)
#     ACGs.drawACC_for_quantiles(aList2,objects,line_width,dot_size,quantiles)
#     ACGs.drawACC_FF_for_quantiles(aList2,titles,line_width,dot_size,quantiles)
#     ACGs.drawFingertips_for_quantiles(aList2,objects,line_width,dot_size,quantiles)
#     ACGs.drawQuantiles(aList2,titles,objects,quantiles,line_width,dot_size,0,4,0.9,0.4)
    
    #ZOOMSTAR
    # zoomstar_conf indicates what is included in graph 
    #0 - 1-span, 2-color intensity , 3- quantiles
    #1 - 1-plot eq interval
    #2 - 1- plot median line with 100% color intensity 2-median with black color
    #3 - 1- axis texts are includes, 0 - not included
    
    
    hues=[] #if different colors are desired for all objects, values from 0 to 1 should be assigned for every object. Empty list is replaced with 0.5 in method
#     hues=getHues(aList2,nre,nrf,objects,titles,type,0,0,quantiles) 


    #QUANTILE ZOOMSTARS
    zoomstar_conf=[3,0,1,1]
    
    #EQUIVALENT INTERVAL ZOOMSTARS
#     zoomstar_conf=[1,1,1,1]
    
    #EQUIVALENT INTERVAL ZOOMSTARS (for interval without extra span)
#     zoomstar_conf=[0,1,1,1]

    #COLOR INTENSITY ZOOMSTARS
#     zoomstar_conf=[2,0,2,1]
    
    #single  zoomstar
#     object_nr=0
#     star.drawZoomStar(aList2,nrf,nre,objects,titles,zoomstar_conf,Fmin,Fmax,object_nr,hues=hues,labels=cat,quantiles=quantiles)
    #all zoomstars
#     cols=5 #how many zoomstars is a row
#     star.drawZoomStar(aList2,nrf,nre,objects,titles,zoomstar_conf,Fmin,Fmax,cols=cols,hues=hues,labels=cat,quantiles=quantiles)

#SHAPE ENCODING
#     new_order=[x for x in range(nre)] #for reordering objects
#     se.drawShapeCoding(aList2,objects,titles,bins,new_order)
    
    
#SCATTER PLOT OF PCA
#     D=[]
#     for i in range(nre):
#         abi=[]
#         for j in range(nrf):
#             abi.append(aList2[i].getHistogram(j).get_median())
#         D.append(abi)
#     #PERFORM PCA AND GET POINTS
#     P=cp.getPCAPoints(D,nre,nrf,titles,objects,type)  
#     star.drawScatterPlot(aList2,nrf,nre,objects,titles,type,Fmin,Fmax,P,quantiles,hues=hues,labels=cat)


#MINIMUM SPANNING TREE
#     MST.MST(aList2,nre,nrf,objects,quantiles)

#predict missing values
#     i=0 #object to be predicted
#     qu=[0] #index of quantiles to be predicted. One quantile at the time gives best results
#     f=0 #feature to be predicted
#     k=3 #how many closest neighbours
#     
#     q=quantiles #quantiles used for corrections
#     FS=[1 for x in range(nrf)] #features to be used for predictions
#     
#     #if mode is 1, corrections are used
#     #t is prediction
#     k,t=KNN.predict(aList2, i, q,qu,FS, nre, nrf, f, k,mode=0)


if __name__ == '__main__':
    
    Fmin=[]
    Fmax=[]
    cat=[]
    bins=[]
    nMode=0
    
#    OILS DATASET
    objects=['Linsead','Perilla','Cotton','Sesame','Camellia','Olive','Beef','Hog']
    titles = ['Specific Gravity','Freezing Point','Iodine Value','Saponification value','Major Fatty Acids']
    type=[0,0,0,0,0]
    macro_mode=4
    X = np.loadtxt("data/oils.txt", dtype='str', delimiter='\n')
#     cat=[[],[],[],[],["L", "Ln", "O", "P", "M","S","A","C","Lu"]]
#     bins=[10,10,10,10,9]
    
    #CITY DATASET
#     objects=['Amsterdam','Athens','Bahrain','Bombay','Cairo','Calcutta','Colombo','Copenhagen','Dubai','Frankfurt','Geneva','Hong_Kong','Kuala_Lumpur','Lisbon','London','Madras','Madrid','Manila','Mauritius','Mexico_City','Moscow','Munich','Nairobi','New_Delhi','New_York','Paris','Rome','San_Francisco','Seoul','Singapore','Stockholm','Sydney','Tehran','Tokyo','Toronto','Vienna','Zurich']
#     titles=['jan','feb','mar','apr','may','jun','jul','aug','sept','oct','nov','dec']
#     type=[0,0,0,0,0,0,0,0,0,0,0,0]
#     macro_mode=4
#     X = np.loadtxt("data/city.txt", dtype='str', delimiter='\n')

# #    HARDWOOD DATASET
    objects=['ACER_EAST','ACER_WEST','ALNUS_EAST','ALNUS_WEST','FRAXINUS_EAST','FRAXINUS_WEST','JUGLANS_EAST','JUGLANS_WEST','QUERCUS_EAST','QUERCUS_WEST']
    titles=['ANNT','JANT','JULT','ANNP','JANP','JULP','GDC5','MITM']
    type=[0,0,0,0,0,0,0,0,0]
    macro_mode=0
    X = np.loadtxt("data/hardwood.txt", dtype='str', delimiter='\n')
# 
# #    STATES DATASET
#     objects=['Alabama','Arizona','Arkansas','California','Colorado','Connecticut','Delaware','Florida','Georgia','Idaho','Illinois','Indiana','Iowa','Kansas','Kentucky','Louisiana','Maine','Maryland','Massachusetts','Michigan','Minnesota','Mississippi','Missouri','Montana','Nebraska','Nevada','New_Hampshire','New_Jersey','New_Mexico','New_York','North_Carolina','North_Dakota','Ohio','Oklahoma','Oregon','Pennsylvania','Rhode_Island','South_Carolina','South_Dakota','Tennessee','Texas','Utah','Vermont','Virginia','Washington','West_Virginia','Wisconsin','Wyoming']
#     titles=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sept','Oct','Nov','Dec']
#     type=[0,0,0,0,0,0,0,0,0,0,0,0]
#     macro_mode=0
#     nMode=1
#     X = np.loadtxt("data/states.txt", dtype='str', delimiter='\n')
    
    
    callMethod(X,objects,titles,type,macro_mode,Fmin=Fmin,Fmax=Fmax,nMode=nMode,cat=cat,bins=bins)