'''
Class for HCC.
Two methods are available:
(mode=4) default method by Ichino based on concept sizes in feature space for interval valued data
(mode=0) modified version of previous for histogram valued data that considers histogram shape in more detail than just span/area in feature space
@author: Kadri Umbleja
'''

from main.Common import SymbolicObject, calcBack


#mode=4 rectangle/concept size method by Ichino
#mode=0 histogram version by Umbleja
def HCC(aList2,nrf,nre,till,mode=4,Fmin=[],Fmax=[],l=[],CT=[],r=0):
    aList=[]
    #make copy of symbolic objects as they will be destroyed durning HCC
    for i in range(nre):
        aList.append(SymbolicObject(aList2[i].hist,aList2[i].name,aList2[i].types,id=aList2[i].id))
    
    
    include_CT=0
    if(len(CT)>0): #make 
        include_CT=1  
    
    #for normalization if needed
    feature_len=[]
    if(len(Fmin)==0):
        for i in range(nrf):
            Fmin.append(0)
    if(len(Fmax)==0):
        for i in range(nrf):
            Fmax.append(1)
    
    if(len(l)==0):
        for i in range(nrf):
            l.append(1)
    
    cnt=0        
    for i in range(nrf):
        feature_len.append(Fmax[i]-Fmin[i])
        if(l[i]==1):
            cnt+=1
    
    BG=[0 for x in range(nrf)]
    for i in range(len(aList)):
        for f in range(nrf):
            BG[f]+=aList[i].getHistogram(f).get_quantile(0.5)
    for f in range(nrf):
        BG[f]/=len(aList)

    #dissimilarity matrix
    diss=[ [ 0 for y in range( nre ) ] for x in range( nre ) ]
    
    tss=[0 for x in range(nrf)]
    TSSG=0
    for i in range(len(aList)):
        for f in range(nrf):
            tss[f]+=(BG[f]-aList[i].getHistogram(f).get_quantile(0.5))**2
    for f in range(nrf):
        TSSG+=tss[f]
    
    notes=[]
    while(len(aList)>till):
        
#         calcDunn(aList,mode,nrf)

        curMin=float("inf")
        minj=-1
        mini=-1
        minW =None
        curr_size=len(aList)

        for i in range(curr_size):
            for j in range(curr_size):
                if(i>j):
                    if(mode==0): #histogram method
                        b = aList[i].combine(aList[j])
                        sum = b.calcCompactness()
                        if(curMin >sum):
                            curMin=sum
                            minW=b
                            mini=i
                            minj=j
                    elif(mode==4): #ichino method
                        b = aList[i].combineR(aList[j])
                        
                        sum=b.get_avg_Span(l)
                        if(curMin >sum):
                            curMin=sum
                            minW=b
                            mini=i
                            minj=j

#         print(minW.getHistRounded(3))
        print(str(curMin)+","+minW.getName())


        if(include_CT==1):
            rr=""
            for i in range(len(CT)):
                if(i>0):
                    rr+=", "
                rr+=""+str(round(calcBack(minW.getHistogram(CT[i]).get_nonZeroS(),Fmin,Fmax,CT[i]),3))+"<=F"+str(CT[i]+1)+"<="+str(round(calcBack(minW.getHistogram(CT[i]).get_nonZeroE(),Fmin,Fmax,CT[i]),3))
            notes.append(rr)
            
        #update dissimilarity matrix
        for i in range(len(aList[mini].list)):
            for j in range(len(aList[minj].list)):
                k = aList[mini].list[i].id;
                m = aList[minj].list[j].id;
                diss[k][m]=curMin
                diss[m][k]=curMin
        #update objects
        aList.append(minW)
        aList.pop(mini)
        aList.pop(minj)
        
#         if(len(aList)<=10):
#             calcVRC(aList,nrf,BG,nre,TSSG)
    
    if(r==0): #return dissimilarity matrix with names
        aNew=[]
        for i in range(len(aList)):
            aNew.append(aList[i].getName())
        return diss,aNew
    elif(r==2): #return dissimilarity matrix with node descriptions
        return diss,notes
    elif(r==1): #return the remaining list
        return aList

def calcDunn(aList,mode,nrf):
    
    mx=0
    for i in range(len(aList)): #for every cluster
        if(mode==0):
            s2=aList[i].get_avg_compactness()
        elif(mode==4): #ichino method 
            s2=aList[i].get_avg_compactnessR()
        if(s2>mx):
            mx=s2
    mn=float("inf")
        
    for i in range(len(aList)): #for every cluster
        for j in range(len(aList)): #for every cluster
            
            if(j>i):
                s2=0
                if(mode==0): #histogram method
                    b = aList[i].combine(aList[j])
                    s2 = b.calcCompactness()
                elif(mode==4): #ichino method
                    b = aList[i].combineR(aList[j])
                    s2=b.get_avg_compactnessR()

                if(s2<mn):
                    mn=s2
    print(str(len(aList))+" "+str(mn)+" "+str(mx), end=" ")
    if(mx>0):
        print(mn/mx)
    else:
        print()

def calcVRC(aList,nrf,BG,nre,TSSG):
#     WGSS=0

    if(len(aList)<2):
        return
    WGSS2=0
#     BGSS=0
    BGSS2=0
    for i in range(len(aList)):
#         s=0
        s2=0
#         G=[0 for x in range(nrf)]
        G2=[0 for x in range(nrf)]
#         for f in range(nrf):
#             G[f]=aList[i].getHistogram(f).get_quantile(0.5)
#             s+=(BG[f]-G[f])**2
        for j in range(len(aList[i].list)):
            for f in range(nrf):
                G2[f]+=aList[i].list[j].getHistogram(f).get_quantile(0.5)
#                 WGSS+=(G[f]-aList[i].list[j].getHistogram(f).get_quantile(0.5))**2
        for f in range(nrf):
            G2[f]/=len(aList[i].list)
            s2+=(BG[f]-G2[f])**2
        for j in range(len(aList[i].list)):
            for f in range(nrf):
                WGSS2+=(G2[f]-aList[i].list[j].getHistogram(f).get_quantile(0.5))**2
#         BGSS+=len(aList[i].list)*s
        BGSS2+=len(aList[i].list)*s2
#     CHI=((nre-len(aList))*BGSS)/((len(aList)-1)*WGSS)
    CHI2=((nre-len(aList))*BGSS2)/((len(aList)-1)*WGSS2)
    print(str(len(aList))+" "+str(CHI2))