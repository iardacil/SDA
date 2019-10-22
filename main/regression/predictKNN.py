#predicting values based by K closest neighbours


def findMostSimilarE(aList,i,k,l,nre,nrf,quantiles):
    nrq=len(quantiles)
    similar=[]
    dissM=[]
    
    
    for j in range(nre):
        if(i!=j):
            #calculate distance
            diss=0
            cnt=0
            for f in range(nrf):
                if l[f]==1:
                    for q in range(nrq):
                        diss+=abs(aList[i].getHistogram(f).get_quantile(quantiles[q])-aList[j].getHistogram(f).get_quantile(quantiles[q]))
                        cnt+=1
            diss=diss/cnt

            if(diss<=k):
                similar.append(j)
                dissM.append(diss)
    return similar,dissM

def findMostSimilar(aList,i,k,l,nre,nrf,quantiles):
    nrq=len(quantiles)
    similar=[-1 for x in range(k)]
    dissM=[float("inf") for x in range(k)]
    
    
    for j in range(nre):
        if(i!=j):
            #calculate distance
            diss=0
            cnt=0
            for f in range(nrf):
                if l[f]==1:
                    for q in range(nrq):
                        diss+=abs(aList[i].getHistogram(f).get_quantile(quantiles[q])-aList[j].getHistogram(f).get_quantile(quantiles[q]))
                        cnt+=1
            diss=diss/cnt
            for p in range(k):
                if(diss<dissM[p]):
                    pp=k-1
                    while(pp>p):
                        similar[pp]=similar[pp-1]
                        dissM[pp]=dissM[pp-1]
                        pp-=1
                    similar[p]=j
                    dissM[p]=diss
                    break
    return similar,dissM

def predictQ(aList,similar,f,q,quantiles,corrections=[],dissM=[]):
    if(len(corrections)==0):
        corrections=[0 for x in range(len(quantiles))]
    if(len(dissM)==0):
        dissM=[0 for x in range(len(similar))]
    value=[]
    for i in range(len(similar)):
#         value.append((aList[similar[i]].getHistogram(f).get_quantile(quantiles[q])+corrections[i])*dissM[i])
        value.append(aList[similar[i]].getHistogram(f).get_quantile(quantiles[q]))
    return value

def findCorrections(aList,i,similar,l,nrf,quantiles):
    corrections=[0 for x in range(len(similar))]
    for k in range(len(similar)):
        cnt=0
        for f in range(nrf):
            if l[f]==1:
                for q in range(len(quantiles)):
                    corrections[k]+=aList[i].getHistogram(f).get_quantile(quantiles[q])-aList[similar[k]].getHistogram(f).get_quantile(quantiles[q])
                    cnt+=1
        corrections[k]/=cnt
    return corrections

def predict(aList,i,quantiles,quantiles2,l,nre,nrf,f,k,mode=0):
    nrq=len(quantiles)

    if(k==-1):
        e=0.02
        similar,dissM=findMostSimilarE(aList,i,e,l,nre,nrf,quantiles2)
        k=len(similar)
#         print(similar, end=" ")
#         print(k, end=" ")
    elif(k==-2):
        e=1
        similar,dissM=findMostSimilarE(aList,i,e,l,nre,nrf,quantiles2)
        k=len(similar)
    else:
        similar,dissM=findMostSimilar(aList,i,k,l,nre,nrf,quantiles2)
#     sM=0
#     diffM=dissM[k-1]-dissM[0]
#     minM=dissM[0]
# 
#     for ii in range(k):
#         if(diffM!=0):
#             dissM[ii]=1-0.5*(dissM[ii]-minM)/diffM
#         else:
#             dissM[ii]=1


    if(mode==1):
        corrections=findCorrections(aList,i,similar,l,nrf,quantiles)
    else:
        corrections=[0 for x in range(k)]

    pv=[]
    t=-1
    if(k!=0):
        t=0
        for q in range(nrq):
            pv.append(predictQ(aList,similar,f,q,quantiles,corrections))


            t+=abs((sum(pv[q])/k)-(aList[i].getHistogram(f).get_quantile(quantiles[q])))
    return k,t
        
        
                        
                        
                        