'''
Class for Microscopib Hierarchical Conceptual Clustering
@author: Kadri Umbleja
'''

from main.Common import calcToC, calcBack, SymbolicObject

#microscopic similiarities
#r=0 returns diss and structure
#r1 returns the list
def MHCC(aList2,nre,nrf,objects,quantiles,till=1,r=0,point=-1,fs=[],CT=[],Fmin=[],Fmax=[],noteQ=[],qRmin=[],qRmax=[],clusters=[],CC=-1,minQ1=[],maxQ1=[]):
    nrq=len(quantiles)
    
    #feature selection, use only selected features or if not specificed, then all
    if(len(fs)==0):
        fs=[1 for y in range(nrf)]
    
    #produce specific notes on nodes for dendogram
    include_CT=0
    if(len(CT)>0 or len(qRmin)>0): #make 
        include_CT=1
    
    pp=[] #placeholder for quantile points
    kk=[] #placeholder for names
    nn=[] #element index holder for dissimilarity matrix
    pp2=[] #placeholder for original obects
    
    if(point==-1):
        for q in range(nrq):
            if(quantiles[q]==0.5):
                point=q
    
    #dissimilarity matrix
    diss=[ [ 0 for y in range( nre ) ] for x in range( nre ) ]
    #initialise placeholders
    for i in range(nre):
        kk.append(objects[i])
        nn.append([i])
    
    #min/max values for quantiles
    minQ=[[ float("inf") for y in range(nrq ) ] for y in range( nrf )]
    maxQ=[[ float("-inf") for y in range( nrq) ]  for y in range( nrf )]
    
    BG=[[0 for x in range(nrf)] for y in range(nrq)]
    
    #initialize quantile values and min/max values
    for i in range(nre):
        ll=[]
        for f in range(nrf):
            a=[]
            for j in range(nrq):
                val=aList2[i].getHistogram(f).get_quantile(quantiles[j])
                if(val<minQ[f][j]):
                    minQ[f][j]=val
                if(val>maxQ[f][j]):
                    maxQ[f][j]=val    
                b=[val]
                a.append([val])
#                 if(j==point):
                BG[j][f]+=val   
            ll.append(a)
        pp.append(ll)
        pp2.append(ll)
    
        
    if(len(minQ1)>0):
        minQ=minQ1
    if(len(maxQ1)>0):
        maxQ=maxQ1
        
    

    
    for j in range(nrq):
        for f in range(nrf):
            BG[j][f]/=nre

    
#     tss=[0 for x in range(nrf)]
#     TSSG=0
#     for i in range(len(pp)):
#         for f in range(nrf):
#             tss[f]+=(BG[f]-pp[i][f][point][0])**2
#     for f in range(nrf):
#         TSSG+=tss[f]
        
    notes=[]
    #clustering
    while (len(pp)>till):
        
#         calcDunn(pp,nrf,nrq,maxQ,minQ)
        


        mm=float("inf")
        ind1=-1
        ind2=-1
        
        #over all elements left
        for i in range(len(pp)):
            for j in range(len(pp)):
                if(i>j):
                    d=0
                    nr=0
                    #over all features
                    for f in range(nrf):
                        if(fs[f]==1):
                            #over all quantiles
                            for n in range(nrq):
                                if(minQ[f][n]!=maxQ[f][n]):                
                                    mn=min(pp[i][f][n]+pp[j][f][n])
                                    mx=max(pp[i][f][n]+pp[j][f][n])
                                    d+=(mx-mn)/(maxQ[f][n]-minQ[f][n])
                                    #for normalization later
                                    nr+=1
                    #check if best result          
                    if(d/(nr)<mm):
                        mm=d/(nr)
                        ind1=i
                        ind2=j
        #merge features
        for f in range(nrf):
            for n in range(len(quantiles)):
                pp[ind2][f][n]+=pp[ind1][f][n]
        #update dissimilarity matrix
        for i in range(len(nn[ind1])):
            for j in range(len(nn[ind2])):
                diss[nn[ind1][i]][nn[ind2][j]]=mm
                diss[nn[ind2][j]][nn[ind1][i]]=mm
        
        
        
        if(include_CT==1):
            if(len(CT)>0):
                rr=""
                if(len(noteQ)==0):
                    for i in range(len(CT)):
                        if(i>0):
                            rr+=", "
                        rr+=""+str(round(calcBack(min(pp[ind2][CT[i]][0]),Fmin,Fmax,CT[i]),3))+"<=F"+str(CT[i]+1)+"<="+str(round(calcBack(max(pp[ind2][CT[i]][nrq-1]),Fmin,Fmax,CT[i]),3))
                else:
                    for i in range(len(CT)):
                        if(i>0):
                            rr+=", "
                        rr+=""+str(round(calcBack(min(pp[ind2][CT[i]][noteQ[i]]),Fmin,Fmax,CT[i]),3))+"<=F"+str(CT[i]+1)+" "+str(round(quantiles[noteQ[i]]*100))+"%<="+str(round(calcBack(max(pp[ind2][CT[i]][noteQ[i]]),Fmin,Fmax,CT[i]),3))
                notes.append(rr)
            else: #the supervised case
                rr=""
                inc=-1
                print(len(clusters))
                print(len(nn))
                print(len(nn[ind1]))
                print(nn[ind1][i])
                if(clusters[nn[ind1][i]]!=CC):
                    for f in range(nrf):
                        if(inc!=-1 and inc>0):
                            rr+="\n"
                        inc=0
                        for n in range(len(quantiles)):
                            if(min(pp[ind2][f][n])>qRmax[f][n] or max(pp[ind2][f][n])<qRmin[f][n]):
                                if(inc==0):
                                    rr+="F"+str(f+1)+" "+str(round(quantiles[n]*100))+"%"
                                    inc+=1
                                else:
                                    rr+=""+str(round(quantiles[n]*100))+"%"
                notes.append(rr)
        
        
        #new name for concept
        kk[ind2]="("+str(kk[ind2])+"-"+str(kk[ind1])+")"

        print("merging "+str(kk[ind2])+" with "+str(mm))
        nn[ind2]+=nn[ind1]
        #remove not needed 
        pp.pop(ind1)
        kk.pop(ind1)
        nn.pop(ind1)
        
        
        
        
#         calcVRC(pp,pp2,nn,point,nrf,BG,nre,nrq,minQ,maxQ)
    
        
#         print(str(len(pp))+" "+str(CHI))
#     for i in range(len(pp)):
#         i=3
#         print(kk[i])
#         for j in range(nrf):
#             print("F"+str(j+1),end="&")
#             for q in range(nrq):
#              
#                 if(q<nrq-1):
# #                     print("["+str(calcToC(min(pp[i][j][q])))+"--"+str(calcToC(max(pp[i][j][q])))+"]",end="&")
#                     print("["+str(calcToC(min(pp[i][j][q])))+"]",end="&")
#                 else:
# #                     print("["+str(calcToC(min(pp[i][j][q])))+"--"+str(calcToC(max(pp[i][j][q])))+"]\\\\")
#                     print("["+str(calcToC(min(pp[i][j][q])))+"]\\\\")
#         print()
#         print()


    if(r==0):
        return diss,kk
    elif(r==2 or r==3):
        return diss,notes
    elif(r==1):
        ANew=[]
        
        for i in range(len(pp)):
            hist=""
            for f in range(nrf):
                hist+="{"
                for n in range(len(quantiles)-1):
                    hist+="["+str(min(pp[i][f][n]))+","+str(max(pp[i][f][n+1]))+"]"+str(quantiles[n+1]-quantiles[n])
                    if(n<len(quantiles)-2):
                        hist+=";"
                hist+="}"
#                 print(str(calcToC(min(pp[i][f][2])))+" "+str(calcToC(max(pp[i][f][4]))))
                if(f<(nrf-1)):
                    hist+=","            
             
            ANew.append(SymbolicObject(hist,kk[i],aList2[0].types,[0 for x in range(nrf)],[1 for x in range(nrf)],id=i))
        return ANew
    elif(r==11):
        return pp,kk



def calcDunn(pp,nrf,nrq,maxQ,minQ):
    mx=0
    for i in range(len(pp)): #for every cluster
        s2=0
        
        for f in range(nrf): #for every feature
#             s2+=max(pp[i][f][nrq-1])-min(pp[i][f][0])
#         if(s2/nrf>mx):
#             mx=s2/nrf
            for q in range(nrq):
                s2+=(max(pp[i][f][q])-min(pp[i][f][q]))/(maxQ[f][q]-minQ[f][q])   
        if(s2/(nrf*nrq)>mx):
            mx=s2/(nrf*nrq)
    mn=float("inf")
    for i in range(len(pp)): #for every cluster
        for j in range(len(pp)): #for every cluster
            s2=0
            if(j>i):
                for f in range(nrf): #for every feature
#                     s2+=max(pp[i][f][nrq-1]+pp[j][f][nrq-1])-min(pp[i][f][0]+pp[j][f][0])
#                 if(s2/nrf<mn):
#                     mn=s2/nrf
                    for q in range(nrq):
                        if(max(pp[i][f][q])<min(pp[j][f][q])):
                            s2+=(min(pp[j][f][q])-max(pp[i][f][q]))/(maxQ[f][q]-minQ[f][q])
                        elif(max(pp[j][f][q])<min(pp[i][f][q])):
                            s2+=(min(pp[i][f][q])-max(pp[j][f][q]))/(maxQ[f][q]-minQ[f][q])
                
#                         mn2=min(pp[i][f][q]+pp[j][f][q])
#                         mx2=max(pp[i][f][q]+pp[j][f][q])
#                         s2+=(mx2-mn2)/(maxQ[f][q]-minQ[f][q])         
                if(s2/(nrf*nrq)<mn):
                    mn=s2/(nrf*nrq)
    print(str(len(pp))+" "+str(mn)+" "+str(mx), end=" ")
    if(mx>0):
        print(mn/mx)
    else:
        print()
         

def calcVRC(pp,pp2,nn,point,nrf,BG,nre,nrq,minQ,maxQ):
    
    if(len(pp)==1):
        return

    CHI2=0
    WGSS2=0
    BGSS2=0
    WGSS=[0 for y in range(nrq)]
    BGSS=[0 for y in range(nrq)]
    CHI=[0 for y in range(nrq)]
    for i in range(len(pp)):
        s2=0
        G2=[[0 for x in range(nrf)] for y in range(nrq)]
        for j in range(len(nn[i])):
            for f in range(nrf):
                for q in range(nrq):
                    G2[q][f]+=pp2[nn[i][j]][f][q][0]

        s1=[0 for y in range(nrq)]
        for q in range(nrq):
            for f in range(nrf):
                G2[q][f]/=len(nn[i])
                s1[q]+=((BG[q][f]-G2[q][f]))**2
        
        for j in range(len(nn[i])):
            for q in range(nrq):
                for f in range(nrf):
                    WGSS[q]+=((G2[q][f]-pp2[nn[i][j]][f][q][0]))**2
#         BGSS+=len(aList[i].list)*s
        for q in range(nrq):
            BGSS[q]+=len(nn[i])*s1[q]
    
    for q in range(nrq):
        if(WGSS[q]!=0):
            CHI[q]=((nre-len(pp))*BGSS[q])/((len(pp)-1)*WGSS[q])
#     CHI=((nre-len(aList))*BGSS)/((len(aList)-1)*WGSS)
    
    CHI2=max(CHI)
    ind=-1
    for q in range(nrq):
        if(CHI[q]==CHI2):
            ind=q
#     CHI2=sum(CHI)/len(CHI)
    print(str(len(pp))+" "+str(CHI2)+" "+str(ind))
#     print(str(len(pp))+" "+str(BGSS)+" "+str(WGSS)+" "+str(CHI)+" "+str(CHI2))
    

# def calcVRC(pp,pp2,nn,point,nrf,BG,nre,TSSG):
#     WGSS=0
#     BGSS=0
#     
#     WGSS2=0
#     BGSS2=0
#     for i in range(len(pp)):
#         s2=0
#         G2=[0 for x in range(nrf)]
#         for j in range(len(nn[i])):
#             for f in range(nrf):
#                 G2[f]+=pp2[nn[i][j]][f][point][0]
# #                 WGSS+=(G[f]-aList[i].list[j].getHistogram(f).get_quantile(0.5))**2
#         for f in range(nrf):
#             G2[f]/=len(nn[i])
#             s2+=(BG[f]-G2[f])**2
#         for j in range(len(nn[i])):
#             for f in range(nrf):
#                 WGSS2+=(G2[f]-pp2[nn[i][j]][f][point][0])**2
# #         BGSS+=len(aList[i].list)*s
#         BGSS2+=len(nn[i])*s2
# #     CHI=((nre-len(aList))*BGSS)/((len(aList)-1)*WGSS)
#     CHI2=((nre-len(pp))*BGSS2)/((len(pp)-1)*WGSS2)
#     print(str(len(pp))+" "+str(BGSS2)+" "+str(WGSS2)+" "+str(CHI2))
#     

#     for i in range(len(pp)): #for every cluster
#         s2=0
#         s3=0
#         LG=[0 for x in range(nrf)]
#         
#         for f in range(nrf): #for every feature
#             for j in range(len(nn[i])):
#                 LG[f]+=pp2[nn[i][j]][f][point][0]
#             LG[f]/=len(nn[i])
#         
# 
#         for f in range(nrf): #for every featur              
#             
#             for q in range(nrq):
#                 s2+=(max(pp[i][f][q])-min(pp[i][f][q]))  
#                 s3+=abs(LG[f]-GBC[f])
#         WGSS+=(s2/(len(nn[i])*nrf*nrq))
#         BGSS+=((s3*len(nn[i]))/(nrf*nrq))
#         print(str(kk[i])+" "+str(s2/(len(nn[i])*nrf*nrq))+" "+str((s3*len(nn[i]))/(nrf*nrq)))
#         CHI=(BGSS*(nre-len(pp)))/(WGSS*(len(pp)-1))
    
    