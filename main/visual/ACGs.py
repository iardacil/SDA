'''
ACGs are drawn in this cass

'''
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from main.Common import Histogram,SymbolicObject,normalize
import matplotlib.patches as mpatches

width=1

#help functions for plotting
def addPoint(i,v1,ax,colors,dot_size,j=None):
    if j is None:
        j=i
    v=[v1]
    y=[(i)*width]
    ax.plot(y, v, 'ro',color=colors[j],markersize=dot_size)
    
    
def addLineByPoints(v,y,ax,color,line_width,dot_size):
    ax.plot(y, v, color=color,linewidth=line_width,markersize=dot_size)

def addLineByPoints2(v,y,ax,color,line_width,dot_size):
    ax.plot(y, v,  'ro-',color=color,linewidth=line_width,markersize=dot_size)
    
def addLineUniliteral(j,i,v1,v2,ax,colors,line_width,dot_size):
    v=[v1,v2]
    y=[(i)*width,(i+1)*width]
    ax.plot(y, v, 'ro-',color=colors[j],linewidth=line_width,markersize=dot_size)

def transformSD(aList,nre,nrf,quantiles=[0,1]):
    nrq=len(quantiles)
    data=[ [ [ 0 for z in range( nrq ) ] for y in range( nrf ) ] for x in range( nre ) ]
       
        #parse data to Nxdxn
    for i in range(nre):
        for j in range(nrf):
            for k in range (nrq):
                data[i][j][k]=aList[i].getHistogram(j).get_quantile(quantiles[k])
    return data

#MAIN VISUALIZATION FOR SINGLE POINT DATA          
def drawSimple_by_points(data,objects,line_width,dot_size,l=[]):
    

    fig, axes = plt.subplots(1, 1, sharey=False,figsize=(16, 8)) 
    nrf=len(data[0])
    nre=len(data)
    #diffrence between objects
    m=nrf+3 
    
    #color shades are produced
    colors=[0 for x in range(nrf)]
    cmap2 = plt.cm.get_cmap('hsv')
    start=0
    end=1
    for i in range(nrf):
        colors[i] = cmap2(start+((end-start)/(nrf))*(i+1))
    
    if(len(l)==0):
        l=[x for x in range(nre)]
    #lines are drawn    
    for ll in range(nre):
        i=l[ll]
        so_far=0;
        for j in range(nrf):
                if(j==0):
                    so_far=data[i][0] 
                else:
                    addLineUniliteral(j,(i)*m+j,so_far,so_far+(data[i][j] ),axes,colors,line_width,dot_size)
                    so_far=so_far+(data[i][j])
        
        #add name    
        axes.text((i+1)*m-(round(m/2)), so_far+0.2, objects[i],{'ha': 'left', 'va': 'bottom'}, rotation=90)
    
    

    axes.set_xlim([0,nre*m])
    axes.set_ylim([0, nrf])
    
    plt.savefig('graphs/SACG.jpg', format='jpg', dpi=300)
    plt.show()
    plt.close()
    
 
    
    
    
#METHODS FOR THREE WAY DATA
def drawSimple_for_quantiles(aList,objects,line_width,dot_size,quantiles,l=[]):
    
    
    fig, axes = plt.subplots(1, 1, sharey=False,figsize=(16, 8)) 
    
    nrf=aList[0].nrf
    nre=len(aList)
    nrq=len(quantiles)
    data=transformSD(aList,nre,nrf,quantiles)
    
    #diffrence between objects
    m=nrf+3 
    
    #color shades are produced
    colors=[0 for x in range(nrf)]
    cmap2 = plt.cm.get_cmap('hsv')
    start=0
    end=1
    for i in range(nrf):
        colors[i] = cmap2(start+((end-start)/(nrf))*(i+1))
    
    #for adequate limits
    max=0
    if(len(l)==0):
        l=[x for x in range(nre)]
    #lines are drawn    
    for ll in range(nre):
        i=l[ll]
        so_far=0;
        for j in range(nrf):
            if(i>=0):
                if(j==0):
                    so_far=data[i][j][0] 
                else:
                    addLineUniliteral(j,(i)*m+j,so_far,so_far+(data[i][j][nrq-1] -data[i][j][0] ),axes,colors,line_width,dot_size)
                    so_far=so_far+(data[i][j][nrq-1] -data[i][j][0] )
        if(so_far>max):
            max=so_far
        
        #add name    
        axes.text((i+1)*m-(round(nrf/2)), so_far+0.1, objects[i],{'ha': 'left', 'va': 'bottom'}, rotation=90)
    
    

    axes.set_xlim([0,nre*m])
    axes.set_ylim([0, max*1.1])
    
    plt.savefig('graphs/SACG.jpg', format='jpg', dpi=300)
    plt.show()
    plt.close()
    
def drawMinMax_for_quantiles(aList,objects,line_width,dot_size,quantiles,l=[]):
  
    fig, axes = plt.subplots(1, 1, sharey=False,figsize=(16, 8)) 
    nrf=aList[0].nrf
    nre=len(aList)
    nrq=len(quantiles)
    data=transformSD(aList,nre,nrf,quantiles)
    #diffrence between objects
    m1=2
    m2=4
    m=nrf+m1+nrf+m2
    
    #color shades are produced, one for min, other to max
    colors=[0 for x in range(nrf)]
    cmap = plt.cm.get_cmap('Blues')
    start=0.5
    end=1
    for i in range(nrf):
        colors[i] = cmap(start+((end-start)/(nrf+1))*(i+1))
        
    colors2=[0 for x in range(nrf)]
    cmap = plt.cm.get_cmap('Reds')
    start=0.5
    end=1
    for i in range(nrf):
        colors2[i] = cmap(start+((end-start)/(nrf+1))*(i+1))
    
    if(len(l)==0):
        l=[x for x in range(nre)]
    #lines are drawn    
    for ll in range(nre):
        i=l[ll]
        so_far=0;
        for j in range(nrf):
            if(j==0):
                so_far=data[i][j][0]
            else:
                addLineUniliteral(j,(i)*m+j,so_far,so_far+(data[i][j][0]),axes,colors,line_width,dot_size)
                so_far=so_far+(data[i][j][0])
                
        so_far=0;
        for j in range(nrf):
            if(j==0):
                so_far=data[i][j][nrq-1]  
            else:
                addLineUniliteral(j,(i)*m+nrf+m1+j,so_far,so_far+(data[i][j][nrq-1] ),axes,colors2,line_width,dot_size)
                so_far=so_far+(data[i][j][nrq-1] )
                    
        #add name    
        axes.text((i+1)*m-(nrf), so_far+0.1, objects[i],{'ha': 'left', 'va': 'bottom'}, rotation=90)


    axes.set_xlim([0,nre*m])
    axes.set_ylim([0, nrf])
    
        
        
    
    plt.savefig('graphs/MinMaxACG.jpg', format='jpg', dpi=300)
    plt.show()
    plt.close()
    
    
def drawUniliteral_for_quantiles(aList,objects,line_width,dot_size,quantiles,l=[]):
    fig, axes = plt.subplots(1, 1, sharey=False,figsize=(16, 8))
    nrf=aList[0].nrf
    nre=len(aList)
    nrq=len(quantiles)
    data=transformSD(aList,nre,nrf,quantiles)
    #diffrence between objects
    m1=3
    m=2*nrf+m1
     
     #color shades are produced
    colors=[0 for x in range(nrf)]
    cmap2 = plt.cm.get_cmap('hsv')
    start=0
    end=1
    for i in range(nrf):
        colors[i] = cmap2(start+((end-start)/(nrf))*(i+1))
    
    
    #for adequate limits
    max=0 
    

    if(len(l)==0):
        l=[x for x in range(nre)]
    #lines are drawn    
    for ll in range(nre):
        i=l[ll]
        so_far=0;
        for j in range(nrf):
            if(j==0):
                so_far=data[i][j][0]
                addLineUniliteral(j,(i*m)+j+1,so_far,data[i][j][0],axes,colors,line_width,dot_size)
                so_far=data[i][j][0]
            else:
                addLineUniliteral(j,(i*m)+j*2,so_far,so_far+data[i][j][0],axes,colors,line_width,dot_size)
                addPoint((i*m)+j*2,so_far,axes,colors,dot_size,j-1)
                so_far=so_far+data[i][j][0]
                addLineUniliteral(j,(i*m)+j*2+1,so_far,so_far+data[i][j][nrq-1]-data[i][j][0],axes,colors,line_width,dot_size)
                so_far=so_far+data[i][j][nrq-1]-data[i][j][0]
        if(so_far>max):
            max=so_far
 
         
        axes.text((i+1)*m-(round(nrf/2)), so_far+0.1, objects[i],{'ha': 'left', 'va': 'bottom'}, rotation=90)
     
 
    axes.set_xlim([0,nre*m])
    axes.set_ylim([0, max*1.1])
    
    
    plt.savefig('graphs/uniliteral.jpg', format='jpg', dpi=300)
    plt.show()
    plt.close()
    
def drawBiliteral_for_quantiles(aList,objects,line_width,dot_size,quantiles,l=[]):
    fig, axes = plt.subplots(1, 1, sharey=False,figsize=(16, 8))

    nrf=aList[0].nrf
    nre=len(aList)
    nrq=len(quantiles)
    data=transformSD(aList,nre,nrf,quantiles)
    #diffrence between objects
    m1=3
    m=2*nrf+m1
     
     #color shades are produced
    colors=[0 for x in range(nrf)]
    cmap2 = plt.cm.get_cmap('hsv')
    start=0
    end=1
    for i in range(nrf):
        colors[i] = cmap2(start+((end-start)/(nrf))*(i+1))
    
    
    #for adequate limits
    max=0 
    #lines are drawn
    if(len(l)==0):
        l=[x for x in range(nre)]
    #lines are drawn    
    for ll in range(nre):
        i=l[ll]
        so_far=0;
        for j in range(nrf):
            if(j==0):
                so_far=data[i][j][0]
                addLineUniliteral(j,(i*m)+j+1,so_far,data[i][j][nrq-1],axes,colors,line_width,dot_size)
                so_far+=data[i][j][nrq-1]
            else:
                addLineUniliteral(j,(i*m)+j*2,so_far,so_far+data[i][j][0],axes,colors,line_width,dot_size)
                addPoint((i*m)+j*2,so_far,axes,colors,dot_size,j-1)
                so_far=so_far+data[i][j][0]
                addLineUniliteral(j,(i*m)+j*2+1,so_far,so_far+data[i][j][nrq-1],axes,colors,line_width,dot_size)
                so_far=so_far+data[i][j][nrq-1]
        if(so_far>max):
            max=so_far
 
         
        axes.text((i+1)*m-(round(nrf/2)), so_far+0.1, objects[i],{'ha': 'left', 'va': 'bottom'}, rotation=90)
     
 
    axes.set_xlim([0,nre*m])
    axes.set_ylim([0, max*1.1])
    
    
    plt.savefig('graphs/biliteral.jpg', format='jpg', dpi=300)
    plt.show()
    plt.close()
    

#l is custom order
def drawQVACG_for_quantiles(aList,objects,line_width,dot_size,quantiles,l=[]):

    fig, axes = plt.subplots(1, 1, sharey=False,figsize=(16, 8)) 
    nrf=aList[0].nrf
    nre=len(aList)
    nrq=len(quantiles)
    data=transformSD(aList,nre,nrf,quantiles)
    
    #diffrence between features
    m1=1
    #diffrence between objects
    m2=3
    m=(m1+nrq)*nrf+m2
    
    #for adequate limits
    max=0 
     #color shades are produced
    colors=[0 for x in range(nrf)]
    cmap2 = plt.cm.get_cmap('hsv')
    start=0
    end=1
    for i in range(nrf):
        colors[i] = cmap2(start+((end-start)/(nrf))*(i+1))
    
    
    if(len(l)==0):
        l=[x for x in range(nre)]
    #lines are drawn    
    for ll in range(nre):
        i=l[ll]
          
        for k in range(nrq):
            so_far2=[0 for x in range(nrq)]
            for j in range(nrf-1):
                if(j==0):
                    v1=so_far2[k]+data[i][j][k]
                    so_far2[k]=v1
                    v2=so_far2[k]+data[i][j+1][k]
                    so_far2[k]=v2
                else:
                    v1=so_far2[k]
                    v2=so_far2[k]+data[i][j+1][k]
                    so_far2[k]=v2
                addLineUniliteral(j,(i*m)+nrf*k+j,v1,v2,axes,colors,line_width,dot_size)
            addPoint((i*m)+nrf*(k+1),so_far2[k],axes,colors,dot_size,nrf-1)
        
            if(so_far2[nrq-1]>max):
                max=so_far2[nrq-1]
               

        axes.text(i*m+round(m/2), so_far2[nrq-1]+0.2, objects[i],{'ha': 'left', 'va': 'bottom'}, rotation=90)
    
    axes.set_xlim([0,nre*m])
    axes.set_ylim([0, max*1.1])
    
    
    plt.savefig('graphs/QV_ACG.jpg', format='jpg', dpi=300)
    plt.show()
    plt.close()

def drawSQVACG_for_quantiles(aList,objects,line_width,dot_size,quantiles,l=[]):
    
    fig, axes = plt.subplots(1, 1, sharey=False,figsize=(16, 8)) 
    nrf=aList[0].nrf
    nre=len(aList)
    nrq=len(quantiles)
    data=transformSD(aList,nre,nrf,quantiles)
    
    
    m=(nrf+1)*nrq

    
    colors=[0 for x in range(nrf)]
    colors2=[0 for x in range(nrq)]
    cmap = plt.cm.get_cmap('Greys')
    cmap2 = plt.cm.get_cmap('hsv')
    start=0.3
    end=1
    for i in range(nrf):
        colors[i] = cmap2(start+((end-start)/(nrf))*(i+1))
    for i in range(nrq):
        colors2[i] = cmap2((1/(nrq+1))*(i+1))
    
    
    pp=[[0 for y in range(nre*(nrq-1))] for x in range(nrf)]


    gmx=0
    if(len(l)==0):
        l=[x for x in range(nre)]
    #lines are drawn    
    for ll in range(nre):
        e=l[ll]
        so_far=0
        for j in range(nrf):
            sf2=0
            for k in range(nrq-1):
                v1=so_far
                v2=so_far+(data[e][j][k+1]-data[e][j][k])
                addLineUniliteral(j,(ll)*m+(nrq-1)*j+k,v1,v2,axes,colors,line_width,dot_size)
                so_far=v2
                sf2+=(data[e][j][k+1]-data[e][j][k])
                pp[j][(nrq-1)*e+k]=so_far
            mx=so_far 
        gmx=max(gmx,so_far) 
        
        axes.text((ll+1)*m-(m/2), mx+0.1, objects[e],{'ha': 'left', 'va': 'bottom'}, rotation=90)  
        e=e+1
    
    axes.set_xlim([0,m*(nre)+6])
    axes.set_ylim([0, gmx+2])
    
    plt.savefig('graphs/SQVACG.png', format='png', dpi=300)
    plt.show()
    plt.close()
    return pp

def drawTACG_for_quantiles(aList,objects,line_width,dot_size,quantiles,l=[]):
   
    fig, axes = plt.subplots(1, 1, sharey=False,figsize=(16, 8)) 
    nrf=aList[0].nrf
    nre=len(aList)
    nrq=len(quantiles)
    data=transformSD(aList,nre,nrf,quantiles)

    m=(nrf+1)*nrq

    
    colors=[0 for x in range(nrf)]
    colors2=[0 for x in range(nrq)]
    cmap = plt.cm.get_cmap('Greys')
    cmap2 = plt.cm.get_cmap('hsv')
    start=0
    end=1
    for i in range(nrf):
        colors[i] = cmap2(start+((end-start)/(nrf))*(i+1))
    for i in range(nrq):
        colors2[i] = cmap2((1/(nrq+1))*(i+1))
    
    gmx=0
    c=0

    if(len(l)==0):
        l=[x for x in range(nre)]
    for ll in range(nre):
        e=l[ll]
        
#         axes.set_title(objects[e],rotation='vertical')
        mx=0

        first=0
        so_far2=0
        for j in range(nrf):
            for i in range(nrq):
                if(i==0):
                    if(first!=0):
                        v1=so_far2
                        v2=so_far2+data[e][j][i] 
                        addLineUniliteral(j,(ll)*m+nrq*j+i,v1,v2,axes,colors,line_width,dot_size)

                else:
                    v1=so_far2+data[e][j][i-1] 
                    v2=so_far2+data[e][j][i] 
                    addLineUniliteral(j,(ll)*m+nrq*j+i,v1,v2,axes,colors,line_width,dot_size)
            so_far2=v2
            c+=1
            first+=1

        mx=so_far2  

        gmx=max(gmx,so_far2) 

        axes.text((ll+1)*m-(m/2), mx+0.1, objects[e],{'ha': 'left', 'va': 'bottom'}, rotation=90)  
        e=e+1
    
    axes.set_xlim([0,m*(nre)+6])
    axes.set_ylim([0, gmx+2])
    
    plt.savefig('graphs/TACG.png', format='png', dpi=300)
    plt.show()
    plt.close()
    
def drawACC_for_quantiles(aList,objects,line_width,dot_size,quantiles,l=[]):

    fig, axes = plt.subplots(1, 1, sharey=False,figsize=(16, 8)) 
    nrf=aList[0].nrf
    nre=len(aList)
    nrq=len(quantiles)
    data=transformSD(aList,nre,nrf,quantiles)

    m=nrf*nrq+5
        
    colors=[0 for x in range(nrq)]
    colors2=[0 for x in range(nrf)]
    cmap = plt.cm.get_cmap('Greys')
    cmap2 = plt.cm.get_cmap('hsv')
    start=0.3
    end=1
    for i in range(nrq):
        colors[i] = cmap(start+((end-start)/(nrq+1))*(i+1))
    for i in range(nrf):
        colors2[i] = cmap2((1/(nrf+1))*(i+1))
        
    points=[[[0 for x in range(nrq)] for y in range(nrf)] for z in range(nre)]
    if(len(l)==0):
        l=[x for x in range(nre)]
    for ll in range(nre):
        e=l[ll]
        
        sf=[[0 for y in range(nrf-1)] for x in range(nrq)]
        mx=0
        for j in range(nrq):
            so_far2=[0 for x in range(nrq)]
            for i in range(nrf-1):
                if(i==0):
                    v1=so_far2[j]+data[e][i][j] 
                    so_far2[j]=v1
                    points[e][0][j]=v1
                    if(j>0):
                        addLineByPoints([points[e][0][j-1],points[e][0][j]],[(e)*m+nrf*(j-1)+i,(e)*m+nrf*j+i],axes,colors2[i],line_width,dot_size) 
                    v2=so_far2[j]+data[e][i+1][j] 
                    sf[j][i]=v2
                    so_far2[j]=v2
                else:
                    v1=so_far2[j]
                    v2=so_far2[j]+data[e][i+1][j] 
                    so_far2[j]=v2
                    sf[j][i]=v2
                addLineUniliteral(j,(e)*m+nrf*j+i,v1,v2,axes,colors,line_width,dot_size)
                
                
                if(j>0):
                    addLineByPoints([sf[j-1][i],sf[j][i]],[(e)*m+nrf*(j-1)+i+1,(e)*m+nrf*j+i+1],axes,colors2[i+1],line_width,dot_size)             
            
#             axes.text((j)*m+e*nrq, so_far2[e]+0.1, objects[e],{'ha': 'left', 'va': 'bottom'}, rotation=90) 
            mx=max(mx,so_far2[j])

            points[e][i+1][j]=sf[j][i]
            addPoint((e)*m+nrf*(j)+(nrf-1),so_far2[j],axes,colors,dot_size,j)
        
        axes.text((e+1)*m-(m/2), mx+0.1, objects[e],{'ha': 'left', 'va': 'bottom'}, rotation=90)  
        e=e+1
        
    
    axes.set_xlim([0,m*(nre+1)])
    axes.set_ylim([0, 10])
    
    plt.savefig('graphs/ACC.png', format='png', dpi=300)
    plt.show()
    plt.close()

def drawACC_FF_for_quantiles(aList,titles,line_width,dot_size,quantiles,l=[]):

    fig, axes = plt.subplots(1, 1, sharey=False,figsize=(16, 8)) 
    nrf=len(titles)
    nre=len(aList)
    nrq=len(quantiles)
    data=transformSD(aList,nre,nrf,quantiles)

    m=(nre+1)*nrq
    
    colors=[0 for x in range(nre)]
    colors2=[0 for x in range(nrq)]
    cmap = plt.cm.get_cmap('Greys')
    cmap2 = plt.cm.get_cmap('hsv')
    start=0.3
    end=1
    for i in range(nre):
        colors[i] = cmap(start+((end-start)/(nre+1))*(i+1))
    for i in range(nrq):
        colors2[i] = cmap2((1/(nrq+1))*(i+1))
  
   
    if(len(l)==0):
        l=[x for x in range(nre)]
    #lines are drawn    
    for j in range(nrf):
#         axes.set_title(objects[e],rotation='vertical')
        mx=0

        sf=[[0 for y in range(nrq-1)] for x in range(nre)]
        for ll in range(nre):
            e=l[ll]
            so_far2=[0 for x in range(nre)]
            
            for i in range(nrq-1):
                
                if(i==0):
                    v1=so_far2[e]+data[e][j][i]
                    so_far2[e]=v1
                    v2=so_far2[e]+data[e][j][i+1]
                    so_far2[e]=v2
                    sf[e][i]=v2
                else:
                    v1=so_far2[e]
                    v2=so_far2[e]+data[e][j][i+1]
                    so_far2[e]=v2
                    sf[e][i]=v2
                addLineUniliteral(e,(j)*m+nrq*e+i,v1,v2,axes,colors,line_width,dot_size)
                if(e>0):
                    addLineByPoints([sf[e-1][i],sf[e][i]],[(j)*m+nrq*(e-1)+i,(j)*m+nrq*e+i],axes,colors2[i],line_width,dot_size)             
            
#             axes.text((j)*m+e*nrq, so_far2[e]+0.1, objects[e],{'ha': 'left', 'va': 'bottom'}, rotation=90) 
            mx=max(mx,so_far2[e])
            addPoint((j)*m+nrq*(e)+(nrq-1),so_far2[j],axes,colors,dot_size,j)
            
        
        axes.text((j+1)*m-(m/2), mx+0.1, titles[j],{'ha': 'left', 'va': 'bottom'}, rotation=90)  

    
    axes.set_xlim([0,m*(nrf)+6])
    axes.set_ylim([0, nrq])
    
    plt.savefig('graphs/ACC_FF.png', format='png', dpi=300)
    plt.show()
    plt.close()
    
    
def drawFingertips_for_quantiles(aList,objects,line_width,dot_size,quantiles,l=[]):
    fig, axes = plt.subplots(1, 1, sharey=False,figsize=(16, 8)) 
    nrf=aList[0].nrf
    nre=len(aList)
    nrq=len(quantiles)
    data=transformSD(aList,nre,nrf,quantiles)
    m=(nrf+1)*nrq

    
    colors=[0 for x in range(nrf)]
    colors2=[0 for x in range(nrq)]
    cmap = plt.cm.get_cmap('Greys')
    cmap2 = plt.cm.get_cmap('hsv')
    start=0.3
    end=1
    for i in range(nrf):
        colors[i] = cmap(start+((end-start)/(nrf))*(i+1))
    for i in range(nrq):
        colors2[i] = cmap2((1/(nrq))*(i+1))
    
    
    points=[[[0 for x in range(nrq)] for y in range(nrf)] for z in range(nre)]
    

    if(len(l)==0):
        l=[x for x in range(nre)]
    for ll in range(nre):
        e=l[ll]
#         axes.set_title(objects[e],rotation='vertical')
        mx=0

        sf=[[0 for y in range(nrq-1)] for x in range(nrf)]
        for j in range(nrf):
            so_far2=[0 for x in range(nrf)]
            for i in range(nrq-1):
                if(i==0):
                    v1=so_far2[j]+data[e][j][i] 
                    so_far2[j]=v1
                    points[e][j][0]=v1
                    if(j>0):
                        addLineByPoints([points[e][j-1][0],points[e][j][0]],[(e)*m+nrq*(j-1)+i,(e)*m+nrq*j+i],axes,colors2[i],line_width,dot_size) 
                    v2=so_far2[j]+data[e][j][i+1]
                    so_far2[j]=v2
                    sf[j][i]=v2
                else:
                    v1=so_far2[j]
                    v2=so_far2[j]+data[e][j][i+1]
                    so_far2[j]=v2
                    sf[j][i]=v2
                addLineUniliteral(j,(e)*m+nrq*j+i,v1,v2,axes,colors,line_width,dot_size)
                if(j>0):
                    addLineByPoints([sf[j-1][i],sf[j][i]],[(e)*m+nrq*(j-1)+i+1,(e)*m+nrq*j+i+1],axes,colors2[i+1],line_width,dot_size)
                points[e][j][i+1]=sf[j][i]
            mx=max(mx,so_far2[j])
            
            addPoint((e)*m+nrq*(j)+(nrq-1),so_far2[j],axes,colors,dot_size,j)
            
        
        axes.text((e+1)*m-(m/2), mx+0.1, objects[e],{'ha': 'left', 'va': 'bottom'}, rotation=90)  
        e=e+1

        
    
    axes.set_xlim([0,m*(nre)+6])
    axes.set_ylim([0, nrq])
    
    plt.savefig('graphs/fingertips.png', format='png', dpi=300)
    plt.show()
    plt.close()
    

#plots quantiles in 2D plot
#f1 and f2 are features plotted
#L_x,L_y are coordinates for legend
#l may contain custom list of specific objects to be plotted - if empty all are plotted 
def drawQuantiles(aList,titles,objects,quantiles,line_width,dot_size,f1=0,f2=1,L_x=0.95,L_y=1,l=[]):
    if(len(l)==0):
        l=[x for x in range(len(aList))]
    
    fig, axes = plt.subplots(1, 1, sharey=False,figsize=(20, 20)) 

    nre=len(objects)
    nrl=len(l)
    nrq=len(quantiles)
    
    colors=[0 for x in range(nre)]
    cmap = plt.cm.get_cmap('Spectral')
    start=0
    end=1
    for i in range(nre):
        colors[i] = cmap(start+((end-start)/(nre+1))*(i+1))

    for ll in range(nrl):
        i=l[ll]
        for j in range(nrq-1):
            x1=aList[i].getHistogram(f1).get_quantile(quantiles[j])
            x2=aList[i].getHistogram(f1).get_quantile(quantiles[j+1])
            y1=aList[i].getHistogram(f2).get_quantile(quantiles[j])
            y2=aList[i].getHistogram(f2).get_quantile(quantiles[j+1])
            v=[x1,x2]
            y=[y1,y2]
            addLineByPoints2(v,y,axes,colors[i],line_width,dot_size)
    axes.set_xlim([0, 1])
    axes.set_ylim([0, 1])
    axes.set_xlabel(titles[f1])
    axes.set_ylabel(titles[f2])
    

 
    handles=[]   
    for ll in range(nrl):
        e=l[ll]
        color = colors[e]
        handles.append(mpatches.Patch(color=color, label=objects[e]))
     
    chartBox = axes.get_position()
    axes.set_position([chartBox.x0, chartBox.y0, chartBox.width*1, chartBox.height*1])
    
    axes.legend(handles=handles, loc='upper left', bbox_to_anchor=(1,1))
    axes.legend(handles=handles,loc='upper center', bbox_to_anchor= (L_x, L_y), ncol=1, 
            borderaxespad=0, frameon=False)
    plt.tight_layout(pad=7)
    
    plt.savefig('graphs/quantiles.png', format='png', dpi=300)
    plt.show()
    plt.close()