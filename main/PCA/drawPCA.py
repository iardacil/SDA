import matplotlib.patches as patches
import matplotlib.lines as lines    
import matplotlib.pyplot as plt

import numpy as np


#just one graph
def plot_resultsMAIN(pca,newX,n,nrf,ncomp=2,npoint=0):
    conf=0.5
    fig, ax = plt.subplots(figsize=(16*conf, 10*conf)) 
    n_samples=len(n)
    
    colors=[0 for x in range(n_samples)]
    cmap = plt.cm.get_cmap('Spectral')
    start=0
    end=1
    for i in range(n_samples):
        colors[i] = cmap(start+((end-start)/((n_samples)+1))*(i+1))
        
    cnt=int(len(newX)/n_samples)
    for i in range(0,n_samples):
        i1=i*cnt
        i2=(i+1)*cnt
        ax.scatter(newX[i1:i2,0],newX[i1:i2,1],color=colors[i])


    
    if(npoint==0):
        npoint= 2**nrf
    
    for i in range(0,len(n)):
        i1=i*npoint
        i2=(i+1)*npoint
        P1=newX[:,0]
        P2=newX[:,1]
        min1=min(P1[i1:i2])
        max1=max(P1[i1:i2])
        min2=min(P2[i1:i2])
        max2=max(P2[i1:i2])
        ax.add_patch(
        patches.Rectangle(
            (min1, min2),
            (max1-min1),
            (max2-min2),
            color=colors[i],
            fill=False      # remove background
            )
        )

    for i, txt in enumerate(n):
        ax.annotate(txt, (newX[i*npoint,0],newX[i*npoint,1]))   
        
    pca.explained_variance_ratio_[0]
    ax.set_xlabel(round(pca.explained_variance_ratio_[0]*100,2))
    ax.set_ylabel(round(pca.explained_variance_ratio_[1]*100,2))
    
    ax.set_title('First {0} Principal Vectors'.format(ncomp))
    ax.set_ylim([-6,6])
    ax.set_xlim([-6,6])
    
    plt.savefig('graphs/PCA_oneplot.png', format='png', dpi=300)
    plt.show()
    plt.close()

#draws 4 plots
def plot_results(X_scaled, pca,newX,n_samples,n_features,n,titles,backRes=[],Y=[],ncomp=2,imp1=0,imp2=1,config=[0,0,0],pperelem=0):
    conf=2   
    fig, ax = plt.subplots(2, 2, figsize=(16*conf, 16*conf))
     
    colors=[0 for x in range(n_samples)]
    cmap = plt.cm.get_cmap('Spectral')
    start=0
    end=1
    for i in range(n_samples):
        colors[i] = cmap(start+((end-start)/((n_samples)+1))*(i+1))
 
 
    if(pperelem==0):
        pperelem=1
    npoint= 2**n_features
 
     
##original
 
    minX=0
    minY=0
    maxX=0
    maxY=0
    
    cnt=int(len(X_scaled)/n_samples)
    for i in range(0,n_samples):
        i1=i*cnt
        i2=(i+1)*cnt
        ax[0,0].scatter(X_scaled[i1:i2,imp1],X_scaled[i1:i2,imp2],color=colors[i])
#     ax[0, 0].scatter(X_scaled[:,imp1],X_scaled[:,imp2])
    
 
    if(config[0]==1):
 
        for i in range(0,n_samples):
            i1=i*npoint
            i2=(i+1)*npoint
            P1=X_scaled[:,imp1]
            P2=X_scaled[:,imp2]
            min1=min(P1[i1:i2])
            max1=max(P1[i1:i2])
            min2=min(P2[i1:i2])
            max2=max(P2[i1:i2])
            ax[0, 0].add_patch(
            patches.Rectangle(
                (min1, min2),
                (max1-min1),
                (max2-min2),
                color=colors[i],
                fill=False      # remove background
                )
            )
            ax[0, 0].annotate(n[i], (max1, max2))
            minX=min([minX,min1,max1])
            maxX=max([maxX,min1,max1])
            minY=min([minY,min2,max2])
            maxY=max([maxY,min2,max2])
             
            l1 = lines.Line2D([min1,max1],[min2,max2], lw=1, color=colors[i], axes=ax[0, 0])
            ax[0, 0].add_line(l1)
        ax[0, 0].set_xlim([minX-abs(minX*0.1),maxX+abs(maxX*0.1)])
        ax[0, 0].set_ylim([minY-abs(minY*0.1),maxY+abs(maxY*0.1)])
             
    else:
        for i, txt in enumerate(n):
            ax[0, 0].annotate(txt, (X_scaled[i*pperelem,imp1],X_scaled[i*pperelem,imp2]))
  
#PCs       
#     ax[0, 1].scatter(newX[:,0],newX[:,1])
    cnt=int(len(newX)/n_samples)
    for i in range(0,n_samples):
        i1=i*cnt
        i2=(i+1)*cnt
        ax[0,1].scatter(newX[i1:i2,0],newX[i1:i2,1],color=colors[i])
     
    if(config[1]==1):
        minX=0
        minY=0
        maxX=0
        maxY=0
        for i in range(0,n_samples):
            i1=i*npoint
            i2=(i+1)*npoint
            P1=backRes[:,0]
            P2=backRes[:,1]
            min1=min(P1[i1:i2])
            max1=max(P1[i1:i2])
            min2=min(P2[i1:i2])
            max2=max(P2[i1:i2])
            ax[0, 1].add_patch(
            patches.Rectangle(
                (min1, min2),
                (max1-min1),
                (max2-min2),
                color=colors[i],
                fill=False      # remove background
                )
            )
            ax[0, 1].annotate(n[i], (max1, max2))
            minX=min([minX,min1,max1])
            maxX=max([maxX,min1,max1])
            minY=min([minY,min2,max2])
            maxY=max([maxY,min2,max2])
            l1 = lines.Line2D([min1,max1],[min2,max2], lw=1, color=colors[i], axes=ax[0, 1])
            ax[0, 1].add_line(l1)
             
        ax[0, 1].set_xlim([minX-abs(minX*0.1),maxX+abs(maxX*0.1)])
        ax[0, 1].set_ylim([minY-abs(minY*0.1),maxY+abs(maxY*0.1)])
 
    else:
        for i, txt in enumerate(n):
            ax[0, 1].annotate(txt, (newX[i*pperelem,0],newX[i*pperelem,1]))  
#      
  
      
    index = np.arange(2, n_features+2)
    bar_width=1/n_features
  
 
#back to original
    for i in range(0,n_samples):
        i1=i*npoint
        i2=(i+1)*npoint
        ax[1,0].scatter(Y[i1:i2,imp1],Y[i1:i2,imp2],color=colors[i])
    if(config[2]==1):
             
        minX=0
        minY=0
        maxX=0
        maxY=0
        for i in range(0,n_samples):
            i1=i*npoint
            i2=(i+1)*npoint
            P1=Y[:,imp1]
            P2=Y[:,imp2]
            min1=min(P1[i1:i2])
            max1=max(P1[i1:i2])
            min2=min(P2[i1:i2])
            max2=max(P2[i1:i2])
             
            ax[1, 0].add_patch(
            patches.Rectangle(
                (min1, min2),
                (max1-min1),
                (max2-min2),
                color=colors[i],
                fill=False      # remove background
                )
            )
            ax[1, 0].annotate(n[i], (max1, max2))
            minX=min([minX,min1,max1])
            maxX=max([maxX,min1,max1])
            minY=min([minY,min2,max2])
            maxY=max([maxY,min2,max2])
            l1 = lines.Line2D([min1,max1],[min2,max2], lw=1, color=colors[i], axes=ax[1, 0])
            ax[1, 0].add_line(l1)
 
        ax[1, 0].set_xlim([minX-abs(minX*0.1),maxX+abs(maxX*0.1)])
        ax[1, 0].set_ylim([minY-abs(minY*0.1),maxY+abs(maxY*0.1)])
 
    else:
        for i, txt in enumerate(n):
            ax[1, 0].annotate(txt, (Y[i*pperelem,imp1],Y[i*pperelem,imp2]))   
  
 
#     l1 = lines.Line2D([pca.components_[0,0]*4, pca.components_[1,0]*4], [pca.components_[0,0]*4, pca.components_[1,0]*4], transform=ax[0, 0].transFigure, figure=ax[0, 0])                              
#     l2 = lines.Line2D([pca.components_[0,1]*4, pca.components_[1,1]*4], [pca.components_[0,1]*4, pca.components_[1,1]*4], transform=ax[0, 0].transFigure, figure=ax[0, 0]) 
#
    
#     ax[0, 0].lines.extend([l1, l2])
     
#     l1 = lines.Line2D([pca.components_[0,0]*4, pca.components_[1,0]*4], [pca.components_[0,0]*4, pca.components_[1,0]*4],
#                     lw=1, color='black', axes=ax[0, 0])   
    l2 = lines.Line2D([pca.components_[0,0]*-8,pca.components_[0,0]*8],[pca.components_[1,0]*-8,pca.components_[1,0]*8],
                    lw=1, color='black', axes=ax[0, 0])
    l1 = lines.Line2D([pca.components_[0,1]*-8,pca.components_[0,1]*8],[pca.components_[1,1]*-8,pca.components_[1,1]*8],
                    lw=1, color='black', axes=ax[0, 0])
    l3 = lines.Line2D([-8,8],[0,0],
                    lw=1, color='grey', axes=ax[0, 0])
    l4 = lines.Line2D([0,0],[-4,4],
                    lw=1, color='grey', axes=ax[0, 0])
 
    
#     ax[0, 0].add_line(l1)
    ax[0, 0].add_line(l1)
    ax[0, 0].add_line(l2)
    ax[0, 0].add_line(l3)
    ax[0, 0].add_line(l4)
     
    l3 = lines.Line2D([-8,8],[0,0], lw=1, color='grey', axes=ax[0, 1])
    l4 = lines.Line2D([0,0],[-8,8],  lw=1, color='grey', axes=ax[0, 1])
    l5 = lines.Line2D([0,0],[-8,8],  lw=1, color='grey', axes=ax[1, 0])
    l6 = lines.Line2D([-8,8],[0,0],  lw=1, color='grey', axes=ax[1, 0])
    ax[0, 1].add_line(l3)
    ax[0, 1].add_line(l4)
    ax[1, 0].add_line(l5)
    ax[1, 0].add_line(l6)
#     plt.show()
 
#     print(pca.explained_variance_)
     
#     ax[0, 1].plot(np.arange(1, ncomp+1), pca.explained_variance_ratio_)
    rect =ax[1, 1].bar(index, pca.explained_variance_ratio_, bar_width,color='b')
    ax[1, 1].set_xlim(1, n_features+2)
    ax[1, 1].set_ylim(0, None)
 
     
    labels=['PC %s' %i for i in range(1,(n_features+1))]
    plt.xticks(index, labels)
     
          
    def autolabel(rects):
        i=0
        for rect in rects:
            height = rect.get_height()
            ax[1, 1].text(rect.get_x() + rect.get_width()/2., 1.05*height,
                    '%.2f' % pca.explained_variance_ratio_[i],
                    ha='center', va='bottom')
            i=i+1
   
    autolabel(rect)
  
  
#     ax[0, 0].xaxis.set_major_formatter(plt.NullFormatter())
#     ax[0, 1].xaxis.set_major_formatter(plt.NullFormatter())
#     ax[1, 1].xaxis.set_major_formatter(plt.NullFormatter())
#     ax[1, 0].xaxis.set_major_formatter(plt.NullFormatter())
      
    ax[0, 0].set_title('Input Data')
    ax[0, 1].set_title('First {0} Principal Vectors'.format(ncomp))
    ax[1, 0].set_title('Reconstructed Data ({0} components)'.format(ncomp))
    ax[1, 1].set_title('PCA variance ratio')
    ax[1, 1].set_xlabel('principal vector')
    ax[1, 1].set_ylabel('proportion of total variance')
    ax[0, 0].set_xlabel('{0}'.format(titles[imp1]))
    ax[0, 0].set_ylabel('{0}'.format(titles[imp2]))
    ax[1, 0].set_xlabel('{0}'.format(titles[imp1]))
    ax[1, 0].set_ylabel('{0}'.format(titles[imp2]))
 
      
    fig.suptitle('PCA', fontsize=16)
     
    plt.savefig('graphs/PCA.png', format='png', dpi=300)
    plt.show()
    plt.close()
 
#one graph with quantiles
def plot_resultsQuantile(PCA_res,n,objects,ratio):
    conf=1
    fig, ax = plt.subplots(figsize=(16*conf, 10*conf))   
     
    colors=[0 for x in range(len(objects))]
    cmap = plt.cm.get_cmap('Spectral')
    start=0
    end=1
    for i in range(len(objects)):
        colors[i] = cmap(start+((end-start)/((len(objects))+1))*(i+1))
 
 
    minX=0
    minY=0
    maxX=0
    maxY=0
    for i in range(len(objects)):
     
        if(n>1):
            for j in range(n-1):
                ax.arrow(PCA_res[n*i+j][0], PCA_res[n*i+j][1], (PCA_res[n*i+j+1][0]-PCA_res[n*i+j][0]), (PCA_res[n*i+j+1][1]-PCA_res[n*i+j][1]), length_includes_head=True, head_width=.005,color=colors[i])
#                 l1 = lines.Line2D([PCA_res[n*i+j][0],PCA_res[n*i+j+1][0]],[PCA_res[n*i+j][1],PCA_res[n*i+j+1][1]], lw=1, color=colors[i], axes=ax)
#                 ax.add_line(l1)
        else:
            ax.plot(PCA_res[i][0],PCA_res[i][1],'ro')
         
 
        min1=PCA_res[n*i][0]
        max1=PCA_res[n*i+n-1][0]
        min2=PCA_res[n*i][1]
        max2=PCA_res[n*i+n-1][1]
        ax.add_patch(
        patches.Rectangle(
            (min1, min2),
            (max1-min1),
            (max2-min2),
            fill=False,color=colors[i]      # remove background
            )
        )
        minX=min([minX,min1,max1])
        maxX=max([maxX,min1,max1])
        minY=min([minY,min2,max2])
        maxY=max([maxY,min2,max2])
 
 
 
 
#ANNOTATE
        ax.annotate(objects[i], (PCA_res[n*i+n-1][0],PCA_res[n*i+n-1][1]))
 
    ax.set_xlim([minX-abs(minX*0.1),maxX+abs(maxX*0.1)])
    ax.set_ylim([minY-abs(minY*0.1),maxY+abs(maxY*0.1)])
    ax.set_xlabel('%.2f' % ratio[0])
    ax.set_ylabel('%.2f' % ratio[1])
    plt.savefig('graphs/PCA_quantile.png', format='png', dpi=300)
    plt.show()
    plt.close()