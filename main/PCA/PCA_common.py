#get all points
def itterateAll(A,so_far,nr,nrf,res):
    lcnr=nr
    l_so_far=so_far.copy()
    l_so_far1=so_far.copy()
    if(lcnr<nrf-1):
        l_so_far.append(A.getHistogram(lcnr).get_nonZeroS())
        itterateAll(A,l_so_far,lcnr+1,nrf,res)
        l_so_far1.append(A.getHistogram(lcnr).get_nonZeroE())
        itterateAll(A,l_so_far1,lcnr+1,nrf,res)
#         itterateAll(A,l_so_far1.append(A.getHistogram(lcnr).get_nonZeroE()),lcnr+1,nrf,res)
    else:
        l_so_far.append(A.getHistogram(lcnr).get_nonZeroS())
        res.append(l_so_far)
        l_so_far1.append(A.getHistogram(lcnr).get_nonZeroE())
        res.append(l_so_far1)
 
#get minimum points       
def itterateMin(A,so_far,nr,nrf,res):
    lcnr=nr
    l_so_far=so_far.copy()
    if(lcnr<nrf-1):
        l_so_far.append(A.getHistogram(lcnr).get_nonZeroS())
        itterateMin(A,l_so_far,lcnr+1,nrf,res)
    else:
        l_so_far.append(A.getHistogram(lcnr).get_nonZeroS())
        res.append(l_so_far)

#getmaxpoints      
def itterateMax(A,so_far,nr,nrf,res):
    lcnr=nr
    l_so_far=so_far.copy()
    if(lcnr<nrf-1):
        l_so_far.append(A.getHistogram(lcnr).get_nonZeroE())
        itterateMax(A,l_so_far,lcnr+1,nrf,res)
    else:
        l_so_far.append(A.getHistogram(lcnr).get_nonZeroE())
        res.append(l_so_far)