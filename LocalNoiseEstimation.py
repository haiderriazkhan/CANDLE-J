'''
Created on Aug 18, 2014


'''

import array
from array import zeros
import math
import time

# Function farras
def farras():
  
    af = zeros('f', 20)
    af[0] = 0.0
    af[1] = -0.01122679215254
  
    af[2] = 0.0
    af[3] = 0.01122679215254
  
    af[4] = -0.08838834764832
    af[5] = 0.08838834764832
  
    af[6] = 0.08838834764832
    af[7] = 0.08838834764832
  
    af[8] = 0.69587998903400
    af[9] = -0.69587998903400
  
    af[10] = 0.69587998903400
    af[11] = 0.69587998903400
  
    af[12] = 0.08838834764832
    af[13] = -0.08838834764832
  
    af[14] = -0.08838834764832
    af[15] = -0.08838834764832
  
    af[16] = 0.01122679215254 
    af[17] = 0.0
  
    af[18] = 0.01122679215254 
    af[19] = 0.0
    
    return af

# Function cshift3D
def cshift3D(x, m, d, size):
    N1 = size[0]
    N2 = size[1]
    N3 = size[2]
    if d == 1:
        n = xrange(N1) #(1)
        n = [(val - m) % N1 for val in n] #(2)
        out = array.array('f')
        for zloc in xrange(N3):
            for xloc in n:
                for yloc in xrange(N2):
                    out.append(x[int(zloc*N1*N2 + xloc*N2 + yloc)])
    
    elif d == 2:
        n = xrange(N2) #(1)
        n = [(val - m) % N2 for val in n] #(2)
        out = array.array('f')
        for zloc in xrange(N3):
            for xloc in xrange(N1):
                for yloc in n:
                    out.append(x[int(zloc*N1*N2 + xloc*N2 + yloc)])
    
    elif d == 3:
        n = xrange(N3) #(1)
        n = [(val - m) % N3 for val in n] #(2)
        out = array.array('f')
        for zloc in n:
            for xloc in xrange(N1):
                for yloc in xrange(N2):
                    out.append(x[int(zloc*N1*N2 + xloc*N2 + yloc)]) 
    
    return out 


# Function permute
def permute(A,order):
    print "To Do"
    
    
# Function upfirdn
def upfirdn(xin,h,p,q):
    print "To Do"    


# Function ipermute
def ipermute(B, order):
    print "To Do"




# Function afb3D_A
def afb3D_A(x, af, d):
    
    lpf = array.array('f', [af[i] for i in xrange(0,len(af),2)])
    hpf = array.array('f', [af[i] for i in xrange(1,len(af),2)])
    
    p = array.array('f', [(((d - 1) + val) % 3) + 1 for val in xrange(3)])
    
    # TO DO LATER
    (x, size) = permute(x, p)
    L = 5
    # TO DO LATER
    x = cshift3D(x, -L, 1, size);
    N1 = size[0]
    N2 = size[1]
    N3 = size[2]
    lo = zeros('f', (L+N1/2)*N2*N3)
    hi = zeros('f', (L+N1/2)*N2*N3)
    
    for k in xrange(N3):
        xTemp = array.array('f')
        for xloc in xrange(N1):
            for yloc in xrange(N2):
                xTemp.extend(x[k*N1*N2 + xloc*N2 + yloc]) # Corresponds to the kth slice
        
        
        loTemp = upfirdn(xTemp, lpf, 1, 2)    
        for xloc in xrange(N1):
            for yloc in xrange(N2):
                lo[k*N1*N2 + xloc*N2 + yloc] = loTemp[xloc*N2 + yloc]
    
    loTemp = lo
    for zloc in xrange(N3):
        for xloc in xrange(L):
            for yloc in xrange(N2):
                loTemp[zloc*N1*N2 + xloc*N2 + yloc] += lo[zloc*N1*N2 + (xloc + (N1/2))*N2 + yloc]
    
                
    
    lo = zeros('f', (N1/2)*N2*N3)
    for zloc in xrange(N3):
        for xloc in xrange(N1/2):
            for yloc in xrange(N2):
                lo[zloc*N1*N2 + xloc*N2 + yloc] = loTemp[zloc*N1*N2 + xloc*N2 + yloc]
    

    for k in xrange(N3):
        xTemp = array.array('f')
        for xloc in xrange(N1):
            for yloc in xrange(N2):
                xTemp.extend(x[k*N1*N2 + xloc*N2 + yloc])
        
    
        hiTemp = upfirdn(xTemp, hpf, 1, 2)    
        for xloc in xrange(N1):
            for yloc in xrange(N2):
                hi[k*N1*N2 + xloc*N2 + yloc] = hiTemp[xloc*N2 + yloc]
                
    
    hiTemp = hi
    for zloc in xrange(N3):
        for xloc in xrange(L):
            for yloc in xrange(N2):
                hiTemp[zloc*N1*N2 + xloc*N2 + yloc] += hi[zloc*N1*N2 + (xloc + (N1/2))*N2 + yloc]
    


    hi = zeros('f', (N1/2)*N2*N3)
    for zloc in xrange(N3):
        for xloc in xrange(N1/2):
            for yloc in xrange(N2):
                hi[zloc*N1*N2 + xloc*N2 + yloc] = hiTemp[zloc*N1*N2 + xloc*N2 + yloc]
 
    
    
    lo = ipermute(lo, p)
    hi = ipermute(hi, p)
    
    return (lo, hi) 

# Function afb3D
def afb3D(x, af1, af2, af3):
    
    # filter along dimension 1
    (L, H) = afb3D_A(x, af1, 1);
    
    # filter along dimension 2
    (LL , LH) = afb3D_A(L, af2, 2)
    (HL , HH) = afb3D_A(H, af2, 2);
    
    # filter along dimension 3
    (LLL , LLH) = afb3D_A(LL, af3, 3)
    (LHL , LHH) = afb3D_A(LH, af3, 3);
    (HLL , HLH) = afb3D_A(HL, af3, 3);
    (HHL , HHH) = afb3D_A(HH, af3, 3);
    
    lo = LLL
    hi = []
    hi.append(LLH)
    hi.append(LHL)
    hi.append(LHH)
    hi.append(HLL)
    hi.append(HLH)
    hi.append(HHL)
    hi.append(HHH)
    
    return (lo, hi) 




# Function dwt3D    
def dwt3D(x, J, af):
  
    w = [afb3D(x, af, af, af) for idx in xrange(J)]
    w.extend(x)
    
    return w




    
# Function estimate
def estimate(ima,x,y,z):
    
    minimum = reduce(min, ima)
    
    if (minimum < 0):
        ima = map(lambda x: x - minimum, ima)
     
    p1 =  math.pow(2 , math.ceil(math.log(x, 2)))
    p2 =  math.pow(2 , math.ceil(math.log(y, 2)))         
    p3 =  math.pow(2 , math.ceil(math.log(z, 2)))        
                
    # zeros pading 
    
    if p1 == x and p2 ==y and p3 == z:
        pad1 = array.array('f', ima)
        
    else:
        p = int(p1*p2*p3)
        pad1 = zeros('f' , p) #(1)
        for zloc in xrange(z): # slices
            for xloc in xrange(x): # rows
                for yloc in xrange(y): #columns
                    pad1[int(zloc*p1*p2 + xloc*p2 + yloc)] = ima[int(zloc*x*y + xloc*y + yloc)]
    
               
    
    af = farras()