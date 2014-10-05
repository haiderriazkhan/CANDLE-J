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


# Function permute and ipermute
def permute(A, rows, cols, slices, dims, per_or_iper):
    
    num_elements = len(A)
    B = zeros('f', num_elements)
    A_dim = [rows, cols, slices] # Creating a list that has the dimensions of input
    B_dim = [0]*3
  
    if per_or_iper: # Permute
        B_dim[0] = A_dim[dims[0]] # Creating a list that has dimensions of output
        B_dim[1] = A_dim[dims[1]]
        B_dim[2] = A_dim[dims[2]]
        rows_final = rows
        cols_final = cols
        slices_final = slices
    else: # Inverse permute
        B_dim[dims[0]] = 0
        B_dim[dims[1]] = 1
        B_dim[dims[2]] = 2  
        rows_final = B_dim[0]
        cols_final = B_dim[1]
        slices_final = B_dim[2]
  
    for ii in xrange(rows_final):
        for jj in xrange(cols_final):
            for kk in xrange(slices_final):
                ind = [ii,jj,kk]
                ind_val = [0]*3
                ind_val[0] = ind[dims[0]] # ind[2] - slices
                ind_val[1] = ind[dims[1]] # ind[0] - rows
                ind_val[2] = ind[dims[2]] # ind[1] - columns
        
            if per_or_iper: # Permute
                B[int(kk*B_dim[0]*B_dim[1] + ii*B_dim[1] + jj)] = A[int(kk*rows*cols + ii*cols + jj)]
            else: # Inverse permute
                B[int(kk*rows*cols + ii*cols + jj)] = A[int(kk*B_dim[0]*B_dim[1] + ii*B_dim[1] + jj)]
        
   
    return B
    
    
# Function upfirdn
def upfirdn(xin,h,p,q):
    print "To Do"    


def medfilt2(A, rows, cols, rows_neigh, cols_neigh):
  
    out = zeros('f', rows*cols)
    rows_offset = int(rows_neigh/2)
    cols_offset = int(cols_neigh/2)
    total_elements = rows_neigh*cols_neigh
    middle_index = int(total_elements/2)
  
    # For each pixel in the slice
    for ii in xrange(rows_offset,rows-rows_offset):
        for jj in xrange(cols_offset,cols-cols_offset):
      
            # Grabs the neighbourhood of pixels centered
            # at (ii,jj)
            neigh_pix = []
            for yy in xrange(ii-rows_offset,ii+rows_offset+1):
                for xx in xrange(jj-cols_offset,jj+cols_offset+1):
                    neigh_pix.append(A[int(yy*cols + xx)])
          
            # Sort the values
            sorted_values = sorted(neigh_pix)
            out[int(ii*cols+jj)] = sorted_values[middle_index]
  
    return out




# Function afb3D_A
def afb3D_A(x, af, d):
    
    lpf = array.array('f', [af[i] for i in xrange(0,len(af),2)])
    hpf = array.array('f', [af[i] for i in xrange(1,len(af),2)])
    
    p = array.array('f', [(((d - 1) + val) % 3) + 1 for val in xrange(3)])
    
    # Problem - Getting size 
    (x, size) = permute(x, p , 1)
    L = 5
    
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
 
    
    # Need double checking
    lo = permute(lo,(N1/2), N2 , N3, p , 0)
    hi = permute(hi, (N1/2), N2 , N3 , p , 0)
    
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
    
    # What is going on here?
    w = [afb3D(x, af, af, af) for idx in xrange(J)]
    w.extend(x)
    
    return w




    
# Function estimate
def estimate(ima,x,y,z , ps):
    
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
    w1 = dwt3D(pad1,1,af)
    # Variable is unused. Big fucking problem.
    tmp = w1[6]
    dims_half = [round(x/2), round(y/2), round(z/2)]
    num_elements_half = dims_half[0]*dims_half[1]*dims_half[2]
    tmp = array.zeros('f', num_elements_half)
    tmp_abs = array.zeros('f', num_elements_half)
    
    counter = 0
    for ii in xrange(0,x,2):
        for jj in xrange(0,y,2):
            for kk in xrange(0,z,2):
                tmp[counter] = ima[int(kk*x*y + ii*y + jj)]
                tmp_abs[counter] = abs(tmp[counter])
                counter += 1
                
                
    Sig = sorted(tmp_abs)[int(num_elements_half/2)] / 0.6745
    
    NsigMAP = []
    img2d = zeros('f', dims_half[0]*dims_half[1])
    
    for k in dims_half[2]:
        # Get the kth slice
        for ii in xrange(dims_half[0]):
            for jj in xrange(dims_half[1]):
                
                img2d[int(ii*dims_half[1] + jj)] = tmp[int(k*dims_half[0]*dims_half[1] + ii +jj)]
                # Median filter slice
                img2dp = medfilt2(img2d, dims_half[0], dims_half[1], 2*ps + 1, 2*ps + 1)
                # Median filter absolute values
                img2dp_abs = array.array('f' ,[abs(val) / 0.6745 for val in img2dp])
                
                LocalSTD = medfilt2(img2dp_abs, dims_half[0], dims_half[1], 2*ps + 1, 2*ps + 1)
                NsigMAP.append(LocalSTD)
    # Interpolation code begins
    ### Output matrix
    HHH = array.zeros('f', x*y*z)
    
    # 2D slice
    output_slc = array.zeros('f', x*y)
    
    for k in xrange(dims_half[2]):
        # Get the kth slice
        slc = NsigMAP[k]
        
        # interpolates in 2D
        for xx in xrange(dims_half[0]):
            for yy in xrange(dims_half[1]):
                output_slc[int(xx*y + yy)] = slc[int(xx*dims_half[1] + yy)] # Original
                output_slc[int((xx+1)*y + yy)] = slc[int(xx*dims_half[1] + yy)] # Down
                output_slc[int(xx*y + (yy+1))] = slc[int(xx*dims_half[1] + yy)] # Right
                output_slc[int((xx+1)*y + (yy+1))] = slc[int(xx*dims_half[1] + yy)] # Diagonal
        
        # interpolates in 3D
        # k = 0 --> [0,2) 0, 1
        # k = 1 --> [2,4) 2, 3
        # k = 2 --> [4,6) 4, 5
        for kk in xrange(2*k,2*k+2):
            for xx in xrange(x):
                for yy in xrange(y):
                    HHH[int(kk*x*y + xx*y + yy)] = output_slc[int(xx*y + yy)]
    
    
    for val in xrange(x*y*z):
        if HHH[val] < Sig:
            HHH[val] = Sig
    