'''
Created on Aug 18, 2014


'''

import array
from array import zeros
import math
import time
from com.sun.jmx.snmp.IPAcl import JJTParserState

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
        B_dim[0] = A_dim[int(dims[0])] # Creating a list that has dimensions of output
        B_dim[1] = A_dim[int(dims[1])]
        B_dim[2] = A_dim[int(dims[2])]
        rows_final = rows
        cols_final = cols
        slices_final = slices
    else: # Inverse permute
        B_dim[int(dims[0])] = A_dim[0]
        B_dim[int(dims[1])] = A_dim[1]
        B_dim[int(dims[2])] = A_dim[2]  
        rows_final = B_dim[0]
        cols_final = B_dim[1]
        slices_final = B_dim[2]
        
        
  
    for kk in xrange(slices_final):
        for ii in xrange(rows_final):
            for jj in xrange(cols_final):
                ind = [ii,jj,kk]
                ind_val = [0]*3
                ind_val[0] = ind[int(dims[0])] # ind[2] - slices
                ind_val[1] = ind[int(dims[1])] # ind[0] - rows
                ind_val[2] = ind[int(dims[2])] # ind[1] - columns
        
                if per_or_iper: # Permute
                    
                    B[int(ind_val[2]*B_dim[0]*B_dim[1] + ind_val[0]*B_dim[1] + ind_val[1])] = A[int(kk*rows*cols + ii*cols + jj)]
                else: # Inverse permute
                     
                    #B[int(kk*rows*cols + ii*cols + jj)] = A[int(ind_val[2]*B_dim[0]*B_dim[1] + ind_val[0]*B_dim[1] + ind_val[1])]
                    B[int(kk*B_dim[0]*B_dim[1] + ii*B_dim[1] + jj)] = A[int(ind_val[2]*rows*cols + ind_val[0]*cols + ind_val[1])]
   
    return B


def upsample(x , p , nr , nc):
    print 'useless'
    
# h is the column
# img is a 2D image
def convn(h, img, rows, cols):
  
    # Allocate space for output array
    hflip = h[::-1]
    outrows = (rows+len(h)-1)
    out = zeros('f', outrows*cols)    

    # Counter used to index into out
    cnt = 0
    
    # For each valid row location
    for y in range(-len(h)+1, rows):
        
        
        # For each valid column location
        for x in range(cols):
            s = 0 # Stores the summation
            counter = 0 # Access elements in the column vector
            
            
            
            for z in range(y,y+len(h)): # For each valid location within the column vector
                
                if z < 0: # If we are too far above, skip
                    counter += 1
                    continue
                elif z >= rows: # If we are too far below, skip
                    counter += 1
                    continue
                else: # Do the point by point multiplication and sum
                    s += hflip[counter]*img[z*cols + x]
                    counter += 1
        
            # Store sum in corresponding output location
            #print 'y:', y
            #print 'x:', x
            #print 'z:', z
            out[cnt] = s
            cnt += 1
    
    return (out, outrows)
    
def downsample(y , N1 , N2):
    
    counter = 0
    
    if N1 % 2 == 0:
    
        rows = N1/2
    
    else:
        rows = (N1/2) + 1
        
    out = zeros('f', rows*N2)
    
    for i in xrange(0,N1,2):
        for j in xrange(N2):
            out[counter] = y[int((i*N2)+j)]
            counter = counter + 1
            
    return (out , rows)       
            

    
    
# Function upfirdn
def upfirdn(xin,h, N1 , N2):
    
    (yout , r ) = convn(h,xin , N1 , N2)
    (yout, rout) = downsample(yout, r , N2)
    
    return yout
    
    

def medfilt2(A, rows, cols, rows_neigh, cols_neigh):
  
    out = zeros('f', int(rows*cols))
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
def afb3D_A(x, af, d, xx, yy, zz):
    
    
    hpf = array.array('f', [af[i] for i in xrange(1,len(af),2)])
    
    p = array.array('f', [(((d - 1) + val) % 3) for val in xrange(3)])
    
    
    size = zeros('i' , 3 )
    
    
   
    x = permute(x, xx, yy, zz , p , 1)
    L = 5
    
    if d==1:
        size[0] = xx
        xx = xx/2
        size[1] = yy
        size[2] = zz
    elif d==2:
        size[0] = yy
        yy = yy/2
        size[1] = zz
        size[2] = xx
    elif d==3:
        size[0] = zz
        zz = zz/2
        size[1] = xx
        size[2] = yy
    
    
    
    
    x = cshift3D(x, -L, 1, size);
    N1 = size[0]
    N2 = size[1]
    N3 = size[2]
    
    hi = zeros('f', (L+N1/2)*N2*N3)
     
    for k in xrange(N3):
        xTemp = zeros('f' , N1*N2)
        counter = 0
        for xloc in xrange(N1):
            for yloc in xrange(N2):
                xTemp[counter] = x[int(k*N1*N2 + xloc*N2 + yloc)]
                counter = counter + 1
        
    
        hiTemp = upfirdn(xTemp, hpf, N1, N2)    
        for xloc in xrange(L+N1/2):
            for yloc in xrange(N2):
                hi[int(k*(L+N1/2)*N2 + xloc*N2 + yloc)] = hiTemp[int(xloc*N2 + yloc)]
                
    
    hiTemp = hi
    for zloc in xrange(N3):
        for xloc in xrange(L):
            for yloc in xrange(N2):
                hiTemp[int(zloc*(L+N1/2)*N2 + xloc*N2 + yloc)] += hi[int(zloc*(L+N1/2)*N2 + (xloc + (N1/2))*N2 + yloc)]
    


    hi = zeros('f', (N1/2)*N2*N3)
    for zloc in xrange(N3):
        for xloc in xrange(N1/2):
            for yloc in xrange(N2):
                hi[int(zloc*(N1/2)*N2 + xloc*N2 + yloc)] = hiTemp[int(zloc*(L+N1/2)*N2 + xloc*N2 + yloc)]
 
    

    hi = permute(hi, (N1/2), N2 , N3 , p , 0)
    return (hi , xx , yy , zz) 

# Function afb3D
def afb3D(x, af1, af2, af3, x1, y1, z1):
    
    # filter along dimension 1
    (H, x2 , y2 , z2) = afb3D_A(x, af1, 1, x1, y1, z1);
    
    # filter along dimension 2
    (HH , x3 , y3 , z3) = afb3D_A(H, af2, 2, x2, y2, z2);
    
    
    (HHH , xu, yu , zu) = afb3D_A(HH, af3, 3, x3, y3, z3);
    
    
    hi = HHH
   
    
    return hi 




# Function dwt3D    
def dwt3D(x, af, xx, yy, zz):
    
    
    w = afb3D(x, af, af, af, xx, yy, zz) 
    
    
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
    xx = x
    yy = y
    zz = z
    if p1 == x and p2 ==y and p3 == z:
        pad1 = ima
        
    else:
        xx = p1
        yy = p2
        zz = p3
        p = int(p1*p2*p3)
        pad1 = zeros('f' , p) #(1)
        for zloc in xrange(z): # slices
            for xloc in xrange(x): # rows
                for yloc in xrange(y): #columns
                    pad1[int(zloc*p1*p2 + xloc*p2 + yloc)] = ima[int(zloc*x*y + xloc*y + yloc)]
    
               
    
    af = farras()
    
    # Simplified
    w = dwt3D(pad1,af, xx, yy, zz)
    print 'dwt3D done'
    # Variable is unused. Big  problem.
    locnois = w
    
    dims_half = [round(x/2), round(y/2), round(z/2)]
    num_elements_half = dims_half[0]*dims_half[1]*dims_half[2]
    
    tmp = zeros('f', int(num_elements_half))
    tmp_abs = zeros('f', int(num_elements_half))
    
    counter = 0
    for kk in xrange(z/2):
        for ii in xrange(x/2):
            for jj in xrange(y/2):
                
                tmp[counter] = locnois[int(kk*(x/2)*(y/2) + ii*(y/2) + jj)]
                tmp_abs[counter] = abs(tmp[counter])
                counter += 1
    
    
    SigArr = sorted(tmp_abs)
    
    if num_elements_half % 2 == 0:
        Sig = (SigArr[int(num_elements_half/2)]+SigArr[int((num_elements_half/2) - 1)])/(2*0.6745)
    else:
        Sig =  SigArr[int(num_elements_half/2)] / 0.6745  
    
    print Sig
    print 'Sig done'
    
    #NsigMAP = []
    #img2d = zeros('f', int(dims_half[0]*dims_half[1]))
    
    #for k in xrange(dims_half[2]):
        # Get the kth slice
        #for ii in xrange(dims_half[0]):
            #for jj in xrange(dims_half[1]):
                
                #img2d[int(ii*dims_half[1] + jj)] = tmp[int(k*dims_half[0]*dims_half[1] + ii +jj)]
                # absolute values
                #img2dp_abs = array.array('f' ,[abs(val) / 0.6745 for val in img2d])
        
        #print k        
        #LocalSTD = medfilt2(img2dp_abs, dims_half[0], dims_half[1], 2*ps + 1, 2*ps + 1)
        #print LocalSTD        
        #NsigMAP.append(LocalSTD)
                
    #print 'medfilt2 done'            
    # Interpolation code begins
    ### Output matrix
    HHH = array.zeros('f', int(x*y*z))
    
    # 2D slice
    #output_slc = array.zeros('f', int(x*y))
    
    #for k in xrange(dims_half[2]):
        # Get the kth slice
        
        #slc = NsigMAP[k]
        
        # interpolates in 2D
        #for xx in xrange(dims_half[0]):
            #for yy in xrange(dims_half[1]):
                #output_slc[int(xx*y + yy)] = slc[int(xx*dims_half[1] + yy)] # Original
                #output_slc[int((xx+1)*y + yy)] = slc[int(xx*dims_half[1] + yy)] # Down
                #output_slc[int(xx*y + (yy+1))] = slc[int(xx*dims_half[1] + yy)] # Right
                #output_slc[int((xx+1)*y + (yy+1))] = slc[int(xx*dims_half[1] + yy)] # Diagonal
        
        # interpolates in 3D
        # k = 0 --> [0,2) 0, 1
        # k = 1 --> [2,4) 2, 3
        # k = 2 --> [4,6) 4, 5
        #for kk in xrange(2*k,2*k+2):
            #for xx in xrange(x):
                #for yy in xrange(y):
                    #HHH[int(kk*x*y + xx*y + yy)] = output_slc[int(xx*y + yy)]
    
    
    for val in xrange(x*y*z):
        #if HHH[val] < Sig:
            #HHH[val] = Sig
        HHH[val] = Sig    
    
    return HHH