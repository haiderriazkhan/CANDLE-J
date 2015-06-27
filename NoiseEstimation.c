/* Haider Khan - haiderriazkhan@hotmail.com                               */
/* Ruthazer Lab, Montreal Neurological Institute.                         */
/* McGill University                                                      */
/*                                                                        */
/* The algorithms in this program were developed with the guidance of     */
/* Raymond Phan (ray@bublcam.com).                                        */
/*                                                                        */
/* Copyright (C) 2015 Haider Riaz Khan                                    */

/*                                                                        */
/***************************************************************************
 *  The Wavelet-based local noise estimation is described in:              *
 *                                                                         *
 *  P. Coupe, Martin Munz, Jose V.Manjon, Edward Ruthazer, D.Louis Collins.*
 *  A CANDLE for a deeper in-vivo insight. Medical Image Analysis,         *
 *  16(4):849-64 (2012).                                                   *
 ***************************************************************************/


#include "math.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>




// If x is a power of 2, the method simply retunrs x.
// If x is not a power of 2, the method returns the next greatest power of 2.
int lp2 (unsigned int x) {
    if (x == 0) return 0;
    unsigned int result = 1;
    while ((result < x) && (result != 0))
        result <<= 1;
    return (int) result;
}




// Returns the minimum value of the input array
float minimum(int size, float A[static size] ){
    
    float min = A[0];
    
    for (int i =1; i < size; i++) {
        
        if (A[i] < min) {
            
            min = A[i];
            
        }
        
    }
    
    return min;
    
    
}


// Sorts the array using the quick sort algorithm
void quick_sort (float *a, int n) {
    int i, j;
    float p, t;
    if (n < 2)
        return;
    p = a[n / 2];
    for (i = 0, j = n - 1;; i++, j--) {
        while (a[i] < p)
            i++;
        while (p < a[j])
            j--;
        if (i >= j)
            break;
        t = a[i];
        a[i] = a[j];
        a[j] = t;
    }
    quick_sort(a, i);
    quick_sort(a + i, n - i);
}


// Returns the median value of the input array
float median(int size , float A[static size]){
    
    float *Sorted, med;
    
    Sorted = (float*)malloc(sizeof(float) * size );
    
    memcpy(Sorted , A, sizeof(float) * size ) ;
    
    quick_sort(Sorted , size);
    
    if (size % 2) {
        
        med = Sorted[size/2] ;
        
    }else{
        
        med = ( Sorted[size/2] + Sorted[(size/2) - 1] ) / 2 ;
        
        
    }
    
    free(Sorted);
    return med;
    
}





// 3D circular shift (along dimension 1)
float * cshift3D(float*x , int N1 , int N2 , int N3 ){
    int i, j , k, counter;
    float *y;
    
    int n[N1];
    
    y = (float*)malloc(sizeof(float) * ( N1 * N2 * N3 ) );
    
    for (i = 0; i < N1; i++) {
        
        n[i] = (i + 5 ) % N1;
    }
    
    counter = 0;
    for (k=0; k < N3; k++){
        for (i=0; i < N1; i++){
            for (j=0; j < N2; j++){
                y[counter] = x[(k*N1*N2) + (n[i]*N2) + j];
                counter++;
                
            }
            
        }
        
    }
    
    free(x);
    return y;
    
}





// Rearranges dimension of A as specified
float * permute(float*A , int rows ,int cols ,int slices ,int p[static 3] ,int per_or_iper){
    
    int ii, jj , kk, rows_final , cols_final , slices_final;
    float *B;
    
    int A_dim[3] = {rows, cols , slices};
    
    int B_dim[3];
    
    B = (float*)malloc(sizeof(float) * ( rows*cols*slices ) );
    
    if (per_or_iper) {
        
        B_dim[0] = A_dim[p[0]] ;
        B_dim[1] = A_dim[p[1]] ;
        B_dim[2] = A_dim[p[2]] ;
        
        rows_final = rows;
        cols_final = cols;
        slices_final = slices;
        
        
        
    }else{
        
        B_dim[p[0]] = A_dim[0];
        B_dim[p[1]] = A_dim[1];
        B_dim[p[2]] = A_dim[2];
        
        rows_final = B_dim[0];
        cols_final = B_dim[1];
        slices_final = B_dim[2];
    
    }
    
    
    
    
    
    int ind_val[3];
    int ind[3];
    for (kk=0; kk<slices_final; kk++){
        ind[2] = kk;
        for (ii=0; ii<rows_final; ii++){
            ind[0] = ii;
            for (jj=0; jj<cols_final; jj++){
                ind[1] = jj;
                ind_val[0] = ind[p[0]];
                ind_val[1] = ind[p[1]];
                ind_val[2] = ind[p[2]];
    
                if (per_or_iper) {
                    
                    B[ind_val[2]*B_dim[0]*B_dim[1] + ind_val[0]*B_dim[1] + ind_val[1]] = A[kk*rows*cols + ii*cols + jj];
                    
                }else{
                
                    B[kk*B_dim[0]*B_dim[1] + ii*B_dim[1] + jj] = A[ind_val[2]*rows*cols + ind_val[0]*cols + ind_val[1]];
                    
                
                }
                
            }
            
        }
        
    }
    
    
    
    free(A);
    return B;
    
}

// N-D convolution algorithm similar to the one provided by MatLab
// The pre-computed convolution kernel is initialized within the method
// Also downsamples the output at the end
float * convn(float *const xin, int rows, int cols){
    int outrows, temprows, count, counter, x , y , z , i , j;
    float s;
    float *out , *temporary;
    
    float hpf[10] = {0 , 0 , -0.08838834764832 , -0.08838834764832 , 0.69587998903400 , -0.69587998903400, 0.08838834764832, 0.08838834764832, 0.01122679215254 , -0.01122679215254} ;
    
    temprows = rows + 9;
    
    temporary = (float*)malloc(sizeof(float) * ( temprows * cols ) );
    
    
    
    count = 0;
    
    
    for (y = -9 ; y < rows; y++) {
        for (x = 0; x < cols; x++) {
            
            s = 0;
            counter = 0;
            
            for (z = y; z < (y+10); z++) {
                
                if (z < 0) {
                    counter++;
                    continue;
                }
                else if(z >= rows){
                    counter++;
                    continue;
                }else{
                    
                    s += ( hpf[counter++] * xin[z*cols+x] ) ;
                    
                
                }
                
                
            }
            
            temporary[count++] = s;
        }
        
    }
    
    count = 0;
    
   
        
    
    outrows = (temprows/2) + 1 ;
    
    
    out = (float*)malloc(sizeof(float) * ( outrows * cols ) );
    
    for (i = 0; i < temprows; i += 2) {
        for (j = 0; j < cols; j++) {
            out[count++] = temporary[(i*cols) + j];
        }
    }
    
    free(temporary);
    return out;
    
    
}





// 3D Filter Bank
float * afb3D_A(float*x , int d, int xx, int yy, int zz){

    int L, i, N1, N2, N3 , zloc, xloc, yloc, counter;
    
    float *perm,  *cshif, *iperm, *hi , *xTemp;
    
    
    int p[3];
    for(i = 0; i < 3; i++ ){
        
        p[i] = ((d-1)+i) % 3 ;
        
        
    }
    
    
    
    
    
    perm = permute(x , xx , yy , zz , p , 1);
    
    
    

    L = 5;
    
    
    
    if (d == 1) {
        N1 = xx;
        N2 = yy;
        N3 = zz;
        
    }
    else if(d == 2){
        N1 = yy;
        N2 = zz;
        N3 = xx;
    
    
    }else{
        
        N1 = zz;
        N2 = xx;
        N3 = yy;
    
    }
    
    
    
    cshif = cshift3D(perm , N1, N2 , N3);
    
    
    hi = (float*)malloc(sizeof(float) * ( (L+(N1/2) )* N2 * N3 ) );
    xTemp = (float*)malloc(sizeof(float) * ( N1 * N2 ) );
    float *hiTemp;
    
    for (zloc = 0; zloc < N3; zloc++) {
        counter = 0;
        for (xloc = 0; xloc <  N1; xloc++) {
            for (yloc = 0; yloc < N2; yloc++) {
                xTemp[counter++] = cshif[zloc*N1*N2 + xloc*N2 + yloc];
            }
        }
        hiTemp = convn(xTemp, N1, N2);
        for (xloc = 0; xloc <  (L+(N1/2)); xloc++) {
            for (yloc = 0; yloc < N2; yloc++) {
                hi[zloc*(L+N1/2)*N2 + xloc*N2 + yloc] = hiTemp[xloc*N2 + yloc];
            }
        }
        
        free(hiTemp);
        
    }
    
    
    free(xTemp);
    free(cshif);
    
    
    
    
    hiTemp = (float*)malloc(sizeof(float) * ( (L+(N1/2) )* N2 * N3) );
    memcpy(hiTemp , hi , sizeof(float) * ( (L+(N1/2) )* N2 * N3 ) );
           
    
    
    for (zloc = 0; zloc < N3; zloc++){
        for (xloc = 0; xloc <  L; xloc++){
            for (yloc = 0; yloc < N2; yloc++){
                hiTemp[zloc*(L+N1/2)*N2 + xloc*N2 + yloc] += hi[zloc*(L+N1/2)*N2 + (xloc + (N1/2))*N2 + yloc];
           }
        }
    }
    
           
    hi = (float*)realloc(hi , sizeof(float) * ( (N1/2) * N2 * N3) );
    for (zloc = 0; zloc < N3; zloc++){
        for (xloc = 0; xloc <  (N1/2); xloc++){
            for (yloc = 0; yloc < N2; yloc++){
                hi[zloc*(N1/2)*N2 + xloc*N2 + yloc] = hiTemp[zloc*(L+N1/2)*N2 + xloc*N2 + yloc];
            }
        }
    }
    
    free(hiTemp);
    
    iperm = permute(hi, (N1/2), N2 , N3 , p , 0);
    
    
    
    return iperm;
    
    

}



// Pads the 2-D image by circular repitition of the floating point numbers
void pad2d(float *arr , float *newarr , int x , int y , int axis){
    
    int xx, yy , rowcoun, colcoun, i , row , j , col;
    
    
    xx = x + 2*axis;
    yy = y + 2*axis;
    
    rowcoun = 0;
    
    for(i=0; i < xx; i++){
        colcoun = 0;
        if(i <= (axis - 1)){
            
            row = axis - rowcoun;
            row = row % x;
            
            if (row){
                
                row -= 1;
                row = (x-1) - row;
                
            }
            rowcoun++;
            
            
        }
        else if(i <= (x+axis-1) ){
            
            row = i - axis;
        
        }else{
            
            row = rowcoun  - axis + 1 ;
            row = row % x ;
            
            if(row){
                
                row -= 1;
            
            }else{
                
                row = x - 1;
            }
        
            rowcoun++;
        }
    
        for(j=0; j < yy; j++){
            
            
            
            if(j <= (axis-1) ){
                
                col = axis - colcoun ;
                col = col % y ;
                
                
                if(col){
                    
                    col -= 1;
                    col = (y-1) - col;
                    
                }
                
                colcoun++;
            }
            else if(j <= (y+axis-1)){
                
                col = j - axis;
                
            
            }else{
                
                col = colcoun  - axis + 1;
                col = col % y ;
                if(col){
                    
                    col -= 1;
                    
                }else{
                    
                    col = y - 1;
                }
                
                colcoun++;
            }
            
            
            newarr[yy*i + j] = arr[row*y + col];
            
        
        }
        
        
    }
    
    
}

// 2-D Median Filter
// Array A is already padded
// Array B is the output
void medfilt2d(float *A, float *B, int rows, int cols, int axis){
    
    int winelem , ii , jj , xx , yy , inc, ind;
    
    
    
    winelem =  ( (2*axis) + 1) * ( (2*axis) + 1) ;
    
    float window[winelem];
    
    
    for(ii = 0; ii < rows; ii++){
        
        for (jj=0; jj < cols; jj++) {
            
            inc = 0;
            
            for (xx = 0; xx < ( (2*axis) + 1); xx++) {
                for (yy = 0; yy < ( (2*axis) + 1); yy++) {
                    
                    ind = ( (ii+xx)* ((2*axis) + cols ) + (jj + yy)) ;
                    window[inc++] = A[ind];
                }
            }
            
            B[ii*cols + jj] = median(winelem , window);
            
        }
        
    
    }
    

}



// Parent method. 
void estimate(float*ima , int x , int y , int z , int ps, float **HHH){
    
    int size , xx , yy , zz , i, j, k, p1 , p2 , p3, counter , z_half , x_half, y_half , padx , pady, indexx, zi, xi, yi, val ;
    float minim, Sig;
    
    float *filt1, *filt2, *filt3, *padarray, *temp, *NsigMAP , *img2d , *img2dp;
    
    size = x*y*z;
    
    
    
    minim = minimum(size , ima);
    
    if (minim < 0) {
        
        for(i=0; i<size; i++){
            
            ima[i] = ima[i] - minim;
        }
        
    }
    
    // Assign values to p1, p2 and p3
    
    
    p1 = lp2( (unsigned int)x ) ;
    
    p2 = lp2( (unsigned int)y ) ;
    
    p3 = lp2( (unsigned int)z ) ;
    
    
    
    
    
    // Make the image dimesnions powers of 2
    if(p1 == x & p2==y & p3 == z){
        
        padarray = (float*)malloc(sizeof(float) * (p1*p2*p3));
        memcpy(padarray , ima , sizeof(float) * (p1*p2*p3));
    
    
    }else{
    
        
        
        padarray = (float*)calloc((p1*p2*p3) , sizeof(float));
        
        for(k = 0; k < z; k++ ){
            
            for(i=0; i < x; i++){
                
                for(j=0; j < y; j++){
                    
                    padarray[(k*p1*p2)+(i*p2)+ j] = ima[(k*x*y) + (i*y) + j];
                }
            }
        }
    }
    
    
    
    // Filter along dimension 1
    
    
    
    filt1 = afb3D_A(padarray, 1 , p1, p2, p3);
    
    
    p1 = p1/2;
    
    
    
    // Filter along dimension 2
    
    
    filt2 = afb3D_A(filt1, 2 , p1, p2, p3);
    
    p2 = p2/2 ;
    
    
    
    // Filter along dimension 3
    
    filt3 = afb3D_A(filt2, 3 , p1, p2, p3);
    
    p3 = p3/2;
    
    
    
    
    z_half = z/2 + (z % 2 != 0);
    x_half = x/2 + (x % 2 != 0);
    y_half = y/2 + (y % 2 != 0);
    
    // Remove Regions corresponding to zero padding
    temp = (float*)malloc(sizeof(float) * ( z_half*x_half*y_half ) );
    
    
    
    counter = 0;
    for(k=0; k < z_half; k++){
        for (i = 0; i < x_half; i++){
            for (j=0; j < y_half; j++){
                
                temp[counter++] = fabsf( filt3[k*p1*p2 + i*p2 + j] ) / 0.6745;
                
                
            }
        
        }
    }
    
    
    free(filt3);
    Sig = median( (z_half*x_half*y_half) , temp );
    printf("Sig:\n");
    printf("%.6f", Sig);
    padx = x_half + 2*ps;
    pady = y_half + 2*ps;
    
    
    img2d = (float*)malloc(sizeof(float) * ( x_half*y_half ) );
    img2dp = (float*)malloc(sizeof(float) * ( padx * pady ) );
    
    
    NsigMAP = (float*)malloc(sizeof(float) * ( z_half * x_half * y_half ) ) ;
    
    
    for(k=0; k < z_half; k++){
        counter = 0;
        for(i=0; i < x_half; i++ ){
            for(j=0; j < y_half; j++){
                
                // Get the kth Slice
                img2d[counter++] = temp[k*x_half*y_half + i*y_half + j];
                
            }
        
        }
        
        
        
        memset(img2dp,0, sizeof(float) * (padx*pady) );
        
        
        // Pad the 2d image
        pad2d(img2d , img2dp , x_half , y_half, ps);
        
        // Apply 2d median filter
        medfilt2d(img2dp , img2d,  x_half , y_half , ps);
        
        indexx = k*x_half*y_half ;
        
        memcpy( &NsigMAP[indexx] , img2d, sizeof(float) * (x_half*y_half) );
        
    }
    
    free(temp);
    free(img2d);
    free(img2dp);
    
    *HHH = (float*)malloc(sizeof(float) * ( z * x * y ) );
    
    
    // 3-d interpolation
    
    for (k=0; k < z; k++) {
        
        zi = k/2;
        
        for (i=0; i < x; i++) {
            
            xi = i/2;
            
            for (j=0; j < y; j++) {
                
                yi = j/2;
                
                if (k && i && j) {
                    
                    (*HHH)[k*x*y + y*i + j] = NsigMAP[zi*x_half*y_half + xi*y_half + yi] ;
                    
                }else{
                    
                    (*HHH)[k*x*y + y*i + j] = NAN;
                    
                }
                
            }
        }
    }
    
    
    
    free(NsigMAP);
    
    for (val = 0; val < (x*y*z); val++) {
        
        if ((*HHH)[val] < Sig) {
            (*HHH)[val] = Sig;
        }
    }
    
    
    
    
}




