/* Pierrick Coupe - pierrick.coupe@gmail.com                               */
/* Jose V. Manjon - jmanjon@fis.upv.es                                     */
/* Brain Imaging Center, Montreal Neurological Institute.                  */
/* Mc Gill University                                                      */
/*                                                                         */
/* Copyright (C) 2010 Pierrick Coupe and Jose V. Manjon                    */

/*                          Details on ONLM filter                        */
/***************************************************************************
 *  The ONLM filter is described in:                                       *
 *                                                                         *
 *  P. Coupe, P. Yger, S. Prima, P. Hellier, C. Kervrann, C. Barillot.     *
 *  An Optimized Blockwise Non Local Means Denoising Filter for 3D Magnetic*
 *  Resonance Images. IEEE Transactions on Medical Imaging, 27(4):425-441, *
 *  Avril 2008                                                             *
 ***************************************************************************/


#include "math.h"
//#include "mex.h"
#include <stdlib.h>
//#include "matrix.h"

// Include the JNI
#include <jni.h>
#include <stdio.h>

/* Multithreading stuff*/
#ifdef _WIN32
#include <windows.h>
#include <process.h>
#else
#include <pthread.h>
#endif

typedef struct {
    int rows;
    int cols;
    int slices;
    float * in_image;
    float * ref_image;
    float * means_image;
    float * var_image;
    float * estimate;
    float * mask_image;
    unsigned short * label;
    int ini;
    int fin;
    int radioB;
    int radioS;
    float *sigma;
    float beta;
} myargument;


float max;

/* Function which compute the weighted average for one block */
void Average_block(float *ima,int x,int y,int z,int neighborhoodsize,float *average, float weight, int sx,int sy,int sz)
{
    int x_pos,y_pos,z_pos;
    bool is_outside;
    int a,b,c,ns,sxy,count;
    
    ns=2*neighborhoodsize+1;
    sxy=sx*sy;
    
    count = 0;
    
    for (c = 0; c<ns;c++)
    {
        z_pos = z+c-neighborhoodsize;
        for (b = 0; b<ns;b++)
        {
            y_pos = y+b-neighborhoodsize;
            for (a = 0; a<ns;a++)
            {
                x_pos = x+a-neighborhoodsize;
                
                is_outside = false;
                if ((x_pos < 0) || (x_pos > sx-1)) is_outside = true;
                if ((y_pos < 0) || (y_pos > sy-1)) is_outside = true;
                if ((z_pos < 0) || (z_pos > sz-1)) is_outside = true;
                
                
                if (is_outside)
                    average[count] = average[count] + ima[z*(sxy)+(y*sx)+x]*weight;
                else
                    average[count] = average[count] + ima[z_pos*(sxy)+(y_pos*sx)+x_pos]*weight;
                
                count++;
            }
        }
    }
}

/* Function which computes the value assigned to each voxel */
void Value_block(float *Estimate, unsigned short *Label,int x,int y,int z,int neighborhoodsize,float *average, float global_sum, int sx,int sy,int sz)
{
    int x_pos,y_pos,z_pos;
    int ret;
    bool is_outside;
    float value = 0.0;
    unsigned short label = 0;
    int count=0 ;
    int a,b,c,ns,sxy;
    
    extern bool rician;
    
    ns=2*neighborhoodsize+1;
    sxy=sx*sy;
    
    for (c = 0; c<ns;c++)
    {
        z_pos = z+c-neighborhoodsize;
        for (b = 0; b<ns;b++)
        {
            y_pos = y+b-neighborhoodsize;
            for (a = 0; a<ns;a++)
            {
                x_pos = x+a-neighborhoodsize;
                
                is_outside = false;
                if ((x_pos < 0) || (x_pos > sx-1)) is_outside = true;
                if ((y_pos < 0) || (y_pos > sy-1)) is_outside = true;
                if ((z_pos < 0) || (z_pos > sz-1)) is_outside = true;
                
                if (!is_outside)
                {
                    value = Estimate[z_pos*(sxy)+(y_pos*sx)+x_pos];
                    
                    value = value + (average[count]/global_sum);
                    
                    label = Label[(x_pos + y_pos*sx + z_pos*sxy)];
                    Estimate[z_pos*(sxy)+(y_pos*sx)+x_pos] = value;
                    Label[(x_pos + y_pos*sx + z_pos *sxy)] = label +1;
                }
                count++;
            }
        }
    }
}

// nx - columns
// ny - rows
// nz - slices
float distance(float* ima,int x,int y,int z,int nx,int ny,int nz,int f,int sx,int sy,int sz)
{
    float d,acu,distancetotal;
    int i,j,k,ni1,nj1,ni2,nj2,nk1,nk2;
    
    distancetotal=0;
    
    for(k=-f;k<=f;k++)
    {
        nk1=z+k;
        nk2=nz+k;
        if(nk1<0) nk1=-nk1;
        if(nk2<0) nk2=-nk2;
        if(nk1>=sz) nk1=2*sz-nk1-1;
        if(nk2>=sz) nk2=2*sz-nk2-1;
        
        for(j=-f;j<=f;j++)
        {
            nj1=y+j;
            nj2=ny+j;
            if(nj1<0) nj1=-nj1;
            if(nj2<0) nj2=-nj2;
            if(nj1>=sy) nj1=2*sy-nj1-1;
            if(nj2>=sy) nj2=2*sy-nj2-1;
            
            for(i=-f;i<=f;i++)
            {
                ni1=x+i;
                ni2=nx+i;
                if(ni1<0) ni1=-ni1;
                if(ni2<0) ni2=-ni2;
                if(ni1>=sx) ni1=2*sx-ni1-1;
                if(ni2>=sx) ni2=2*sx-ni2-1;
                
                distancetotal = distancetotal + ((ima[nk1*(sx*sy)+(nj1*sx)+ni1]-ima[nk2*(sx*sy)+(nj2*sx)+ni2])*(ima[nk1*(sx*sy)+(nj1*sx)+ni1]-ima[nk2*(sx*sy)+(nj2*sx)+ni2]));
            }
        }
    }
    
    acu=(2*f+1)*(2*f+1)*(2*f+1);
    d=distancetotal/acu;
    
    return d;
    
}

void* ThreadFunc( void* pArguments )
{
    float *Estimate,*ima,*means,*variances,*average,*sigma,*ref,*mask,epsilon,totalweight,wmax,t1,d,w;
    float beta;
    unsigned short *Label;
    int rows,cols,slices,ini,fin,v,f,init,i,j,k,rc,ii,jj,kk,ni,nj,nk,Ndims;
    
    myargument arg;
    arg=*(myargument *) pArguments;
    
    rows=arg.rows;
    cols=arg.cols;
    slices=arg.slices;
    ini=arg.ini;
    fin=arg.fin;
    ima=arg.in_image;
    ref=arg.ref_image;
    means=arg.means_image;
    variances=arg.var_image;
    mask = arg.mask_image;
    Estimate=arg.estimate;
    Label=arg.label;
    v=arg.radioB;
    f=arg.radioS;
    sigma=arg.sigma;
    beta= arg.beta;
    
    epsilon = 0.0000001;
    init = 0;
    rc=rows*cols;
    
    Ndims = (2*f+1)*(2*f+1)*(2*f+1);
    
    average=(float*)malloc(Ndims*sizeof(float));
    
    for(k=ini;k<fin;k+=2)
    {
        for(j=0;j<rows;j+=2) // j is rows
        {
            for(i=0;i<cols;i+=2) // i is columns
            {
                for (init=0 ; init < Ndims; init++) average[init]=0.0;
                
                totalweight=0.0;
                
                wmax=0.0;
                
                
                if(mask[k*rc+(j*cols)+i]>0)
                {
                    for(kk=-v;kk<=v;kk++)
                    {
                        for(jj=-v;jj<=v;jj++)
                        {
                            for(ii=-v;ii<=v;ii++)
                            {
                                ni=i+ii; // columns
                                nj=j+jj; // rows
                                nk=k+kk; // slices
                                
                                if(ii==0 && jj==0 && kk==0) continue;
                                
                                if(ni>=0 && nj>=0 && nk>=0 && ni<cols && nj<rows && nk<slices)
                                {
                                    t1 = ((2 * means[k*rc+(j*cols)+i] * means[nk*rc+(nj*cols)+ni]+ epsilon) / ( means[k*rc+(j*cols)+i]*means[k*rc+(j*cols)+i] + means[nk*rc+(nj*cols)+ni]*means[nk*rc+(nj*cols)+ni] + epsilon) )  * ((2 * sqrt(variances[k*rc+(j*cols)+i]) * sqrt(variances[nk*rc+(nj*cols)+ni]) + epsilon) / (variances[k*rc+(j*cols)+i] + variances[nk*rc+(nj*cols)+ni] + epsilon));
                                    if(t1 > 0.9)
                                    {
                                        
                                        d=distance(ref,i,j,k,ni,nj,nk,f,cols,rows,slices);
                                        if (d<=0) d = epsilon;
                                        
                                        w = exp(-d/(beta*sigma[k*(rc)+(j*cols)+i]*sigma[k*(rc)+(j*cols)+i]));
                                        
                                        if(w>wmax) wmax = w;
                                        
                                        
                                        Average_block(ima,ni,nj,nk,f,average,w,cols,rows,slices);
                                        totalweight = totalweight + w;
                                        
                                    }
                                    
                                    
                                }
                            }
                        }
                        
                    }
                    
                    if(wmax==0.0) wmax=1.0;
                    
                    Average_block(ima,i,j,k,f,average,wmax,cols,rows,slices);
                    
                    totalweight = totalweight + wmax;
                    
                    
                    if(totalweight != 0.0)
                        Value_block(Estimate,Label,i,j,k,f,average,totalweight,cols,rows,slices);
                }
                else
                {
                    wmax=1.0;
                    
                    Average_block(ref,i,j,k,f,average,wmax,cols,rows,slices);
                    
                    totalweight = totalweight + wmax;
                    
                    
                    if(totalweight != 0.0)
                        Value_block(Estimate,Label,i,j,k,f,average,totalweight,cols,rows,slices);
                }
                
                
            }
        }
    }
    
#ifdef _WIN32
_endthreadex(0);
return 0;
#else
pthread_exit(0);
return 0;
#endif

}


//void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
// First parameter - image
// Second parameter - v (search radius)
// Third parameter - f (patch radius)
// Fourth parameter - sigma (MAP)
// Fifth parameter - beta
// Sixth parameter - median filtered image
// Seventh parameter - mask
// Eighth parameter - rows
// Ninth parameter - columns
// Tenth parameter - slices
JNIEXPORT jfloatArray JNICALL Java_CANDLEWrapper_adaptONLMCollaborative(JNIEnv * env, jobject jobj, jfloatArray jima, jint jv, jint jf, jfloatArray jsigma, jfloat jbeta, jfloatArray jref,
        jfloatArray jmask, jint jrows, jint jcols, jint jnumSlices)
{
    /*Declarations*/
    //mxArray *xData;
    jfloat *xData; // Stores the output
    
    float *ima, *fima,*average, *sigma,*ref, *mask;
    
    //mxArray *Mxmeans, *Mxvariances, *MxEstimate, *MxLabel;
    jfloat *Mxmeans, *Mxvariances, *MxEstimate, *MxLabel;
    
    float *means, *variances, *Estimate;
    
    // Is unsigned short supported?
    unsigned short *Label;
    
    //mxArray *pv;
    jfloat *pv;
    float w,totalweight,wmax,d,mean,var,t1,t2,epsilon,label,estimate;
    float beta;
    int Ndims,i,j,k,ii,jj,kk,ni,nj,nk,v,f,ndim,indice,init,Nthreads,ini,fin;
    //const int  *dims;
    
    // Defined earlier - structure
    myargument *ThreadArgs;
    
#ifdef _WIN32
HANDLE *ThreadList; // Handles to the worker threads
#else
pthread_t * ThreadList;
#endif

/*Copy input pointer x*/


/*Get matrix x*/
//ima = (float*)mxGetPr(prhs[0]);
ima = (*env)->GetFloatArrayElements(env, jima, 0);

//ndim = mxGetNumberOfDimensions(prhs[0]);
//dims= mxGetDimensions(prhs[0]);

/*Copy input parameters*/

/*Get the Integer*/
//v = (int)(mxGetScalar(prhs[1]));
v = (int)jv;

//f = (int)(mxGetScalar(prhs[2]));
f = (int)jf;

//sigma =  (float*)mxGetPr(prhs[3]);
sigma = (*env)->GetFloatArrayElements(env, jsigma, 0);

//beta = (float)(mxGetScalar(prhs[4]));
beta = (float)jbeta;

//ref =  (float*)mxGetPr(prhs[5]);
// Corrected
ref = (*env)->GetFloatArrayElements(env, jref, 0);

//mask =  (float*)mxGetPr(prhs[6]);
mask = (*env)->GetFloatArrayElements(env, jmask, 0);

// ndim is undefined
Ndims = pow((2*f+1),ndim);

int rows = (int)jrows;
int cols = (int)jcols;
int slices = (int)jslices;
int totalElements = rows*cols*slices;

/*Allocate memory and assign output pointer*/

//plhs[0] = mxCreateNumericArray(ndim,dims,mxGetClassID(prhs[0]), mxREAL);
// What is happening here?
xData = env->NewFloatArray(totalElements);

//Mxmeans = mxCreateNumericArray(ndim,dims,mxGetClassID(prhs[0]), mxREAL);
Mxmeans = env->NewFloatArray(totalElements);

//Mxvariances = mxCreateNumericArray(ndim,dims,mxGetClassID(prhs[0]), mxREAL);
Mxvariances = env->NewFloatArray(totalElements);

//MxEstimate = mxCreateNumericArray(ndim,dims,mxGetClassID(prhs[0]), mxREAL);
MxEstimate = env->NewFloatArray(totalElements);

//MxLabel = mxCreateNumericArray(ndim,dims,mxUINT16_CLASS, mxREAL);
MxLabel = env->NewFloatArray(totalElements);

average=(float*)malloc(Ndims*sizeof(float));

/*Get a pointer to the data space in our newly allocated memory*/
//fima = (float*)mxGetPr(plhs[0]);
fima = (*env)->GetFloatArrayElements(env, xData, 0);

//means = (float*)mxGetPr(Mxmeans);
means = (*env)->GetFloatArrayElements(env, Mxmeans, 0);

//variances = (float*)mxGetPr(Mxvariances);
variances = (*env)->GetFloatArrayElements(env, Mxvariances, 0);
//Estimate = (float*)mxGetPr(MxEstimate);
// Corrected
Estimate = (*env)->GetFloatArrayElements(env, MxEstimate, 0);

//Label = (unsigned short*)mxGetPr(MxLabel);
Label = (*env)->GetIntArrayElements(env, MxLabel, 0);


// for (i = 0; i < dims[2] *dims[1] * dims[0];i++)
for (i = 0; i <totalElements;i++)
{
    Estimate[i] = 0.0;
    Label[i] = 0.0;
    fima[i] = 0.0;
}

max=0;
// for(k=0;k<dims[2];k++) // For all slices
for(k=0;k<slices;k++) // For all slices
{
    //for(i=0;i<dims[1];i++)
    for(i=0;i<cols;i++)
    {
        //for(j=0;j<dims[0];j++)
        for(j=0;j<rows;j++)
        {
            mean=0;
            indice=0;
            for(ii=-1;ii<=1;ii++)
            {
                for(jj=-1;jj<=1;jj++)
                {
                    for(kk=-1;kk<=1;kk++)
                    {
                        ni=i+ii;
                        nj=j+jj;
                        nk=k+kk;
                        
                        if(ni<0) ni=-ni;
                        if(nj<0) nj=-nj;
                        if(nk<0) nk=-nk;
                        if(ni>=cols) ni=2*cols-ni-1;
                        if(nj>=rows) nj=2*rows-nj-1;
                        if(nk>=slices) nk=2*slices-nk-1;
                        
                        
                        // mean = mean + ref[nk*(dims[0]*dims[1])+(ni*dims[0])+nj];
                        mean = mean + ref[nk*(rows*cols)+(nj*cols)+ni];
                        indice=indice+1;
                        
                    }
                }
            }
            mean=mean/indice;
            //means[k*(dims[0]*dims[1])+(i*dims[0])+j]=mean;
            means[k*(rows*cols)+(j*cols)+i]=mean;
            if (mean > max) max =mean;
        }
    }
}

//for(k=0;k<dims[2];k++)
for(k=0;k<slices;k++)
{
    //for(i=0;i<dims[1];i++)
    for(i=0;i<cols;i++)
    {
        //for(j=0;j<dims[0];j++)
        for(j=0;j<rows;j++)
        {
            var=0;
            indice=0;
            for(ii=-1;ii<=1;ii++)
            {
                for(jj=-1;jj<=1;jj++)
                {
                    for(kk=-1;kk<=1;kk++)
                    {
                        ni=i+ii;
                        nj=j+jj;
                        nk=k+kk;
                        if(ni>=0 && nj>=0 && nk>0 && ni<cols && nj<rows && nk<dims[2])
                        {
                            //var = var + (ref[nk*(dims[0]*dims[1])+(ni*dims[0])+nj]-means[k*(dims[0]*dims[1])+(i*dims[0])+j])*(ref[nk*(dims[0]*dims[1])+(ni*dims[0])+nj]-means[k*(dims[0]*dims[1])+(i*dims[0])+j]);
                            // Corrected
                            var = var + (ref[nk*(rows*cols)+(nj*cols)+ni]-means[k*(rows*cols)+(j*cols)+i])*(ref[nk*(rows*cols)+(nj*cols)+ni]-means[k*(rows*cols)+(j*cols)+i]);
                            indice=indice+1;
                        }
                    }
                }
            }
            var=var/(indice-1);
            //variances[k*(dims[0]*dims[1])+(i*dims[0])+j]=var;
            variances[k*(rows*cols)+(j*cols)+i]=var;
        }
    }
}

Nthreads=8;


#ifdef _WIN32

// Reserve room for handles of threads in ThreadList
ThreadList = (HANDLE*)malloc(Nthreads* sizeof( HANDLE ));
ThreadArgs = (myargument*) malloc( Nthreads*sizeof(myargument));

for (i=0; i<Nthreads; i++)
{
    /* Make Thread Structure*/
    //ini=(i*dims[2])/Nthreads;
    ini=(i*slices)/Nthreads;
    //        fin=((i+1)*dims[2])/Nthreads;
    fin=((i+1)*slices)/Nthreads;
    //        ThreadArgs[i].cols=dims[0];
    ThreadArgs[i].cols=cols;
    //        ThreadArgs[i].rows=dims[1];
    ThreadArgs[i].rows=rows;
    //        ThreadArgs[i].slices=dims[2];
    ThreadArgs[i].slices=slices;
    ThreadArgs[i].in_image=ima;
    ThreadArgs[i].ref_image=ref;
    ThreadArgs[i].mask_image=mask;
    ThreadArgs[i].var_image=variances;
    ThreadArgs[i].means_image=means;
    ThreadArgs[i].estimate=Estimate;
    ThreadArgs[i].label=Label;
    ThreadArgs[i].ini=ini;
    ThreadArgs[i].fin=fin;
    ThreadArgs[i].radioB=v;
    ThreadArgs[i].radioS=f;
    ThreadArgs[i].sigma=sigma;
    ThreadArgs[i].beta=beta;
    ThreadList[i] = (HANDLE)_beginthreadex( NULL, 0, &ThreadFunc, &ThreadArgs[i] , 0, NULL );
}

for (i=0; i<Nthreads; i++) { WaitForSingleObject(ThreadList[i], INFINITE); }
for (i=0; i<Nthreads; i++) { CloseHandle( ThreadList[i] ); }

#else

/* Reserve room for handles of threads in ThreadList*/
ThreadList = (pthread_t *) calloc(Nthreads,sizeof(pthread_t));
ThreadArgs = (myargument*) calloc( Nthreads,sizeof(myargument));

for (i=0; i<Nthreads; i++)
{
    /* Make Thread Structure*/
    //        ini=(i*dims[2])/Nthreads;
    //        fin=((i+1)*dims[2])/Nthreads;
    ini=(i*slices)/Nthreads;
    fin=((i+1)*slices)/Nthreads;
    //        ThreadArgs[i].cols=dims[0];
    ThreadArgs[i].cols=cols;
    //        ThreadArgs[i].rows=dims[1];
    ThreadArgs[i].rows=rows;
    //        ThreadArgs[i].slices=dims[2];
    ThreadArgs[i].slices=slices;
    ThreadArgs[i].in_image=ima;
    ThreadArgs[i].ref_image=ref;
    ThreadArgs[i].mask_image=mask;
    ThreadArgs[i].var_image=variances;
    ThreadArgs[i].means_image=means;
    ThreadArgs[i].estimate=Estimate;
    ThreadArgs[i].label=Label;
    ThreadArgs[i].ini=ini;
    ThreadArgs[i].fin=fin;
    ThreadArgs[i].radioB=v;
    ThreadArgs[i].radioS=f;
    ThreadArgs[i].sigma=sigma;
    ThreadArgs[i].beta=beta;
}

for (i=0; i<Nthreads; i++)
{
    
    if(pthread_create(&ThreadList[i], NULL, ThreadFunc,&ThreadArgs[i]))
    {
        printf("Threads cannot be created\n");
        exit(1);
    }
    
}



for (i=0; i<Nthreads; i++)
{
    pthread_join(ThreadList[i],NULL);
}

#endif

free(ThreadArgs);
free(ThreadList);


label = 0.0;
estimate = 0.0;


/* Aggregation of the estimators (i.e. means computation) */
//for (k = 0; k < dims[2]; k++ )
for (k = 0; k < slices; k++ )
{
    //        for (i = 0; i < dims[1]; i++ )
    for (i = 0; i < cols; i++ )
    {
        //        for (j = 0; j < dims[0]; j++ )
        for (j = 0; j < rows; j++ )
        {
            //label = Label[k*(rows*dims[1])+(i*dims[0])+j];
            label = Label[k*(rows*cols)+(j*cols)+i];
            if (label == 0.0)
            {
                //                    fima[k*(dims[0]*dims[1])+(i*dims[1])+j] = ima[k*(dims[0]*dims[1])+(i*dims[0])+j];
                // Corrected ?
                fima[k*(rows*cols)+(j*rows)+i] = ima[k*(rows*cols)+(j*cols)+i];
                
            }
            else
            {
                //estimate = Estimate[k*(dims[0]*dims[1])+(i*dims[0])+j];
                estimate = Estimate[k*(rows*cols)+(j*cols)+i];
                estimate = (estimate/label);
                //fima[k*(dims[0]*dims[1])+(i*dims[0])+j]=estimate;
                fima[k*(rows*cols)+(j*cols)+i]=estimate;
                
            }
        }
    }
}

    
(*env)->ReleaseFloatArrayElements(env, jima, ima, 0);
(*env)->ReleaseFloatArrayElements(env, jsigma, sigma, 0);
(*env)->ReleaseFloatArrayElements(env, jref, ref, 0);
(*env)->ReleaseFloatArrayElements(env, jmask, mask, 0);
    
(*env)->ReleaseFloatArrayElements(env, Mxmeans, means, 0);
(*env)->ReleaseFloatArrayElements(env, Mxvariances, variances, 0);
(*env)->ReleaseFloatArrayElements(env, MxEstimate, Estimate, 0);
(*env)->ReleaseFloatArrayElements(env, MxLabel, Label, 0);

    
    
// Do not understand
return xData;

}




//public class CANDLEWrapper {
//   native void adaptONLMCollaborative(float[] jima, int v, int f,
//                                        float[] sigma, float beta,
//                                        float[] ref, float[] mask,
//                                        int rows, int cols, int numSlices); /* (1) */
//   static {
//       System.loadLibrary("adapt_onlm_collaborative"); /* (2) */
//   }
//  static public void main(String argv[]) {
//      CANDLEWrapper candle = new CANDLEWrapper();
//        // Do parameter set up here
//        // Read in the images - convert to 1D arrays
//        // Call method
//      float[] fima = candle.adaptONLMCollaborative(jima, v, f, sigma, beta,
//								        ref, mask, rows, cols, numSlices);
//        // Convert to a 3D stack after for fima
//
//  }
//}

/*
