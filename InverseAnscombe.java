// Haider Khan - haiderriazkhan@hotmail.com
//
// Optimal inverse Anscombe Transform (AT) 
//
// References
// [1] M. Makitalo and A. Foi, "On the inversion of the Anscombe transformation in low-count Poisson image denoising", 
// Proc. Int. Workshop on Local and Non-Local Approx. in Image Process., LNLA 2009, Tuusula, Finland, pp. 26-32, August 2009.
// [2] M. Makitalo and A. Foi, "Optimal inversion of the Anscombe transformation in low-count Poisson image denoising",
// TIP 2010.
//
// -------------------------------------------------------------------------------





import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import ij.IJ;


public class InverseAnscombe {
	
// Load pre-computed correspondence table	
public static Object[] loadOVSTable(){
	
	double[] firstColumn = new double[30000];
	double[] secondColumn = new double[30000];
	int index = 0;

	try {
	      BufferedReader in = new BufferedReader(new FileReader("plugins/CANDLE-J/Table.dat"));
	      String str ;

	      while ((str = in.readLine()) != null) {
	             String[] values = str.split(",");
	             firstColumn[index] = Double.parseDouble(values[0]);
	             secondColumn[index++] = Double.parseDouble(values[1]);           
	      }

	       in.close();
	    } 
	catch (IOException e) 
	{
		IJ.log("Error: Could not load pre-computed correspondence table.");
	}
	
	return new Object[]{firstColumn, secondColumn};
	
}


// Binary Search algorithm for a sorted array of doubles and a floating point search value
public static int binary_search(double[] x, float value){
	
	Integer middle = null;
	int first = 0;
	int last = x.length - 1;
	
	while(first <= last){
		
		middle = (first + last) / 2 ;
		
		if(x[middle] < value){
			
			first = middle + 1;
		}else if(x[middle] == value){
			
			break;
		}else{
			
			last = middle - 1;
		}
		
		if( first > last){
			middle = last;
		    break;
		}
		
	}
	
	return middle;
	
}

// The method implements Matlab's 1D interpolation function
public static float[] interp1(double[] x, double[] y, float[]xq){
	
	double x1,x2,y1,y2,m,interp;
	int index;
	int counter = 0;
	
	float [] vq = new float[xq.length];
	
	for(float pt : xq){
		
		if(Float.isNaN(pt)){
			
			vq[counter++] = Float.NaN;
			continue;
		}
		
		if(pt < x[0]){
			x1 = x[0];
		    x2 = x[1];
		    y1 = y[0];
		    y2 = y[1];
		}else if(pt > x[x.length-1]){
			x1 = x[x.length-2];
		    x2 = x[x.length-1];
		    y1 = y[x.length-2];
		    y2 = y[x.length-1];
			
		}else{
			
			index = binary_search(x, pt);
			
			x1 = x[index];
		    x2 = x[index+1];
		    y1 = y[index];
		    y2 = y[index+1];
			
		}
		
		m = (y2 - y1) / (x2 - x1);
	    interp = y1 + m*(pt - x1);
        
        if(interp < 0){
            vq[counter++] = 0;
            continue;
        }
        
	    vq[counter++] = (float)interp;
		
	}
	
	return vq;
}
	
	
// Parent method that implements the Inverse Anscombe Transform	
public static float[] OVST(float[] fimg){
		
		int len = fimg.length;
		
		float [] asympinvfimg = new float[len];
		
		// Estimation of the asympotic Inverse Anscombe Transform (Classical one)
		for(int i = 0; i < len; i++){
			asympinvfimg[i] =  (float) ( Math.pow( (fimg[i] / 2) , 2 ) - (0.125) );	
		}
		
		Object[] tuple = loadOVSTable();
		
		double[] Efz = (double[]) tuple[0];
		double[] Ez = (double[]) tuple[1];
		
		// Estimation of the optimal Inverse Anscombe Transform (see paper [2]).
		float[] optimalinvfimg = interp1(Efz, Ez , fimg);
		  
		double max_Efz = 100.0025;
		double min_const = 2*Math.sqrt(3/8) ;
		
		// Areas under and above the validity range are filled 
		for(int i=0; i < len; i++){
			
			if( fimg[i] > max_Efz ){
				optimalinvfimg[i] = asympinvfimg[i];
			}else if( fimg[i] < min_const ){
				optimalinvfimg[i] = 0;
			}
			
		}
		
		return optimalinvfimg;
    
		
	}

}


