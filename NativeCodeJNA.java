
import com.sun.jna.Library;
import com.sun.jna.Memory;
import com.sun.jna.Native;
import com.sun.jna.Pointer;
import com.sun.jna.ptr.PointerByReference;
import ij.IJ;


public class NativeCodeJNA {
	
	// Interface that declares native methods of the ONLM C program
	public interface ONLMInterface extends Library {
        public void ONLM(Pointer p1, int x1, int x2, Pointer p2, float x3, Pointer p3, int x4, int x5, int x6, PointerByReference valsRef);
        public void cleanup(Pointer s);
    }
	
	// Interface that declares native methods of the NoiseEstimation C program
	public interface NoiseEstimationInterface extends Library {
        public void estimate(Pointer p1, int x1, int x2, int x3, int x4, PointerByReference valsRef);
    }
	
	// Method that calls the two (NoiseEstimation and ONLM) native programs.
	static public float[] NativeCall(float[] img , float[] medfiltimg, int searchrad , int patchrad, float beta, int x, int y, int z ){
		
		
		// Load Native Library: libNoiseEstimation.dylib
    	NoiseEstimationInterface noisest = (NoiseEstimationInterface) Native.loadLibrary("NoiseEstimation", NoiseEstimationInterface.class);
    	
    	int size = x*y*z;
    	
		// Allocate memory for native pointers
    	Pointer ima = new Memory(size * Native.getNativeSize(Float.TYPE));
    	Pointer ref = new Memory(size * Native.getNativeSize(Float.TYPE));
    	
    	
    	for (int dloop=0; dloop<(size); dloop++) {
    		
    		ima.setFloat(dloop * Native.getNativeSize(Float.TYPE), img[dloop] );
    		ref.setFloat(dloop * Native.getNativeSize(Float.TYPE), medfiltimg[dloop] );
    		
    	}
    	
    	// Pointer to a pointer : This lets us read the output array from the first native method
    	PointerByReference PMAP = new PointerByReference();
        
        IJ.log("Wavelet-based local noise estimation");
    	
        long startTime = System.currentTimeMillis();
        
    	// Call to the first native method
        noisest.estimate(ima, x, y, z, 2*searchrad, PMAP);
    	
        long endTime = System.currentTimeMillis();
        
        IJ.log("Elapsed time: " + ((endTime - startTime)/1000) );
        
    	Pointer MAP = PMAP.getValue();
        
        
    	// Load Native Library: libONLM_Mod.dylib
    	ONLMInterface onlm = (ONLMInterface) Native.loadLibrary("ONLM_Mod", ONLMInterface.class);
        
    	
    	// Pointer to a pointer : This lets us read the output array from the second native method
    	PointerByReference Pfima = new PointerByReference();
        
        IJ.log("Denoising: 3D Optimized Non-local Means Filter");
        
        startTime = System.currentTimeMillis();
		
    	// Call to the second native method
    	onlm.ONLM(ima, searchrad, patchrad, MAP, beta, ref, x, y, z, Pfima);
    	
        endTime = System.currentTimeMillis();
        
        IJ.log("Elapsed time: " + ((endTime - startTime)/1000) );
    	
    	Pointer fima = Pfima.getValue();
    	
    	// Get final result
    	for (int dloop=0; dloop<size; dloop++) {
    		
    		medfiltimg[dloop] = fima.getFloat(dloop * Native.getNativeSize(Float.TYPE));
    			
    	}
    	
    	// free memory
    	onlm.cleanup(fima);
    	
    	return medfiltimg;
    	
		
	}
	

}
