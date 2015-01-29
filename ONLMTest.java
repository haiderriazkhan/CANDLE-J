package JNApackage;



import com.sun.jna.Memory;
import com.sun.jna.Native;
import com.sun.jna.Library;
import com.sun.jna.Pointer;
import com.sun.jna.ptr.PointerByReference;

public class ONLMTest {
	
	
	public interface CTest extends Library {
        public float ONLM(Pointer p1, int x1, int x2, Pointer p2, float x3, Pointer p3, Pointer p4, int x4, int x5, int x6, PointerByReference valsRef);
        public void cleanup(Pointer s);
        public void Hello();
    }
	
	static public float[] ONLMInputs(float[] img , int searchrad , int patchrad, float[] MAP, float beta , float[] fimgMed, float[] mask , int x, int y, int z ) {
		System.setProperty("jna.library.path", "/Users/haiderriaz/Desktop/JNA-C/ONLMTest");
    	CTest ctest = (CTest) Native.loadLibrary("ONLM_Mod", CTest.class);
    	int size = x*y*z;
    	
    	Pointer ima = new Memory(size * Native.getNativeSize(Float.TYPE));
    	Pointer sigma = new Memory(size * Native.getNativeSize(Float.TYPE));
    	Pointer ref = new Memory(size * Native.getNativeSize(Float.TYPE));
    	Pointer maskk = new Memory(size * Native.getNativeSize(Float.TYPE));
    	
    	for (int dloop=0; dloop<(size); dloop++) {
    		
    		ima.setFloat(dloop * Native.getNativeSize(Float.TYPE), img[dloop] );
    		
    		sigma.setFloat(dloop * Native.getNativeSize(Float.TYPE), MAP[dloop] );
    		ref.setFloat(dloop * Native.getNativeSize(Float.TYPE), fimgMed[dloop] );
    		maskk.setFloat(dloop * Native.getNativeSize(Float.TYPE), mask[dloop] );
    		
    	}
    	PointerByReference Pfima = new PointerByReference();
    	
    	ctest.ONLM(ima, searchrad, patchrad, sigma, beta, ref, maskk, x, y, z, Pfima);
    	
    	Pointer fima = Pfima.getValue();
    	
    	for (int dloop=0; dloop<size; dloop++) {
    		
    		MAP[dloop] = fima.getFloat(dloop * Native.getNativeSize(Float.TYPE));
    		
    		
    	}
    	
    	ctest.cleanup(fima);
    	return MAP;
		
	}
	
	
	
	
	
	
	static public void main(String argv[]) {
		System.setProperty("jna.library.path", "/Users/haiderriaz/Desktop/JNA-C/ONLMTest");
    	CTest ctest = (CTest) Native.loadLibrary("ONLM_Mod", CTest.class);
    	int size = 1000;
    	Pointer ima = new Memory(size * Native.getNativeSize(Float.TYPE));
    	Pointer sigma = new Memory(size * Native.getNativeSize(Float.TYPE));
    	Pointer ref = new Memory(size * Native.getNativeSize(Float.TYPE));
    	Pointer maskk = new Memory(size * Native.getNativeSize(Float.TYPE));
    	
    	for (int dloop=0; dloop<size; dloop++) {
    		
    		ima.setFloat(dloop * Native.getNativeSize(Float.TYPE), (float)dloop );
    		sigma.setFloat(dloop * Native.getNativeSize(Float.TYPE), (float)dloop );
    		ref.setFloat(dloop * Native.getNativeSize(Float.TYPE), (float)dloop );
    		maskk.setFloat(dloop * Native.getNativeSize(Float.TYPE), (float)dloop );
    		
    	}
    	
    	
    	
    	PointerByReference Pfima = new PointerByReference();
    	System.out.println("Here1");
    	ctest.ONLM(ima, 6, 6, sigma, (float)1, ref, maskk, 10, 10, 10, Pfima);
    	System.out.println("Here2");
    	Pointer fima = Pfima.getValue();
    	
    	
    	for (int dloop=0; dloop<size; dloop++) {
    		float print = fima.getFloat(dloop * Native.getNativeSize(Float.TYPE));
    		System.out.println("\tvalue: " + print);
    		
    	}
    	
    	ctest.cleanup(fima);
    	
    }
	

}
