package algorithms;

import com.github.fommil.netlib.BLAS;
import java.lang.reflect.*;
import java.util.Random;
import junit.framework.TestCase;

/**
 * from http://jeshua.me/blog/NetlibJavaJNI
 * Code to do some simple speed comparisons between netlib java JNI and Java implementations.
 *  by author Jeshua Bratman 
 * 
 * below are edits to work with version 1.1.2
 */
@SuppressWarnings({"rawtypes"})
public class NetlibTest extends TestCase {
    
  public void test() {
    
      //com.github.fommil.netlib.ARPACK
      
    try{
      //get java blas object and make it accessible
      //Class javaBlasClass = Class.forName("org.netlib.blas.JBLAS");
      Class javaBlasClass = com.github.fommil.netlib.BLAS.class;
      Field javaBlas = javaBlasClass.getDeclaredField("INSTANCE");
      Field jInstance = javaBlas.getClass().getDeclaredField("modifiers");
      jInstance.setAccessible(true);
      jInstance.setInt(javaBlas,javaBlas.getModifiers() & ~Modifier.FINAL);
      javaBlas.setAccessible(true);

      //get native blas object and make it accessible
      //Class nativeBlasClass = Class.forName("org.netlib.blas.NativeBLAS");
      Class nativeBlasClass = com.github.fommil.netlib.NativeRefBLAS.class;
      //Class nativeBlasClass = com.github.fommil.netlib.NativeSystemBLAS.class;
      Field nativeBlas = nativeBlasClass.getSuperclass().getSuperclass().
          getDeclaredField("INSTANCE");
      Field nInstance = nativeBlas.getClass().getDeclaredField("modifiers");
      nInstance.setAccessible(true);
      nInstance.setInt(nativeBlas,nativeBlas.getModifiers() & ~Modifier.FINAL);
      nativeBlas.setAccessible(true);

      //get blas current object and make it accessible
      //Field blasCurrent = Class.forName("org.netlib.blas.BLAS").getDeclaredField("current");
      //Field blasCurrent = com.github.fommil.netlib.BLAS.class.getDeclaredField("current");
      //Field bInstance = blasCurrent.getClass().getDeclaredField("modifiers");
      //bInstance.setAccessible(true);
      //bInstance.setInt(blasCurrent, blasCurrent.getModifiers() & ~Modifier.FINAL);
      //blasCurrent.setAccessible(true);
      Field blasCurrent = nativeBlas;

      //multiply various sized matrices
      int[] szs = {100, 200, 500, 800, 1000, 2000, 3000, 4000};
      for(int k = 0;k<szs.length;k++){
        int sz = szs[k];

        double[] mat1 = new double[sz*sz];
        double[] mat2 = new double[sz*sz];
        double[] matOut = new double[sz*sz];
        Random rand = new Random();
        for(int i=0;i<sz*sz;i++){
          mat1[i] = rand.nextDouble();
          mat2[i] = rand.nextDouble();
          matOut[i] = 0;
        }				
        //SET TO JBLAS
        blasCurrent.set(null, javaBlas.get(null));

        long time1 = System.nanoTime();//bean.getCurrentThreadUserTime();	
        BLAS.getInstance().dgemm("N", "N",
            sz, sz, sz,
            1.0, mat1, sz, mat2, sz,
            0.0, matOut, sz);		
        double elapsed = (System.nanoTime() - time1)/1000000000d;
        System.out.printf("JBLAS:\t\t %.4f seconds to multiply %dx%d matrices.\n",elapsed,sz,sz);							

        //SET TO NativeBLAS
        blasCurrent.set(null, nativeBlas.get(null));

        time1 = System.nanoTime();
        BLAS.getInstance().dgemm("N", "N",
            sz, sz, sz,
            1.0, mat1, sz, mat2, sz,
            0.0, matOut, sz);		
        elapsed = (System.nanoTime() - time1)/1000000000d;
        System.out.printf("NativeBLAS:\t %.4f seconds to multiply %dx%d matrices.\n",elapsed,sz,sz);
      }
    }catch(Exception e){e.printStackTrace();}
  }
  
  /*
  [junit] JBLAS:               0.0735 seconds to multiply 100x100 matrices.
    [junit] NativeBLAS:  0.0009 seconds to multiply 100x100 matrices.
    [junit] JBLAS:               0.0024 seconds to multiply 200x200 matrices.
    [junit] NativeBLAS:  0.0028 seconds to multiply 200x200 matrices.
    [junit] JBLAS:               0.0310 seconds to multiply 500x500 matrices.
    [junit] NativeBLAS:  0.0305 seconds to multiply 500x500 matrices.
    [junit] JBLAS:               0.1825 seconds to multiply 800x800 matrices.
    [junit] NativeBLAS:  0.1577 seconds to multiply 800x800 matrices.
    [junit] JBLAS:               0.3085 seconds to multiply 1000x1000 matrices.
    [junit] NativeBLAS:  0.3614 seconds to multiply 1000x1000 matrices.
    [junit] JBLAS:               2.9021 seconds to multiply 2000x2000 matrices.
    [junit] NativeBLAS:  2.8109 seconds to multiply 2000x2000 matrices.
    [junit] JBLAS:               8.8778 seconds to multiply 3000x3000 matrices.
    [junit] NativeBLAS:  9.0282 seconds to multiply 3000x3000 matrices.
    [junit] JBLAS:               23.2959 seconds to multiply 4000x4000 matrices.
    [junit] NativeBLAS:  25.0183 seconds to multiply 4000x4000 matrices.
  */
}
