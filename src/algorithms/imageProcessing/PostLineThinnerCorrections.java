package algorithms.imageProcessing;

import algorithms.util.PairInt;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.Set;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class PostLineThinnerCorrections {
    
    private Logger log = Logger.getLogger(this.getClass().getName());

    public void correctForArtifacts(GreyscaleImage input) {
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        Set<PairInt> points = imageProcessor.readNonZeroPixels(input);
        
        //TODO: reduce the number of patterns here if possible
        // and make sure that true corners aren't drastically reduced to less
        // usable smaller corners 
        
        //correctForHoleArtifacts3(input);
        //correctForHoleArtifacts2(input);

        /*
        correctForHoleArtifacts1(points);
        
        correctForHoleArtifacts1_2(points);
                 
        correctForHoleArtifacts1_3(points);
        
        correctForHoleArtifacts1_4(points);
        
        correctForZigZag0(points);
        
        correctForZigZag0Alt(points);
                
        correctForZigZag1(points);
     
        correctForZigZag2(points);
        
        correctForZigZag1Alt(points);
        
        correctForZigZag3(points);
                
        correctForZigZag5(points);
        
        correctForZigZag6(points);
        
        correctForWs(points);
        
        // TODO: revisit, not sure this is always an artifact:
        correctForLine0(points);
        
        // better edge extraction at the expense of unsharpening true corners:
        correctForLs(input);
        correctForLs2(points);
        
        correctForZigZag1(points);
        
        correctForSpurs(points);
        
        correctForZigZag7(points);
*/
        imageProcessor.writeAsBinaryToImage(input, points);
    }
    
    private void correctForZigZag0(GreyscaleImage input) {
       
        /*
        looking for pattern
       
           0  0         2
           0  #  #      1
           #* #< 0      0
        #     0  0     -1
        
       -1  0  1  2
        
        and removing the topmost left #'s
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(0, -1));
        zeroes.add(new PairInt(0, -2));
        zeroes.add(new PairInt(1, 1));
        zeroes.add(new PairInt(1, -2));
        zeroes.add(new PairInt(2, 0));
        zeroes.add(new PairInt(2, 1));
        
        ones.add(new PairInt(-1, 1));
        ones.add(new PairInt(1, -1));
        ones.add(new PairInt(1, 0));
        ones.add(new PairInt(2, -1));
        
        changeToZeroes.add(new PairInt(1, 0));
        
        int startValue = 1;
            
        replacePattern(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
        rotate90ThreeTimes(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
    }
    
    private void correctForZigZag0Alt(GreyscaleImage input) {
       
        /*
        keep
       
           0  0  #  0   2
        0  0  #  0  0   1
        0  #* #< 0  0   0
        0  #  0  0     -1
        0  #  0        -2
        
       -1  0  1  2  3       
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
     
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-1, 2));
        zeroes.add(new PairInt(-1, 1));
        zeroes.add(new PairInt(-1, 0));
        zeroes.add(new PairInt(-1, -1));
        zeroes.add(new PairInt(0, -1));
        zeroes.add(new PairInt(0, -2));
        zeroes.add(new PairInt(1, 1));
        zeroes.add(new PairInt(1, 2));
        zeroes.add(new PairInt(1, -2));
        zeroes.add(new PairInt(2, 1));
        zeroes.add(new PairInt(2, 0));
        zeroes.add(new PairInt(2, -1));
        zeroes.add(new PairInt(3, 0));
        zeroes.add(new PairInt(3, -1));
        zeroes.add(new PairInt(3, -2));
       
        ones.add(new PairInt(0, 2));
        ones.add(new PairInt(0, 1));
        ones.add(new PairInt(1, 0));
        ones.add(new PairInt(1, -1));
        ones.add(new PairInt(2, -2));
        
        changeToZeroes.add(new PairInt(1, 0));
        
        int startValue = 1;
            
        replacePattern(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
        rotate90ThreeTimes(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
    }
    
    private void correctForZigZag5(GreyscaleImage input) {
       
        /*
        keep
            
        0  #  0        2
        0  #  0  0     1
        0  #  #*<0     0
        0  0  #  0    -1
              0  #    -2
        
       -2 -1  0  1         
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
     
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-2, 1));
        zeroes.add(new PairInt(-2, 0));
        zeroes.add(new PairInt(-2, -1));
        zeroes.add(new PairInt(-2, -2));
        zeroes.add(new PairInt(-1, 1));
        zeroes.add(new PairInt(0, 2));
        zeroes.add(new PairInt(0, -1));
        zeroes.add(new PairInt(0, -2));
        zeroes.add(new PairInt(1, 1));
        zeroes.add(new PairInt(1, 0));
        zeroes.add(new PairInt(1, -1));
        
        ones.add(new PairInt(-1, 0));
        ones.add(new PairInt(-1, -1));
        ones.add(new PairInt(-1, -2));
        ones.add(new PairInt(0, 1));
        ones.add(new PairInt(1, 2));
        
        changeToZeroes.add(new PairInt(0, 0));
        
        int startValue = 1;
            
        replacePattern(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
        rotate90ThreeTimes(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
    }
    
    private void correctForZigZag6(GreyscaleImage input) {
       
        /*       
           0  0  0  0   2
           0  #< #  #   1
        #  #* #  0      0
           0  0  0     -1
        
       -1  0  1  2        
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(0, 1));
        zeroes.add(new PairInt(0, -1));
        zeroes.add(new PairInt(0, -2));
        zeroes.add(new PairInt(1, 1));
        zeroes.add(new PairInt(1, -2));
        zeroes.add(new PairInt(2, 1));
        zeroes.add(new PairInt(2, 0));
        zeroes.add(new PairInt(2, -2));
        zeroes.add(new PairInt(3, -2));
       
        ones.add(new PairInt(-1, 0));
        ones.add(new PairInt(1, 0));
        ones.add(new PairInt(1, -1));
        ones.add(new PairInt(2, -1));
        ones.add(new PairInt(3, -1));
        
        changeToZeroes.add(new PairInt(1, -1));
        
        int startValue = 1;
            
        replacePattern(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
        rotate90ThreeTimes(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
    }
    
    /**
     * using this one unfortunately reduces the sharpness of real corners,
     * but it allows clear point ordering when extracting edges from the
     * line thinned image.
     * 
     * @param input 
     */
    private void correctForZigZag7(GreyscaleImage input) {
       
        /*  
           #
        0  #  0         1
        #  #*<0         0
        #  0  0        -1
        #              -2
        
       -1  0  1  2 
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
    
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(0, 1));
        zeroes.add(new PairInt(1, 1));
        zeroes.add(new PairInt(1, 0));
        zeroes.add(new PairInt(1, -1));
        zeroes.add(new PairInt(-1, -1));
        
        ones.add(new PairInt(0, -1));
        ones.add(new PairInt(-1, 0));
        ones.add(new PairInt(-1, 1));
        ones.add(new PairInt(-1, 2));
        ones.add(new PairInt(0, -2));
        
        changeToZeroes.add(new PairInt(0, 0));
        
        int startValue = 1;
            
        replacePattern(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
        rotate90ThreeTimes(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
    }
    
    /**
     * using this one unfortunately reduces the sharpness of real corners,
     * but it allows clear point ordering when extracting edges from the
     * line thinned image.
     * 
     * @param input 
     */
    private void correctForLs(GreyscaleImage input) {
       
        /*  
        
        0  #  0         1
        0  #* 0         0
        0  #< #        -1
        0  0  0        -2
        
       -1  0  1  2 
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
    
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-1, 2));
        zeroes.add(new PairInt(-1, 1));
        zeroes.add(new PairInt(-1, 0));
        zeroes.add(new PairInt(-1, -1));
        zeroes.add(new PairInt(0, 2));
        zeroes.add(new PairInt(1, 2));
        zeroes.add(new PairInt(1, 0));
        zeroes.add(new PairInt(1, -1));
        
        ones.add(new PairInt(0, 1));
        ones.add(new PairInt(0, -1));
        ones.add(new PairInt(1, 1));
        
        changeToZeroes.add(new PairInt(0, 1));
        
        int startValue = 1;
            
        replacePattern(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
        rotate90ThreeTimes(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
    }
    
    private void correctForLs2(GreyscaleImage input) {
       
        /*  
        0  #  0  0      1
        0  #*<#  #      0
        0  0  0        -1
        
       -1  0  1  2 
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
   
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-1, -1));
        zeroes.add(new PairInt(-1, 0));
        zeroes.add(new PairInt(-1, 1));
        zeroes.add(new PairInt(0, 1));
        zeroes.add(new PairInt(1, 1));
        zeroes.add(new PairInt(1, -1));
        zeroes.add(new PairInt(2, -1));
        
        ones.add(new PairInt(0, -1));
        ones.add(new PairInt(1, 0));
        ones.add(new PairInt(2, 0));
        
        changeToZeroes.add(new PairInt(0, 0));
        
        int startValue = 1;
            
        replacePattern(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
        rotate90ThreeTimes(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
    }
    
    private void correctForWs(GreyscaleImage input) {
       
        /*        
        
        0  0  #         1
        0  #*<#         0
        #  #  0        -1
                       -2
       -1  0  1  2 
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
    
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-1, 0));
        zeroes.add(new PairInt(-1, -1));
        zeroes.add(new PairInt(0, -1));
        zeroes.add(new PairInt(1, 1));
        
        ones.add(new PairInt(-1, 1));
        ones.add(new PairInt(0, 1));
        ones.add(new PairInt(1, 0));
        ones.add(new PairInt(1, -1));
        
        changeToZeroes.add(new PairInt(0, 0));
        
        int startValue = 1;
            
        replacePattern(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
        rotate90ThreeTimes(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
    }
    
    private void correctForZigZag1(GreyscaleImage input) {
        
        /*
        looking for pattern
       
                 #      3
           0     #      2
        0  0  #  0      1
        0  #* #< 0      0
        #     0         -1
        
       -1  0  1  2 
        
        and removing the topmost left #'s
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
        
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-1, 0));
        zeroes.add(new PairInt(-1, -1));
        zeroes.add(new PairInt(0, -1));
        zeroes.add(new PairInt(0, -2));
        zeroes.add(new PairInt(2, 0));
        zeroes.add(new PairInt(2, -1));
        zeroes.add(new PairInt(1, 1));
        
        ones.add(new PairInt(-1, 1));
        ones.add(new PairInt(1, 0));
        ones.add(new PairInt(1, -1));
        ones.add(new PairInt(2, -2));
        ones.add(new PairInt(2, -3));
        
        changeToZeroes.add(new PairInt(1, 0));
        
        int startValue = 1;
        
        replacePattern(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
        rotate90ThreeTimes(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
    }
    
     private void correctForZigZag2(GreyscaleImage input) {
        
        /*
        looking for pattern
       
           0  0  0  #      1
           0  #* # 0       0
           0  #  0         -1
           #  0            -2
        
       -2 -1  0  1  2
        
        and removing the topmost left #'s
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-1, 1));
        zeroes.add(new PairInt(-1, 0));
        zeroes.add(new PairInt(-1, -1));
        zeroes.add(new PairInt(0, 2));
        zeroes.add(new PairInt(0, -1));
        zeroes.add(new PairInt(1, 1));
        zeroes.add(new PairInt(1, -1));
        zeroes.add(new PairInt(2, 0));
      
        ones.add(new PairInt(-1, 2));
        ones.add(new PairInt(0, 1));
        ones.add(new PairInt(1, 0));
        ones.add(new PairInt(2, -1));
        
        changeToZeroes.add(new PairInt(0, 0));
        
        int startValue = 1;
        
        replacePattern(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
        rotate90ThreeTimes(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
    }
    
    private void correctForZigZag1Alt(GreyscaleImage input) {
        
        /*
        looking for pattern
       
           0  0  #  #       1
           0  #* #< 0       0
           0  #  0  0      -1
           #  0            -2
        
       -2 -1  0  1  2        
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(1, 1));
        zeroes.add(new PairInt(0, 2));
        zeroes.add(new PairInt(-1, 1));
        zeroes.add(new PairInt(-1, 0));
        zeroes.add(new PairInt(-1, -1));
        zeroes.add(new PairInt(0, -1));
        zeroes.add(new PairInt(2, 1));
        zeroes.add(new PairInt(2, 0));
      
        ones.add(new PairInt(0, -1)); // NOT YET REV
        ones.add(new PairInt(1, 0));
        ones.add(new PairInt(1, 1));
        ones.add(new PairInt(-1, -2));
        ones.add(new PairInt(2, 1));
        
        changeToZeroes.add(new PairInt(1, 0));
        
        int startValue = 1;
        
        replacePattern(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
        rotate90ThreeTimes(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
    }
    
    private void correctForZigZag3(GreyscaleImage input) {
        
        /*
        looking for pattern
                 0  #      3
              0  #         2
           0  #< #  0      1
           0  #* 0         0
           #  0  0        -1
        
       -2 -1  0  1  2
        
        and removing the topmost left #'s
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-1, 0));
        zeroes.add(new PairInt(-1, -1));
        zeroes.add(new PairInt(0, 1));
        zeroes.add(new PairInt(0, -2));
        zeroes.add(new PairInt(1, 1));
        zeroes.add(new PairInt(1, 0));
        zeroes.add(new PairInt(1, -3));
        zeroes.add(new PairInt(2, -1));
      
        ones.add(new PairInt(-1, 1));
        ones.add(new PairInt(0, -1));
        ones.add(new PairInt(1, -1));
        ones.add(new PairInt(1, -2));
        ones.add(new PairInt(2, -3));
        
        changeToZeroes.add(new PairInt(0, -1));
        
        int startValue = 1;
        
        replacePattern(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
        rotate90ThreeTimes(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
    }
    
    /**
     * using this one unfortunately reduces the sharpness of real corners,
     * but it allows clear point ordering when extracting edges from the
     * line thinned image.
     * 
     * @param input 
     */
    private void correctForRemaining(GreyscaleImage input) {
       
        /*  
        
        #               1
        0  #* 0         0
        0  #< #        -1
        0  0  0        -2
        
       -1  0  1  2 
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
    
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-1, 2));
        zeroes.add(new PairInt(-1, 1));
        zeroes.add(new PairInt(-1, 0));
        zeroes.add(new PairInt(0, 2));
        zeroes.add(new PairInt(1, 2));
        zeroes.add(new PairInt(1, 0));
        
        ones.add(new PairInt(-1, -1));
        ones.add(new PairInt(0, 1));
        ones.add(new PairInt(1, 1));
        
        changeToZeroes.add(new PairInt(0, 1));
        
        int startValue = 1;
            
        replacePattern(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
        rotate90ThreeTimes(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
    }
    
    private void rotate90ThreeTimes(GreyscaleImage input, 
        final LinkedHashSet<PairInt> zeroes, final LinkedHashSet<PairInt> ones, 
        LinkedHashSet<PairInt> changeToZeroes, 
        final LinkedHashSet<PairInt> changeToOnes, final int startCenterValue) {
        
        // ----- change the sign of x to handle other direction -----
        for (PairInt p : zeroes) {
            p.setX(-1 * p.getX());
        }
        
        for (PairInt p : ones) {
            p.setX(-1 * p.getX());
        }
          
        for (PairInt p : changeToZeroes) {
            p.setX(-1 * p.getX());
        }
        
        for (PairInt p : changeToOnes) {
            p.setX(-1 * p.getX());
        }
                    
        replacePattern(input, zeroes, ones, changeToZeroes, changeToOnes,
            startCenterValue);
             
        // ----- change the sign of y to handle other direction -----
        for (PairInt p : zeroes) {
            p.setY(-1 * p.getY());
        }
        
        for (PairInt p : ones) {
            p.setY(-1 * p.getY());
        }
        
        for (PairInt p : changeToZeroes) {
            p.setY(-1 * p.getY());
        }
        
        for (PairInt p : changeToOnes) {
            p.setY(-1 * p.getY());
        }
        
        replacePattern(input, zeroes, ones, changeToZeroes, changeToOnes,
            startCenterValue);
        
        // ----- change the sign of x to handle another direction -----
        for (PairInt p : zeroes) {
            p.setX(-1 * p.getX());
        }
        
        for (PairInt p : ones) {
            p.setX(-1 * p.getX());
        }
        
        for (PairInt p : changeToZeroes) {
            p.setX(-1 * p.getX());
        }
        
        for (PairInt p : changeToOnes) {
            p.setX(-1 * p.getX());
        }
                    
        replacePattern(input, zeroes, ones, changeToZeroes, changeToOnes,
            startCenterValue);
    }
    
    private void replacePattern(GreyscaleImage input, 
        final LinkedHashSet<PairInt> zeroes, final LinkedHashSet<PairInt> ones, 
        final LinkedHashSet<PairInt> changeToZeroes, final LinkedHashSet<PairInt> changeToOnes, 
        final int startCenterValue) {
        
        int w = input.getWidth();
        int h = input.getHeight();
        
        for (int col = 0; col < w; col++) {
            
            for (int row = 0; row < h; row++) {
                
                int v = input.getValue(col, row);
                   
                if (v != startCenterValue) {
                    continue;
                }
                
                boolean foundPattern = true;
                
                for (PairInt p : zeroes) {
                    int x = col + p.getX();
                    int y = row + p.getY();
                    if ((x < 0) || (y < 0) || (x > (w - 1)) || (y > (h - 1))) {
                        //TODO: revisit this
                        foundPattern = false;
                        break;
                    }
                    int vz = input.getValue(x, y);
                    if (vz != 0) {
                        foundPattern = false;
                        break;
                    }
                }
                
                if (!foundPattern) {
                    continue;
                }
                
                for (PairInt p : ones) {
                    int x = col + p.getX();
                    int y = row + p.getY();
                    if ((x < 0) || (y < 0) || (x > (w - 1)) || (y > (h - 1))) {
                        foundPattern = false;
                        break;
                    }
                    int vz = input.getValue(x, y);
                    if (vz != 1) {
                        foundPattern = false;
                        break;
                    }
                }
                
                if (!foundPattern) {
                    continue;
                }
                
                for (PairInt p : changeToZeroes) {
                    int x = col + p.getX();
                    int y = row + p.getY();
                    if ((x < 0) || (y < 0) || (x > (w - 1)) || (y > (h - 1))) {
                        continue;
                    }
                    input.setValue(x, y, 0);
                }
                
                for (PairInt p : changeToOnes) {
                    int x = col + p.getX();
                    int y = row + p.getY();
                    if ((x < 0) || (y < 0) || (x > (w - 1)) || (y > (h - 1))) {
                        continue;
                    }
                    input.setValue(x, y, 1);
                }
            }
        }
    }
    
    private void debugPrint(GreyscaleImage input, int xStart, int xStop,
        int yStart, int yStop) {
        
        StringBuilder sb = new StringBuilder();
                    
        for (int row = yStart; row <= yStop; row++) {
            sb.append(String.format("%3d: ", row));
            for (int col = xStart; col <= xStop; col++) {
                sb.append(String.format(" %3d ", input.getValue(col, row)));
            }
            sb.append(String.format("\n"));
        }
        
        System.out.println(sb.toString());
    }

    /**
     * removes a hole artifact in inclined lines.  note that this should
     * probably be adjusted for gaussian convolution combined radius
     * if used outside of the gradientXY image produced by the
     * CannyEdgeFilter.
     * @param input 
     */
    private void correctForHoleArtifacts3(GreyscaleImage input) {
        
        /*
        looking for pattern
        
            0    0    0    0    0    0           3
            0    0    0    0    1    1    1      2
            0    0    0    1    0    1    0      1
            0    0    1    0*   1    1    0      0
            0    1    0    1    1    0    0     -1
            0    1    1    1    0    0    0     -2
            0    1    0    0    0    0    0     -3
        
           -3   -2   -1    0    1    2    3     
        
        and removing the topmost left #'s
        
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
        
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-1, -1));
        zeroes.add(new PairInt(-1, 1)); 
        zeroes.add(new PairInt(1, -1)); 
        zeroes.add(new PairInt(0, -2)); 
        zeroes.add(new PairInt(1, 2)); 
        zeroes.add(new PairInt(2, 2)); 
        zeroes.add(new PairInt(2, 1)); 
        zeroes.add(new PairInt(0, -3)); 
        zeroes.add(new PairInt(1, -3)); 
        zeroes.add(new PairInt(-1, -2)); 
        zeroes.add(new PairInt(-1, -3)); 
        zeroes.add(new PairInt(2, -3));
        zeroes.add(new PairInt(-2, -3));
        zeroes.add(new PairInt(-2, -2));
        zeroes.add(new PairInt(-2, -1));
        zeroes.add(new PairInt(-2, 0));
        zeroes.add(new PairInt(-1, 3));
        zeroes.add(new PairInt(0, 3));
        zeroes.add(new PairInt(1, 3));
        zeroes.add(new PairInt(2, 3));
        zeroes.add(new PairInt(3, 3));
        zeroes.add(new PairInt(3, 2));
        zeroes.add(new PairInt(3, 1));
        zeroes.add(new PairInt(3, 0));
        zeroes.add(new PairInt(3, -1));
        zeroes.add(new PairInt(-3, -3));
        zeroes.add(new PairInt(-3, 2));
        zeroes.add(new PairInt(-3, -1));
        zeroes.add(new PairInt(-3, 0));
        zeroes.add(new PairInt(-3, -1));
        zeroes.add(new PairInt(-3, -2));
        zeroes.add(new PairInt(-3, -3));
        
        ones.add(new PairInt(0, 1));
        ones.add(new PairInt(0, 2));
        ones.add(new PairInt(0, -1));
        ones.add(new PairInt(-1, 0));
        ones.add(new PairInt(-1, 2));
        ones.add(new PairInt(1, -2));
        ones.add(new PairInt(1, 0));
        ones.add(new PairInt(1, 1));
        ones.add(new PairInt(2, -2));
        ones.add(new PairInt(2, -1));
        ones.add(new PairInt(2, 0));
        ones.add(new PairInt(3, -2));
        ones.add(new PairInt(-2, 1));
        ones.add(new PairInt(-2, 2));
        ones.add(new PairInt(-2, 3));
    
        changeToZeroes.add(new PairInt(-2, 2));
        changeToZeroes.add(new PairInt(-2, 1));
        changeToZeroes.add(new PairInt(-1, 0));
        changeToZeroes.add(new PairInt(0, -1));
        changeToZeroes.add(new PairInt(1, -2));
        changeToZeroes.add(new PairInt(2, -2));
        changeToZeroes.add(new PairInt(2, 0));
        changeToZeroes.add(new PairInt(1, 1));
        changeToZeroes.add(new PairInt(0, 2));
              
        int centralValue = 0;
        
        replacePattern(input, zeroes, ones, changeToZeroes, changeToOnes, 
            centralValue);
        
        rotate90ThreeTimes(input, zeroes, ones, changeToZeroes, changeToOnes, 
            centralValue);
       
    }
    
    /**
     * removes a hole artifact in inclined lines.  note that this should
     * probably be adjusted for gaussian convolution combined radius
     * if used outside of the gradientXY image produced by the
     * CannyEdgeFilter.
     * @param input 
     */
    private void correctForHoleArtifacts2(GreyscaleImage input) {
        
        /*
        looking for pattern
        
         0    0    0    0    0    0     3
         0    0    0    1    1    1     2
         0    0    1    0    1    0     1
         0    1    0*   1    1    0     0
         0    1    1    1    0    0    -1
         0    1    0    0    0    0    -2
        
        -2   -1    0    1    2    3                
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
      
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(0, 2));
        zeroes.add(new PairInt(0, -2));
        zeroes.add(new PairInt(-1, -1));
        zeroes.add(new PairInt(-1, -2));
        zeroes.add(new PairInt(-1, -3));
        zeroes.add(new PairInt(0, -3));
        zeroes.add(new PairInt(1, -3));
        zeroes.add(new PairInt(1, -1));
        zeroes.add(new PairInt(1, 2));
        zeroes.add(new PairInt(2, -3));
        zeroes.add(new PairInt(2, 1));
        zeroes.add(new PairInt(2, 2));
        zeroes.add(new PairInt(-2, -3));
        zeroes.add(new PairInt(-2, -2));
        zeroes.add(new PairInt(-2, -1));
        zeroes.add(new PairInt(-2, 0));
        zeroes.add(new PairInt(-2, 1));
        zeroes.add(new PairInt(-2, 2));
        zeroes.add(new PairInt(3, -3));
        zeroes.add(new PairInt(3, -1));
        zeroes.add(new PairInt(3, 0));
        zeroes.add(new PairInt(3, 1));
        zeroes.add(new PairInt(3, 2));
        
        ones.add(new PairInt(0, -1));
        ones.add(new PairInt(0, 1));
        ones.add(new PairInt(-1, 0));
        ones.add(new PairInt(-1, 1));
        ones.add(new PairInt(-1, 2));
        ones.add(new PairInt(1, -2));
        ones.add(new PairInt(1, 0));
        ones.add(new PairInt(1, 1));
        ones.add(new PairInt(2, -2));
        ones.add(new PairInt(2, -1));
        ones.add(new PairInt(2, 0));
        ones.add(new PairInt(3, -2));
    
        changeToZeroes.add(new PairInt(-1, 0));
        changeToZeroes.add(new PairInt(-1, 1));
        changeToZeroes.add(new PairInt(0, -1));
        changeToZeroes.add(new PairInt(1, -2));
        changeToZeroes.add(new PairInt(1, 1));
        changeToZeroes.add(new PairInt(2, -2));
        changeToZeroes.add(new PairInt(2, 0));
        
        int centralValue = 0;
        
        replacePattern(input, zeroes, ones, changeToZeroes, changeToOnes, 
            centralValue);
        
        rotate90ThreeTimes(input, zeroes, ones, changeToZeroes, changeToOnes, 
            centralValue);
       
    }
    
    private void correctForLine0(GreyscaleImage input) {
        
        /*
        looking for pattern
       
        0  #  0         2
        0  #  0  0      1
        0  0* #  0      0
        0  #  0  0     -1
        0  #  0
        
       -1  0  1  2
        
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-1, 2));
        zeroes.add(new PairInt(-1, 1));
        zeroes.add(new PairInt(-1, 0));
        zeroes.add(new PairInt(-1, -1));
        zeroes.add(new PairInt(-1, -2));
        zeroes.add(new PairInt(1, 2));
        zeroes.add(new PairInt(1, 1));
        zeroes.add(new PairInt(1, -1));
        zeroes.add(new PairInt(1, -2));
        zeroes.add(new PairInt(2, 1));
        zeroes.add(new PairInt(2, 0));
        zeroes.add(new PairInt(2, -1));
        
        ones.add(new PairInt(0, 2));
        ones.add(new PairInt(0, 1));
        ones.add(new PairInt(0, -1));
        ones.add(new PairInt(0, -2));
        ones.add(new PairInt(1, 0));
        
        changeToZeroes.add(new PairInt(1, 0));
        changeToOnes.add(new PairInt(0, 0));
        
        int startValue = 0;
        
        replacePattern(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
        rotate90ThreeTimes(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
    }
    
    /**
     * possibly decreases sharpness of a true diagonal edge at the expense of
     * making a better line width for the edge extractor.
     * @param input 
     */
    private void correctForSpurs(GreyscaleImage input) {
        
        /*
        looking for pattern
                        
        #  #  0         1
        0  #*<0         0
        0  0  0         -1
        
       -1  0  1  2
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
      
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-1, 1));
        zeroes.add(new PairInt(-1, 0));
        zeroes.add(new PairInt(0, 1));
        zeroes.add(new PairInt(1, 1));
        zeroes.add(new PairInt(1, -1));
        
        ones.add(new PairInt(-1, -1));
        ones.add(new PairInt(0, -1));
        
        changeToZeroes.add(new PairInt(0, 0));
        
        int startValue = 1;
        
        replacePattern(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
        rotate90ThreeTimes(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
    }
    
    private void correctForHoleArtifacts1(Set<PairInt> points, int imageWidth,
        int imageHeight) {
             
        /* 
                      1               2
                 1    0    1          1     
                      1*              0
                                     -1
        
           -2   -1    0    1    2
        */  
        
        ErosionFilter erosionFilter = new ErosionFilter();
        
        Set<PairInt> ones = new HashSet<PairInt>();
        Set<PairInt> zeroes = new HashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        ones.add(new PairInt(-1, -1));
        ones.add(new PairInt(0, -2));
        ones.add(new PairInt(1, -1));
        
        zeroes.add(new PairInt(0, -1));
        
        int w = imageWidth;
        int h = imageHeight;
                
        for (PairInt p : points) {
            
            // test for the pattern of ones and zeroes in the neighbors,
            // then make a temporary set of center to zero and test if each of
            // the four sorrounding can be deleted
            
            int col = p.getX();
            int row = p.getY();
            
            boolean foundPattern = true;

            for (PairInt p2 : ones) {
                int x = col + p2.getX();
                int y = row + p2.getY();
                if ((x < 0) || (y < 0) || (x > (w - 1)) || (y > (h - 1))) {
                    foundPattern = false;
                    break;
                }
                if (!points.contains(p2)) {
                    foundPattern = false;
                    break;
                }
            }
            
            if (!foundPattern) {
                continue;
            }
            
            for (PairInt p2 : zeroes) {
                int x = col + p2.getX();
                int y = row + p2.getY();
                if ((x < 0) || (y < 0) || (x > (w - 1)) || (y > (h - 1))) {
                    foundPattern = false;
                    break;
                }
                if (points.contains(p2)) {
                    foundPattern = false;
                    break;
                }
            }
            
            if (!foundPattern) {
                continue;
            }
            
            Set<PairInt> tmpPointsRemoved = new HashSet<PairInt>();
            Set<PairInt> tmpPointsAdded = new HashSet<PairInt>();
            
            // change the central 0 to a 1
            tmpPointsAdded.add(new PairInt(0, -1));
                        
            // test if can set the surrounding 1's to 0's without disconnecting
            // lines
            for (PairInt p2 : ones) {
                
                int x = col + p2.getX();
                int y = row + p2.getY();
                
                PairInt p3 = new PairInt(x, y);
                
                // adds to tmpPointsRemoved
                boolean nullable = erosionFilter.process(p3, points, 
                    tmpPointsAdded, tmpPointsRemoved, w, h);
            }
            
            for (PairInt p2 : tmpPointsRemoved) {
                points.remove(p2);
            }
            for (PairInt p2 : tmpPointsAdded) {
                points.add(p2);
            }
        }
        
    }

    private void correctForHoleArtifacts1_2(GreyscaleImage input) {
        
        /*     
             0  0              2
          0  0  #  #  #        1
             #  0  0  #        0
             #  #  #  0  0     -1
                   0  0        -2
        
         -3 -2 -1  0  1  2
        */  

        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();

        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-3, -1));
        zeroes.add(new PairInt(-2, -1));
        zeroes.add(new PairInt(-2, -2));
        zeroes.add(new PairInt(-1, 0));
        zeroes.add(new PairInt(-1, -2));
        zeroes.add(new PairInt(0, 2));
        zeroes.add(new PairInt(1, 2));
        zeroes.add(new PairInt(1, 1));
        zeroes.add(new PairInt(2, 1));
        
        ones.add(new PairInt(-2, 1));
        ones.add(new PairInt(-2, 0));
        ones.add(new PairInt(-1, 1));
        ones.add(new PairInt(-1, -1));
        ones.add(new PairInt(0, 1));
        ones.add(new PairInt(0, -1));
        ones.add(new PairInt(1, 0));
        ones.add(new PairInt(1, -1));
        
        changeToZeroes.add(new PairInt(-2, 0));
        changeToZeroes.add(new PairInt(-1, -1));
        changeToZeroes.add(new PairInt(-1, 1));
        changeToZeroes.add(new PairInt(0, 1));
        changeToZeroes.add(new PairInt(0, -1));
        changeToZeroes.add(new PairInt(1, 0));
        
        changeToOnes.add(new PairInt(-1, 0));
        changeToOnes.add(new PairInt(0, 0));
        
        int startValue = 0;
            
        replacePattern(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
        rotate90ThreeTimes(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
       
    }
    
    private void correctForHoleArtifacts1_3(GreyscaleImage input) {
        
        /*     
             0  0              2
          0  0  #  #           1
          0  #  0  0< #  0     0
                #  #  0  0     -1
                   0  0        -2
        
         -3 -2 -1  0  1  2
        */  

        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();

        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-3, 0));
        zeroes.add(new PairInt(-3, -1));
        zeroes.add(new PairInt(-2, -1));
        zeroes.add(new PairInt(-2, -2));
        zeroes.add(new PairInt(-1, 0));
        zeroes.add(new PairInt(-1, -2));
        zeroes.add(new PairInt(0, 2));
        zeroes.add(new PairInt(1, 2));
        zeroes.add(new PairInt(1, 1));
        zeroes.add(new PairInt(2, 1));
        zeroes.add(new PairInt(2, 0));
        
        ones.add(new PairInt(-2, 0));
        ones.add(new PairInt(-1, 1));
        ones.add(new PairInt(-1, -1));
        ones.add(new PairInt(0, 1));
        ones.add(new PairInt(0, -1));
        ones.add(new PairInt(1, 0));
        
        changeToZeroes.add(new PairInt(-2, 0));
        changeToZeroes.add(new PairInt(-1, -1));
        changeToZeroes.add(new PairInt(-1, 1));
        changeToZeroes.add(new PairInt(0, -1));
        changeToZeroes.add(new PairInt(0, 1));
        changeToZeroes.add(new PairInt(1, 0));
        
        changeToOnes.add(new PairInt(-2, 1));
        changeToOnes.add(new PairInt(-1, 0));
        changeToOnes.add(new PairInt(0, 0));
        changeToOnes.add(new PairInt(1, -1));
        
        int startValue = 0;
            
        replacePattern(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
        rotate90ThreeTimes(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
       
    }

    private void correctForHoleArtifacts1_4(GreyscaleImage input) {
        
        /*     
                   0  0        3
                #  #  0  0     2
                #  0  #  0     1
             0  #  0< #        0
             0  0  #           -1
                0  0           -2
        
         -3 -2 -1  0  1  2
        */

        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();

        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-2, 1));
        zeroes.add(new PairInt(-2, 0));
        zeroes.add(new PairInt(-1, 2));
        zeroes.add(new PairInt(-1, 1));
        zeroes.add(new PairInt(0, 2));
        zeroes.add(new PairInt(0, -1));
        zeroes.add(new PairInt(0, -3));
        zeroes.add(new PairInt(1, -2));
        zeroes.add(new PairInt(1, -3));
        zeroes.add(new PairInt(2, -1));
        zeroes.add(new PairInt(2, -2));

        ones.add(new PairInt(-1, 0));
        ones.add(new PairInt(-1, -1));
        ones.add(new PairInt(-1, -2));
        ones.add(new PairInt(0, 1));
        ones.add(new PairInt(0, -2));
        ones.add(new PairInt(1, 0));
        ones.add(new PairInt(1, -1));
     
        changeToZeroes.add(new PairInt(-1, 0));
        changeToZeroes.add(new PairInt(-1,-1));
        changeToZeroes.add(new PairInt(0, 1));
        changeToZeroes.add(new PairInt(0, -2));
        changeToZeroes.add(new PairInt(1, 0));
        changeToZeroes.add(new PairInt(1, -1));
        
        changeToOnes.add(new PairInt(0, -1));
        changeToOnes.add(new PairInt(0, 0));
        changeToOnes.add(new PairInt(1, 1));
        
        int startValue = 0;
            
        replacePattern(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
        rotate90ThreeTimes(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
       
    }

    protected GreyscaleImage sumOver8Neighborhood(GreyscaleImage img) {
        
        GreyscaleImage summed = img.copyImage();
        
        int[] dxs = new int[]{-1, -1,  0,  1, 1, 1, 0, -1};
        int[] dys = new int[]{ 0, -1, -1, -1, 0, 1, 1,  1};
        
        int w = img.getWidth();
        int h = img.getHeight();
        
        // for each pixel, sum it's neighbors
        for (int col = 0; col < w; col++) {
            for (int row = 0; row < h; row++) {
                
                int sum = 0;
                
                for (int nIdx = 0; nIdx < dxs.length; nIdx++) {
                    
                    int x = dxs[nIdx] + col;
                    int y = dys[nIdx] + row;
                    
                    if ((x<0) || (y<0) || (x>(w-1)) || (y>(h-1))) {
                        continue;
                    }
                    int v = img.getValue(x, y);
                    
                    sum += v;                    
                }
                summed.setValue(col, row, sum);
            }
        }
        
        return summed;
    }

    private void correctForMinorOffsetsByIntensity(GreyscaleImage input, 
        GreyscaleImage summed) {
     
        int[] dxs = new int[]{-1, -1,  0,  1, 1, 1, 0, -1};
        int[] dys = new int[]{ 0, -1, -1, -1, 0, 1, 1,  1};
       
        ErosionFilter erosionFilter = new ErosionFilter();
        
        int w = input.getWidth();
        int h = input.getHeight();
            
        for (int col = 1; col < (w - 1); col++) {
            for (int row = 1; row < (h - 1); row++) {
                
                int v = input.getValue(col, row);
                if (v == 0) {
                    continue;
                }
                
                int vSum = summed.getValue(col, row);
                
                if (vSum == 8) {
                    continue;
                }
                                
                if (erosionFilter.doesDisconnect(input, col, row)) {
                    continue;
                }
                
                int maxSum = vSum;
                int maxIdx = -1;
                
                for (int nIdx = 0; nIdx < dxs.length; nIdx++) {
                    
                    int x = dxs[nIdx] + col;
                    int y = dys[nIdx] + row;
                    
                    if ((x<0) || (y<0) || (x>(w-1)) || (y>(h-1))) {
                        continue;
                    }
                    
                    // only compare the neighbors which are swappable
                    if (input.getValue(x, y) == 0) {
                        int sum = summed.getValue(x, y);
                        if (sum > maxSum) {
                            maxIdx = nIdx;
                            maxSum = sum;
                        }
                    }
                }
                                
                if (maxIdx > -1) {
                    int x = dxs[maxIdx] + col;
                    int y = dys[maxIdx] + row;
                    input.setValue(x, y, 1);
                    input.setValue(col, row, 0);
                }
            }    
        }
    }

    protected GreyscaleImage addOnePixelBorders(GreyscaleImage input) {
        
        int w = input.getWidth();
        int h = input.getHeight();
        
        GreyscaleImage output = new GreyscaleImage(w + 2, h + 2);
        
        //TODO: make a method internal to GreyscaleImage that uses
        //   System.arrays.copy
        for (int col = 0; col < w; col++) {
            for (int row = 0; row < h; row++) {
                output.setValue(col + 1, row + 1, input.getValue(col, row));
            }
        }
        
        return output;
    }

    protected GreyscaleImage removeOnePixelBorders(GreyscaleImage input) {
        
        int w = input.getWidth();
        int h = input.getHeight();
        
        GreyscaleImage output = new GreyscaleImage(w - 2, h - 2);
        
        //TODO: make a method internal to GreyscaleImage that uses
        //   System.arrays.copy
        for (int col = 0; col < (w - 2); col++) {
            for (int row = 0; row < (h - 2); row++) {
                output.setValue(col, row, input.getValue(col + 1, row + 1));
            }
        }
        
        return output;
    }

    /**
     * return true if at least one pixel is found on the border of the images
     * to have a non-zero value (value > 0 || value < 0).
     * @param input
     * @return 
     */
    protected boolean hasAtLeastOneBorderPixel(GreyscaleImage input) {
        
        int lastCol = input.getWidth() - 1;
        int lastRow = input.getHeight() - 1;
        
        for (int i = 0; i <= lastCol; i++) {
            if (input.getValue(i, 0) != 0) {
                return true;
            }
        }
        
        for (int i = 0; i <= lastCol; i++) {
            if (input.getValue(i, lastRow) != 0) {
                return true;
            }
        }
        
        for (int i = 0; i <= lastRow; i++) {
            if (input.getValue(0, i) != 0) {
                return true;
            }
        }
        
        for (int i = 0; i <= lastRow; i++) {
            if (input.getValue(lastCol, i) != 0) {
                return true;
            }
        }
        
        return false;
    }
    
}
