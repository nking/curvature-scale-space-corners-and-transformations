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
        
        int w = input.getWidth();
        int h = input.getHeight();
        
        Set<PairInt> points = imageProcessor.readNonZeroPixels(input);
        
        //TODO: reduce the number of patterns here if possible
        // and make sure that true corners aren't drastically reduced to less
        // usable smaller corners 
        
        //correctForHoleArtifacts3(input);
        //correctForHoleArtifacts2(input);

        correctForHoleArtifacts1(points, w, h);
        
        correctForHoleArtifacts1_2(points, w, h);

        correctForHoleArtifacts1_3(points, w, h);
        
        correctForHoleArtifacts1_4(points, w, h);
        
        correctForZigZag0(points, w, h);
        
        correctForZigZag0Alt(points, w, h);
                
        correctForZigZag1(points, w, h);
     
        correctForZigZag2(points, w, h);
        
        correctForZigZag1Alt(points, w, h);
        
        correctForZigZag3(points, w, h);
                
        correctForZigZag5(points, w, h);
        
        correctForZigZag6(points, w, h);
        
        correctForWs(points, w, h);
        
        // TODO: revisit, not sure this is always an artifact:
        correctForLine0(points, w, h);
        
        // better edge extraction at the expense of unsharpening true corners:
        correctForLs(points, w, h);
        correctForLs2(points, w, h);
        
        correctForZigZag1(points, w, h);
        
        correctForSpurs(points, w, h);
        
        correctForZigZag7(points, w, h);

        imageProcessor.writeAsBinaryToImage(input, points);
    }
    
    private void correctForZigZag0(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
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
            
        replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
    }
    
    private void correctForZigZag0Alt(Set<PairInt> points, int imageWidth,
        int imageHeight) {
       
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
            
        replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
    }
    
    private void correctForZigZag5(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
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
            
        replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
    }
    
    private void correctForZigZag6(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
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
            
        replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
    }
    
    /**
     * using this one unfortunately reduces the sharpness of real corners,
     * but it allows clear point ordering when extracting edges from the
     * line thinned image.
     * 
     * @param input 
     */
    private void correctForZigZag7(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
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
            
        replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
    }
    
    /**
     * using this one unfortunately reduces the sharpness of real corners,
     * but it allows clear point ordering when extracting edges from the
     * line thinned image.
     * 
     * @param input 
     */
    private void correctForLs(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
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
            
        replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
    }
    
    private void correctForLs2(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
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
            
        replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
    }
    
    private void correctForWs(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
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
            
        replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
    }
    
    private void correctForZigZag1(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
        
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
        
        replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
    }
    
    private void correctForZigZag2(Set<PairInt> points, int imageWidth, 
         int imageHeight) {
        
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
        
        replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
    }
    
    private void correctForZigZag1Alt(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
        
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
        
        replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
    }
    
    private void correctForZigZag3(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
        
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
        
        replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
    }
    
    /**
     * using this one unfortunately reduces the sharpness of real corners,
     * but it allows clear point ordering when extracting edges from the
     * line thinned image.
     * 
     * @param input 
     */
    private void correctForRemaining(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
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
            
        replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
    }
    
    private void rotate90ThreeTimes(
        Set<PairInt> points, int imageWidth, int imageHeight,
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
                    
        replacePattern(
            points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes,
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
        
        replacePattern(
            points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes,
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
                    
        replacePattern(
            points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes,
            startCenterValue);
    }
    
    private void replacePattern(
        Set<PairInt> points, int imageWidth, int imageHeight,
        final LinkedHashSet<PairInt> zeroes, final LinkedHashSet<PairInt> ones, 
        final LinkedHashSet<PairInt> changeToZeroes, final LinkedHashSet<PairInt> changeToOnes, 
        final int startCenterValue) {
        
        int w = imageWidth;
        int h = imageHeight;

        Set<PairInt> tmpPointsRemoved = new HashSet<PairInt>();
        Set<PairInt> tmpPointsAdded = new HashSet<PairInt>();

        for (PairInt p : points) {

            boolean isNotPresent = tmpPointsRemoved.contains(p) ||
                (!points.contains(p) && !tmpPointsAdded.contains(p));

            if (startCenterValue == 0) {
                // skip if point is in set
                if (!isNotPresent) {
                    continue;
                }
            } else if (isNotPresent) {
                // skip if point is not in set
                continue;
            }
            
            int col = p.getX();
            int row = p.getY();
            
            boolean foundPattern = true;
            
            for (PairInt p2 : zeroes) {
                int x = col + p2.getX();
                int y = row + p2.getY();
                if ((x < 0) || (y < 0) || (x > (w - 1)) || (y > (h - 1))) {
                    //TODO: revisit this
                    foundPattern = false;
                    break;
                }
                PairInt p3 = new PairInt(x, y);
                if (!tmpPointsRemoved.contains(p3)
                    && (points.contains(p3) || tmpPointsAdded.contains(p3))) {
                    foundPattern = false;
                    break;
                }
            }
                
            if (!foundPattern) {
                continue;
            }

            for (PairInt p2 : ones) {
                int x = col + p2.getX();
                int y = row + p2.getY();
                if ((x < 0) || (y < 0) || (x > (w - 1)) || (y > (h - 1))) {
                    foundPattern = false;
                    break;
                }
                PairInt p3 = new PairInt(x, y);
                if (tmpPointsRemoved.contains(p3) ||
                    (!points.contains(p3) && !tmpPointsAdded.contains(p3))) {
                    foundPattern = false;
                    break;
                }
            }
                
            if (!foundPattern) {
                continue;
            }
            
            for (PairInt p2 : changeToZeroes) {
                int x = col + p2.getX();
                int y = row + p2.getY();
                if ((x < 0) || (y < 0) || (x > (w - 1)) || (y > (h - 1))) {
                    continue;
                }
                tmpPointsRemoved.add(new PairInt(x, y));
            }

            for (PairInt p2 : changeToOnes) {
                int x = col + p2.getX();
                int y = row + p2.getY();
                if ((x < 0) || (y < 0) || (x > (w - 1)) || (y > (h - 1))) {
                    continue;
                }
                tmpPointsAdded.add(new PairInt(x, y));
            }
        }
        
        for (PairInt p2 : tmpPointsRemoved) {
            points.remove(p2);
        }
        for (PairInt p2 : tmpPointsAdded) {
            points.add(p2);
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
    private void correctForHoleArtifacts3(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
        
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
              
        int startValue = 0;
        
        replacePattern(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
       
    }
    
    /**
     * removes a hole artifact in inclined lines.  note that this should
     * probably be adjusted for gaussian convolution combined radius
     * if used outside of the gradientXY image produced by the
     * CannyEdgeFilter.
     * @param input 
     */
    private void correctForHoleArtifacts2(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
        
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
        
        int startValue = 0;
        
        replacePattern(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
       
    }
    
    private void correctForLine0(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
        
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
        
        replacePattern(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
    }
    
    /**
     * possibly decreases sharpness of a true diagonal edge at the expense of
     * making a better line width for the edge extractor.
     * @param input 
     */
    private void correctForSpurs(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
        
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
        
        replacePattern(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
    }
    
    protected void correctForHoleArtifacts1(Set<PairInt> points, int imageWidth,
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
        
        Set<PairInt> tmpPointsRemoved = new HashSet<PairInt>();
        Set<PairInt> tmpPointsAdded = new HashSet<PairInt>();
                    
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
                PairInt p3 = new PairInt(x, y);
                if (tmpPointsRemoved.contains(p3) ||
                    (!points.contains(p3) && !tmpPointsAdded.contains(p3))) {
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
                PairInt p3 = new PairInt(x, y);
                if (!tmpPointsRemoved.contains(p3) 
                    && (points.contains(p3) || tmpPointsAdded.contains(p3))) {
                    foundPattern = false;
                    break;
                }
            }
            
            if (!foundPattern) {
                continue;
            }
            
            // change the central 0 to a 1
            tmpPointsAdded.add(new PairInt(col, row - 1));
                        
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
        }
        
        for (PairInt p2 : tmpPointsRemoved) {
            points.remove(p2);
        }
        for (PairInt p2 : tmpPointsAdded) {
            points.add(p2);
        }
    }

    private void correctForHoleArtifacts1_2(Set<PairInt> points, int imageWidth,
        int imageHeight) {
        
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
        
        replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
        rotate90ThreeTimes(
            points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
       
    }
    
    private void correctForHoleArtifacts1_3(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
        
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
            
        replacePattern(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
       
    }

    private void correctForHoleArtifacts1_4(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
        
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
            
        replacePattern(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
       
    }

}
