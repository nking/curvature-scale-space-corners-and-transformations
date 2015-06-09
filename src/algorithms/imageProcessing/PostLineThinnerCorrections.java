package algorithms.imageProcessing;

import algorithms.CountingSort;
import algorithms.util.PairInt;
import algorithms.util.ResourceFinder;
import java.io.IOException;
import java.util.Arrays;
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
    
    protected static final int[] fourNeighborsX = new int[]{0,  1,  0,  -1};
    protected static final int[] fourNeighborsY = new int[]{-1, 0,  1,   0};
    
    protected static final int[] eightNeighborsX = 
        new int[]{-1, -1,  0,  1, 1, 1, 0, -1};
    protected static final int[] eightNeighborsY = 
        new int[]{ 0, -1, -1, -1, 0, 1, 1,  1};
    
    public void correctForArtifacts(GreyscaleImage input) {
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        int w = input.getWidth();
        int h = input.getHeight();
        
        Set<PairInt> points = imageProcessor.readNonZeroPixels(input);
        
        //TODO: reduce the number of patterns here if possible
        // and make sure that true corners aren't drastically reduced to less
        // usable smaller corners 
        
try {
String dirPath = ResourceFinder.findDirectory("bin");
ImageIOHelper.writeOutputImage(dirPath + "/nonZero.png", input);
} catch (IOException e){}

        correctForHoleArtifacts00_0(points, w, h);
        correctForHoleArtifacts00_1(points, w, h);
        correctForHoleArtifacts00_2(points, w, h);
        correctForHoleArtifacts00_3(points, w, h);
        correctForHoleArtifacts00_4(points, w, h);
       
        /*
        correctForHoleArtifacts0(points, w, h);
        correctForHoleArtifacts0_1(points, w, h);
        correctForHoleArtifacts0_2(points, w, h);
        correctForHoleArtifacts0_3(points, w, h);
        correctForHoleArtifacts0_4(points, w, h);
        correctForHoleArtifacts0_5(points, w, h);

        correctForHoleArtifacts1(points, w, h);
        correctForHoleArtifacts1_1(points, w, h);
       
        correctForHoleArtifacts1_2(points, w, h);
        correctForHoleArtifacts1_2_1(points, w, h);

        correctForHoleArtifacts1_3(points, w, h);
        
        correctForHoleArtifacts1_4(points, w, h);
        
        correctForZigZag00(points, w, h);
        correctForZigZag00_1(points, w, h);
        correctForZigZag00_2(points, w, h);
        correctForZigZag00_3(points, w, h);
        correctForZigZag00_5(points, w, h);
        correctForZigZag00_6(points, w, h);
        correctForZigZag00_7(points, w, h);
        correctForZigZag00_8(points, w, h);
        correctForZigZag00_9(points, w, h);
        correctForZigZag00_10(points, w, h);
        correctForZigZag00_11(points, w, h);
        correctForZigZag00_12(points, w, h);
        correctForZigZag00_13(points, w, h);
        
        correctForZigZag00_4(points, w, h);
        
        correctForTs(points, w, h);
        
        correctForWs(points, w, h);
        
        // TODO: revisit, not sure this is always an artifact:
        correctForLine0(points, w, h);
        
        // better edge extraction at the expense of unsharpening true corners:
        correctForLs(points, w, h);
        correctForLs2(points, w, h);
                
        correctForZigZag00_14(points, w, h);
        correctForZigZag00_15(points, w, h);
        correctForZigZag00_16(points, w, h);
        
        //correctForSpurs(points, w, h);
        */
        imageProcessor.writeAsBinaryToImage(input, points);

try {
String dirPath = ResourceFinder.findDirectory("bin");
ImageIOHelper.writeOutputImage(dirPath + "/nonZero2.png", input);
} catch (IOException e){}

    }
    
    private void correctForZigZag00(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*       
           0  0     #   2
        0  0  #  #      1
        0  #* #  0      0
           #  0  0     -1
        #              -2
       -1  0  1  2  3        
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-1, 0)); zeroes.add(new PairInt(-1, -1));
        zeroes.add(new PairInt(0, -1)); zeroes.add(new PairInt(0, -2));
        zeroes.add(new PairInt(1, 1)); zeroes.add(new PairInt(1, -2));
        zeroes.add(new PairInt(2, 1)); zeroes.add(new PairInt(2, 0));
        
        ones.add(new PairInt(-1, 2));
        ones.add(new PairInt(0, 1));
        ones.add(new PairInt(1, 0)); ones.add(new PairInt(1, -1));
        ones.add(new PairInt(2, -1));
        ones.add(new PairInt(3, -2));
        
        changeToZeroes.add(new PairInt(0, 0));
        changeToZeroes.add(new PairInt(1, -1));
                    
        replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
    }
    
    private void correctForZigZag00_1(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*       
                   0  0  #      2
                0  0  #         1
                0  #* #  0      0
                   #  0  0     -1
                #  0  0        -2
        
           -2  -1  0  1  2  3        
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-1, 0)); zeroes.add(new PairInt(-1, -1));
        zeroes.add(new PairInt(0, 2)); zeroes.add(new PairInt(0, -1)); zeroes.add(new PairInt(0, -2));
        zeroes.add(new PairInt(1, 2)); zeroes.add(new PairInt(1, 1)); zeroes.add(new PairInt(1, -2));
        zeroes.add(new PairInt(2, 1)); zeroes.add(new PairInt(2, 0));
        
        ones.add(new PairInt(-1, 2));
        ones.add(new PairInt(0, 1));
        ones.add(new PairInt(1, 0)); ones.add(new PairInt(1, -1));
        ones.add(new PairInt(2, -2));
        
        changeToZeroes.add(new PairInt(1, 0));
                    
        replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
    }
    
    private void correctForZigZag00_2(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*       
                   0  0     #   2
                0  0  #  #      1
                0  #* #< 0      0
                #     0  0     -1
             #                 -2
        
            -2 -1  0  1  2  3        
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-1, 0)); zeroes.add(new PairInt(-1, -1));
        zeroes.add(new PairInt(0, -1)); zeroes.add(new PairInt(0, -2));
        zeroes.add(new PairInt(1, 1)); zeroes.add(new PairInt(1, -2));
        zeroes.add(new PairInt(2, 1)); zeroes.add(new PairInt(2, 0));
        
        ones.add(new PairInt(-2, 2));
        ones.add(new PairInt(-1, 1));
        ones.add(new PairInt(1, 0)); ones.add(new PairInt(1, -1));
        ones.add(new PairInt(2, -1));
        ones.add(new PairInt(3, -2));
        
        changeToZeroes.add(new PairInt(1, 0));
                    
        replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
    }
    
    private void correctForZigZag00_3(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*       
             0  0  0  0         2
             #  #  #< 0         1
             0  0  #* #         0
                0  0  0  #     -1
        
            -2 -1  0  1  2     
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-2, 0)); zeroes.add(new PairInt(-2, -2));
        zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, 0)); zeroes.add(new PairInt(-1, -2));
        zeroes.add(new PairInt(0, 1)); zeroes.add(new PairInt(0, -2));
        zeroes.add(new PairInt(1, 1)); zeroes.add(new PairInt(1, -1)); zeroes.add(new PairInt(1, -2));
        
        ones.add(new PairInt(-2, -1));
        ones.add(new PairInt(-1, -1));
        ones.add(new PairInt(0, -1));
        ones.add(new PairInt(1, 0));
        ones.add(new PairInt(2, 1));
        
        changeToZeroes.add(new PairInt(0, -1));
                    
        replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
    }
    
    private void correctForZigZag00_5(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*       
          #  0  0  0  0         2
             #  #  #< 0         1
             0  0  #* #         0
                0  0  0  #     -1
        
         -3 -2 -1  0  1  2     
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        ones.add(new PairInt(-3, -2));
        ones.add(new PairInt(-2, -1));
        ones.add(new PairInt(-1, -1));
        ones.add(new PairInt(0, -1));
        ones.add(new PairInt(1, 0));
        
        zeroes.add(new PairInt(-2, 0)); zeroes.add(new PairInt(-2, -2));
        zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, 0)); zeroes.add(new PairInt(-1, -2));
        zeroes.add(new PairInt(0, 1)); zeroes.add(new PairInt(0, -2));
        zeroes.add(new PairInt(1, 1)); zeroes.add(new PairInt(1, -1)); zeroes.add(new PairInt(1, -2));
        
        changeToZeroes.add(new PairInt(0, -1));
                    
        replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
    }
    
    private void correctForZigZag00_6(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*       
             #  0  0  0         2
             0  #  #  0         1
             0  0  #*<#  #      0
                0  0  0  #     -1
        
         -3 -2 -1  0  1  2     
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        ones.add(new PairInt(-2, -2));
        ones.add(new PairInt(-1, -1));
        ones.add(new PairInt(0, -1));
        ones.add(new PairInt(1, 0));
        ones.add(new PairInt(2, 1)); ones.add(new PairInt(2, 0));
        
        zeroes.add(new PairInt(-2, 0)); zeroes.add(new PairInt(-2, -1));
        zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, 0)); zeroes.add(new PairInt(-1, -2));
        zeroes.add(new PairInt(0, 1)); zeroes.add(new PairInt(0, -2));
        zeroes.add(new PairInt(1, 1)); zeroes.add(new PairInt(1, -1)); zeroes.add(new PairInt(1, -2));
        
        changeToZeroes.add(new PairInt(0, 0));
                    
        replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
    }
    
    private void correctForZigZag00_7(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*       
             #  0  0  0         2
             0  #  #  0         1
             0  0  #*<#  0      0
                0  0  0  #     -1
        
         -3 -2 -1  0  1  2     
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        ones.add(new PairInt(-2, -2));
        ones.add(new PairInt(-1, -1));
        ones.add(new PairInt(0, -1));
        ones.add(new PairInt(1, 0));
        ones.add(new PairInt(2, 1));
        
        zeroes.add(new PairInt(-2, 0)); zeroes.add(new PairInt(-2, -1));
        zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, 0)); zeroes.add(new PairInt(-1, -2));
        zeroes.add(new PairInt(0, 1)); zeroes.add(new PairInt(0, -2));
        zeroes.add(new PairInt(1, 1)); zeroes.add(new PairInt(1, -1)); zeroes.add(new PairInt(1, -2));
        zeroes.add(new PairInt(2, 0));
        
        changeToZeroes.add(new PairInt(0, 0));
                    
        replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
    }
    
    private void correctForZigZag00_8(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*       
          #  0  0  0  0         2
             #  #  #< 0  0      1
             0  0  #* #  0      0
                0  0  #  #     -1
        
         -3 -2 -1  0  1  2     
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        ones.add(new PairInt(-3, -2));
        ones.add(new PairInt(-2, -1));
        ones.add(new PairInt(-1, -1));
        ones.add(new PairInt(0, -1));
        ones.add(new PairInt(1, 1)); ones.add(new PairInt(1, 0));
        ones.add(new PairInt(2, 1));
        
        zeroes.add(new PairInt(-2, 0)); zeroes.add(new PairInt(-2, -2));
        zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, 0)); zeroes.add(new PairInt(-1, -2));
        zeroes.add(new PairInt(0, 1)); zeroes.add(new PairInt(0, -2));
        zeroes.add(new PairInt(1, -1)); zeroes.add(new PairInt(1, -2));
        zeroes.add(new PairInt(2, 0)); zeroes.add(new PairInt(2, -1));
        
        changeToZeroes.add(new PairInt(0, -1));
                    
        replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
    }
    
    private void correctForZigZag00_9(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*       
             #  0  0  0         2
             0  #  #  0  0      1
             0  0  #*<#  #      0
                0  0  0  0     -1
        
         -3 -2 -1  0  1  2     
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        ones.add(new PairInt(-2, -2));
        ones.add(new PairInt(-1, -1));
        ones.add(new PairInt(0, -1));
        ones.add(new PairInt(1, 0));
        ones.add(new PairInt(2, 0));
        
        zeroes.add(new PairInt(-2, 0)); zeroes.add(new PairInt(-2, -1));
        zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, 0)); zeroes.add(new PairInt(-1, -2));
        zeroes.add(new PairInt(0, 1)); zeroes.add(new PairInt(0, -2));
        zeroes.add(new PairInt(1, 1)); zeroes.add(new PairInt(1, -1)); zeroes.add(new PairInt(1, -2));
        zeroes.add(new PairInt(2, 1)); zeroes.add(new PairInt(2, -1));
        
        changeToZeroes.add(new PairInt(0, 0));
                    
        replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
    }
    
    private void correctForZigZag00_10(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*       
                0  0            2
             0  0  #  #         1
             0  #< #* 0         0
             #  #  0  0        -1
        
         -3 -2 -1  0  1  2     
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        ones.add(new PairInt(-2, 1));
        ones.add(new PairInt(-1, 1)); ones.add(new PairInt(-1, 0));
        ones.add(new PairInt(0, -1));
        ones.add(new PairInt(1, -1));
        
        zeroes.add(new PairInt(-2, 0)); zeroes.add(new PairInt(-2, -1));
        zeroes.add(new PairInt(-1, -1)); zeroes.add(new PairInt(-1, -2));
        zeroes.add(new PairInt(0, 1)); zeroes.add(new PairInt(0, -2));
        zeroes.add(new PairInt(1, 1)); zeroes.add(new PairInt(1, 0));
        
        changeToZeroes.add(new PairInt(-1, 0));
                    
        replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
    }
    
    private void correctForZigZag00_11(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*       
                                2
                0  0  0         1
                #  #* 0         0
                0  #  #        -1
                0  0  0        -2
         -3 -2 -1  0  1  2     
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        ones.add(new PairInt(-1, 0));
        ones.add(new PairInt(0, 1));
        ones.add(new PairInt(1, 1));
        
        zeroes.add(new PairInt(-1, 2)); zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, -1));
        zeroes.add(new PairInt(0, 2)); zeroes.add(new PairInt(0, -1));
        zeroes.add(new PairInt(1, 2)); zeroes.add(new PairInt(1, 0)); zeroes.add(new PairInt(1, -1));
        
        changeToZeroes.add(new PairInt(0, 0));
                    
        replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
    }
    
    private void correctForZigZag00_12(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*       
                0  0  0         2
             0  0  #  #         1
             0  #  #* 0         0
             0  #  0  0        -1
                               -2
         -3 -2 -1  0  1  2     
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        ones.add(new PairInt(-1, 1)); ones.add(new PairInt(-1, 0));
        ones.add(new PairInt(0, -1));
        ones.add(new PairInt(1, -1));
        
        zeroes.add(new PairInt(-2, 1)); zeroes.add(new PairInt(-2, 0)); zeroes.add(new PairInt(-2, -1));
        zeroes.add(new PairInt(-1, -1)); zeroes.add(new PairInt(-1, -2));
        zeroes.add(new PairInt(0, 1)); zeroes.add(new PairInt(0, -2));
        zeroes.add(new PairInt(1, 1)); zeroes.add(new PairInt(1, 0)); zeroes.add(new PairInt(1, -2));
        
        changeToZeroes.add(new PairInt(0, 0));
                    
        replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
    }
    
    private void correctForZigZag00_13(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*       
                0  0  #         1
                0  #* #         0
                0  #  0        -1
        
         -3 -2 -1  0  1  2     
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        ones.add(new PairInt(0, 1));
        ones.add(new PairInt(1, 0)); ones.add(new PairInt(1, -1));
        
        zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, 0)); zeroes.add(new PairInt(-1, -1));
        zeroes.add(new PairInt(0, -1));
        zeroes.add(new PairInt(1, 1));
        
        changeToZeroes.add(new PairInt(0, 0));
                    
        replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
    }
    
    private void correctForZigZag00_14(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*       
                0  #  0         1
                0  #*<#  0      0
                0  #  #< 0     -1
                #  0  0        -2
         -3 -2 -1  0  1  2     
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        ones.add(new PairInt(-1, 2));
        ones.add(new PairInt(0, 1)); ones.add(new PairInt(0, -1));
        ones.add(new PairInt(1, 1)); ones.add(new PairInt(1, 0));
        
        zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, 0)); zeroes.add(new PairInt(-1, -1));
        zeroes.add(new PairInt(0, 2));
        zeroes.add(new PairInt(1, 2)); zeroes.add(new PairInt(1, -1));
        zeroes.add(new PairInt(2, 1)); zeroes.add(new PairInt(2, 0));
        
        changeToZeroes.add(new PairInt(0, 0));
        changeToZeroes.add(new PairInt(1, 1));
                    
        replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
    }
    
    private void correctForZigZag00_15(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*       
                #  0            2
                0  #  0         1
                0  #*<#  0      0
                0  #< #  0     -1
                   0  0  #     -2
         -3 -2 -1  0  1  2     
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        ones.add(new PairInt(-1, -2));
        ones.add(new PairInt(0, 1)); ones.add(new PairInt(0, -1));
        ones.add(new PairInt(1, 1)); ones.add(new PairInt(1, 0));
        ones.add(new PairInt(2, 2));
        
        zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, 0)); zeroes.add(new PairInt(-1, -1));
        zeroes.add(new PairInt(0, 2)); zeroes.add(new PairInt(0, -2));
        zeroes.add(new PairInt(1, 2)); zeroes.add(new PairInt(1, -1));
        zeroes.add(new PairInt(2, 1)); zeroes.add(new PairInt(2, 0));
        
        changeToZeroes.add(new PairInt(0, 1));
        changeToZeroes.add(new PairInt(0, 0));
                    
        replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
    }
    
    private void correctForZigZag00_16(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*       
                                2
                0  0  0  #      1
                0  #* #  0      0
             0  #  #< #< 0     -1
             #  0  0  0  0     -2
         -3 -2 -1  0  1  2     
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        ones.add(new PairInt(-2, 2));
        ones.add(new PairInt(-1, 1));
        ones.add(new PairInt(0, 1));
        ones.add(new PairInt(1, 1)); ones.add(new PairInt(1, 0));
        ones.add(new PairInt(2, -1));
        
        zeroes.add(new PairInt(-2, 1));
        zeroes.add(new PairInt(-1, 2)); zeroes.add(new PairInt(-1, 0)); zeroes.add(new PairInt(-1, -1));
        zeroes.add(new PairInt(0, 2)); zeroes.add(new PairInt(0, -1));
        zeroes.add(new PairInt(1, 2)); zeroes.add(new PairInt(1, -1));
        zeroes.add(new PairInt(2, 2)); zeroes.add(new PairInt(2, 1)); zeroes.add(new PairInt(2, 0));
        
        changeToZeroes.add(new PairInt(0, 1));
        changeToZeroes.add(new PairInt(1, 1));
                    
        replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
    }
    
    private void correctForTs(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*       
          0  0  0  0            2
          #  #  #  0  0  0      1
          0  #< 0  #* #  #      0
          0  0  0  0  #< 0     -1
                   0  0  0     -2
         -3 -2 -1  0  1  2     
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        ones.add(new PairInt(-3, -1));
        ones.add(new PairInt(-2, 0));ones.add(new PairInt(-2, -1));
        ones.add(new PairInt(-1, -1));
        ones.add(new PairInt(1, 1)); ones.add(new PairInt(1, 0));
        ones.add(new PairInt(2, 0));
        
        zeroes.add(new PairInt(-3, 1)); zeroes.add(new PairInt(-3, 0)); zeroes.add(new PairInt(-3, -2));
        zeroes.add(new PairInt(-2, 1)); zeroes.add(new PairInt(-2, -2));
        zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, 0)); zeroes.add(new PairInt(-1, -2));
        zeroes.add(new PairInt(0, 2)); zeroes.add(new PairInt(0, 1)); zeroes.add(new PairInt(0, -1)); zeroes.add(new PairInt(0, -2));
        zeroes.add(new PairInt(1, 2)); zeroes.add(new PairInt(1, -1));
        zeroes.add(new PairInt(2, 2)); zeroes.add(new PairInt(2, 1)); zeroes.add(new PairInt(2, -1));
        
        changeToZeroes.add(new PairInt(-2, 0));
        changeToZeroes.add(new PairInt(1, 1));
                    
        replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
    }
    
    private void correctForZigZag00_4(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*       
                0  0  0  #      1
                0  #* #  0      0
                0  #  0  0     -1
                #  0  0        -2
        
           -2  -1  0  1  2  3        
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, 0)); zeroes.add(new PairInt(-1, -1));
        zeroes.add(new PairInt(0, 2)); zeroes.add(new PairInt(0, -1));
        zeroes.add(new PairInt(1, 2)); zeroes.add(new PairInt(1, 1)); zeroes.add(new PairInt(1, -1));
        zeroes.add(new PairInt(2, 1)); zeroes.add(new PairInt(2, 0));
        
        ones.add(new PairInt(-1, 2));
        ones.add(new PairInt(0, 1));
        ones.add(new PairInt(1, 0));
        ones.add(new PairInt(2, -1));
        
        changeToZeroes.add(new PairInt(0, 0));
                    
        replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
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
                    
        replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
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
                    
        replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
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
                    
        replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
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
                    
        replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
    }
    
    private void reverseXs(
        final Set<PairInt> zeroes, final Set<PairInt> ones, 
        Set<PairInt> changeToZeroes, final Set<PairInt> changeToOnes) {
        
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
    }
    
    private void reverseXs(final Set<PairInt> zeroes, final Set<PairInt> ones) {
        
        // ----- change the sign of x to handle other direction -----
        for (PairInt p : zeroes) {
            p.setX(-1 * p.getX());
        }
        
        for (PairInt p : ones) {
            p.setX(-1 * p.getX());
        }
        
    }
    
    private void reverseYs(final Set<PairInt> zeroes, final Set<PairInt> ones, 
        Set<PairInt> changeToZeroes, final Set<PairInt> changeToOnes) {
        
        // ----- change the sign of y  -----
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
        
    }
    
    private void reverseYs(final Set<PairInt> zeroes, final Set<PairInt> ones) {
        
        // ----- change the sign of y  -----
        for (PairInt p : zeroes) {
            p.setY(-1 * p.getY());
        }
        
        for (PairInt p : ones) {
            p.setY(-1 * p.getY());
        }
        
    }
    
    private void rotate90ThreeTimes(
        Set<PairInt> points, int imageWidth, int imageHeight,
        final LinkedHashSet<PairInt> zeroes, final LinkedHashSet<PairInt> ones, 
        LinkedHashSet<PairInt> changeToZeroes, 
        final LinkedHashSet<PairInt> changeToOnes) {
        
        // ----- change the sign of x to handle other direction -----
        reverseXs(zeroes, ones, changeToZeroes, changeToOnes);
        
        replacePattern(
            points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
             
        // ----- change the sign of y to handle other direction -----
        reverseYs(zeroes, ones, changeToZeroes, changeToOnes);
        
        replacePattern(
            points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        // ----- change the sign of x to handle another direction -----
        reverseXs(zeroes, ones, changeToZeroes, changeToOnes);
                    
        replacePattern(
            points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
    }
    
    private void replacePattern(
        Set<PairInt> points, int imageWidth, int imageHeight,
        final LinkedHashSet<PairInt> zeroes, final LinkedHashSet<PairInt> ones, 
        final LinkedHashSet<PairInt> changeToZeroes, 
        final LinkedHashSet<PairInt> changeToOnes) {
        
        int w = imageWidth;
        int h = imageHeight;

        Set<PairInt> tmpPointsRemoved = new HashSet<PairInt>();
        Set<PairInt> tmpPointsAdded = new HashSet<PairInt>();

        for (PairInt p : points) {
            
            boolean isNotPresent = tmpPointsRemoved.contains(p);
            if (isNotPresent) {
                continue;
            }
            
            int col = p.getX();
            int row = p.getY();           
            
            boolean foundPattern = true;
//43,429    
if (col==43 && row==429) {
    int z = 1;
}            
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
        
            0    0    0    0    0    0           3 4
            0    0    0    0    1    1    1      2 3
            0    0    0    1    0    1    0      1 2
            0    0    1    0*   1    1    0      0 1
            0    1    0    1**  1    0    0     -1 0
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
        zeroes.add(new PairInt(-1, -1 - 1));
        zeroes.add(new PairInt(-1, 1 - 1));
        zeroes.add(new PairInt(1, -1 - 1));
        zeroes.add(new PairInt(0, -2 - 1));
        zeroes.add(new PairInt(1, 2 - 1));
        zeroes.add(new PairInt(2, 2 - 1)); 
        zeroes.add(new PairInt(2, 1 - 1));
        zeroes.add(new PairInt(0, -3 - 1));
        zeroes.add(new PairInt(1, -3 - 1));
        zeroes.add(new PairInt(-1, -2 - 1));
        zeroes.add(new PairInt(-1, -3 - 1));
        zeroes.add(new PairInt(2, -3 - 1));
        zeroes.add(new PairInt(-2, -3 - 1));
        zeroes.add(new PairInt(-2, -2 - 1));
        zeroes.add(new PairInt(-2, -1 - 1));
        zeroes.add(new PairInt(-2, 0 - 1));
        zeroes.add(new PairInt(-1, 3 - 1));
        zeroes.add(new PairInt(0, 3 - 1));
        zeroes.add(new PairInt(1, 3 - 1));
        zeroes.add(new PairInt(2, 3 - 1));
        zeroes.add(new PairInt(3, 3 - 1));
        zeroes.add(new PairInt(3, 2 - 1));
        zeroes.add(new PairInt(3, 1 - 1));
        zeroes.add(new PairInt(3, 0 - 1));
        zeroes.add(new PairInt(3, -1 - 1));
        zeroes.add(new PairInt(-3, -3 - 1));
        zeroes.add(new PairInt(-3, 2 - 1));
        zeroes.add(new PairInt(-3, -1 - 1));
        zeroes.add(new PairInt(-3, 0 - 1));
        zeroes.add(new PairInt(-3, -1 - 1));
        zeroes.add(new PairInt(-3, -2 - 1));
        zeroes.add(new PairInt(-3, -3 - 1));
        
        ones.add(new PairInt(0, 1 - 1));
        ones.add(new PairInt(0, 2 - 1));
        ones.add(new PairInt(0, -1 - 1));
        ones.add(new PairInt(-1, 0 - 1));
        ones.add(new PairInt(-1, 2 - 1));
        ones.add(new PairInt(1, -2 - 1));
        ones.add(new PairInt(1, 0 - 1));
        ones.add(new PairInt(1, 1 - 1));
        ones.add(new PairInt(2, -2 - 1));
        ones.add(new PairInt(2, -1 - 1));
        ones.add(new PairInt(2, 0 - 1));
        ones.add(new PairInt(3, -2 - 1));
        ones.add(new PairInt(-2, 1 - 1));
        ones.add(new PairInt(-2, 2 - 1));
        ones.add(new PairInt(-2, 3 - 1));
    
        changeToZeroes.add(new PairInt(-2, 2 - 1));
        changeToZeroes.add(new PairInt(-2, 1 - 1));
        changeToZeroes.add(new PairInt(-1, 0 - 1));
        changeToZeroes.add(new PairInt(0, -1 - 1));
        changeToZeroes.add(new PairInt(1, -2 - 1));
        changeToZeroes.add(new PairInt(2, -2 - 1));
        changeToZeroes.add(new PairInt(2, 0 - 1));
        changeToZeroes.add(new PairInt(1, 1 - 1));
        changeToZeroes.add(new PairInt(0, 2 - 1));
                      
        replacePattern(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
       
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
        
         0    0    0    0    0    0     3   4
         0    0    0    1    1    1     2   3
         0    0    1    0    1    0     1   2
         0    1    0    1    1    0     0   1
         0    1    1*    1    0    0    -1  0
         0    1    0    0    0    0    -2   -1
        
        -2   -1    0    1    2    3                
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
      
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(0, 2 - 1));
        zeroes.add(new PairInt(0, -2 - 1));
        zeroes.add(new PairInt(-1, -1 - 1));
        zeroes.add(new PairInt(-1, -2 - 1));
        zeroes.add(new PairInt(-1, -3 - 1));
        zeroes.add(new PairInt(0, -3 - 1));
        zeroes.add(new PairInt(1, -3 - 1));
        zeroes.add(new PairInt(1, -1 - 1));
        zeroes.add(new PairInt(1, 2 - 1));
        zeroes.add(new PairInt(2, -3 - 1));
        zeroes.add(new PairInt(2, 1 - 1));
        zeroes.add(new PairInt(2, 2 - 1));
        zeroes.add(new PairInt(-2, -3 - 1));
        zeroes.add(new PairInt(-2, -2 - 1));
        zeroes.add(new PairInt(-2, -1 - 1));
        zeroes.add(new PairInt(-2, 0 - 1));
        zeroes.add(new PairInt(-2, 1 - 1));
        zeroes.add(new PairInt(-2, 2 - 1));
        zeroes.add(new PairInt(3, -3 - 1));
        zeroes.add(new PairInt(3, -1 - 1));
        zeroes.add(new PairInt(3, 0 - 1));
        zeroes.add(new PairInt(3, 1 - 1));
        zeroes.add(new PairInt(3, 2 - 1));
        
        ones.add(new PairInt(0, -1 - 1));
        ones.add(new PairInt(0, 1 - 1));
        ones.add(new PairInt(-1, 0 - 1));
        ones.add(new PairInt(-1, 1 - 1));
        ones.add(new PairInt(-1, 2 - 1));
        ones.add(new PairInt(1, -2 - 1));
        ones.add(new PairInt(1, 0 - 1));
        ones.add(new PairInt(1, 1 - 1));
        ones.add(new PairInt(2, -2 - 1));
        ones.add(new PairInt(2, -1 - 1));
        ones.add(new PairInt(2, 0 - 1));
        ones.add(new PairInt(3, -2 - 1));
    
        changeToZeroes.add(new PairInt(-1, 0 - 1));
        changeToZeroes.add(new PairInt(-1, 1 - 1));
        changeToZeroes.add(new PairInt(0, -1 - 1));
        changeToZeroes.add(new PairInt(1, -2 - 1));
        changeToZeroes.add(new PairInt(1, 1 - 1));
        changeToZeroes.add(new PairInt(2, -2 - 1));
        changeToZeroes.add(new PairInt(2, 0 - 1));
                
        replacePattern(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
       
    }
    
    private void correctForLine0(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
        
        /*
        looking for pattern
       
        0  #  0         2
        0  #  0  0      1
        0  0  #  0      0
        0  #* 0  0     -1
        0  #  0
        
       -1  0  1  2
        
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-1, 2 - 1));
        zeroes.add(new PairInt(-1, 1 - 1));
        zeroes.add(new PairInt(-1, 0 - 1));
        zeroes.add(new PairInt(-1, -1 - 1));
        zeroes.add(new PairInt(-1, -2 - 1));
        zeroes.add(new PairInt(1, 2 - 1));
        zeroes.add(new PairInt(1, 1 - 1));
        zeroes.add(new PairInt(1, -1 - 1));
        zeroes.add(new PairInt(1, -2 - 1));
        zeroes.add(new PairInt(2, 1 - 1));
        zeroes.add(new PairInt(2, 0 - 1));
        zeroes.add(new PairInt(2, -1 - 1));
        
        ones.add(new PairInt(0, 2 - 1));
        ones.add(new PairInt(0, 1 - 1));
        ones.add(new PairInt(0, -1 - 1));
        ones.add(new PairInt(0, -2 - 1));
        ones.add(new PairInt(1, 0 - 1));
        
        changeToZeroes.add(new PairInt(1, 0 - 1));
        changeToOnes.add(new PairInt(0, 0 - 1));
                
        replacePattern(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
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
                
        replacePattern(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
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
             #  #  #* 0  0     -1
                   0  0        -2
        
         -3 -2 -1  0  1  2
        */  

        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();

        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-3, -1 - 1));
        zeroes.add(new PairInt(-2, -1 - 1));
        zeroes.add(new PairInt(-2, -2 - 1));
        zeroes.add(new PairInt(-1, 0 - 1));
        zeroes.add(new PairInt(-1, -2 - 1));
        zeroes.add(new PairInt(0, 2 - 1));
        zeroes.add(new PairInt(1, 2 - 1));
        zeroes.add(new PairInt(1, 1 - 1));
        zeroes.add(new PairInt(2, 1 - 1));
        
        ones.add(new PairInt(-2, 1 - 1));
        ones.add(new PairInt(-2, 0 - 1));
        ones.add(new PairInt(-1, 1 - 1));
        ones.add(new PairInt(-1, -1 - 1));
        ones.add(new PairInt(0, 1 - 1));
        ones.add(new PairInt(0, -1 - 1));
        ones.add(new PairInt(1, 0 - 1));
        ones.add(new PairInt(1, -1 - 1));
        
        changeToZeroes.add(new PairInt(-2, 0 - 1));
        changeToZeroes.add(new PairInt(-1, -1 - 1));
        changeToZeroes.add(new PairInt(-1, 1 - 1));
        changeToZeroes.add(new PairInt(0, 1 - 1));
        changeToZeroes.add(new PairInt(0, -1 - 1));
        changeToZeroes.add(new PairInt(1, 0 - 1));
        
        changeToOnes.add(new PairInt(-1, 0 - 1));
        changeToOnes.add(new PairInt(0, 0 - 1));
                
        replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        rotate90ThreeTimes(
            points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
       
    }
    
    private void correctForHoleArtifacts1_2_1(Set<PairInt> points, int imageWidth,
        int imageHeight) {
        
        /*     
                               
                #  #  #        2
             #  0  0  #        1
             #  #  #* #        0
        
         -3 -2 -1  0  1  2
        */  

        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();

        // y's are inverted here because sketch above is top left is (0,0)
        ones.add(new PairInt(-2, 0));
        ones.add(new PairInt(-2, -1));
        ones.add(new PairInt(-1, 0));
        ones.add(new PairInt(-1, -2));
        ones.add(new PairInt(0, -2));
        ones.add(new PairInt(1, 0));
        ones.add(new PairInt(1, -1));
        ones.add(new PairInt(1, -2));
        
        zeroes.add(new PairInt(-1, -1));
        zeroes.add(new PairInt(0, -1));        
        
        int nRotations = 3;
        
        replaceAndRotateOnesIfNullable(points, imageWidth, imageHeight,
            zeroes, ones, nRotations);
       
    }
    
    private void correctForHoleArtifacts1_3(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
        
        /*     
             0  0              2
          0  0  #  #           1
          0  #  0  0< #  0     0
                #  #* 0  0     -1
                   0  0        -2
        
         -3 -2 -1  0  1  2
        */  

        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();

        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-3, 0 - 1));
        zeroes.add(new PairInt(-3, -1 - 1));
        zeroes.add(new PairInt(-2, -1 - 1));
        zeroes.add(new PairInt(-2, -2 - 1));
        zeroes.add(new PairInt(-1, 0 - 1));
        zeroes.add(new PairInt(-1, -2 - 1));
        zeroes.add(new PairInt(0, 2 - 1));
        zeroes.add(new PairInt(1, 2 - 1));
        zeroes.add(new PairInt(1, 1 - 1));
        zeroes.add(new PairInt(2, 1 - 1));
        zeroes.add(new PairInt(2, 0 - 1));
        
        ones.add(new PairInt(-2, 0 - 1));
        ones.add(new PairInt(-1, 1 - 1));
        ones.add(new PairInt(-1, -1 - 1));
        ones.add(new PairInt(0, 1 - 1));
        ones.add(new PairInt(0, -1 - 1));
        ones.add(new PairInt(1, 0 - 1));
        
        changeToZeroes.add(new PairInt(-2, 0 - 1));
        changeToZeroes.add(new PairInt(-1, -1 - 1));
        changeToZeroes.add(new PairInt(-1, 1 - 1));
        changeToZeroes.add(new PairInt(0, -1 - 1));
        changeToZeroes.add(new PairInt(0, 1 - 1));
        changeToZeroes.add(new PairInt(1, 0 - 1));
        
        changeToOnes.add(new PairInt(-2, 1 - 1));
        changeToOnes.add(new PairInt(-1, 0 - 1));
        changeToOnes.add(new PairInt(0, 0 - 1));
        changeToOnes.add(new PairInt(1, -1 - 1));
                    
        replacePattern(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
       
    }

    private void correctForHoleArtifacts1_4(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
        
        /*     
                   0  0        3
                #  #  0  0     2
                #  0  #  0     1
             0  #  0< #        0
             0  0  #*          -1
                0  0           -2
        
         -3 -2 -1  0  1  2
        */

        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();

        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-2, 1 - 1));
        zeroes.add(new PairInt(-2, 0 - 1));
        zeroes.add(new PairInt(-1, 2 - 1));
        zeroes.add(new PairInt(-1, 1 - 1));
        zeroes.add(new PairInt(0, 2 - 1));
        zeroes.add(new PairInt(0, -1 - 1));
        zeroes.add(new PairInt(0, -3 - 1));
        zeroes.add(new PairInt(1, -2 - 1));
        zeroes.add(new PairInt(1, -3 - 1));
        zeroes.add(new PairInt(2, -1 - 1));
        zeroes.add(new PairInt(2, -2 - 1));

        ones.add(new PairInt(-1, 0 - 1));
        ones.add(new PairInt(-1, -1 - 1));
        ones.add(new PairInt(-1, -2 - 1));
        ones.add(new PairInt(0, 1 - 1));
        ones.add(new PairInt(0, -2 - 1));
        ones.add(new PairInt(1, 0 - 1));
        ones.add(new PairInt(1, -1 - 1));
     
        changeToZeroes.add(new PairInt(-1, 0 - 1));
        changeToZeroes.add(new PairInt(-1, -1 - 1));
        changeToZeroes.add(new PairInt(0, 1 - 1));
        changeToZeroes.add(new PairInt(0, -2 - 1));
        changeToZeroes.add(new PairInt(1, 0 - 1));
        changeToZeroes.add(new PairInt(1, -1 - 1));
        
        changeToOnes.add(new PairInt(0, -1 - 1));
        changeToOnes.add(new PairInt(0, 0 - 1));
        changeToOnes.add(new PairInt(1, 1 - 1));
                    
        replacePattern(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
       
    }

    private void correctForHoleArtifacts0(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
        
        /*
                #  #  #      2
                #  0  #      1
                #  #* #      0
                               
         -3 -2 -1  0  1  2
        
        if the pattern is found, 
            -- set the center 0's to '1'
            -- for each pixel in the open squares, test which 
               values can be nulled and set them
        */
                
        Set<PairInt> ones = new HashSet<PairInt>();
        Set<PairInt> zeroes = new HashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        ones.add(new PairInt(-1, 0));
        ones.add(new PairInt(-1, -1));
        ones.add(new PairInt(-1, -2));
        ones.add(new PairInt(0, -2));
        ones.add(new PairInt(1, 0));
        ones.add(new PairInt(1, -1));
        ones.add(new PairInt(1, -2));
                
        zeroes.add(new PairInt(0, -1));
        
        int nRotations = 0;
        
        replaceAndRotateOnesIfNullable(points, imageWidth, imageHeight,
            zeroes, ones, nRotations);
    }
    
    private void correctForHoleArtifacts0_4(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
        
        /*
                #  #         2
                #  0  #      1
                #  #* #      0
                               
         -3 -2 -1  0  1  2
        
        if the pattern is found, 
            -- set the center 0's to '1'
            -- for each pixel in the open squares, test which 
               values can be nulled and set them
        */
                
        Set<PairInt> ones = new HashSet<PairInt>();
        Set<PairInt> zeroes = new HashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        ones.add(new PairInt(-1, 0));
        ones.add(new PairInt(-1, -1));
        ones.add(new PairInt(-1, -2));
        ones.add(new PairInt(0, -2));
        ones.add(new PairInt(1, 0));
        ones.add(new PairInt(1, -1));
                
        zeroes.add(new PairInt(0, -1));
        
        int nRotations = 3;
        
        replaceAndRotateOnesIfNullable(points, imageWidth, imageHeight,
            zeroes, ones, nRotations);
    }
    
    private void correctForHoleArtifacts0_1(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
        
        /*
                #  #  #      2
             #  #  0  #      1
          #  #  0  #* #      0
          #  0  #  #        -1
          #  #  #           -2
         -3 -2 -1  0  1  2
        
        if the pattern is found, 
            -- set the center 0's to '1'
            -- for each pixel in the open squares, test which 
               values can be nulled and set them
        */
                
        Set<PairInt> ones = new HashSet<PairInt>();
        Set<PairInt> zeroes = new HashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        ones.add(new PairInt(-3, 2));  ones.add(new PairInt(-3, 1)); ones.add(new PairInt(-3, 0));
        ones.add(new PairInt(-2, 2));  ones.add(new PairInt(-2, 0)); ones.add(new PairInt(-2, -1));
        ones.add(new PairInt(-1, 2));  ones.add(new PairInt(-1, 1)); ones.add(new PairInt(-1, -1)); ones.add(new PairInt(-1, -2));
        ones.add(new PairInt(0, 1));  ones.add(new PairInt(0, -2));
        ones.add(new PairInt(1, 0));  ones.add(new PairInt(1, -1)); ones.add(new PairInt(1, -2));
        
        zeroes.add(new PairInt(-2, 1));  
        zeroes.add(new PairInt(-1, 0));
        zeroes.add(new PairInt(0, -1));
        
        int nRotations = 3;
        
        replaceAndRotateOnesIfNullable(points, imageWidth, imageHeight,
            zeroes, ones, nRotations);
    }

    private void correctForHoleArtifacts0_2(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
        
        /*
                #  #  #      2
             #  #  0  #      1
             #  0  #* #      0
             #  #  #        -1
                            -2
         -3 -2 -1  0  1  2
        
        if the pattern is found, 
            -- set the center 0's to '1'
            -- for each pixel in the open squares, test which 
               values can be nulled and set them
        */
                
        Set<PairInt> ones = new HashSet<PairInt>();
        Set<PairInt> zeroes = new HashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        ones.add(new PairInt(-2, 1));  ones.add(new PairInt(-2, 0)); ones.add(new PairInt(-2, -1));
        ones.add(new PairInt(-1, 1)); ones.add(new PairInt(-1, -1)); ones.add(new PairInt(-1, -2));
        ones.add(new PairInt(0, 1));  ones.add(new PairInt(0, -2));
        ones.add(new PairInt(1, 0));  ones.add(new PairInt(1, -1)); ones.add(new PairInt(1, -2));
        
        zeroes.add(new PairInt(-1, 0));
        zeroes.add(new PairInt(0, -1));
        
        int nRotations = 3;
        
        replaceAndRotateOnesIfNullable(points, imageWidth, imageHeight,
            zeroes, ones, nRotations);
    }

    private void replaceAndRotateOnesIfNullable(Set<PairInt> points, 
        int imageWidth, int imageHeight,
        Set<PairInt> zeroes, Set<PairInt> ones, int nRotations) {
        
        ErosionFilter erosionFilter = new ErosionFilter();
        
        int w = imageWidth;
        int h = imageHeight;
        
        for (int nRot = 0; nRot <= nRotations; nRot++) {
       
            switch(nRot) {
                case 1:
                    reverseXs(zeroes, ones);
                    break;
                case 2:
                    reverseYs(zeroes, ones);
                    break;
                case 3:
                    reverseXs(zeroes, ones);
                    break;
                default:
                    break;
            }
            
            Set<PairInt> tmpPointsRemoved = new HashSet<PairInt>();
            Set<PairInt> tmpPointsAdded = new HashSet<PairInt>();

            for (PairInt p : points) {

                int col = p.getX();
                int row = p.getY();

                // make sure the current point hasn't been added to tmpPointsRemoved
                boolean isNotPresent = tmpPointsRemoved.contains(p);
                if (isNotPresent) {
                    continue;
                }

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
                
                //change the zeroes to ones.  bounds have been checked
                for (PairInt p2 : zeroes) {
                    int x = col + p2.getX();
                    int y = row + p2.getY();
                    PairInt p3 = new PairInt(x, y);
                    tmpPointsAdded.add(p3);
                }

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
    }

    private void correctForHoleArtifacts0_3(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
        
        /*
                   #  #      2
                #  0  #      1
             #  0  #*        0
                #           -1
                            -2
         -3 -2 -1  0  1  2
        
        if the pattern is found, 
            -- set the center 0's to '1'
            -- for each pixel in the open squares, test which 
               values can be nulled and set them
        */
                
        Set<PairInt> ones = new HashSet<PairInt>();
        Set<PairInt> zeroes = new HashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        ones.add(new PairInt(-2, 0)); 
        ones.add(new PairInt(-1, 1)); ones.add(new PairInt(-1, -1));
        ones.add(new PairInt(0, -2));
        ones.add(new PairInt(1, -1)); ones.add(new PairInt(1, -2));
        
        zeroes.add(new PairInt(-1, 0));
        zeroes.add(new PairInt(0, -1));
        
        int nRotations = 3;
        
        replaceAndRotateOnesIfNullable(points, imageWidth, imageHeight,
            zeroes, ones, nRotations);
    }

    private void correctForHoleArtifacts1_1(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
        
        /*     
             0  0  0  0  #     2
             0  #  #  #        1
             0  #  #* 0        0
             0  #  0           -1
             #
         -3 -2 -1  0  1  2
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
        zeroes.add(new PairInt(-1, -2));
        zeroes.add(new PairInt(0, 1));
        zeroes.add(new PairInt(0, -2));
        zeroes.add(new PairInt(1, 0));
        zeroes.add(new PairInt(1, -2));
         
        ones.add(new PairInt(-2, 2));
        ones.add(new PairInt(-1, 1));
        ones.add(new PairInt(-1, 0));
        ones.add(new PairInt(-1, -1));
        ones.add(new PairInt(0, -1));
        ones.add(new PairInt(1, -1));
        ones.add(new PairInt(2, -2));
        
        changeToZeroes.add(new PairInt(-1, -1));
        changeToZeroes.add(new PairInt(-1, 0));
        changeToZeroes.add(new PairInt(0, -1));
                    
        replacePattern(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
       
    }

    private void correctForHoleArtifacts00_0(Set<PairInt> points, 
        int imageWidth, int imageHeight) {
        
        /*   
                            5                        5
             # # 0          4         # #            4
           # 0 # # 0        3       #   # .          3
           # # 0 # #        2       # # @ # #        2
           0 # # 0 #        1         . #   #        1
             0 # #*         0           # #*         0
                           -1                       -1
          -3-2-1 0 1 2             -3-2-1 0 1 2
        */ 
       
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
        
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-3, -1));
        zeroes.add(new PairInt(-2, 0)); zeroes.add(new PairInt(-2, -3));
        zeroes.add(new PairInt(-1, -2));
        zeroes.add(new PairInt(0, -1)); zeroes.add(new PairInt(0, -4));
        zeroes.add(new PairInt(1, -3)); 
        
        ones.add(new PairInt(-3, -2)); ones.add(new PairInt(-3, -3));
        ones.add(new PairInt(-2, -1)); ones.add(new PairInt(-2, -2)); ones.add(new PairInt(-2, -4));
        ones.add(new PairInt(-1, 0)); ones.add(new PairInt(-1, -1)); ones.add(new PairInt(-1, -3)); ones.add(new PairInt(-1, -4));
        ones.add(new PairInt(0, -2)); ones.add(new PairInt(0, -3));
        ones.add(new PairInt(1, -1)); ones.add(new PairInt(1, -2));
        
        changeToZeroes.add(new PairInt(-2, -1));
        changeToZeroes.add(new PairInt(0, -3));
        
        changeToOnes.add(new PairInt(-1, -2));
        
        replacePattern(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
    }

    private void correctForHoleArtifacts00_1(Set<PairInt> points, 
        int imageWidth, int imageHeight) {
                
        /*   
             # 0 0 0        4                        4
           0 # # # 0 0      3         # . .          3
           0 # 0 # # 0      2         . @ . .        2
           # # # 0 # 0      1         # . @ .        1
             0 # #*# 0      0           . . #        0
             0 0 0 #       -1                       -1
          -3-2-1 0 1 2             -3-2-1 0 1 2
        */ 
       
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
        
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-3, -2)); zeroes.add(new PairInt(-3, -3));
        zeroes.add(new PairInt(-2, 1)); zeroes.add(new PairInt(-2, 0));
        zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, -2)); zeroes.add(new PairInt(-1, -4));
        zeroes.add(new PairInt(0, 1)); zeroes.add(new PairInt(0, -1)); zeroes.add(new PairInt(0, -4));
        zeroes.add(new PairInt(1, -3)); zeroes.add(new PairInt(1, -4));
        zeroes.add(new PairInt(2, 0)); zeroes.add(new PairInt(2, -1));
        zeroes.add(new PairInt(2, -2)); zeroes.add(new PairInt(2, -3));
        
        ones.add(new PairInt(-3, -1));
        ones.add(new PairInt(-2, -1)); ones.add(new PairInt(-2, -2)); ones.add(new PairInt(-2, -3)); ones.add(new PairInt(-2, -4));
        ones.add(new PairInt(-1, 0)); ones.add(new PairInt(-1, -1)); ones.add(new PairInt(-1, -3));
        ones.add(new PairInt(0, -2)); ones.add(new PairInt(0, -3));
        ones.add(new PairInt(1, 1)); ones.add(new PairInt(1, 0)); ones.add(new PairInt(1, -1)); ones.add(new PairInt(1, -2));
        
        changeToZeroes.add(new PairInt(-2, -2));
        changeToZeroes.add(new PairInt(-1, 0)); changeToZeroes.add(new PairInt(-1, -1)); changeToZeroes.add(new PairInt(-1, -3));
        changeToZeroes.add(new PairInt(0, 0)); changeToZeroes.add(new PairInt(0, -2)); changeToZeroes.add(new PairInt(0, -3));
        changeToZeroes.add(new PairInt(1, -1)); changeToZeroes.add(new PairInt(1, -2));
        
        changeToOnes.add(new PairInt(-1, -2));
        changeToOnes.add(new PairInt(0, -1));
        
        replacePattern(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
    }
    
    private void correctForHoleArtifacts00_2(Set<PairInt> points, 
        int imageWidth, int imageHeight) {
                
        /*   
           # 0 0 0 0        4                        4
           0 # # # 0 #      3         # . .   #      3
           0 # 0 # # #      2         . @ . . #      2
           0 # # 0 #        1         . . @ #        1
             0 # #*         0           .            0
                           -1                       -1
          -3-2-1 0 1 2             -3-2-1 0 1 2
        */ 
       
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
        
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-3, -1)); zeroes.add(new PairInt(-3, -2)); zeroes.add(new PairInt(-3, -3));
        zeroes.add(new PairInt(-2, 0)); zeroes.add(new PairInt(-2, -4));
        zeroes.add(new PairInt(-1, -2)); zeroes.add(new PairInt(-1, -4));
        zeroes.add(new PairInt(0, -1)); zeroes.add(new PairInt(0, -4));
        zeroes.add(new PairInt(1, -3)); zeroes.add(new PairInt(1, -4));
        
        ones.add(new PairInt(-3, -4));
        ones.add(new PairInt(-2, -1)); ones.add(new PairInt(-2, -2)); ones.add(new PairInt(-2, -3));
        ones.add(new PairInt(-1, 0)); ones.add(new PairInt(-1, -1)); ones.add(new PairInt(-1, -3));
        ones.add(new PairInt(0, -2)); ones.add(new PairInt(0, -3));
        ones.add(new PairInt(1, -1)); ones.add(new PairInt(1, -2));
        ones.add(new PairInt(2, -2)); ones.add(new PairInt(2, -3));
        
        changeToZeroes.add(new PairInt(-2, -1)); changeToZeroes.add(new PairInt(-2, -2));
        changeToZeroes.add(new PairInt(-1, 0)); changeToZeroes.add(new PairInt(-1, -1)); changeToZeroes.add(new PairInt(-1, -3));
        changeToZeroes.add(new PairInt(0, -2)); changeToZeroes.add(new PairInt(0, -3));
        changeToZeroes.add(new PairInt(1, -2));
        
        changeToOnes.add(new PairInt(-1, -2));
        changeToOnes.add(new PairInt(0, -1));
        
        replacePattern(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
    }
    
    private void correctForHoleArtifacts00_3(Set<PairInt> points, 
        int imageWidth, int imageHeight) {

        /*   
               0 0 0        4                        4
             # # # 0 0      3         # . .          3
           0 # 0 # # 0      2         . @ . .        2
           0 # # 0 # 0      1         . . @ .        1
             0 # #*#        0           . . #        0
             0 0 0         -1                       -1
          -3-2-1 0 1 2             -3-2-1 0 1 2
        */ 
       
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
        
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-3, -1)); zeroes.add(new PairInt(-3, -2));
        zeroes.add(new PairInt(-2, 1)); zeroes.add(new PairInt(-2, 0));
        zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, -2)); zeroes.add(new PairInt(-1, -4));
        zeroes.add(new PairInt(0, 1)); zeroes.add(new PairInt(0, -1)); zeroes.add(new PairInt(0, -4));
        zeroes.add(new PairInt(1, -3)); zeroes.add(new PairInt(1, -4));
        zeroes.add(new PairInt(2, -1)); zeroes.add(new PairInt(2, -2)); zeroes.add(new PairInt(2, -3));
        
        ones.add(new PairInt(-2, -1)); ones.add(new PairInt(-2, -2)); ones.add(new PairInt(-2, -3));
        ones.add(new PairInt(-1, 0)); ones.add(new PairInt(-1, -1)); ones.add(new PairInt(-1, -3));
        ones.add(new PairInt(0, -2)); ones.add(new PairInt(0, -3));
        ones.add(new PairInt(1, 0)); ones.add(new PairInt(1, -1)); ones.add(new PairInt(1, -2));
        
        changeToZeroes.add(new PairInt(-2, -1)); changeToZeroes.add(new PairInt(-2, -2));
        changeToZeroes.add(new PairInt(-1, 0)); changeToZeroes.add(new PairInt(-1, -1)); changeToZeroes.add(new PairInt(-1, -3));
        changeToZeroes.add(new PairInt(0, 0)); changeToZeroes.add(new PairInt(0, -2)); changeToZeroes.add(new PairInt(0, -3));
        changeToZeroes.add(new PairInt(1, -1)); changeToZeroes.add(new PairInt(1, -2));
        
        changeToOnes.add(new PairInt(-1, -2));
        changeToOnes.add(new PairInt(0, -1));

        replacePattern(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
    }
    
    private void correctForHoleArtifacts00_4(Set<PairInt> points, 
        int imageWidth, int imageHeight) {

        /*   
               0 0 0 0 0    4                        4
             0 0 # # # 0    3             . . .      3
             0 # # 0 # #    2           . . @ # #    2
             0 # 0 # # 0    1           . @ . .      1
             # # #*# 0      0         # # .*.        0
               0 0 0       -1                       -1
          -3-2-1 0 1 2 3           -3-2-1 0 1 2
        */ 
       
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
        
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-2, -1)); zeroes.add(new PairInt(-2, -2)); zeroes.add(new PairInt(-2, -3));
        zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, -3)); zeroes.add(new PairInt(-1, -4));
        zeroes.add(new PairInt(0, 1)); zeroes.add(new PairInt(0, -1)); zeroes.add(new PairInt(0, -4));
        zeroes.add(new PairInt(1, 1)); zeroes.add(new PairInt(1, -2)); zeroes.add(new PairInt(1, -4));
        zeroes.add(new PairInt(2, 0)); zeroes.add(new PairInt(2, -4));
        zeroes.add(new PairInt(3, -1)); zeroes.add(new PairInt(3, -3)); zeroes.add(new PairInt(3, -4));
        
        ones.add(new PairInt(-2, 0));
        ones.add(new PairInt(-1, 0)); ones.add(new PairInt(-1, -1)); ones.add(new PairInt(-1, -2));
        ones.add(new PairInt(0, -2)); ones.add(new PairInt(0, -3));
        ones.add(new PairInt(1, 0)); ones.add(new PairInt(1, -1)); ones.add(new PairInt(1, -3));
        ones.add(new PairInt(2, -1)); ones.add(new PairInt(2, -2)); ones.add(new PairInt(2, -3));
        ones.add(new PairInt(3, -2));
        
        changeToZeroes.add(new PairInt(-1, -1)); changeToZeroes.add(new PairInt(-1, -2));
        changeToZeroes.add(new PairInt(0, 0)); changeToZeroes.add(new PairInt(0, -2)); changeToZeroes.add(new PairInt(0, -3));
        changeToZeroes.add(new PairInt(1, 0)); changeToZeroes.add(new PairInt(1, -1)); changeToZeroes.add(new PairInt(1, -3));
        changeToZeroes.add(new PairInt(2, -1)); changeToZeroes.add(new PairInt(2, -3));
        
        changeToOnes.add(new PairInt(0, -1));  
        changeToOnes.add(new PairInt(1, -2));
        
        replacePattern(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
    }
    
    private void correctForHoleArtifacts00_10(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
        /*
        TODO: redo this one to check for connections for the diagonal corners
        of the surrounding 8 pixels.
        Since this artifact occurs on diagonal lines mainly, we're looking to
        remove all pixels in the square except the connecting diagonals.
        (along w/ setting the center to '1'.
        */
        int w = imageWidth;
        int h = imageHeight;
        
        int[] neighborsCount = new int[8];
        int[] neighborsIdx = new int[8];
        
        /*
                   #         2
                #  0  #      1
                   #*        0
                               
         -3 -2 -1  0  1  2
        */
                
        Set<PairInt> ones = new HashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        ones.add(new PairInt(-1, -1));
        ones.add(new PairInt(0, -2));
        ones.add(new PairInt(1, -1));
        
        ErosionFilter erosionFilter = new ErosionFilter();
                        
        Set<PairInt> tmpPointsRemoved = new HashSet<PairInt>();
        Set<PairInt> tmpPointsAdded = new HashSet<PairInt>();

        for (PairInt p : points) {

            int col = p.getX();
            int row = p.getY();

            // make sure the current point hasn't been added to tmpPointsRemoved
            boolean isNotPresent = tmpPointsRemoved.contains(p);
            if (isNotPresent) {
                continue;
            }

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
            
            /*
            set the center value to zero.
            visit each of the 8 neighbors and count the number of neighbors
            it has.
            sort those by count.
            visit the smallest count neighbors (w/ count > 0) and test
               if nullable.
            */
            
            // set the center point
            tmpPointsAdded.add(new PairInt(col, row - 1));

            Arrays.fill(neighborsCount, 0);
            Arrays.fill(neighborsIdx, 0);
            for (int n1Idx = 0; n1Idx < eightNeighborsX.length; ++n1Idx) {
                int x1 = col + eightNeighborsX[n1Idx];
                int y1 = row + eightNeighborsY[n1Idx];                
                neighborsIdx[n1Idx] = n1Idx;
                if ((x1 < 0) || (y1 < 0) || (x1 > (w - 1)) || (y1 > (h - 1))) {
                    neighborsCount[n1Idx] = 0;
                    continue;
                }
                int count = 0;
                // count the neighbors of the neighbor
                for (int n2Idx = 0; n2Idx < eightNeighborsX.length; ++n2Idx) {
                    int x2 = x1 + eightNeighborsX[n2Idx];
                    int y2 = y1 + eightNeighborsY[n2Idx];                
                    if ((x2 < 0) || (y2 < 0) || (x2 > (w - 1)) || (y2 > (h - 1))) {
                        neighborsCount[n2Idx] = 0;
                    }
                    PairInt p3 = new PairInt(x2, y2);
                    if (!tmpPointsRemoved.contains(p3) && 
                        (tmpPointsAdded.contains(p3) || points.contains(p3))) {
                        count++;
                    }
                }
                neighborsCount[n1Idx] = count;
            }
            CountingSort.sort(neighborsCount, neighborsIdx, 8);
            
            for (int ii = 0; ii < neighborsIdx.length; ++ii) {
                int nIdx = neighborsIdx[ii];
                if (neighborsCount[nIdx] == 0) {
                    continue;
                }
                
                int x2 = col + eightNeighborsX[nIdx];
                int y2 = row + eightNeighborsY[nIdx];
                
                PairInt p3 = new PairInt(x2, y2);
                
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
}
