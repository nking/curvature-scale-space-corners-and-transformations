package algorithms.imageProcessing;

import algorithms.CountingSort;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.util.CornerArray;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
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
    
    protected GreyscaleImage edgeGuideImage = null;
    
    public void setEdgeGuideImage(final GreyscaleImage gradientXY) {
        this.edgeGuideImage = gradientXY;
    }
    
    public void correctForArtifacts(GreyscaleImage input) {
           
        if (edgeGuideImage != null) {
            if (input.getXRelativeOffset() != edgeGuideImage.getXRelativeOffset()) {
                throw new IllegalStateException(
                    "input and edgeGuideImage must have same x offsets");
            }
            if (input.getYRelativeOffset() != edgeGuideImage.getYRelativeOffset()) {
                throw new IllegalStateException(
                    "input and edgeGuideImage must have same y offsets");
            }
            if (input.getWidth() != edgeGuideImage.getWidth()) {
                throw new IllegalStateException(
                    "input and edgeGuideImage must have same widths");
            }
            if (input.getHeight() != edgeGuideImage.getHeight()) {
                throw new IllegalStateException(
                    "input and edgeGuideImage must have same heights");
            }
        }
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        final int w = input.getWidth();
        final int h = input.getHeight();
        
        Set<PairInt> points = imageProcessor.readNonZeroPixels(input);
          
        for (int i = 0; i < 1; i++) {
            straightenLines(points, w, h);
        }
        
        correctForArtifacts(points, w, h);
       
        for (int i = 0; i < 2; i++) {
            straightenLines(points, w, h);
        }
        
        correctForArtifacts(points, w, h);

        imageProcessor.writeAsBinaryToImage(input, points);

    }
 
    public void correctForArtifacts(Set<PairInt> points, int w, int h) {
     
        correctForHoleArtifacts1_2_2(points, w, h);
        correctForHoleArtifacts1_2_3(points, w, h);
        
        correctForSingleHole_01(points, w, h);
        correctForSingleHole_02(points, w, h);
        correctForSingleHole_03(points, w, h);
        correctForSingleHole_04(points, w, h);
        correctForSingleHole_05(points, w, h);
        correctForSingleHole_06(points, w, h);
        correctForSingleHole_07(points, w, h);
        correctForSingleHole_08(points, w, h);
        correctForSingleHole_09(points, w, h);
        correctForSingleHole_10(points, w, h);
        correctForSingleHole_11(points, w, h);
        correctForSingleHole_12(points, w, h);

        correctForZigZag000_00(points, w, h);
        correctForZigZag000_01(points, w, h);
        correctForZigZag000_02(points, w, h);
        correctForZigZag000_03(points, w, h);
        
        correctForZigZag000_04(points, w, h);
        correctForZigZag000_05(points, w, h);
        correctForZigZag000_06(points, w, h);
        correctForZigZag000_07(points, w, h);
        correctForZigZag000_08(points, w, h);
        correctForZigZag000_09(points, w, h);

        //correctForHoleArtifacts1_1(points, w, h);       
        correctForHoleArtifacts1_2(points, w, h);
        correctForHoleArtifacts1_2_1(points, w, h);
        correctForHoleArtifacts1_2_4(points, w, h);
        correctForHoleArtifacts1_2_5(points, w, h);
        correctForHoleArtifacts1_2_6(points, w, h);
        correctForHoleArtifacts1_2_7(points, w, h);
        correctForHoleArtifacts1_3(points, w, h);
        
        correctForZigZag00(points, w, h);
        correctForZigZag00_0(points, w, h);
        correctForZigZag00_0_00(points, w, h);
        correctForZigZag00_0_01(points, w, h);
        correctForZigZag00_0_0(points, w, h);
        correctForZigZag00_0_1(points, w, h);
        correctForZigZag00_0_2(points, w, h);
        correctForZigZag00_1(points, w, h);
        correctForZigZag00_2(points, w, h);
        correctForZigZag00_3(points, w, h);
        correctForZigZag00_5(points, w, h);
        correctForZigZag00_6(points, w, h);
        correctForZigZag00_7(points, w, h);
        correctForZigZag00_8(points, w, h);
        correctForZigZag00_9(points, w, h);
        correctForZigZag00_10(points, w, h);
                
        correctForZigZag00_4(points, w, h);
        
        correctForZigZag00_11(points, w, h);
        correctForZigZag00_12(points, w, h);
        correctForZigZag00_13(points, w, h);
      
        correctForTs(points, w, h);
        correctForTs_0(points, w, h);
        correctForTs_0_1(points, w, h);
        correctForTs_0_2(points, w, h);
        correctForTs_0_3(points, w, h);
        correctForTs_1(points, w, h);
        correctForTs_1_1(points, w, h);
                
        // TODO: revisit, not sure this is always an artifact:
        correctForLine0(points, w, h);
        
        // better edge extraction at the expense of unsharpening true corners:
        correctForLs_0(points, w, h);
        correctForLs_1(points, w, h);        
        correctForLs2_0(points, w, h);
        correctForLs2_1(points, w, h);
        
        correctForZigZag00_14(points, w, h);
        correctForZigZag00_15(points, w, h);
        correctForZigZag00_16(points, w, h);
            
        correctForHoleArtifacts00_10(points, w, h);
        
        correctForZigZag00_0_1(points, w, h);
       
        correctForCorrectionCreatedSquares(points, w, h);
      
        correctForZigZag000_04(points, w, h);
        correctForZigZag000_07(points, w, h);
        
        correctForLs_0(points, w, h);
        correctForLs_1(points, w, h);
        correctForLs_2(points, w, h); 
        correctForLs_3(points, w, h);
        correctForLs2_0(points, w, h);
        correctForLs2_1(points, w, h);
        
        correctForTs_0(points, w, h);
        correctForTs_0_1(points, w, h);
        correctForTs_0_2(points, w, h);
        correctForTs_0_3(points, w, h);
        correctForTs_1(points, w, h);
        correctForTs_1_1(points, w, h);
        
        correctForZigZag00_17(points, w, h);
        
        correctForExtCorner(points, w, h);
        
        //correctForRemaining(points, w, h);
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
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    private void correctForZigZag00_0_00(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
                   
        /*       
            0  0  0         2
            #  #  0         1
            0  #* 0         0
            0  #  #        -1
            0  0  0        -2
           -1  0  1  2  3        
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-1, 2)); zeroes.add(new PairInt(-1, 1)); 
        zeroes.add(new PairInt(-1, 0)); zeroes.add(new PairInt(-1, -2));
        zeroes.add(new PairInt(0, 2)); zeroes.add(new PairInt(0, -2));        
        zeroes.add(new PairInt(1, 2)); zeroes.add(new PairInt(1, 0)); 
        zeroes.add(new PairInt(1, -1)); zeroes.add(new PairInt(1, -2));
                
        ones.add(new PairInt(-1, -1));
        ones.add(new PairInt(0, 1)); ones.add(new PairInt(0, -1));
        ones.add(new PairInt(1, 1));
        
        changeToZeroes.add(new PairInt(0, -1));
        changeToZeroes.add(new PairInt(0, 1));
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        reverseXs(zeroes, ones, changeToZeroes, changeToOnes);

        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    private void correctForZigZag00_0_01(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
        
        /*       
            0  0  0  #  0     1
            0  #  #* #  0     0
            0  #  0  0  0    -1
           -2 -1  0  1  2   
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-2, 1)); zeroes.add(new PairInt(-2, 0)); zeroes.add(new PairInt(-2, -1));
        zeroes.add(new PairInt(-1, -1));
        zeroes.add(new PairInt(0, 1)); zeroes.add(new PairInt(0, -1));        
        zeroes.add(new PairInt(1, 1)); 
        zeroes.add(new PairInt(2, 1)); zeroes.add(new PairInt(2, 0)); zeroes.add(new PairInt(2, -1));
                
        ones.add(new PairInt(-1, 1)); ones.add(new PairInt(-1, 0));
        ones.add(new PairInt(1, 0));  ones.add(new PairInt(1, -1));
        
        changeToZeroes.add(new PairInt(-1, 0));
        changeToZeroes.add(new PairInt(1, 0));
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        reverseXs(zeroes, ones, changeToZeroes, changeToOnes);

        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    private void correctForZigZag00_0(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*       
               0  0     #   2
               0  #  #      1
            0  #* #  0      0
               #  0  0     -1
               #           -2
           -1  0  1  2  3        
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-1, 0));
        zeroes.add(new PairInt(0, -1)); zeroes.add(new PairInt(0, -2));
        zeroes.add(new PairInt(1, 1)); zeroes.add(new PairInt(1, -2));
        zeroes.add(new PairInt(2, 1)); zeroes.add(new PairInt(2, 0));
        
        ones.add(new PairInt(0, 2)); ones.add(new PairInt(0, 1));
        ones.add(new PairInt(1, 0)); ones.add(new PairInt(1, -1));
        ones.add(new PairInt(2, -1));
        ones.add(new PairInt(3, -2));
        
        changeToZeroes.add(new PairInt(1, 0));
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    private void correctForZigZag00_0_0(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*       
               0  0         2
               0  #  #  #   1
            0  #* #  0      0
               #  0  0     -1
               #           -2
           -1  0  1  2  3        
        */
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-1, 0));
        zeroes.add(new PairInt(0, -1)); zeroes.add(new PairInt(0, -2));
        zeroes.add(new PairInt(1, 1)); zeroes.add(new PairInt(1, -2));
        zeroes.add(new PairInt(2, 1)); zeroes.add(new PairInt(2, 0));
        
        ones.add(new PairInt(0, 2)); ones.add(new PairInt(0, 1));
        ones.add(new PairInt(1, 0)); ones.add(new PairInt(1, -1));
        ones.add(new PairInt(2, -1));
        ones.add(new PairInt(3, -1));
        
        changeToZeroes.add(new PairInt(1, 0));
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    private void correctForZigZag00_0_1(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*       
               0            2
               0  #  #      1
            0  #* #  0      0
            0  #  0  0     -1
               0           -2
           -1  0  1  2  3        
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, 0));
        zeroes.add(new PairInt(0, 2));  zeroes.add(new PairInt(0, -1)); zeroes.add(new PairInt(0, -2));
        zeroes.add(new PairInt(1, 1));
        zeroes.add(new PairInt(2, 1)); zeroes.add(new PairInt(2, 0));
        
        ones.add(new PairInt(0, 1));
        ones.add(new PairInt(1, 0)); ones.add(new PairInt(1, -1));
        ones.add(new PairInt(2, -1));
        
        changeToZeroes.add(new PairInt(1, 0));
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    private void correctForZigZag00_0_2(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*       
               0  0  0      2
               0  #< #      1
            #  #* #  0      0
            0  #  0        -1
                           -2
           -1  0  1  2  3        
        */
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-1, 1));
        zeroes.add(new PairInt(0, -1)); zeroes.add(new PairInt(0, -2));
        zeroes.add(new PairInt(1, 1)); zeroes.add(new PairInt(1, -2));
        zeroes.add(new PairInt(2, 0)); zeroes.add(new PairInt(2, -2));
        
        ones.add(new PairInt(-1, 0));
        ones.add(new PairInt(0, 1));
        ones.add(new PairInt(1, 0)); ones.add(new PairInt(1, -1));
        ones.add(new PairInt(2, -1));
        
        changeToZeroes.add(new PairInt(1, -1));
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    /**
     * note, this one diminishes corner to make an unambiguous thinned section
     * @param points
     * @param imageWidth
     * @param imageHeight 
     */
    private void correctForRemaining(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*       
            0  0  0         2
            0  #< #         1
            0  #* 0         0
                           -1
                           -2
           -1  0  1  2  3        
        */
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-1, 0)); zeroes.add(new PairInt(-1, -1)); zeroes.add(new PairInt(-1, -2));
        zeroes.add(new PairInt(0, -2));
        zeroes.add(new PairInt(1, 0)); zeroes.add(new PairInt(1, -2));
        
        ones.add(new PairInt(0, -1)); 
        ones.add(new PairInt(1, -1));
        
        changeToZeroes.add(new PairInt(0, -1));
                    
        int nCorrections = 0;
       
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
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
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
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
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
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
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    private void correctForZigZag00_11(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*       
                0  0  0         1
                0  #* #  #      0
                0  0  #  0     -1
            -2 -1  0  1  2  3  
        */
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, 0)); zeroes.add(new PairInt(-1, -1));
        zeroes.add(new PairInt(0, 1)); zeroes.add(new PairInt(0, -1));
        zeroes.add(new PairInt(1, -1));
        zeroes.add(new PairInt(2, 1));
        
        ones.add(new PairInt(1, 1)); ones.add(new PairInt(1, 0));
        ones.add(new PairInt(2, 0));
        
        changeToZeroes.add(new PairInt(0, 0));
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    private void correctForZigZag00_12(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
    
        /*       
                0  #  0         1
                0  #* #  #      0
                #  #  0        -1
                0  0  0        -2
            -2 -1  0  1  2  3  
        */
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-1, 2)); zeroes.add(new PairInt(-1, 0)); zeroes.add(new PairInt(-1, -1));
        zeroes.add(new PairInt(0, 2));
        zeroes.add(new PairInt(1, 2)); zeroes.add(new PairInt(1, 1)); zeroes.add(new PairInt(1, -1));
        
        ones.add(new PairInt(-1, 1));
        ones.add(new PairInt(0, 1)); ones.add(new PairInt(0, -1));
        ones.add(new PairInt(1, 0));
        ones.add(new PairInt(2, 0));
 
        changeToZeroes.add(new PairInt(0, 1));
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    private void correctForZigZag00_13(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
    
        /*      0  0  0         3
                0  #  0         2
                0  #  #  #      1
             #  #  #* 0         0
             0  #  0  0        -1
             0  0  0           -2
            -2 -1  0  1  2  3  
        */
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-2, 2)); zeroes.add(new PairInt(-2, 1));
        zeroes.add(new PairInt(-1, 2)); zeroes.add(new PairInt(-1, -1)); 
        zeroes.add(new PairInt(-1, -2)); zeroes.add(new PairInt(-1, -3));
        zeroes.add(new PairInt(0, 2)); zeroes.add(new PairInt(0, 1)); zeroes.add(new PairInt(0, -3));
        zeroes.add(new PairInt(1, 1)); zeroes.add(new PairInt(1, 0)); 
        zeroes.add(new PairInt(1, -2)); zeroes.add(new PairInt(1, -3));
        
        ones.add(new PairInt(-2, 0));
        ones.add(new PairInt(-1, 1)); ones.add(new PairInt(-1, 0));
        ones.add(new PairInt(0, -1)); ones.add(new PairInt(0, -2));
        ones.add(new PairInt(1, -1));
        ones.add(new PairInt(2, -1));
 
        changeToZeroes.add(new PairInt(-1, 1));
        changeToZeroes.add(new PairInt(0, -2));
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
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
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
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
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
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
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
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
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
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
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
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
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    private void correctForZigZag000_07(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*       
                0  0  0         2
                #  #< 0         1
                0  #  #*        0
                0  0  0        -1
            -3 -2 -1  0  1  2     
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        ones.add(new PairInt(-2, -1));
        ones.add(new PairInt(-1, 0)); ones.add(new PairInt(-1, -1));
        
        zeroes.add(new PairInt(-2, 1)); zeroes.add(new PairInt(-2, 0)); zeroes.add(new PairInt(-2, -2));
        zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, -2));
        zeroes.add(new PairInt(0, 1)); zeroes.add(new PairInt(0, -1)); zeroes.add(new PairInt(0, -2));
        
        changeToZeroes.add(new PairInt(-1, -1));
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    private void correctForZigZag000_06(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*       
           0  0  0        2
           #  #  0  0     1
           0  #< #* 0     0
           0  0  #  0    -1
                         -2
       -3 -2 -1  0  1  2     
        */ 
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        ones.add(new PairInt(-2, -1));
        ones.add(new PairInt(-1, 0)); ones.add(new PairInt(-1, -1));
        ones.add(new PairInt(0, 1));
        
        zeroes.add(new PairInt(-2, 1)); zeroes.add(new PairInt(-2, 0)); zeroes.add(new PairInt(-2, -2));
        zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, -2));
        zeroes.add(new PairInt(0, -1)); zeroes.add(new PairInt(0, -2));
        zeroes.add(new PairInt(1, 1)); zeroes.add(new PairInt(1, 0)); zeroes.add(new PairInt(1, -1));
        
        changeToZeroes.add(new PairInt(-1, 0));
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    private void correctForZigZag000_08(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*       
                #  0  0         1
                #  #*<0         0
                0  #  0        -1
        
         -3 -2 -1  0  1  2     
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        ones.add(new PairInt(-1, 0)); ones.add(new PairInt(-1, -1));
        ones.add(new PairInt(0, 1));
        
        zeroes.add(new PairInt(-1, 1));
        zeroes.add(new PairInt(0, -1));
        zeroes.add(new PairInt(1, 1)); zeroes.add(new PairInt(1, 0)); zeroes.add(new PairInt(1, -1));
        
        changeToZeroes.add(new PairInt(0, 0));
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
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
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
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
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
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
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    private void correctForZigZag00_17(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*       
             0  0  0            2
             0  #  0            1
             #  #  #* #         0
             0  0  #  0        -1
                0  0  0        -2
        
         -3 -2 -1  0  1  2     
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        ones.add(new PairInt(-2, 0));
        ones.add(new PairInt(-1, 0)); ones.add(new PairInt(-1, -1));
        ones.add(new PairInt(0, 1));
        ones.add(new PairInt(1, 0));
        
        zeroes.add(new PairInt(-2, 1)); zeroes.add(new PairInt(-2, -1)); zeroes.add(new PairInt(-2, -2));
        zeroes.add(new PairInt(-1, 2)); zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, -2));
        zeroes.add(new PairInt(0, 2)); zeroes.add(new PairInt(0, -1)); zeroes.add(new PairInt(0, -2));
        zeroes.add(new PairInt(1, 2)); zeroes.add(new PairInt(1, 1));
        
        changeToZeroes.add(new PairInt(-1, -1));
        changeToZeroes.add(new PairInt(0, 1));
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
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
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    private void correctForTs_0(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*       
             #  0  0  0         1
             0  #  #* #         0
                0  #  0        -1
                0  0  0        -2
            -2 -1  0  1  2     
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        ones.add(new PairInt(-2, -1));
        ones.add(new PairInt(-1, 0));
        ones.add(new PairInt(0, 1));
        ones.add(new PairInt(1, 0));
        
        zeroes.add(new PairInt(-2, 0));
        zeroes.add(new PairInt(-1, 2)); zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, -1));
        zeroes.add(new PairInt(0, 2)); zeroes.add(new PairInt(0, -1));
        zeroes.add(new PairInt(1, 2)); zeroes.add(new PairInt(1, 1)); zeroes.add(new PairInt(1, -1));
        
        changeToZeroes.add(new PairInt(0, 0));
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    private void correctForTs_0_1(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*       
                0  0  0  0      1
             #  #  #* #  #      0
                0  #  0        -1
                0  0  0        -2
            -2 -1  0  1  2     
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        ones.add(new PairInt(-2, 0));
        ones.add(new PairInt(-1, 0));
        ones.add(new PairInt(0, 1));
        ones.add(new PairInt(1, 0));
        ones.add(new PairInt(2, 0));
        
        zeroes.add(new PairInt(-1, 2)); zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, -1));
        zeroes.add(new PairInt(0, 2)); zeroes.add(new PairInt(0, -1));
        zeroes.add(new PairInt(1, 2)); zeroes.add(new PairInt(1, 1)); zeroes.add(new PairInt(1, -1));
        zeroes.add(new PairInt(2, -1));
        
        changeToZeroes.add(new PairInt(0, 1));
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    private void correctForTs_0_2(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*       
                0  0  0         1
             0  #  #* #         0
             #  0  #  0        -1
                0  0  0        -2
            -2 -1  0  1  2     
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        ones.add(new PairInt(-2, 1));
        ones.add(new PairInt(-1, 0));
        ones.add(new PairInt(0, 1));
        ones.add(new PairInt(1, 0));
        
        zeroes.add(new PairInt(-2, 0));
        zeroes.add(new PairInt(-1, 2)); zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, -1));
        zeroes.add(new PairInt(0, 2)); zeroes.add(new PairInt(0, -1));
        zeroes.add(new PairInt(1, 2)); zeroes.add(new PairInt(1, 1)); zeroes.add(new PairInt(1, -1));
        
        changeToZeroes.add(new PairInt(0, 1));
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    private void correctForTs_0_3(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*       
             #  0  0  0         1
             0  #  #* #         0
                0  #  0        -1
                0  #  0        -2
            -2 -1  0  1  2     
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        ones.add(new PairInt(-2, -1));
        ones.add(new PairInt(-1, 0));
        ones.add(new PairInt(0, 2)); ones.add(new PairInt(0, 1));
        ones.add(new PairInt(1, 0));
        
        zeroes.add(new PairInt(-2, 0));
        zeroes.add(new PairInt(-1, 2)); zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, -1));
        zeroes.add(new PairInt(0, -1));
        zeroes.add(new PairInt(1, 2)); zeroes.add(new PairInt(1, 1)); zeroes.add(new PairInt(1, -1));
        
        changeToZeroes.add(new PairInt(0, 0));
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    private void correctForTs_1(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*       
                                2
                0  #  0  0      1
                0  #* #  0      0
                0  #  0  0     -1
                #  0           -2
            -2 -1  0  1  2     
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        ones.add(new PairInt(-1, 2));
        ones.add(new PairInt(0, 1)); ones.add(new PairInt(0, -1));
        ones.add(new PairInt(1, 0));
        
        zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, 0)); zeroes.add(new PairInt(-1, -1));
        zeroes.add(new PairInt(0, 2));
        zeroes.add(new PairInt(1, 1)); zeroes.add(new PairInt(1, -1));
        zeroes.add(new PairInt(2, 1)); zeroes.add(new PairInt(2, 0)); zeroes.add(new PairInt(2, -1));
        
        changeToZeroes.add(new PairInt(0, 0));
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    private void correctForTs_1_1(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*       
                0  #            2
                0  #  0  0      1
                0  #* #  0      0
                0  #  0  0     -1
                0  #           -2
            -2 -1  0  1  2     
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        ones.add(new PairInt(0, 2)); ones.add(new PairInt(0, 1)); 
        ones.add(new PairInt(0, -1)); ones.add(new PairInt(0, -2));
        ones.add(new PairInt(1, 0));
        
        zeroes.add(new PairInt(-1, 2)); zeroes.add(new PairInt(-1, 1)); 
        zeroes.add(new PairInt(-1, 0)); zeroes.add(new PairInt(-1, -1)); zeroes.add(new PairInt(-1, -2));
        zeroes.add(new PairInt(1, 1)); zeroes.add(new PairInt(1, -1));
        zeroes.add(new PairInt(2, 1)); zeroes.add(new PairInt(2, 0)); zeroes.add(new PairInt(2, -1));
        
        changeToZeroes.add(new PairInt(1, 0));
        
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        // ----- change the sign of x to handle other direction -----
        reverseXs(zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += replacePattern(
            points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
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
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    /**
     * using this one unfortunately reduces the sharpness of real corners,
     * but it allows clear point ordering when extracting edges from the
     * line thinned image.
     */
    protected void correctForLs_0(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*          
        0  #  0         1
        0  #* 0         0
        0  #< #        -1
        0  0  0  #     -2
        
       -1  0  1  2 
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
    
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-1, 2)); zeroes.add(new PairInt(-1, 1)); 
        zeroes.add(new PairInt(-1, 0)); zeroes.add(new PairInt(-1, -1)); 
        zeroes.add(new PairInt(0, 2));
        zeroes.add(new PairInt(1, 2)); zeroes.add(new PairInt(1, 0)); zeroes.add(new PairInt(1, -1));
        
        ones.add(new PairInt(0, 1)); ones.add(new PairInt(0, -1));
        ones.add(new PairInt(1, 1));
        ones.add(new PairInt(2, 2));
        
        changeToZeroes.add(new PairInt(0, 1));
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    /**
     * using this one unfortunately reduces the sharpness of real corners,
     * but it allows clear point ordering when extracting edges from the
     * line thinned image.
     */
    protected void correctForLs_1(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*          
        0  #  0         1
        0  #* 0         0
        0  #< #        -1
        0  0  #        -2
        
       -1  0  1  2 
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
    
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-1, 2)); zeroes.add(new PairInt(-1, 1)); 
        zeroes.add(new PairInt(-1, 0)); zeroes.add(new PairInt(-1, -1)); 
        zeroes.add(new PairInt(0, 2));
        zeroes.add(new PairInt(1, 2)); zeroes.add(new PairInt(1, 0)); zeroes.add(new PairInt(1, -1));
        
        ones.add(new PairInt(0, 1)); ones.add(new PairInt(0, -1));
        ones.add(new PairInt(1, 1));
        ones.add(new PairInt(1, 2));
        
        changeToZeroes.add(new PairInt(0, 1));
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    protected void correctForLs_2(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*          
             0  #  0         1
             0  #* 0  0      0
             0  #  #< 0     -1
             #  0  0  0     -2
          #                 -3
         -2 -1  0  1  2 
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
    
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, 0)); zeroes.add(new PairInt(-1, -1)); 
        zeroes.add(new PairInt(0, 2));
        zeroes.add(new PairInt(1, 2)); zeroes.add(new PairInt(1, 0)); zeroes.add(new PairInt(1, -1));
        zeroes.add(new PairInt(2, 2)); zeroes.add(new PairInt(2, 1)); zeroes.add(new PairInt(2, 0)); 
        
        ones.add(new PairInt(-2, 3));
        ones.add(new PairInt(-1, 2));
        ones.add(new PairInt(0, 1)); ones.add(new PairInt(0, -1));
        ones.add(new PairInt(1, 1));
        
        changeToZeroes.add(new PairInt(1, 1));
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    protected void correctForLs_3(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*          
            0  0  0         1
            0  #* 0         0
            0  #< #  #     -1
            #  0  0  0     -2
         #                 -3       
        -2 -1  0  1  2 
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
    
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, 0)); zeroes.add(new PairInt(-1, -1)); 
        zeroes.add(new PairInt(0, 2)); zeroes.add(new PairInt(0, -1));
        zeroes.add(new PairInt(1, 2)); zeroes.add(new PairInt(1, 0)); zeroes.add(new PairInt(1, -1));
        zeroes.add(new PairInt(2, 2));
        
        ones.add(new PairInt(-2, 3));
        ones.add(new PairInt(-1, 2));
        ones.add(new PairInt(0, 1));
        ones.add(new PairInt(1, 1));
        ones.add(new PairInt(2, 1));
        
        changeToZeroes.add(new PairInt(0, 0));
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    protected void correctForLs2_0(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*        
        #               2
        0  #  0  0      1
        0  #*<#  #      0
        0  0  0  0     -1
        
       -1  0  1  2  3
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
   
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-1, -1)); zeroes.add(new PairInt(-1, 0)); zeroes.add(new PairInt(-1, 1));
        zeroes.add(new PairInt(0, 1));
        zeroes.add(new PairInt(1, 1)); zeroes.add(new PairInt(1, -1));
        zeroes.add(new PairInt(2, 1)); zeroes.add(new PairInt(2, -1));
        
        ones.add(new PairInt(-1, -2));
        ones.add(new PairInt(0, -1));
        ones.add(new PairInt(1, 0));
        ones.add(new PairInt(2, 0));
        
        changeToZeroes.add(new PairInt(0, 0));
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    protected void correctForLs2_1(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*        
           #            2
        0  #  0  0      1
        0  #*<#  #      0
        0  0  0  0     -1
        
       -1  0  1  2  3
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
   
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-1, -1)); zeroes.add(new PairInt(-1, 0)); zeroes.add(new PairInt(-1, 1));
        zeroes.add(new PairInt(0, 1));
        zeroes.add(new PairInt(1, 1)); zeroes.add(new PairInt(1, -1));
        zeroes.add(new PairInt(2, 1)); zeroes.add(new PairInt(2, -1));
        
        ones.add(new PairInt(0, -2));
        ones.add(new PairInt(0, -1));
        ones.add(new PairInt(1, 0));
        ones.add(new PairInt(2, 0));
        
        changeToZeroes.add(new PairInt(0, 0));
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    protected void correctForExtCorner(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*        
           0  0  #      2
        0  0  #  0      1
        0  #*<#  0      0
        0  0  0  0     -1
        
       -1  0  1  2  3
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
   
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, 0)); zeroes.add(new PairInt(-1, 1));
        zeroes.add(new PairInt(0, 1)); zeroes.add(new PairInt(0, -1)); zeroes.add(new PairInt(0, -2));
        zeroes.add(new PairInt(1, 1)); zeroes.add(new PairInt(1, -2));
        zeroes.add(new PairInt(2, 1)); zeroes.add(new PairInt(2, 0)); zeroes.add(new PairInt(2, -1));
        
        ones.add(new PairInt(1, 0)); ones.add(new PairInt(1, -1));
        ones.add(new PairInt(2, -2));
        
        changeToZeroes.add(new PairInt(1, 0));
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    private void correctForZigZag000_09(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*        
            #  0  0         1
            #  #*<0         0
            0  #  #        -1
           -1  0  1  2 
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
    
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-1, 1));
        zeroes.add(new PairInt(0, -1));
        zeroes.add(new PairInt(1, 0)); zeroes.add(new PairInt(1, -1));
        
        ones.add(new PairInt(-1, 0)); ones.add(new PairInt(-1, -1));
        ones.add(new PairInt(0, 1));
        ones.add(new PairInt(1, -1));
        
        changeToZeroes.add(new PairInt(0, 0));
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    protected void reverseXs(
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
    
    protected void reverseXs(final Set<PairInt> zeroes, final Set<PairInt> ones) {
        
        // ----- change the sign of x to handle other direction -----
        for (PairInt p : zeroes) {
            p.setX(-1 * p.getX());
        }
        
        for (PairInt p : ones) {
            p.setX(-1 * p.getX());
        }
        
    }
    
    protected void reverseYs(final Set<PairInt> zeroes, final Set<PairInt> ones, 
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
    
    protected void reverseYs(final Set<PairInt> zeroes, final Set<PairInt> ones) {
        
        // ----- change the sign of y  -----
        for (PairInt p : zeroes) {
            p.setY(-1 * p.getY());
        }
        
        for (PairInt p : ones) {
            p.setY(-1 * p.getY());
        }
        
    }
    
    private int rotate90ThreeTimes(
        Set<PairInt> points, int imageWidth, int imageHeight,
        final LinkedHashSet<PairInt> zeroes, final LinkedHashSet<PairInt> ones, 
        LinkedHashSet<PairInt> changeToZeroes, 
        final LinkedHashSet<PairInt> changeToOnes) {
        
        // ----- change the sign of x to handle other direction -----
        reverseXs(zeroes, ones, changeToZeroes, changeToOnes);
        
        int nCorrections = 0;
        
        nCorrections += replacePattern(
            points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
             
        // ----- change the sign of y to handle other direction -----
        reverseYs(zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += replacePattern(
            points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        // ----- change the sign of x to handle another direction -----
        reverseXs(zeroes, ones, changeToZeroes, changeToOnes);
                    
        nCorrections += replacePattern(
            points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        return nCorrections;
    }
    
    private int replacePattern(
        Set<PairInt> points, int imageWidth, int imageHeight,
        final LinkedHashSet<PairInt> zeroes, final LinkedHashSet<PairInt> ones, 
        final LinkedHashSet<PairInt> changeToZeroes, 
        final LinkedHashSet<PairInt> changeToOnes) {
        
        int w = imageWidth;
        int h = imageHeight;

        Set<PairInt> tmpPointsRemoved = new HashSet<PairInt>();
        Set<PairInt> tmpPointsAdded = new HashSet<PairInt>();
        
        for (PairInt p : points) {

            boolean isRemoved = tmpPointsRemoved.contains(p);
            if (isRemoved) {
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
                PairInt p3 = new PairInt(x, y);
                tmpPointsAdded.remove(p3);              
                tmpPointsRemoved.add(p3);
            }

            for (PairInt p2 : changeToOnes) {
                int x = col + p2.getX();
                int y = row + p2.getY();
                if ((x < 0) || (y < 0) || (x > (w - 1)) || (y > (h - 1))) {
                    continue;
                }
                PairInt p3 = new PairInt(x, y);                 
                tmpPointsRemoved.remove(p3);
                tmpPointsAdded.add(p3);
            }
        }
        
        int nCorrections = tmpPointsRemoved.size() + tmpPointsAdded.size();
        
        for (PairInt p2 : tmpPointsRemoved) {
            points.remove(p2);
        }
        for (PairInt p2 : tmpPointsAdded) {
            points.add(p2);
        }
        
        return nCorrections;
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
                      
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(
            points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
       
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
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
                
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(
            points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
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
                
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(
            points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
       
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
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
        
        int nCorrections = replaceAndRotateOnesIfNullable(points, imageWidth, 
            imageHeight, zeroes, ones, nRotations);
       
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    private void correctForHoleArtifacts1_2_4(Set<PairInt> points, int imageWidth,
        int imageHeight) {
        
        /*     
                #  #  #        2
             #  0  0  #        1
             #  #  #*          0
        
         -3 -2 -1  0  1  2
        */

        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();

        // y's are inverted here because sketch above is top left is (0, 0)
        ones.add(new PairInt(-2, 0)); ones.add(new PairInt(-2, -1));
        ones.add(new PairInt(-1, 0)); ones.add(new PairInt(-1, -2));
        ones.add(new PairInt(0, -2));
        ones.add(new PairInt(1, -1)); ones.add(new PairInt(1, -2));
        
        zeroes.add(new PairInt(-1, -1));
        zeroes.add(new PairInt(0, -1));
        
        int nRotations = 1;
        
        int nCorrections = replaceAndRotateOnesIfNullable(points, imageWidth, 
            imageHeight, zeroes, ones, nRotations);
       
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    private void correctForHoleArtifacts1_2_5(Set<PairInt> points, int imageWidth,
        int imageHeight) {
        
        /*     
                #  #  #        2
             #  0  0  #        1
                #  #*          0
        
         -3 -2 -1  0  1  2
        */

        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();

        // y's are inverted here because sketch above is top left is (0, 0)
        ones.add(new PairInt(-2, -1));
        ones.add(new PairInt(-1, 0)); ones.add(new PairInt(-1, -2));
        ones.add(new PairInt(0, -2));
        ones.add(new PairInt(1, -1)); ones.add(new PairInt(1, -2));
        
        zeroes.add(new PairInt(-1, -1));
        zeroes.add(new PairInt(0, -1));
        
        int nRotations = 1;
        
        int nCorrections = replaceAndRotateOnesIfNullable(points, imageWidth, 
            imageHeight, zeroes, ones, nRotations);
       
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    private void correctForHoleArtifacts1_2_6(Set<PairInt> points, int imageWidth,
        int imageHeight) {
        
        /*     
                #  #  #        3
                #  0  #        2
                #  0  #        1
                   #* #        0
        
         -3 -2 -1  0  1  2
        */

        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();

        // y's are inverted here because sketch above is top left is (0, 0)
        ones.add(new PairInt(-1, -1)); ones.add(new PairInt(-1, -2)); ones.add(new PairInt(-1, -3));
        ones.add(new PairInt(0, -3));
        ones.add(new PairInt(1, 0)); ones.add(new PairInt(1, -1));
        ones.add(new PairInt(1, -2)); ones.add(new PairInt(1, -3));
        
        zeroes.add(new PairInt(0, -1)); zeroes.add(new PairInt(0, -2));
        
        int nRotations = 1;
        
        int nCorrections = replaceAndRotateOnesIfNullable(points, imageWidth, 
            imageHeight, zeroes, ones, nRotations);
       
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    private void correctForHoleArtifacts1_2_7(Set<PairInt> points, int imageWidth,
        int imageHeight) {
        
        /*     
                   #  #        3
                #  0  #        2
                #  0  #        1
                   #* #        0
        
         -3 -2 -1  0  1  2
        */

        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();

        // y's are inverted here because sketch above is top left is (0, 0)
        ones.add(new PairInt(-1, -1)); ones.add(new PairInt(-1, -2));
        ones.add(new PairInt(0, -3));
        ones.add(new PairInt(1, 0)); ones.add(new PairInt(1, -1));
        ones.add(new PairInt(1, -2)); ones.add(new PairInt(1, -3));
        
        zeroes.add(new PairInt(0, -1)); zeroes.add(new PairInt(0, -2));
        
        int nRotations = 1;
        
        int nCorrections = replaceAndRotateOnesIfNullable(points, imageWidth, 
            imageHeight, zeroes, ones, nRotations);
       
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    private void correctForHoleArtifacts1_2_2(Set<PairInt> points, int imageWidth,
        int imageHeight) {
        
        /*   0  0  0  0  0     4                        4
             0  #  #  #  0     3         .  .  .        3 
             0  #  0  #  #     2         .     #  #     2
             0  #  0  #  0     1         .  @  .        1
             0  #  #* #  0     0         #  .* .        0
                #  0  0  0     -1        #              -1
        
         -3 -2 -1  0  1  2        -3 -2 -1  0  1  2
        */
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();

        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-2, 0)); zeroes.add(new PairInt(-2, -1));
        zeroes.add(new PairInt(-2, -2));  zeroes.add(new PairInt(-2, -3)); zeroes.add(new PairInt(-2, -4));
        zeroes.add(new PairInt(-1, -4));
        zeroes.add(new PairInt(0, 1)); zeroes.add(new PairInt(0, -1)); 
        zeroes.add(new PairInt(0, -2)); zeroes.add(new PairInt(0, -4));
        zeroes.add(new PairInt(1, 1)); zeroes.add(new PairInt(1, -4));
        zeroes.add(new PairInt(2, 1)); zeroes.add(new PairInt(2, 0)); zeroes.add(new PairInt(2, -1));
        zeroes.add(new PairInt(2, -3)); zeroes.add(new PairInt(2, -4));
        
        ones.add(new PairInt(-1, 1)); ones.add(new PairInt(-1, 0)); ones.add(new PairInt(-1, -1));
        ones.add(new PairInt(-1, -2)); ones.add(new PairInt(-1, -3));
        ones.add(new PairInt(0, -3));
        ones.add(new PairInt(1, 0)); ones.add(new PairInt(1, -1)); ones.add(new PairInt(1, -2));
        ones.add(new PairInt(1, -3));
        ones.add(new PairInt(2, -2));
        
        changeToOnes.add(new PairInt(0, -1));
        
        changeToZeroes.add(new PairInt(-1, -1)); changeToZeroes.add(new PairInt(-1, -2));
        changeToZeroes.add(new PairInt(-1, -3));
        changeToZeroes.add(new PairInt(0, 0)); changeToZeroes.add(new PairInt(0, -3));
        changeToZeroes.add(new PairInt(1, 0)); changeToZeroes.add(new PairInt(1, -1));
        changeToZeroes.add(new PairInt(1, -3));
                    
        int nCorrections = 0;
      
        nCorrections += replacePattern(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
       
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
       
    }
    
    private void correctForHoleArtifacts1_2_3(Set<PairInt> points, int imageWidth,
        int imageHeight) {
        
        /*                                               
           0  0  0  0  0  0  0     3                                 3
           0  #  #  #  #  #  0     2            .  .  .  .  .        2
           0  #  0  0  0  #  0     1            .     @  @  .        1
           #  #  #  #  #* #  0     0         #  #  #  .  .* .        0
              0  #  0  0  0  0     -1           0  #                 -1

          -4 -3 -2 -1  0  1  2        
        */
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();

        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-4, -1)); zeroes.add(new PairInt(-4, -2)); zeroes.add(new PairInt(-4, -3));
        zeroes.add(new PairInt(-3, 1)); zeroes.add(new PairInt(-3, -3));
        zeroes.add(new PairInt(-2, -1)); zeroes.add(new PairInt(-2, -3));
        zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, -1)); zeroes.add(new PairInt(-1, -3));
        zeroes.add(new PairInt(0, 1)); zeroes.add(new PairInt(0, -1)); zeroes.add(new PairInt(0, -3)); 
        zeroes.add(new PairInt(1, 1)); zeroes.add(new PairInt(1, -3));
        zeroes.add(new PairInt(2, 1)); zeroes.add(new PairInt(2, 0)); zeroes.add(new PairInt(2, -1));
        zeroes.add(new PairInt(2, -2)); zeroes.add(new PairInt(2, -3));
        
        ones.add(new PairInt(-4, 0)); 
        ones.add(new PairInt(-3, 0)); ones.add(new PairInt(-3, -1)); ones.add(new PairInt(-3, -2));
        ones.add(new PairInt(-2, 1)); ones.add(new PairInt(-2, 0)); ones.add(new PairInt(-2, -2));
        ones.add(new PairInt(-1, 0)); ones.add(new PairInt(-1, -2));
        ones.add(new PairInt(0, -2));
        ones.add(new PairInt(1, 0)); ones.add(new PairInt(1, -1)); ones.add(new PairInt(1, -2));
      
        changeToOnes.add(new PairInt(-1, -1));
        changeToOnes.add(new PairInt(0, -1));
        
        changeToZeroes.add(new PairInt(-3, -1)); changeToZeroes.add(new PairInt(-3, -2));
        changeToZeroes.add(new PairInt(-2, -2));
        changeToZeroes.add(new PairInt(-1, 0)); changeToZeroes.add(new PairInt(-1, -2));
        changeToZeroes.add(new PairInt(0, 0)); changeToZeroes.add(new PairInt(0, -2));
        changeToZeroes.add(new PairInt(1, 0)); changeToZeroes.add(new PairInt(1, -1));
        changeToZeroes.add(new PairInt(1, -2));
        
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
       
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
       
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
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
       
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }

    private int replaceAndRotateOnesIfNullable(Set<PairInt> points, 
        int imageWidth, int imageHeight,
        Set<PairInt> zeroes, Set<PairInt> ones, int nRotations) {
        
        ErosionFilter erosionFilter = new ErosionFilter();
        erosionFilter.overrideEndOfLineCheck();
        
        int w = imageWidth;
        int h = imageHeight;
        
        int[] neighborsCount = new int[ones.size()];
        PairInt[] neighbors = new PairInt[neighborsCount.length];
        
        int nCorrected = 0;
        
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
                if (tmpPointsRemoved.contains(p)) {
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
                    tmpPointsRemoved.remove(p3);
                    tmpPointsAdded.add(p3);
                }
               
                Arrays.fill(neighborsCount, 0);
                int nIdx = 0;
                for (PairInt p2 : ones) {
                    int x1 = col + p2.getX();
                    int y1 = row + p2.getY();
                    PairInt p3 = new PairInt(x1, y1);
                    neighbors[nIdx] = p3;
                    int count = 0;
                    for (int n2Idx = 0; n2Idx < eightNeighborsX.length; ++n2Idx) {
                        int x2 = x1 + eightNeighborsX[n2Idx];
                        int y2 = y1 + eightNeighborsY[n2Idx];
                        if ((x2 < 0) || (y2 < 0) || (x2 > (w - 1)) || (y2 > (h - 1))) {
                            continue;
                        }
                        PairInt p4 = new PairInt(x2, y2);
                        if (!tmpPointsRemoved.contains(p4)
                            && (tmpPointsAdded.contains(p4) || 
                            points.contains(p4))) {
                            count++;
                        }
                    }
                    neighborsCount[nIdx] = count;
                    nIdx++;
                }
                
                CountingSort.sort(neighborsCount, neighbors, 8);

                for (nIdx = 0; nIdx < neighbors.length; ++nIdx) {
                    
                    PairInt p3 = neighbors[nIdx];
                    
                    // test if can set to 0's without disconnecting lines

                    // adds to tmpPointsRemoved
                    boolean nullable = erosionFilter.process(p3, points, 
                        tmpPointsAdded, tmpPointsRemoved, w, h);
                }
            }

            nCorrected += (tmpPointsRemoved.size() + tmpPointsAdded.size());
            
            for (PairInt p2 : tmpPointsRemoved) {
                points.remove(p2);
            }
            for (PairInt p2 : tmpPointsAdded) {
                points.add(p2);
            }
        }
        
        return nCorrected;
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
           
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
       
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }

    private void correctForZigZag000_00(Set<PairInt> points, 
        int imageWidth, int imageHeight) {
        
        /* diagonal lines on both sides of a narrow cavity are sometimes
           zig-zag and can be thinned in the direction of widening the
           cavity.
                                3                          3
                   0 #          2             0 #          2
                 # # 0 #        1           # . 0 #        1
                 0 # #*0 #      0             # .*0        0
                   0 #         -1               #         -1
        
        -6-5-4-3-2-1 0 1 2 3       -6-5-4-3-2-1 0 1 2 3
        */
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
        
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-2, 0));
        zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, -2));
        zeroes.add(new PairInt(0, -1));
        zeroes.add(new PairInt(1, 0));
        
        ones.add(new PairInt(-2, -1));
        ones.add(new PairInt(-1, 0)); ones.add(new PairInt(-1, -1));
        ones.add(new PairInt(0, 1)); ones.add(new PairInt(0, -2));
        ones.add(new PairInt(1, -1));
        ones.add(new PairInt(2, 0));
        
        changeToZeroes.add(new PairInt(-1, -1)); 
        changeToZeroes.add(new PairInt(0, 0));
                
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        // ----- change the sign of x to handle other direction -----
        reverseXs(zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += replacePattern(
            points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }

    private void correctForZigZag000_01(Set<PairInt> points, 
        int imageWidth, int imageHeight) {
        
        /* diagonal lines on both sides of a narrow cavity are sometimes
           zig-zag and can be thinned in the direction of widening the
           cavity.
                                3                          3
                     # #        2               #          2
                   # 0 # #      1             # 0 .        1
                     #*0 # #    0               #*0 .      0
                       #       -1                 #       -1
        
        -6-5-4-3-2-1 0 1 2 3       -6-5-4-3-2-1 0 1 2 3
        */
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
        
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(0, -1));
        zeroes.add(new PairInt(1, 0));
        
        ones.add(new PairInt(-1, -1));
        ones.add(new PairInt(0, -2));
        ones.add(new PairInt(1, 1)); ones.add(new PairInt(1, -1)); ones.add(new PairInt(1, -2));
        ones.add(new PairInt(2, 0)); ones.add(new PairInt(2, -1));
        ones.add(new PairInt(3, 0));
        
        changeToZeroes.add(new PairInt(1, -1)); 
        changeToZeroes.add(new PairInt(2, 0));
                
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        // ----- change the sign of x to handle other direction -----
        reverseXs(zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += replacePattern(
            points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    private void correctForZigZag000_02(Set<PairInt> points, 
        int imageWidth, int imageHeight) {
        
        /* diagonal lines on both sides of a narrow cavity are sometimes
           zig-zag and can be thinned in the direction of widening the
           cavity.
                                3                          3
                     # #        2               #          2
                   # 0 # #      1             # 0 .        1
                     #*0 #      0               #*0        0
                       #       -1                 #       -1
        
        -6-5-4-3-2-1 0 1 2 3       -6-5-4-3-2-1 0 1 2 3
        */
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
        
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(0, -1));
        zeroes.add(new PairInt(1, 0));
        
        ones.add(new PairInt(-1, -1));
        ones.add(new PairInt(0, -2));
        ones.add(new PairInt(1, 1)); ones.add(new PairInt(1, -1)); ones.add(new PairInt(1, -2));
        ones.add(new PairInt(2, 0)); ones.add(new PairInt(2, -1));
        
        changeToZeroes.add(new PairInt(1, -1)); 
                
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        // ----- change the sign of x to handle other direction -----
        reverseXs(zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += replacePattern(
            points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    private void correctForZigZag000_04(Set<PairInt> points, 
        int imageWidth, int imageHeight) {
        
        /*
                                3                           3
               # 0 0 0          2         # 0 0 0          2
               0 # # 0 0        1         0 # #<0 0        1
                 0 # #*0        0           0 # #<0        0
                     0 0       -1               0 0       -1
        
        -6-5-4-3-2-1 0 1 2 3       -6-5-4-3-2-1 0 1 2 3
        */
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
        
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-3, -1));
        zeroes.add(new PairInt(-2, 0)); zeroes.add(new PairInt(-2, -2));
        zeroes.add(new PairInt(-1, -2));
        zeroes.add(new PairInt(0, 1)); zeroes.add(new PairInt(0, -1)); zeroes.add(new PairInt(0, -2));
        zeroes.add(new PairInt(1, 1)); zeroes.add(new PairInt(1, 0)); zeroes.add(new PairInt(1, -1));
        
        ones.add(new PairInt(-3, -2));
        ones.add(new PairInt(-2, -1));
        ones.add(new PairInt(-1, 0)); ones.add(new PairInt(-1, -1));
        
        changeToZeroes.add(new PairInt(0, 0));
        changeToZeroes.add(new PairInt(-1, -1));
                
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        // ----- change the sign of x to handle other direction -----
        reverseXs(zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += replacePattern(
            points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    private void correctForZigZag000_03(Set<PairInt> points, 
        int imageWidth, int imageHeight) {
        
        /* diagonal lines on both sides of a narrow cavity are sometimes
           zig-zag and can be thinned in the direction of widening the
           cavity.
                                3                          3
                   0 # #        2             0 #          2
                 # # 0 # #      1           # . 0 #        1
                 # # #*0 #      0             # .*0        0
                   # #         -1               #         -1
        
        -6-5-4-3-2-1 0 1 2 3       -6-5-4-3-2-1 0 1 2 3
        */
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
        
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-1, -2));
        zeroes.add(new PairInt(0, -1));
        zeroes.add(new PairInt(1, 0));
        
        ones.add(new PairInt(-2, -1)); ones.add(new PairInt(-2, 0));
        ones.add(new PairInt(-1, 1)); ones.add(new PairInt(-1, 0)); ones.add(new PairInt(-1, -1));
        ones.add(new PairInt(0, 1)); ones.add(new PairInt(0, -2));
        ones.add(new PairInt(1, -1)); ones.add(new PairInt(1, -2));
        ones.add(new PairInt(2, 0)); ones.add(new PairInt(2, -1));
        
        changeToZeroes.add(new PairInt(-1, -1)); 
        changeToZeroes.add(new PairInt(0, 0));
                
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        // ----- change the sign of x to handle other direction -----
        reverseXs(zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += replacePattern(
            points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }

    private void correctForZigZag000_05(Set<PairInt> points, 
        int imageWidth, int imageHeight) {
        
        /* diagonal lines on both sides of a narrow cavity are sometimes
           zig-zag and can be thinned in the direction of widening the
           cavity.
                                3                          3
                 # 0 0 0        2           #              2
               0 #.# 0 0 0      1           #.#            1
               0 0 #.#*0 0      0             #.#*         0
                 0 0 #         -1               #         -1
        
        -6-5-4-3-2-1 0 1 2 3       -6-5-4-3-2-1 0 1 2 3
        */
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
        
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-3, 0)); zeroes.add(new PairInt(-3, -1));
        zeroes.add(new PairInt(-2, 1)); zeroes.add(new PairInt(-2, 0));
        zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, -2));
        zeroes.add(new PairInt(0, -1)); zeroes.add(new PairInt(0, -2));
        zeroes.add(new PairInt(1, 0)); zeroes.add(new PairInt(1, -1)); zeroes.add(new PairInt(1, -2));
        zeroes.add(new PairInt(2, 0)); zeroes.add(new PairInt(2, -1));
        
        ones.add(new PairInt(-2, -1)); ones.add(new PairInt(-2, -2));
        ones.add(new PairInt(-1, 0)); ones.add(new PairInt(-1, -1));
        ones.add(new PairInt(0, 1)); ones.add(new PairInt(0, 0));
        
        changeToZeroes.add(new PairInt(-2, -1)); 
        changeToZeroes.add(new PairInt(-1, 0));
         
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        // ----- change the sign of x to handle other direction -----
        reverseXs(zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += replacePattern(
            points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    public static int correctForHoleArtifacts00_10(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
        
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
        erosionFilter.overrideEndOfLineCheck();
                        
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
            
            int centerRow = row - 1;
            
            // set the center point
            PairInt p3 = new PairInt(col, centerRow);
            tmpPointsRemoved.remove(p3);
            tmpPointsAdded.add(p3);

            Arrays.fill(neighborsCount, 0);
            Arrays.fill(neighborsIdx, 0);
            for (int n1Idx = 0; n1Idx < eightNeighborsX.length; ++n1Idx) {
                int x1 = col + eightNeighborsX[n1Idx];
                int y1 = centerRow + eightNeighborsY[n1Idx];                
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
                        continue;
                    }
                    PairInt p4 = new PairInt(x2, y2);
                    if (!tmpPointsRemoved.contains(p4) && 
                        (tmpPointsAdded.contains(p4) || points.contains(p4))) {
                        count++;
                    }
                }
                neighborsCount[n1Idx] = count;
            }
            CountingSort.sort(neighborsCount, neighborsIdx, 8);

            for (int ii = 0; ii < neighborsIdx.length; ++ii) {
                
                if (neighborsCount[ii] == 0) {
                    continue;
                }
                
                int nIdx = neighborsIdx[ii];
                
                int x2 = col + eightNeighborsX[nIdx];
                int y2 = centerRow + eightNeighborsY[nIdx];
                
                PairInt p4 = new PairInt(x2, y2);

                // adds to tmpPointsRemoved
                boolean nullable = erosionFilter.process(p4, points, 
                    tmpPointsAdded, tmpPointsRemoved, w, h);
               
            }
        }
        
        int nCorrections = tmpPointsRemoved.size() + tmpPointsAdded.size();

        for (PairInt p2 : tmpPointsRemoved) {
            points.remove(p2);
        }
        for (PairInt p2 : tmpPointsAdded) {
            points.add(p2);
        }
        
        Logger.getLogger(PostLineThinnerCorrections.class.getName())
            .fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
        
        return nCorrections;
    }

    private void correctForCorrectionCreatedSquares(Set<PairInt> points, 
        int imageWidth, int imageHeight) {
        
        int w = imageWidth;
        int h = imageHeight;
        
        /*
        Find a filled in 2x2 square.
       
                   #  #      1
                   #* #      0
                               
         -3 -2 -1  0  1  2
        
         -- expand the region to the surrounding pixels and for each,
           null them if possible
        */
                
        Set<PairInt> ones = new HashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        ones.add(new PairInt(0, -1));
        ones.add(new PairInt(1, 0));
        ones.add(new PairInt(1, -1));
        
        ErosionFilter erosionFilter = new ErosionFilter();
        erosionFilter.overrideEndOfLineCheck();
                        
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
            visit each of the surrounding pixels and the center 4 in order
            of smallest number of neighbors, and null the ones which are 
            nullable without disconnecting lines.
            
                    +  +  +  +   2
                    +  #  #  +   1
                    +  #* #  +   0
                    +  +  +  +  -1
        
             -3 -2 -1  0  1  2
            */
            
            int nn = 0;
            int[] neighborsCount = new int[16];
            PairInt[] neighbors = new PairInt[16];
            
            for (int ii = -1; ii <= 2; ++ii) {
                for (int jj = -2; jj <= 1; ++jj) {
                    
                    int x1 = col + ii;
                    int y1 = row + jj;
                    
                    if ((x1 < 0) || (y1 < 0) || (x1 > (w - 1)) || (y1 > (h - 1))) {
                        continue;
                    }
                    PairInt p3 = new PairInt(x1, y1);
                    if (tmpPointsRemoved.contains(p3) || 
                        (!tmpPointsAdded.contains(p3) && !points.contains(p3))) {
                        continue;
                    }

                    neighbors[nn] = p3;
                    
                    int count = 0;
                    // count the neighbors of the neighbor
                    for (int n2Idx = 0; n2Idx < eightNeighborsX.length; ++n2Idx) {
                        int x2 = x1 + eightNeighborsX[n2Idx];
                        int y2 = y1 + eightNeighborsY[n2Idx];                
                        if ((x2 < 0) || (y2 < 0) || (x2 > (w - 1)) || (y2 > (h - 1))) {
                            continue;
                        }
                        PairInt p4 = new PairInt(x2, y2);
                        if (!tmpPointsRemoved.contains(p4) && 
                            (tmpPointsAdded.contains(p4) || points.contains(p4))) {
                            count++;
                        }
                    }
                    
                    neighborsCount[nn] = count;
                    nn++;
                }
            }
            
            neighbors = Arrays.copyOf(neighbors, nn);
            neighborsCount = Arrays.copyOf(neighborsCount, nn);
            
            CountingSort.sort(neighborsCount, neighbors, 8);

            for (int ii = 0; ii < neighbors.length; ++ii) {
                if (neighborsCount[ii] == 0) {
                    continue;
                }
                PairInt p3 = neighbors[ii]; 
                // adds to tmpPointsRemoved

                boolean nullable = erosionFilter.process(p3, points, 
                    tmpPointsAdded, tmpPointsRemoved, w, h);
            }
        }

        int nCorrections = tmpPointsRemoved.size() + tmpPointsAdded.size();
        
        for (PairInt p2 : tmpPointsRemoved) {
            points.remove(p2);
        }
        for (PairInt p2 : tmpPointsAdded) {
            points.add(p2);
        }
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }

    protected void straightenLines(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
        
        if (edgeGuideImage == null) {
            return;
        }
       
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        curveHelper.straightenLines(points, edgeGuideImage);

    }
    
    private void correctForSingleHole_01(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
        
        /*
          0  0  0  0  0         3                           3
          0  0  #  #  0         2           .  .            2
          0  #  0  #  #         1        .  0  .  #         1
          0  #  #  #* 0         0        .  #  #*           0
             #  0  0  0        -1        #                 -1
         -3 -2 -1  0  1  2  3        -3 -2 -1  0  1  2  3  
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-3, 0)); zeroes.add(new PairInt(-3, -1));
        zeroes.add(new PairInt(-3, -2)); zeroes.add(new PairInt(-3, -3));
        zeroes.add(new PairInt(-2, -2)); zeroes.add(new PairInt(-2, -3));
        zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, -1)); zeroes.add(new PairInt(-1, -3));
        zeroes.add(new PairInt(0, 1)); zeroes.add(new PairInt(0, -3));
        zeroes.add(new PairInt(1, 0)); zeroes.add(new PairInt(1, 1)); zeroes.add(new PairInt(1, -2)); zeroes.add(new PairInt(1, -3));
        
        ones.add(new PairInt(-2, 1)); ones.add(new PairInt(-2, 0)); ones.add(new PairInt(-2, -1));
        ones.add(new PairInt(-1, 0)); ones.add(new PairInt(-1, -2));
        ones.add(new PairInt(0, -1)); ones.add(new PairInt(0, -2));
        ones.add(new PairInt(1, -1));
        
        changeToZeroes.add(new PairInt(-2, -1));
        changeToZeroes.add(new PairInt(-2, 0));
        changeToZeroes.add(new PairInt(-1, -2));
        changeToZeroes.add(new PairInt(0, -2));
        changeToZeroes.add(new PairInt(0, -1));
        
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    private void correctForSingleHole_02(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
        
        /*    
          0  0  0  0  0         3                           3
          0  0  #  0  0         2           .               2
          0  #  0  #  0         1        .  0  .            1
          0  #  #  #* #         0        .  #  #* #         0
             #  0  0  0        -1        #                 -1
         -3 -2 -1  0  1  2           -3 -2 -1  0  1  2        
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-3, 0)); zeroes.add(new PairInt(-3, -1)); 
        zeroes.add(new PairInt(-3, -2)); zeroes.add(new PairInt(-3, -3));
        zeroes.add(new PairInt(-2, -2)); zeroes.add(new PairInt(-2, -3));
        zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, -1)); zeroes.add(new PairInt(-1, -3));
        zeroes.add(new PairInt(0, 1)); zeroes.add(new PairInt(0, -2)); zeroes.add(new PairInt(0, -3));
        zeroes.add(new PairInt(1, 1)); zeroes.add(new PairInt(1, -1)); 
        zeroes.add(new PairInt(1, -2)); zeroes.add(new PairInt(1, -3));
        
        ones.add(new PairInt(-2, 1)); ones.add(new PairInt(-2, 0)); ones.add(new PairInt(-2, -1));
        ones.add(new PairInt(-1, 0)); ones.add(new PairInt(-1, -2));
        ones.add(new PairInt(0, -1));
        ones.add(new PairInt(1, 0));
        
        changeToZeroes.add(new PairInt(-2, -1)); changeToZeroes.add(new PairInt(-2, 0));
        changeToZeroes.add(new PairInt(-1, -2));
        changeToZeroes.add(new PairInt(0, -1));

        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    private void correctForSingleHole_03(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
        
        /*    
          0  0  0  0  0         3                           3
          0  0  #  #  0         2           .  .            2
          0  #  0  #  #         1        .  0  .  #         1
          0  #  #  #* 0         0        .  .  #*           0
          0  0  #  0  0        -1           #              -1
         -3 -2 -1  0  1  2  3        -3 -2 -1  0  1  2  3  
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-3, 1)); zeroes.add(new PairInt(-3, 0)); zeroes.add(new PairInt(-3, -1));
        zeroes.add(new PairInt(-3, -2)); zeroes.add(new PairInt(-3, -3));
        zeroes.add(new PairInt(-2, 1)); zeroes.add(new PairInt(-2, -2)); zeroes.add(new PairInt(-2, -3));
        zeroes.add(new PairInt(-1, -1)); zeroes.add(new PairInt(-1, -3));
        zeroes.add(new PairInt(0, 1)); zeroes.add(new PairInt(0, -3));
        zeroes.add(new PairInt(1, 1)); zeroes.add(new PairInt(1, 0)); 
        zeroes.add(new PairInt(1, -2)); zeroes.add(new PairInt(1, -3));
        
        ones.add(new PairInt(-2, 0)); ones.add(new PairInt(-2, -1));
        ones.add(new PairInt(-1, 1)); ones.add(new PairInt(-1, 0)); ones.add(new PairInt(-1, -2));
        ones.add(new PairInt(0, -1)); ones.add(new PairInt(0, -2));
        ones.add(new PairInt(1, -1));
        
        changeToZeroes.add(new PairInt(-2, 0)); changeToZeroes.add(new PairInt(-2, -1));
        changeToZeroes.add(new PairInt(-1, 0)); changeToZeroes.add(new PairInt(-1, -2));
        changeToZeroes.add(new PairInt(0, -1)); changeToZeroes.add(new PairInt(0, -2));
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    private void correctForSingleHole_04(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
        
        /*    
          0  0  0  0  0         3                           3
          0  0  #  #  0         2           .  .            2
          0  #  0  #  0         1        .  0  .            1
          0  #  #  #* #         0        .  .  #*           0
          0  0  #  0  0        -1           #              -1
         -3 -2 -1  0  1  2  3        -3 -2 -1  0  1  2  3  
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-3, 1)); zeroes.add(new PairInt(-3, 0)); zeroes.add(new PairInt(-3, -1));
        zeroes.add(new PairInt(-3, -2)); zeroes.add(new PairInt(-3, -3));
        zeroes.add(new PairInt(-2, 1)); zeroes.add(new PairInt(-2, -2)); zeroes.add(new PairInt(-2, -3));
        zeroes.add(new PairInt(-1, -1)); zeroes.add(new PairInt(-1, -3));
        zeroes.add(new PairInt(0, 1)); zeroes.add(new PairInt(0, -3));
        zeroes.add(new PairInt(1, 1)); zeroes.add(new PairInt(1, -1)); 
        zeroes.add(new PairInt(1, -2)); zeroes.add(new PairInt(1, -3));
        
        ones.add(new PairInt(-2, 0)); ones.add(new PairInt(-2, -1));
        ones.add(new PairInt(-1, 1)); ones.add(new PairInt(-1, 0)); ones.add(new PairInt(-1, -2));
        ones.add(new PairInt(0, -1)); ones.add(new PairInt(0, -2));
        ones.add(new PairInt(1, 0));
        
        changeToZeroes.add(new PairInt(-2, 0)); changeToZeroes.add(new PairInt(-2, -1));
        changeToZeroes.add(new PairInt(-1, 0)); changeToZeroes.add(new PairInt(-1, -2));
        changeToZeroes.add(new PairInt(0, -1)); changeToZeroes.add(new PairInt(0, -2));
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    private void correctForSingleHole_05(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
        
        /*    
          0  0  0  0  0         3                           3
          0  0  #  #  0         2           .  .            2
          0  #  0  #  #         1        .  0  .  #         1
          #  #  #  #* 0         0        #  #  #*           0
          0  0  0  0  0        -1                          -1
         -3 -2 -1  0  1  2  3        -3 -2 -1  0  1  2  3  
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-3, 1)); zeroes.add(new PairInt(-3, -1));
        zeroes.add(new PairInt(-3, -2)); zeroes.add(new PairInt(-3, -3));
        zeroes.add(new PairInt(-2, 1)); zeroes.add(new PairInt(-2, -2)); zeroes.add(new PairInt(-2, -3));
        zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, -1)); zeroes.add(new PairInt(-1, -3));
        zeroes.add(new PairInt(0, 1)); zeroes.add(new PairInt(0, -3));
        zeroes.add(new PairInt(1, 1)); zeroes.add(new PairInt(1, 0)); 
        zeroes.add(new PairInt(1, -2)); zeroes.add(new PairInt(1, -3));
        
        ones.add(new PairInt(-3, 0));
        ones.add(new PairInt(-2, 0)); ones.add(new PairInt(-2, -1));
        ones.add(new PairInt(-1, 0)); ones.add(new PairInt(-1, -2));
        ones.add(new PairInt(0, -1)); ones.add(new PairInt(0, -2));
        ones.add(new PairInt(1, -1));
        
        changeToZeroes.add(new PairInt(-2, -1));
        changeToZeroes.add(new PairInt(-1, -2));
        changeToZeroes.add(new PairInt(0, -1)); changeToZeroes.add(new PairInt(0, -2));
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    private void correctForSingleHole_06(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
        
        /*    
          0  0  0  0  #         3                 #         3
          0  #  #  #            2        .  .  #            2
          0  #  0  #  0         1        .  @  .            1
          0  #  #  #* 0         0        #  .  .*           0
          0  #  0  0  0        -1        #                 -1
         -3 -2 -1  0  1  2  3        -3 -2 -1  0  1  2  3  
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-3, 1)); zeroes.add(new PairInt(-3, 0)); zeroes.add(new PairInt(-3, -1));
        zeroes.add(new PairInt(-3, -2)); zeroes.add(new PairInt(-3, -3));
        zeroes.add(new PairInt(-2, -3));
        zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, -1)); zeroes.add(new PairInt(-1, -3));
        zeroes.add(new PairInt(0, 1)); zeroes.add(new PairInt(0, -3));
        zeroes.add(new PairInt(1, 1)); zeroes.add(new PairInt(1, 0)); 
        zeroes.add(new PairInt(1, -1));
        
        ones.add(new PairInt(-2, 1)); ones.add(new PairInt(-2, 0)); 
        ones.add(new PairInt(-2, -1)); ones.add(new PairInt(-2, -2));
        ones.add(new PairInt(-1, 0)); ones.add(new PairInt(-1, -2));
        ones.add(new PairInt(0, -1)); ones.add(new PairInt(0, -2));
        ones.add(new PairInt(1, -3));
        
        changeToZeroes.add(new PairInt(-2, -1)); changeToZeroes.add(new PairInt(-2, -2));
        changeToZeroes.add(new PairInt(-1, 0)); changeToZeroes.add(new PairInt(-1, -2));
        changeToZeroes.add(new PairInt(0, 0)); changeToZeroes.add(new PairInt(0, -1));
                    
        changeToOnes.add(new PairInt(-1, -1));
        
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    private void correctForSingleHole_07(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
        
        /*  
                 0  0  0  0  #      4                    #      4
              0  0  #  #  #         3           .  .  #         3
              0  #  #  0  #  0      2        .  .  @  .         2
              0  #  0  #  #  0      1        .  @  .  .         1
              0  #  #  #* 0  0      0        .  #  .*           0
              0  #  0  0  0        -1        #                 -1
             -3 -2 -1  0  1  2  3        -3 -2 -1  0  1  2  3  
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-3, 1)); zeroes.add(new PairInt(-3, 0)); zeroes.add(new PairInt(-3, -1));
        zeroes.add(new PairInt(-3, -2)); zeroes.add(new PairInt(-3, -3));
        zeroes.add(new PairInt(-2, -3)); zeroes.add(new PairInt(-2, -4));
        zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, -1)); zeroes.add(new PairInt(-1, -4));
        zeroes.add(new PairInt(0, 1)); zeroes.add(new PairInt(0, -2)); zeroes.add(new PairInt(0, -4));
        zeroes.add(new PairInt(1, 1)); zeroes.add(new PairInt(1, 0)); 
        zeroes.add(new PairInt(1, -4));
        zeroes.add(new PairInt(2, 0)); zeroes.add(new PairInt(2, -1)); zeroes.add(new PairInt(2, -2));
        
        ones.add(new PairInt(-2, 1)); ones.add(new PairInt(-2, 0)); 
        ones.add(new PairInt(-2, -1)); ones.add(new PairInt(-2, -2));
        ones.add(new PairInt(-1, 0)); ones.add(new PairInt(-1, -2)); ones.add(new PairInt(-1, -3));
        ones.add(new PairInt(0, -1)); ones.add(new PairInt(0, -3));
        ones.add(new PairInt(1, -1)); ones.add(new PairInt(1, -2)); ones.add(new PairInt(1, -3));
        ones.add(new PairInt(2, -4));
        
        changeToZeroes.add(new PairInt(-2, 0)); changeToZeroes.add(new PairInt(-2, -1)); changeToZeroes.add(new PairInt(-2, -2));
        changeToZeroes.add(new PairInt(-1, -2)); changeToZeroes.add(new PairInt(-1, -3));
        changeToZeroes.add(new PairInt(0, 0)); changeToZeroes.add(new PairInt(0, -1)); changeToZeroes.add(new PairInt(0, -3));
        changeToZeroes.add(new PairInt(1, -1)); changeToZeroes.add(new PairInt(1, -2));
        
        changeToOnes.add(new PairInt(-1, -1));
        changeToOnes.add(new PairInt(0, -2));
        
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    private void correctForSingleHole_08(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
        
        /*  
                 0  0  0  0  #      4                    #      4
              0  0  #  #  #         3           .  .  #         3
           0  0  #  #  0  #  0      2        .  .  @  .         2
           0  #  #  0  #  #  0      1     .  .  @  .  .         1
           0  #  0  #  #* 0  0      0     .  @  .  .*           0
           0  #  #  #  0  0        -1     #  .  .              -1
           0  #  0  0              -2     #                    -2
          -4 -3 -2 -1  0  1  2  3        -3 -2 -1  0  1  2  3  
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-4, 2)); zeroes.add(new PairInt(-4, 1)); zeroes.add(new PairInt(-4, 0));
        zeroes.add(new PairInt(-4, -1)); zeroes.add(new PairInt(-4, -2));
        zeroes.add(new PairInt(-3, -2)); zeroes.add(new PairInt(-3, -3));
        zeroes.add(new PairInt(-2, 2)); zeroes.add(new PairInt(-2, 0));
        zeroes.add(new PairInt(-2, -3)); zeroes.add(new PairInt(-2, -4));
        zeroes.add(new PairInt(-1, 2)); zeroes.add(new PairInt(-1, -1)); zeroes.add(new PairInt(-1, -4));
        zeroes.add(new PairInt(0, 1)); zeroes.add(new PairInt(0, -2)); zeroes.add(new PairInt(0, -4));
        zeroes.add(new PairInt(1, 1)); zeroes.add(new PairInt(1, 0)); zeroes.add(new PairInt(1, -4));
        zeroes.add(new PairInt(2, 0)); zeroes.add(new PairInt(2, -1)); zeroes.add(new PairInt(2, -2));
        
        ones.add(new PairInt(-3, 2)); ones.add(new PairInt(-3, 1)); 
        ones.add(new PairInt(-3, 0)); ones.add(new PairInt(-3, -1));
        ones.add(new PairInt(-2, 1)); ones.add(new PairInt(-2, -1)); ones.add(new PairInt(-2, -2));
        ones.add(new PairInt(-1, 1)); ones.add(new PairInt(-1, 0)); 
        ones.add(new PairInt(-1, -2)); ones.add(new PairInt(-1, -3));
        ones.add(new PairInt(0, -1)); ones.add(new PairInt(0, -3));
        ones.add(new PairInt(1, -1)); ones.add(new PairInt(1, -2)); ones.add(new PairInt(1, -3));
        ones.add(new PairInt(2, -4));
        
        changeToZeroes.add(new PairInt(-3, 0)); changeToZeroes.add(new PairInt(-3, -1));
        changeToZeroes.add(new PairInt(-2, 1)); changeToZeroes.add(new PairInt(-2, -1)); changeToZeroes.add(new PairInt(-2, -2));
        changeToZeroes.add(new PairInt(-1, 1)); changeToZeroes.add(new PairInt(-1, 0)); 
        changeToZeroes.add(new PairInt(-1, -2)); changeToZeroes.add(new PairInt(-1, -3));
        changeToZeroes.add(new PairInt(0, 0)); changeToZeroes.add(new PairInt(0, -1)); changeToZeroes.add(new PairInt(0, -3));
        changeToZeroes.add(new PairInt(1, -1)); changeToZeroes.add(new PairInt(1, -2));
        
        changeToOnes.add(new PairInt(-2, 0));
        changeToOnes.add(new PairInt(-1, -1));
        changeToOnes.add(new PairInt(0, -2));
        
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    private void correctForSingleHole_09(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
        
        /*  
                                    4                           4
                 0  0  0  #  0      3                 #         3
                 0  #  #  #  0      2           .  #  #         2
                 0  #  0  #  0      1           #     #         1 
                 #  #  #* #  #      0        #  .  .* .  #      0
                    0  0  0        -1                          -1
                                   -2                          -2
          -4 -3 -2 -1  0  1  2  3        -3 -2 -1  0  1  2  3  
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-2, -1)); zeroes.add(new PairInt(-2, -2)); zeroes.add(new PairInt(-2, -3));
        zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, -4));
        zeroes.add(new PairInt(0, 1)); zeroes.add(new PairInt(0, -1)); zeroes.add(new PairInt(0, -3));
        zeroes.add(new PairInt(1, 1));
        zeroes.add(new PairInt(2, -1)); zeroes.add(new PairInt(2, -2)); zeroes.add(new PairInt(2, -3));
        
        ones.add(new PairInt(-2, 0));
        ones.add(new PairInt(-1, 0)); ones.add(new PairInt(-1, -1)); ones.add(new PairInt(-1, -2)); 
        ones.add(new PairInt(0, -2));
        ones.add(new PairInt(1, 0)); ones.add(new PairInt(1, -1));
        ones.add(new PairInt(1, -2)); ones.add(new PairInt(1, -3));
        ones.add(new PairInt(2, 0));
        
        changeToZeroes.add(new PairInt(-1, 0)); changeToZeroes.add(new PairInt(-1, -2)); 
        changeToZeroes.add(new PairInt(0, 0));
        changeToZeroes.add(new PairInt(1, 0));
        
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    private void correctForSingleHole_10(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
        
        /*
                     0 0        3                          3
                   0 # 0 0      2           # 0 .          2
                 # # 0 # 0 0    1           # # 0 .        1
                   # #*0 # 0    0             # #*0 .      0
                     # # #     -1               #         -1
        
        -6-5-4-3-2-1 0 1 2 3       -6-5-4-3-2-1 0 1 2 3
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-1, -2));
        zeroes.add(new PairInt(0, -1)); zeroes.add(new PairInt(0, -3));
        zeroes.add(new PairInt(1, 0)); zeroes.add(new PairInt(1, -2)); zeroes.add(new PairInt(1, -3));
        zeroes.add(new PairInt(2, -1)); zeroes.add(new PairInt(2, -2));
        zeroes.add(new PairInt(3, 0)); zeroes.add(new PairInt(3, -1));
        
        ones.add(new PairInt(-2, -1));
        ones.add(new PairInt(-1, 0)); ones.add(new PairInt(-1, -1)); 
        ones.add(new PairInt(0, 1)); ones.add(new PairInt(0, -2));
        ones.add(new PairInt(1, 1)); ones.add(new PairInt(1, -1));
        ones.add(new PairInt(2, 1)); ones.add(new PairInt(2, 0));
        
        changeToZeroes.add(new PairInt(0, -2));
        changeToZeroes.add(new PairInt(1, -1));
        changeToZeroes.add(new PairInt(2, 0));
        
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    private void correctForSingleHole_11(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
        
        /*
                 0  0  0  #         4                 #         4
              0  0  #  #  #  0      3           .  .  #         3
              0  #  #  0  #  0      2        .  .  @  .         2
              0  #  0  #  #  0      1        .  @  .  .         1
              0  #  #  #* 0  0      0        #  .  .*           0
                 #  0  0  0        -1        #                 -1
             -3 -2 -1  0  1  2  3        -3 -2 -1  0  1  2  3
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-3, 0)); zeroes.add(new PairInt(-3, -1));
        zeroes.add(new PairInt(-3, -2)); zeroes.add(new PairInt(-3, -3));
        zeroes.add(new PairInt(-2, -3)); zeroes.add(new PairInt(-2, -4));
        zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, -1)); zeroes.add(new PairInt(-1, -4));
        zeroes.add(new PairInt(0, 1)); zeroes.add(new PairInt(0, -2)); zeroes.add(new PairInt(0, -4));
        zeroes.add(new PairInt(1, 1)); zeroes.add(new PairInt(1, 0)); 
        zeroes.add(new PairInt(2, 0)); zeroes.add(new PairInt(2, -1));
        zeroes.add(new PairInt(2, -2)); zeroes.add(new PairInt(2, -3));
        
        ones.add(new PairInt(-2, 1)); ones.add(new PairInt(-2, 0)); ones.add(new PairInt(-2, -1));
        ones.add(new PairInt(-2, -2)); 
        ones.add(new PairInt(-1, 0)); ones.add(new PairInt(-1, -2)); ones.add(new PairInt(-1, -3)); 
        ones.add(new PairInt(0, -1)); ones.add(new PairInt(0, -3));
        ones.add(new PairInt(1, -1)); ones.add(new PairInt(1, -2));
        ones.add(new PairInt(1, -3)); ones.add(new PairInt(1, -4));
        
        changeToZeroes.add(new PairInt(-2, -1)); changeToZeroes.add(new PairInt(-2, -2));
        changeToZeroes.add(new PairInt(-1, 0)); changeToZeroes.add(new PairInt(-1, -2));
        changeToZeroes.add(new PairInt(-1, -3));
        changeToZeroes.add(new PairInt(0, 0)); changeToZeroes.add(new PairInt(0, -1));
        changeToZeroes.add(new PairInt(0, -3));
        changeToZeroes.add(new PairInt(1, -1)); changeToZeroes.add(new PairInt(1, -2));
        
        changeToOnes.add(new PairInt(-1, -1)); changeToOnes.add(new PairInt(0, -2));
        
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    private void correctForSingleHole_12(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
        
        /*
                 0  0  0  0  #      4                   #       4
              0  0  #  #  #  #      3           .  .  # .       3
              0  #  #  0  #  0      2        .  .  @  .         2
              0  #  0  #  #  0      1        .  @  .  .         1
              0  #  #  #* 0  0      0        #  .  .*           0
                 #  0  0  0        -1        #                 -1
             -3 -2 -1  0  1  2  3        -3 -2 -1  0  1  2  3
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
  
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-3, 0)); zeroes.add(new PairInt(-3, -1));
        zeroes.add(new PairInt(-3, -2)); zeroes.add(new PairInt(-3, -3));
        zeroes.add(new PairInt(-2, -3)); zeroes.add(new PairInt(-2, -4));
        zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, -1)); zeroes.add(new PairInt(-1, -4));
        zeroes.add(new PairInt(0, 1)); zeroes.add(new PairInt(0, -2)); zeroes.add(new PairInt(0, -4));
        zeroes.add(new PairInt(1, 1)); zeroes.add(new PairInt(1, 0)); zeroes.add(new PairInt(1, -4)); 
        zeroes.add(new PairInt(2, 0)); zeroes.add(new PairInt(2, -1));
        zeroes.add(new PairInt(2, -2));
        
        ones.add(new PairInt(-2, 1)); ones.add(new PairInt(-2, 0)); ones.add(new PairInt(-2, -1));
        ones.add(new PairInt(-2, -2));
        ones.add(new PairInt(-1, 0)); ones.add(new PairInt(-1, -2)); ones.add(new PairInt(-1, -3)); 
        ones.add(new PairInt(0, -1)); ones.add(new PairInt(0, -3));
        ones.add(new PairInt(1, -1)); ones.add(new PairInt(1, -2)); 
        ones.add(new PairInt(1, -3)); 
        ones.add(new PairInt(2, -3)); ones.add(new PairInt(2, -4));
        
        /*
                 0  0  0  0  #      4                    #      4
              0  0  #  #  #  #      3           .  .  #  .      3
              0  #  #  0  #  0      2        .  .  @  .         2
              0  #  0  #  #  0      1        .  @  .  .         1
              0  #  #  #* 0  0      0        #  .  .*           0
                 #  0  0  0        -1        #                 -1
             -3 -2 -1  0  1  2  3        -3 -2 -1  0  1  2  3
        */
        changeToZeroes.add(new PairInt(-2, -1)); changeToZeroes.add(new PairInt(-2, -2));
        changeToZeroes.add(new PairInt(-1, 0)); changeToZeroes.add(new PairInt(-1, -2));
        changeToZeroes.add(new PairInt(-1, -3));
        changeToZeroes.add(new PairInt(0, 0)); changeToZeroes.add(new PairInt(0, -1));
        changeToZeroes.add(new PairInt(0, -3));
        changeToZeroes.add(new PairInt(1, -1)); changeToZeroes.add(new PairInt(1, -2));
        changeToZeroes.add(new PairInt(2, -3));
        
        changeToOnes.add(new PairInt(-1, -1)); changeToOnes.add(new PairInt(0, -2));
        
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    /**
     * given an ordered list of corners, remove any corners that appear to be 
     * the artifact of a single pixel staircases in a straight line that is 
     * nearly vertical or nearly horizontal. Note that ordered list of corners
     * means that edgeCorners data are internally ordered by increasing values
     * of idx.
     * runtime complexity is O(N_corners).
     * @param edgeCorners 
     * @param curveIsClosed 
     */
    public static void removeSingleStairsAliasArtifacts(CornerArray edgeCorners, 
        boolean curveIsClosed) {
        
        int nMaxIter = edgeCorners.getN();
        int nIter = 0;
        int nRemoved = 0;
        while ((nIter == 0) || ((nIter < nMaxIter) && (nRemoved > 0))) {
            
            nRemoved = removeSingleStairAliasArtifacts(edgeCorners, 
                curveIsClosed);
            
            nIter++;
        }
    }
    
    /**
     * given an ordered list of corners, remove any corners that appear to be 
     * the artifact of a single pixel staircase in a straight line that is 
     * nearly vertical or nearly horizontal. Note that ordered list of corners
     * means that edgeCorners data are internally ordered by increasing values
     * of idx.
     * runtime complexity is O(N_corners).
     * @param edgeCorners 
     * @param curveIsClosed 
     */
    protected static int removeSingleStairAliasArtifacts(CornerArray edgeCorners, 
        boolean curveIsClosed) {
        
        if (edgeCorners.getN() < 3) {
            return 0;
        }
        
        /*
        look for pattern __-- embedded between 2 corners that have the same 
        y's as the 1 pixel step
        then do same for vertical pattern.
        */
        
        List<Integer> remove = new ArrayList<Integer>();
        
        for (int i = 1; i < edgeCorners.getN() - 1; ++i) {
            int xa = Math.round(edgeCorners.getX(i - 1));
            int ya = Math.round(edgeCorners.getY(i - 1));
            
            int xb = Math.round(edgeCorners.getX(i));
            int yb = Math.round(edgeCorners.getY(i));
            
            int xc = Math.round(edgeCorners.getX(i + 1));
            int yc = Math.round(edgeCorners.getY(i + 1));
            
            int diffYac = Math.abs(ya - yc);
            int diffYab = Math.abs(ya - yb);
            int diffYbc = Math.abs(yb - yc);
                        
            if ((diffYac == 1) && (
                ((diffYab == 0) && (diffYbc == 1)) ||
                ((diffYab == 1) && (diffYbc == 0)))) {
                
                if (((xa < xb) && (xb < xc)) || ((xa > xb) && (xb > xc))) {
                    remove.add(Integer.valueOf(i));
                    ++i;
                    continue;
                }
            }
        
            int diffXac = Math.abs(xa - xc);
            int diffXab = Math.abs(xa - xb);
            int diffXbc = Math.abs(xb - xc);
            
            if ((diffXac == 1) && (
                ((diffXab == 0) && (diffXbc == 1)) ||
                ((diffXab == 1) && (diffXbc == 0)))) {
                
                if (((ya < yb) && (yb < yc)) || ((ya > yb) && (yb > yc))) {
                    remove.add(Integer.valueOf(i));
                    ++i;
                }
            }
        }
        
        if (curveIsClosed) {
            
            int n = edgeCorners.getN();
            
            boolean contains = remove.contains(Integer.valueOf(n - 1)) || 
                remove.contains(Integer.valueOf(0)) ||
                remove.contains(Integer.valueOf(1));
            
            if (!contains) {
                
                int xa = Math.round(edgeCorners.getX(n - 1));
                int ya = Math.round(edgeCorners.getY(n - 1));

                int xb = Math.round(edgeCorners.getX(0));
                int yb = Math.round(edgeCorners.getY(0));

                int xc = Math.round(edgeCorners.getX(1));
                int yc = Math.round(edgeCorners.getY(1));
                
                int diffYac = Math.abs(ya - yc);
                int diffYab = Math.abs(ya - yb);
                int diffYbc = Math.abs(yb - yc);

                if ((diffYac == 1) && (((diffYab == 0) && (diffYbc == 1))
                    || ((diffYab == 1) && (diffYbc == 0)))) {
                    
                    if (((xa < xb) && (xb < xc)) || ((xa > xb) && (xb > xc))) {
                        remove.add(Integer.valueOf(0));
                    }
                    
                } else {

                    int diffXac = Math.abs(xa - xc);
                    int diffXab = Math.abs(xa - xb);
                    int diffXbc = Math.abs(xb - xc);

                    if ((diffXac == 1) && (((diffXab == 0) && (diffXbc == 1))
                        || ((diffXab == 1) && (diffXbc == 0)))) {
                        
                        if (((ya < yb) && (yb < yc)) || ((ya > yb) && (yb > yc))) {
                            remove.add(Integer.valueOf(0));
                        }
                    }
                }
            }
        }
        
        for (int i = (remove.size() - 1); i > -1; --i) {
            int idx = remove.get(i).intValue();
            edgeCorners.removeRange(idx, idx);
        }
        
        return remove.size();
    }
    
    /**
     * corrects for the curve artifacts introduced after phase congruency edge 
     * thinning followed by phase angle step concatenation.
     * <pre>
     * this pattern embedded in many different ways in a curve can be corrected
     * by removing two oppossing points around the empty center and filling
     * the center as long as connecting points are not disconnected
     * (avoiding breaking a curve into two pieces).
     *       #
     *      # #
     *       #
     * </pre>
     * @param points
     * @param imageWidth
     * @param imageHeight 
     */
    public void correctForHolePattern100(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
                
        // first try to change the 2 horizontal points to zeros and center to one
        // then repeat check with the 2 vertical points if pattern still exists
        
        PairInt[][] neighborCoordOffsets
            = AbstractLineThinner.createCoordinatePointsForEightNeighbors(0, 0);
        
        int nCorrections = 0;
        int nIter = 0;
        
        Set<PairInt> editedPoints = new HashSet<PairInt>(points);
        
        while ((nIter == 0) || (nCorrections > 0)) {
            
            nCorrections = 0;
            
            for (PairInt p : points) {
                int x = p.getX();
                int y = p.getY();

                //TODO: should simplify this to use transformations by 90.

                // test for pattern where the current (x,y) is at bottom of diamond
                if (editedPoints.contains(new PairInt(x - 1, y + 1)) && 
                    editedPoints.contains(new PairInt(x + 1, y + 1)) && 
                    editedPoints.contains(new PairInt(x, y + 2)) && 
                    !editedPoints.contains(new PairInt(x, y + 1)) 
                    ){

                    /*
                        0  #  0     2    
                        #  0  #     1   
                        0  *  0     0  x,y is at (0, 0) here
                       -1  0  1     
                    */
                    // see if setting x,y and nulling x-1 and x+1 disconnects
                    // any points

                    PairInt pCenter = new PairInt(x, y + 1);

                    editedPoints.add(pCenter);

                    PairInt pLeft = new PairInt(x - 1, y + 1);

                    // test for left disconnects
                    boolean disconnects = ImageSegmentation.doesDisconnect(editedPoints, 
                        neighborCoordOffsets, x - 1, y + 1, imageWidth, imageHeight);

                    if (!disconnects) {

                        editedPoints.remove(pLeft);

                        // test for right point
                        disconnects = ImageSegmentation.doesDisconnect(editedPoints, 
                            neighborCoordOffsets, x + 1, y + 1, imageWidth, imageHeight);

                        if (!disconnects) {
                            // tests passed, so finish changes and continue
                            editedPoints.remove(new PairInt(x + 1, y + 1));
                            ++nCorrections;
                            continue;                        
                        } else {
                            // undo remove pLeft to reset for next tests
                            editedPoints.add(pLeft);
                        }
                    }

                    // test above
                    PairInt pAbove = new PairInt(x, y + 2);

                    disconnects = ImageSegmentation.doesDisconnect(editedPoints, 
                        neighborCoordOffsets, x, y + 2, imageWidth, imageHeight);

                    if (disconnects) {
                        // tests failed, so undo add pCenter and continue
                        editedPoints.remove(pCenter);
                        continue;
                    }

                    // test for pBelow
                    editedPoints.remove(pAbove);

                    PairInt pBelow = p;

                    disconnects = ImageSegmentation.doesDisconnect(editedPoints, 
                        neighborCoordOffsets, x, y, imageWidth, imageHeight);

                    if (!disconnects) {
                        // pattern succeeded, so remove last point and continue
                        editedPoints.remove(pBelow);
                        ++nCorrections;
                        continue;
                    }

                    // else tests didn't pass so undo temporary changes
                    editedPoints.remove(pCenter);
                    editedPoints.add(pAbove);
                }

                // test for pattern where the current (x,y) is at top of diamond
                if (editedPoints.contains(new PairInt(x - 1, y - 1)) && 
                    editedPoints.contains(new PairInt(x + 1, y - 1)) && 
                    editedPoints.contains(new PairInt(x, y - 2)) && 
                    !editedPoints.contains(new PairInt(x, y - 1)) 
                    ){

                    /*
                        0  *  0     0    x,y is at (0, 0) here
                        #  0  #    -1   
                        0  #  0    -2  
                       -1  0  1     
                    */

                    // see if setting x,y and nulling x-1 and x+1 disconnects
                    // any points

                    PairInt pCenter = new PairInt(x, y - 1);

                    editedPoints.add(pCenter);

                    PairInt pLeft = new PairInt(x - 1, y - 1);

                    // test for left disconnects
                    boolean disconnects = ImageSegmentation.doesDisconnect(editedPoints, 
                        neighborCoordOffsets, x - 1, y - 1, imageWidth, imageHeight);

                    if (!disconnects) {

                        editedPoints.remove(pLeft);

                        // test for right point
                        disconnects = ImageSegmentation.doesDisconnect(editedPoints, 
                            neighborCoordOffsets, x + 1, y - 1, imageWidth, imageHeight);

                        if (!disconnects) {
                            // tests passed, so finish changes and continue
                            editedPoints.remove(new PairInt(x + 1, y - 1));
                            ++nCorrections;
                            continue;                        
                        } else {
                            // undo remove pLeft to reset for next tests
                            editedPoints.add(pLeft);
                        }
                    }

                    // test below
                    PairInt pBelow = new PairInt(x, y - 2);

                    disconnects = ImageSegmentation.doesDisconnect(editedPoints, 
                        neighborCoordOffsets, x, y - 2, imageWidth, imageHeight);

                    if (disconnects) {
                        // tests failed, so undo add pCenter and continue
                        editedPoints.remove(pCenter);
                        continue;
                    }

                    // test for pAbove
                    editedPoints.remove(pBelow);

                    PairInt pAbove = p;

                    disconnects = ImageSegmentation.doesDisconnect(editedPoints, 
                        neighborCoordOffsets, x, y, imageWidth, imageHeight);

                    if (!disconnects) {
                        // pattern succeeded, so remove last point and continue
                        editedPoints.remove(pAbove);
                        ++nCorrections;
                        continue;
                    }

                    // else tests didn't pass so undo temporary changes
                    editedPoints.remove(pCenter);
                    editedPoints.add(pBelow);
                }

                // test for pattern where the current (x,y) is at left of diamond
                if (editedPoints.contains(new PairInt(x + 1, y + 1)) && 
                    editedPoints.contains(new PairInt(x + 1, y - 1)) && 
                    editedPoints.contains(new PairInt(x + 2, y)) && 
                    !editedPoints.contains(new PairInt(x + 1, y)) 
                    ){

                    /*
                              0  #  0     1    
                              *  0  #     0  x,y is at (0, 0) here
                              0  #  0    -1  
                          -1  0  1  2 
                    */
                    // see if setting x,y and nulling x-1 and x+1 disconnects
                    // any points

                    PairInt pCenter = new PairInt(x, y + 1);

                    editedPoints.add(pCenter);

                    // test for above disconnects
                    PairInt pAbove = new PairInt(x + 1, y + 1);

                    boolean disconnects = ImageSegmentation.doesDisconnect(editedPoints, 
                        neighborCoordOffsets, x + 1, y + 1, imageWidth, imageHeight);

                    if (!disconnects) {

                        editedPoints.remove(pAbove);

                        // test for below disconnects
                        disconnects = ImageSegmentation.doesDisconnect(editedPoints, 
                            neighborCoordOffsets, x + 1, y - 1, imageWidth, imageHeight);

                        if (!disconnects) {
                            // tests passed, so finish changes and continue
                            editedPoints.remove(new PairInt(x + 1, y - 1));
                            ++nCorrections;
                            continue;                        
                        } else {
                            // undo remove pAbove to reset for next tests
                            editedPoints.add(pAbove);
                        }
                    }
                    /*
                              0  #  0     1    
                              *  0  #     0  x,y is at (0, 0) here
                              0  #  0    -1  
                          -1  0  1  2 
                    */

                    // test right
                    PairInt pRight = new PairInt(x + 2, y);

                    disconnects = ImageSegmentation.doesDisconnect(editedPoints, 
                        neighborCoordOffsets, x + 2, y, imageWidth, imageHeight);

                    if (disconnects) {
                        // tests failed, so undo add pCenter and continue
                        editedPoints.remove(pCenter);
                        continue;
                    }

                    // test for pLeft
                    editedPoints.remove(pRight);

                    PairInt pLeft = p;

                    disconnects = ImageSegmentation.doesDisconnect(editedPoints, 
                        neighborCoordOffsets, x, y, imageWidth, imageHeight);

                    if (!disconnects) {
                        // pattern succeeded, so remove last point and continue
                        editedPoints.remove(pLeft);
                        ++nCorrections;
                        continue;
                    }

                    // else tests didn't pass so undo temporary changes
                    editedPoints.remove(pCenter);
                    editedPoints.add(pRight);
                }

                // test for pattern where the current (x,y) is at right of diamond
                if (editedPoints.contains(new PairInt(x - 1, y + 1)) && 
                    editedPoints.contains(new PairInt(x - 1, y - 1)) && 
                    editedPoints.contains(new PairInt(x - 2, y)) && 
                    !editedPoints.contains(new PairInt(x - 1, y)) 
                    ){

                    /*
                              0  #  0     1    
                              #  0  *     0  x,y is at (0, 0) here
                              0  #  0    -1  
                             -2 -1  0  1  2 
                    */
                    // see if setting x,y and nulling x-1 and x+1 disconnects
                    // any points

                    PairInt pCenter = new PairInt(x - 1, y);

                    editedPoints.add(pCenter);

                    // test for above disconnects
                    PairInt pAbove = new PairInt(x - 1, y + 1);

                    boolean disconnects = ImageSegmentation.doesDisconnect(editedPoints, 
                        neighborCoordOffsets, x - 1, y + 1, imageWidth, imageHeight);

                    if (!disconnects) {

                        editedPoints.remove(pAbove);

                        // test for below disconnects
                        disconnects = ImageSegmentation.doesDisconnect(editedPoints, 
                            neighborCoordOffsets, x - 1, y - 1, imageWidth, imageHeight);

                        if (!disconnects) {
                            // tests passed, so finish changes and continue
                            editedPoints.remove(new PairInt(x - 1, y - 1));
                            ++nCorrections;
                            continue;                        
                        } else {
                            // undo remove pAbove to reset for next tests
                            editedPoints.add(pAbove);
                        }
                    }
                    /*
                              0  #  0     1    
                              #  0  *     0  x,y is at (0, 0) here
                              0  #  0    -1  
                             -2 -1  0  1  2 
                    */ 

                    // test left
                    PairInt pLeft = new PairInt(x - 2, y);

                    disconnects = ImageSegmentation.doesDisconnect(editedPoints, 
                        neighborCoordOffsets, x - 2, y, imageWidth, imageHeight);

                    if (disconnects) {
                        // tests failed, so undo add pCenter and continue
                        editedPoints.remove(pCenter);
                        continue;
                    }

                    // test for pRight
                    editedPoints.remove(pLeft);

                    PairInt pRight = p;

                    disconnects = ImageSegmentation.doesDisconnect(editedPoints, 
                        neighborCoordOffsets, x, y, imageWidth, imageHeight);

                    if (!disconnects) {
                        // pattern succeeded, so remove last point and continue
                        editedPoints.remove(pRight);
                        ++nCorrections;
                        continue;
                    }

                    // else tests didn't pass so undo temporary changes
                    editedPoints.remove(pCenter);
                    editedPoints.add(pLeft);
                }
            }
            ++nIter;
        }
    
        points.clear();
        points.addAll(editedPoints);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    /**
     * corrects for the curve artifacts introduced after phase congruency edge 
     * thinning followed by phase angle step concatenation.
     * <pre>
     * the pattern is a single pixel displaced on a vertical or horizontal 
     * line
     *        #
     *     ### ###
     *     
     * </pre>
     * @param points
     * @param imageWidth
     * @param imageHeight 
     */
    public void correctForLineHatHoriz(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
                
        /*
                    .  .  .         1                           2
                 .  .  #  .  .      0                           1
                 #  #  .  #  #     -1                           0
                 .  .  .  .  .     -2                          -1
             -3 -2 -1  0  1  2  3        -3 -2 -1  0  1  2  3
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-2, 2)); zeroes.add(new PairInt(-2, 0));
        zeroes.add(new PairInt(-1, 2)); zeroes.add(new PairInt(-1, 0)); zeroes.add(new PairInt(-1, -1));
        zeroes.add(new PairInt(0, 2)); zeroes.add(new PairInt(0, 1)); zeroes.add(new PairInt(0, -1));
        zeroes.add(new PairInt(1, 2)); zeroes.add(new PairInt(1, 0)); zeroes.add(new PairInt(1, -1));
        zeroes.add(new PairInt(2, 2)); zeroes.add(new PairInt(2, 0));
        
        ones.add(new PairInt(-2, 1));  ones.add(new PairInt(-1, 1));
        ones.add(new PairInt(1, 1)); ones.add(new PairInt(2, 1));
        
        changeToZeroes.add(new PairInt(0, 0));
        
        changeToOnes.add(new PairInt(0, 1));
        
        int nCorrections = 0;
       
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        // ----- change the sign of y to handle other direction -----
        reverseYs(zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += replacePattern(
            points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }

    /**
     * corrects for the curve artifacts introduced after phase congruency edge 
     * thinning followed by phase angle step concatenation.
     * <pre>
     * the pattern is a single pixel displaced on a vertical or horizontal 
     * line
     *     
     *     #
     *     #
     *      #
     *     #
     *     #
     * </pre>
     * @param points
     * @param imageWidth
     * @param imageHeight 
     */
    public void correctForLineHatVert(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
                
        /*
                    .  #  .         2
                 .  .  #  .         1                          
                 .  #  .  .         0                          
                 .  .  #  .        -1                          
                    .  #  .        -2                         
          -3 -2 -1  0  1  2  3       
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, 0)); zeroes.add(new PairInt(-1, -1));
        zeroes.add(new PairInt(0, 2)); zeroes.add(new PairInt(0, 1)); 
        zeroes.add(new PairInt(0, -1)); zeroes.add(new PairInt(0, -2));
        zeroes.add(new PairInt(1, 0));
        zeroes.add(new PairInt(2, 2)); zeroes.add(new PairInt(2, 1));
        zeroes.add(new PairInt(2, 0)); zeroes.add(new PairInt(2, -1));
        zeroes.add(new PairInt(2, -2));
        
        ones.add(new PairInt(1, 2)); ones.add(new PairInt(1, 1)); ones.add(new PairInt(0, 0));
        ones.add(new PairInt(1, -1));  ones.add(new PairInt(1, -2));
        
        changeToZeroes.add(new PairInt(0, 0));
        
        changeToOnes.add(new PairInt(1, 0));
        
        int nCorrections = 0;

        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        // ----- change the sign of y to handle other direction -----
        reverseXs(zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += replacePattern(
            points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }

    /**
     * note, line should already be thinned to use this method.  it's intended
     * to remove the pixels with 3 neighbors.
     * @param lines
     * @param points
     * @param imageWidth
     * @param imageHeight
     */
    public void thinLineStaircases(Map<Set<PairInt>, PairInt> lines, 
        Set<PairInt> points, int imageWidth, int imageHeight) {
                
        /*
        looking for stair case patterns in each line that are thicker than 
        1 pixel.
        if found, determines longest connected segment to both sides and 
        trims the pixel from that side.
        
        for example:
                    would trim this pixel instead of pixel below it
                    \/
            # # # # # 
                    # # #
        */
        
        int[] dxs = new int[]{ 0,  1,  1,  1,  0, -1, -1, -1};
        int[] dys = new int[]{ 1,  1,  0, -1, -1, -1,  0,  1};
  
        List<Set<PairInt>> lineList = new ArrayList<Set<PairInt>>(lines.keySet());
        for (int i = 0; i < lineList.size(); ++i) {
            Set<PairInt> line = lineList.get(i);            
            Set<PairInt> rm = new HashSet<PairInt>();
            
            /*
             7 0 1
             6 # 2
             5 4 3
            */
            boolean[] present = new boolean[8];
            for (PairInt p : line) {                
                if (rm.contains(p)) {
                    continue;
                }
                int x = p.getX();
                int y = p.getY();
                int c = 0;
                for (int k = 0; k < dxs.length; ++k) {
                    int x2 = x + dxs[k];
                    int y2 = y + dys[k];
                    PairInt p2 = new PairInt(x2, y2);
                    if (line.contains(p2) && !rm.contains(p2)) {
                        ++c;
                        present[k] = true;
                    } else {
                        present[k] = false;
                    }
                }                
                if (c < 3) {
                    continue;
                }
                /*
                 7 0 1
                 6 # 2
                 5 4 3
                */
                //  #  #  #
                //        *#  #
                if (present[7] && present[0] && present[2] 
                    && !present[1] && present[5] && present[4] && !present[3]
                    && !present[6]
                    ) {
                    int left = 0;
                    int right = 0;
                    int y2 = y + 1;
                    for (int x2 = x; x2 > -1; --x2) {
                        PairInt p2 = new PairInt(x2, y2);
                        if (line.contains(p2) && !rm.contains(p2)) {
                            ++left;
                        } else {
                            break;
                        }
                    }
                    y2 = y;
                    for (int x2 = (x + 1); x2 < Integer.MAX_VALUE; ++x2) {
                        PairInt p2 = new PairInt(x2, y2);
                        if (line.contains(p2) && !rm.contains(p2)) {
                            ++right;
                        } else {
                            break;
                        }
                    }
                    if (left > right) {
                        rm.add(new PairInt(x, y + 1));
                    } else {
                        rm.add(p);
                    }
                } else if (present[6] && present[0] && present[1]
                    && !present[7] && !present[2]
                    && !present[5] && !present[4] && !present[3]
                    ) {
                    /*
                     7 0 1
                     6 # 2
                     5 4 3
                     */
                    //        #  #
                    //  #  #  #*
                    int left = 0;
                    int right = 0;
                    int y2 = y;
                    for (int x2 = (x - 1); x2 > -1; --x2) {
                        PairInt p2 = new PairInt(x2, y2);
                        if (line.contains(p2) && !rm.contains(p2)) {
                            ++left;
                        } else {
                            break;
                        }
                    }
                    y2 = y + 1;
                    for (int x2 = x; x2 < Integer.MAX_VALUE; ++x2) {
                        PairInt p2 = new PairInt(x2, y2);
                        if (line.contains(p2) && !rm.contains(p2)) {
                            ++right;
                        } else {
                            break;
                        }
                    }
                    if (right > left) {
                        rm.add(new PairInt(x, y + 1));
                    } else {
                        rm.add(p);
                    }
                } else if (present[5] && present[4] && present[2]
                    && !present[0] && !present[1] && !present[3] && !present[6] && !present[7]
                    ) {
                    /*
                     7 0 1
                     6 # 2
                     5 4 3
                     */
                    //        *#  #
                    //  #  #  #
                    
                    int left = 0;
                    int right = 0;
                    int y2 = y - 1;
                    for (int x2 = x; x2 > -1; --x2) {
                        PairInt p2 = new PairInt(x2, y2);
                        if (line.contains(p2) && !rm.contains(p2)) {
                            ++left;
                        } else {
                            break;
                        }
                    }
                    y2 = y;
                    for (int x2 = (x + 1); x2 < Integer.MAX_VALUE; ++x2) {
                        PairInt p2 = new PairInt(x2, y2);
                        if (line.contains(p2) && !rm.contains(p2)) {
                            ++right;
                        } else {
                            break;
                        }
                    }
                    if (left > right) {
                        rm.add(new PairInt(x, y - 1));
                    } else {
                        rm.add(p);
                    }
                } else if (present[6] && present[4] && present[3]
                    && !present[0] && !present[1] && !present[2]
                    && !present[5] && !present[7]
                    ) {
                    /*
                     7 0 1
                     6 # 2
                     5 4 3
                     */
                    //  #  #  #*
                    //        #  #
                    int left = 0;
                    int right = 0;
                    int y2 = y;
                    for (int x2 = x - 1; x2 > -1; --x2) {
                        PairInt p2 = new PairInt(x2, y2);
                        if (line.contains(p2) && !rm.contains(p2)) {
                            ++left;
                        } else {
                            break;
                        }
                    }
                    y2 = y - 1;
                    for (int x2 = x; x2 < Integer.MAX_VALUE; ++x2) {
                        PairInt p2 = new PairInt(x2, y2);
                        if (line.contains(p2) && !rm.contains(p2)) {
                            ++right;
                        } else {
                            break;
                        }
                    }
                    if (right > left) {
                        rm.add(new PairInt(x, y - 1));
                    } else {
                        rm.add(p);
                    }
                } else if (present[0] && present[2] && present[3]
                    && !present[1] && !present[4] && !present[5]
                    && !present[6] && !present[7]
                    ) {
                    /*
                     7 0 1
                     6 # 2
                     5 4 3
                     */
                    //     #
                    //     #
                    //    *#  #
                    //        #
                    int up = 0;
                    int down = 0;
                    int x2 = x;
                    for (int y2 = (y + 1); y2 < Integer.MAX_VALUE; ++y2) {
                        PairInt p2 = new PairInt(x2, y2);
                        if (line.contains(p2) && !rm.contains(p2)) {
                            ++up;
                        } else {
                            break;
                        }
                    }
                    x2 = x + 1;
                    for (int y2 = y; y2 > -1; --y2) {
                        PairInt p2 = new PairInt(x2, y2);
                        if (line.contains(p2) && !rm.contains(p2)) {
                            ++down;
                        } else {
                            break;
                        }
                    }                    
                    if (down > up) {
                        rm.add(new PairInt(x + 1, y));
                    } else {
                        rm.add(p);
                    }
                } else if (present[5] && present[6] && present[0]
                    && !present[1] && !present[2] && !present[3]
                    && !present[4] && !present[7]
                    ) {
                    /*
                     7 0 1
                     6 # 2
                     5 4 3
                     */
                    //     #
                    //     #
                    //  # *#   
                    //  #  
                    int up = 0;
                    int down = 0;
                    int x2 = x - 1;
                    for (int y2 = y; y2 > -1; --y2) {
                        PairInt p2 = new PairInt(x2, y2);
                        if (line.contains(p2) && !rm.contains(p2)) {
                            ++down;
                        } else {
                            break;
                        }
                    }
                    x2 = x;
                    for (int y2 = (y + 1); y2 < Integer.MAX_VALUE; ++y2) {
                        PairInt p2 = new PairInt(x2, y2);
                        if (line.contains(p2) && !rm.contains(p2)) {
                            ++up;
                        } else {
                            break;
                        }
                    }                    
                    if (down > up) {
                        rm.add(new PairInt(x - 1, y));
                    } else {
                        rm.add(p);
                    }
                } else if (present[7] && present[6] && present[4]
                    && !present[0] && !present[1] && !present[2]
                    && !present[3] && !present[5]
                    ) {
                    /*
                     7 0 1
                     6 # 2
                     5 4 3
                     */
                    //     #
                    //     #
                    //     #  #*
                    //        #
                    int up = 0;
                    int down = 0;
                    int x2 = x - 1;
                    for (int y2 = y; y2 < Integer.MAX_VALUE; ++y2) {
                        PairInt p2 = new PairInt(x2, y2);
                        if (line.contains(p2) && !rm.contains(p2)) {
                            ++up;
                        } else {
                            break;
                        }
                    }
                    x2 = x;
                    for (int y2 = (y - 1); y2 > -1; --y2) {
                        PairInt p2 = new PairInt(x2, y2);
                        if (line.contains(p2) && !rm.contains(p2)) {
                            ++down;
                        } else {
                            break;
                        }
                    }                    
                    if (up > down) {
                        rm.add(new PairInt(x - 1, y));
                    } else {
                        rm.add(p);
                    }
                } else if (present[4] && present[2] && present[1]
                    && !present[0] && !present[3] && !present[5]
                    && !present[6] && !present[7]) {
                    /*
                     7 0 1
                     6 # 2
                     5 4 3
                     */
                    //     #
                    //     #
                    // *#  #   
                    //  #  
                    int up = 0;
                    int down = 0;
                    int x2 = x;
                    for (int y2 = (y - 1); y2 > -1; --y2) {
                        PairInt p2 = new PairInt(x2, y2);
                        if (line.contains(p2) && !rm.contains(p2)) {
                            ++down;
                        } else {
                            break;
                        }
                    }
                    x2 = x + 1;
                    for (int y2 = y; y2 < Integer.MAX_VALUE; ++y2) {
                        PairInt p2 = new PairInt(x2, y2);
                        if (line.contains(p2) && !rm.contains(p2)) {
                            ++up;
                        } else {
                            break;
                        }
                    }                    
                    if (up > down) {
                        rm.add(new PairInt(x + 1, y));
                    } else {
                        rm.add(p);
                    }
                }
            }
            
            for (PairInt p : rm) {
                line.remove(p);
                points.remove(p);
            }
        }        
    }
    
    /**
     * corrects for the curve artifacts introduced after phase congruency edge 
     * thinning followed by phase angle step concatenation.
     * <pre>
     * the pattern is a single pixel displaced on a vertical or horizontal 
     * line
     *        #
     *     #######
     *     
     * </pre>
     * @param points
     * @param imageWidth
     * @param imageHeight 
     */
    public void correctForLineSpurHoriz(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
                
        /*
                    .  .  .         1                           2
                 .  .  #  .  .      0                           1
                 #  #  #  #  #     -1                           0
                 .  .  .  .  .     -2                          -1
             -3 -2 -1  0  1  2  3        -3 -2 -1  0  1  2  3
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-2, 2)); zeroes.add(new PairInt(-2, 0));
        zeroes.add(new PairInt(-1, 2)); zeroes.add(new PairInt(-1, 0)); zeroes.add(new PairInt(-1, -1));
        zeroes.add(new PairInt(0, 2)); zeroes.add(new PairInt(0, -1));
        zeroes.add(new PairInt(1, 2)); zeroes.add(new PairInt(1, 0)); zeroes.add(new PairInt(1, -1));
        zeroes.add(new PairInt(2, 2)); zeroes.add(new PairInt(2, 0));
        
        ones.add(new PairInt(-2, 1));  ones.add(new PairInt(-1, 1));
        ones.add(new PairInt(0, 1)); 
        ones.add(new PairInt(1, 1)); ones.add(new PairInt(2, 1));
        
        changeToZeroes.add(new PairInt(0, 0));
                
        int nCorrections = 0;
       
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        // ----- change the sign of y to handle other direction -----
        reverseYs(zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += replacePattern(
            points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    /**
     * corrects for the curve artifacts introduced after phase congruency edge 
     * thinning followed by phase angle step concatenation.
     * <pre>
     * the pattern is a single pixel displaced on a vertical or horizontal 
     * line
     *     
     *     #
     *     #
     *     # #
     *     #
     *     #
     * </pre>
     * @param points
     * @param imageWidth
     * @param imageHeight 
     */
    public void correctForLineSpurVert(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
                
        /*
                    .  #  .         2
                 .  .  #  .         1                          
                 .  #  #  .         0                          
                 .  .  #  .        -1                          
                    .  #  .        -2                         
          -3 -2 -1  0  1  2  3       
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, 0)); zeroes.add(new PairInt(-1, -1));
        zeroes.add(new PairInt(0, 2)); zeroes.add(new PairInt(0, 1)); 
        zeroes.add(new PairInt(0, -1)); zeroes.add(new PairInt(0, -2));
        zeroes.add(new PairInt(2, 2)); zeroes.add(new PairInt(2, 1));
        zeroes.add(new PairInt(2, 0)); zeroes.add(new PairInt(2, -1));
        zeroes.add(new PairInt(2, -2));
        
        ones.add(new PairInt(1, 0));
        ones.add(new PairInt(1, 2)); ones.add(new PairInt(1, 1)); ones.add(new PairInt(0, 0));
        ones.add(new PairInt(1, -1));  ones.add(new PairInt(1, -2));
        
        changeToZeroes.add(new PairInt(0, 0));
                
        int nCorrections = 0;

        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        // ----- change the sign of y to handle other direction -----
        reverseXs(zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += replacePattern(
            points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }

    /**
     * corrects for the curve artifacts introduced after phase congruency edge 
     * thinning followed by phase angle step concatenation.
     * <pre>
     * the pattern is a 2 pixel pixel spur displaced on a vertical or horizontal 
     * line
     *        #
     *        #
     *     #######
     *     
     * </pre>
     * @param points
     * @param imageWidth
     * @param imageHeight 
     */
    public void correctForLine2SpurVert(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
                
        /*
                 .  .  .  .  .      3
                 .  .  .  .  .      2
                 .  .  #  .  .      1                           2
                 .  .  #  .  .      0                           1
                 #  #  #  #  #     -1                           0
                 .  .  .  .  .     -2                          -1
             -3 -2 -1  0  1  2  3        -3 -2 -1  0  1  2  3
        */
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-2, 2)); zeroes.add(new PairInt(-2, 0));
        zeroes.add(new PairInt(-2, -1)); zeroes.add(new PairInt(-2, -2));
        zeroes.add(new PairInt(-2, -3));
        zeroes.add(new PairInt(-1, 2)); zeroes.add(new PairInt(-1, 0)); zeroes.add(new PairInt(-1, -1));
        zeroes.add(new PairInt(-1, -2)); zeroes.add(new PairInt(-1, -3));
        zeroes.add(new PairInt(0, 2)); zeroes.add(new PairInt(0, -2)); zeroes.add(new PairInt(0, -3));
        zeroes.add(new PairInt(1, 2)); zeroes.add(new PairInt(1, 0)); zeroes.add(new PairInt(1, -1));
        zeroes.add(new PairInt(1, -2)); zeroes.add(new PairInt(1, -3));
        zeroes.add(new PairInt(2, 2)); zeroes.add(new PairInt(2, 0));
        zeroes.add(new PairInt(2, -1)); zeroes.add(new PairInt(2, -2));
        zeroes.add(new PairInt(2, -3));
        
        ones.add(new PairInt(-2, 1));  ones.add(new PairInt(-1, 1));
        ones.add(new PairInt(0, 1)); ones.add(new PairInt(0, 0));  ones.add(new PairInt(0, -1));
        ones.add(new PairInt(1, 1)); ones.add(new PairInt(2, 1));
        
        changeToZeroes.add(new PairInt(0, 0));
        changeToZeroes.add(new PairInt(0, -1));
                
        int nCorrections = 0;
       
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        // ----- change the sign of y to handle other direction -----
        reverseYs(zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += replacePattern(
            points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    /**
     * corrects for the curve artifacts introduced after phase congruency edge 
     * thinning followed by phase angle step concatenation.
     * <pre>
     * the pattern is a 2 pixel spur displaced on a vertical or horizontal 
     * line
     *     
     *     #
     *     #
     *     # # #
     *     #
     *     #
     * </pre>
     * @param points
     * @param imageWidth
     * @param imageHeight 
     */
    public void correctForLine2SpurHoriz(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
        
        //TODO: may want to consider extending the vertical line out 3 pix on each side
                
        /*
           .  .  .  .  #  .         2
           .  .  .  .  #  .         1                          
           .  .  #  #  #  .         0                          
           .  .  .  .  #  .        -1                          
           .  .  .  .  #  .        -2                         
          -3 -2 -1  0  1  2  3       
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        for (int i = 0; i < 5; ++i) {
            zeroes.add(new PairInt(-3, i - 2));
            zeroes.add(new PairInt(-2, i - 2));
            zeroes.add(new PairInt(2, i - 2));
            ones.add(new PairInt(1, i - 2));
        }
        zeroes.add(new PairInt(-1, 2)); zeroes.add(new PairInt(-1, 1)); 
        zeroes.add(new PairInt(-1, -1)); zeroes.add(new PairInt(-1, -2));
        zeroes.add(new PairInt(0, 2)); zeroes.add(new PairInt(0, 1)); 
        zeroes.add(new PairInt(0, -1)); zeroes.add(new PairInt(0, -2));
        
        ones.add(new PairInt(-1, 0)); ones.add(new PairInt(0, 0));
        
        changeToZeroes.add(new PairInt(-1, 0)); changeToZeroes.add(new PairInt(0, 0));
                
        int nCorrections = 0;

        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        // ----- change the sign of y to handle other direction -----
        reverseXs(zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += replacePattern(
            points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }

}
