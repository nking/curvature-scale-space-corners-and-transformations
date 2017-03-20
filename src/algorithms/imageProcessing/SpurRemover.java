package algorithms.imageProcessing;

import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import java.util.HashSet;
import java.util.Set;

/**
 * and incomplete class to remove single pixels spur patterns from curves, 
 * with the assumption that the given points are curves.  It also corrects
 * a few features to make the main use of the results easier, that is,
 * assembly of the border points of a curve into an ordered sequence around
 * the perimeter.
 * 
 * @author nichole
 */
public class SpurRemover {
    
    private boolean debug = false;
    
    public void setToDebug() {
        debug = true;
    }
    
    /**
     * an incomplete set of patterns are applied to the curve to remove spurs
     * that are longer than a pixel and are 1 pixel wide.  this method is not
     * complete. 
     * 
     * @param img 
     */
    public void remove(GreyscaleImage img) {
      
        Set<PairInt> points = new HashSet<PairInt>();
        for (int i = 0; i < img.getNPixels(); ++i) {
            if (img.getValue(i) > 0) {
                points.add(new PairInt(img.getCol(i), img.getRow(i)));
            }
        }
        
        Set<PairInt> cp = new HashSet<PairInt>(points);
        
        remove(points, img.getWidth(), img.getHeight());
        
        cp.removeAll(points);
        
        for (PairInt p : cp) {
            img.setValue(p.getX(), p.getY(), 0);
        }
    }
    
    /**
     * an incomplete set of patterns are applied to the curve to remove spurs
     * that are longer than a pixel and are 1 pixel wide.  this method is not
     * complete.
     * 
     * @param points 
     */
    public void remove(Set<PairInt> points, int imageWidth, int imageHeight) {
        
        int nChanged = 0;
        int nMaxIter = 10;
        int nIter = 0;

        while ((nIter == 0) || ((nChanged > 0) && (nIter < nMaxIter))) {
        
//if (debug) {        
Image img3 = new Image(imageWidth, imageHeight);
for (PairInt p : points) {
    img3.setRGB(p.getX(), p.getY(), 255, 0, 0);
}
MiscDebug.writeImage(img3, "spur_removal_" + nIter + "_" 
+ MiscDebug.getCurrentTimeFormatted());
//}

            nChanged = 0;
            
            nChanged += pattern0(points, imageWidth, imageHeight);

            nChanged += pattern00(points, imageWidth, imageHeight);

            nChanged += pattern1(points, imageWidth, imageHeight);

            nChanged += pattern2(points, imageWidth, imageHeight);

            nChanged += pattern3(points, imageWidth, imageHeight);

            nChanged += pattern4(points, imageWidth, imageHeight);
            
            nChanged += pattern5(points, imageWidth, imageHeight);
            
            nChanged += pattern6(points, imageWidth, imageHeight);
            
            nChanged += pattern7(points, imageWidth, imageHeight);
            
            nChanged += pattern8(points, imageWidth, imageHeight);
            
            nChanged += pattern9(points, imageWidth, imageHeight);
            
            nChanged += pattern10(points, imageWidth, imageHeight);
            
            nChanged += pattern11(points, imageWidth, imageHeight);
            
            nChanged += pattern13(points, imageWidth, imageHeight);
           
            nChanged += curveCorrection12(points, imageWidth, imageHeight);
            
            ++nIter;
        }
        
if (debug) {        
Image img3 = new Image(imageWidth, imageHeight);
for (PairInt p : points) {
    img3.setRGB(p.getX(), p.getY(), 255, 0, 0);
}
MiscDebug.writeImage(img3, "spur_removal_" + nIter + "_" 
+ MiscDebug.getCurrentTimeFormatted());
}

    }
    
    public boolean isASpur(int xCoord, int yCoord, Set<PairInt> curve, 
        int imageWidth, int imageHeight) {
        
        boolean hasPattern = testPatternSwapDirections(xCoord, yCoord, curve, 
            imageWidth, imageHeight, getPattern0());
        
        if (hasPattern) {
            return true;
        }

        hasPattern = testPatternSwapDirections(xCoord, yCoord, curve, 
            imageWidth, imageHeight, getPattern00());
        
        if (hasPattern) {
            return true;
        }
        
        hasPattern = testPatternSwapDirections(xCoord, yCoord, curve, 
            imageWidth, imageHeight, getPattern1());
        
        if (hasPattern) {
            return true;
        }

        hasPattern = testPatternSwapDirections(xCoord, yCoord, curve, 
            imageWidth, imageHeight, getPattern2());
        
        if (hasPattern) {
            return true;
        }

        hasPattern = testPatternSwapDirections(xCoord, yCoord, curve, 
            imageWidth, imageHeight, getPattern3());
        
        if (hasPattern) {
            return true;
        }

        hasPattern = testPatternSwapDirections(xCoord, yCoord, curve, 
            imageWidth, imageHeight, getPattern4());
        
        if (hasPattern) {
            return true;
        }
            
        hasPattern = testPatternSwapDirections(xCoord, yCoord, curve, 
            imageWidth, imageHeight, getPattern5());
        
        if (hasPattern) {
            return true;
        }
            
        hasPattern = testPatternSwapDirections(xCoord, yCoord, curve, 
            imageWidth, imageHeight, getPattern6());
        
        if (hasPattern) {
            return true;
        }
            
        hasPattern = testPatternSwapDirections(xCoord, yCoord, curve, 
            imageWidth, imageHeight, getPattern7());
        
        if (hasPattern) {
            return true;
        }
            
        hasPattern = testPatternSwapDirections(xCoord, yCoord, curve, 
            imageWidth, imageHeight, getPattern8());
        
        if (hasPattern) {
            return true;
        }
            
        hasPattern = testPatternSwapDirections(xCoord, yCoord, curve, 
            imageWidth, imageHeight, getPattern9());
        
        if (hasPattern) {
            return true;
        }
            
        hasPattern = testPatternSwapDirections(xCoord, yCoord, curve, 
            imageWidth, imageHeight, getPattern10());
        
        if (hasPattern) {
            return true;
        }
            
        hasPattern = testPatternSwapDirections(xCoord, yCoord, curve, 
            imageWidth, imageHeight, getPattern11());
        
        if (hasPattern) {
            return true;
        }
        
        hasPattern = testPatternSwapDirections(xCoord, yCoord, curve, 
            imageWidth, imageHeight, getPattern13());
        
        if (hasPattern) {
            return true;
        }
            
        hasPattern = testPatternSwapDirections(xCoord, yCoord, curve, 
            imageWidth, imageHeight, getPattern12());
        
        return hasPattern;
    }
    
    public static class PatternReplacement {
        Set<PairInt> ones;
        Set<PairInt> zeroes;
        Set<PairInt> changeToZeroes;
    }
    
    private PatternReplacement getPattern0() {
        
        /*
          0  #  0        1
          0 >#  0        0
          0  0  0       -1
         -1  0  1  2  3
         */
        PatternReplacement pr = new PatternReplacement();
        pr.ones = new HashSet<PairInt>();
        pr.zeroes = new HashSet<PairInt>();
        pr.changeToZeroes = new HashSet<PairInt>();

        pr.zeroes.add(new PairInt(-1, -1)); pr.zeroes.add(new PairInt(-1, 0)); pr.zeroes.add(new PairInt(-1, 1));
        pr.zeroes.add(new PairInt(0, -1));
        pr.zeroes.add(new PairInt(1, -1)); pr.zeroes.add(new PairInt(1, 0)); pr.zeroes.add(new PairInt(1, 1));
        
        pr.ones.add(new PairInt(0, 1)); 
        
        pr.changeToZeroes.add(new PairInt(0, 0));
        
        return pr;
    }

    private int pattern0(Set<PairInt> curve, int imageWidth, int imageHeight) {
        
        PatternReplacement pr = getPattern0();
        
        int nChanged = 0;
        
        /*
        invoke pattern and then reversed for Y
        */

        nChanged += replacePattern(curve, imageWidth, imageHeight, pr);

        // ----- change the sign of x to handle other direction -----
        reverseYs(pr);

        nChanged += replacePattern(curve, imageWidth, imageHeight, pr);
        
        return nChanged;
    }
    
    private PatternReplacement getPattern00() {
         
        /*
          0  0  0        1
          0 >#  #        0
          0  0  0       -1
         -1  0  1  2  3
         */
        PatternReplacement pr = new PatternReplacement();
        pr.ones = new HashSet<PairInt>();
        pr.zeroes = new HashSet<PairInt>();
        pr.changeToZeroes = new HashSet<PairInt>();

        pr.zeroes.add(new PairInt(-1, -1)); pr.zeroes.add(new PairInt(-1, 0)); pr.zeroes.add(new PairInt(-1, 1));
        pr.zeroes.add(new PairInt(0, -1)); pr.zeroes.add(new PairInt(0, 1));
        pr.zeroes.add(new PairInt(1, -1)); pr.zeroes.add(new PairInt(1, 1));
        
        pr.ones.add(new PairInt(1, 0)); 
        
        pr.changeToZeroes.add(new PairInt(0, 0));
        
        return pr;
    }
    
    private int pattern00(Set<PairInt> curve, int imageWidth, int imageHeight) {
       
        int nChanged = 0;
        
        PatternReplacement pr = getPattern00();
        
        nChanged += replacePattern(curve, imageWidth, imageHeight, pr);
        
        // ----- change the sign of x to handle other direction -----
        reverseXs(pr);

        nChanged += replacePattern(curve, imageWidth, imageHeight, pr);
            
        return nChanged;
    }
    
    private PatternReplacement getPattern1() {
        
        /*
          #  #  0        1
          0 >#  0        0
          0  0  0       -1
         -1  0  1  2  3
         */
        PatternReplacement pr = new PatternReplacement();
        pr.ones = new HashSet<PairInt>();
        pr.zeroes = new HashSet<PairInt>();
        pr.changeToZeroes = new HashSet<PairInt>();

        pr.zeroes.add(new PairInt(-1, -1)); pr.zeroes.add(new PairInt(-1, 0));
        pr.zeroes.add(new PairInt(0, -1));
        pr.zeroes.add(new PairInt(1, -1)); pr.zeroes.add(new PairInt(1, 0)); pr.zeroes.add(new PairInt(1, 1));
        
        pr.ones.add(new PairInt(-1, 1));
        pr.ones.add(new PairInt(0, 1)); 
        
        pr.changeToZeroes.add(new PairInt(0, 0));
        
        return pr;
    }
    
    private int pattern1(Set<PairInt> curve, int imageWidth, int imageHeight) {
        
        int nChanged = 0;
        
        PatternReplacement pr = getPattern1();
        
        nChanged = replacePatternSwapDirections(curve, imageWidth, imageHeight, 
            pr);
       
        return nChanged;
    }
    
    private PatternReplacement getPattern2() {
         
        /*
          #              2
          0  #  0        1
          0 >#  0        0
          0  0  0       -1
         -1  0  1  2  3
         */
        PatternReplacement pr = new PatternReplacement();
        pr.ones = new HashSet<PairInt>();
        pr.zeroes = new HashSet<PairInt>();
        pr.changeToZeroes = new HashSet<PairInt>();

        pr.zeroes.add(new PairInt(-1, -1)); pr.zeroes.add(new PairInt(-1, 0)); pr.zeroes.add(new PairInt(-1, 1));
        pr.zeroes.add(new PairInt(0, -1));
        pr.zeroes.add(new PairInt(1, -1)); pr.zeroes.add(new PairInt(1, 0)); pr.zeroes.add(new PairInt(1, 1));
        
        pr.ones.add(new PairInt(-1, 2));
        pr.ones.add(new PairInt(0, 1)); 
        
        pr.changeToZeroes.add(new PairInt(0, 0));
        
        return pr;
    }
    
    private int pattern2(Set<PairInt> curve, int imageWidth, int imageHeight) {
       
        int nChanged = 0;
                
        PatternReplacement pr = getPattern2();
        
        nChanged = replacePatternSwapDirections(curve, imageWidth, imageHeight, 
            pr);
           
        return nChanged;
    }
    
    private PatternReplacement getPattern3() {
          
        /*
          0  0  #        1
          0 >#  #        0
          0  0  0       -1
         -1  0  1  2  3
         */
        PatternReplacement pr = new PatternReplacement();
        pr.ones = new HashSet<PairInt>();
        pr.zeroes = new HashSet<PairInt>();
        pr.changeToZeroes = new HashSet<PairInt>();

        pr.zeroes.add(new PairInt(-1, -1)); pr.zeroes.add(new PairInt(-1, 0)); pr.zeroes.add(new PairInt(-1, 1));
        pr.zeroes.add(new PairInt(0, -1)); pr.zeroes.add(new PairInt(0, 1));
        pr.zeroes.add(new PairInt(1, -1));
        
        pr.ones.add(new PairInt(1, 0)); pr.ones.add(new PairInt(1, 1)); 
        
        pr.changeToZeroes.add(new PairInt(0, 0));
        
        return pr;
    }
    
    private int pattern3(Set<PairInt> curve, int imageWidth, int imageHeight) {
      
        int nChanged = 0;
        
        PatternReplacement pr = getPattern3();
        
        nChanged = replacePatternSwapDirections(curve, imageWidth, imageHeight, 
            pr);
        
        return nChanged;
    }
    
    private PatternReplacement getPattern4() {
          
        /*
          0  0  0  #     1
          0 >#  #        0
          0  0  0       -1
         -1  0  1  2  3
         */
        PatternReplacement pr = new PatternReplacement();
        pr.ones = new HashSet<PairInt>();
        pr.zeroes = new HashSet<PairInt>();
        pr.changeToZeroes = new HashSet<PairInt>();

        pr.zeroes.add(new PairInt(-1, -1)); pr.zeroes.add(new PairInt(-1, 0)); pr.zeroes.add(new PairInt(-1, 1));
        pr.zeroes.add(new PairInt(0, -1)); pr.zeroes.add(new PairInt(0, 1));
        pr.zeroes.add(new PairInt(1, -1)); pr.zeroes.add(new PairInt(1, 1));
        
        pr.ones.add(new PairInt(1, 0)); 
        pr.ones.add(new PairInt(2, 1)); 
        
        pr.changeToZeroes.add(new PairInt(0, 0));
        
        return pr;
    }
    
    private int pattern4(Set<PairInt> curve, int imageWidth, int imageHeight) {
      
        int nChanged = 0;
                   
        PatternReplacement pr = getPattern4();
        
        nChanged = replacePatternSwapDirections(curve, imageWidth, imageHeight, 
            pr);
        
        return nChanged;
    }
    
    private PatternReplacement getPattern5() {
          
        /*
             0  0  0     2
          0  0  #  #     1
          0 >#  0        0
          0  0  0       -1
         -1  0  1  2  3
         */
        PatternReplacement pr = new PatternReplacement();
        pr.ones = new HashSet<PairInt>();
        pr.zeroes = new HashSet<PairInt>();
        pr.changeToZeroes = new HashSet<PairInt>();

        pr.zeroes.add(new PairInt(-1, -1)); pr.zeroes.add(new PairInt(-1, 0)); pr.zeroes.add(new PairInt(-1, 1));
        pr.zeroes.add(new PairInt(0, -1)); pr.zeroes.add(new PairInt(0, 1)); pr.zeroes.add(new PairInt(0, 2));
        pr.zeroes.add(new PairInt(1, -1)); pr.zeroes.add(new PairInt(1, 0)); pr.zeroes.add(new PairInt(1, 2));
        pr.zeroes.add(new PairInt(2, 2));
        
        pr.ones.add(new PairInt(1, 1)); 
        pr.ones.add(new PairInt(2, 1)); 
        
        pr.changeToZeroes.add(new PairInt(0, 0));
        
        return pr;
    }
    
    private int pattern5(Set<PairInt> curve, int imageWidth, int imageHeight) {
      
        int nChanged = 0;
                        
        PatternReplacement pr = getPattern5();
        
        nChanged = replacePatternSwapDirections(curve, imageWidth, imageHeight, 
            pr);
        
        return nChanged;
    }
    
    private PatternReplacement getPattern6() {
         /*
                         
             #           2
          0  #  0  0     1
          0  0 >#  0     0
          0  0  0  0    -1
         -2 -1  0  1  2  3
         */
        PatternReplacement pr = new PatternReplacement();
        pr.ones = new HashSet<PairInt>();
        pr.zeroes = new HashSet<PairInt>();
        pr.changeToZeroes = new HashSet<PairInt>();

        pr.zeroes.add(new PairInt(-2, -1)); pr.zeroes.add(new PairInt(-2, 0)); pr.zeroes.add(new PairInt(-2, 1));
        pr.zeroes.add(new PairInt(-1, -1)); pr.zeroes.add(new PairInt(-1, 0));
        pr.zeroes.add(new PairInt(0, -1)); pr.zeroes.add(new PairInt(0, 1));
        pr.zeroes.add(new PairInt(1, -1)); pr.zeroes.add(new PairInt(1, 0)); pr.zeroes.add(new PairInt(1, 1));
        
        pr.ones.add(new PairInt(-1, 1)); pr.ones.add(new PairInt(-1, 2));
        
        pr.changeToZeroes.add(new PairInt(0, 0));
        
        return pr;
    }
    
    private int pattern6(Set<PairInt> curve, int imageWidth, int imageHeight) {
        
        PatternReplacement pr = getPattern6();
       
        int nChanged = 0;
                        
        nChanged = replacePatternSwapDirections(curve, imageWidth, imageHeight, 
            pr);
        
        return nChanged;
    }
    
    private PatternReplacement getPattern7() {
        
        /*
             0  0        2
          0  0  #  #     1
          0 >#  0  0     0
          0  0  0       -1
         -1  0  1  2  3
         */
        PatternReplacement pr = new PatternReplacement();
        pr.ones = new HashSet<PairInt>();
        pr.zeroes = new HashSet<PairInt>();
        pr.changeToZeroes = new HashSet<PairInt>();

        pr.zeroes.add(new PairInt(-1, -1)); pr.zeroes.add(new PairInt(-1, 0)); pr.zeroes.add(new PairInt(-1, 1));
        pr.zeroes.add(new PairInt(0, -1)); pr.zeroes.add(new PairInt(0, 1)); pr.zeroes.add(new PairInt(0, 2));
        pr.zeroes.add(new PairInt(1, -1)); pr.zeroes.add(new PairInt(1, 0)); pr.zeroes.add(new PairInt(1, 2));
        pr.zeroes.add(new PairInt(2, 0));
        
        pr.ones.add(new PairInt(1, 1)); 
        pr.ones.add(new PairInt(2, 1)); 
        
        pr.changeToZeroes.add(new PairInt(0, 0));
        
        return pr;
    }
    
    private int pattern7(Set<PairInt> curve, int imageWidth, int imageHeight) {
        
        PatternReplacement pr = getPattern7();
        
        int nChanged = 0;
                        
        nChanged = replacePatternSwapDirections(curve, imageWidth, imageHeight, 
            pr);
        
        return nChanged;
    }
    
    private PatternReplacement getPattern8() {
        /*
                         2
             #  0        1
          0  #  0  0     0
          0  0 >#  0    -1
          0  0  0  0    -2
         -1  0  1  2  3
         */
        PatternReplacement pr = new PatternReplacement();
        pr.ones = new HashSet<PairInt>();
        pr.zeroes = new HashSet<PairInt>();
        pr.changeToZeroes = new HashSet<PairInt>();

        pr.zeroes.add(new PairInt(-1, -1)); pr.zeroes.add(new PairInt(-1, 0)); pr.zeroes.add(new PairInt(-1, 1));
        pr.zeroes.add(new PairInt(0, -2)); pr.zeroes.add(new PairInt(0, -1));
        pr.zeroes.add(new PairInt(1, -2)); pr.zeroes.add(new PairInt(1, 0)); pr.zeroes.add(new PairInt(1, 1));
        pr.zeroes.add(new PairInt(2, -2)); pr.zeroes.add(new PairInt(2, -1)); pr.zeroes.add(new PairInt(2, 0));
        
        pr.ones.add(new PairInt(0, 1)); 
        pr.ones.add(new PairInt(1, -1)); 
        
        pr.changeToZeroes.add(new PairInt(1, -1));
        
        return pr;
    }
    
    private int pattern8(Set<PairInt> curve, int imageWidth, int imageHeight) {
        
        PatternReplacement pr = getPattern8();
        
        int nChanged = 0;
                        
        nChanged = replacePatternSwapDirections(curve, imageWidth, imageHeight, 
            pr);
        
        return nChanged;
    }
    
    private PatternReplacement getPattern9() {
         
        /*
                         2
          #  0  0        1
          0  #  0  0     0
          #  0 >#  0    -1
          0  0  0  0    -2
         -1  0  1  2  3
         */
        PatternReplacement pr = new PatternReplacement();
        pr.ones = new HashSet<PairInt>();
        pr.zeroes = new HashSet<PairInt>();
        pr.changeToZeroes = new HashSet<PairInt>();

        pr.zeroes.add(new PairInt(-1, -2)); pr.zeroes.add(new PairInt(-1, 0));
        pr.zeroes.add(new PairInt(0, -2)); pr.zeroes.add(new PairInt(0, -1)); pr.zeroes.add(new PairInt(0, 1));
        pr.zeroes.add(new PairInt(1, -2)); pr.zeroes.add(new PairInt(1, 0)); pr.zeroes.add(new PairInt(1, 1));
        pr.zeroes.add(new PairInt(2, -2)); pr.zeroes.add(new PairInt(2, -1)); pr.zeroes.add(new PairInt(2, 0));
        
        pr.ones.add(new PairInt(-1, -1)); pr.ones.add(new PairInt(-1, 1));
        pr.ones.add(new PairInt(1, -1)); 
        
        pr.changeToZeroes.add(new PairInt(1, -1));
        
        return pr;
    }
    
    private int pattern9(Set<PairInt> curve, int imageWidth, int imageHeight) {
       
        int nChanged = 0;
                  
        PatternReplacement pr = getPattern9();
        
        nChanged = replacePatternSwapDirections(curve, imageWidth, imageHeight, 
            pr);
        
        return nChanged;
    }
    
    private PatternReplacement getPattern10() {
         
        /*
                         2
             #  0  #        1
             0  #  0  0     0
             0  0  #  0    -1
             0  0  0  0    -2
            -1  0  1  2  3
         */
        PatternReplacement pr = new PatternReplacement();
        pr.ones = new HashSet<PairInt>();
        pr.zeroes = new HashSet<PairInt>();
        pr.changeToZeroes = new HashSet<PairInt>();

        pr.zeroes.add(new PairInt(-1, -2)); pr.zeroes.add(new PairInt(-1, 0));
        pr.zeroes.add(new PairInt(0, -2)); pr.zeroes.add(new PairInt(0, -1)); pr.zeroes.add(new PairInt(0, 1));
        pr.zeroes.add(new PairInt(1, -2)); pr.zeroes.add(new PairInt(1, -1)); pr.zeroes.add(new PairInt(1, 0));
        pr.zeroes.add(new PairInt(2, -2)); pr.zeroes.add(new PairInt(2, -1)); pr.zeroes.add(new PairInt(2, 0));
        
        pr.ones.add(new PairInt(-1, 1));
        pr.ones.add(new PairInt(1, -1)); pr.ones.add(new PairInt(1, 1)); 
        
        pr.changeToZeroes.add(new PairInt(1, -1));
        
        return pr;
    }
    
    private int pattern10(Set<PairInt> curve, int imageWidth, int imageHeight) {
       
        int nChanged = 0;
                    
        PatternReplacement pr = getPattern10();
        
        nChanged = replacePatternSwapDirections(curve, imageWidth, imageHeight, 
            pr);
        
        return nChanged;
    }
    
    private PatternReplacement getPattern11() {
         
        /*
             0  0  #        2
             #  #  0  0     1
             0  0 >#  0     0
                0  0  0    -1
            -2 -1  0  1  2  3
         */
        PatternReplacement pr = new PatternReplacement();
        pr.ones = new HashSet<PairInt>();
        pr.zeroes = new HashSet<PairInt>();
        pr.changeToZeroes = new HashSet<PairInt>();

        pr.zeroes.add(new PairInt(-2, 0)); pr.zeroes.add(new PairInt(-2, 2));
        pr.zeroes.add(new PairInt(-1, -1)); pr.zeroes.add(new PairInt(-1, 0)); pr.zeroes.add(new PairInt(-1, 2));
        pr.zeroes.add(new PairInt(0, -1)); pr.zeroes.add(new PairInt(0, 1));
        pr.zeroes.add(new PairInt(1, -1)); pr.zeroes.add(new PairInt(1, 0)); pr.zeroes.add(new PairInt(1, 1));
        
        pr.ones.add(new PairInt(-2, 1));
        pr.ones.add(new PairInt(-1, 1)); 
        pr.ones.add(new PairInt(0, 2)); 
        
        pr.changeToZeroes.add(new PairInt(0, 0));
        
        return pr;
    }
    
    private int pattern11(Set<PairInt> curve, int imageWidth, int imageHeight) {
       
        int nChanged = 0;
                    
        PatternReplacement pr = getPattern11();
        
        nChanged = replacePatternSwapDirections(curve, imageWidth, imageHeight, 
            pr);
        
        return nChanged;
    }
    
    private PatternReplacement getPattern12() {
        
        /*
                0  #  0     1
                0 >#  #     0
                0  0  0    -1
               -1  0  1  2  3
         */
        PatternReplacement pr = new PatternReplacement();
        pr.ones = new HashSet<PairInt>();
        pr.zeroes = new HashSet<PairInt>();
        pr.changeToZeroes = new HashSet<PairInt>();

        pr.zeroes.add(new PairInt(-1, -1)); pr.zeroes.add(new PairInt(-1, 0)); pr.zeroes.add(new PairInt(-1, 1));
        pr.zeroes.add(new PairInt(0, -1));
        pr.zeroes.add(new PairInt(1, -1)); pr.zeroes.add(new PairInt(1, 1));
        
        pr.ones.add(new PairInt(0, 1));
        pr.ones.add(new PairInt(1, 0)); 
        
        pr.changeToZeroes.add(new PairInt(0, 0));
        
        return pr;
    }
    
    private int curveCorrection12(Set<PairInt> curve, int imageWidth, int imageHeight) {
       
        int nChanged = 0;
                    
        PatternReplacement pr = getPattern12();
        
        nChanged = replacePatternSwapDirections(curve, imageWidth, imageHeight, 
            pr);
        
        return nChanged;
    }
    
    private PatternReplacement getPattern13() {
    
        /*
                            2
                0  0  0     1
                0 >#  0     0
                0  0  #    -1
            -2 -1  0  1  2  3
         */
        PatternReplacement pr = new PatternReplacement();
        pr.ones = new HashSet<PairInt>();
        pr.zeroes = new HashSet<PairInt>();
        pr.changeToZeroes = new HashSet<PairInt>();

        pr.zeroes.add(new PairInt(-1, -1)); pr.zeroes.add(new PairInt(-1, 0)); pr.zeroes.add(new PairInt(-1, 1));
        pr.zeroes.add(new PairInt(0, -1)); pr.zeroes.add(new PairInt(0, 1));
        pr.zeroes.add(new PairInt(1, 0)); pr.zeroes.add(new PairInt(1, 1));
        
        pr.ones.add(new PairInt(1, -1)); 
        
        pr.changeToZeroes.add(new PairInt(0, 0));
        
        return pr;
    }
    
    private int pattern13(Set<PairInt> curve, int imageWidth, int imageHeight) {
       
        int nChanged = 0;
                    
        PatternReplacement pr = getPattern13();
        
        nChanged = replacePatternSwapDirections(curve, imageWidth, imageHeight, 
            pr);
        
        return nChanged;
    }
    
    private int replacePattern(Set<PairInt> points, int imageWidth, 
        int imageHeight, final PatternReplacement pr) {
        
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
            
            for (PairInt p2 : pr.zeroes) {
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

            for (PairInt p2 : pr.ones) {
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
            
            for (PairInt p2 : pr.changeToZeroes) {
                int x = col + p2.getX();
                int y = row + p2.getY();
                if ((x < 0) || (y < 0) || (x > (w - 1)) || (y > (h - 1))) {
                    continue;
                }
                PairInt p3 = new PairInt(x, y);
                tmpPointsAdded.remove(p3);                
                tmpPointsRemoved.add(p3);
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

    private int replacePatternSwapDirections(Set<PairInt> curve, 
        int imageWidth, int imageHeight, PatternReplacement pr) {
        
        int nChanged = replacePattern(curve, imageWidth, imageHeight, pr);
        
        // ----- change the sign of x to handle other direction -----
        reverseXs(pr);
        
        nChanged += replacePattern(curve, imageWidth, imageHeight, pr);
        
        // ----- change the sign of y to handle other direction -----
        reverseYs(pr);
        
        nChanged = replacePattern(curve, imageWidth, imageHeight, pr);
        
        // ----- change the sign of x to handle other direction -----
        reverseXs(pr);
        
        nChanged += replacePattern(curve, imageWidth, imageHeight, pr);
        
        // ----- change the sign of y to handle other direction -----
        reverseYs(pr);
        
        return nChanged;
    }
    
    protected void reverseXs(final PatternReplacement pr) {
        
        // ----- change the sign of x to handle other direction -----
        for (PairInt p : pr.zeroes) {
            p.setX(-1 * p.getX());
        }
        
        for (PairInt p : pr.ones) {
            p.setX(-1 * p.getX());
        }
          
        for (PairInt p : pr.changeToZeroes) {
            p.setX(-1 * p.getX());
        }
    }
    
    protected void reverseYs(final PatternReplacement pr) {
        
        // ----- change the sign of y  -----
        for (PairInt p : pr.zeroes) {
            p.setY(-1 * p.getY());
        }
        
        for (PairInt p : pr.ones) {
            p.setY(-1 * p.getY());
        }
        
        for (PairInt p : pr.changeToZeroes) {
            p.setY(-1 * p.getY());
        }
    }
    
    private boolean testPatternSwapDirections(
        final int xCoord, final int yCoord, Set<PairInt> curve, 
        int imageWidth, int imageHeight, PatternReplacement pr) {
        
        boolean hasPattern = testPattern(xCoord, yCoord, curve, imageWidth, 
            imageHeight, pr);
        
        if (hasPattern) {
            return hasPattern;
        }
        
        // ----- change the sign of x to handle other direction -----
        reverseXs(pr);
        
        hasPattern = testPattern(xCoord, yCoord, curve, imageWidth, 
            imageHeight, pr);
        
        if (hasPattern) {
            return hasPattern;
        }
        
        // ----- change the sign of y to handle other direction -----
        reverseYs(pr);
        
        hasPattern = testPattern(xCoord, yCoord, curve, imageWidth, 
            imageHeight, pr);
        
        if (hasPattern) {
            return hasPattern;
        }
        
        // ----- change the sign of x to handle other direction -----
        reverseXs(pr);
        
        hasPattern = testPattern(xCoord, yCoord, curve, imageWidth, 
            imageHeight, pr);
        
        if (hasPattern) {
            return hasPattern;
        }
        
        // ----- change the sign of y to handle other direction -----
        reverseYs(pr);
        
        return hasPattern;
    }
    
    private boolean testPattern(int xCoord, int yCoord, Set<PairInt> points, 
        int imageWidth, int imageHeight, final PatternReplacement pr) {
        
        int w = imageWidth;
        int h = imageHeight;

        Set<PairInt> tmpPointsRemoved = new HashSet<PairInt>();

        PairInt p = new PairInt(xCoord, yCoord);
            
        if (!points.contains(p)) {
            return false;
        }
            
        boolean foundPattern = true;

        for (PairInt p2 : pr.zeroes) {
            int x = xCoord + p2.getX();
            int y = yCoord + p2.getY();
            if ((x < 0) || (y < 0) || (x > (w - 1)) || (y > (h - 1))) {
                //TODO: revisit this
                foundPattern = false;
                break;
            }
            PairInt p3 = new PairInt(x, y);
            if (!tmpPointsRemoved.contains(p3) && points.contains(p3)) {
                foundPattern = false;
                break;
            }
        }
               
        if (!foundPattern) {
            return false;
        }

        for (PairInt p2 : pr.ones) {
            int x = xCoord + p2.getX();
            int y = yCoord + p2.getY();
            if ((x < 0) || (y < 0) || (x > (w - 1)) || (y > (h - 1))) {
                foundPattern = false;
                break;
            }
            PairInt p3 = new PairInt(x, y);
            if (tmpPointsRemoved.contains(p3) || !points.contains(p3)) {
                foundPattern = false;
                break;
            }
        }
          
        if (!foundPattern) {
            return false;
        }
            
        for (PairInt p2 : pr.changeToZeroes) {
            int x = xCoord + p2.getX();
            int y = yCoord + p2.getY();
            if ((x < 0) || (y < 0) || (x > (w - 1)) || (y > (h - 1))) {
                continue;
            }
            PairInt p3 = new PairInt(x, y);
            tmpPointsRemoved.add(p3);
        }
        
        int nCorrections = tmpPointsRemoved.size();
        
        return (nCorrections > 0);
    }
}
