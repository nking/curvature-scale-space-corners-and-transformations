package algorithms.imageProcessing;

import algorithms.util.PairInt;
import java.util.HashSet;
import java.util.Set;

/**
 * and incomplete class to remove single pixels spur patterns from curves, 
 * with the assumption that the given points are curves.
 * 
 * @author nichole
 */
class SpurRemover {
    
    /**
     * an incomplete set of patterns are applied to the curve to remove spurs
     * that are longer than a pixel and are 1 pixel wide.  this method is not
     * complete.
     * 
     * @param points 
     */
    public void remove(Set<PairInt> curve, int imageWidth, int imageHeight) {
        
        pattern1(curve, imageWidth, imageHeight);
        
        pattern2(curve, imageWidth, imageHeight);
       
        pattern3(curve, imageWidth, imageHeight);
        
        pattern4(curve, imageWidth, imageHeight);
       
        pattern5(curve, imageWidth, imageHeight);
        
        pattern6(curve, imageWidth, imageHeight);
        
        pattern7(curve, imageWidth, imageHeight);
        
        pattern8(curve, imageWidth, imageHeight);
    }

    private void pattern1(Set<PairInt> curve, int imageWidth, int imageHeight) {
        
        /*
        
               #           2
            0  0  #        1
            0 >#  #        0
            0  0  0  #    -1
        
           -1  0  1  2  3
        */
        
        Set<PairInt> ones = new HashSet<PairInt>();
        Set<PairInt> zeroes = new HashSet<PairInt>();
        Set<PairInt> changeToZeroes = new HashSet<PairInt>();
        
        zeroes.add(new PairInt(-1, -1)); zeroes.add(new PairInt(-1, 0)); zeroes.add(new PairInt(-1,  1)); 
        zeroes.add(new PairInt(0, -1)); zeroes.add(new PairInt(0, 1));
        zeroes.add(new PairInt(1, -1));
        
        ones.add(new PairInt(0, 2)); 
        ones.add(new PairInt(1, 0)); ones.add(new PairInt(1, 1));
        ones.add(new PairInt(2, -1));
        
        changeToZeroes.add(new PairInt(0, 0));
        
        int nChanged = replacePatternSwapDirections(curve, 
            imageWidth, imageHeight, zeroes, ones, changeToZeroes);
        
    }
    
    private void pattern2(Set<PairInt> curve, int imageWidth, int imageHeight) {
        
        /*
        
                  #        2
            0  0  #        1
            0 >#  #        0
            0  0  0  #    -1
           -1  0  1  2  3
        */
        
        Set<PairInt> ones = new HashSet<PairInt>();
        Set<PairInt> zeroes = new HashSet<PairInt>();
        Set<PairInt> changeToZeroes = new HashSet<PairInt>();
        
        zeroes.add(new PairInt(-1, -1)); zeroes.add(new PairInt(-1, 0)); zeroes.add(new PairInt(-1,  1)); 
        zeroes.add(new PairInt(0, -1)); zeroes.add(new PairInt(0, 1));
        zeroes.add(new PairInt(1, -1));
        
        ones.add(new PairInt(1, 0)); ones.add(new PairInt(1, 1)); ones.add(new PairInt(1, 2)); 
        ones.add(new PairInt(2, -1));
        
        changeToZeroes.add(new PairInt(0, 0));
        
        int nChanged = replacePatternSwapDirections(curve, 
            imageWidth, imageHeight, zeroes, ones, changeToZeroes);
    }
    
    private void pattern3(Set<PairInt> curve, int imageWidth, int imageHeight) {
        
        /*
                     #      3
            0  0  0  #      2
            0 >#  #  #      1
            0 >#  0         0
            0  0  0        -1
           -1  0  1  2  3
        */
        
        Set<PairInt> ones = new HashSet<PairInt>();
        Set<PairInt> zeroes = new HashSet<PairInt>();
        Set<PairInt> changeToZeroes = new HashSet<PairInt>();
        
        zeroes.add(new PairInt(-1, -1)); zeroes.add(new PairInt(-1, 0)); 
        zeroes.add(new PairInt(-1,  1)); zeroes.add(new PairInt(-1, 2));
        zeroes.add(new PairInt(0, -1)); zeroes.add(new PairInt(0, 2));
        zeroes.add(new PairInt(1, -1)); zeroes.add(new PairInt(1, 0)); zeroes.add(new PairInt(1, 2));
        
        ones.add(new PairInt(0, 1)); 
        ones.add(new PairInt(1, 1)); 
        ones.add(new PairInt(2, 1)); ones.add(new PairInt(2, 2)); ones.add(new PairInt(2, 3)); 
        
        changeToZeroes.add(new PairInt(0, 0));
        changeToZeroes.add(new PairInt(0, 1));
    
        int nChanged = replacePatternSwapDirections(curve, imageWidth, 
            imageHeight, zeroes, ones, changeToZeroes);
    }
    
    private void pattern4(Set<PairInt> curve, int imageWidth, int imageHeight) {
        
        /*
                            3
            0  0  0  0      2
            0 >#  #  0      1
            0  #  0  0      0
            #  0  #        -1
           -1  0  1  2  3
        */
        
        Set<PairInt> ones = new HashSet<PairInt>();
        Set<PairInt> zeroes = new HashSet<PairInt>();
        Set<PairInt> changeToZeroes = new HashSet<PairInt>();
        
        zeroes.add(new PairInt(-1, 0)); zeroes.add(new PairInt(-1,  1)); zeroes.add(new PairInt(-1, 2));
        zeroes.add(new PairInt(0, -1)); zeroes.add(new PairInt(0, 2));
        zeroes.add(new PairInt(1, 0)); zeroes.add(new PairInt(1, 2));
        zeroes.add(new PairInt(2, 0)); zeroes.add(new PairInt(2, 1)); zeroes.add(new PairInt(2, 2));
        
        ones.add(new PairInt(-1, -1));
        ones.add(new PairInt(0, 1)); 
        ones.add(new PairInt(1, -1)); ones.add(new PairInt(1, 1)); 
        
        changeToZeroes.add(new PairInt(0, 1));
        changeToZeroes.add(new PairInt(1, 1));
    
        int nChanged = replacePatternSwapDirections(curve, imageWidth, 
            imageHeight, zeroes, ones, changeToZeroes);
    }
    
    private void pattern5(Set<PairInt> curve, int imageWidth, int imageHeight) {
        
        /*
                            3
            0  0  0  0      2
            0 ># >#  0      1
            #  #  0  0      0
            #  0  #        -1
           -1  0  1  2  3
        */
        
        Set<PairInt> ones = new HashSet<PairInt>();
        Set<PairInt> zeroes = new HashSet<PairInt>();
        Set<PairInt> changeToZeroes = new HashSet<PairInt>();
        
        zeroes.add(new PairInt(-1,  1)); zeroes.add(new PairInt(-1, 2));
        zeroes.add(new PairInt(0, -1)); zeroes.add(new PairInt(0, 2));
        zeroes.add(new PairInt(1, 0)); zeroes.add(new PairInt(1, 2));
        zeroes.add(new PairInt(2, 0)); zeroes.add(new PairInt(2, 1)); zeroes.add(new PairInt(2, 2));
        
        ones.add(new PairInt(-1, -1)); ones.add(new PairInt(-1, 0));
        ones.add(new PairInt(0, 1)); 
        ones.add(new PairInt(1, -1)); ones.add(new PairInt(1, 1)); 
        
        changeToZeroes.add(new PairInt(0, 1));
        changeToZeroes.add(new PairInt(1, 1));
    
        int nChanged = replacePatternSwapDirections(curve, imageWidth, 
            imageHeight, zeroes, ones, changeToZeroes);
    }
    
    private void pattern6(Set<PairInt> curve, int imageWidth, int imageHeight) {
        
        /*
                            3
            0  0  0  0      2
            0 ># >#  0      1
            0  #  0  0      0
               #           -1
           -1  0  1  2  3
        */
        
        Set<PairInt> ones = new HashSet<PairInt>();
        Set<PairInt> zeroes = new HashSet<PairInt>();
        Set<PairInt> changeToZeroes = new HashSet<PairInt>();
        
        zeroes.add(new PairInt(-1,  0)); zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, 2));
        zeroes.add(new PairInt(0, 2));
        zeroes.add(new PairInt(1, 0)); zeroes.add(new PairInt(1, 2));
        zeroes.add(new PairInt(2, 0)); zeroes.add(new PairInt(2, 1)); zeroes.add(new PairInt(2, 2));
        
        ones.add(new PairInt(0, -1)); ones.add(new PairInt(0, 1));
        ones.add(new PairInt(1, 1)); 
        
        changeToZeroes.add(new PairInt(0, 1));
        changeToZeroes.add(new PairInt(1, 1));
    
        int nChanged = replacePatternSwapDirections(curve, imageWidth, 
            imageHeight, zeroes, ones, changeToZeroes);
    }
    
    private void pattern7(Set<PairInt> curve, int imageWidth, int imageHeight) {
        
        /*
                            3
            0  0  0         2
            0 >#  #         1
            0 >#  0  #      0
            0  0  0        -1
           -1  0  1  2  3
        */
        
        Set<PairInt> ones = new HashSet<PairInt>();
        Set<PairInt> zeroes = new HashSet<PairInt>();
        Set<PairInt> changeToZeroes = new HashSet<PairInt>();
        
        zeroes.add(new PairInt(-1,  -1)); zeroes.add(new PairInt(-1,  0)); 
        zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, 2));
        zeroes.add(new PairInt(0, -1)); zeroes.add(new PairInt(0, 2));
        zeroes.add(new PairInt(1, -1)); zeroes.add(new PairInt(1, 0)); zeroes.add(new PairInt(1, 2));
        
        ones.add(new PairInt(0, 1));
        ones.add(new PairInt(1, 1));
        ones.add(new PairInt(2, 0));
        
        changeToZeroes.add(new PairInt(0, 0));
        changeToZeroes.add(new PairInt(0, 1));
    
        int nChanged = replacePatternSwapDirections(curve, imageWidth, 
            imageHeight, zeroes, ones, changeToZeroes);
    }
    
    private void pattern8(Set<PairInt> curve, int imageWidth, int imageHeight) {
        
        /*
                            3
            0  0  0  #      2
            0 >#  #         1
            0 >#  0         0
            0  0  0        -1
           -1  0  1  2  3
        */
        
        Set<PairInt> ones = new HashSet<PairInt>();
        Set<PairInt> zeroes = new HashSet<PairInt>();
        Set<PairInt> changeToZeroes = new HashSet<PairInt>();
        
        zeroes.add(new PairInt(-1,  -1)); zeroes.add(new PairInt(-1,  0)); 
        zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, 2));
        zeroes.add(new PairInt(0, -1)); zeroes.add(new PairInt(0, 2));
        zeroes.add(new PairInt(1, -1)); zeroes.add(new PairInt(1, 0)); zeroes.add(new PairInt(1, 2));
        
        ones.add(new PairInt(0, 1));
        ones.add(new PairInt(1, 1));
        ones.add(new PairInt(2, 2));
        
        changeToZeroes.add(new PairInt(0, 0));
        changeToZeroes.add(new PairInt(0, 1));
    
        int nChanged = replacePatternSwapDirections(curve, imageWidth, 
            imageHeight, zeroes, ones, changeToZeroes);
    }
    
    private int replacePattern(Set<PairInt> points, int imageWidth, 
        int imageHeight, final Set<PairInt> zeroes, final Set<PairInt> ones, 
        final Set<PairInt> changeToZeroes) {
        
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
            
if ((Math.abs(col - 87) < 7) && (Math.abs(row - 89) < 4)) {
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
        int imageWidth, int imageHeight, Set<PairInt> zeroes, 
        Set<PairInt> ones, Set<PairInt> changeToZeroes) {
        
        int nChanged = replacePattern(curve, imageWidth, imageHeight, zeroes, 
            ones, changeToZeroes);
        
        // ----- change the sign of x to handle other direction -----
        reverseXs(zeroes, ones, changeToZeroes);
        
        nChanged += replacePattern(curve, imageWidth, imageHeight, zeroes, 
            ones, changeToZeroes);
        
        // ----- change the sign of y to handle other direction -----
        reverseYs(zeroes, ones, changeToZeroes);
        
        nChanged = replacePattern(curve, imageWidth, imageHeight, zeroes, 
            ones, changeToZeroes);
        
        // ----- change the sign of x to handle other direction -----
        reverseXs(zeroes, ones, changeToZeroes);
        
        nChanged += replacePattern(curve, imageWidth, imageHeight, zeroes, 
            ones, changeToZeroes);
        
        // ----- change the sign of y to handle other direction -----
        reverseYs(zeroes, ones, changeToZeroes);
        
        return nChanged;
    }
    
    protected void reverseXs(final Set<PairInt> zeroes, final Set<PairInt> ones, 
        Set<PairInt> changeToZeroes) {
        
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
    }
    
    protected void reverseYs(final Set<PairInt> zeroes, final Set<PairInt> ones, 
        Set<PairInt> changeToZeroes) {
        
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
        
    }
    
}
