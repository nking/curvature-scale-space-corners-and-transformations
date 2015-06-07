package algorithms.imageProcessing;

import algorithms.util.PairInt;
import java.util.HashSet;
import java.util.Set;
import java.util.logging.Logger;

/**
 * The erosion filter iterates over all pixels in the image, checking whether
 * the current non-null pixel can be nulled without disconnecting any adjacent
 * pixels. 
 * 
 * TODO: the erosion filter doesn't know the overall shape so corners and lines
 * could be improved in the final result by adding templates to the filter
 * to find and handle potential lines and corners in a non-resolution dependent 
 * way.
 * 
 * TODO: It's all boolean rules, so it looks like a clever way of using boolean 
 * logic on more than 1 pixel as the current being tested for nullability,
 * at the same time should reduce common comparisons and make a faster filter.
 * 
 * @author nichole
 */
public class ErosionFilter extends AbstractLineThinner {
    
    protected boolean useLineDrawingMode = false;
    
    protected boolean debug = false;
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    /**
     * for images which are already line drawings, that is images such as
     * maps with only lines, or for block images, use this to avoid a gap filling
     * stage that fills single pixel gaps surrounded by non-zero pixels.  
     * (Else, the filter applies such a gap filling algorithm to help avoid 
     * creating bubbles in thick lines).
     */
    @Override
    public void useLineDrawingMode() {
        useLineDrawingMode = true;
    }

    @Override
    public void setDebug(boolean setToDebug) {
        debug = setToDebug;
    }
    
    @Override
    public void applyFilter(final GreyscaleImage input) {
        
        GreyscaleImage output = input.copyImage();
        
        if (!useLineDrawingMode) {
            prefillGapsSurroundedByNeighbors(output);
        }
        
        int count = 1;
        
        int maxIter = 100;
        int nIter = 0;
        
        GreyscaleImage input2;
        
        while ((count > 0) && (nIter < maxIter)) {
            
            count = 0;
            
            input2 = output.copyImage();
            
            // alternate the directions of approach to 'erode' from both 'sides'
            if ((nIter & 1) == 0) {
                for (int i = 0; i < input2.getWidth(); i++) {
                    for (int j = 0; j < input2.getHeight(); j++) {
                        boolean nulled = process(input2, output, i, j);
                        if (nulled) {
                            count++;
                        }
                    }
                }
            } else {
                for (int i = (input2.getWidth() - 1); i > -1; i--) {
                    for (int j = (input2.getHeight() - 1); j > -1; j--) {
                        boolean nulled = process(input2, output, i, j);
                        if (nulled) {
                            count++;
                        }
                    }
                }
            }                                    
            
            log.fine("nulled " + count + " pixels");
            
            nIter++;
        }
        
        input.resetTo(output);
    }

    /**
     * instead of dilation, just fill in any gaps that have 4 neighbors already
     * filled in.  This should not act upon many pixels.
     * 
     * @param input
     */
    protected void prefillGapsSurroundedByNeighbors(final GreyscaleImage input) {
        
        int w = input.getWidth();
        int h = input.getHeight();
        
        int count = 1;
               
        for (int i = 0; i < w; i++) {

            for (int j = 0; j < h; j++) {

                if (input.getValue(i, j) > 0) {
                    continue;
                }
                       
                if (!isASmallerSearchRegion(i, j, w, h) &&
                    hasImmediateFourNeighbors(input, i, j)) {
                      
                    float avg = input.getValue(i, j + 1) + input.getValue(i + 1, j)
                        + input.getValue(i, j - 1) + input.getValue(i - 1, j);
                    
                    avg /= 4.f;
                    
                    input.setValue(i, j, (int)avg);
                    
                    count++;
                }
            }
        }

        log.fine("filled in  " + count + " pixels");
    }

    /**
     * within the surrounding 8 pixels, return false if the point is the 
     * endpoint of a line with width of 1 pixel, or
     * return true if nulling the center point
     * would not disconnect the othe points or return true if the point.
     * 
     * TODO: compute the runtime complexity of this.
     * 
     * @param input
     * @param row
     * @param col
     * @return 
     */
    private boolean canBeNulled(final GreyscaleImage input, int col, int row) {
        
        /*
        handle edge cases, that is searches constrained to regions smaller than
        the full 3x3 neighborhood region, that is regions of sizes
           2X2 or 2x3 or 3X2
        */
        
        int w = input.getWidth();
        int h = input.getHeight();
       
        boolean isASmallerSearchRegion = isASmallerSearchRegion(
            col, row, w, h);
        
        if (isASmallerSearchRegion) {
            
            boolean doesDisconnect = isASmallerSearchRegionAndDoesDisconnect(
                input, col, row);

            return !doesDisconnect;
        }
        
        /*
        check that this isn't the endpoint of a 1 pixel width line.
        */
        if (isTheEndOfAOnePixelLine(input, col, row)) {
            return false;
        }
        
        /*
        else, try the 8 pixel neighborhood with C being the current (col, row)
        and @ being the pixel tested for disconnection by setting C to 0. 
       
        . is a 0
        ? is an AND between possible '1's.  ("At Least One")
        _ is a pixel that can be ignored

             (1)       (2)     (3)     (4)     (5)     (6)     (7)     (8)
             @ . ?    _ @ _   ? . @   ? . _   ? ? ?   ? ? ?   ? ? ?   _ . ?
             . C ?    . C .   ? C .   ? C @   ? C .   . C .   . C ?   @ C ?
             ? ? ?    ? ? ?   ? ? ?   ? . _   ? . @   _ @ _   @ . ?   _ . ?
        */

        boolean doesDisconnect = doesDisconnect(input, col, row);
        
        return !doesDisconnect;
       
    }
    
    /**
     * this assumes that the invoker has checked that (col, row) has 1 pixel
     * to each side and 1 above and below.
     * @param input
     * @param col
     * @param row
     * @return 
     */
    protected boolean hasImmediateFourNeighbors(final GreyscaleImage input, 
        int col, int row) {
        
        /*   @
           @ C @
             @
        */
        if ((input.getValue(col, row + 1) > 0) &&
            (input.getValue(col + 1, row) > 0) &&
            (input.getValue(col, row - 1) > 0) &&
            (input.getValue(col - 1, row) > 0)) {
            return true;
        }
        
        return false;
    }
   
    /**
     * this assumes that the invoker has checked that (col, row) has 1 pixel
     * to each side and 1 above and below.
    
     * @return 
     */
    protected boolean hasImmediateFourNeighbors(PairInt p, Set<PairInt> points, 
        Set<PairInt> overridePointsAdded, Set<PairInt> overridePointsRemoved) {
        
        /*   @
           @ C @
             @
        */        
        for (int nIdx = 0; nIdx < fourNeighborsX.length; nIdx++) {
            int x = p.getX() + fourNeighborsX[nIdx];
            int y = p.getY() + fourNeighborsY[nIdx];
            PairInt p2 = new PairInt(x, y);
            
            if (overridePointsRemoved.contains(p2)) {
                return false;
            }
            if (!overridePointsAdded.contains(p2) && !points.contains(p2)) {
                return false;
            }
        }
        
        return true;
    }

    /**
     * for the full 8 neighbor region, determine whether nulling the pixel
     * at (col, row) would disconnect the remaining line.  Note that the
     * boolean logic is embedded in the comments.  One should be able to
     * combine the rules for multiple pixel tests to reduce the redundant
     * comparisons for the regions in common.
     * 
     * Note, that row and col are expected to be at least 1 pixel distant
     * from the image borders.
     * 
     * @param input
     * @param col
     * @param row
     * @return 
     */
    protected boolean doesDisconnect(final GreyscaleImage input, int col, 
        int row) {
        
        /*
            6   7  8 
           11 *C* 12
           15  16 17
        */
                
        //6  Ʌ  8  Ʌ  ¬((7)  V  (11  Ʌ  16  Ʌ  12))
        if (
            ((input.getValue(col - 1, row + 1) > 0) 
            && (input.getValue(col + 1, row + 1) > 0)) && 
            !((input.getValue(col, row + 1) > 0) 
            || 
            ((input.getValue(col - 1, row) > 0)
            && (input.getValue(col, row - 1) > 0) 
            && (input.getValue(col + 1, row) > 0)))) {
            return true;
        } else 
            //6  Ʌ  12  Ʌ  ¬((7)  V  (11  Ʌ  16))
            if (
            ((input.getValue(col - 1, row + 1) > 0) 
                && (input.getValue(col + 1, row) > 0)) && 
                !((input.getValue(col, row + 1) > 0)
                || 
                ((input.getValue(col - 1, row) > 0) 
                && (input.getValue(col, row - 1) > 0)))) {
            return true;
        } else 
            //6  Ʌ  17  Ʌ  ¬((7  Ʌ  12) V (11  Ʌ  16))
            if (
            ((input.getValue(col - 1, row + 1) > 0) 
                && (input.getValue(col + 1, row - 1) > 0)) && 
                !(((input.getValue(col, row + 1) > 0)
                && (input.getValue(col + 1, row) > 0)) 
                || ((input.getValue(col - 1, row) > 0) 
                && (input.getValue(col, row - 1) > 0)))) {
            return true;
        } else 
            //6  Ʌ  16  Ʌ  ¬((7  Ʌ  12) V (11))
            if (
            ((input.getValue(col - 1, row + 1) > 0) 
                && (input.getValue(col, row - 1) > 0)) && 
                !(((input.getValue(col, row + 1) > 0)
                && (input.getValue(col + 1, row) > 0)) 
                || (input.getValue(col - 1, row) > 0))) {
            return true;
        } else 
            //6  Ʌ  15  Ʌ  ¬((7  Ʌ  12 Ʌ 16) V (11))
            if (
            ((input.getValue(col - 1, row + 1) > 0) 
                && (input.getValue(col - 1, row - 1) > 0)) && 
                !(((input.getValue(col, row + 1) > 0) 
                && (input.getValue(col + 1, row) > 0) 
                && (input.getValue(col, row - 1) > 0)) 
                || (input.getValue(col - 1, row) > 0))) {
            return true;
        } else 
            //7  Ʌ  17  Ʌ  ¬((12)  V  (11  Ʌ  16))
            if (
            ((input.getValue(col, row + 1) > 0) 
                && (input.getValue(col + 1, row - 1) > 0)) && 
                !(((input.getValue(col - 1, row) > 0)
                && (input.getValue(col, row - 1) > 0)) 
                || (input.getValue(col + 1, row) > 0))) {
            return true;
        } else 
            //7  Ʌ  16  Ʌ  ¬(12  V  11)
            if (
            ((input.getValue(col, row + 1) > 0) 
                && (input.getValue(col, row - 1) > 0)) && 
                !((input.getValue(col + 1, row) > 0)
                || (input.getValue(col - 1, row) > 0))) {
            return true;
        } else 
            //7  Ʌ  15  Ʌ  ¬((11)  V  (12 Ʌ 16))
            if (
            ((input.getValue(col, row + 1) > 0) 
                && (input.getValue(col - 1, row - 1) > 0)) && 
                !(((input.getValue(col + 1, row) > 0) 
                && (input.getValue(col, row - 1) > 0)) 
                || (input.getValue(col - 1, row) > 0))) {
            return true;
        } else 
            //8  Ʌ  17  Ʌ  ¬((12)  V  (7 Ʌ 11  Ʌ  16))
            if (
            ((input.getValue(col + 1, row + 1) > 0) 
                && (input.getValue(col + 1, row - 1) > 0)) && 
                !(((input.getValue(col, row + 1) > 0) 
                && (input.getValue(col - 1, row) > 0) 
                && (input.getValue(col, row - 1) > 0)) 
                || (input.getValue(col + 1, row) > 0))) {
            return true;
        } else 
            //8  Ʌ  16  Ʌ  ¬((12)  V  (7 Ʌ 11))
            if (
            ((input.getValue(col + 1, row + 1) > 0) 
                && (input.getValue(col, row - 1) > 0)) && 
                !(((input.getValue(col, row + 1) > 0) 
                && (input.getValue(col - 1, row) > 0)) 
                || (input.getValue(col + 1, row) > 0))) {
            return true;
        } else 
            //8  Ʌ  15  Ʌ  ¬((12 Ʌ  16)  V  (7 Ʌ 11))
            if (
            ((input.getValue(col + 1, row + 1) > 0) 
                && (input.getValue(col - 1, row - 1) > 0)) && 
                !(((input.getValue(col + 1, row) > 0) 
                && (input.getValue(col, row - 1) > 0)) || 
                ((input.getValue(col, row + 1) > 0) 
                && (input.getValue(col - 1, row) > 0)))) {
            return true;
        } else 
            //8  Ʌ  11  Ʌ  ¬((7)  V  (12 Ʌ 16))
            if (
            ((input.getValue(col + 1, row + 1) > 0) 
                && (input.getValue(col - 1, row) > 0)) && 
                !(((input.getValue(col + 1, row) > 0) 
                && (input.getValue(col, row - 1) > 0)) 
                || (input.getValue(col, row + 1) > 0))) {
            return true;
        } else 
            //12  Ʌ  15  Ʌ  ¬((16)  V  (7 Ʌ 11))
            if (
            ((input.getValue(col + 1, row) > 0) 
                && (input.getValue(col - 1, row - 1) > 0)) && 
                !(((input.getValue(col, row + 1) > 0) 
                && (input.getValue(col - 1, row) > 0)) 
                || (input.getValue(col, row - 1) > 0))) {
            return true;
        } else 
            //12  Ʌ  11  Ʌ  ¬((16)  V  (7))
            if (
            ((input.getValue(col + 1, row) > 0) 
                && (input.getValue(col - 1, row) > 0)) && 
                !((input.getValue(col, row + 1) > 0) 
                || (input.getValue(col, row - 1) > 0))) {
            return true;
        } else 
            //17  Ʌ  15  Ʌ  ¬((16)  V  (7 Ʌ 11 Ʌ 12))
            if (
            ((input.getValue(col + 1, row - 1) > 0) 
                && (input.getValue(col - 1, row - 1) > 0)) && 
                !(((input.getValue(col, row + 1) > 0) 
                && (input.getValue(col - 1, row) > 0) 
                && (input.getValue(col + 1, row) > 0)) 
                || (input.getValue(col, row - 1) > 0))) {
            return true;
        } else 
            //17  Ʌ  11  Ʌ  ¬((16)  V  (7 Ʌ 12))
            if (
            ((input.getValue(col + 1, row - 1) > 0) 
                && (input.getValue(col - 1, row) > 0)) && 
                !(((input.getValue(col, row + 1) > 0) 
                && (input.getValue(col + 1, row) > 0)) 
                || (input.getValue(col, row - 1) > 0))) {
            return true;
        }
        
        return false;
    }
    
    private boolean hasAZeroInNeighbors(final GreyscaleImage input, int col, 
        int row) {
         
        int ks = 3;
        int h = ks >> 1;
        
        // search for a zero in the neighboring 8 pixels
        for (int kCol = 0; kCol < ks; kCol++) {

            int x = kCol - h;
            int imgX = col + x;

            if ((imgX < 0) || (imgX > (input.getWidth() - 1))) {
                continue;
            }

            for (int kRow = 0; kRow < ks; kRow++) {

                int y = kRow - h;
                int imgY = row + y;

                if ((imgY < 0) || (imgY > (input.getHeight() - 1))) {
                    continue;
                }

                if (input.getValue(imgX, imgY) == 0) {
                    return true;
                }
            }
        }
        
        return false;
    }

    /**
     * return true if this endpoint (col, row) is the end of a line that has 
     * a width of 1 pixel.
     * 
     * @param input
     * @param col
     * @param row
     * @return 
     */
    protected boolean isTheEndOfAOnePixelLine(final GreyscaleImage input, 
        int col, int row) {
        
        //TODO: this method could be simplified
        
        /*
          . @ .
          . C .
          . . .
        */
        if ((input.getValue(col, row + 1) > 0) && 
            (input.getValue(col + 1, row + 1) == 0) && 
            (input.getValue(col + 1, row) == 0) && 
            (input.getValue(col + 1, row - 1) == 0) && 
            (input.getValue(col, row - 1) == 0) && 
            (input.getValue(col - 1, row - 1) == 0) && 
            (input.getValue(col - 1, row) == 0) && 
            (input.getValue(col - 1, row + 1) == 0)
        ) {
            /* further check that line thickness is 1 outside of the 8 neighbors
            
            # # #
            . @ .
            . C .
            . . .
            */
            if ((row + 2) > (input.getHeight() - 1)) {
                return true;
            }
            // only 1 of the '#' can be non-zero
            int count = 0;
            if (input.getValue(col - 1, row + 2) > 0) {
                count++;
            }
            if (input.getValue(col, row + 2) > 0) {
                count++;
            }
            if (input.getValue(col + 1, row + 2) > 0) {
                count++;
            }
            if (count > 1) {
                return false;
            }
            
            return true;
        }
        /*
          . . @
          . C .
          . . .
        */
        if ((input.getValue(col + 1, row + 1) > 0) && 
            (input.getValue(col + 1, row) == 0) && 
            (input.getValue(col + 1, row - 1) == 0) &&
            (input.getValue(col, row - 1) == 0) && 
            (input.getValue(col - 1, row - 1) == 0) && 
            (input.getValue(col - 1, row) == 0) && 
            (input.getValue(col - 1, row + 1) == 0) && 
            (input.getValue(col, row + 1) == 0)
            
        ) {
            /* further check that line thickness is 1 outside of the 8 neighbors
            
              # # #
            . . @ #
            . C . #
            . . .
            */
            if ((row + 2) > (input.getHeight() - 1)) {
                /*
                  . . @ #
                  . C . #
                  . . 
                */
                if ((col + 2) > (input.getWidth() - 1)) {
                    return true;
                }
                int count = 0;
                if (input.getValue(col + 2, row + 1) > 0) {
                    count++;
                }
                if (input.getValue(col + 2, row) > 0) {
                    count++;
                }
                if (count > 1) {
                    return false;
                }
                return true;
            } else if ((col + 2) > (input.getWidth() - 1)) {
                /*     # #
                     . . @
                     . C .
                     . . .
                */
                int count = 0;
                if (input.getValue(col, row + 2) > 0) {
                    count++;
                }
                if (input.getValue(col + 1, row + 2) > 0) {
                    count++;
                }
                if (count > 1) {
                    return false;
                }
                return true;
            }
            
            /*
              # # #
            . . @ #
            . C . #
            . . .
            */
            
            // only 1 of the '#' can be non-zero
            int count = 0;
            if (input.getValue(col, row + 2) > 0) {
                count++;
            }
            if (input.getValue(col + 1, row + 2) > 0) {
                count++;
            }
            if (input.getValue(col + 2, row + 2) > 0) {
                count++;
            }
            if (input.getValue(col + 2, row + 1) > 0) {
                count++;
            }
            if (input.getValue(col + 2, row) > 0) {
                count++;
            }
            if (count > 1) {
                return false;
            }
            
            return true;
        }
        /*
          . . .
          . C @
          . . .
        */
        if ((input.getValue(col + 1, row) > 0) && 
            (input.getValue(col + 1, row + 1) == 0) && 
            (input.getValue(col, row + 1) == 0) && 
            (input.getValue(col - 1, row + 1) == 0) && 
            (input.getValue(col - 1, row) == 0) && 
            (input.getValue(col - 1, row - 1) == 0) && 
            (input.getValue(col, row - 1) == 0) && 
            (input.getValue(col + 1, row - 1) == 0)
        ) {
            /*
               . . . #
               . C @ #
               . . . #
            */
            if ((col + 2) > (input.getWidth() - 1)) {
                return true;
            }
            
            // only 1 of the '#' can be non-zero
            int count = 0;
            if (input.getValue(col + 2, row + 1) > 0) {
                count++;
            }
            if (input.getValue(col + 2, row) > 0) {
                count++;
            }
            if (input.getValue(col + 2, row - 1) > 0) {
                count++;
            }
            
            if (count > 1) {
                return false;
            }
            return true;
        }
        
        /*
          . . .
          . C .
          . . @
        */
        if ((input.getValue(col + 1, row - 1) > 0) && 
            (input.getValue(col + 1, row) == 0) && 
            (input.getValue(col + 1, row + 1) == 0) && 
            (input.getValue(col, row + 1) == 0) && 
            (input.getValue(col - 1, row + 1) == 0) && 
            (input.getValue(col - 1, row) == 0) && 
            (input.getValue(col - 1, row - 1) == 0) && 
            (input.getValue(col, row - 1) == 0)
        ) {
            /*
               . . .
               . C . #
               . . @ #
                 # # #
            */
            
            if ((row - 2) < 0) {
                if ((col + 2) > (input.getWidth() - 1)) {
                    return true;
                }
                /*
                     . . .
                     . C . #
                     . . @ #
                 */
                int count = 0;
                if (input.getValue(col + 2, row) > 0) {
                    count++;
                }
                if (input.getValue(col + 2, row - 1) > 0) {
                    count++;
                }

                if (count > 1) {
                    return false;
                }
                return true;
            } else if ((col + 2) > (input.getWidth() - 1)) {
                /*
                     . . .
                     . C .
                     . . @
                       # #
                 */
                int count = 0;
                if (input.getValue(col, row - 2) > 0) {
                    count++;
                }
                if (input.getValue(col + 1, row - 2) > 0) {
                    count++;
                }

                if (count > 1) {
                    return false;
                }
                return true;
            }
            /*
                  . . .
                  . C . #
                  . . @ #
                    # # #
             */
            int count = 0;
            if (input.getValue(col + 2, row) > 0) {
                count++;
            }
            if (input.getValue(col + 2, row - 1) > 0) {
                count++;
            }
            if (input.getValue(col + 2, row - 2) > 0) {
                count++;
            }
            if (input.getValue(col + 1, row - 2) > 0) {
                count++;
            }
            if (input.getValue(col, row - 2) > 0) {
                count++;
            }
            if (count > 1) {
                return false;
            }
            return true;
        }
        
        /*
          . . .
          . C .
          . @ .
        */
        if ((input.getValue(col, row - 1) > 0) && 
            (input.getValue(col + 1, row - 1) == 0) && 
            (input.getValue(col + 1, row) == 0) && 
            (input.getValue(col + 1, row + 1) == 0) && 
            (input.getValue(col, row + 1) == 0) && 
            (input.getValue(col - 1, row + 1) == 0) && 
            (input.getValue(col - 1, row) == 0) && 
            (input.getValue(col - 1, row - 1) == 0)
        ) {
            /*
               . . . 
               . C . 
               . @ . 
               # # #
            */
            if ((row - 2) < 0) {
                return true;
            }
            
            // only 1 of the '#' can be non-zero
            int count = 0;
            if (input.getValue(col - 1, row - 2) > 0) {
                count++;
            }
            if (input.getValue(col, row - 2) > 0) {
                count++;
            }
            if (input.getValue(col + 1, row - 2) > 0) {
                count++;
            }
            
            if (count > 1) {
                return false;
            }
            return true;
        }
        
        /*
          . . .
          . C .
          @ . .
        */
        if ((input.getValue(col - 1, row - 1) > 0) && 
            (input.getValue(col, row - 1) == 0) && 
            (input.getValue(col + 1, row - 1) == 0) && 
            (input.getValue(col + 1, row) == 0) && 
            (input.getValue(col + 1, row + 1) == 0) && 
            (input.getValue(col, row + 1) == 0) && 
            (input.getValue(col - 1, row + 1) == 0) && 
            (input.getValue(col - 1, row) == 0)
        ) {
            /*
                 . . .
               # . C .
               # @ . .
               # # # 
            */
            if ((col - 2) < 0) {
                if ((row - 2) < 0) {
                    return true;
                }
                /*
                      . . .
                      . C .
                      @ . .
                      # # 
                */
                // only 1 of the '#' can be non-zero
                int count = 0;
                if (input.getValue(col - 1, row - 2) > 0) {
                    count++;
                }
                if (input.getValue(col, row - 2) > 0) {
                    count++;
                }
                if (count > 1) {
                    return false;
                }
                return true;
            } else if ((row - 2) < 0) {
                /*
                       . . .
                     # . C .
                     # @ . .
                 */
                int count = 0;
                if (input.getValue(col - 2, row) > 0) {
                    count++;
                }
                if (input.getValue(col - 2, row - 1) > 0) {
                    count++;
                }
                if (count > 1) {
                    return false;
                }
                return true;
            }
            /*
                 . . .
               # . C .
               # @ . .
               # # # 
            */
            int count = 0;
            if (input.getValue(col - 2, row) > 0) {
                count++;
            }
            if (input.getValue(col - 2, row - 1) > 0) {
                count++;
            }
            if (input.getValue(col - 2, row - 2) > 0) {
                count++;
            }
            if (input.getValue(col - 1, row - 2) > 0) {
                count++;
            }
            if (input.getValue(col, row - 2) > 0) {
                count++;
            }

            if (count > 1) {
                return false;
            }
            return true;
        }
        
        /*
          . . .
          @ C .
          . . .
        */
        if ((input.getValue(col - 1, row) > 0) && 
            (input.getValue(col - 1, row + 1) == 0) && 
            (input.getValue(col, row + 1) == 0) && 
            (input.getValue(col + 1, row + 1) == 0) && 
            (input.getValue(col + 1, row) == 0) && 
            (input.getValue(col + 1, row - 1) == 0) && 
            (input.getValue(col, row - 1) == 0) && 
            (input.getValue(col - 1, row - 1) == 0)
        ) {
            /*
                # . . .
                # @ C .
                # . . .
            */
            if ((col - 2) < 0) {
                return true;
            }
            
            // only 1 of the '#' can be non-zero
            int count = 0;
            if (input.getValue(col - 2, row - 1) > 0) {
                count++;
            }
            if (input.getValue(col - 2, row) > 0) {
                count++;
            }
            if (input.getValue(col - 2, row + 1) > 0) {
                count++;
            }
            
            if (count > 1) {
                return false;
            }
            return true;
        }
        
        /*
          @ . .
          . C .
          . . .
        */
        if ((input.getValue(col - 1, row + 1) > 0) && 
            (input.getValue(col, row + 1) == 0) && 
            (input.getValue(col + 1, row + 1) == 0) && 
            (input.getValue(col + 1, row) == 0) && 
            (input.getValue(col + 1, row - 1) == 0) && 
            (input.getValue(col, row - 1) == 0) && 
            (input.getValue(col - 1, row - 1) == 0) && 
            (input.getValue(col - 1, row) == 0)
        ) {
            if ((col - 2) < 0) {
                if ((row + 2) > (input.getHeight() - 1)) {
                    return true;
                }
                /*
                     # #
                     @ . .
                     . C .
                     . . .
                */
                // only 1 of the '#' can be non-zero
                int count = 0;
                if (input.getValue(col - 1, row + 2) > 0) {
                    count++;
                }
                if (input.getValue(col, row + 2) > 0) {
                    count++;
                }

                if (count > 1) {
                    return false;
                }
                return true;
            } else if ((row + 2) > (input.getHeight() - 1)) {
                /*
                    # @ . .
                    # . C .
                    . . .
                 */
                int count = 0;
                if (input.getValue(col - 2, row) > 0) {
                    count++;
                }
                if (input.getValue(col - 2, row + 1) > 0) {
                    count++;
                }

                if (count > 1) {
                    return false;
                }
                return true;
            }
            /*
               # # #
               # @ . .
               # . C .
                 . . .
            */
            int count = 0;
            if (input.getValue(col - 2, row) > 0) {
                count++;
            }
            if (input.getValue(col - 2, row + 1) > 0) {
                count++;
            }
            if (input.getValue(col - 2, row + 2) > 0) {
                count++;
            }
            if (input.getValue(col - 1, row + 2) > 0) {
                count++;
            }
            if (input.getValue(col, row + 2) > 0) {
                count++;
            }

            if (count > 1) {
                return false;
            }
            return true;
        }
        
        return false;
    }

    /**
     * return true if pixel at (i, j) is nullable, that is does not disconnect
     * other currently connected points.  the method does set the pixel value
     * to 0 if it can and return true if so.
     * @param input
     * @param output
     * @param i
     * @param j
     * @return 
     */
    boolean process(final GreyscaleImage input, final GreyscaleImage 
        output, int i, int j) {
        
        if (output.getValue(i, j) == 0) {
            return false;
        }
        
        int w = input.getWidth();
        int h = input.getHeight();

        if (!isASmallerSearchRegion(i, j, w, h)
            && hasImmediateFourNeighbors(output, i, j)) {

            return false;
        }

        // choosing input here helps to thin a line from one direction
        // at a time rather than continue thinning from the same
        // direction
        if (hasAZeroInNeighbors(input, i, j)) {

            //if (!doesDisconnectOrIsEndPoint(output, i, j)) {
            if (canBeNulled(output, i, j)) {
                
                output.setValue(i, j, 0);

                return true;
            }
        }
        
        return false;
    }

    /**
     * return true if pixel at (i, j) is nullable, that is does not disconnect
     * other currently connected points. note that if point is not present
     * in points (or is in overridePointsRemoved) then it returns an immediate
     * false.  the method does set the pixel value
     * to 0 if it can and return true if so.
     * 
     * @return 
     */
    boolean process(PairInt p, Set<PairInt> points, Set<PairInt> 
        overridePointsAdded, Set<PairInt> overridePointsRemoved, 
        int imageWidth, int imageHeight) {
                
        // if point does not exist already in point set, return false
        if (overridePointsRemoved.contains(p)) {
            return false;
        } else if (!points.contains(p) && !overridePointsAdded.contains(p)) {
            return false;
        }
        
        if (!isASmallerSearchRegion(p.getX(), p.getY(), imageWidth, imageHeight)
            && hasImmediateFourNeighbors(p, points, overridePointsAdded, 
            overridePointsRemoved)) {

            return false;
        }

        // choosing input here helps to thin a line from one direction
        // at a time rather than continue thinning from the same
        // direction
        if (hasAZeroInNeighbors(p, points, overridePointsAdded, 
            overridePointsRemoved, imageWidth, imageHeight)) {

            if (canBeNulled(p, points, overridePointsAdded, 
                overridePointsRemoved, imageWidth, imageHeight)) {
                
                overridePointsRemoved.add(p);

                return true;
            }
        }
        
        return false;
    }
    
    private boolean hasAZeroInNeighbors(PairInt p, Set<PairInt> points, 
        Set<PairInt> overridePointsAdded, Set<PairInt> overridePointsRemoved,
        int imageWidth, int imageHeight) {
        
        for (int nIdx = 0; nIdx < fourNeighborsX.length; nIdx++) {
            int x = p.getX() + fourNeighborsX[nIdx];
            int y = p.getY() + fourNeighborsY[nIdx];
            
            if ((x < 0) || (x > (imageWidth - 1)) || (y < 0) || 
                (y > (imageHeight - 1))) {
                continue;
            }

            PairInt p2 = new PairInt(x, y);
            
            if (overridePointsRemoved.contains(p2)) {
                return true;
            }
            if (!overridePointsAdded.contains(p2) && !points.contains(p2)) {
                return true;
            }
        }
        
        return false;
                
    }

    private boolean canBeNulled(PairInt p, Set<PairInt> points, 
        Set<PairInt> overridePointsAdded, Set<PairInt> overridePointsRemoved, 
        int imageWidth, int imageHeight) {
        
        /*
        handle edge cases, that is searches constrained to regions smaller than
        the full 3x3 neighborhood region, that is regions of sizes
           2X2 or 2x3 or 3X2
        */
        
        boolean isASmallerSearchRegion = isASmallerSearchRegion(
            p.getX(), p.getY(), imageWidth, imageHeight);
        
        if (isASmallerSearchRegion) {
            
            boolean doesDisconnect = isASmallerSearchRegionAndDoesDisconnect(
                p, points, overridePointsAdded, overridePointsRemoved,
                imageWidth, imageHeight);

            return !doesDisconnect;
        }
        
        /*
        check that this isn't the endpoint of a 1 pixel width line.
        */
        if (isTheEndOfAOnePixelLine(p, points, overridePointsAdded, 
            overridePointsRemoved, imageWidth, imageHeight)) {
            return false;
        }
        
        /*
        else, try the 8 pixel neighborhood with C being the current (col, row)
        and @ being the pixel tested for disconnection by setting C to 0. 
       
        . is a 0
        ? is an AND between possible '1's.  ("At Least One")
        _ is a pixel that can be ignored

             (1)       (2)     (3)     (4)     (5)     (6)     (7)     (8)
             @ . ?    _ @ _   ? . @   ? . _   ? ? ?   ? ? ?   ? ? ?   _ . ?
             . C ?    . C .   ? C .   ? C @   ? C .   . C .   . C ?   @ C ?
             ? ? ?    ? ? ?   ? ? ?   ? . _   ? . @   _ @ _   @ . ?   _ . ?
        */

        boolean doesDisconnect = doesDisconnect(p, points, overridePointsAdded, 
            overridePointsRemoved, imageWidth, imageHeight);
        
        return !doesDisconnect;
    }

    private boolean isTheEndOfAOnePixelLine(PairInt p, Set<PairInt> points, 
        Set<PairInt> overridePointsAdded, Set<PairInt> overridePointsRemoved, 
        int imageWidth, int imageHeight) {
        
        //TODO: this method could be simplified
        
        int col = p.getX();
        int row = p.getY();
        
        /*
        coordinates of the 8 neighbors as already created PairInts without 
        bound checks.
        indexes are found as +1 of the difference relative to center,
        for example, a point for (col-1, row-1) is found as neighborCoords[0][0]
        */
        PairInt[][] neighborCoords = createCoordinatePointsForEightNeighbors(
            col, row);
        
        /*
          . @ .
          . C .
          . . .
        */
        if (
            (!overridePointsRemoved.contains(neighborCoords[1][2]) &&
            (overridePointsAdded.contains(neighborCoords[1][2]) ||
            points.contains(neighborCoords[1][2])))
            && 
            (overridePointsRemoved.contains(neighborCoords[2][2]) ||
            (!overridePointsAdded.contains(neighborCoords[2][2]) &&
            !points.contains(neighborCoords[2][2])))
            && 
            (overridePointsRemoved.contains(neighborCoords[2][1]) ||
            (!overridePointsAdded.contains(neighborCoords[2][1]) &&
            !points.contains(neighborCoords[2][1])))
            && 
            (overridePointsRemoved.contains(neighborCoords[2][0]) ||
            (!overridePointsAdded.contains(neighborCoords[2][0]) &&
            !points.contains(neighborCoords[2][0])))
            && 
            (overridePointsRemoved.contains(neighborCoords[1][0]) ||
            (!overridePointsAdded.contains(neighborCoords[1][0]) &&
            !points.contains(neighborCoords[1][0])))
            && 
            (overridePointsRemoved.contains(neighborCoords[0][0]) ||
            (!overridePointsAdded.contains(neighborCoords[0][0]) &&
            !points.contains(neighborCoords[0][0])))
            && 
            (overridePointsRemoved.contains(neighborCoords[0][1]) ||
            (!overridePointsAdded.contains(neighborCoords[0][1]) &&
            !points.contains(neighborCoords[0][1])))
            && 
            (overridePointsRemoved.contains(neighborCoords[0][2]) ||
            (!overridePointsAdded.contains(neighborCoords[0][2]) &&
            !points.contains(neighborCoords[0][2])))
        ) {
            /* further check that line thickness is 1 outside of the 8 neighbors
            
            # # #
            . @ .
            . C .
            . . .
            */
            if ((row + 2) > (imageHeight - 1)) {
                return true;
            }
            // only 1 of the '#' can be non-zero
            int count = 0;
            PairInt p2 = new PairInt(col - 1, row + 2);
            if (
                (!overridePointsRemoved.contains(p2) &&
                (overridePointsAdded.contains(p2) ||
                points.contains(p2)))
                ) {
                count++;
            }
            p2 = new PairInt(col, row + 2);
            if (
                (!overridePointsRemoved.contains(p2) &&
                (overridePointsAdded.contains(p2) ||
                points.contains(p2)))
                ) {
                count++;
            }
            p2 = new PairInt(col + 1, row + 2);
            if (
                (!overridePointsRemoved.contains(p2) &&
                (overridePointsAdded.contains(p2) ||
                points.contains(p2)))
                ) {
                count++;
            }
            if (count > 1) {
                return false;
            }
            
            return true;
        }
        /*
          . . @
          . C .
          . . .
        */
        if (
            (!overridePointsRemoved.contains(neighborCoords[2][2]) &&
            (overridePointsAdded.contains(neighborCoords[2][2]) ||
            points.contains(neighborCoords[2][2])))
            && 
            (overridePointsRemoved.contains(neighborCoords[2][1]) ||
            (!overridePointsAdded.contains(neighborCoords[2][1]) &&
            !points.contains(neighborCoords[2][1])))
            && 
            (overridePointsRemoved.contains(neighborCoords[2][0]) ||
            (!overridePointsAdded.contains(neighborCoords[2][0]) &&
            !points.contains(neighborCoords[2][0])))
            &&
            (overridePointsRemoved.contains(neighborCoords[1][0]) ||
            (!overridePointsAdded.contains(neighborCoords[1][0]) &&
            !points.contains(neighborCoords[1][0])))
            && 
            (overridePointsRemoved.contains(neighborCoords[0][0]) ||
            (!overridePointsAdded.contains(neighborCoords[0][0]) &&
            !points.contains(neighborCoords[0][0])))
            && 
            (overridePointsRemoved.contains(neighborCoords[0][1]) ||
            (!overridePointsAdded.contains(neighborCoords[0][1]) &&
            !points.contains(neighborCoords[0][1])))
            && 
            (overridePointsRemoved.contains(neighborCoords[0][2]) ||
            (!overridePointsAdded.contains(neighborCoords[0][2]) &&
            !points.contains(neighborCoords[0][2])))
            && 
            (overridePointsRemoved.contains(neighborCoords[1][2]) ||
            (!overridePointsAdded.contains(neighborCoords[1][2]) &&
            !points.contains(neighborCoords[1][2])))
            
        ) {
            /* further check that line thickness is 1 outside of the 8 neighbors
            
              # # #
            . . @ #
            . C . #
            . . .
            */
            if ((row + 2) > (imageHeight - 1)) {
                /*
                  . . @ #
                  . C . #
                  . . 
                */
                if ((col + 2) > (imageWidth - 1)) {
                    return true;
                }
                int count = 0;
                PairInt p2 = new PairInt(col + 2, row + 1);
                if (
                    (!overridePointsRemoved.contains(p2) &&
                    (overridePointsAdded.contains(p2) ||
                    points.contains(p2)))
                    ) {
                    count++;
                }
                p2 = new PairInt(col + 2, row);
                if (
                    (!overridePointsRemoved.contains(p2) &&
                    (overridePointsAdded.contains(p2) ||
                    points.contains(p2)))
                    ) {
                    count++;
                }
                if (count > 1) {
                    return false;
                }
                return true;
            } else if ((col + 2) > (imageWidth - 1)) {
                /*     # #
                     . . @
                     . C .
                     . . .
                */
                int count = 0;
                PairInt p2 = new PairInt(col, row + 2);
                if (
                    (!overridePointsRemoved.contains(p2) &&
                    (overridePointsAdded.contains(p2) ||
                    points.contains(p2)))
                    ) {
                    count++;
                }
                p2 = new PairInt(col + 1, row + 2);
                if (
                    (!overridePointsRemoved.contains(p2) &&
                    (overridePointsAdded.contains(p2) ||
                    points.contains(p2)))
                    ) {
                    count++;
                }
                if (count > 1) {
                    return false;
                }
                return true;
            }
            
            /*
              # # #
            . . @ #
            . C . #
            . . .
            */
            
            // only 1 of the '#' can be non-zero
            int count = 0;
            PairInt p2 = new PairInt(col, row + 2);
            if (
                (!overridePointsRemoved.contains(p2) &&
                (overridePointsAdded.contains(p2) ||
                points.contains(p2)))
                ) {
                count++;
            }
            p2 = new PairInt(col + 1, row + 2);
            if (
                (!overridePointsRemoved.contains(p2) &&
                (overridePointsAdded.contains(p2) ||
                points.contains(p2)))
                ) {
                count++;
            }
            p2 = new PairInt(col + 2, row + 2);
            if (
                (!overridePointsRemoved.contains(p2) &&
                (overridePointsAdded.contains(p2) ||
                points.contains(p2)))
                ) {
                count++;
            }
            p2 = new PairInt(col + 2, row + 1);
            if (
                (!overridePointsRemoved.contains(p2) &&
                (overridePointsAdded.contains(p2) ||
                points.contains(p2)))
                ) {
                count++;
            }
            p2 = new PairInt(col + 2, row);
            if (
                (!overridePointsRemoved.contains(p2) &&
                (overridePointsAdded.contains(p2) ||
                points.contains(p2)))
                ) {
                count++;
            }
            if (count > 1) {
                return false;
            }
            
            return true;
        }
        /*
          . . .
          . C @
          . . .
        */
        if (
            (!overridePointsRemoved.contains(neighborCoords[2][1]) &&
            (overridePointsAdded.contains(neighborCoords[2][1]) ||
            points.contains(neighborCoords[2][1])))
            && 
            (overridePointsRemoved.contains(neighborCoords[2][2]) ||
            (!overridePointsAdded.contains(neighborCoords[2][2]) &&
            !points.contains(neighborCoords[2][2])))
            && 
            (overridePointsRemoved.contains(neighborCoords[1][2]) ||
            (!overridePointsAdded.contains(neighborCoords[1][2]) &&
            !points.contains(neighborCoords[1][2])))
            && 
            (overridePointsRemoved.contains(neighborCoords[0][2]) ||
            (!overridePointsAdded.contains(neighborCoords[0][2]) &&
            !points.contains(neighborCoords[0][2])))
            && 
            (overridePointsRemoved.contains(neighborCoords[0][1]) ||
            (!overridePointsAdded.contains(neighborCoords[0][1]) &&
            !points.contains(neighborCoords[0][1])))
            && 
            (overridePointsRemoved.contains(neighborCoords[0][0]) ||
            (!overridePointsAdded.contains(neighborCoords[0][0]) &&
            !points.contains(neighborCoords[0][0])))
            && 
            (overridePointsRemoved.contains(neighborCoords[1][0]) ||
            (!overridePointsAdded.contains(neighborCoords[1][0]) &&
            !points.contains(neighborCoords[1][0])))
            && 
            (overridePointsRemoved.contains(neighborCoords[2][0]) ||
            (!overridePointsAdded.contains(neighborCoords[2][0]) &&
            !points.contains(neighborCoords[2][0])))
        ) {
            /*
               . . . #
               . C @ #
               . . . #
            */
            if ((col + 2) > (imageWidth - 1)) {
                return true;
            }
            
            // only 1 of the '#' can be non-zero
            int count = 0;
            PairInt p2 = new PairInt(col + 2, row + 1);
            if (
                (!overridePointsRemoved.contains(p2) &&
                (overridePointsAdded.contains(p2) ||
                points.contains(p2)))
                ) {
                count++;
            }
            p2 = new PairInt(col + 2, row);
            if (
                (!overridePointsRemoved.contains(p2) &&
                (overridePointsAdded.contains(p2) ||
                points.contains(p2)))
                ) {
                count++;
            }
            p2 = new PairInt(col + 2, row - 1);
            if (
                (!overridePointsRemoved.contains(p2) &&
                (overridePointsAdded.contains(p2) ||
                points.contains(p2)))
                ) {
                count++;
            }
            
            if (count > 1) {
                return false;
            }
            return true;
        }
        
        /*
          . . .
          . C .
          . . @
        */
        if (
            (!overridePointsRemoved.contains(neighborCoords[2][0]) &&
            (overridePointsAdded.contains(neighborCoords[2][0]) ||
            points.contains(neighborCoords[2][0])))
            && 
            (overridePointsRemoved.contains(neighborCoords[2][1]) ||
            (!overridePointsAdded.contains(neighborCoords[2][1]) &&
            !points.contains(neighborCoords[2][1])))
            && 
            (overridePointsRemoved.contains(neighborCoords[2][2]) ||
            (!overridePointsAdded.contains(neighborCoords[2][2]) &&
            !points.contains(neighborCoords[2][2])))
            && 
            (overridePointsRemoved.contains(neighborCoords[1][2]) ||
            (!overridePointsAdded.contains(neighborCoords[1][2]) &&
            !points.contains(neighborCoords[1][2])))
            && 
            (overridePointsRemoved.contains(neighborCoords[0][2]) ||
            (!overridePointsAdded.contains(neighborCoords[0][2]) &&
            !points.contains(neighborCoords[0][2])))
            && 
            (overridePointsRemoved.contains(neighborCoords[0][1]) ||
            (!overridePointsAdded.contains(neighborCoords[0][1]) &&
            !points.contains(neighborCoords[0][1])))
            && 
            (overridePointsRemoved.contains(neighborCoords[0][0]) ||
            (!overridePointsAdded.contains(neighborCoords[0][0]) &&
            !points.contains(neighborCoords[0][0])))
            && 
            (overridePointsRemoved.contains(neighborCoords[1][0]) ||
            (!overridePointsAdded.contains(neighborCoords[1][0]) &&
            !points.contains(neighborCoords[1][0])))
        ) {
            /*
               . . .
               . C . #
               . . @ #
                 # # #
            */
            
            if ((row - 2) < 0) {
                if ((col + 2) > (imageWidth - 1)) {
                    return true;
                }
                /*
                     . . .
                     . C . #
                     . . @ #
                 */
                int count = 0;
                PairInt p2 = new PairInt(col + 2, row);
                if (
                    (!overridePointsRemoved.contains(p2) &&
                    (overridePointsAdded.contains(p2) ||
                    points.contains(p2)))
                    ) {
                    count++;
                }
                p2 = new PairInt(col + 2, row - 1);
                if (
                    (!overridePointsRemoved.contains(p2) &&
                    (overridePointsAdded.contains(p2) ||
                    points.contains(p2)))
                    ) {
                    count++;
                }

                if (count > 1) {
                    return false;
                }
                return true;
            } else if ((col + 2) > (imageWidth - 1)) {
                /*
                     . . .
                     . C .
                     . . @
                       # #
                 */
                int count = 0;
                PairInt p2 = new PairInt(col, row - 2);
                if (
                    (!overridePointsRemoved.contains(p2) &&
                    (overridePointsAdded.contains(p2) ||
                    points.contains(p2)))
                    ) {
                    count++;
                }
                p2 = new PairInt(col + 1, row - 2);
                if (
                    (!overridePointsRemoved.contains(p2) &&
                    (overridePointsAdded.contains(p2) ||
                    points.contains(p2)))
                    ) {
                    count++;
                }

                if (count > 1) {
                    return false;
                }
                return true;
            }
            /*
                  . . .
                  . C . #
                  . . @ #
                    # # #
             */
            int count = 0;
            PairInt p2 = new PairInt(col + 2, row);
            if (
                (!overridePointsRemoved.contains(p2) &&
                (overridePointsAdded.contains(p2) ||
                points.contains(p2)))
                ) {
                count++;
            }
            p2 = new PairInt(col + 2, row - 1);
            if (
                (!overridePointsRemoved.contains(p2) &&
                (overridePointsAdded.contains(p2) ||
                points.contains(p2)))
                ) {
                count++;
            }
            p2 = new PairInt(col + 2, row - 2);
            if (
                (!overridePointsRemoved.contains(p2) &&
                (overridePointsAdded.contains(p2) ||
                points.contains(p2)))
                ) {
                count++;
            }
            p2 = new PairInt(col + 1, row - 2);
            if (
                (!overridePointsRemoved.contains(p2) &&
                (overridePointsAdded.contains(p2) ||
                points.contains(p2)))
                ) {
                count++;
            }
            p2 = new PairInt(col, row - 2);
            if (
                (!overridePointsRemoved.contains(p2) &&
                (overridePointsAdded.contains(p2) ||
                points.contains(p2)))
                ) {
                count++;
            }
            if (count > 1) {
                return false;
            }
            return true;
        }
        
        /*
          . . .
          . C .
          . @ .
        */
        if (
            (!overridePointsRemoved.contains(neighborCoords[1][0]) &&
            (overridePointsAdded.contains(neighborCoords[1][0]) ||
            points.contains(neighborCoords[1][0])))
            && 
            (overridePointsRemoved.contains(neighborCoords[2][0]) ||
            (!overridePointsAdded.contains(neighborCoords[2][0]) &&
            !points.contains(neighborCoords[2][0])))
            && 
            (overridePointsRemoved.contains(neighborCoords[2][1]) ||
            (!overridePointsAdded.contains(neighborCoords[2][1]) &&
            !points.contains(neighborCoords[2][1])))
            && 
            (overridePointsRemoved.contains(neighborCoords[2][2]) ||
            (!overridePointsAdded.contains(neighborCoords[2][2]) &&
            !points.contains(neighborCoords[2][2])))
            && 
            (overridePointsRemoved.contains(neighborCoords[1][2]) ||
            (!overridePointsAdded.contains(neighborCoords[1][2]) &&
            !points.contains(neighborCoords[1][2])))
            && 
            (overridePointsRemoved.contains(neighborCoords[0][2]) ||
            (!overridePointsAdded.contains(neighborCoords[0][2]) &&
            !points.contains(neighborCoords[0][2])))
            && 
            (overridePointsRemoved.contains(neighborCoords[0][1]) ||
            (!overridePointsAdded.contains(neighborCoords[0][1]) &&
            !points.contains(neighborCoords[0][1])))
            && 
            (overridePointsRemoved.contains(neighborCoords[0][0]) ||
            (!overridePointsAdded.contains(neighborCoords[0][0]) &&
            !points.contains(neighborCoords[0][0])))
        ) {
            /*
               . . . 
               . C . 
               . @ . 
               # # #
            */
            if ((row - 2) < 0) {
                return true;
            }
            
            // only 1 of the '#' can be non-zero
            int count = 0;
            PairInt p2 = new PairInt(col - 1, row - 2);
            if (
                (!overridePointsRemoved.contains(p2) &&
                (overridePointsAdded.contains(p2) ||
                points.contains(p2)))
                ) {
                count++;
            }
            p2 = new PairInt(col, row - 2);
            if (
                (!overridePointsRemoved.contains(p2) &&
                (overridePointsAdded.contains(p2) ||
                points.contains(p2)))
                ) {
                count++;
            }
            p2 = new PairInt(col + 1, row - 2);
            if (
                (!overridePointsRemoved.contains(p2) &&
                (overridePointsAdded.contains(p2) ||
                points.contains(p2)))
                ) {
                count++;
            }
            
            if (count > 1) {
                return false;
            }
            return true;
        }
        
        /*
          . . .
          . C .
          @ . .
        */
        if (
            (!overridePointsRemoved.contains(neighborCoords[0][0]) &&
            (overridePointsAdded.contains(neighborCoords[0][0]) ||
            points.contains(neighborCoords[0][0])))
            && 
            (overridePointsRemoved.contains(neighborCoords[1][0]) ||
            (!overridePointsAdded.contains(neighborCoords[1][0]) &&
            !points.contains(neighborCoords[1][0])))
            && 
            (overridePointsRemoved.contains(neighborCoords[2][0]) ||
            (!overridePointsAdded.contains(neighborCoords[2][0]) &&
            !points.contains(neighborCoords[2][0])))
            && 
            (overridePointsRemoved.contains(neighborCoords[2][1]) ||
            (!overridePointsAdded.contains(neighborCoords[2][1]) &&
            !points.contains(neighborCoords[2][1])))
            && 
            (overridePointsRemoved.contains(neighborCoords[2][2]) ||
            (!overridePointsAdded.contains(neighborCoords[2][2]) &&
            !points.contains(neighborCoords[2][2])))
            && 
            (overridePointsRemoved.contains(neighborCoords[1][2]) ||
            (!overridePointsAdded.contains(neighborCoords[1][2]) &&
            !points.contains(neighborCoords[1][2])))
            && 
            (overridePointsRemoved.contains(neighborCoords[0][2]) ||
            (!overridePointsAdded.contains(neighborCoords[0][2]) &&
            !points.contains(neighborCoords[0][2])))
            && 
            (overridePointsRemoved.contains(neighborCoords[0][1]) ||
            (!overridePointsAdded.contains(neighborCoords[0][1]) &&
            !points.contains(neighborCoords[0][1])))
        ) {
            /*
                 . . .
               # . C .
               # @ . .
               # # # 
            */
            if ((col - 2) < 0) {
                if ((row - 2) < 0) {
                    return true;
                }
                /*
                      . . .
                      . C .
                      @ . .
                      # # 
                */
                // only 1 of the '#' can be non-zero
                int count = 0;
                PairInt p2 = new PairInt(col - 1, row - 2);
                if (
                    (!overridePointsRemoved.contains(p2) &&
                    (overridePointsAdded.contains(p2) ||
                    points.contains(p2)))
                    ) {
                    count++;
                }
                p2 = new PairInt(col, row - 2);
                if (
                    (!overridePointsRemoved.contains(p2) &&
                    (overridePointsAdded.contains(p2) ||
                    points.contains(p2)))
                    ) {
                    count++;
                }
                if (count > 1) {
                    return false;
                }
                return true;
            } else if ((row - 2) < 0) {
                /*
                       . . .
                     # . C .
                     # @ . .
                 */
                int count = 0;
                PairInt p2 = new PairInt(col - 2, row);
                if (
                    (!overridePointsRemoved.contains(p2) &&
                    (overridePointsAdded.contains(p2) ||
                    points.contains(p2)))
                    ) {
                    count++;
                }
                p2 = new PairInt(col - 2, row - 1);
                if (
                    (!overridePointsRemoved.contains(p2) &&
                    (overridePointsAdded.contains(p2) ||
                    points.contains(p2)))
                    ) {
                    count++;
                }
                if (count > 1) {
                    return false;
                }
                return true;
            }
            /*
                 . . .
               # . C .
               # @ . .
               # # # 
            */
            int count = 0;
            PairInt p2 = new PairInt(col - 2, row);
            if (
                (!overridePointsRemoved.contains(p2) &&
                (overridePointsAdded.contains(p2) ||
                points.contains(p2)))
                ) {
                count++;
            }
            p2 = new PairInt(col - 2, row - 1);
            if (
                (!overridePointsRemoved.contains(p2) &&
                (overridePointsAdded.contains(p2) ||
                points.contains(p2)))
                ) {
                count++;
            }
            p2 = new PairInt(col - 2, row - 2);
            if (
                (!overridePointsRemoved.contains(p2) &&
                (overridePointsAdded.contains(p2) ||
                points.contains(p2)))
                ) {
                count++;
            }
            p2 = new PairInt(col - 1, row - 2);
            if (
                (!overridePointsRemoved.contains(p2) &&
                (overridePointsAdded.contains(p2) ||
                points.contains(p2)))
                ) {
                count++;
            }
            p2 = new PairInt(col, row - 2);
            if (
                (!overridePointsRemoved.contains(p2) &&
                (overridePointsAdded.contains(p2) ||
                points.contains(p2)))
                ) {
                count++;
            }

            if (count > 1) {
                return false;
            }
            return true;
        }
        
        /*
          . . .
          @ C .
          . . .
        */
        if (
            (!overridePointsRemoved.contains(neighborCoords[0][1]) &&
            (overridePointsAdded.contains(neighborCoords[0][1]) ||
            points.contains(neighborCoords[0][1])))
            && 
            (overridePointsRemoved.contains(neighborCoords[0][2]) ||
            (!overridePointsAdded.contains(neighborCoords[0][2]) &&
            !points.contains(neighborCoords[0][2])))
            && 
            (overridePointsRemoved.contains(neighborCoords[1][2]) ||
            (!overridePointsAdded.contains(neighborCoords[1][2]) &&
            !points.contains(neighborCoords[1][2])))
            && 
            (overridePointsRemoved.contains(neighborCoords[2][2]) ||
            (!overridePointsAdded.contains(neighborCoords[2][2]) &&
            !points.contains(neighborCoords[2][2])))
            && 
            (overridePointsRemoved.contains(neighborCoords[2][1]) ||
            (!overridePointsAdded.contains(neighborCoords[2][1]) &&
            !points.contains(neighborCoords[2][1])))
            && 
            (overridePointsRemoved.contains(neighborCoords[2][0]) ||
            (!overridePointsAdded.contains(neighborCoords[2][0]) &&
            !points.contains(neighborCoords[2][0])))
            && 
            (overridePointsRemoved.contains(neighborCoords[1][0]) ||
            (!overridePointsAdded.contains(neighborCoords[1][0]) &&
            !points.contains(neighborCoords[1][0])))
            && 
            (overridePointsRemoved.contains(neighborCoords[0][0]) ||
            (!overridePointsAdded.contains(neighborCoords[0][0]) &&
            !points.contains(neighborCoords[0][0])))
        ) {
            /*
                # . . .
                # @ C .
                # . . .
            */
            if ((col - 2) < 0) {
                return true;
            }
            
            // only 1 of the '#' can be non-zero
            int count = 0;
            PairInt p2 = new PairInt(col - 2, row - 1);
            if (
                (!overridePointsRemoved.contains(p2) &&
                (overridePointsAdded.contains(p2) ||
                points.contains(p2)))
                ) {
                count++;
            }
            p2 = new PairInt(col - 2, row);
            if (
                (!overridePointsRemoved.contains(p2) &&
                (overridePointsAdded.contains(p2) ||
                points.contains(p2)))
                ) {
                count++;
            }
            p2 = new PairInt(col - 2, row + 1);
            if (
                (!overridePointsRemoved.contains(p2) &&
                (overridePointsAdded.contains(p2) ||
                points.contains(p2)))
                ) {
                count++;
            }
            
            if (count > 1) {
                return false;
            }
            return true;
        }
        
        /*
          @ . .
          . C .
          . . .
        */
        if (
            (!overridePointsRemoved.contains(neighborCoords[0][2]) &&
            (overridePointsAdded.contains(neighborCoords[0][2]) ||
            points.contains(neighborCoords[0][2])))
            && 
            (overridePointsRemoved.contains(neighborCoords[1][2]) ||
            (!overridePointsAdded.contains(neighborCoords[1][2]) &&
            !points.contains(neighborCoords[1][2])))
            && 
            (overridePointsRemoved.contains(neighborCoords[2][2]) ||
            (!overridePointsAdded.contains(neighborCoords[2][2]) &&
            !points.contains(neighborCoords[2][2])))
            && 
            (overridePointsRemoved.contains(neighborCoords[2][1]) ||
            (!overridePointsAdded.contains(neighborCoords[2][1]) &&
            !points.contains(neighborCoords[2][1])))
            && 
            (overridePointsRemoved.contains(neighborCoords[2][0]) ||
            (!overridePointsAdded.contains(neighborCoords[2][0]) &&
            !points.contains(neighborCoords[2][0])))
            && 
            (overridePointsRemoved.contains(neighborCoords[1][0]) ||
            (!overridePointsAdded.contains(neighborCoords[1][0]) &&
            !points.contains(neighborCoords[1][0])))
            && 
            (overridePointsRemoved.contains(neighborCoords[0][0]) ||
            (!overridePointsAdded.contains(neighborCoords[0][0]) &&
            !points.contains(neighborCoords[0][0])))
            && 
            (overridePointsRemoved.contains(neighborCoords[0][1]) ||
            (!overridePointsAdded.contains(neighborCoords[0][1]) &&
            !points.contains(neighborCoords[0][1])))
        ) {
            if ((col - 2) < 0) {
                if ((row + 2) > (imageHeight - 1)) {
                    return true;
                }
                /*
                     # #
                     @ . .
                     . C .
                     . . .
                */
                // only 1 of the '#' can be non-zero
                int count = 0;
                PairInt p2 = new PairInt(col - 1, row + 2);
                if (
                    (!overridePointsRemoved.contains(p2) &&
                    (overridePointsAdded.contains(p2) ||
                    points.contains(p2)))
                    ) {
                    count++;
                }
                p2 = new PairInt(col, row + 2);
                if (
                    (!overridePointsRemoved.contains(p2) &&
                    (overridePointsAdded.contains(p2) ||
                    points.contains(p2)))
                    ) {
                    count++;
                }

                if (count > 1) {
                    return false;
                }
                return true;
            } else if ((row + 2) > (imageHeight - 1)) {
                /*
                    # @ . .
                    # . C .
                    . . .
                 */
                int count = 0;
                PairInt p2 = new PairInt(col - 2, row);
                if (
                    (!overridePointsRemoved.contains(p2) &&
                    (overridePointsAdded.contains(p2) ||
                    points.contains(p2)))
                    ) {
                    count++;
                }
                p2 = new PairInt(col - 2, row + 1);
                if (
                    (!overridePointsRemoved.contains(p2) &&
                    (overridePointsAdded.contains(p2) ||
                    points.contains(p2)))
                    ) {
                    count++;
                }

                if (count > 1) {
                    return false;
                }
                return true;
            }
            /*
               # # #
               # @ . .
               # . C .
                 . . .
            */
            int count = 0;
            PairInt p2 = new PairInt(col - 2, row);
            if (
                (!overridePointsRemoved.contains(p2) &&
                (overridePointsAdded.contains(p2) ||
                points.contains(p2)))
                ) {
                count++;
            }
            p2 = new PairInt(col - 2, row + 1);
            if (
                (!overridePointsRemoved.contains(p2) &&
                (overridePointsAdded.contains(p2) ||
                points.contains(p2)))
                ) {
                count++;
            }
            p2 = new PairInt(col - 2, row + 2);
            if (
                (!overridePointsRemoved.contains(p2) &&
                (overridePointsAdded.contains(p2) ||
                points.contains(p2)))
                ) {
                count++;
            }
            p2 = new PairInt(col - 1, row + 2);
            if (
                (!overridePointsRemoved.contains(p2) &&
                (overridePointsAdded.contains(p2) ||
                points.contains(p2)))
                ) {
                count++;
            }
            p2 = new PairInt(col, row + 2);
            if (
                (!overridePointsRemoved.contains(p2) &&
                (overridePointsAdded.contains(p2) ||
                points.contains(p2)))
                ) {
                count++;
            }

            if (count > 1) {
                return false;
            }
            return true;
        }
        
        return false;
    }

    /**
     * for the full 8 neighbor region, determine whether nulling the pixel
     * at (col, row) would disconnect the remaining line.  Note that the
     * boolean logic is embedded in the comments.  One should be able to
     * combine the rules for multiple pixel tests to reduce the redundant
     * comparisons for the regions in common.
     * 
     * Note, that row and col are expected to be at least 1 pixel distant
     * from the image borders.
    
     * @return 
     */
    private boolean doesDisconnect(PairInt p, Set<PairInt> points, 
        Set<PairInt> overridePointsAdded, Set<PairInt> overridePointsRemoved, 
        int imageWidth, int imageHeight) {
       
        int col = p.getX();
        int row = p.getY();
        
        /*
        coordinates of the 8 neighbors as already created PairInts without 
        bound checks.
        indexes are found as +1 of the difference relative to center,
        for example, a point for (col-1, row-1) is found as neighborCoords[0][0]
        */
        PairInt[][] neighborCoords = createCoordinatePointsForEightNeighbors(
            col, row);
        
        /*
            6   7  8 
           11 *C* 12
           15  16 17
        */
                
        //6  Ʌ  8  Ʌ  ¬((7)  V  (11  Ʌ  16  Ʌ  12))
        if (
            (
               (!overridePointsRemoved.contains(neighborCoords[0][2]) &&
               (overridePointsAdded.contains(neighborCoords[0][2]) ||
               points.contains(neighborCoords[0][2])))
               && 
               (!overridePointsRemoved.contains(neighborCoords[2][2]) &&
               (overridePointsAdded.contains(neighborCoords[2][2]) ||
               points.contains(neighborCoords[2][2])))
            ) && 
            !(
               (!overridePointsRemoved.contains(neighborCoords[1][2]) &&
               (overridePointsAdded.contains(neighborCoords[1][2]) ||
               points.contains(neighborCoords[1][2])))
               || 
               (
                   (!overridePointsRemoved.contains(neighborCoords[0][1]) &&
                   (overridePointsAdded.contains(neighborCoords[0][1]) ||
                   points.contains(neighborCoords[0][1])))
                   && 
                   (!overridePointsRemoved.contains(neighborCoords[1][0]) &&
                   (overridePointsAdded.contains(neighborCoords[1][0]) ||
                   points.contains(neighborCoords[1][0])))
                   && 
                   (!overridePointsRemoved.contains(neighborCoords[2][1]) &&
                   (overridePointsAdded.contains(neighborCoords[2][1]) ||
                   points.contains(neighborCoords[2][1])))
            ))) {
            return true;
        } else 
            //6  Ʌ  12  Ʌ  ¬((7)  V  (11  Ʌ  16))
            if (
                (
                    (!overridePointsRemoved.contains(neighborCoords[0][2]) &&
                    (overridePointsAdded.contains(neighborCoords[0][2]) ||
                    points.contains(neighborCoords[0][2])))
                    && 
                    (!overridePointsRemoved.contains(neighborCoords[2][1]) &&
                    (overridePointsAdded.contains(neighborCoords[2][1]) ||
                    points.contains(neighborCoords[2][1])))
                ) 
                && 
                !(
                    (!overridePointsRemoved.contains(neighborCoords[1][2]) &&
                    (overridePointsAdded.contains(neighborCoords[1][2]) ||
                    points.contains(neighborCoords[1][2])))
                    || 
                    (
                        (!overridePointsRemoved.contains(neighborCoords[0][1]) &&
                        (overridePointsAdded.contains(neighborCoords[0][1]) ||
                        points.contains(neighborCoords[0][1])))
                        && 
                        (!overridePointsRemoved.contains(neighborCoords[1][0]) &&
                        (overridePointsAdded.contains(neighborCoords[1][0]) ||
                        points.contains(neighborCoords[1][0])))
                ))) {
            return true;
        } else 
            //6  Ʌ  17  Ʌ  ¬((7  Ʌ  12) V (11  Ʌ  16))
            if (
                (
                    (!overridePointsRemoved.contains(neighborCoords[0][2]) &&
                   (overridePointsAdded.contains(neighborCoords[0][2]) ||
                   points.contains(neighborCoords[0][2])))
                    && 
                    (!overridePointsRemoved.contains(neighborCoords[2][0]) &&
                   (overridePointsAdded.contains(neighborCoords[2][0]) ||
                   points.contains(neighborCoords[2][0])))
                ) 
                && 
                !(
                    (
                        (!overridePointsRemoved.contains(neighborCoords[1][2]) &&
                        (overridePointsAdded.contains(neighborCoords[1][2]) ||
                        points.contains(neighborCoords[1][2])))
                        && 
                        (!overridePointsRemoved.contains(neighborCoords[2][1]) &&
                        (overridePointsAdded.contains(neighborCoords[2][1]) ||
                        points.contains(neighborCoords[2][1])))
                    ) 
                    || 
                    (
                        (!overridePointsRemoved.contains(neighborCoords[0][1]) &&
                        (overridePointsAdded.contains(neighborCoords[0][1]) ||
                        points.contains(neighborCoords[0][1])))
                        && 
                        (!overridePointsRemoved.contains(neighborCoords[1][0]) &&
                        (overridePointsAdded.contains(neighborCoords[1][0]) ||
                        points.contains(neighborCoords[1][0])))
                ))) {
            return true;
        } else 
            //6  Ʌ  16  Ʌ  ¬((7  Ʌ  12) V (11))
            if (
                (
                    (!overridePointsRemoved.contains(neighborCoords[0][2]) &&
                    (overridePointsAdded.contains(neighborCoords[0][2]) ||
                    points.contains(neighborCoords[0][2])))
                    && 
                    (!overridePointsRemoved.contains(neighborCoords[1][0]) &&
                    (overridePointsAdded.contains(neighborCoords[1][0]) ||
                    points.contains(neighborCoords[1][0])))
                ) 
                && 
                !(
                    (
                        (!overridePointsRemoved.contains(neighborCoords[1][2]) &&
                        (overridePointsAdded.contains(neighborCoords[1][2]) ||
                        points.contains(neighborCoords[1][2])))
                        && 
                        (!overridePointsRemoved.contains(neighborCoords[2][1]) &&
                        (overridePointsAdded.contains(neighborCoords[2][1]) ||
                        points.contains(neighborCoords[2][1])))
                    ) 
                    || 
                    (!overridePointsRemoved.contains(neighborCoords[0][1]) &&
                    (overridePointsAdded.contains(neighborCoords[0][1]) ||
                    points.contains(neighborCoords[0][1])))
                )) {
            return true;
        } else 
            //6  Ʌ  15  Ʌ  ¬((7  Ʌ  12 Ʌ 16) V (11))
            if (
                (
                    (!overridePointsRemoved.contains(neighborCoords[0][2]) &&
                    (overridePointsAdded.contains(neighborCoords[0][2]) ||
                    points.contains(neighborCoords[0][2])))
                    && 
                    (!overridePointsRemoved.contains(neighborCoords[0][0]) &&
                    (overridePointsAdded.contains(neighborCoords[0][0]) ||
                    points.contains(neighborCoords[0][0])))
                ) 
                && 
                !(
                    (
                        (!overridePointsRemoved.contains(neighborCoords[1][2]) &&
                        (overridePointsAdded.contains(neighborCoords[1][2]) ||
                        points.contains(neighborCoords[1][2])))
                        && 
                        (!overridePointsRemoved.contains(neighborCoords[2][1]) &&
                        (overridePointsAdded.contains(neighborCoords[2][1]) ||
                        points.contains(neighborCoords[2][1])))
                        && 
                        (!overridePointsRemoved.contains(neighborCoords[1][0]) &&
                        (overridePointsAdded.contains(neighborCoords[1][0]) ||
                        points.contains(neighborCoords[1][0])))
                    ) 
                    || 
                    (!overridePointsRemoved.contains(neighborCoords[0][1]) &&
                    (overridePointsAdded.contains(neighborCoords[0][1]) ||
                    points.contains(neighborCoords[0][1])))
                )) {
            return true;
        } else 
            //7  Ʌ  17  Ʌ  ¬((12)  V  (11  Ʌ  16))
            if (
                (
                    (!overridePointsRemoved.contains(neighborCoords[1][2]) &&
                    (overridePointsAdded.contains(neighborCoords[1][2]) ||
                    points.contains(neighborCoords[1][2])))
                    && 
                    (!overridePointsRemoved.contains(neighborCoords[2][0]) &&
                    (overridePointsAdded.contains(neighborCoords[2][0]) ||
                    points.contains(neighborCoords[2][0])))
                ) 
                && 
                !(
                    (
                        (!overridePointsRemoved.contains(neighborCoords[0][1]) &&
                        (overridePointsAdded.contains(neighborCoords[0][1]) ||
                        points.contains(neighborCoords[0][1])))
                        && 
                        (!overridePointsRemoved.contains(neighborCoords[1][0]) &&
                        (overridePointsAdded.contains(neighborCoords[1][0]) ||
                        points.contains(neighborCoords[1][0])))
                    ) 
                    || 
                    (!overridePointsRemoved.contains(neighborCoords[2][1]) &&
                    (overridePointsAdded.contains(neighborCoords[2][1]) ||
                    points.contains(neighborCoords[2][1])))
                )) {
            return true;
        } else 
            //7  Ʌ  16  Ʌ  ¬(12  V  11)
            if (
                (
                    (!overridePointsRemoved.contains(neighborCoords[1][2]) &&
                    (overridePointsAdded.contains(neighborCoords[1][2]) ||
                    points.contains(neighborCoords[1][2])))
                    && 
                    (!overridePointsRemoved.contains(neighborCoords[1][0]) &&
                    (overridePointsAdded.contains(neighborCoords[1][0]) ||
                    points.contains(neighborCoords[1][0])))
                ) 
                && 
                !(
                    (!overridePointsRemoved.contains(neighborCoords[2][1]) &&
                    (overridePointsAdded.contains(neighborCoords[2][1]) ||
                    points.contains(neighborCoords[2][1])))
                    || 
                    (!overridePointsRemoved.contains(neighborCoords[0][1]) &&
                    (overridePointsAdded.contains(neighborCoords[0][1]) ||
                    points.contains(neighborCoords[0][1])))
                )) {
            return true;
        } else 
            //7  Ʌ  15  Ʌ  ¬((11)  V  (12 Ʌ 16))
            if (
                (
                    (!overridePointsRemoved.contains(neighborCoords[1][2]) &&
                    (overridePointsAdded.contains(neighborCoords[1][2]) ||
                    points.contains(neighborCoords[1][2])))
                    && 
                    (!overridePointsRemoved.contains(neighborCoords[0][0]) &&
                    (overridePointsAdded.contains(neighborCoords[0][0]) ||
                    points.contains(neighborCoords[0][0])))
                ) 
                && 
                !(
                    (
                        (!overridePointsRemoved.contains(neighborCoords[2][1]) &&
                        (overridePointsAdded.contains(neighborCoords[2][1]) ||
                        points.contains(neighborCoords[2][1])))
                        && 
                        (!overridePointsRemoved.contains(neighborCoords[1][0]) &&
                        (overridePointsAdded.contains(neighborCoords[1][0]) ||
                        points.contains(neighborCoords[1][0])))
                    ) 
                    || 
                    (!overridePointsRemoved.contains(neighborCoords[0][1]) &&
                    (overridePointsAdded.contains(neighborCoords[0][1]) ||
                    points.contains(neighborCoords[0][1])))
                )) {
            return true;
        } else 
            //8  Ʌ  17  Ʌ  ¬((12)  V  (7 Ʌ 11  Ʌ  16))
            if (
                (
                    (!overridePointsRemoved.contains(neighborCoords[2][2]) &&
                    (overridePointsAdded.contains(neighborCoords[2][2]) ||
                    points.contains(neighborCoords[2][2])))
                    && 
                    (!overridePointsRemoved.contains(neighborCoords[2][0]) &&
                    (overridePointsAdded.contains(neighborCoords[2][0]) ||
                    points.contains(neighborCoords[2][0])))
                ) 
                && 
                !(
                    (
                        (!overridePointsRemoved.contains(neighborCoords[1][2]) &&
                        (overridePointsAdded.contains(neighborCoords[1][2]) ||
                        points.contains(neighborCoords[1][2])))
                        && 
                        (!overridePointsRemoved.contains(neighborCoords[0][1]) &&
                        (overridePointsAdded.contains(neighborCoords[0][1]) ||
                        points.contains(neighborCoords[0][1])))
                        && 
                        (!overridePointsRemoved.contains(neighborCoords[1][0]) &&
                        (overridePointsAdded.contains(neighborCoords[1][0]) ||
                        points.contains(neighborCoords[1][0])))
                    ) 
                    || 
                    (!overridePointsRemoved.contains(neighborCoords[2][1]) &&
                    (overridePointsAdded.contains(neighborCoords[2][1]) ||
                    points.contains(neighborCoords[2][1])))
                )) {
            return true;
        } else 
            //8  Ʌ  16  Ʌ  ¬((12)  V  (7 Ʌ 11))
            if (
                (
                    (!overridePointsRemoved.contains(neighborCoords[2][2]) &&
                    (overridePointsAdded.contains(neighborCoords[2][2]) ||
                    points.contains(neighborCoords[2][2])))
                    && 
                    (!overridePointsRemoved.contains(neighborCoords[1][0]) &&
                    (overridePointsAdded.contains(neighborCoords[1][0]) ||
                    points.contains(neighborCoords[1][0])))
                ) 
                && 
                !(
                    (
                        (!overridePointsRemoved.contains(neighborCoords[1][2]) &&
                        (overridePointsAdded.contains(neighborCoords[1][2]) ||
                        points.contains(neighborCoords[1][2])))
                        && 
                        (!overridePointsRemoved.contains(neighborCoords[0][1])
                        && (overridePointsAdded.contains(neighborCoords[0][1])
                        || points.contains(neighborCoords[0][1])))
                    ) 
                    || 
                    (!overridePointsRemoved.contains(neighborCoords[2][1]) &&
                    (overridePointsAdded.contains(neighborCoords[2][1]) ||
                    points.contains(neighborCoords[2][1])))
                )) {
            return true;
        } else 
            //8  Ʌ  15  Ʌ  ¬((12 Ʌ  16)  V  (7 Ʌ 11))
            if (
                (
                    (!overridePointsRemoved.contains(neighborCoords[2][2]) &&
                    (overridePointsAdded.contains(neighborCoords[2][2]) ||
                    points.contains(neighborCoords[2][2])))
                    && 
                    (!overridePointsRemoved.contains(neighborCoords[0][0]) &&
                    (overridePointsAdded.contains(neighborCoords[0][0]) ||
                    points.contains(neighborCoords[0][0])))
                ) 
                && 
                !(
                    (
                        (!overridePointsRemoved.contains(neighborCoords[2][1]) &&
                        (overridePointsAdded.contains(neighborCoords[2][1]) ||
                        points.contains(neighborCoords[2][1])))
                        && 
                        (!overridePointsRemoved.contains(neighborCoords[1][0]) &&
                        (overridePointsAdded.contains(neighborCoords[1][0]) ||
                        points.contains(neighborCoords[1][0])))
                    ) 
                    || 
                    (
                        (!overridePointsRemoved.contains(neighborCoords[1][2]) &&
                        (overridePointsAdded.contains(neighborCoords[1][2]) ||
                        points.contains(neighborCoords[1][2])))
                        && 
                        (!overridePointsRemoved.contains(neighborCoords[0][1]) &&
                        (overridePointsAdded.contains(neighborCoords[0][1]) ||
                        points.contains(neighborCoords[0][1])))
                ))) {
            return true;
        } else 
            //8  Ʌ  11  Ʌ  ¬((7)  V  (12 Ʌ 16))
            if (
                (
                    (!overridePointsRemoved.contains(neighborCoords[2][2]) &&
                    (overridePointsAdded.contains(neighborCoords[2][2]) ||
                    points.contains(neighborCoords[2][2])))
                    && 
                    (!overridePointsRemoved.contains(neighborCoords[0][1]) &&
                    (overridePointsAdded.contains(neighborCoords[0][1]) ||
                    points.contains(neighborCoords[0][1])))
                ) 
                && 
                !(
                    (
                        (!overridePointsRemoved.contains(neighborCoords[2][1]) &&
                        (overridePointsAdded.contains(neighborCoords[2][1]) ||
                        points.contains(neighborCoords[2][1])))
                        && 
                        (!overridePointsRemoved.contains(neighborCoords[1][0]) &&
                        (overridePointsAdded.contains(neighborCoords[1][0]) ||
                        points.contains(neighborCoords[1][0])))
                    ) 
                    || 
                    (!overridePointsRemoved.contains(neighborCoords[1][2]) &&
                    (overridePointsAdded.contains(neighborCoords[1][2]) ||
                    points.contains(neighborCoords[1][2])))
                )) {
            return true;
        } else 
            //12  Ʌ  15  Ʌ  ¬((16)  V  (7 Ʌ 11))
            if (
                (
                    (!overridePointsRemoved.contains(neighborCoords[2][1]) &&
                    (overridePointsAdded.contains(neighborCoords[2][1]) ||
                    points.contains(neighborCoords[2][1])))
                    && 
                    (!overridePointsRemoved.contains(neighborCoords[0][0]) &&
                    (overridePointsAdded.contains(neighborCoords[0][0]) ||
                    points.contains(neighborCoords[0][0])))
                ) 
                && 
                !(
                    (
                        (!overridePointsRemoved.contains(neighborCoords[1][2]) &&
                        (overridePointsAdded.contains(neighborCoords[1][2]) ||
                        points.contains(neighborCoords[1][2])))
                        && 
                        (!overridePointsRemoved.contains(neighborCoords[0][1]) &&
                        (overridePointsAdded.contains(neighborCoords[0][1]) ||
                        points.contains(neighborCoords[0][1])))
                    ) 
                    || 
                    (!overridePointsRemoved.contains(neighborCoords[1][0]) &&
                    (overridePointsAdded.contains(neighborCoords[1][0]) ||
                    points.contains(neighborCoords[1][0])))
                )) {
            return true;
        } else 
            //12  Ʌ  11  Ʌ  ¬((16)  V  (7))
            if (
                (
                    (!overridePointsRemoved.contains(neighborCoords[2][1]) &&
                    (overridePointsAdded.contains(neighborCoords[2][1]) ||
                    points.contains(neighborCoords[2][1])))
                    && 
                    (!overridePointsRemoved.contains(neighborCoords[0][1]) &&
                    (overridePointsAdded.contains(neighborCoords[0][1]) ||
                    points.contains(neighborCoords[0][1])))
                ) 
                && 
                !(
                    (!overridePointsRemoved.contains(neighborCoords[1][2]) &&
                    (overridePointsAdded.contains(neighborCoords[1][2]) ||
                    points.contains(neighborCoords[1][2])))
                    || 
                    (!overridePointsRemoved.contains(neighborCoords[1][0]) &&
                    (overridePointsAdded.contains(neighborCoords[1][0]) ||
                    points.contains(neighborCoords[1][0])))
                )) {
            return true;
        } else 
            //17  Ʌ  15  Ʌ  ¬((16)  V  (7 Ʌ 11 Ʌ 12))
            if (
                (
                    (!overridePointsRemoved.contains(neighborCoords[2][0]) &&
                    (overridePointsAdded.contains(neighborCoords[2][0]) ||
                    points.contains(neighborCoords[2][0])))
                    && 
                    (!overridePointsRemoved.contains(neighborCoords[0][0]) &&
                    (overridePointsAdded.contains(neighborCoords[0][0]) ||
                    points.contains(neighborCoords[0][0])))
                ) 
                && 
                !(
                    (
                        (!overridePointsRemoved.contains(neighborCoords[1][2]) &&
                        (overridePointsAdded.contains(neighborCoords[1][2]) ||
                        points.contains(neighborCoords[1][2])))
                        && 
                        (!overridePointsRemoved.contains(neighborCoords[0][1]) &&
                        (overridePointsAdded.contains(neighborCoords[0][1]) ||
                        points.contains(neighborCoords[0][1])))
                        && 
                        (!overridePointsRemoved.contains(neighborCoords[2][1]) &&
                        (overridePointsAdded.contains(neighborCoords[2][1]) ||
                        points.contains(neighborCoords[2][1])))
                    ) 
                    || 
                    (!overridePointsRemoved.contains(neighborCoords[1][0]) &&
                    (overridePointsAdded.contains(neighborCoords[1][0]) ||
                    points.contains(neighborCoords[1][0])))
                )) {
            return true;
        } else 
            //17  Ʌ  11  Ʌ  ¬((16)  V  (7 Ʌ 12))
            if (
                (
                    (!overridePointsRemoved.contains(neighborCoords[2][0]) &&
                    (overridePointsAdded.contains(neighborCoords[2][0]) ||
                    points.contains(neighborCoords[2][0])))
                    && 
                    (!overridePointsRemoved.contains(neighborCoords[0][1]) &&
                    (overridePointsAdded.contains(neighborCoords[0][1]) ||
                    points.contains(neighborCoords[0][1])))
                ) 
                && 
                !(
                    (
                        (!overridePointsRemoved.contains(neighborCoords[1][2]) &&
                        (overridePointsAdded.contains(neighborCoords[1][2]) ||
                        points.contains(neighborCoords[1][2])))
                        && 
                        (!overridePointsRemoved.contains(neighborCoords[2][1]) &&
                        (overridePointsAdded.contains(neighborCoords[2][1]) ||
                        points.contains(neighborCoords[2][1])))
                    ) 
                    || 
                    (!overridePointsRemoved.contains(neighborCoords[1][0]) &&
                    (overridePointsAdded.contains(neighborCoords[1][0]) ||
                    points.contains(neighborCoords[1][0])))
                )) {
            return true;
        }
        
        return false;
    }

}
