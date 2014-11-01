package algorithms.imageProcessing;

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
        
        int count = 1;
               
        for (int i = 0; i < input.getWidth(); i++) {

            for (int j = 0; j < input.getHeight(); j++) {

                if (input.getValue(i, j) > 0) {
                    continue;
                }
                       
                if (!isASmallerSearchRegion(input, i, j) &&
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
       
        boolean isASmallerSearchRegion = isASmallerSearchRegion(
            input, col, row);
        
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
     * for the full 8 neighbor region, determine whether nulling the pixel
     * at (col, row) would disconnect the remaining line.  Note that the
     * boolean logic is embedded in the comments.  One should be able to
     * combine the rules for multiple pixel tests to reduce the redundant
     * comparisons for the regions in common.
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
     * other currently connected points
     * @param input
     * @param output
     * @param i
     * @param j
     * @return 
     */
    private boolean process(final GreyscaleImage input, final GreyscaleImage 
        output, int i, int j) {
        
        if (output.getValue(i, j) == 0) {
            return false;
        }

        if (!isASmallerSearchRegion(output, i, j)
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

}
