package algorithms.imageProcessing;

import algorithms.util.PairInt;
import java.util.HashSet;
import java.util.Set;

/**
 *
 * @author nichole
 */
public abstract class AbstractLineThinner implements ILineThinner {
    
    protected static final int[] fourNeighborsX = new int[]{0,  1,  0,  -1};
    protected static final int[] fourNeighborsY = new int[]{-1, 0,  1,   0};
    
    protected static final int[] eightNeighborsX = 
        new int[]{-1, -1,  0,  1, 1, 1, 0, -1};
    protected static final int[] eightNeighborsY = 
        new int[]{ 0, -1, -1, -1, 0, 1, 1,  1};

    protected boolean useLineDrawingMode = false;
    
    protected boolean debug = false;
    
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
    
    protected boolean isASmallerSearchRegionAndDoesDisconnect(
        final GreyscaleImage input, int col, int row) {
        
        if ((col - 1) < 0) {
            /*
            | @ @
            | C @
            | @ @
            */
            if ((col + 1) > (input.getWidth() - 1)) {
                /*
                | @ |
                | C |
                | @ |
                */
                // region is 1 column thin, so don't delete the pix no matter
                // what its position
                return true;
            }
            if ((row - 1) < 0) {
                /*
                | @ @
                | C @
                */
                if ((row + 1) > (input.getHeight() - 1)) {
                    /*
                    | C @
                    */
                    return true;
                }
                if ((input.getValue(col, row + 1) > 0) || 
                    (input.getValue(col + 1, row + 1) > 0) || 
                    (input.getValue(col + 1, row) > 0)) {
                    return true;
                }
                return false;
            } else if ((row + 1) > (input.getHeight() - 1)) {
                /*
                | C @
                | @ @
                */
                if ((input.getValue(col, row - 1) > 0) || (input.getValue(col + 1, row - 1) > 0) || (input.getValue(col + 1, row) > 0)) {
                    return true;
                }
                return false;
            }
            /*
            | @ @
            | C @
            | @ @
             */
            /*
            | # _
            | C .     '#' is disconnected if All .'s are 0s
            | ? ?     AND AtLeastOne (? is 1)
            |
            | _ #
            | C .     '#' is disconnected if All .'s are 0s
            | ? ?     AND AtLeastOne (? is 1)
             */
            if (((input.getValue(col, row + 1) > 0) || (input.getValue(col + 1, row + 1) > 0)) && (input.getValue(col + 1, row) == 0)) {
                if ((input.getValue(col, row - 1) > 0) || (input.getValue(col + 1, row - 1) > 0)) {
                    return true;
                }
            }
            /*
            | . _
            | C #   '#' is disconnected if All .'s are 0s
            | . _
             */
            if (input.getValue(col + 1, row) > 0) {
                if ((input.getValue(col, row + 1) == 0) || (input.getValue(col, row - 1) == 0)) {
                    return true;
                }
            }
            /*
            | ? ?
            | C .     '#' is disconnected if All .'s are 0s
            | _ #     AND AtLeastOne (? is 1)
            |
            | ? ?
            | C .     '#' is disconnected if All .'s are 0s
            | # _     AND AtLeastOne (? is 1)
             */
            if (((input.getValue(col, row - 1) > 0) || (input.getValue(col + 1, row - 1) > 0)) && (input.getValue(col + 1, row) == 0)) {
                if ((input.getValue(col, row + 1) > 0) || (input.getValue(col + 1, row + 1) > 0)) {
                    return true;
                }
            }
            return false;
        } else if ((col + 1) > (input.getWidth() - 1)) {
            /*
            @ @ |
            @ C |
            @ @ |
             */
            // already tested that col-1 exists
            if ((row - 1) < 0) {
                /*
                @ @ |
                @ C |
                 */
                if ((row + 1) > (input.getHeight() - 1)) {
                    /*
                    @ C |
                     */
                    return true;
                }
                if ((input.getValue(col - 1, row) > 0) || (input.getValue(col - 1, row + 1) > 0) || (input.getValue(col, row + 1) > 0)) {
                    return true;
                }
                return false;
            } else if ((row + 1) > (input.getHeight() - 1)) {
                /*
                @ C |
                @ @ |
                 */
                if ((input.getValue(col - 1, row) > 0) || (input.getValue(col - 1, row - 1) > 0) || (input.getValue(col, row - 1) > 0)) {
                    return true;
                }
                return false;
            }
            /*
            @ @ |
            @ C |
            @ @ |
             */
            /*
            - # |
            . C | '#' is disconnected if All .'s are 0s
            ? ? |   AND AtLeastOne (? is 1)
            # _ |
            . C | '#' is disconnected if All .'s are 0s
            ? ? |   AND AtLeastOne (? is 1)
             */
            if (((input.getValue(col, row + 1) > 0) || (input.getValue(col - 1, row + 1) > 0)) && (input.getValue(col - 1, row) == 0)) {
                if ((input.getValue(col, row - 1) > 0) || (input.getValue(col - 1, row - 1) > 0)) {
                    return true;
                }
            }
            /*
            _ . |
            # C |  '#' is disconnected if All .'s are 0s
            _ . |
             */
            if (input.getValue(col - 1, row) > 0) {
                if ((input.getValue(col, row + 1) == 0) || (input.getValue(col, row - 1) == 0)) {
                    return true;
                }
            }
            /*
            ? ? |
            . C | '#' is disconnected if All .'s are 0s
            - # |   AND AtLeastOne (? is 1)
            ? ? |
            . C | '#' is disconnected if All .'s are 0s
            # _ |   AND AtLeastOne (? is 1)
             */
            if (((input.getValue(col, row - 1) > 0) || (input.getValue(col - 1, row - 1) > 0)) && (input.getValue(col - 1, row) == 0)) {
                if ((input.getValue(col, row + 1) > 0) || (input.getValue(col - 1, row + 1) > 0)) {
                    return true;
                }
            }
            return false;
        } else if ((row - 1) < 0) {
            /*
            @ @ @
            @ C @
             */
            // have already tested that column 0 and column 2 are present
            /*
            ?  .  #    '#' is disconnected if All .'s are 0s
            ?  C  _        AND AtLeastOne (? is 1)
            ?  .  _    '#' is disconnected if All .'s are 0s
            ?  C  #        AND AtLeastOne (? is 1)
             */
            if (((input.getValue(col + 1, row) > 0) || (input.getValue(col + 1, row + 1) > 0)) && (input.getValue(col, row + 1) == 0)) {
                if ((input.getValue(col - 1, row + 1) > 0) || (input.getValue(col - 1, row) > 0)) {
                    return true;
                }
            }
            /*
            -  #  _
            .  C  .    '#' is disconnected if All .'s are 0s
             */
            if (input.getValue(col, row + 1) > 0) {
                if ((input.getValue(col - 1, row) == 0) || (input.getValue(col + 1, row) == 0)) {
                    return true;
                }
            }
            /*
            #  .  ?    '#' is disconnected if All .'s are 0s
            _  C  ?       AND AtLeastOne (? is 1)
            _  .  ?    '#' is disconnected if All .'s are 0s
            #  C  ?        AND AtLeastOne (? is 1)
             */
            if (((input.getValue(col - 1, row) > 0) || (input.getValue(col - 1, row + 1) > 0)) && (input.getValue(col, row + 1) == 0)) {
                if ((input.getValue(col + 1, row + 1) > 0) || (input.getValue(col + 1, row) > 0)) {
                    return true;
                }
            }
            return false;
        } else if ((row + 1) > (input.getHeight() - 1)) {
            /*
            @ C @
            @ @ @
             */
            // have already tested that column 0 and column 2 are present
            // and that the bottom row is present
            /*
            ?  C  _    '#' is disconnected if All .'s are 0s
            ?  .  #       AND AtLeastOne (? is 1)
            ?  C  #    '#' is disconnected if All .'s are 0s
            ?  .  _       AND AtLeastOne (? is 1)
             */
            if (((input.getValue(col + 1, row - 1) > 0) || (input.getValue(col + 1, row) > 0)) && (input.getValue(col, row - 1) == 0)) {
                if ((input.getValue(col - 1, row) > 0) || (input.getValue(col - 1, row - 1) > 0)) {
                    return true;
                }
            }
            /*
            .  C  .
            -  #  _  '#' is disconnected if All .'s are 0s
             */
            if (input.getValue(col, row - 1) > 0) {
                if ((input.getValue(col - 1, row) == 0) || (input.getValue(col + 1, row) == 0)) {
                    return true;
                }
            }
            /*
            -  C  ?    '#' is disconnected if All .'s are 0s
            #  .  ?       AND AtLeastOne (? is 1)
            #  C  ?    '#' is disconnected if All .'s are 0s
            _  .  ?      AND AtLeastOne (? is 1)
             */
            if (((input.getValue(col - 1, row - 1) > 0) || (input.getValue(col - 1, row) > 0)) && (input.getValue(col, row - 1) == 0)) {
                if ((input.getValue(col + 1, row) > 0) || (input.getValue(col + 1, row - 1) > 0)) {
                    return true;
                }
            }
            return false;
        }
        return false;
    }
    
    public boolean isASmallerSearchRegionAndDoesDisconnect(PairInt p, 
        Set<PairInt> points, Set<PairInt> overridePointsAdded, 
        Set<PairInt> overridePointsRemoved, int imageWidth, int imageHeight) {
        
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
        
        if ((col - 1) < 0) {
            /*
             | @ @
             | C @
             | @ @
             */
            if ((col + 1) > (imageWidth - 1)) {
                /*
                 | @ |
                 | C |
                 | @ |
                 */
                // region is 1 column thin, so don't delete the pix no matter
                // what its position
                return true;
            }
            if ((row - 1) < 0) {
                /*
                 | @ @
                 | C @
                 */
                if ((row + 1) > (imageHeight - 1)) {
                    /*
                     | C @
                     */
                    return true;
                }
                
                if (
                    (!overridePointsRemoved.contains(neighborCoords[1][2]) &&
                    (overridePointsAdded.contains(neighborCoords[1][2]) ||
                    points.contains(neighborCoords[1][2])))
                    || 
                    (!overridePointsRemoved.contains(neighborCoords[2][2]) &&
                    (overridePointsAdded.contains(neighborCoords[2][2]) ||
                    points.contains(neighborCoords[2][2])))
                    ||
                    (!overridePointsRemoved.contains(neighborCoords[2][1]) &&
                    (overridePointsAdded.contains(neighborCoords[2][1]) ||
                    points.contains(neighborCoords[2][1])))) {
                    
                    return true;
                }
                
                return false;
            } else if ((row + 1) > (imageHeight - 1)) {
                /*
                 | C @
                 | @ @
                 */
                if (
                    (!overridePointsRemoved.contains(neighborCoords[1][0]) &&
                    (overridePointsAdded.contains(neighborCoords[1][0]) ||
                    points.contains(neighborCoords[1][0])))
                    || 
                    (!overridePointsRemoved.contains(neighborCoords[2][0]) &&
                    (overridePointsAdded.contains(neighborCoords[2][0]) ||
                    points.contains(neighborCoords[2][0])))
                    ||
                    (!overridePointsRemoved.contains(neighborCoords[2][1]) &&
                    (overridePointsAdded.contains(neighborCoords[2][1]) ||
                    points.contains(neighborCoords[2][1])))) {
                    
                    return true;
                }
                
                return false;
            }
            /*
             | @ @
             | C @
             | @ @
             */
            /*
             | # _
             | C .     '#' is disconnected if All .'s are 0s
             | ? ?     AND AtLeastOne (? is 1)
             |
             | _ #
             | C .     '#' is disconnected if All .'s are 0s
             | ? ?     AND AtLeastOne (? is 1)
             */
            if (
                (
                    (!overridePointsRemoved.contains(neighborCoords[1][2]) &&
                    (overridePointsAdded.contains(neighborCoords[1][2]) ||
                    points.contains(neighborCoords[1][2])))
                    || 
                    (!overridePointsRemoved.contains(neighborCoords[2][2]) &&
                    (overridePointsAdded.contains(neighborCoords[2][2]) ||
                    points.contains(neighborCoords[2][2])))
                ) 
                && 
                    (overridePointsRemoved.contains(neighborCoords[2][1]) ||
                    (!overridePointsAdded.contains(neighborCoords[2][1]) &&
                    !points.contains(neighborCoords[2][1])))
                ) {
                
                if (
                    (!overridePointsRemoved.contains(neighborCoords[1][0]) &&
                    (overridePointsAdded.contains(neighborCoords[1][0]) ||
                    points.contains(neighborCoords[1][0])))
                    || 
                    (!overridePointsRemoved.contains(neighborCoords[2][0]) &&
                    (overridePointsAdded.contains(neighborCoords[2][0]) ||
                    points.contains(neighborCoords[2][0])))
                    ) {
                    
                    return true;
                }
            }
            /*
             | . _
             | C #   '#' is disconnected if All .'s are 0s
             | . _
             */
            if (
                (!overridePointsRemoved.contains(neighborCoords[2][1]) &&
                (overridePointsAdded.contains(neighborCoords[2][1]) ||
                points.contains(neighborCoords[2][1])))
                ) {
                
                if (
                    (overridePointsRemoved.contains(neighborCoords[1][2]) ||
                    (!overridePointsAdded.contains(neighborCoords[1][2]) &&
                    !points.contains(neighborCoords[1][2])))
                    || 
                    (overridePointsRemoved.contains(neighborCoords[1][0]) ||
                    (!overridePointsAdded.contains(neighborCoords[1][0]) &&
                    !points.contains(neighborCoords[1][0])))
                    ) {
                    
                    return true;
                }
            }
            /*
             | ? ?
             | C .     '#' is disconnected if All .'s are 0s
             | _ #     AND AtLeastOne (? is 1)
             |
             | ? ?
             | C .     '#' is disconnected if All .'s are 0s
             | # _     AND AtLeastOne (? is 1)
             */
            if (
                (
                    (!overridePointsRemoved.contains(neighborCoords[1][0]) &&
                    (overridePointsAdded.contains(neighborCoords[1][0]) ||
                    points.contains(neighborCoords[1][0])))
                    || 
                    (!overridePointsRemoved.contains(neighborCoords[2][0]) &&
                    (overridePointsAdded.contains(neighborCoords[2][0]) ||
                    points.contains(neighborCoords[2][0])))
                ) 
                && 
                    (overridePointsRemoved.contains(neighborCoords[2][1]) ||
                    (!overridePointsAdded.contains(neighborCoords[2][1]) &&
                    !points.contains(neighborCoords[2][1])))
                ) {
                
                if (
                    (!overridePointsRemoved.contains(neighborCoords[1][2]) &&
                    (overridePointsAdded.contains(neighborCoords[1][2]) ||
                    points.contains(neighborCoords[1][2])))
                    || 
                    (!overridePointsRemoved.contains(neighborCoords[2][2]) &&
                    (overridePointsAdded.contains(neighborCoords[2][2]) ||
                    points.contains(neighborCoords[2][2])))
                    ) {
                    
                    return true;
                }
            }
            return false;
            
        } else if ((col + 1) > (imageWidth - 1)) {
            /*
             @ @ |
             @ C |
             @ @ |
             */
            // already tested that col-1 exists
            if ((row - 1) < 0) {
                /*
                 @ @ |
                 @ C |
                 */
                if ((row + 1) > (imageHeight - 1)) {
                    /*
                     @ C |
                     */
                    return true;
                }
                if (
                    (!overridePointsRemoved.contains(neighborCoords[0][1]) &&
                    (overridePointsAdded.contains(neighborCoords[0][1]) ||
                    points.contains(neighborCoords[0][1])))
                    || 
                    (!overridePointsRemoved.contains(neighborCoords[0][2]) &&
                    (overridePointsAdded.contains(neighborCoords[0][2]) ||
                    points.contains(neighborCoords[0][2])))
                    ||
                    (!overridePointsRemoved.contains(neighborCoords[1][2]) &&
                    (overridePointsAdded.contains(neighborCoords[1][2]) ||
                    points.contains(neighborCoords[1][2])))) {
                    
                    return true;
                }
                return false;
            } else if ((row + 1) > (imageHeight - 1)) {
                /*
                 @ C |
                 @ @ |
                 */
                if (
                    (!overridePointsRemoved.contains(neighborCoords[0][1]) &&
                    (overridePointsAdded.contains(neighborCoords[0][1]) ||
                    points.contains(neighborCoords[0][1])))
                    || 
                    (!overridePointsRemoved.contains(neighborCoords[0][0]) &&
                    (overridePointsAdded.contains(neighborCoords[0][0]) ||
                    points.contains(neighborCoords[0][0])))
                    ||
                    (!overridePointsRemoved.contains(neighborCoords[1][0]) &&
                    (overridePointsAdded.contains(neighborCoords[1][0]) ||
                    points.contains(neighborCoords[1][0])))) {
                    
                    return true;
                }
                return false;
            }
            /*
             @ @ |
             @ C |
             @ @ |
             */
            /*
             - # |
             . C | '#' is disconnected if All .'s are 0s
             ? ? |   AND AtLeastOne (? is 1)
             # _ |
             . C | '#' is disconnected if All .'s are 0s
             ? ? |   AND AtLeastOne (? is 1)
             */
            if (
                (
                    (!overridePointsRemoved.contains(neighborCoords[1][2]) &&
                    (overridePointsAdded.contains(neighborCoords[1][2]) ||
                    points.contains(neighborCoords[1][2])))
                    || 
                    (!overridePointsRemoved.contains(neighborCoords[0][2]) &&
                    (overridePointsAdded.contains(neighborCoords[0][2]) ||
                    points.contains(neighborCoords[0][2])))
                ) 
                && 
                (overridePointsRemoved.contains(neighborCoords[0][1]) ||
                (!overridePointsAdded.contains(neighborCoords[0][1]) &&
                !points.contains(neighborCoords[0][1])))
                ) {
                
                if (
                    (!overridePointsRemoved.contains(neighborCoords[1][0]) &&
                    (overridePointsAdded.contains(neighborCoords[1][0]) ||
                    points.contains(neighborCoords[1][0])))
                    || 
                    (!overridePointsRemoved.contains(neighborCoords[0][0]) &&
                    (overridePointsAdded.contains(neighborCoords[0][0]) ||
                    points.contains(neighborCoords[0][0])))
                    ) {
                    
                    return true;
                }
            }
            /*
             _ . |
             # C |  '#' is disconnected if All .'s are 0s
             _ . |
             */
            if (
                (!overridePointsRemoved.contains(neighborCoords[0][1]) &&
                (overridePointsAdded.contains(neighborCoords[0][1]) ||
                points.contains(neighborCoords[0][1])))
                ) {
                if (
                    (overridePointsRemoved.contains(neighborCoords[1][2]) ||
                    (!overridePointsAdded.contains(neighborCoords[1][2]) &&
                    !points.contains(neighborCoords[1][2])))
                    || 
                    (overridePointsRemoved.contains(neighborCoords[1][0]) ||
                    (!overridePointsAdded.contains(neighborCoords[1][0]) &&
                    !points.contains(neighborCoords[1][0])))
                    ) {
                    return true;
                }
            }
            /*
             ? ? |
             . C | '#' is disconnected if All .'s are 0s
             - # |   AND AtLeastOne (? is 1)
             ? ? |
             . C | '#' is disconnected if All .'s are 0s
             # _ |   AND AtLeastOne (? is 1)
             */
            if (
                (
                    (!overridePointsRemoved.contains(neighborCoords[1][0])
                    && (overridePointsAdded.contains(neighborCoords[1][0])
                    || points.contains(neighborCoords[1][0])))
                    || 
                    (!overridePointsRemoved.contains(neighborCoords[0][0])
                    && (overridePointsAdded.contains(neighborCoords[0][0])
                    || points.contains(neighborCoords[0][0])))
                ) 
                && 
                (
                    (overridePointsRemoved.contains(neighborCoords[0][1]) ||
                    (!overridePointsAdded.contains(neighborCoords[0][1]) &&
                    !points.contains(neighborCoords[0][1])))
                )
                ) {
                
                if (
                    (!overridePointsRemoved.contains(neighborCoords[1][2])
                    && (overridePointsAdded.contains(neighborCoords[1][2])
                    || points.contains(neighborCoords[1][2])))
                    || 
                    (!overridePointsRemoved.contains(neighborCoords[0][2])
                    && (overridePointsAdded.contains(neighborCoords[0][2])
                    || points.contains(neighborCoords[0][2])))
                    ) {
                    
                    return true;
                }
            }
            return false;
        } else if ((row - 1) < 0) {
            /*
             @ @ @
             @ C @
             */
            // have already tested that column 0 and column 2 are present
            /*
             ?  .  #    '#' is disconnected if All .'s are 0s
             ?  C  _        AND AtLeastOne (? is 1)
             ?  .  _    '#' is disconnected if All .'s are 0s
             ?  C  #        AND AtLeastOne (? is 1)
             */
            if (
                (
                    (!overridePointsRemoved.contains(neighborCoords[2][1]) &&
                    (overridePointsAdded.contains(neighborCoords[2][1]) ||
                    points.contains(neighborCoords[2][1])))
                    || 
                    (!overridePointsRemoved.contains(neighborCoords[2][2]) &&
                    (overridePointsAdded.contains(neighborCoords[2][2]) ||
                    points.contains(neighborCoords[2][2])))
                ) 
                && 
                (overridePointsRemoved.contains(neighborCoords[1][2]) ||
                (!overridePointsAdded.contains(neighborCoords[1][2]) &&
                !points.contains(neighborCoords[1][2])))
                ) {
                
                if (
                    (!overridePointsRemoved.contains(neighborCoords[0][2]) &&
                    (overridePointsAdded.contains(neighborCoords[0][2]) ||
                    points.contains(neighborCoords[0][2])))
                    || 
                    (!overridePointsRemoved.contains(neighborCoords[0][1]) &&
                    (overridePointsAdded.contains(neighborCoords[0][1]) ||
                    points.contains(neighborCoords[0][1])))
                    ) {
                    return true;
                }
            }
            /*
             -  #  _
             .  C  .    '#' is disconnected if All .'s are 0s
             */
            if (
                (!overridePointsRemoved.contains(neighborCoords[1][2]) &&
                (overridePointsAdded.contains(neighborCoords[1][2]) ||
                points.contains(neighborCoords[1][2])))
                ) {
                if (
                    (overridePointsRemoved.contains(neighborCoords[0][1]) ||
                    (!overridePointsAdded.contains(neighborCoords[0][1]) &&
                    !points.contains(neighborCoords[0][1])))
                    || 
                    (overridePointsRemoved.contains(neighborCoords[2][1]) ||
                    (!overridePointsAdded.contains(neighborCoords[2][1]) &&
                    !points.contains(neighborCoords[2][1])))
                    ) {
                    return true;
                }
            }
            /*
             #  .  ?    '#' is disconnected if All .'s are 0s
             _  C  ?       AND AtLeastOne (? is 1)
             _  .  ?    '#' is disconnected if All .'s are 0s
             #  C  ?        AND AtLeastOne (? is 1)
             */
            if (
                (
                    (!overridePointsRemoved.contains(neighborCoords[0][1]) &&
                    (overridePointsAdded.contains(neighborCoords[0][1]) ||
                    points.contains(neighborCoords[0][1])))
                    || 
                    (!overridePointsRemoved.contains(neighborCoords[0][2]) &&
                    (overridePointsAdded.contains(neighborCoords[0][2]) ||
                    points.contains(neighborCoords[0][2])))
                ) 
                && 
                    (overridePointsRemoved.contains(neighborCoords[1][2]) ||
                    (!overridePointsAdded.contains(neighborCoords[1][2]) &&
                    !points.contains(neighborCoords[1][2])))
                ) {
                
                if (
                    (!overridePointsRemoved.contains(neighborCoords[2][2]) &&
                    (overridePointsAdded.contains(neighborCoords[2][2]) ||
                    points.contains(neighborCoords[2][2])))
                    || 
                    (!overridePointsRemoved.contains(neighborCoords[2][1]) &&
                    (overridePointsAdded.contains(neighborCoords[2][1]) ||
                    points.contains(neighborCoords[2][1])))
                    ) {
                    
                    return true;
                }
            }
            return false;
        } else if ((row + 1) > (imageHeight - 1)) {
            /*
             @ C @
             @ @ @
             */
            // have already tested that column 0 and column 2 are present
            // and that the bottom row is present
            /*
             ?  C  _    '#' is disconnected if All .'s are 0s
             ?  .  #       AND AtLeastOne (? is 1)
             ?  C  #    '#' is disconnected if All .'s are 0s
             ?  .  _       AND AtLeastOne (? is 1)
             */
            if (
                (
                    (!overridePointsRemoved.contains(neighborCoords[2][0]) &&
                    (overridePointsAdded.contains(neighborCoords[2][0]) ||
                    points.contains(neighborCoords[2][0])))
                    || 
                    (!overridePointsRemoved.contains(neighborCoords[2][1]) &&
                    (overridePointsAdded.contains(neighborCoords[2][1]) ||
                    points.contains(neighborCoords[2][1])))
                ) && 
                    (overridePointsRemoved.contains(neighborCoords[1][0]) ||
                    (!overridePointsAdded.contains(neighborCoords[1][0]) &&
                    !points.contains(neighborCoords[1][0])))
                ) {
                
                if (
                    (!overridePointsRemoved.contains(neighborCoords[0][1]) &&
                    (overridePointsAdded.contains(neighborCoords[0][1]) ||
                    points.contains(neighborCoords[0][1])))
                    || 
                    (!overridePointsRemoved.contains(neighborCoords[0][0]) &&
                    (overridePointsAdded.contains(neighborCoords[0][0]) ||
                    points.contains(neighborCoords[0][0])))
                    ) {
                    
                    return true;
                }
            }
            /*
             .  C  .
             -  #  _  '#' is disconnected if All .'s are 0s
             */
            if (
                (!overridePointsRemoved.contains(neighborCoords[1][0]) &&
                (overridePointsAdded.contains(neighborCoords[1][0]) ||
                points.contains(neighborCoords[1][0])))
                ) {
                
                if (
                    (overridePointsRemoved.contains(neighborCoords[0][1]) ||
                    (!overridePointsAdded.contains(neighborCoords[0][1]) &&
                    !points.contains(neighborCoords[0][1])))
                    || 
                    (overridePointsRemoved.contains(neighborCoords[2][1]) ||
                    (!overridePointsAdded.contains(neighborCoords[2][1]) &&
                    !points.contains(neighborCoords[2][1])))
                    ) {
                    return true;
                }
            }
            /*
             -  C  ?    '#' is disconnected if All .'s are 0s
             #  .  ?       AND AtLeastOne (? is 1)
             #  C  ?    '#' is disconnected if All .'s are 0s
             _  .  ?      AND AtLeastOne (? is 1)
             */
            if (
                (
                    (!overridePointsRemoved.contains(neighborCoords[0][0]) &&
                    (overridePointsAdded.contains(neighborCoords[0][0]) ||
                    points.contains(neighborCoords[0][0])))
                    || 
                    (!overridePointsRemoved.contains(neighborCoords[0][1]) &&
                    (overridePointsAdded.contains(neighborCoords[0][1]) ||
                    points.contains(neighborCoords[0][1])))
                ) 
                && 
                    (overridePointsRemoved.contains(neighborCoords[1][0]) ||
                    (!overridePointsAdded.contains(neighborCoords[1][0]) &&
                    !points.contains(neighborCoords[1][0])))
                ) {
                
                if (
                    (!overridePointsRemoved.contains(neighborCoords[2][1]) &&
                    (overridePointsAdded.contains(neighborCoords[2][1]) ||
                    points.contains(neighborCoords[2][1])))
                    || 
                    (!overridePointsRemoved.contains(neighborCoords[2][0]) &&
                    (overridePointsAdded.contains(neighborCoords[2][0]) ||
                    points.contains(neighborCoords[2][0])))
                    ) {
                    return true;
                }
            }
            return false;
        }
        return false;
    }

    protected boolean isASmallerSearchRegion(int col, int row, int imageWidth,
        int imageHeight) {
        
        if ((col - 1) < 0) {
            return true;
        } else if ((col + 1) > (imageWidth - 1)) {
            return true;
        } else if ((row - 1) < 0) {
            return true;
        } else if ((row + 1) > (imageHeight - 1)) {
            return true;
        }
        return false;
    }
    
    public PairInt[][] createCoordinatePointsForEightNeighbors(
        int col, int row) {
        
        PairInt[][] set = new PairInt[3][];
        
        for (int i = 0; i < 3; ++i) {
            set[i] = new PairInt[3];
            int x = col + i - 1;
            for (int j = 0; j < 3; ++j) {
                int y = row + j - 1;
                set[i][j] = new PairInt(x, y);
            }
        }
        
        return set;
    }

    public void rotateBy90(PairInt[][] neighborCoords) {
        
        /*
            6   7  8     +1  2      transformed by 90 rot:     15  11  6
           11 *C* 12     0   1                                 16  C*  7
           15  16 17     -1  0                                 17  12  8
        
           -1  0   1
            0  1   2
         */
        
        PairInt tmp6 = neighborCoords[0][2];
        PairInt tmp7 = neighborCoords[1][2];
        PairInt tmp8 = neighborCoords[2][2];
        
        neighborCoords[0][2] = neighborCoords[0][0];
        neighborCoords[1][2] = neighborCoords[0][1];
        neighborCoords[2][2] = tmp6;
        
        PairInt tmp12 = neighborCoords[2][1];
        
        neighborCoords[0][1] = neighborCoords[1][0];
        neighborCoords[2][1] = tmp7;
        
        neighborCoords[0][0] = neighborCoords[2][0];
        neighborCoords[1][0] = tmp12;
        neighborCoords[2][0] = tmp8;
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
    protected boolean doesDisconnect(PairInt p, Set<PairInt> points, 
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
            6   7  8     +1  2      transformed by 90 rot:     15  11  6
           11 *C* 12     0   1                                 16  C*  7
           15  16 17     -1  0                                 17  12  8
        
           -1  0   1
            0  1   2
        
        disconnects:
           -- if (6 || 7 || 8) && (11 && 17) && !(15) && !(16)
           -- if (6 || 7 || 8) && (12 && 15) && !(16) && !(17)
           -- if (6 || 7 || 8) && (15 || 16 || 17) && !(11) && !(12)
           -- if !(6) && !(7) && (8) && (11) && !(12)
           -- if (6) && !(7) && !(8) && !(11) && (12)
        does not disconnect
           -- if (6 || 7 || 8) && !(15) && !(16) && !(17) && !(11) && !(12)
        
        then rotate 90 and test, then rotate 90 and test, then rotate 90 and test
        */
        
        for (int nRot = 0; nRot < 4; nRot++) {
        
            if (nRot > 0) {
                rotateBy90(neighborCoords);
            }
            
            if (//if (6 || 7 || 8)
            (!overridePointsRemoved.contains(neighborCoords[0][2]) &&
            (overridePointsAdded.contains(neighborCoords[0][2]) ||
            points.contains(neighborCoords[0][2])))
            ||
            (!overridePointsRemoved.contains(neighborCoords[1][2]) &&
            (overridePointsAdded.contains(neighborCoords[1][2]) ||
            points.contains(neighborCoords[1][2])))
            ||
            (!overridePointsRemoved.contains(neighborCoords[2][2]) &&
            (overridePointsAdded.contains(neighborCoords[2][2]) ||
            points.contains(neighborCoords[2][2])))
            ) {
                if ( //if (6 || 7 || 8) && (11 && 17) && !(15) && !(16)
                (
                    (!overridePointsRemoved.contains(neighborCoords[0][1]) &&
                    (overridePointsAdded.contains(neighborCoords[0][1]) ||
                    points.contains(neighborCoords[0][1])))
                    &&
                    (!overridePointsRemoved.contains(neighborCoords[2][0]) &&
                    (overridePointsAdded.contains(neighborCoords[2][0]) ||
                    points.contains(neighborCoords[2][0])))
                )
                &&
                (
                    !(!overridePointsRemoved.contains(neighborCoords[0][0]) &&
                    (overridePointsAdded.contains(neighborCoords[0][0]) ||
                    points.contains(neighborCoords[0][0])))
                    &&
                    !(!overridePointsRemoved.contains(neighborCoords[1][0]) &&
                    (overridePointsAdded.contains(neighborCoords[1][0]) ||
                    points.contains(neighborCoords[1][0])))
                )
                ) {
                    return true;
                } else if (//if (6 || 7 || 8) && (12 && 15) && !(16) && !(17)
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
                (
                    !(!overridePointsRemoved.contains(neighborCoords[1][0]) &&
                    (overridePointsAdded.contains(neighborCoords[1][0]) ||
                    points.contains(neighborCoords[1][0])))
                    &&
                    !(!overridePointsRemoved.contains(neighborCoords[2][0]) &&
                    (overridePointsAdded.contains(neighborCoords[2][0]) ||
                    points.contains(neighborCoords[2][0])))
                )
                ) {
                    return true;
                } else if (//if (6 || 7 || 8) && (15 || 16 || 17) && !(11) && !(12)
                (
                    (!overridePointsRemoved.contains(neighborCoords[0][0]) &&
                    (overridePointsAdded.contains(neighborCoords[0][0]) ||
                    points.contains(neighborCoords[0][0])))
                    ||
                    (!overridePointsRemoved.contains(neighborCoords[1][0]) &&
                    (overridePointsAdded.contains(neighborCoords[1][0]) ||
                    points.contains(neighborCoords[1][0])))
                    ||
                    (!overridePointsRemoved.contains(neighborCoords[2][0]) &&
                    (overridePointsAdded.contains(neighborCoords[2][0]) ||
                    points.contains(neighborCoords[2][0])))
                )
                &&
                (
                    !(!overridePointsRemoved.contains(neighborCoords[0][1]) &&
                    (overridePointsAdded.contains(neighborCoords[0][1]) ||
                    points.contains(neighborCoords[0][1])))
                    &&
                    !(!overridePointsRemoved.contains(neighborCoords[2][1]) &&
                    (overridePointsAdded.contains(neighborCoords[2][1]) ||
                    points.contains(neighborCoords[2][1])))
                )
                ) {
                    return true;
                } else if (//if (6 || 7 || 8) && !(15) && !(16) && !(17) && !(11) && !(12)
                !(!overridePointsRemoved.contains(neighborCoords[0][0]) &&
                (overridePointsAdded.contains(neighborCoords[0][0]) ||
                points.contains(neighborCoords[0][0])))
                &&
                !(!overridePointsRemoved.contains(neighborCoords[1][0]) &&
                (overridePointsAdded.contains(neighborCoords[1][0]) ||
                points.contains(neighborCoords[1][0])))
                &&
                !(!overridePointsRemoved.contains(neighborCoords[2][0]) &&
                (overridePointsAdded.contains(neighborCoords[2][0]) ||
                points.contains(neighborCoords[2][0])))
                &&
                !(!overridePointsRemoved.contains(neighborCoords[0][1]) &&
                (overridePointsAdded.contains(neighborCoords[0][1]) ||
                points.contains(neighborCoords[0][1])))
                &&
                !(!overridePointsRemoved.contains(neighborCoords[2][1]) &&
                (overridePointsAdded.contains(neighborCoords[2][1]) ||
                points.contains(neighborCoords[2][1])))
                ) {
                    return false;
                }
            } else if (//if !(6) && !(7) && (8) && (11) && !(12)
            !(!overridePointsRemoved.contains(neighborCoords[0][2]) &&
            (overridePointsAdded.contains(neighborCoords[0][2]) ||
            points.contains(neighborCoords[0][2])))
            &&
            !(!overridePointsRemoved.contains(neighborCoords[1][2]) &&
            (overridePointsAdded.contains(neighborCoords[1][2]) ||
            points.contains(neighborCoords[1][2])))
            &&
            (!overridePointsRemoved.contains(neighborCoords[2][2]) &&
            (overridePointsAdded.contains(neighborCoords[2][2]) ||
            points.contains(neighborCoords[2][2])))
            &&
            (!overridePointsRemoved.contains(neighborCoords[0][1]) &&
            (overridePointsAdded.contains(neighborCoords[0][1]) ||
            points.contains(neighborCoords[0][1])))
            &&
            !(!overridePointsRemoved.contains(neighborCoords[2][1]) &&
            (overridePointsAdded.contains(neighborCoords[2][1]) ||
            points.contains(neighborCoords[2][1])))
            ) {
                return true;
            } else if (//if (6) && !(7) && !(8) && !(11) && (12)
            (!overridePointsRemoved.contains(neighborCoords[0][2]) &&
            (overridePointsAdded.contains(neighborCoords[0][2]) ||
            points.contains(neighborCoords[0][2])))
            &&
            !(!overridePointsRemoved.contains(neighborCoords[1][2]) &&
            (overridePointsAdded.contains(neighborCoords[1][2]) ||
            points.contains(neighborCoords[1][2])))
            &&
            !(!overridePointsRemoved.contains(neighborCoords[2][2]) &&
            (overridePointsAdded.contains(neighborCoords[2][2]) ||
            points.contains(neighborCoords[2][2])))
            &&
            !(!overridePointsRemoved.contains(neighborCoords[0][1]) &&
            (overridePointsAdded.contains(neighborCoords[0][1]) ||
            points.contains(neighborCoords[0][1])))
            &&
            (!overridePointsRemoved.contains(neighborCoords[2][1]) &&
            (overridePointsAdded.contains(neighborCoords[2][1]) ||
            points.contains(neighborCoords[2][1])))
            ) {
                return true;
            }            
        }
        
        return false;
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

}
