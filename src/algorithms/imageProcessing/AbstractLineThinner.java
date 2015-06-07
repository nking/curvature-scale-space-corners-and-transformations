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
    
    protected PairInt[][] createCoordinatePointsForEightNeighbors(
        int col, int row) {
        
        PairInt[][] set = new PairInt[3][];
        
        for (int i = 0; i < 3; i++) {
            set[i] = new PairInt[3];
            int x = col + i - 1;
            for (int j = 0; j < 3; j++) {
                if (i == 1 && j == 1) {
                    continue;
                }
                int y = row + j - 1;
                set[i][j] = new PairInt(x, y);
            }
        }
        
        return set;
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
    
}
