package algorithms.imageProcessing;

/**
 *
 * @author nichole
 */
public abstract class AbstractLineThinner implements ILineThinner {

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

    protected boolean isASmallerSearchRegion(final GreyscaleImage input, 
        int col, int row) {
        
        if ((col - 1) < 0) {
            return true;
        } else if ((col + 1) > (input.getWidth() - 1)) {
            return true;
        } else if ((row - 1) < 0) {
            return true;
        } else if ((row + 1) > (input.getHeight() - 1)) {
            return true;
        }
        return false;
    }
    
}
