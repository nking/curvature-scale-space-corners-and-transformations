package algorithms.imageProcessing;

import algorithms.util.PairInt;
import java.util.Set;

/**
 *
 * @author nichole
 */
class UntraversableLobeRemover {

    public boolean applyFilter(Set<PairInt> borderPixels) {
        
        /*
        Example:
                         -3
                 #       -2
              -  #  -    -1    The center pattern can be found, then if each
           #  #  -  -     0    of the 3 has a neighbor that isn't adjacent
              -  #  -     1    to the others, it is a junction that cannot
                    #          be traversed.
        
        0 -1  0  1  2


        To further locate the untraversable section,

        */
        
        return false;
    }
    
}
