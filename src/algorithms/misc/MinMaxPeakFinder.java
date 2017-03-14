package algorithms.misc;

import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import java.util.Arrays;

/**
 * find peaks in sequential data by finding the minima above  threshold and
 * the maxima and then the maxima which are a gain factor above either
 * adjacent minima.
 * The success of the method depends upon a reasonable lowThreshold and
 * gain.
 * 
 * NOTE: methods such as those in MedianSmooth could be used to find a
 * windowed mean and max (and hence, standard deviation, but they are
 * dependent upon the size of the window.  A peak that is a wide gradual
 * plateau and a narrow window might miss the peak.
 * 
 * @author nichole
 */
public class MinMaxPeakFinder {
    
    public int[] findPeaks(float[] values) {
        
        float[] a = Arrays.copyOf(values, values.length);
        
        Arrays.sort(a);
        
        float mean3Percent = 0;
        int end = Math.round(0.03f * values.length);
        for (int i = 0; i < end; ++i) {
            mean3Percent += a[i];
        }
        mean3Percent /= (float)end;
        
        return findPeaks(values, mean3Percent, 2.5f);
    }
    
    public int[] findPeaks(float[] values, float lowThreshold, 
        float factorAboveMin) {
        
        int[] minMaxIdxs = findMinimaMaxima(values, lowThreshold);
        
        if (minMaxIdxs.length == 0) {
            return null;
        } else if (minMaxIdxs.length == 1) {
            return minMaxIdxs;
        } 
        
        TIntList peaks = new TIntArrayList(minMaxIdxs.length/2);

        // choose candidates from minMaxIndexes that are
        //     >= factorAboveMin for one adjacent minima
        for (int ii = 0; ii < minMaxIdxs.length; ii++) {

            int idx = minMaxIdxs[ii];

            if (idx > -1) {
                // this is a maxima

                boolean found = false;
                
                // compare to preceding minimum
                for (int iii = (ii - 1); iii > -1; iii--) {
                    int idx2 = minMaxIdxs[iii];
                    if (idx2 < 0) {
                        float compare = values[-1*idx2];
                        if (compare < lowThreshold) {
                            // avoids divide by very small number sometimes
                            compare = lowThreshold;
                        }
                        if (values[idx] >= lowThreshold && 
                            values[idx] >= factorAboveMin * compare) {
                            
                            peaks.add(idx);
                            found = true;
                        }
                        break;
                    }
                }
                
                if (found && (ii > 0)) {
                    continue;
                }

                //compare to proceeding minimum
                for (int iii = (ii + 1); iii < minMaxIdxs.length; iii++) {
                    int idx2 = minMaxIdxs[iii];
                    if (idx2 < 0) {
                        float compare = values[-1*idx2];
                        if (compare < lowThreshold) {
                            // avoids divide by very small number sometimes
                            compare = lowThreshold;
                        }
                        if (values[idx] >= lowThreshold 
                            && values[idx] >= factorAboveMin * compare) {
                            
                            peaks.add(idx);
                        }
                        
                        break;
                    }
                }
            }
        }

        return peaks.toArray(new int[peaks.size()]);
    }
    
    /**
     * find the minima above lowThreshold and the maxima in values
     * and return their indexes.  the negative values are -1*index for a minima
     * while positive indexes are the indexes of maxima.
     * 
     * @param values
     * @param lowThreshold
     * @return 
     */
    public int[] findMinimaMaxima(float[] values, float lowThreshold) {
        
        TIntList minMaxIdxs = new TIntArrayList();
        
        float lastK = values[0];
        boolean incr = true;
        for (int ii = 1; ii < values.length; ii++) {

            float currentK = values[ii];

            if ((currentK < lastK) && incr) {
                if (values[ii - 1] > lowThreshold) {
                    minMaxIdxs.add(ii - 1);
                }
                incr = false;
            } else if ((currentK > lastK) && !incr) {
                // values below outputLowThreshold[0] are handled by
                // callers.  TODO: redesign the caller and this method
                // to not need to understand peculiarities of the data.
                minMaxIdxs.add(-1*(ii - 1));
                incr = true;
            }

            lastK = currentK;
        }

        if (incr) {
            // add the last point
             minMaxIdxs.add(values.length - 1);
        }
        
        return minMaxIdxs.toArray(new int[minMaxIdxs.size()]);
    }
}
