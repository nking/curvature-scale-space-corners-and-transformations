package algorithms.imageProcessing.features;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.PixelHelper;
import gnu.trove.iterator.TLongIterator;
import gnu.trove.set.TLongSet;
import java.util.Arrays;
import java.util.Collection;

/**
 *
 * @author nichole
 */
public class HOGUtil {
    
    private static float eps = 0.000001f;
    
    /**
     * CAVEAT: small amount of testing done, not yet throughly tested.
     *
     * calculate the intersection of histA and histB which have already
     * been normalized to the same scale.
     * A result of 0 is maximally dissimilar and a result of 1 is maximally similar.
     *
     * The orientations are needed to compare the correct rotated bins to one another.
     * Internally, orientation of 90 leads to no shift for rotation,
     * and orientation near 0 results in rotation of nBins/2, etc...
     *
     * Note that an orientation of 90 is a unit vector from x,y=0,0 to
     * x,y=0,1.
     * 
     * Note also that you may want to try the rotation of oppossite direction.
     *
     * @param histA
     * @param orientationA
     * @param histB
     * @param orientationB
     * @return
     */
    public static float intersection(int[] histA, int orientationA, int[] histB,
        int orientationB) {

        if ((histA.length != histB.length)) {
            throw new IllegalArgumentException(
                "histA and histB must be same dimensions");
        }

        if (orientationA < 0 || orientationA > 180 || orientationB < 0 ||
            orientationB > 180) {
            throw new IllegalArgumentException("orientations must be in range 0 to 180,"
                + "  inclusive,  or!=" + orientationA + " orB=" + orientationB);
        }
        if (orientationA == 180) {
            orientationA = 0;
        }
        if (orientationB == 180) {
            orientationB = 0;
        }

        int nBins = histA.length;

        int binWidth = 180/nBins;

        int shiftA = (orientationA - 90)/binWidth;
        int shiftB = (orientationB - 90)/binWidth;

        /*
        histograms are already normalized

        K(a,b) =
            (summation_over_i_from_1_to_n( min(a_i, b_i))
             /
            (min(summation_over_i(a_i), summation_over_i(b_i))
        */

        float sum = 0;
        float sumA = 0;
        float sumB = 0;
        for (int j = 0; j < nBins; ++j) {

            int idxA = j + shiftA;
            if (idxA < 0) {
                idxA += nBins;
            } else if (idxA > (nBins - 1 )) {
                idxA -= nBins;
            }

            int idxB = j + shiftB;
            if (idxB < 0) {
                idxB += nBins;
            } else if (idxB > (nBins - 1 )) {
                idxB -= nBins;
            }
            
            float yA = histA[idxA];
            float yB = histB[idxB];

            sum += Math.min(yA, yB);
            sumA += yA;
            sumB += yB;

            //System.out.println(" " + yA + " -- " + yB + " sum="+sum + ", " + sumA + "," + sumB);
        }

        float d = eps + Math.min(sumA, sumB);
        float sim = sum/d;

        //System.out.format("  (%d) hA=%s\n", orientationA, Arrays.toString(histA));
        //System.out.format("  (%d) hB=%s\n", orientationB, Arrays.toString(histB));
        //System.out.println("->inters=" + sim);
        
        return sim;
    }
    
    /**
     * 
     * calculate the intersection of histA and histB which have already
     * been normalized to the same scale.
     * A result of 0 is maximally dissimilar and a result of 1 is maximally similar.
     * 
     * The orientations are needed to compare the correct rotated bins to one another.
     * Internally, orientation of 90 leads to no shift for rotation,
     * and orientation near 0 results in rotation of nBins/2, etc...
     * 
     * Note that an orientation of 90 is a unit vector from x,y=0,0 to
     * x,y=0,1.
     * 
     * @param histA
     * @param histB
     * @return 
     */
    public static float intersection(int[] histA, int[] histB) {
        
        //return HOGUtil.intersection(histA, orientationA, histB, orientationB);
    
        if ((histA.length != histB.length)) {
            throw new IllegalArgumentException(
                "histA and histB must be same dimensions");
        }
        
        int nBins = histA.length;
        
        int binWidth = 256/nBins;
        
        /*
        histograms are already normalized
        
        K(a,b) = 
            (summation_over_i_from_1_to_n( min(a_i, b_i))
             /
            (min(summation_over_i(a_i), summation_over_i(b_i))
        */
                
        float sum = 0;
        float sumA = 0;
        float sumB = 0;
        for (int j = 0; j < nBins; ++j) {
            
            float yA = histA[j];
            float yB = histB[j];
            
            sum += Math.min(yA, yB);
            sumA += yA;
            sumB += yB;
            
            //System.out.println(" " + yA + " -- " + yB + " sum="+sum + ", " + sumA + "," + sumB);
        }
        
        float d = eps +  Math.min(sumA, sumB);
        
        float sim = sum/d;
        
        return sim;
    }
    
    /**
     * CAVEAT: small amount of testing done, not yet throughly tested.
     * 
     * calculate the intersection of histA and histB which have already
     * been normalized to the same scale.
     * A result of 0 is maximally dissimilar and a result of 1 is maximally similar.
     * 
     The orientations are needed to compare the correct rotated bins to one another.
     * Internally, orientation of 90 leads to no shift for rotation,
     * and orientation near 0 results in rotation of nBins/2, etc...
     *
     * Note that an orientation of 90 is a unit vector from x,y=0,0 to
     * x,y=0,1.
     *
     * @param histA
     * @param orientationA
     * @param histB
     * @param orientationB
     * @return difference, error
     */
    public static float[] diff(int[] histA, int orientationA, int[] histB,
        int orientationB) {

        if ((histA.length != histB.length)) {
            throw new IllegalArgumentException(
                "histA and histB must be same dimensions");
        }

        if (orientationA < 0 || orientationA > 180 || orientationB < 0 ||
            orientationB > 180) {
            throw new IllegalArgumentException("orientations must be in range 0 to 180,"
                + "  inclusive,  or!=" + orientationA + " orB=" + orientationB);
        }
        if (orientationA == 180) {
            orientationA = 0;
        }
        if (orientationB == 180) {
            orientationB = 0;
        }

        int nBins = histA.length;

        int binWidth = 180/nBins;

        int shiftA = (orientationA - 90)/binWidth;
        int shiftB = (orientationB - 90)/binWidth;

        /*
        histograms are already normalized

        K(a,b) =
            (summation_over_i_from_1_to_n( min(a_i, b_i))
             /
            (min(summation_over_i(a_i), summation_over_i(b_i))
        */

        double sumDiff = 0;
        
        double err = 0;
        for (int j = 0; j < nBins; ++j) {

            int idxA = j + shiftA;
            if (idxA < 0) {
                idxA += nBins;
            } else if (idxA > (nBins - 1 )) {
                idxA -= nBins;
            }

            int idxB = j + shiftB;
            if (idxB < 0) {
                idxB += nBins;
            } else if (idxB > (nBins - 1 )) {
                idxB -= nBins;
            }

            float yA = histA[idxA];
            float yB = histB[idxB];

            float maxValue = Math.max(yA, yB) + eps;

            float diff = Math.abs((yA - yB)/maxValue);
            
            //sumDiff += (diff * diff);
            sumDiff += diff;

            //      already squared
            err += (diff/maxValue);
        }
        
        sumDiff /= (double)nBins;
        
        //sumDiff = Math.sqrt(sumDiff);
        
        err /= (double)nBins;
        err = Math.sqrt(err);
        
        return new float[]{(float)sumDiff, (float)err};
    }
    
    /**
     * CAVEAT: small amount of testing done, not yet throughly tested.
     * 
     * calculate the intersection of histA and histB which have already
     * been normalized to the same scale.
     * A result of 0 is maximally dissimilar and a result of 1 is maximally similar.
     * 
     * The orientations are needed to compare the correct rotated bins to one another.
     * Internally, orientation of 90 leads to no shift for rotation,
     * and orientation near 0 results in rotation of nBins/2, etc...
     * 
     * Note that an orientation of 90 is a unit vector from x,y=0,0 to
     * x,y=0,1.
     * 
     * @param histA
     * @param histB
     * @return 
     */
    public static float[] diff(int[] histA, int[] histB) {

        if ((histA.length != histB.length)) {
            throw new IllegalArgumentException(
                "histA and histB must be same dimensions");
        }

        int nBins = histA.length;

        int binWidth = 180/nBins;

        double sumDiff = 0;
        double err = 0;
                        
        for (int j = 0; j < nBins; ++j) {
            
            float yA = histA[j];
            float yB = histB[j];
            
            float maxValue = Math.max(yA, yB) + eps;

            float diff = Math.abs((yA - yB)/maxValue);
            
            //sumDiff += (diff * diff);
            sumDiff += diff;

            //      already squared
            err += (diff/maxValue);           
        }
        
        sumDiff /= (double)nBins;

        //sumDiff = Math.sqrt(sumDiff);

        err /= (double)nBins;
        err = Math.sqrt(err);
        
        return new float[]{(float)sumDiff, (float)err};
    }
    
    /**
     * 
     * @param img
     * @param points
     * @param outputMinMaxXY  populated from bounds of points
     * @param outputRefFramePixs populated for subimage referece frame
     * @return 
     */
    public static GreyscaleImage createAndMaskSubImage(GreyscaleImage img, 
        Collection<PairInt> points,
        int[] outputMinMaxXY, TLongSet outputRefFramePixs) {
        
        int maskValue = 0;
        
        return createAndMaskSubImage(img, maskValue, points, 
            outputMinMaxXY, outputRefFramePixs);
    }
    
    /**
     * 
     * @param img
     * @param maskValue
     * @param points
     * @param outputMinMaxXY  populated from bounds of points
     * @param outputRefFramePixs populated for subimage referece frame
     * @return 
     */
    public static GreyscaleImage createAndMaskSubImage(GreyscaleImage img, 
        int maskValue, Collection<PairInt> points,
        int[] outputMinMaxXY, TLongSet outputRefFramePixs) {

        PixelHelper ph = new PixelHelper();
                
        outputMinMaxXY[0] = Integer.MAX_VALUE;
        outputMinMaxXY[1] = Integer.MIN_VALUE;
        outputMinMaxXY[2] = Integer.MAX_VALUE;
        outputMinMaxXY[3] = Integer.MIN_VALUE;
        for (PairInt xy : points) {
            if (xy.getX() < outputMinMaxXY[0]) {
                outputMinMaxXY[0] = xy.getX();
            }
            if (xy.getX() > outputMinMaxXY[1]) {
                outputMinMaxXY[1] = xy.getX();
            }
            if (xy.getY() < outputMinMaxXY[2]) {
                outputMinMaxXY[2] = xy.getY();
            }
            if (xy.getY() > outputMinMaxXY[3]) {
                outputMinMaxXY[3] = xy.getY();
            }
        }
        
        GreyscaleImage img2 = img.subImage2(outputMinMaxXY[0], 
            outputMinMaxXY[1], outputMinMaxXY[2], outputMinMaxXY[3]);
                
        int w0 = img.getWidth();
        int h0 = img.getHeight();
        int w2 = img2.getWidth();
        int h2 = img2.getHeight();
        int xOffset = outputMinMaxXY[0];
        int yOffset = outputMinMaxXY[2];

        //In reference frame of subImage
        for (PairInt xy : points) {
            long pixIdx = ph.toPixelIndex(xy.getX() - xOffset, 
                xy.getY() - yOffset, w2);
            
            assert(pixIdx >= 0);
            assert(pixIdx < (w2 * h2));
            
            outputRefFramePixs.add(pixIdx);
        }
        
        // mask out pixels not in the region
        int c = 0;
        for (int i2 = 0; i2 < w2; ++i2) {
            for (int j2 = 0; j2 < h2; ++j2) {
                long pixIdx = ph.toPixelIndex(i2, j2, w2);
                if (!outputRefFramePixs.contains(pixIdx)) {
                    img2.setValue(i2, j2, maskValue);
                    c++;
                }
            }
        }
        System.out.println("  masked " + c + " out of " + (w2*h2));
        
        //DEBUG
        {
            int ts = MiscDebug.getCurrentTimeFormatted();
            Image img0 = img.copyToColorGreyscale();
            Image imgSub = img2.copyToColorGreyscale();
            
            int[] xy2 = new int[2];            
            for (PairInt p : points) {
                ImageIOHelper.addPointToImage(p.getX(), p.getY(), img0,
                    0, 255, 0, 0);
            }
            TLongIterator iter3 = outputRefFramePixs.iterator();
            while (iter3.hasNext()) {
                long pixIdx = iter3.next();
                ph.toPixelCoords(pixIdx, w2, xy2);
                ImageIOHelper.addPointToImage(xy2[0], xy2[1], imgSub,
                    0, 255, 0, 0);
            };
            MiscDebug.writeImage(img0, "_DGB_IMG__" + ts);
            MiscDebug.writeImage(imgSub, "_DGB_IMGSUB__" + ts);
        }
        
        return img2;
    }
    
    /**
     * 
     * @param img
     * @param points
     * @param minMaxXY given the bounds of points
     * @param refFramePixs given the coordinates already clipped to the sub-image
     *    that will be returned.
     * @return 
     */
    public static GreyscaleImage createAndMaskSubImage2(GreyscaleImage img, 
        int[] minMaxXY, TLongSet refFramePixs) {
        
        int maskValue = 0;
        
        return createAndMaskSubImage2(img, maskValue, minMaxXY, 
            refFramePixs);
    }

    /**
     * 
     * @param img
     * @param maskValue
     * @param points
     * @param minMaxXY given the bounds of points
     * @param refFramePixs given the coordinates already clipped to the sub-image
     *    that will be returned.
     * @return 
     */
    public static GreyscaleImage createAndMaskSubImage2(GreyscaleImage img, 
        int maskValue, int[] minMaxXY, TLongSet refFramePixs) {

        PixelHelper ph = new PixelHelper();
                
        GreyscaleImage img2 = img.subImage2(minMaxXY[0], minMaxXY[1], 
            minMaxXY[2], minMaxXY[3]);
                
        int w2 = img2.getWidth();
        int h2 = img2.getHeight();
        //int xOffset = minMaxXY[0];
        //int yOffset = minMaxXY[2];
        
        // mask out pixels not in the region
        for (int i2 = 0; i2 < w2; ++i2) {
            for (int j2 = 0; j2 < h2; ++j2) {
                long pixIdx = ph.toPixelIndex(i2, j2, w2);
                if (!refFramePixs.contains(pixIdx)) {
                    img2.setValue(i2, j2, maskValue);
                }
            }
        }
        
        return img2;
    }
}
