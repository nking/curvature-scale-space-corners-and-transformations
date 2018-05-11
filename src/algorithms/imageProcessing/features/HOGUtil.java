package algorithms.imageProcessing.features;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.IntegralHistograms;
import algorithms.imageProcessing.SummedAreaTable;
import algorithms.util.PairInt;
import algorithms.util.PixelHelper;
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
     * create a 2D integral histogram image, apply a windows sum of cell size
     * to it, then transform back into a 2D histogram integral image again.
     * @param gXY
     * @param theta
     * @param nAngleBins
     * @param nPixPerCellDimension
     * @return 
     */
    public static int[][] createHOGHistogram(GreyscaleImage gXY, 
        GreyscaleImage theta, int nAngleBins, int nPixPerCellDimension) {
        
        GradientIntegralHistograms gh = new GradientIntegralHistograms();
        
        int[][] histograms = gh.createHistograms(gXY, theta, nAngleBins);

        //apply a windowed sum across the integral image.
        // result is that at each pixel is a histogram holding the sum of histograms
        //    from the surrounding N_PIX_PER_CELL_DIM window. 
        // The result is a 2D histogram integral image
        gh.applyWindowedSum(histograms, gXY.getWidth(), gXY.getHeight(), 
            nPixPerCellDimension);

        return histograms;
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
    
    public static int[][] createHCPTHistogram(GreyscaleImage ptImg, 
        TLongSet regionPixelCoords, int nHistBins, int nPixPerCellDimension) {

        int w2 = ptImg.getWidth();
        int h2 = ptImg.getHeight();
        
        PolarThetaIntegralHistograms gh = new PolarThetaIntegralHistograms();
        
        int[][] histograms = gh.createHistograms(ptImg, 
            regionPixelCoords, nHistBins);

        //apply a windowed avg across the integral image
        gh.applyWindowedSum(histograms, w2, h2, nPixPerCellDimension);
        
        return histograms;
    }
    
    public static int[][] createHGSHistogram(GreyscaleImage gsImg, 
        TLongSet regionPixelCoords, int nHistBins, int nPixPerCellDimension) {

        int w2 = gsImg.getWidth();
        int h2 = gsImg.getHeight();
        
        IntegralHistograms gh = new IntegralHistograms();
            
        int[][] histogramsHGS = gh.create(gsImg, regionPixelCoords, 
            0, 255, nHistBins);
        //apply a windowed avg across the integral image
        gh.applyWindowedSum(histogramsHGS, w2, h2, nPixPerCellDimension);
        
        return histogramsHGS;
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
        //System.out.println("  masked " + c + " out of " + (w2*h2));
        
        //DEBUG
        /*{
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
        }*/
        
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
    
    /**
     * 
     * runtime complexity is O(nBins)
     * 
     * @param histograms the 2D histogram integral image.  the first dimension 
     *    is the pixel index and the 2nd is the histogram bin, e.g.
     *    histograms[pixIdx][binIdx].
     * @param startX
     * @param stopX the last x pixel in the window, inclusive
     * @param startY
     * @param stopY the last y pixel in the window, inclusive
     * @param w image width
     * @param h image height
     * @param output
     * @param outputN an empty 1 dimensional array of size 1 to return the 
     * number of pixels in the cell
     */
    public static void extractWindow(int[][] histograms, int startX, int stopX, 
        int startY, int stopY, int w, int h, 
        int output[], int[] outputN) {

        if (output.length != histograms[0].length) {
            throw new IllegalArgumentException("output.length must == nThetaBins");
        }
        
        if (stopX < startX || stopY < startY) {
            throw new IllegalArgumentException("stopX must be >= startX and "
                + "stopY >= startY");
        }
        
        PixelHelper ph = new PixelHelper();
        
        if (startX == 0 && startY == 0) {
            if (stopX == startX && stopY == startY) {
                outputN[0] = 1;
                System.arraycopy(histograms[(int)ph.toPixelIndex(stopX, stopY, w)], 0, 
                    output, 0, output.length);
            } else if (stopX > startX && stopY > startY) {
                outputN[0] = (stopX + 1) * (stopY + 1);
                System.arraycopy(histograms[(int)ph.toPixelIndex(stopX, stopY, w)], 0, 
                    output, 0, output.length);
            } else if (stopX > startX) {
                //startY==0 && stopY=0
                outputN[0] = (stopX + 1);
                System.arraycopy(histograms[(int)ph.toPixelIndex(stopX, stopY, w)], 0,
                    output, 0, output.length);
            } else if (stopY > startY) {
                outputN[0] = (stopY + 1);
                System.arraycopy(histograms[(int)ph.toPixelIndex(stopX, stopY, w)], 0,
                    output, 0, output.length);
            }
        } else if (startX > 0 && startY > 0) {
            outputN[0] = ((stopX - startX) + 1) * ((stopY - startY) + 1);

            System.arraycopy(histograms[(int)ph.toPixelIndex(stopX, stopY, w)], 0,
                output, 0, output.length);
            subtract(output, histograms[(int)ph.toPixelIndex(startX - 1, stopY, w)]);
            subtract(output, histograms[(int)ph.toPixelIndex(stopX, startY - 1, w)]);
            add(output, histograms[(int)ph.toPixelIndex(startX - 1, startY - 1, w)]);
                
        } else if (startX > 0) {
            //startY == 0
            if (stopX == startX && stopY == startY) {
                outputN[0] = 1;
                System.arraycopy(histograms[(int)ph.toPixelIndex(stopX, stopY, w)], 0, 
                    output, 0, output.length);
            } else {
                outputN[0] = ((stopX - startX) + 1) * ((stopY - startY) + 1);

                System.arraycopy(histograms[(int)ph.toPixelIndex(stopX, stopY, w)], 0,
                    output, 0, output.length);
                subtract(output, histograms[(int)ph.toPixelIndex(startX - 1, stopY, w)]);
            }       
        } else if (startY > 0) {
            //startX == 0
            if (stopX == startX && stopY == startY) {
                outputN[0] = 1;
                System.arraycopy(histograms[(int)ph.toPixelIndex(stopX, stopY, w)], 0, 
                    output, 0, output.length);
            } else {
                outputN[0] = ((stopX - startX) + 1) * ((stopY - startY) + 1);

                System.arraycopy(histograms[(int)ph.toPixelIndex(stopX, stopY, w)], 0,
                    output, 0, output.length);
                subtract(output, histograms[(int)ph.toPixelIndex(stopX, startY - 1, w)]);
            }   
        }
    }
    
    /**
     * extract the sum of histograms in the window inclusively defined as
     * (startX:stopX, startY:stopY).
     * 
     * runtime complexity is O(nBins)
     * 
     * @param histograms the 2D histogram integral image.  the first dimension 
     *    is the pixel index and the 2nd is the histogram bin, e.g.
     *    histograms[pixIdx][binIdx].
     * @param startX
     * @param stopX the last x pixel in the window, inclusive
     * @param startY
     * @param stopY the last y pixel in the window, inclusive
     * @param output
     * @param outputN an empty 1 dimensional array of size 1 to return the 
     * number of pixels in the cell
     */
    public static void extractWindow(int[][] histograms, int startX, int stopX, 
        int startY, int stopY, int w, int h, 
        long[] output, int[] outputN) {

        if (output.length != histograms[0].length) {
            throw new IllegalArgumentException("output.length must == nThetaBins");
        }
        
        if (stopX < startX || stopY < startY) {
            throw new IllegalArgumentException("stopX must be >= startX and "
                + "stopY >= startY");
        }
        
        PixelHelper ph = new PixelHelper();
        
        if (startX == 0 && startY == 0) {
            if (stopX == startX && stopY == startY) {
                outputN[0] = 1;
                copy(histograms[(int)ph.toPixelIndex(stopX, stopY, w)], 
                    output);
            } else if (stopX > startX && stopY > startY) {
                outputN[0] = (stopX + 1) * (stopY + 1);
                copy(histograms[(int)ph.toPixelIndex(stopX, stopY, w)], 
                    output);
            } else if (stopX > startX) {
                //startY==0 && stopY=0
                outputN[0] = (stopX + 1);
                copy(histograms[(int)ph.toPixelIndex(stopX, stopY, w)],
                    output);
            } else if (stopY > startY) {
                outputN[0] = (stopY + 1);
                copy(histograms[(int)ph.toPixelIndex(stopX, stopY, w)],
                    output);
            }
        } else if (startX > 0 && startY > 0) {
            outputN[0] = ((stopX - startX) + 1) * ((stopY - startY) + 1);

            copy(histograms[(int)ph.toPixelIndex(stopX, stopY, w)],
                output);
            subtract(output, histograms[(int)ph.toPixelIndex(startX - 1, stopY, w)]);
            subtract(output, histograms[(int)ph.toPixelIndex(stopX, startY - 1, w)]);
            add(output, histograms[(int)ph.toPixelIndex(startX - 1, startY - 1, w)]);
                
        } else if (startX > 0) {
            //startY == 0
            if (stopX == startX && stopY == startY) {
                outputN[0] = 1;
                copy(histograms[(int)ph.toPixelIndex(stopX, stopY, w)], 
                    output);
            } else {
                outputN[0] = ((stopX - startX) + 1) * ((stopY - startY) + 1);

                copy(histograms[(int)ph.toPixelIndex(stopX, stopY, w)],
                    output);
                subtract(output, histograms[(int)ph.toPixelIndex(startX - 1, stopY, w)]);
            }       
        } else if (startY > 0) {
            //startX == 0
            if (stopX == startX && stopY == startY) {
                outputN[0] = 1;
                copy(histograms[(int)ph.toPixelIndex(stopX, stopY, w)],
                    output);
            } else {
                outputN[0] = ((stopX - startX) + 1) * ((stopY - startY) + 1);

                copy(histograms[(int)ph.toPixelIndex(stopX, stopY, w)],
                    output);
                subtract(output, histograms[(int)ph.toPixelIndex(stopX, startY - 1, w)]);
            }   
        }
    }
    
    /**
     * apply a windowed sum across the 2D integral histogram image,
     * where the window size is N_PIX_PER_CELL_DIM.
     * The result is an image of histograms, where each histogram represents
     * the integration over the surrounding N_PIX_PER_CELL_DIM window.
     * For windows near the image edge, a factor is applied to bring the
     * counts up by factor (N_PIX_PER_CELL_DIM/n_pix_in_window).
     * The 2D histograms array is then made into a 2D integral histogram image
     * again.
     *  
     * @param histograms the 2D histogram integral image.  the first dimension 
     *    is the pixel index and the 2nd is the histogram bin, e.g.
     *    histograms[pixIdx][binIdx]
     * @param w
     * @param h
     * @param N_PIX_PER_CELL_DIM
     */    
    public static void applyWindowedSum(int[][] histograms, int w, int h, 
        int N_PIX_PER_CELL_DIM) {
        
        int[][] img2 = applyWindowedSum0(histograms, w, h, N_PIX_PER_CELL_DIM);
        
        img2 = transformIntoIntegral2DHist(img2, w, h);
        
        for (int i = 0; i < histograms.length; ++i) {
            System.arraycopy(img2[i], 0, histograms[i], 0, img2[i].length);
        }
    }
    
    /**
     * apply a windowed sum across the gradient integral histogram image,
     * where the window size is N_PIX_PER_CELL_DIM.
     * The result is an image of histograms, where each histogram represents
     * the integration over the surrounding N_PIX_PER_CELL_DIM window.
     * The result is NOT a 2D integral histogram image, just a 2D histogram image;
     *  
     * @param histograms the 2D histogram integral image.  the first dimension 
     *    is the pixel index and the 2nd is the histogram bin, e.g.
     *    histograms[pixIdx][binIdx]
     * @param w
     * @param h
     * @param N_PIX_PER_CELL_DIM
     * @return the 2D histogram image smoothed over a window of
     *     N_PIX_PER_CELL_DIM size.  Note that this is not an integral image.
     */    
    static int[][] applyWindowedSum0(int[][] histograms, int w, int h, 
        int N_PIX_PER_CELL_DIM) {
        
        if (N_PIX_PER_CELL_DIM < 1) {
            throw new IllegalArgumentException("N_PIX_PER_CELL_DIM must be >= 1");
        }
        
        int nBins = histograms[0].length;
                
        int[][] img2 = new int[w * h][];
        
        int[] outN = new int[1];
        
        int windowSize = N_PIX_PER_CELL_DIM * N_PIX_PER_CELL_DIM;
        
        // a centered window sum
        int r = N_PIX_PER_CELL_DIM >> 1;
        int r0, r1;
        if (r == 0) {
            r0 = 0;
            r1 = 0;
        } else if ((N_PIX_PER_CELL_DIM & 1) == 1) {
            r0 = -r;
            r1 = r;
        } else {
            r0 = -r;
            r1 = r - 1;
        }
        
        float factor;
                        
        // extract the summed area of each dxd window centered on x,y
        for (int x = 0; x < w; ++x) {
            
            int x2 = x + r0;
            int x3 = x + r1;
            if (x3 < 0) {
                continue;
            } else if (x2 < 0) {
                x2 = 0;
            } else if (x2 >= w) {
                break;
            }
            if (x3 >= w) {
                x3 = w - 1;
            }

            for (int y = 0; y < h; ++y) {
                
                int y2 = y + r0;
                int y3 = y + r1;
                if (y3 < 0) {
                    continue;
                } else if (y2 < 0) {
                    y2 = 0;
                } else if (y2 >= h) {
                    break;
                }
                if (y3 >= h) {
                    y3 = h - 1;
                }
                                
                int pixIdx = (y * w) + x;
                
                img2[pixIdx] = new int[nBins];

                extractWindow(histograms, x2, x3, y2, y3, w, h, img2[pixIdx], outN);
                
                if (outN[0] < windowSize) {
                    factor = (float)windowSize/(float)outN[0];
                    HOGUtil.mult(img2[pixIdx], factor);
                }             
            }
        }
        
        return img2;
    }
    
    /**
     * 
     * @param hist a 2D histogram image of format hist[pixIdx][hitBinIdx].
     * @param imageWidth
     * @param imageHeight
     * @return a 2D histogram integral image of format hist[pixIdx][hitBinIdx].
     */
    public static int[][] transformIntoIntegral2DHist(int[][] hist, int imageWidth,
        int imageHeight) {
        
        int nPix = imageWidth * imageHeight;
        PixelHelper ph = new PixelHelper();
        
        int[][] out = new int[nPix][];
        for (int i = 0; i < nPix; ++i) {
            out[i] = Arrays.copyOf(hist[i], hist[i].length);
        }
        
        int[] tmp;
        int pixIdx;
                
        for (int x = 0; x < imageWidth; ++x) {
            for (int y = 0; y < imageHeight; ++y) {
                
                pixIdx = (int)ph.toPixelIndex(x, y, imageWidth);
                
                tmp = out[pixIdx];
                                
                if (x > 0 && y > 0) {
                    HOGUtil.add(tmp, out[(int)ph.toPixelIndex(x - 1, y, imageWidth)]);
                    HOGUtil.add(tmp, out[(int)ph.toPixelIndex(x, y - 1, imageWidth)]);
                    HOGUtil.subtract(tmp, out[(int)ph.toPixelIndex(x - 1, y - 1, 
                        imageWidth)]);
                } else if (x > 0) {
                    HOGUtil.add(tmp, out[(int)ph.toPixelIndex(x - 1, y, imageWidth)]);
                    
                } else if (y > 0) {
                    HOGUtil.add(tmp, out[(int)ph.toPixelIndex(x, y - 1, imageWidth)]);
                    
                }
            }
        }
        return out;
    }
    
    public static void add(int[] addTo, int[] addFrom) {
        for (int i = 0; i < addTo.length; ++i) {
            addTo[i] += addFrom[i];
        }
    }
    
    public static void mult(int[] a, float factor) {
        for (int i = 0; i < a.length; ++i) {
            a[i] = Math.round(a[i] * factor);
        }
    }
    
    public static void add(long[] addTo, int[] addFrom) {
        for (int i = 0; i < addTo.length; ++i) {
            addTo[i] += addFrom[i];
        }
    }
    
    public static void subtract(int[] subtractFrom, int[] subtract) {
        for (int i = 0; i < subtractFrom.length; ++i) {
            subtractFrom[i] -= subtract[i];
        }
    }
    
    public static void subtract(long[] subtractFrom, int[] subtract) {
        for (int i = 0; i < subtractFrom.length; ++i) {
            subtractFrom[i] -= subtract[i];
        }
    }
    
    public static void divide(int[] a, int d) {
        for (int i = 0; i < a.length; ++i) {
            a[i] /= d;
        }
    }
    
    public static void copy(int[] src, long[] dest) {
        for (int i = 0; i < dest.length; ++i) {
            dest[i] = src[i];
        }
    }
   
    public static GreyscaleImage applyMeanSmoothing(GreyscaleImage img,
        int window) {
        
        SummedAreaTable sat = new SummedAreaTable();
        
        GreyscaleImage imgS = sat.createAbsoluteSummedAreaTable(img);
        
        imgS = sat.applyMeanOfWindowFromSummedAreaTable(imgS, window);
        
        return imgS;
    }

    public static GreyscaleImage applyBiasLevelToMakePositive(GreyscaleImage gradient) {
        int min = gradient.min();
        if (min >= 0) {
            return gradient;
        }
        GreyscaleImage g2;
        int range = gradient.max() - min;
        if (range > gradient.getMaxAllowed()) {
            g2 = gradient.copyToFullRangeIntImage();
        } else {
            g2 = gradient.copyImage();
        } 
        for (int pixIdx = 0; pixIdx < g2.getNPixels(); ++pixIdx) {
            int v = g2.getValue(pixIdx) - min;
            g2.setValue(pixIdx, v);
        }
        return g2;
    }

    public static boolean containsNegativeNumber(int[] a) {
        for (int ai : a) {
            if (ai < 0) {
                return true;
            }
        }
        return false;
    }
   
    
    public static void _printHistograms_xy(int[][] gHists, int w, int h) {
        PixelHelper ph = new PixelHelper();
        for (int row = 0; row < h; ++row) {
            for (int col = 0; col < w; ++col) {
                int pixIdx = (int)ph.toPixelIndex(col, row, w);
                int[] gh0 = gHists[pixIdx];
                if (hasNonZeroes(gh0)) {
                    System.out.format("(%d,%d) %s\n", col, row, Arrays.toString(gh0));
                }
            }
        }
    }
    public static void _printHistograms_xy_4(int[][] gHists, int w, int h) {
        PixelHelper ph = new PixelHelper();
        for (int row = 0; row < h; ++row) {
            for (int col = 0; col < w; ++col) {
                int pixIdx = (int)ph.toPixelIndex(col, row, w);
                int[] gh0 = gHists[pixIdx];
                //if (hasNonZeroes(gh0)) {
                    StringBuilder sb = new StringBuilder(
                        String.format("(%4d,%4d)", col, row));
                    for (int j = 0; j < gh0.length; ++j) {
                        sb.append(String.format(", %4d", gh0[j]));
                    }
                    System.out.println(sb.toString());
                //}
            }
        }
    }
    private static boolean hasNonZeroes(int[] a) {
        for (int b : a) {
            if (b != 0) {
                return true;
            }
        }
        return false;
    }
}
