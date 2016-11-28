package algorithms.imageProcessing;

import algorithms.imageProcessing.segmentation.ColorSpace;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import java.util.List;
import java.util.Set;

/**
 * class to calculate the average and standard deviation of
 * colors within a class of objects.
 * The average and standard deviation of each object is
 * calculated first and then the average and standard
 * deviation of all object average is determined.
 * This class holds this calculation for different color spaces.
 * 
 * @author nichole
 */
public class IntraClassStats {
   
    /**
     * calculate the mean of the images in each color band and the
     * standard deviation of the mean and the SSD of each object's internal
     * standard deviation of the mean.
     * <pre>
     * The result is a two dimensional array of format:
     * [color band][mean of color][stdv of mean][SSD of indiv stdvs]
     * </pre>
     * @param imgs
     * @param shapes
     * @param clrSpace
     * @return two dimensional array of format:
     * [color band][mean of color][stdv of mean][SSD of indiv stdevs]
     */
    public float[][] calculateWithinClassMeans(ImageExt[] imgs,
        List<Set<PairInt>> shapes, ColorSpace clrSpace) {
       
        if (clrSpace.equals(ColorSpace.HSV)) {
            return calcMeanForHSV(imgs, shapes);
        } else if (clrSpace.equals(ColorSpace.CIELAB)) {
            return calcMeanForCIELAB(imgs, shapes);
        }
    }
    
     /**
     * calculate within the class, the differences of the mean colors and the
     * standard deviations of those differences.
     * (NOTE, for some color spaces, you might prefer to use the delteE method).
     * 
     * The result is a two dimensional array of format:
     * [color band][mean diff of color averages][stdv of mean diff][SSD of indiv stdevs]
     * </pre>
     * @param imgs
     * @param shapes
     * @param clrSpace
     * @return two dimensional array of format:
     * [color band][mean diff of color averages][stdv of mean diff][SSD of indiv stdevs]
     */
    public float[][] calculateWithinClassDiffererences(ImageExt[] imgs,
        List<Set<PairInt>> shapes, ColorSpace clrSpace) {
        
        if (clrSpace.equals(ColorSpace.HSV)) {
            return calcDiffForHSV(imgs, shapes);
        } else if (clrSpace.equals(ColorSpace.CIELAB)) {
            return calcDiffForCIELAB(imgs, shapes);
        }
    }
    
    /**
     * calculate within the class, the differences of the mean colors using
     * the DELTA E 2000 formula, and calculate the
     * standard deviations of those differences.
     * (NOTE, that deltaE is only relevant for some color spaces, and for those
     * you can use the difference method instead).
     * 
     * The result is a two dimensional array of format:
     * [color band][mean deltaE of color averages][stdv of mean diff][SSD of indiv stdevs]
     * </pre>
     * @param imgs
     * @param shapes
     * @param clrSpace
     * @return two dimensional array of format:
     * [color band][mean deltaE of color averages][stdv of mean diff][SSD of indiv stdevs]
     */
    public float[][] calculateWithinClassDeltaE(ImageExt[] imgs,
        List<Set<PairInt>> shapes, ColorSpace clrSpace) {
        
        if (clrSpace.equals(ColorSpace.HSV)) {
            throw new IllegalArgumentException(" for HSV, use "
                + "calculateWithinClassDiffererences instead");
        } else if (clrSpace.equals(ColorSpace.CIELAB)) {
            return calcDeltaEForCIELAB(imgs, shapes);
        }
    }
    
    // ----- mean calculations ------
    
    //  [clr Idx][avg, stdv1, stdv2]
    protected float[][] calcMeanForHSV(ImageExt[] imgs, List<Set<PairInt>> shapes) {
       
        int nImgs = imgs.length;
        
        float[] a1 = new float[nImgs];  
        float[] a2 = new float[nImgs];  
        float[] a3 = new float[nImgs];
        float[] s1 = new float[nImgs];  
        float[] s2 = new float[nImgs];  
        float[] s3 = new float[nImgs];

        int count = 0;
        for (int nImg = 0; nImg < nImgs; ++nImg) {
            ImageExt img = imgs[nImg];
            Set<PairInt> shape = shapes.get(nImg);
            GroupPixelHSV gr = new GroupPixelHSV();
            gr.calculateColors(shape, img);
            a1[count] = gr.getAvgH();
            a2[count] = gr.getAvgS();
            a3[count] = gr.getAvgV();
            s1[count] = gr.getStdDevH();
            s2[count] = gr.getStdDevS();
            s3[count] = gr.getStdDevV();
            count++;
        }
        
        float[][] results = new float[3][3];
        results[0] = MiscMath.calcMeanAndStDev(a1, s1);
        results[1] = MiscMath.calcMeanAndStDev(a2, s2);
        results[2] = MiscMath.calcMeanAndStDev(a3, s3);
        
        return results;
    }
   
    //  [clr Idx][avg, stdv1, stdv2]
    protected float[][] calcMeanForCIELAB(ImageExt[] imgs, List<Set<PairInt>> shapes) {
       
        int nImgs = imgs.length;
        
        float[] a1 = new float[nImgs];  
        float[] a2 = new float[nImgs];  
        float[] a3 = new float[nImgs];
        float[] s1 = new float[nImgs];  
        float[] s2 = new float[nImgs];  
        float[] s3 = new float[nImgs];

        int count = 0;
        for (int nImg = 0; nImg < nImgs; ++nImg) {
            ImageExt img = imgs[nImg];
            Set<PairInt> shape = shapes.get(nImg);
            GroupPixelCIELAB gr = new GroupPixelCIELAB(shape, img);
            gr.calculateColors(shape, img, 0, 0);
            a1[count] = gr.getAvgL();
            a2[count] = gr.getAvgA();
            a3[count] = gr.getAvgB();
            s1[count] = gr.getStdDevL();
            s2[count] = gr.getStdDevA();
            s3[count] = gr.getStdDevB();
            count++;
        }
        
        float[][] results = new float[3][3];
        results[0] = MiscMath.calcMeanAndStDev(a1, s1);
        results[1] = MiscMath.calcMeanAndStDev(a2, s2);
        results[2] = MiscMath.calcMeanAndStDev(a3, s3);
        
        return results;
    }
   
    // ------ differences ------  
    
    //  [clr Idx][avg, stdv1, stdv2]
    protected float[][] calcDiffForHSV(ImageExt[] imgs, List<Set<PairInt>> shapes) {
       
        int nImgs = imgs.length;
        
        float[] a1 = new float[nImgs];  
        float[] a2 = new float[nImgs];  
        float[] a3 = new float[nImgs];
        float[] s1 = new float[nImgs];  
        float[] s2 = new float[nImgs];  
        float[] s3 = new float[nImgs];

        int count = 0;
        for (int nImg = 0; nImg < nImgs; ++nImg) {
            ImageExt img = imgs[nImg];
            Set<PairInt> shape = shapes.get(nImg);
            GroupPixelHSV gr = new GroupPixelHSV();
            gr.calculateColors(shape, img);
            a1[count] = gr.getAvgH();
            a2[count] = gr.getAvgS();
            a3[count] = gr.getAvgV();
            s1[count] = gr.getStdDevH();
            s2[count] = gr.getStdDevS();
            s3[count] = gr.getStdDevV();
            count++;
        }
        
        int nTot = nImgs * (nImgs - 1);
        float[] d1 = new float[nTot];
        float[] d2 = new float[nTot];
        float[] d3 = new float[nTot];
        float[] ds1 = new float[nTot];
        float[] ds2 = new float[nTot];
        float[] ds3 = new float[nTot];
        count = 0;
        for (int i = 0; i < nImgs; ++i) {
            for (int j = (i + 1); j < nImgs; ++j) {
                d1[count] = Math.abs(a1[i] - a1[j]);
                ds1[count] = (float)Math.sqrt(s1[i]*s1[i] + s1[j]*s1[j]);
                
                d2[count] = Math.abs(a2[i] - a2[j]);
                ds2[count] = (float)Math.sqrt(s2[i]*s2[i] + s2[j]*s2[j]);
                
                d3[count] = Math.abs(a3[i] - a3[j]);
                ds3[count] = (float)Math.sqrt(s3[i]*s3[i] + s3[j]*s3[j]);
                
                count++;
            }
        }
        assert(count == nTot);
        
        float[][] results = new float[3][3];
        results[0] = MiscMath.calcMeanAndStDev(d1, ds1);
        results[1] = MiscMath.calcMeanAndStDev(d2, ds2);
        results[2] = MiscMath.calcMeanAndStDev(d3, ds3);
        
        return results;
    }
   
    //  [clr Idx][avg, stdv1, stdv2]
    protected float[][] calcDiffForCIELAB(ImageExt[] imgs, List<Set<PairInt>> shapes) {
       
        int nImgs = imgs.length;
        
        float[] a1 = new float[nImgs];  
        float[] a2 = new float[nImgs];  
        float[] a3 = new float[nImgs];
        float[] s1 = new float[nImgs];  
        float[] s2 = new float[nImgs];  
        float[] s3 = new float[nImgs];

        int count = 0;
        for (int nImg = 0; nImg < nImgs; ++nImg) {
            ImageExt img = imgs[nImg];
            Set<PairInt> shape = shapes.get(nImg);
            GroupPixelCIELAB gr = new GroupPixelCIELAB(shape, img);
            gr.calculateColors(shape, img, 0, 0);
            a1[count] = gr.getAvgL();
            a2[count] = gr.getAvgA();
            a3[count] = gr.getAvgB();
            s1[count] = gr.getStdDevL();
            s2[count] = gr.getStdDevA();
            s3[count] = gr.getStdDevB();
            count++;
        }
        
        int nTot = nImgs * (nImgs - 1);
        float[] d1 = new float[nTot];
        float[] d2 = new float[nTot];
        float[] d3 = new float[nTot];
        float[] ds1 = new float[nTot];
        float[] ds2 = new float[nTot];
        float[] ds3 = new float[nTot];
        count = 0;
        for (int i = 0; i < nImgs; ++i) {
            for (int j = (i + 1); j < nImgs; ++j) {
                d1[count] = Math.abs(a1[i] - a1[j]);
                ds1[count] = (float)Math.sqrt(s1[i]*s1[i] + s1[j]*s1[j]);
                
                d2[count] = Math.abs(a2[i] - a2[j]);
                ds2[count] = (float)Math.sqrt(s2[i]*s2[i] + s2[j]*s2[j]);
                
                d3[count] = Math.abs(a3[i] - a3[j]);
                ds3[count] = (float)Math.sqrt(s3[i]*s3[i] + s3[j]*s3[j]);
                
                count++;
            }
        }
        assert(count == nTot);
        
        float[][] results = new float[3][3];
        results[0] = MiscMath.calcMeanAndStDev(d1, ds1);
        results[1] = MiscMath.calcMeanAndStDev(d2, ds2);
        results[2] = MiscMath.calcMeanAndStDev(d3, ds3);
        
        return results;
    }
    
    // ------ deltaEs ------  
    //  [deltaE][avg, stdv1]
    protected float[][] calcDeltaEForCIELAB(ImageExt[] imgs, List<Set<PairInt>> shapes) {
       
        int nImgs = imgs.length;
        
        float[] a1 = new float[nImgs];  
        float[] a2 = new float[nImgs];  
        float[] a3 = new float[nImgs];
        float[] s1 = new float[nImgs];  
        float[] s2 = new float[nImgs];  
        float[] s3 = new float[nImgs];

        int count = 0;
        for (int nImg = 0; nImg < nImgs; ++nImg) {
            ImageExt img = imgs[nImg];
            Set<PairInt> shape = shapes.get(nImg);
            GroupPixelCIELAB gr = new GroupPixelCIELAB(shape, img);
            gr.calculateColors(shape, img, 0, 0);
            a1[count] = gr.getAvgL();
            a2[count] = gr.getAvgA();
            a3[count] = gr.getAvgB();
            s1[count] = gr.getStdDevL();
            s2[count] = gr.getStdDevA();
            s3[count] = gr.getStdDevB();
            count++;
        }
        
        CIEChromaticity cieC = new CIEChromaticity();
        
        // calculating the errors of the means is a complex use of variables
        //   because of the transformation into deltaE space.
        //   a short cut is to assume a confidence level and then apply that
        //   factor to stdev to determine the lower and upper boundaries of
        //   L,A,B then the max of the deltaEs can be used to derive a rough stdev
        //   ... can compare that to the stdev among the means below.
        //   for now, just supplying the later below.
        
        int nTot = nImgs * (nImgs - 1);
        float[] d1 = new float[nTot];
        count = 0;
        for (int i = 0; i < nImgs; ++i) {
            for (int j = (i + 1); j < nImgs; ++j) {
                d1[count] = (float)cieC.calcDeltaECIE2000(
                    a1[i], a2[i], a3[i], a1[j], a2[j], a3[j]);
                count++;
            }
        }
        assert(count == nTot);
        
        float[][] results = new float[1][2];
        results[0] = MiscMath.getAvgAndStDev(d1);
        
        return results;
    }
    
}