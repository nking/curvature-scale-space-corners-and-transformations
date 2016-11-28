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
 * deviation of all object averages is determined.
 * This class holds this calculation for different color spaces.
 * It also holds similar calculations but for the differences between the means.
 * @author nichole
 */
public class IntraClassColorStats {
   
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
        } else if (clrSpace.equals(ColorSpace.CIELAB1931)) {
            return calcMeanForCIELAB1931(imgs, shapes);
        } else if (clrSpace.equals(ColorSpace.CIELUV)) {
            return calcMeanForCIELUV(imgs, shapes);
        } else if (clrSpace.equals(ColorSpace.CIELCH)) {
            return calcMeanForCIELCH(imgs, shapes);
        }
        
        throw new IllegalArgumentException("color space is not implementend");
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
        } else if (clrSpace.equals(ColorSpace.CIELAB1931)) {
            return calcDiffForCIELAB1931(imgs, shapes);
        } else if (clrSpace.equals(ColorSpace.CIELUV)) {
            return calcDiffForCIELUV(imgs, shapes);
        } else if (clrSpace.equals(ColorSpace.CIELCH)) {
            return calcDiffForCIELCH(imgs, shapes);
        }
        throw new IllegalArgumentException("color space is not implementend");
    }
    
    /**
     * calculate within the class, the differences of the mean colors and the
     * standard deviations of those differences.
     * (NOTE, for some color spaces, you might prefer to use the delteE method).
     * 
     * The result is a two dimensional array of format:
     * [color band][mean diff of color averages][stdv of mean diff][SSD of indiv stdevs]
     * 
     * The normalization uses results from using a standard illuminant for daylight, D65.
     * </pre>
     * @param imgs
     * @param shapes
     * @param clrSpace
     * @return two dimensional array of format:
     * [color band][mean diff of color averages][stdv of mean diff][SSD of indiv stdevs]
     */
    public float[] calculateWithinClassNormalizedDiffererences(ImageExt[] imgs,
        List<Set<PairInt>> shapes, ColorSpace clrSpace) {
        
        if (clrSpace.equals(ColorSpace.HSV)) {
            return calcNormDiffForHSV(imgs, shapes);
        } else if (clrSpace.equals(ColorSpace.CIELAB)) {
            return calcNormDiffForCIELAB(imgs, shapes);
        } else if (clrSpace.equals(ColorSpace.CIELAB1931)) {
            return calcNormDiffForCIELAB1931(imgs, shapes);
        } else if (clrSpace.equals(ColorSpace.CIELUV)) {
            return calcNormDiffForCIELUV(imgs, shapes);
        } else if (clrSpace.equals(ColorSpace.CIELCH)) {
            return calcNormDiffForCIELCH(imgs, shapes);
        }
        throw new IllegalArgumentException("color space is not implementend");
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
     * [mean deltaE of color averages, stdv of mean diff]
     */
    public float[] calculateWithinClassDeltaE(ImageExt[] imgs,
        List<Set<PairInt>> shapes, ColorSpace clrSpace) {
        
        if (clrSpace.equals(ColorSpace.HSV)) {
            throw new IllegalArgumentException(" for HSV, use "
                + "calculateWithinClassDiffererences instead");
        } else if (clrSpace.equals(ColorSpace.CIELAB)) {
            return calcDeltaEForCIELAB(imgs, shapes);
        } else if (clrSpace.equals(ColorSpace.CIELAB1931)) {
            return calcDeltaEForCIELAB31(imgs, shapes);
        } else if (clrSpace.equals(ColorSpace.CIELUV)) {
            return calcDeltaEForCIELUV(imgs, shapes);
        }
        
        throw new IllegalArgumentException("color space is not implementend");
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
    
    //  [clr Idx][avg, stdv1, stdv2]
    protected float[][] calcMeanForCIELAB1931(ImageExt[] imgs, List<Set<PairInt>> shapes) {
       
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
            GroupPixelCIELAB1931 gr = new GroupPixelCIELAB1931(shape, img);
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
    
    //  [clr Idx][avg, stdv1, stdv2]
    protected float[][] calcMeanForCIELUV(ImageExt[] imgs, List<Set<PairInt>> shapes) {
       
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
            GroupPixelCIELUV gr = new GroupPixelCIELUV(shape, img);
            gr.calculateColors(shape, img, 0, 0);
            a1[count] = gr.getAvgL();
            a2[count] = gr.getAvgU();
            a3[count] = gr.getAvgV();
            s1[count] = gr.getStdDevL();
            s2[count] = gr.getStdDevU();
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
    protected float[][] calcMeanForCIELCH(ImageExt[] imgs, List<Set<PairInt>> shapes) {

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
            GroupPixelCIELCH gr = new GroupPixelCIELCH(shape, img);
            gr.calculateColors(shape, img, 0, 0);
            a1[count] = gr.getAvgL();
            a2[count] = gr.getAvgC();
            a3[count] = gr.getAvgH();
            s1[count] = gr.getStdDevL();
            s2[count] = gr.getStdDevC();
            s3[count] = gr.getStdDevH();
            count++;
        }
        
        float[][] results = new float[3][3];
        results[0] = MiscMath.calcMeanAndStDev(a1, s1);
        results[1] = MiscMath.calcMeanAndStDev(a2, s2);
        results[2] = MiscMath.calcMeanAndStDevWithWrapAround(a3, 359, s3);
        
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
        
        int nTot = nImgs * (nImgs - 1)/2;
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
        
        int nTot = nImgs * (nImgs - 1)/2;
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
    protected float[][] calcDiffForCIELAB1931(ImageExt[] imgs, List<Set<PairInt>> shapes) {
       
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
            GroupPixelCIELAB1931 gr = new GroupPixelCIELAB1931(shape, img);
            gr.calculateColors(shape, img, 0, 0);
            a1[count] = gr.getAvgL();
            a2[count] = gr.getAvgA();
            a3[count] = gr.getAvgB();
            s1[count] = gr.getStdDevL();
            s2[count] = gr.getStdDevA();
            s3[count] = gr.getStdDevB();
            count++;
        }
        
        int nTot = nImgs * (nImgs - 1)/2;
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
    protected float[][] calcDiffForCIELUV(ImageExt[] imgs, List<Set<PairInt>> shapes) {
       
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
            GroupPixelCIELUV gr = new GroupPixelCIELUV(shape, img);
            gr.calculateColors(shape, img, 0, 0);
            a1[count] = gr.getAvgL();
            a2[count] = gr.getAvgU();
            a3[count] = gr.getAvgV();
            s1[count] = gr.getStdDevL();
            s2[count] = gr.getStdDevU();
            s3[count] = gr.getStdDevV();
            count++;
        }
        
        int nTot = nImgs * (nImgs - 1)/2;
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
    protected float[][] calcDiffForCIELCH(ImageExt[] imgs, List<Set<PairInt>> shapes) {
       
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
            GroupPixelCIELCH gr = new GroupPixelCIELCH(shape, img);
            gr.calculateColors(shape, img, 0, 0);
            a1[count] = gr.getAvgL();
            a2[count] = gr.getAvgC();
            a3[count] = gr.getAvgH();
            s1[count] = gr.getStdDevL();
            s2[count] = gr.getStdDevC();
            s3[count] = gr.getStdDevH();
            count++;
        }
        
        int nTot = nImgs * (nImgs - 1)/2;
        float[] d1 = new float[nTot];
        float[] d2 = new float[nTot];
        float[] d3 = new float[nTot];
        float[] ds1 = new float[nTot];
        float[] ds2 = new float[nTot];
        float[] ds3 = new float[nTot];
        count = 0;
        
        float wrapAround = 359;
        float v1, v2;
        for (int i = 0; i < nImgs; ++i) {
            for (int j = (i + 1); j < nImgs; ++j) {
                
                d1[count] = Math.abs(a1[i] - a1[j]);
                ds1[count] = (float)Math.sqrt(s1[i]*s1[i] + s1[j]*s1[j]);
                
                d2[count] = Math.abs(a2[i] - a2[j]);
                ds2[count] = (float)Math.sqrt(s2[i]*s2[i] + s2[j]*s2[j]);
                
                ds3[count] = (float)Math.sqrt(s3[i]*s3[i] + s3[j]*s3[j]);
                v1 = a3[i];
                v2 = a3[j];
                if (v1 > v2) {
                    if (Math.abs(v1 - v2) < Math.abs(v1 - (v2 + wrapAround))) {
                        d3[count] = Math.abs(v1 - v2);
                    } else {
                        d3[count] = Math.abs(v1 - (v2 + wrapAround));
                    }
                } else {
                    if (Math.abs(v1 - v2) < Math.abs((v1 + wrapAround) - v2)) {
                        d3[count] = Math.abs(v1 - v2);
                    } else {
                        d3[count] = Math.abs((v1 + wrapAround) - v2);
                    }
                }
                
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
    protected float[] calcDeltaEForCIELAB(ImageExt[] imgs, List<Set<PairInt>> shapes) {
       
        int nImgs = imgs.length;
        
        float[] a1 = new float[nImgs];  
        float[] a2 = new float[nImgs];  
        float[] a3 = new float[nImgs];

        int count = 0;
        for (int nImg = 0; nImg < nImgs; ++nImg) {
            ImageExt img = imgs[nImg];
            Set<PairInt> shape = shapes.get(nImg);
            GroupPixelCIELAB gr = new GroupPixelCIELAB(shape, img);
            gr.calculateColors(shape, img, 0, 0);
            a1[count] = gr.getAvgL();
            a2[count] = gr.getAvgA();
            a3[count] = gr.getAvgB();
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
        
        int nTot = nImgs * (nImgs - 1)/2;
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
        
        float[] results = MiscMath.getAvgAndStDev(d1);
        
        return results;
    }
    
    //  [deltaE][avg, stdv1]
    protected float[] calcDeltaEForCIELAB31(ImageExt[] imgs, List<Set<PairInt>> shapes) {
       
        int nImgs = imgs.length;
        
        float[] a1 = new float[nImgs];  
        float[] a2 = new float[nImgs];  
        float[] a3 = new float[nImgs];

        int count = 0;
        for (int nImg = 0; nImg < nImgs; ++nImg) {
            ImageExt img = imgs[nImg];
            Set<PairInt> shape = shapes.get(nImg);
            GroupPixelCIELAB1931 gr = new GroupPixelCIELAB1931(shape, img);
            gr.calculateColors(shape, img, 0, 0);
            a1[count] = gr.getAvgL();
            a2[count] = gr.getAvgA();
            a3[count] = gr.getAvgB();
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
        
        int nTot = nImgs * (nImgs - 1)/2;
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
        
        float[] results = MiscMath.getAvgAndStDev(d1);
        
        return results;
    }
    
    //  [deltaE][avg, stdv1]
    protected float[] calcDeltaEForCIELUV(ImageExt[] imgs, List<Set<PairInt>> shapes) {
       
        int nImgs = imgs.length;
        
        float[] a1 = new float[nImgs];  
        float[] a2 = new float[nImgs];  
        float[] a3 = new float[nImgs];

        int count = 0;
        for (int nImg = 0; nImg < nImgs; ++nImg) {
            ImageExt img = imgs[nImg];
            Set<PairInt> shape = shapes.get(nImg);
            GroupPixelCIELUV gr = new GroupPixelCIELUV(shape, img);
            gr.calculateColors(shape, img, 0, 0);
            a1[count] = gr.getAvgL();
            a2[count] = gr.getAvgU();
            a3[count] = gr.getAvgV();
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
        
        int nTot = nImgs * (nImgs - 1)/2;
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
        
        float[] results = MiscMath.getAvgAndStDev(d1);
        
        return results;
    }
    
    // ------ normalized differences ----
    private float[] calcNormDiffForHSV(ImageExt[] imgs, 
        List<Set<PairInt>> shapes) {
        
        float[][] r = calcDiffForHSV(imgs, shapes);

        // averages are already normalized to "1"        

        float sum1 = (r[0][0] + r[1][0] + r[2][0])/3.f;
        
        float sum2 = (float)Math.sqrt(r[0][1]*r[0][1] + r[1][1]*r[1][1] 
            + r[2][1] * r[2][1])/3.f;
        
        float sum3 = (float)Math.sqrt(r[0][2]*r[0][2] + r[1][2]*r[1][2] 
            + r[2][2] * r[2][2])/3.f;
        
        return new float[]{sum1, sum2, sum3};        
    }
    
    private float[] calcNormDiffForCIELAB(ImageExt[] imgs, 
        List<Set<PairInt>> shapes) {
        
        float[][] r = calcDiffForCIELAB(imgs, shapes);
        
        /*
        * using the standard illuminant of daylight, D65,
        * the range of return values is
        *    L    0 to 28.5
        *    A  -46.9  62.5
        *    B  -45.7  48.0
        */
        
        float d1 = 1.f/28.5f;
        float d2 = 1.f/(62.5f + 46.9f);
        float d3 = 1.f/(48.0f + 45.7f);
        
        r[0][0] *= d1; 
        r[1][0] *= d2;
        r[2][0] *= d3;
        
        float sum1 = (r[0][0] + r[1][0] + r[2][0])/3.f;
        
        float sum2 = (float)Math.sqrt(r[0][1]*r[0][1]*d1*d1 
            + r[1][1]*r[1][1]*d2*d2 
            + r[2][1] * r[2][1]*d2*d2)/3.f;
        
        float sum3 = (float)Math.sqrt(r[0][2]*r[0][2]*d1*d1 
            + r[1][2]*r[1][2]*d2*d2 
            + r[2][2] * r[2][2]*d3*d3)/3.f;
        
        return new float[]{sum1, sum2, sum3};
    }

    private float[] calcNormDiffForCIELAB1931(ImageExt[] imgs, List<Set<PairInt>> shapes) {

        float[][] r = calcDiffForCIELAB1931(imgs, shapes);
        
        /*
        * using the standard illuminant of daylight, D65,
        * the range of return values is
        * L*    0 to 104.5
        * a* -190 to 103
        * b* -113 to 99
        */
        
        float d1 = 1.f/104.5f;
        float d2 = 1.f/(103f + 190f);
        float d3 = 1.f/(99f + 113.f);
        
        r[0][0] *= d1; 
        r[1][0] *= d2;
        r[2][0] *= d3;
        
        float sum1 = (r[0][0] + r[1][0] + r[2][0])/3.f;
        
        float sum2 = (float)Math.sqrt(r[0][1]*r[0][1]*d1*d1 
            + r[1][1]*r[1][1]*d2*d2 
            + r[2][1] * r[2][1]*d2*d2)/3.f;
        
        float sum3 = (float)Math.sqrt(r[0][2]*r[0][2]*d1*d1 
            + r[1][2]*r[1][2]*d2*d2 
            + r[2][2] * r[2][2]*d3*d3)/3.f;
        
        return new float[]{sum1, sum2, sum3};
    }

    private float[] calcNormDiffForCIELUV(ImageExt[] imgs, List<Set<PairInt>> shapes) {
        
        /*
           L       0 to 104.5
           u   -86.9 to 183.8
           v  -141.4 to 112.3
        */
        
        float[][] r = calcDiffForCIELUV(imgs, shapes);
        
        float d1 = 1.f/104.5f;
        float d2 = 1.f/(183.8f + 86.9f);
        float d3 = 1.f/(112.3f + 141.4f);
        
        r[0][0] *= d1; 
        r[1][0] *= d2;
        r[2][0] *= d3;
        
        float sum1 = (r[0][0] + r[1][0] + r[2][0])/3.f;
        
        float sum2 = (float)Math.sqrt(r[0][1]*r[0][1]*d1*d1 
            + r[1][1]*r[1][1]*d2*d2 
            + r[2][1] * r[2][1]*d2*d2)/3.f;
        
        float sum3 = (float)Math.sqrt(r[0][2]*r[0][2]*d1*d1 
            + r[1][2]*r[1][2]*d2*d2 
            + r[2][2] * r[2][2]*d3*d3)/3.f;
        
        return new float[]{sum1, sum2, sum3};        
    }

    private float[] calcNormDiffForCIELCH(ImageExt[] imgs, List<Set<PairInt>> shapes) {
        /*
        the range of return values is
        *   luminosity L*  0 to 104.5
        *   magnitude, C:  0 to 139 
        *   angle,     H:  0 to 359
        */
        
        float[][] r = calcDiffForCIELCH(imgs, shapes);
        
        float d1 = 1.f/104.5f;
        float d2 = 1.f/139.f;
        float d3 = 1.f/359.f;
        
        r[0][0] *= d1; 
        r[1][0] *= d2;
        r[2][0] *= d3;
        
        float sum1 = (r[0][0] + r[1][0] + r[2][0])/3.f;
        
        float sum2 = (float)Math.sqrt(r[0][1]*r[0][1]*d1*d1 
            + r[1][1]*r[1][1]*d2*d2 
            + r[2][1] * r[2][1]*d2*d2)/3.f;
        
        float sum3 = (float)Math.sqrt(r[0][2]*r[0][2]*d1*d1 
            + r[1][2]*r[1][2]*d2*d2 
            + r[2][2] * r[2][2]*d3*d3)/3.f;
        
        return new float[]{sum1, sum2, sum3};
    }
    
}