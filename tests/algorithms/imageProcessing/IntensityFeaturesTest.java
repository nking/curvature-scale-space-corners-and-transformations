package algorithms.imageProcessing;

import algorithms.compGeometry.RotatedOffsets;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.ResourceFinder;
import java.util.Set;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class IntensityFeaturesTest extends TestCase {
    
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    public void test0() throws Exception {
        
        int cellDim = 2;
        int nCellsAcross = 6;
        int nColsHalf = nCellsAcross / 2;
        int range0 = cellDim * nColsHalf;
        float[] output = new float[nCellsAcross * nCellsAcross];
        float[] xT = new float[cellDim * cellDim];
        float[] yT = new float[xT.length];
        
        int count = 0;
        
        int xCenter = 10;
        int yCenter = 10;
        
        int w = 20;
        int h = 20;
        
        int rotation = 45;
        
        int tc = 0;
                
        for (int dx = -range0; dx < range0; dx += cellDim) {
            for (int dy = -range0; dy < range0; dy += cellDim) {
                
                // --- calculate values for the cell ---
                boolean withinBounds = IntensityFeatures.transformCellCoordinates(
                    rotation, xCenter, yCenter, dx, dy, cellDim, w, h, xT, yT);
                
                StringBuilder sb = new StringBuilder();
                for (int i = 0; i < xT.length; ++i) {
                    sb.append(String.format("(%.1f,%.1f) ", xT[i] - xCenter, 
                        yT[i] - yCenter));
                }
                                
                //log.fine(String.format("prev %d %d (%d, %d): %s", count, tc, 
                //    dx, dy, sb.toString()));
                
                tc += xT.length;
                
                count++;
            }
        }
    }
    
    public void estCalculate45DegreeOrientation_greyscale() throws Exception {
        
        /*
               |3|3|2|1|1|
               |3|3|2|1|1|
               |4|4| |0|0|
               |5|5|6|7|7|
               |5|5|6|7|7|
        */
        
        RotatedOffsets rOffsets = RotatedOffsets.getInstance();
        
        int blockHalfWidths = 5;
        
        boolean useNormalizedIntensities = true;
        
        GreyscaleImage img = new GreyscaleImage(25, 25, 
            GreyscaleImage.Type.Bits32FullRangeInt);
                
        boolean[] sign = new boolean[]{true, false};
        /*
               |3|3|2|1|1|
               |3|3|2|1|1|
               |4|4| |0|0|
               |5|5|6|7|7|
               |5|5|6|7|7|
        */
        int xc = 2;
        int yc = 2;
        for (boolean pos : sign) {
            
            for (int i = 1; i < 2; ++i) {
            
                // values are cached so create a new instance for each "img"
                IntensityFeatures features = new IntensityFeatures(blockHalfWidths,
                    useNormalizedIntensities, rOffsets);
        
                int v = i;
                if (!pos) {
                    v *= -1;
                }
                
                img.setValue(xc + 1, yc, v);
                img.setValue(xc + 2, yc, v);
                
                v = i + 1;
                if (v > 7) {
                    v = v - 8;
                }
                if (!pos) {
                    v *= -1;
                }
                img.setValue(xc + 1, yc + 1, v);
                img.setValue(xc + 2, yc + 1, v);
                img.setValue(xc + 1, yc + 2, v);
                img.setValue(xc + 2, yc + 2, v);
                
                v = i + 2;
                if (v > 7) {
                    v = v - 8;
                }
                if (!pos) {
                    v *= -1;
                }
                img.setValue(xc, yc + 1, v);
                img.setValue(xc, yc + 2, v);
                /*
                   |3|3|2|1|1|
                   |3|3|2|1|1|
                   |4|4| |0|0|
                   |5|5|6|7|7|
                   |5|5|6|7|7|
                */
                
                v = i + 3;
                if (v > 7) {
                    v = v - 8;
                }
                if (!pos) {
                    v *= -1;
                }
                img.setValue(xc - 1, yc + 1, v);
                img.setValue(xc - 2, yc + 1, v);
                img.setValue(xc - 1, yc + 2, v);
                img.setValue(xc - 2, yc + 2, v);
                
                v = i + 4;
                if (v > 7) {
                    v = v - 8;
                }
                if (!pos) {
                    v *= -1;
                }
                img.setValue(xc - 1, yc, v);
                img.setValue(xc - 2, yc, v);
                
                v = i + 5;
                if (v > 7) {
                    v = v - 8;
                }
                if (!pos) {
                    v *= -1;
                }
                img.setValue(xc - 1, yc - 1, v);
                img.setValue(xc - 2, yc - 1, v);
                img.setValue(xc - 1, yc - 2, v);
                img.setValue(xc - 2, yc - 2, v);
                /*
                   |3|3|2|1|1|
                   |3|3|2|1|1|
                   |4|4| |0|0|
                   |5|5|6|7|7|
                   |5|5|6|7|7|
                */
                
                v = i + 6;
                if (v > 7) {
                    v = v - 8;
                }
                if (!pos) {
                    v *= -1;
                }
                img.setValue(xc, yc - 1, v);
                img.setValue(xc, yc - 2, v);
                
                v = i + 7;
                if (v > 7) {
                    v = v - 8;
                }
                if (!pos) {
                    v *= -1;
                }
                img.setValue(xc + 1, yc - 1, v);
                img.setValue(xc + 2, yc - 1, v);
                img.setValue(xc + 1, yc - 2, v);
                img.setValue(xc + 2, yc - 2, v);
                
                int orientation = features.calculate45DegreeOrientation(img, 2, 2);
                                             
                //System.out.println("i=" + i + " orientation=" + orientation + " pos=" + pos);
                
                //printDebug(img); 
                
                switch(i) {
                    case 0: {
                        if (pos) {
                            assertEquals(135, orientation);
                        } else {
                            assertEquals(315, orientation);
                        }
                        break;
                    }
                    case 1: {
                        if (pos) {
                            assertEquals(90, orientation);
                        } else {
                            assertEquals(270, orientation);
                        }
                        break;
                    }
                    case 2: {
                        if (pos) {
                            assertEquals(45, orientation);
                        } else {
                            assertEquals(225, orientation);
                        }
                        break;
                    }
                    case 3: {
                        if (pos) {
                            assertEquals(0, orientation);
                        } else {
                            assertEquals(180, orientation);
                        }
                        break;
                    }
                    case 4: {
                        if (pos) {
                            assertEquals(315, orientation);
                        } else {
                            assertEquals(135, orientation);
                        }
                        break;
                    }
                    case 5: {
                        if (pos) {
                            assertEquals(270, orientation);
                        } else {
                            assertEquals(90, orientation);
                        }
                        break;
                    }
                    case 6: {
                        if (pos) {
                            assertEquals(225, orientation);
                        } else {
                            assertEquals(45, orientation);
                        }
                        break;
                    } 
                    default: {
                        if (pos) {
                            assertEquals(180, orientation);
                        } else {
                            assertEquals(0, orientation);
                        }
                        break;
                    }
                }
            }
        }
        
    }
    
    public void testCalculate45DegreeOrientation_gradient() throws Exception {
        
        RotatedOffsets rOffsets = RotatedOffsets.getInstance();
        
        int blockHalfWidths = 5;
        boolean useNormalizedIntensities = true;
        
        String fileName = "valve_gaussian_subimage.png";
        //String fileName = "checkerboard_subimage.jpg";
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        GreyscaleImage img = ImageIOHelper.readImageAsGrayScaleAvgRGB(filePath);
        
        // values are cached so create a new instance for each "img"
        IntensityFeatures featuresGS = new IntensityFeatures(blockHalfWidths,
            useNormalizedIntensities, rOffsets);
        IntensityFeatures featuresGradient = new IntensityFeatures(blockHalfWidths,
            useNormalizedIntensities, rOffsets);
        featuresGradient.calculateGradientWithGreyscale(img);
        MiscDebug.writeImage(featuresGradient.theta, "_theta");
        
        int orientationGS, x, y, expected, orientationGradient;
        img = ImageIOHelper.readImageAsGrayScaleAvgRGB(filePath);
        
        /*
        x=52  y=53
        x=47  y=11
        x=104 y=46
        x=40  y=43
        
        98, 72
        96, 76
        32, 32
        83, 84
        */
        //makeCannyCorners(filePath);
        //make2ndDerivPoints(filePath);
        
        Logger log = Logger.getLogger(this.getClass().getName());
        
        x = 52; y = 53; expected=0;
        orientationGS = featuresGS.calculate45DegreeOrientation(img, x, y);
        log.fine(String.format("(%d,%d) gs exp=%d angle=%d", x, y, expected, orientationGS));
        orientationGradient = featuresGradient.calculateOrientation(x, y);
        log.fine(String.format("(%d,%d) gr exp=%d angle=%d\n", x, y, 
            expected, orientationGradient));
        assertTrue(Math.abs(orientationGradient - expected) <= 20);
        
        x = 47; y = 11; expected=270;
        orientationGS = featuresGS.calculate45DegreeOrientation(img, x, y);
        log.fine(String.format("(%d,%d) gs exp=%d angle=%d", x, y, expected, orientationGS));
        orientationGradient = featuresGradient.calculateOrientation(x, y);
        log.fine(String.format("(%d,%d) gr exp=%d angle=%d\n", x, y, 
            expected, orientationGradient));
        assertTrue(Math.abs(orientationGradient - expected) <= 20);
        
        x = 40; y = 43; expected=135;
        orientationGS = featuresGS.calculate45DegreeOrientation(img, x, y);
        log.fine(String.format("(%d,%d) gs exp=%d angle=%d", x, y, expected, orientationGS));
        orientationGradient = featuresGradient.calculateOrientation(x, y);
        log.fine(String.format("(%d,%d) gr exp=%d angle=%d\n", x, y, 
            expected, orientationGradient));
        assertTrue(Math.abs(orientationGradient - expected) <= 20);
        
        x = 104; y = 46; expected=200;
        orientationGS = featuresGS.calculate45DegreeOrientation(img, x, y);
        log.fine(String.format("(%d,%d) gs exp=%d angle=%d", x, y, expected, orientationGS));
        orientationGradient = featuresGradient.calculateOrientation(x, y);
        log.fine(String.format("(%d,%d) gr exp=%d angle=%d\n", x, y, 
            expected, orientationGradient));
        assertTrue(Math.abs(orientationGradient - expected) <= 20);
        
        x = 98; y = 72; expected=150;
        orientationGS = featuresGS.calculate45DegreeOrientation(img, x, y);
        log.fine(String.format("(%d,%d) gs exp=%d angle=%d", x, y, expected, orientationGS));
        orientationGradient = featuresGradient.calculateOrientation(x, y);
        log.fine(String.format("(%d,%d) gr exp=%d angle=%d\n", x, y, 
            expected, orientationGradient));
        assertTrue(Math.abs(orientationGradient - expected) <= 20);
        
        x = 83; y = 84; expected=220;
        orientationGS = featuresGS.calculate45DegreeOrientation(img, x, y);
        log.fine(String.format("(%d,%d) gs exp=%d angle=%d", x, y, expected, orientationGS));
        orientationGradient = featuresGradient.calculateOrientation(x, y);
        log.fine(String.format("(%d,%d) gr exp=%d angle=%d\n", x, y, 
            expected, orientationGradient));
        assertTrue(Math.abs(orientationGradient - expected) <= 20);
        
    }

    private void printDebug(GreyscaleImage img) {
        
        StringBuilder sb = new StringBuilder();
        for (int row = 0; row < img.getHeight(); ++row) {
            sb.append(String.format("row %4d: ", row));
            for (int col = 0; col < img.getWidth(); ++col) {
                sb.append(String.format("%4d ", img.getValue(col, row)));
            }
            sb.append("\n");
        }
        
        log.fine(sb.toString());
    }

    private void makeCannyCorners(String filePath) throws Exception {
        boolean useBinned = false;
        boolean filterOutImageBoundaryBlobs = true;
        boolean filterOutZeroPixels = false;
        boolean doNotAddPoints = true;
        
        ImageExt img = ImageIOHelper.readImageExt(filePath);

        ImageExt imgCp = img.copyToImageExt();
        
        /*        
        SegmentationType type = SegmentationType.GREYSCALE_CANNY;
        BlobPerimeterCornerHelper img1Helper = new BlobPerimeterCornerHelper(img);
        img1Helper.applySegmentation(type, useBinned);

        img1Helper.getBlobs(type, useBinned, filterOutImageBoundaryBlobs, filterOutZeroPixels);

        List<List<CornerRegion>> cornersList = img1Helper.getPerimeterCorners(type, useBinned);
        
        for (List<CornerRegion> list : cornersList) {
            for (CornerRegion cr : list) {
                int x = cr.getX()[cr.getKMaxIdx()];
                int y = cr.getY()[cr.getKMaxIdx()];
                System.out.println("x=" + x + " y=" + y);
            }
            ImageIOHelper.addCornerRegionsToImage(list, imgCp, 0, 0, 1);
        }
        
        MiscDebug.writeImage(imgCp, "canny_corners");
        */
        
        CurvatureScaleSpaceCornerDetector detector = new CurvatureScaleSpaceCornerDetector(img);
        detector.findCorners();
        Set<CornerRegion> corners = detector.getEdgeCornerRegionsInOriginalReferenceFrame(true);
        for (CornerRegion cr : corners) {
            int x = cr.getX()[cr.getKMaxIdx()];
            int y = cr.getY()[cr.getKMaxIdx()];
            System.out.println("x=" + x + " y=" + y);
        }
        ImageIOHelper.addCornerRegionsToImage(corners, imgCp, 0, 0, 1);
        
        MiscDebug.writeImage(imgCp, "canny_corners");
    }
    
    private void make2ndDerivPoints(String filePath) throws Exception {
       
        ImageExt img = ImageIOHelper.readImageExt(filePath);

        ImageExt imgCp = img.copyToImageExt();
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        Set<PairInt> pixels = imageProcessor.extract2ndDerivPoints(img.copyToGreyscale(),
            2000, true);
        
        ImageIOHelper.addCurveToImage(pixels, imgCp, 1, 255, 0, 0);
        
        MiscDebug.writeImage(imgCp, "_pts_2nd_deriv");
    }
}
