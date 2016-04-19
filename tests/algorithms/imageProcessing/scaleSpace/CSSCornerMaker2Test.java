package algorithms.imageProcessing.scaleSpace;

import algorithms.imageProcessing.EdgeExtractorSimple;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.SIGMA;
import algorithms.imageProcessing.features.FeatureHelper;
import algorithms.imageProcessing.features.PhaseCongruencyDetector;
import algorithms.misc.MiscDebug;
import algorithms.util.CornerArray;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class CSSCornerMaker2Test extends TestCase {
    
    public CSSCornerMaker2Test() {
    }

    public void test0() throws Exception {
        
        String fileName = "small_shapes_for_line_thinners.png";
        //String fileName = "merton_college_I_001.jpg";
        //String fileName = "house.gif";

        String filePath = ResourceFinder.findFileInTestResources(fileName);
        
        GreyscaleImage img = ImageIOHelper.readImage(filePath).copyToGreyscale(); 
        int nCols = img.getWidth();
        int nRows = img.getHeight();
        
        float cutOff = 0.5f;//0.3f;//0.5f;
        int nScale = 5;//5;
        int minWavelength = nScale;//3;
        float mult = 2.1f;
        float sigmaOnf = 0.55f;
        int k = 5;//10;//2;
        float g = 10; 
        float deviationGain = 1.5f;
        int noiseMethod = -1;
        double tLow = 0.1;
        double tHigh = 0.3;
        boolean increaseKIfNeeded = false;
        PhaseCongruencyDetector phaseCDetector = new PhaseCongruencyDetector();
        PhaseCongruencyDetector.PhaseCongruencyProducts products =
            phaseCDetector.phaseCongMono(img, nScale, minWavelength, mult, 
                sigmaOnf, k, increaseKIfNeeded,
                cutOff, g, deviationGain, noiseMethod, tLow, tHigh);
        assertNotNull(products);
        
        Image tmp = new Image(nCols, nRows);
        GreyscaleImage pcImg = new GreyscaleImage(nCols, nRows);
        GreyscaleImage thetaImg = new GreyscaleImage(nCols, nRows);
        for (int i = 0; i < nRows; ++i) {
            for (int j = 0; j < nCols; ++j) {
                if (products.getThinned()[i][j] > 0) {
                    tmp.setRGB(j, i, 255, 255, 255);
                }
                pcImg.setValue(j, i, 
                    (int)Math.round(255. * products.getPhaseCongruency()[i][j]));
                thetaImg.setValue(j, i, 
                    (int)Math.round(products.getOrientation()[i][j]));
            }
        }
        
        EdgeExtractorSimple edgeExtractor = new EdgeExtractorSimple(
            products.getThinned());
        edgeExtractor.extractEdges();
        
        Set<PairInt> junctions = edgeExtractor.getJunctions();
        
        List<PairIntArray> theEdges = edgeExtractor.getEdges();
        // put edges in ref frame of image
        for (int i = 0; i < theEdges.size(); ++i) {
            PairIntArray edge = theEdges.get(i);
            PairIntArray edge2 = new PairIntArray(edge.getN());
            for (int j = 0; j < edge.getN(); ++j) {
                int x = edge.getX(j);
                int y = edge.getY(j);
                edge2.add(y, x);
            }
            theEdges.set(i, edge2);
        }
        //FeatureHelper.filterByLocalizability
        CSSCornerMaker cornerMaker = new CSSCornerMaker(nCols, nRows);
        List<CornerArray> cornerList = cornerMaker.findCornersInScaleSpaceMaps(theEdges);
       
        for (CornerArray corners : cornerList) {

            float blurPix = 2.35f * SIGMA.getValue(corners.getSIGMA());            
            System.out.println("sigma=" + corners.getSIGMA() + " blurPix=" + blurPix);
            
            for (int i = 0; i < corners.getN(); ++i) {
                ImageIOHelper.addPointToImage(corners.getX(i), corners.getY(i), 
                    tmp, 0, 255, 0, 0);
                
                int x = corners.getX(i);
                int y = corners.getY(i);
                
                corners.set(i, y, x, 
                    corners.getCurvature(i),
                    corners.getXFirstDeriv(i),
                    corners.getXSecondDeriv(i),
                    corners.getYFirstDeriv(i),
                    corners.getYSecondDeriv(i), corners.getInt(i));
            }
        }
        
        for (PairInt p : junctions) {
            ImageIOHelper.addPointToImage(p.getY(), p.getX(), 
                tmp, 0, 0, 255, 0);
        }
        
        MiscDebug.writeImage(tmp, "_CORNERS_");
        
        
        tmp = new Image(nCols, nRows);
        ImageIOHelper.addAlternatingColorCurvesToImage(theEdges, tmp, 0);
        
        MiscDebug.writeImage(tmp, "_EDGES_");
        
        cornerList = FeatureHelper.filterByLocalizability(img, pcImg, thetaImg, 
            cornerList);
        
        int s0 = 0;
        int s1 = 0;
        int s2 = 0;
        // assert min number of corners
        for (int ii = 0; ii < cornerList.size(); ++ii) {
            
            CornerArray corners = cornerList.get(ii);
            PairIntArray edge = theEdges.get(ii);
            
            for (int i = 0; i < edge.getN(); ++i) {
                int x = edge.getX(i);
                int y = edge.getY(i);
                if (x==22 && y==88) {
                    s0 = corners.getN();
                    break;
                } else if (x==20 && y==20) {
                    s1 = corners.getN();
                    break;
                } else if (x==61 && y==50) {
                    s2 = corners.getN();
                    break;
                }
            }
        }
        
        assertTrue(s0 >= 4);
        assertTrue(s1 >= 3);
        assertTrue(s2 >= 1);
    }
    
    public void est1() throws Exception {
        
        String fileName = "small_shapes_for_line_thinners.png";

        String filePath = ResourceFinder.findFileInTestResources(fileName);
        
        GreyscaleImage img = ImageIOHelper.readImage(filePath).copyToGreyscale(); 
        int nCols = img.getWidth();
        int nRows = img.getHeight();
        
        float cutOff = 0.5f;//0.3f;//0.5f;
        int nScale = 5;//5;
        int minWavelength = nScale;//3;
        float mult = 2.1f;
        float sigmaOnf = 0.55f;
        int k = 5;//10;//2;
        float g = 10; 
        float deviationGain = 1.5f;
        int noiseMethod = -1;
        double tLow = 0.1;
        double tHigh = 0.3;
        boolean increaseKIfNeeded = false;
        PhaseCongruencyDetector phaseCDetector = new PhaseCongruencyDetector();
        PhaseCongruencyDetector.PhaseCongruencyProducts products =
            phaseCDetector.phaseCongMono(img, nScale, minWavelength, mult, 
                sigmaOnf, k, increaseKIfNeeded,
                cutOff, g, deviationGain, noiseMethod, tLow, tHigh);
        assertNotNull(products);
        
        EdgeExtractorSimple edgeExtractor = new EdgeExtractorSimple(
            products.getThinned());
        edgeExtractor.extractEdges();
                
        List<PairIntArray> theEdges = edgeExtractor.getEdges();
        // put edges in ref frame of image
        for (int i = 0; i < theEdges.size(); ++i) {
            PairIntArray edge = theEdges.get(i);
            PairIntArray edge2 = new PairIntArray(edge.getN());
            for (int j = 0; j < edge.getN(); ++j) {
                int x = edge.getX(j);
                int y = edge.getY(j);
                edge2.add(y, x);
            }
            theEdges.set(i, edge2);
        }
        
        CSSCornerMaker cornerMaker = new CSSCornerMaker(nCols, nRows);
        
        List<Map<SIGMA, ScaleSpaceCurve>> outEdgeScaleSpaceMaps =
            new ArrayList<Map<SIGMA, ScaleSpaceCurve>>();
        
        //NOTE: this method is not efficient because it is expected the use
        // of CornerRegions will be obsolete after the next major refactoring.
        List<CornerArray> cornersList = cornerMaker.findCornersInScaleSpaceMaps(
            theEdges, outEdgeScaleSpaceMaps);
                
        int xTest = 29;
        int yTest = 88;
        
        for (int lIdx = 0; lIdx < outEdgeScaleSpaceMaps.size(); ++lIdx) {
             
            Map<SIGMA, ScaleSpaceCurve> mapOfScaleSpacesForAnEdge = 
                outEdgeScaleSpaceMaps.get(lIdx);
                
            PairIntArray edge = theEdges.get(lIdx);
            
            boolean found = false;
            for (int i = 0; i < edge.getN(); ++i) {
                if ((Math.abs(xTest - edge.getX(i)) < 5) && Math.abs(yTest - edge.getY(i)) < 5) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                continue;
            }
            
            CornerArray ca = cornersList.get(lIdx);
                          
            /*
            need to print i (x,y) sigma k
            */
            for (Map.Entry<SIGMA, ScaleSpaceCurve> entry2 : mapOfScaleSpacesForAnEdge.entrySet()) {
                SIGMA sigma = entry2.getKey();
                ScaleSpaceCurve scaleSpaceCurve = entry2.getValue();
                
                for (int i = 0; i < scaleSpaceCurve.getT().length; ++i) {
                    
                    float curvature = scaleSpaceCurve.getK(i);
                    
                    if (Math.abs(curvature) < 0.002) {
                        continue;
                    }
                    
                    String str = String.format("(%d,%d) k=%.4f sigma=%.3f", 
                        Math.round(scaleSpaceCurve.getX(i)),
                        Math.round(scaleSpaceCurve.getY(i)), curvature, 
                        SIGMA.getValue(sigma));
                    
                    System.out.println(str);
                }
            }
            System.out.flush();
            break;
        }
    }
    
}
