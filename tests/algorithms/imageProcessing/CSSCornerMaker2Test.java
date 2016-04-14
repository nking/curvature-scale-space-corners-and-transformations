package algorithms.imageProcessing;

import algorithms.imageProcessing.features.PhaseCongruencyDetector;
import algorithms.imageProcessing.scaleSpace.CSSCornerMaker;
import algorithms.misc.MiscDebug;
import algorithms.util.CornerArray;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.util.List;
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
        double tLow = 0.025;
        double tHigh = 0.3;
        boolean increaseKIfNeeded = false;
        PhaseCongruencyDetector phaseCDetector = new PhaseCongruencyDetector();
        PhaseCongruencyDetector.PhaseCongruencyProducts products =
            phaseCDetector.phaseCongMono(img, nScale, minWavelength, mult, 
                sigmaOnf, k, increaseKIfNeeded,
                cutOff, g, deviationGain, noiseMethod, tLow, tHigh);
        assertNotNull(products);
        
        Image tmp = new Image(nCols, nRows);
        for (int i = 0; i < nRows; ++i) {
            for (int j = 0; j < nCols; ++j) {
                if (products.getThinned()[i][j] > 0) {
                    tmp.setRGB(j, i, 255, 255, 255);
                }
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
        
        CSSCornerMaker cornerMaker = new CSSCornerMaker(nCols, nRows);
        List<CornerArray> cornerList = cornerMaker.findCornersInScaleSpaceMaps(theEdges);
       
        for (CornerArray corners : cornerList) {
            for (int i = 0; i < corners.getN(); ++i) {
                ImageIOHelper.addPointToImage(corners.getX(i), corners.getY(i), 
                    tmp, 0, 255, 0, 0);
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
    
}
