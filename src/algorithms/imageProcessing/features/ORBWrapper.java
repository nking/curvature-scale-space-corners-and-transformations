
package algorithms.imageProcessing.features;

import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.features.ORB.Descriptors;
import algorithms.util.PairInt;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TIntList;
import java.util.List;

/**
 * a class to perform various tasks related to using ORB.java,
 * for example, extracting keypoints from only part of
 * an image.
 * 
 * @author nichole
 */
public class ORBWrapper {
    
    public static void extractKeypointsFromSubImage(Image img,
        int xLL, int yLL, int xUR, int yUR, int nKeypoints,
        List<PairInt> outputKeypoints,
        TDoubleList outputOrientations,
        Descriptors outputDescriptors,
        float fastThreshold, boolean create2ndDerivPointsAlso) {
         
        int buffer = 25;
        
        int startX = xLL - buffer;
        if (startX < 0) {
            startX = 0;
        }
        int stopX = xUR + buffer;
        if (stopX > img.getWidth()) {
            stopX = img.getWidth();
        }
        
        int startY = yLL - buffer;
        if (startY < 0) {
            startY = 0;
        }
        int stopY = yUR + buffer;
        if (stopY > img.getHeight()) {
            stopY = img.getHeight();
        }
        
        Image subImage = img.copySubImage(startX, stopX, startY, stopY);
        
        ORB orb = new ORB(nKeypoints);
        orb.overrideFastThreshold(fastThreshold);
        if (create2ndDerivPointsAlso) {
            orb.overrideToAlsoCreate2ndDerivKeypoints();
        }
        
        orb.detectAndExtract(subImage);

        List<PairInt> kp = orb.getAllKeyPoints();
        TDoubleList or = orb.getAllOrientations();
        Descriptors desc = orb.getAllDescriptors();
        
        int nd = desc.descriptors.length;
        int z = 1;
        
        for (int i = 0; i < kp.size(); ++i) {
            PairInt p = kp.get(i);
            int r = p.getY();
            int c = p.getX();
            
            int x = c + startX;
            int y = r + startY;
            if (x >= xLL && x <= xUR && y >= yLL && y <= yUR) {
                outputKeypoints.add(new PairInt(x, y));
                outputOrientations.add(or.get(i));
            }
        }
        
        int[][] outD = new int[outputKeypoints.size()][];
        int count = 0;
        for (int i = 0; i < kp.size(); ++i) {
            PairInt p = kp.get(i);
            int r = p.getY();
            int c = p.getX();
            
            int x = c + startX;
            int y = r + startY;
            if (x >= xLL && x <= xUR && y >= yLL && y <= yUR) {
                outD[count] = desc.descriptors[i];
                count++;
            }
        }
        outputDescriptors.descriptors = outD;
    }
    
    public static void extractKeypointsFromSubImage(Image img,
        int xLL, int yLL, int xUR, int yUR, int nKeypoints,
        List<PairInt> outputKeypoints,
        TDoubleList outputOrientations,
        Descriptors outputDescriptorsH,
        Descriptors outputDescriptorsS,
        Descriptors outputDescriptorsV,
        float fastThreshold, boolean create2ndDerivPointsAlso) {
         
        int buffer = 25;
        
        int startX = xLL - buffer;
        if (startX < 0) {
            startX = 0;
        }
        int stopX = xUR + buffer;
        if (stopX > img.getWidth()) {
            stopX = img.getWidth();
        }
        
        int startY = yLL - buffer;
        if (startY < 0) {
            startY = 0;
        }
        int stopY = yUR + buffer;
        if (stopY > img.getHeight()) {
            stopY = img.getHeight();
        }
        
        Image subImage = img.copySubImage(startX, stopX, startY, stopY);
        
        ORB orb = new ORB(nKeypoints);
        orb.overrideFastThreshold(fastThreshold);
        if (create2ndDerivPointsAlso) {
            orb.overrideToAlsoCreate2ndDerivKeypoints();
        }
        orb.overrideToCreateHSVDescriptors();
        
        orb.detectAndExtract(subImage);

        List<PairInt> kp = orb.getAllKeyPoints();
        TDoubleList or = orb.getAllOrientations();
        Descriptors[] descHSV = orb.getAllDescriptorsHSV();
        
        int nd = kp.size();
        System.out.println("nd=" + nd);
        assert(nd == descHSV[0].descriptors.length);
        assert(nd == descHSV[1].descriptors.length);
        assert(nd == descHSV[2].descriptors.length);
        
        for (int i = 0; i < kp.size(); ++i) {
            PairInt p = kp.get(i);
            int r = p.getY();
            int c = p.getX();
            
            int x = c + startX;
            int y = r + startY;
            if (x >= xLL && x <= xUR && y >= yLL && y <= yUR) {
                outputKeypoints.add(new PairInt(x, y));
                outputOrientations.add(or.get(i));
            }
        }
        
        int[][] outH = new int[outputKeypoints.size()][];
        int[][] outS = new int[outputKeypoints.size()][];
        int[][] outV = new int[outputKeypoints.size()][];
        int count = 0;
        for (int i = 0; i < kp.size(); ++i) {
            PairInt p = kp.get(i);
            int r = p.getY();
            int c = p.getX();
            
            int x = c + startX;
            int y = r + startY;
            if (x >= xLL && x <= xUR && y >= yLL && y <= yUR) {
                outH[count] = descHSV[0].descriptors[i];
                outS[count] = descHSV[1].descriptors[i];
                outV[count] = descHSV[2].descriptors[i];
                count++;
            }
        }
        outputDescriptorsH.descriptors = outH;
        outputDescriptorsS.descriptors = outS;
        outputDescriptorsV.descriptors = outV;
    }
    
    public static void extractKeypointsFromSubImage(Image img,
        int xLL, int yLL, int xUR, int yUR, int nKeypoints,
        List<PairInt> outputKeypoints) {
        
        int buffer = 25;
        
        int startX = xLL - buffer;
        if (startX < 0) {
            startX = 0;
        }
        int stopX = xUR + buffer;
        if (stopX > img.getWidth()) {
            stopX = img.getWidth();
        }
        
        int startY = yLL - buffer;
        if (startY < 0) {
            startY = 0;
        }
        int stopY = yUR + buffer;
        if (stopY > img.getHeight()) {
            stopY = img.getHeight();
        }
        
        Image subImage = img.copySubImage(startX, stopX, startY, stopY);
        
        ORB orb = new ORB(nKeypoints);
        orb.overrideToNotCreateDescriptors();
        orb.detectAndExtract(subImage);
        
        TIntList kp0 = orb.getAllKeyPoints0();
        TIntList kp1 = orb.getAllKeyPoints1();
        
        for (int i = 0; i < kp0.size(); ++i) {
            int r = kp0.get(i);
            int c = kp1.get(i);
            
            int x = c + startX;
            int y = r + startY;
            if (x >= xLL && x <= xUR && y >= yLL && y <= yUR) {
                outputKeypoints.add(new PairInt(x, y));
            }
        }
    }
}
