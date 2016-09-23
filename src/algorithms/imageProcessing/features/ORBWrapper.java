
package algorithms.imageProcessing.features;

import algorithms.imageProcessing.Image;
import gnu.trove.list.TIntList;

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
        TIntList outputKeypoints0, TIntList outputKeypoints1,
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
        orb.overrideToNotCreateDescriptors();
        if (create2ndDerivPointsAlso) {
            orb.overrideToAlsoCreate2ndDerivKeypoints();
        }
        
        orb.detectAndExtract(subImage);
        
        TIntList kp0 = orb.getAllKeyPoints0();
        TIntList kp1 = orb.getAllKeyPoints1();
        
        for (int i = 0; i < kp0.size(); ++i) {
            int r = kp0.get(i);
            int c = kp1.get(i);
            
            int x = c + startX;
            int y = r + startY;
            if (x >= xLL && x <= xUR && y >= yLL && y <= yUR) {
                outputKeypoints0.add(y);
                outputKeypoints1.add(x);
            }
        }
    }
    
    public static void extractKeypointsFromSubImage(Image img,
        int xLL, int yLL, int xUR, int yUR, int nKeypoints,
        TIntList outputKeypoints0, TIntList outputKeypoints1) {
        
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
                outputKeypoints0.add(y);
                outputKeypoints1.add(x);
            }
        }
    }
}
