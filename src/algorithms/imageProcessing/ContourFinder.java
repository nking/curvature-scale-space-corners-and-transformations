package algorithms.imageProcessing;

import algorithms.util.PairIntArray;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.logging.Logger;

/**
 * class to identity contours in scale space images.
 * 
 * @author nichole
 */
public class ContourFinder {
    
    private double thresholdFactor = 0.15;
        
    protected Logger log = null;
    
    public ContourFinder() {
        log = Logger.getLogger(this.getClass().getName());
    }
    
    /**
     * find contours in this scale space map for an edge of given edgeNumber.
     * Note that the edgeNumber is not used, but is kept for use with 
     * the indexes in debugging later.
     * 
     * @param scaleSpaceImage
     * @param edgeNumber
     * @return 
     */
    public List<CurvatureScaleSpaceContour> findContours(
        ScaleSpaceCurveImage scaleSpaceImage, int edgeNumber) {
        
        List<CurvatureScaleSpaceContour> contours = new
            ArrayList<CurvatureScaleSpaceContour>();
        
        if ((scaleSpaceImage == null) || 
            (scaleSpaceImage.getImageSigmas().length == 0)) {
            return contours;
        }
        
        ScaleSpaceCurveImage space = scaleSpaceImage.copy();
        
        double lowLimit = space.getImageSigmas()[0]*thresholdFactor;
        //TODO: consider a low limit of sigma=3
        if (lowLimit < 3) {
            lowLimit = 3;
        }
                
        // find the first contour at this height and extract it from the
        // dataset, nulling
        for (int i = 0; i < space.getImageSigmas().length; i++) {
            
            float sigma = space.getImageSigmas()[i];
            
            if (sigma < lowLimit) {
                break;
            }
            
            boolean extract = true;
            
            while (extract) {
                
                CurvatureScaleSpaceContour contour = extractNextContour(
                    scaleSpaceImage, i);
                
                if (contour != null) {
                    
                    contours.add(contour);
                    
                } else {
                    extract = false;
                }
            }
        }
        
        correctForWrappedContours(contours);
            
        return contours;
    }

    private CurvatureScaleSpaceContour extractNextContour(
        ScaleSpaceCurveImage scaleSpaceImage, int sigmaIndex) {
        
        if ((scaleSpaceImage == null) || 
            (scaleSpaceImage.getScaleSpaceImage().length == 0)) {
            return null;
        }
        
        float[] t = scaleSpaceImage.getScaleSpaceImage()[sigmaIndex];
        
        if (t == null || t.length == 0) {
            return null;
        }
                
        for (int i = 0; i < t.length; i++) {
            
            if (t[i] < 0) {
                continue;
            }
            
            CurvatureScaleSpaceContour contour = extractContour(scaleSpaceImage,
                sigmaIndex, i);
            
            if (contour != null) {
                
                return contour;
            }
        }
        
        return null;
    }

    private CurvatureScaleSpaceContour extractContour(ScaleSpaceCurveImage 
        scaleSpaceImage, int sigmaIndex, int tIndex) {
  
        if ((scaleSpaceImage == null) || 
            (scaleSpaceImage.getScaleSpaceImage().length == 0)) {
            return null;
        }
        
        float tPoint = scaleSpaceImage.getScaleSpaceImage()[sigmaIndex][tIndex];
        if (tPoint < 0) {
            return null;
        }
        
        int nToRight = 0;

        for (int i = (tIndex + 1); 
            i < scaleSpaceImage.getScaleSpaceImage()[sigmaIndex].length;
            i++) {

            if (scaleSpaceImage.getScaleSpaceImage()[sigmaIndex][i] >= 0) {
                nToRight++;
            }
        }
        
        float sigma = scaleSpaceImage.getImageSigmas()[sigmaIndex];

        if (nToRight == 0) {
            
            // single peak contour if the value is larger than zero            
            
            float t = scaleSpaceImage.getScaleSpaceImage()[sigmaIndex][tIndex];
            
            if (t < 0) {
                // this has alredy been extracted, so return null
                return null;
            }
            
            CurvatureScaleSpaceContour contour = new CurvatureScaleSpaceContour(
                sigma, t);
            
            contour.setEdgeNumber(scaleSpaceImage.getEdgeNumber());
                        
            float t0 = scaleSpaceImage.getScaleSpaceImage()[sigmaIndex][tIndex];
            int idx0 = Math.round(t0 * scaleSpaceImage.getEdgeSize());
            
            CurvatureScaleSpaceImagePoint point0 = 
                new CurvatureScaleSpaceImagePoint(sigma, t,
                    scaleSpaceImage.getXCoord(sigmaIndex, tIndex),
                    scaleSpaceImage.getYCoord(sigmaIndex, tIndex),
                    idx0);
            
            CurvatureScaleSpaceImagePoint[] peakPoints = 
                new CurvatureScaleSpaceImagePoint[]{point0};                
            
            contour.setPeakDetails(peakPoints);
            
            // for case when there's a single point for the peak:
            removeContourFromImage(scaleSpaceImage, sigmaIndex, tIndex);

            return contour;
        }
        
        /*
        Find the next non-negative value in scaleSpaceImage for sigma
        and determine where it's right branch is if any.
        
        For now, will assume that there is never an embedded contour
        which is starting, that is peaking at this same sigma level.
        
        Will look for the first non-negative to be a single peak
        or the left of a left and right of a peak.
        */
       
        boolean isASinglePeak = false;
        
        if (sigmaIndex == (scaleSpaceImage.getImageSigmas().length - 1)) {
            
            // if there's a -1 to the right it's a single point, else, it
            // may be left and right branch or it may not.  doesn't matter
            // very much because this is the bottom sigma of the image
            float[] t = scaleSpaceImage.getScaleSpaceImage()[sigmaIndex];
            if (tIndex == (t.length - 1)) {
                isASinglePeak = true;
            } else if (t[tIndex + 1] == -1) {
                isASinglePeak = true;
            }
            
        } else {
            // descend one level to see if there are 2 peaks 
            // under the current peak that are left and right of it.
            // if the right is closer than the next point on this same level,
            // the current point is a peak
            int leftIndexBelow = -1;
            int rightIndexBelow = -1;
            float minDiffLeftBelow = Float.MAX_VALUE;
            float minDiffRightBelow = Float.MAX_VALUE;


            float[] t = scaleSpaceImage.getScaleSpaceImage()[sigmaIndex + 1];
            // sometimes for low sigma, the contours are misshapen and have
            // vertical gaps, so iterating now to levels below to assert have
            // the contour values below this peak or partial peak
            if (sigma < 3.5) {
                int t0 = (tIndex > 0) ? (tIndex - 1) : tIndex;
                int t1 = ((tIndex + 1) < t.length) ? tIndex + 1 : tIndex;
                if (t1 > t0) {
                    int si = sigmaIndex + 1;
                    while (si < (scaleSpaceImage.getScaleSpaceImage().length - 1)) {
                        t = scaleSpaceImage.getScaleSpaceImage()[si + 1];
                        int n = 0;
                        for (int j = t0; j < t.length; j++) {
                            float tt = t[j];
                            if (tt > -1) {
                                n++;
                            }
                        }
                        if (n > 1) {
                            break;
                        }
                        si++;
                    }
                }
            }

            int start = (tIndex > 0) ? (tIndex - 1) : tIndex;
            for (int j = start; j < t.length; j++) {
                if (t[j] < 0) {
                    continue;
                }
                float lD = tPoint - t[j];
                float rD = t[j] - tPoint;

                if ((lD >= 0) && (lD < minDiffLeftBelow)) {
                    minDiffLeftBelow = lD;
                    leftIndexBelow = j;
                }

                if ((rD >= 0) && (rD < minDiffRightBelow) &&
                    (j > leftIndexBelow)) {
                    minDiffRightBelow = rD;
                    rightIndexBelow = j;
                }
            }

            boolean isAnEdgePair = false;
            
            if (rightIndexBelow == -1) {
                //TODO: the 2nd conditional should be revised
                if ((nToRight == 1) && (tPoint > 0.9)) {
                    isAnEdgePair = true;
                } else {
                    // this can happen if the contour has already been removed
                    return null;
                }
            }
            
            if (!isAnEdgePair) {
                
                float tNext = 
                    scaleSpaceImage.getScaleSpaceImage()[sigmaIndex][tIndex + 1];
                
                isASinglePeak = (t[rightIndexBelow] < tNext) || (tNext < 0);
            }
        }
     
        if (isASinglePeak) {
    
            // it's a single peak
            
            CurvatureScaleSpaceContour contour = new CurvatureScaleSpaceContour(
                sigma, tPoint);
            
            contour.setEdgeNumber(scaleSpaceImage.getEdgeNumber());
                        
            float t = scaleSpaceImage.getScaleSpaceImage()[sigmaIndex][tIndex];
            
            float t0 = scaleSpaceImage.getScaleSpaceImage()[sigmaIndex][tIndex];
            int idx0 = Math.round(t0 * scaleSpaceImage.getEdgeSize());
            
            CurvatureScaleSpaceImagePoint point0 = 
                new CurvatureScaleSpaceImagePoint(sigma, t,
                    scaleSpaceImage.getXCoord(sigmaIndex, tIndex),
                    scaleSpaceImage.getYCoord(sigmaIndex, tIndex), 
                    idx0);
            
            CurvatureScaleSpaceImagePoint[] peakPoints = 
                new CurvatureScaleSpaceImagePoint[]{point0};                
            
            contour.setPeakDetails(peakPoints);
                        
            // for case when there's a single point for the peak:
            removeContourFromImage(scaleSpaceImage, sigmaIndex, tIndex);
                        
            return contour;
            
        } else {
            
            // else its the left side of a left and right point which are the peak
               
            float tNext = 
                scaleSpaceImage.getScaleSpaceImage()[sigmaIndex][tIndex + 1];
    
            CurvatureScaleSpaceContour contour = new CurvatureScaleSpaceContour(
                sigma, (tPoint + tNext)/2.f);
            
            contour.setEdgeNumber(scaleSpaceImage.getEdgeNumber());

            float t0 = scaleSpaceImage.getScaleSpaceImage()[sigmaIndex][tIndex];
            float t1 = scaleSpaceImage.getScaleSpaceImage()[sigmaIndex][tIndex + 1];
            
            int idx0 = Math.round(t0 * scaleSpaceImage.getEdgeSize());
            int idx1 = Math.round(t1 * scaleSpaceImage.getEdgeSize());
            
            CurvatureScaleSpaceImagePoint point0 = 
                new CurvatureScaleSpaceImagePoint(sigma, t0,
                scaleSpaceImage.getXCoord(sigmaIndex, tIndex),
                scaleSpaceImage.getYCoord(sigmaIndex, tIndex), idx0);
            
            CurvatureScaleSpaceImagePoint point1 = 
                new CurvatureScaleSpaceImagePoint(sigma, t1,
                scaleSpaceImage.getXCoord(sigmaIndex, tIndex + 1),
                scaleSpaceImage.getYCoord(sigmaIndex, tIndex + 1), idx1);
            
            CurvatureScaleSpaceImagePoint[] peakPoints = 
                new CurvatureScaleSpaceImagePoint[]{point0, point1};       
            
            contour.setPeakDetails(peakPoints);
            
            removeContourFromImage(scaleSpaceImage, sigmaIndex, tIndex, 
                tIndex + 1);
                        
            return contour;
        }
        
    }

    /**
     * remove the contour under a peak.  Note that the method does not yet
     * handle complex morphologies, such as embedded contours right under
     * the peak;
     * 
     * @param scaleSpaceImage
     * @param sigmaIndex
     * @param tIndex 
     */
    private void removeContourFromImage(ScaleSpaceCurveImage scaleSpaceImage, 
        int sigmaIndex, int tIndex) {
        
        if ((sigmaIndex > (scaleSpaceImage.getScaleSpaceImage().length - 1)) ||
            (scaleSpaceImage.getScaleSpaceImage()[sigmaIndex] == null) ||
            (scaleSpaceImage.getScaleSpaceImage()[sigmaIndex].length == 0)) {
            return;
        }

        if ((tIndex > 
            (scaleSpaceImage.getScaleSpaceImage()[sigmaIndex].length - 1))
            ) {
            return;
        }
        
        float tLeft = scaleSpaceImage.getScaleSpaceImage()[sigmaIndex][tIndex];
        float tRight = tLeft;
        int leftIndex = -1;
        int rightIndex = -1;
            
        for (int i = sigmaIndex; i < (sigmaIndex + 1); i++) {
            
            float[] t = scaleSpaceImage.getScaleSpaceImage()[i];
            
            float tPoint = t[tIndex];
            
            leftIndex = -1;
            rightIndex = -1;
            float minDiffLeft = Float.MAX_VALUE;
            float minDiffRight = Float.MAX_VALUE;
            for (int j = tIndex; j < t.length; j++) {
                
                if (t[j] < 0) {
                    continue;
                }
                
                float lD = tPoint - t[j];
                float rD = t[j] - tPoint;
                
                if ((lD >= 0) && (lD < minDiffLeft)) {
                    minDiffLeft = lD;
                    leftIndex = j;
                }
                
                if ((rD >= 0) && (rD < minDiffRight) && (j > leftIndex)) {
                    minDiffRight = rD;
                    rightIndex = j;
                }
            }
            
            if (leftIndex > -1) {
                tLeft = t[leftIndex];
                t[leftIndex] = -1;
            }
            if (rightIndex > -1) {
                tRight = t[rightIndex];
                t[rightIndex] = -1;
            }
        
            boolean isEmpty = true;
            for (float tt : t) {
                if (!(tt < 0)) {
                    isEmpty = false;
                    break;
                }
            }
            if (isEmpty) {
                scaleSpaceImage.getScaleSpaceImage()[sigmaIndex] = new float[0];
            }
        }
        
        if (rightIndex == -1) {
            rightIndex = leftIndex;
        }
        
        removeContourFromImage(scaleSpaceImage, sigmaIndex + 1, leftIndex, 
            rightIndex);
    }
   
    private void removeContourFromImage(ScaleSpaceCurveImage scaleSpaceImage, 
        int sigmaIndex, int tLeftIndex, int tRightIndex) {
   
        if ((scaleSpaceImage.getScaleSpaceImage() == null) ||
            (sigmaIndex > (scaleSpaceImage.getScaleSpaceImage().length - 1)) ||
            (scaleSpaceImage.getScaleSpaceImage()[sigmaIndex] == null) ||
            (scaleSpaceImage.getScaleSpaceImage()[sigmaIndex].length == 0)) {
            return;
        }
        
        if ((tLeftIndex > 
            (scaleSpaceImage.getScaleSpaceImage()[sigmaIndex].length - 1)) ||
            (tRightIndex > 
            (scaleSpaceImage.getScaleSpaceImage()[sigmaIndex].length - 1))
            ) {
            return;
        }
        
        float tLeft = scaleSpaceImage.getScaleSpaceImage()[sigmaIndex][tLeftIndex];
        float tRight = scaleSpaceImage.getScaleSpaceImage()[sigmaIndex][tRightIndex];
        
        for (int i = sigmaIndex; i < scaleSpaceImage.getImageSigmas().length; i++) {
            
            float[] t = scaleSpaceImage.getScaleSpaceImage()[i];
                        
            int leftIndex = -1;
            int rightIndex = -1;
            float minDiffLeft = Float.MAX_VALUE;
            float minDiffRight = Float.MAX_VALUE;
            
            for (int j = 0; j < t.length; j++) {
                
                if (t[j] < 0) {
                    continue;
                }
                
                float lD = tLeft - t[j];
                float rD = t[j] - tRight;
                
                if (lD < 0) {
                    lD *= -1;
                }
                if (rD < 0) {
                    rD *= -1;
                }
                
                if (lD < minDiffLeft) {
                    minDiffLeft = lD;
                    leftIndex = j;
                }
                
                if ((rD < minDiffRight) && (j > leftIndex)) {
                    minDiffRight = rD;
                    rightIndex = j;
                }
            }
            
            // TODO: improve correction for wrap around.
            // Also, this only includes peaks found at end, not beginning.
            if ((leftIndex > -1) && (leftIndex == (t.length - 1)) 
                && (t[leftIndex] >= 0.9) &&
                (t[0] < 0.1)) {
                rightIndex = 0;
            } 
            
            if (leftIndex > -1) {
                tLeft = t[leftIndex];
                t[leftIndex] = -1;
            }
            if (rightIndex > -1) {
                tRight = t[rightIndex];
                t[rightIndex] = -1;
            }
            
            boolean isEmpty = true;
            for (float tt : t) {
                if (!(tt < 0)) {
                    isEmpty = false;
                    break;
                }
            }
            if (isEmpty) {
                scaleSpaceImage.getScaleSpaceImage()[sigmaIndex] = new float[0];
            }
        }
    }

    /**
     * looks for contours that may be wrap around contours that started near
     * 1.0 and finished on the other side of zero, and corrects for that.
     * Note that this could be done more correctly before the left and right
     * branches are removed from the scale space image, but a correction at 
     * this stage instead of that earlier stage is simpler and easier to 
     * maintain.
     * @param contours 
     */
    private void correctForWrappedContours(final List<CurvatureScaleSpaceContour> 
        contours) {

        if ((contours == null) || contours.isEmpty()) {
            return;
        }
        
        // roughly, look for features with peaks > 0.9 and < 0.1.
        
        // TODO: what shape would produce the widest possible contour in
        // this space?
        
        List<Integer> rightBorderPeakIndexes = new ArrayList<Integer>();
        List<Integer> leftBorderPeakIndexes = new ArrayList<Integer>();
        
        for (int i = 0; i < contours.size(); i++) {
            
            CurvatureScaleSpaceContour contour = contours.get(i);
            
            if (contour.getPeakScaleFreeLength() > 0.9) {
                rightBorderPeakIndexes.add(Integer.valueOf(i));
            } else if (contour.getPeakScaleFreeLength() < 0.1) {
                leftBorderPeakIndexes.add(Integer.valueOf(i));
            }
        }
        
        int maxIter = rightBorderPeakIndexes.size();
        if (leftBorderPeakIndexes.size() > maxIter) {
            maxIter = leftBorderPeakIndexes.size();
        }
        int nIter = 0;
        int i = 0;
        boolean resort = false;
        
        while ((nIter < maxIter) && !rightBorderPeakIndexes.isEmpty() 
            && !leftBorderPeakIndexes.isEmpty()) {
            
            // indexes are ordered by descending sigma
            
            // for now, make an unsafe assumption that there aren't any other
            // full contours within the 0.1 border regions in between the sigma
            // of contours that wrap around
            
            //if ((leftBorderPeakIndexes.size() == 1) 
            //    && (rightBorderPeakIndexes.size() == 1)) {
                
                int idxLeft = leftBorderPeakIndexes.get(i);
                int idxRight = rightBorderPeakIndexes.get(i);

                CurvatureScaleSpaceContour left = contours.get(idxLeft);                
                CurvatureScaleSpaceContour right = contours.get(idxRight);
                
                boolean leftIsTaller = (left.getPeakSigma() 
                    > right.getPeakSigma());
                
                boolean rightIsTaller = (right.getPeakSigma() 
                    > left.getPeakSigma());
                
                if (leftIsTaller && (left.getPeakDetails().length == 2)) {
                    
                    contours.remove(right);
                    resort = true;
                } else if (rightIsTaller && (right.getPeakDetails().length == 2)) {
                    
                    contours.remove(left);
                    resort = true;
                } else {
                    
                    if (leftIsTaller && (right.getPeakDetails().length == 2)) {
                        // left has 1 peak
                        contours.remove(right);
                        resort = true;
                    } else if (rightIsTaller && (left.getPeakDetails().length == 2)) {
                        // right has 1 peak
                        contours.remove(left);
                        resort = true;
                    } else if (!rightIsTaller && !leftIsTaller && 
                        (right.getPeakDetails().length == 2) &&
                        (left.getPeakDetails().length == 2)) {
                        // do nothing, both should remain
                    } else {
                        // both have a single peak, so avg in sigma and t
                        float sigma = (left.getPeakSigma() 
                            + right.getPeakSigma())/2.f;
                        float scaleFreeLength = (left.getPeakScaleFreeLength()
                            + right.getPeakScaleFreeLength())/2.f;
                        
                        CurvatureScaleSpaceContour contour = 
                            new CurvatureScaleSpaceContour(sigma, scaleFreeLength);
                        
                        contour.setEdgeNumber(left.getEdgeNumber());
                        
                        CurvatureScaleSpaceImagePoint[] peakPoints
                            = new CurvatureScaleSpaceImagePoint[]{
                                left.getPeakDetails()[0], 
                                right.getPeakDetails()[0]};
                       
                        contour.setPeakDetails(peakPoints);
                        
                        contours.set(idxLeft, contour);
                        
                        contours.remove(left);
                    
                        contours.remove(right);
                        
                        resort = true;
                    }
                }
            //}
            
            leftBorderPeakIndexes.clear();
            rightBorderPeakIndexes.clear();
            for (int ii = 0; ii < contours.size(); ii++) {

                CurvatureScaleSpaceContour contour = contours.get(ii);

                if (contour.getPeakScaleFreeLength() > 0.9) {
                    rightBorderPeakIndexes.add(Integer.valueOf(ii));
                } else if (contour.getPeakScaleFreeLength() < 0.1) {
                    leftBorderPeakIndexes.add(Integer.valueOf(ii));
                }
            }
            
            nIter++;
        }
        
       if (resort) {
           Collections.sort(contours, new DescendingSigmaComparator());
       }
    }
   
    public boolean reverseIfClockwise(List<CurvatureScaleSpaceContour> result) {
        
        boolean didReverse = false;
        
        PairIntArray testContour = new PairIntArray();
        for (int j = 0; j < result.size(); j++) {
            CurvatureScaleSpaceContour c = result.get(j);
            CurvatureScaleSpaceImagePoint[] points = c.getPeakDetails();
            for (int jj = 0; jj < points.length; jj++) {
                testContour.add(points[jj].getXCoord(), points[jj].getYCoord());
            }
        }
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        boolean isCW = curveHelper.curveIsOrderedClockwise(testContour);
        if (isCW) {
            didReverse = true;
            for (int j = 0; j < result.size(); j++) {
                CurvatureScaleSpaceContour contour = result.get(j);
                CurvatureScaleSpaceContour reversed = 
                    new CurvatureScaleSpaceContour(contour.getPeakSigma(), 
                        1.0f - contour.getPeakScaleFreeLength());
                CurvatureScaleSpaceImagePoint[] points = contour.getPeakDetails();
                if (points.length > 1) {
                    CurvatureScaleSpaceImagePoint tmp = points[0];
                    points[0] = points[1];
                    points[1] = tmp;
                }
                for (int jj = 0; jj < points.length; jj++) {
                    points[jj].setScaleFreeLength(1.0f - points[jj].getScaleFreeLength());
                }
                reversed.setPeakDetails(points);
                reversed.setEdgeNumber(contour.getEdgeNumber());
                result.set(j, reversed);
            }
        }
        
        return didReverse;
    }
}
