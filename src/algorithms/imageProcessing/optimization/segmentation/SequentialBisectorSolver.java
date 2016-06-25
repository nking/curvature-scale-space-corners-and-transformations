package algorithms.imageProcessing.optimization.segmentation;

import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageSegmentation;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

/**
 * given training datasets, finds the best combination of
 * the 4 parameters of a segmentation algorithm using
 * a pattern of bisecting each parameter's range space sequentially until 
 * the solution converges.
 * It assumes the there is one global solution and no local minima.
 * 
 * @author nichole
 */
public class SequentialBisectorSolver {
        
    private final Parameter tLen = new Parameter(1, 40, 5);
    
    private final Parameter tColor;
    
    private final Parameter tR = new Parameter(0.5f * 3.0f, 0.95f * 3.0f, 0.5f * 3.0f);
    
    private final Parameter tSmallMerge = new Parameter(0.005f, 0.05f, 0.005f);
    
    private double difference = Double.MAX_VALUE;
    
    private final int colorSpace;
    
    private final boolean reduceNoise;
    
    private SData[] trainingData;
    
    public SequentialBisectorSolver(boolean useHSV, boolean useLowNoiseEdges,
        SData[] trainingData) {
        colorSpace = useHSV ? 1 : 0;
        reduceNoise = useLowNoiseEdges;
        this.trainingData = trainingData;
        
        if (useHSV) {
            tColor = new Parameter(0.05f, 0.7f, 0.05f);
        } else {
            tColor = new Parameter(2.5f, 9.f, 0.5f);
        }
    }
    
    /*
    4 params:                    range of values
          int tLen,                   trying as an integer, but might need to use 1 to 0.1*nPixels.  try: 1 to 100   d=10, n=10
          double tColor,              for clrSpace=0   2.3 to 9   (2.5 to 5.5 is more reasonable)                    d=0.5, n=13
                                      for clrSpace=1   0.1 to 0.7                                                    d=0.05, n=12
          double tR,                  for clrSpace=0   0.5*3.0 to 0.95*3.0                                           d=0.5*3, n=18
                                      for clrSpace=1   0.5*3.0 to 0.95*3.0                                           d=0.5*3, n=18
          double tSmallMerge,         0.01 to 0.1                                                                    d=0.01, n=10

       *solve for using clrSpace = 0 then =1

       *solve for this as true and false:
          boolean reduceNoise,
    */
    
    /**
     * 
     * @throws IOException exceptions are thrown for errors in finding or
     * reading the training file data.
     */
    public double solve() throws IOException, Exception {
        
        Parameter[] parameters = new Parameter[]{
            tLen, tColor, tR, tSmallMerge};
        
        int np = 1 << parameters.length;
    
        boolean hasConverged = false;
      
        SegmentationResults[] expected = readTrainingFiles();
        
        List<List<PairIntArray>> edgesList = extractEdges();
        
        double lastDifference = Double.MAX_VALUE;
        
        while (!hasConverged) {
            
            double minDiff = Double.MAX_VALUE;

            boolean[] minDiffIsLow = new boolean[parameters.length];
            
            for (int p0 = 0; p0 < 2; ++p0) {
                int lowIdx0 = parameters[p0].lowIdx;
                int highIdx0 = parameters[p0].highIdx;
                int midIdx0 = (highIdx0 + lowIdx0) >> 1;
                int tIdx0;
                if (p0 == 0) {
                    tIdx0 = (midIdx0 + lowIdx0) >> 1;
                } else {
                    tIdx0 = (midIdx0 + highIdx0) >> 1;
                }
                for (int p1 = 0; p1 < 2; ++p1) {
                    int lowIdx1 = parameters[p1].lowIdx;
                    int highIdx1 = parameters[p1].highIdx;
                    int midIdx1 = (highIdx1 + lowIdx1) >> 1;
                    int tIdx1;
                    if (p1 == 0) {
                        tIdx1 = (midIdx1 + lowIdx1) >> 1;
                    } else {
                        tIdx1 = (midIdx1 + highIdx1) >> 1;
                    }
                    for (int p2 = 0; p2 < 2; ++p2) {
                        int lowIdx2 = parameters[p2].lowIdx;
                        int highIdx2 = parameters[p2].highIdx;
                        int midIdx2 = (highIdx2 + lowIdx2) >> 1;
                        int tIdx2;
                        if (p2 == 0) {
                            tIdx2 = (midIdx2 + lowIdx2) >> 1;
                        } else {
                            tIdx2 = (midIdx2 + highIdx2) >> 1;
                        }
                        for (int p3 = 0; p3 < 2; ++p3) {
                            int lowIdx3 = parameters[p3].lowIdx;
                            int highIdx3 = parameters[p3].highIdx;
                            int midIdx3 = (highIdx3 + lowIdx3) >> 1;
                            int tIdx3;
                            if (p3 == 0) {
                                tIdx3 = (midIdx3 + lowIdx3) >> 1;
                            } else {
                                tIdx3 = (midIdx3 + highIdx3) >> 1;
                            }
                            double diff = invoke(
                                parameters[0].getValue(tIdx0),
                                parameters[1].getValue(tIdx1),
                                parameters[2].getValue(tIdx2),
                                parameters[3].getValue(tIdx3),
                                expected, edgesList);
                            if (diff < minDiff) {
                                minDiff = diff;
                                minDiffIsLow[0] = (p0 == 0);
                                minDiffIsLow[1] = (p1 == 0);
                                minDiffIsLow[2] = (p2 == 0);
                                minDiffIsLow[3] = (p3 == 0);
                            }
                        }
                    }
                }
            }
            
            for (int pIdx = 0; pIdx < parameters.length; ++pIdx) {
                int lowIdx = parameters[pIdx].lowIdx;
                int highIdx = parameters[pIdx].highIdx;
                int midIdx = (highIdx + lowIdx) >> 1;
                if (minDiffIsLow[pIdx]) {
                    int tIdx = (midIdx + lowIdx) >> 1;
                    if (minDiff == 0) {
                        lowIdx = tIdx;
                        highIdx = tIdx;
                    } else if (highIdx == tIdx) {
                        highIdx--;
                    } else {
                        highIdx = tIdx;
                    }
                } else {
                    int tIdx = (midIdx + highIdx) >> 1;
                    if (minDiff == 0) {
                        lowIdx = tIdx;
                        highIdx = tIdx;
                    } else if (lowIdx == tIdx) {
                        lowIdx++;
                    } else {
                        lowIdx = tIdx;
                    }
                }
                parameters[pIdx].lowIdx = lowIdx;
                parameters[pIdx].highIdx = highIdx;
            }
           
            lastDifference = minDiff;
            
            if (minDiff == 0) {
                break;
            }                              
            // check for convergence
            hasConverged = true;
            for (int pIdx = 0; pIdx < parameters.length; ++pIdx) {
                //if (parameters[pIdx].lowIdx != parameters[pIdx].highIdx) {
                if (parameters[pIdx].lowIdx < parameters[pIdx].highIdx) {
                    hasConverged = false;
                    break;
                }
            }
        }
        
        return lastDifference;
    }

    private double[] invoke(int pIdx, int tIdx1, int tIdx2,
        Parameter[] parameters, SegmentationResults[] expected,
        List<List<PairIntArray>> edgesList) throws Exception {
        
        double[] sumDiffs = new double[2];
        
        switch (pIdx) {
            case 0:
                sumDiffs[0] = invoke(parameters[0].getValue(tIdx1),
                    parameters[1].getMidValue(),
                    parameters[2].getMidValue(),
                    parameters[3].getMidValue(), expected, edgesList);
                sumDiffs[1] = invoke(parameters[0].getValue(tIdx2),
                    parameters[1].getMidValue(),
                    parameters[2].getMidValue(),
                    parameters[3].getMidValue(), expected, edgesList);
                break;
            case 1:
                sumDiffs[0] = invoke(
                    parameters[0].getMidValue(),
                    parameters[1].getValue(tIdx1),
                    parameters[2].getMidValue(),
                    parameters[3].getMidValue(), expected, edgesList);
                sumDiffs[1] = invoke(
                    parameters[0].getMidValue(),
                    parameters[1].getValue(tIdx2),
                    parameters[2].getMidValue(),
                    parameters[3].getMidValue(), expected, edgesList);
                break;
            case 2:
                sumDiffs[0] = invoke(
                    parameters[0].getMidValue(),
                    parameters[1].getMidValue(),
                    parameters[2].getValue(tIdx1),
                    parameters[3].getMidValue(), expected, edgesList);
                sumDiffs[1] = invoke(
                    parameters[0].getMidValue(),
                    parameters[1].getMidValue(),
                    parameters[2].getValue(tIdx2),
                    parameters[3].getMidValue(), expected, edgesList);
                break;
            default:
                sumDiffs[0] = invoke(
                    parameters[0].getMidValue(),
                    parameters[1].getMidValue(),
                    parameters[2].getMidValue(),
                    parameters[3].getValue(tIdx1), expected, edgesList);
                sumDiffs[1] = invoke(
                    parameters[0].getMidValue(),
                    parameters[1].getMidValue(),
                    parameters[2].getMidValue(),
                    parameters[3].getValue(tIdx2), expected, edgesList);
                break;
        }
        
        return sumDiffs;
    }
    
    private SegmentationResults[] readTrainingFiles() throws IOException {
        
        BerkeleySegmentationFileReader reader = new BerkeleySegmentationFileReader();
        
        SegmentationResults[] output = new SegmentationResults[trainingData.length];
        
        for (int i = 0; i < trainingData.length; ++i) {
            
            //String rootName = trainingData[i].imgFileName.split("\\.")[0];
            String segFilePath = trainingData[i].dirPath + "/" + trainingData[i].segFileName;
            
            List<Set<PairInt>> set = reader.readFile(segFilePath);
            
            output[i] = new SegmentationResults(set);
        }
        
        return output;
    }
    
    private List<List<PairIntArray>> extractEdges() throws Exception {
        
        List<List<PairIntArray>> output = new ArrayList<List<PairIntArray>>();
        
        ImageSegmentation imageSegmentation = new ImageSegmentation();
                
        for (int i = 0; i < trainingData.length; ++i) {
            
            String rootName = trainingData[i].imgFileName.split("\\.")[0];
            String imgFilePath = trainingData[i].dirPath + "/" + trainingData[i].imgFileName;        
            
            ImageExt img = ImageIOHelper.readImageExt(imgFilePath);
        
            List<PairIntArray> edges = imageSegmentation.extractEdges(img, 
                reduceNoise, rootName);
            
            output.add(edges);
        }
        
        return output;
    }
    
    private double invoke(float tLenValue, float tColorValue, float tRValue,
        float tSmallMergeValue, SegmentationResults[] expected,
        List<List<PairIntArray>> edgesList) throws Exception {
                
        ImageSegmentation imageSegmentation = new ImageSegmentation();
        
        double sumDifference = 0;
        
        for (int i = 0; i < trainingData.length; ++i) {
            
            String rootName = trainingData[i].imgFileName.split("\\.")[0];
            String imgFilePath = trainingData[i].dirPath + "/" + trainingData[i].imgFileName;        
            
            ImageExt img = ImageIOHelper.readImageExt(imgFilePath);
            
            List<Set<PairInt>> results = 
                imageSegmentation.createColorEdgeSegmentation(img, 
                    edgesList.get(i),
                    colorSpace, Math.round(tLenValue), tColorValue, tRValue, 
                    reduceNoise, tSmallMergeValue, rootName);
            
            SegmentationResults sr0 = new SegmentationResults(results);
            
            SegmentationResults exp = expected[i];
            
            sumDifference += sr0.evaluate(exp);
        }
        
        return sumDifference;
    }

    public float[] getParameters() {
        
        float[] output = new float[4];
        output[0] = tLen.getMidValue();
        output[1] = tColor.getMidValue();
        output[2] = tR.getMidValue();
        output[3] = tSmallMerge.getMidValue();
        
        return output;
    }

    private static class Parameter {
        final float vFirst;
        final float vLast;
        final float deltaV;
        int lowIdx = 0;
        int highIdx;
        
        public Parameter(float startRange, float stopRange, float delta) {
            vFirst = startRange;
            vLast = stopRange;
            deltaV = delta;
            highIdx = Math.round((stopRange - startRange)/delta);
            //0.1 to 0.7     d=0.05, n=13
            // 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7
            //  0   1    2    3   4    5    6   7   8   9   10   11   12
        }
        
        public float getMidValue() {
            int m = (lowIdx + highIdx) >> 1;
            float v = (m * deltaV) + vFirst;
            return v;
        }
        
        public float getValue(int idx) {
            float v = (idx * deltaV) + vFirst;
            return v;
        }
    }
}
