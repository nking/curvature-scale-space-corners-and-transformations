package algorithms.imageProcessing.optimization.segmentation;

import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageSegmentation;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
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
public class DownhillSimplexSolver {
        
    private final Parameter tLen = new Parameter(1, 40, 5);
    
    private final Parameter tColor;
    
    private final Parameter tR = new Parameter(0.5f * 3.0f, 0.95f * 3.0f, 0.5f * 3.0f);
    
    private final Parameter tSmallMerge = new Parameter(0.005f, 0.05f, 0.005f);
    
    private double difference = Double.MAX_VALUE;
    
    private final int colorSpace;
    
    private final boolean reduceNoise;
    
    private SData[] trainingData;
  
    public DownhillSimplexSolver(boolean useHSV, 
        boolean useLowNoiseEdges, SData[] trainingData) {
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
    public SFit solve() throws IOException, Exception {
        
        //Parameter[] parameters = new Parameter[]{tLen, tColor, 
        //    tR, tSmallMerge};
     
        SegmentationResults[] expected = readTrainingFiles();
        
        List<List<PairIntArray>> edgesList = extractEdges();
        
        SFit[] yFits = createStarterPoints(expected,
            extractEdges());
        
        int tLenCurrent;
        double tColorCurrent;
        double tRCurrent;
        double tSmallMergeCurrent;
              
        /*
        solve for the best parameters for each then average results
        or solve for best parameters in total.
        */
        float alpha = 1;   // > 0
        float gamma = 2;   // > 1
        float beta = 0.9f; // 0 < beta < 1

        boolean go = true;

        int nMaxIter = 100;
        int nIter = 0;
        
        SFit prevBest = null;

        while (go && (nIter < nMaxIter)) {
            
            yFits = condense(yFits);
            
            if (yFits.length < 3) {
                break;
            }
            
            sortFromMinToMax(yFits, 0, 2);

            int minSimplexIndex = 0;
            int maxSimplexIndex = 2;
            int midSimplexIndex = 1;

            // determine center for all points excepting the worse fit
            double sumTLen = 0.0;
            double sumTColor = 0.0;
            double sumTR = 0.0;
            double sumTSmallMerge = 0.0;
            for (int i = 0; i < yFits.length - 1; i++) {
                sumTLen += yFits[i].tLenF;
                sumTColor += yFits[i].tColorF;
                sumTR += yFits[i].tRF;
                sumTSmallMerge += yFits[i].tSmallMergeF;
            }
            tLenCurrent = (int)Math.round(
                sumTLen / (yFits.length - 1));
            tColorCurrent = sumTColor / (yFits.length - 1);
            tRCurrent = sumTR / (yFits.length - 1);
            tSmallMergeCurrent = sumTSmallMerge / (yFits.length - 1);
            
            // "Reflection"
            double tLenReflect = tLenCurrent - 
                (alpha * (yFits[maxSimplexIndex].tLenF - 
                tLenCurrent));
            double tColorReflect = tColorCurrent - 
                (alpha * (yFits[maxSimplexIndex].tColorF - 
                tColorCurrent));
            double tRReflect = tRCurrent - 
                (alpha * (yFits[maxSimplexIndex].tRF - 
                tRCurrent));
            double tSmallMergeReflect = tSmallMergeCurrent - 
                (alpha * (yFits[maxSimplexIndex].tSmallMergeF - 
                tSmallMergeCurrent));
            SFit yFitReflected = invoke(tLenReflect, tColorReflect, 
                tRReflect, tSmallMergeReflect, expected, edgesList);
            
            if ((yFitReflected.costF < 
                yFits[minSimplexIndex].costF)
                && 
                ((tLenReflect >= tLen.vFirst) 
                    && (tLenReflect <= tLen.vLast) 
                && (tColorReflect >= tColor.vFirst) && 
                    (tColorReflect <= tColor.vLast)
                && (tRReflect >= tR.vFirst) 
                    && (tRReflect <= tR.vLast)
                && (tSmallMergeReflect >= tSmallMerge.vFirst) 
                    && (tSmallMergeReflect <= tSmallMerge.vLast)
                )
            ) {

                // "Expansion"
                double tLenExpansion = tLenReflect - 
                    (gamma * (tLenCurrent - tLenReflect));
                double tColorExpansion = tColorReflect - 
                    (gamma * (tColorCurrent - tColorReflect));
                double tRExpansion = tRReflect - 
                    (gamma * (tRCurrent - tRReflect));
                double tSmallMergeExpansion = tSmallMergeReflect - 
                    (gamma * (tSmallMergeCurrent - tSmallMergeReflect));
                SFit yFitExpansion = invoke(
                    tLenExpansion, tColorExpansion, 
                    tRExpansion, tSmallMergeExpansion, 
                    expected, edgesList);
                
                if ((yFitExpansion.costF < 
                yFits[minSimplexIndex].costF)
                && 
                ((tLenExpansion >= tLen.vFirst) 
                    && (tLenExpansion <= tLen.vLast) 
                && (tColorExpansion >= tColor.vFirst) && 
                    (tColorExpansion <= tColor.vLast)
                && (tRExpansion >= tR.vFirst) 
                    && (tRExpansion <= tR.vLast)
                && (tSmallMergeExpansion >= tSmallMerge.vFirst)
                    && (tSmallMergeExpansion <= tSmallMerge.vLast)
                )
            ) {
                    yFits[maxSimplexIndex] = yFitExpansion;
                } else {
                    yFits[maxSimplexIndex] = yFitReflected;
                }

            } else if (
                (yFitReflected.costF > 
                yFits[midSimplexIndex].costF)
                && 
                ((tLenReflect >= tLen.vFirst) 
                    && (tLenReflect <= tLen.vLast) 
                && (tColorReflect >= tColor.vFirst) && 
                    (tColorReflect <= tColor.vLast)
                && (tRReflect >= tR.vFirst) 
                    && (tRReflect <= tR.vLast)
                && (tSmallMergeReflect >= tSmallMerge.vFirst) 
                    && (tSmallMergeReflect <= tSmallMerge.vLast)
                )
            ) {

                if ((yFitReflected.costF 
                    <= yFits[maxSimplexIndex].costF)
                    && 
                    ((tLenReflect >= tLen.vFirst) 
                        && (tLenReflect <= tLen.vLast) 
                    && (tColorReflect >= tColor.vFirst) && 
                        (tColorReflect <= tColor.vLast)
                    && (tRReflect >= tR.vFirst) 
                        && (tRReflect <= tR.vLast)
                    && (tSmallMergeReflect >= tSmallMerge.vFirst) 
                        && (tSmallMergeReflect <= tSmallMerge.vLast)
                )
                ) {

                    yFits[maxSimplexIndex] = yFitReflected;
                }

                // "Contraction"
                double tLenContraction =     
                    (beta * yFits[maxSimplexIndex].tLenF)     
                    + (1 - beta) * tLenCurrent;
                double tColorContraction =     
                    (beta * yFits[maxSimplexIndex].tColorF)     
                    + (1 - beta) * tColorCurrent;
                double tRContraction =     
                    (beta * yFits[maxSimplexIndex].tRF)     
                    + (1 - beta) * tRCurrent;
                double tSmallMergeContraction =     
                    (beta * yFits[maxSimplexIndex].tSmallMergeF)     
                    + (1 - beta) * tSmallMergeCurrent;
                SFit yFitContraction = invoke(
                    tLenContraction, tColorContraction, 
                    tRContraction, tSmallMergeContraction, 
                    expected, edgesList);
                if (yFitContraction.costF > yFits[maxSimplexIndex].costF
                    && 
                    ((tLenContraction >= tLen.vFirst) 
                    && (tLenContraction <= tLen.vLast) 
                    && (tColorContraction >= tColor.vFirst) && 
                        (tColorContraction <= tColor.vLast)
                    && (tRContraction >= tR.vFirst) 
                        && (tRContraction <= tR.vLast)
                    && (tSmallMergeContraction >= tSmallMerge.vFirst)
                        && (tSmallMergeContraction <= tSmallMerge.vLast)
                    )
                ) {
                    double tLenTmp = (yFits[midSimplexIndex].tLenF 
                        + yFits[minSimplexIndex].tLenF)/2;
                    double tColorTmp = (yFits[midSimplexIndex].tColorF 
                        + yFits[minSimplexIndex].tColorF)/2;
                    double tRTmp = (yFits[midSimplexIndex].tRF 
                        + yFits[minSimplexIndex].tRF)/2;
                    double tSmallMergeTmp = 
                        (yFits[midSimplexIndex].tSmallMergeF 
                        + yFits[minSimplexIndex].tSmallMergeF)/2;
                    yFits[midSimplexIndex] = 
                        invoke(tLenTmp, tColorTmp, 
                            tRTmp, tSmallMergeTmp, expected, edgesList);
                } else {
                    yFits[maxSimplexIndex] = yFitContraction;
                }

            } else if ((tLenReflect >= tLen.vFirst)
                && (tLenReflect <= tLen.vLast)
                && (tColorReflect >= tColor.vFirst)
                && (tColorReflect <= tColor.vLast)
                && (tRReflect >= tR.vFirst)
                && (tRReflect <= tR.vLast)
                && (tSmallMergeReflect >= tSmallMerge.vFirst)
                && (tSmallMergeReflect <= tSmallMerge.vLast))
             {
                yFits[maxSimplexIndex] = yFitReflected;
            }

            if ((tLenCurrent < tLen.vFirst)
                && (tLenCurrent > tLen.vLast)
                && (tColorCurrent < tColor.vFirst)
                && (tColorCurrent > tColor.vLast)
                && (tRCurrent < tR.vFirst)
                && (tRCurrent > tR.vLast)
                && (tSmallMergeCurrent < tSmallMerge.vFirst)
                && (tSmallMergeCurrent > tSmallMerge.vLast))
             {
                go = false;
            }
            
            if (hasConverged(prevBest, yFits[0])) {
                go = false;
            } else {
                prevBest = yFits[0].copy();
            }
            
            if ((nIter % 1) == 0) {
                System.out.println("nIter=" + nIter + " "
                    + yFits[0].toString());
            }
            
            nIter++;
        }

        System.out.println("nIter=" + nIter + " bestFit=" +
            yFits[0].toString());
        
        return yFits[0];
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
    
    private SFit invoke(double tLenValue, double tColorValue, 
        double tRValue, double tSmallMergeValue, 
        SegmentationResults[] expected,
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
                    colorSpace, (int)Math.round(tLenValue), tColorValue, tRValue, 
                    reduceNoise, tSmallMergeValue, rootName);
            
            SegmentationResults sr0 = new SegmentationResults(results);
            
            SegmentationResults exp = expected[i];
            
            sumDifference += sr0.calculateDifference(exp);
        }
   
        SFit sFit = new SFit();
        sFit.tLenF = (int) Math.round(tLenValue);
        sFit.tColorF = tColorValue;
        sFit.tRF = tRValue;
        sFit.tSmallMergeF = tSmallMergeValue;
        sFit.colorSpaceF = this.colorSpace;
        sFit.reduceNoiseF = this.reduceNoise;
        sFit.costF = sumDifference;
      
        return sFit;
    }

    public float[] getParameters() {
        
        float[] output = new float[4];
        output[0] = tLen.getMidValue();
        output[1] = tColor.getMidValue();
        output[2] = tR.getMidValue();
        output[3] = tSmallMerge.getMidValue();
        
        return output;
    }

    private void sortFromMinToMax(SFit[] yFits, int p, int r) {

        if (p < r) {

            int q = partition(yFits, p, r);

            sortFromMinToMax(yFits, p, q - 1);

            sortFromMinToMax(yFits, q + 1, r);
        }
    }
    
    int partition(SFit[] yFits, int p, int r) {

        double xxp = yFits[r].costF;

        int i = p - 1;

        for (int j = p; j < r ; j++ ) {
            
            if (yFits[j].costF <= xxp) {

                i++;
                SFit swap = yFits[i];
                yFits[i] = yFits[j];
                yFits[j] = swap;
            }
        }
        i++;
        SFit swap = yFits[i];
        yFits[i] = yFits[r];
        yFits[r] = swap;

        return i;
    }

    private SFit[] createStarterPoints(
        SegmentationResults[] expected, 
        List<List<PairIntArray>> edgesList) throws Exception {
        
        int nStarterPoints = 16;
        
        SFit[] sFits = new SFit[nStarterPoints];
        
        int count = 0;
        int[] i0s = new int[]{
            tLen.lowIdx + 
                (int)(0.25*(tLen.highIdx - tLen.lowIdx)),
            tLen.lowIdx + 
                (int)(0.75*(tLen.highIdx - tLen.lowIdx))
        };
        int[] i1s = new int[]{
            tColor.lowIdx + 
                (int)(0.25*(tColor.highIdx - tColor.lowIdx)),
            tColor.lowIdx + 
                (int)(0.75*(tColor.highIdx - tColor.lowIdx))
        };
        int[] i2s = new int[]{
            tR.lowIdx + 
                (int)(0.25*(tR.highIdx - tR.lowIdx)),
            tR.lowIdx + 
                (int)(0.75*(tR.highIdx - tR.lowIdx))
        };
        int[] i3s = new int[]{
            tSmallMerge.lowIdx + 
                (int)(0.25*(tSmallMerge.highIdx - tSmallMerge.lowIdx)),
            tSmallMerge.lowIdx + 
                (int)(0.75*(tSmallMerge.highIdx - tSmallMerge.lowIdx))
        };
        for (int i0 : i0s) {
            for (int i1 : i1s) {
                for (int i2 : i2s) {
                    for (int i3 : i3s) {
                        sFits[count] =
                            invoke(tLen.getValue(i0),
                            tColor.getValue(i1),
                            tR.getValue(i2),
                            tSmallMerge.getValue(i3),
                            expected, edgesList);
                        count++;
                    }
                }
            }
        }
        
        return sFits;
    }

    private SFit[] condense(SFit[] yFits) {
        SFit[] output = new SFit[yFits.length];
        int count = 0;
        for (SFit sFit : yFits) {
            if (sFit != null) {
                output[count] = sFit;
                count++;
            }
        }
        if (count < yFits.length) {
            output = Arrays.copyOf(output, count);
        }
        return output;
    }

    private boolean hasConverged(SFit prevBest, SFit sFit) {

        if (prevBest == null) {
            return false;
        }
        
        // quick check for exactly the same
        if (prevBest.tLenF == sFit.tLenF) {
            if (prevBest.tColorF == sFit.tColorF) {
                if (prevBest.tRF == sFit.tRF) {
                    if (prevBest.tSmallMergeF == sFit.tSmallMergeF) {
                        if (prevBest.costF == sFit.costF) {
                            return true;
                        }
                    }
                }
            }
        }
        
        return false;
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
      
    public class SFit {
        int tLenF;
        double tColorF;
        double tRF;
        double tSmallMergeF;
        int colorSpaceF;
        double costF;
        boolean reduceNoiseF;
        public SFit copy() {
            SFit cp = new SFit();
            cp.tLenF = tLenF;
            cp.tColorF = tColorF;
            cp.tRF = tRF;
            cp.tSmallMergeF = tSmallMergeF;
            cp.colorSpaceF = colorSpaceF;
            cp.reduceNoiseF = reduceNoiseF;
            cp.costF = costF;
            return cp;
        }
        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append("cost=").append(Double.toString(costF))
                .append(" tLen=").append(Integer.toString(tLenF))
                .append(" tColor=").append(Double.toString(tColorF))
                .append(" tSmallMerge=")
                .append(Double.toString(tSmallMergeF))
                .append(" tR=")
                .append(Double.toString(tRF))
                .append(" reduceNoiseF=")
                .append(Boolean.toString(reduceNoiseF))
                .append(" colorSpace=")
                .append(Double.toString(colorSpaceF))
            ;
            return sb.toString();
        }
    }
}
