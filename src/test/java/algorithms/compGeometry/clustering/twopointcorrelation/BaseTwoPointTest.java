package algorithms.compGeometry.clustering.twopointcorrelation;

import algorithms.compGeometry.clustering.twopointcorrelation.RandomClusterAndBackgroundGenerator.CLUSTER_SEPARATION;
import algorithms.util.ResourceFinder;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.logging.Logger;

import junit.framework.TestCase;

/**
 * utility class for unit tests which need to create random background points and
 * points in clusters.  It's meant to be extended.
 * TODO: change to use composition instead of inheritance?
 *
 * @author nichole
 */
public class BaseTwoPointTest extends TestCase {

    protected RandomClusterAndBackgroundGenerator generator = null;

    protected static Logger log0 = Logger.getLogger(BaseTwoPointTest.class.getName());
    
    @Override
    protected void setUp() throws Exception {

        super.setUp();

        generator = new RandomClusterAndBackgroundGenerator();
    }

    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }

    protected int getExpectedNumberOfClusters() {
        return generator.getExpectedNumberOfClusters();
    }
    protected int getTotalNumberOfPoints() {
        return generator.getTotalNumberOfPoints();
    }

    public void testNeededForJunitRuntime() {}

    protected DoubleAxisIndexer createIndexerWithRandomPoints() throws NoSuchAlgorithmException {
        return generator.createIndexerWithRandomPoints();
    }

    protected DoubleAxisIndexer createIndexerWithRandomPoints(float xmin,
        float xmax, float ymin, float ymax) throws NoSuchAlgorithmException {

        return generator.createIndexerWithRandomPoints(xmin, xmax, ymin, ymax);
    }

    protected DoubleAxisIndexer createIndexerWithRandomPoints(SecureRandom sr, float xmin, float xmax, float ymin, float ymax,
        int[] nClusters, int nBackgroundPoints, CLUSTER_SEPARATION clusterSeparation) {

        return generator.createIndexerWithRandomPoints(sr, xmin, xmax, ymin, ymax,
            nClusters, nBackgroundPoints, clusterSeparation);
    }

    protected DoubleAxisIndexer createIndexerWithRandomPoints(SecureRandom sr, float xmin, float xmax, float ymin, float ymax,
        int numberOfClusters, int minimumNumberOfPointsPerCluster, int maximumNumberOfPointsPerCluster,
        float backgroundPointFractionToClusters) {

        return generator.createIndexerWithRandomPoints(sr, xmin, xmax, ymin, ymax,
            numberOfClusters, minimumNumberOfPointsPerCluster, maximumNumberOfPointsPerCluster,
            backgroundPointFractionToClusters);
    }
    
    protected DoubleAxisIndexer createIndexerWithRandomPointsAroundCenterWithDSquared(
        SecureRandom sr, int numberOfClusterPoints,
        float xmin, float xmax, float ymin, float ymax, float maximumRadius) {
       
        return generator.createIndexerWithRandomPointsAroundCenterWithDSquared(
            sr, numberOfClusterPoints,
            xmin, xmax, ymin, ymax, maximumRadius);
    }

    protected void createRandomPointsAroundCenter(SecureRandom sr, float maxRadius,
        int numberOfPoints, float xc, float yc, float[] x, float[] y, int xyStartOffset) {

        generator.createRandomPointsAroundCenter(sr, maxRadius, numberOfPoints,
            xc, yc, x, y, xyStartOffset);
    }

    protected void createRandomSeparatedPoints(SecureRandom sr, float xmin, float xmax, float ymin, float ymax,
        int numberOfPointsToCreate, float[] xPoints, float[] yPoints, float minSeparationBetweenPoints) {

        generator.createRandomSeparatedPoints(sr, xmin, xmax, ymin, ymax,
            numberOfPointsToCreate, xPoints, yPoints, minSeparationBetweenPoints);
    }

    protected boolean separationBetweenExistingPointsIsLargerThanMin(float[] x, float[] y, int nXY,
        float xp, float yp, float minimumSeparation) {

        return generator.separationBetweenExistingPointsIsLargerThanMin(x, y, nXY,
            xp, yp, minimumSeparation);
    }

    protected void createRandomPointsInRectangle(SecureRandom sr, int nBackgroundPoints,
        float xmin, float xmax, float ymin, float ymax,
        float[] x, float[] y,  int xyStartOffset) {

        generator.createRandomPointsInRectangle(sr, nBackgroundPoints,
            xmin, xmax, ymin, ymax, x, y, xyStartOffset);
    }

    protected void createPoints(int numberOfBackgroundPoints, int[] numberOfClusterPoints,
        CLUSTER_SEPARATION clusterSeparation,
        float xmin, float xmax, float ymin, float ymax, SecureRandom sr, boolean useRandomForErrors) {

        generator.createPoints(numberOfBackgroundPoints, numberOfClusterPoints,
            clusterSeparation, xmin, xmax, ymin, ymax, sr, useRandomForErrors);
    }

     /**
     *      |
     *      |
     * -----|.....  <---- angle is w.r.t y=0, x=xc.  increases in CW order
     *      |
     *      |
     *
     * @param xc
     * @param yc
     * @param radius
     * @param angleInDegreesFromYEQ0XGT0  angle in degrees, CW from point y=0, x=xc
     * @return
     */
    public static float[] calculateXAndYFromXcYcAndRadius(float xc, float yc, float radius, double angleInDegreesFromYEQ0XGT0) {

       return RandomClusterAndBackgroundGenerator.calculateXAndYFromXcYcAndRadius(xc, yc, radius, angleInDegreesFromYEQ0XGT0);
    }

    /**
     *   *  |  *
     *      |
     * -----|.....  <---- angle is w.r.t y=0, x=xc.  increases in CCW order
     *  *   |
     *      |  *
     *
     * @param xc
     * @param yc
     * @param radius
     * @param angleInDegreesFromYEQ0XGT0  angle in degrees, CCW from point y=0, x=xc
     * @return
     */
    static float[] calculateXAndYFromXcYcAndRadiusCCW(float xc, float yc, float radius, double angleInDegreesFromYEQ0XGT0) {

        return RandomClusterAndBackgroundGenerator.calculateXAndYFromXcYcAndRadiusCCW(xc, yc, radius, angleInDegreesFromYEQ0XGT0);
    }
    
    static void printHistogramToStandardOut(float[] xp, int[] yp, float[] xpe, float[] ype) {
        StringBuilder xsb = new StringBuilder();
        StringBuilder ysb = new StringBuilder();
        StringBuilder xesb = new StringBuilder();
        StringBuilder yesb = new StringBuilder();
        for (int z = 0; z < xp.length; z++) {
            if (z > 0) {
                xsb.append("f, ");
                ysb.append("f, ");
                xesb.append("f, ");
                yesb.append("f, ");
            }
            xsb.append(xp[z]);
            ysb.append(yp[z]);
            xesb.append(xpe[z]);
            yesb.append(ype[z]);
        }
        log0.fine("float[] x = new float[]{"  + xsb.append("f").toString() + "};");
        log0.fine("float[] y = new float[]{"  + ysb.append("f").toString() + "};");
        log0.fine("float[] xe = new float[]{" + xesb.append("f").toString() + "};");
        log0.fine("float[] ye = new float[]{" + yesb.append("f").toString() + "};");
    }

    static void writeIndexerToTmpData(DoubleAxisIndexer indexer, int count) throws IOException {
        // write to tmpdata if need to use in tests improve fits, histogram etc
        String str = String.valueOf(count);
        while (str.length() < 3) {
            str = "0" + str;
        }
        String fileNamePostfix = "_clusters_" + str + ".dat";
        String fileName = CreateClusterDataTest.indexerFileNamePrefix + fileNamePostfix;
        String filePath = ResourceFinder.getAFilePathInTmpData(fileName);
        CreateClusterDataTest.writeIndexer(filePath, indexer);
    }

    static void writeVoidDensitiesToTestResources(String fileName, float[] values, float[] valueErrors) 
        throws Exception {

        String filePath = ResourceFinder.getAFilePathInTestResources(fileName);

        FileWriter writer = null;
        BufferedWriter out = null;

        try {
            writer = new FileWriter(new File(filePath));
            out = new BufferedWriter(writer);
            
            out.write(Integer.toString(values.length));
            out.write("\n");

            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < values.length; i++) {
                
                sb.append(values[i]).append("\t").append(valueErrors[i]).append("\n");
            }
            
            out.write(sb.toString());
            
            out.flush();
            
        } finally {
            if (writer != null) {
                writer.close();
            }
            if (out != null) {
                out.close();
            }
        }
    }
}
