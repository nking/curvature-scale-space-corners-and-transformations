package algorithms.util;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.net.URL;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;

import algorithms.compGeometry.clustering.twopointcorrelation.AxisIndexer;
import algorithms.compGeometry.clustering.twopointcorrelation.RandomClusterAndBackgroundGenerator;
import algorithms.compGeometry.clustering.twopointcorrelation.TwoPointCorrelation;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class InputFileReaderTest extends TestCase {

    /**
     *
     */
    protected boolean runForCoverage = false;
    
    /**
     *
     * @throws Exception
     */
    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }

    /**
     *
     * @throws Exception
     */
    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }
    
    /**
     *
     * @throws Exception
     */
    public void testRead() throws Exception {
        
        int n = 1000;
        
        String filePath = writeFile(n);
        
        InputFileReader reader = new InputFileReader();
        
        reader.read(filePath);
        
        assertNotNull(reader.getX());
        
        assertTrue(reader.getX().length == n);
        
        assertNotNull(reader.getY());
        
        assertTrue(reader.getY().length == n);
        

        assertNotNull(reader.getXErrors());
        
        assertTrue(reader.getXErrors().length == n);
        
        assertNotNull(reader.getYErrors());
        
        assertTrue(reader.getYErrors().length == n);
       
    }
    
    /**
     *
     * @throws Exception
     */
    public void testMainRunner() throws Exception {
        
        if (!runForCoverage) {
            return;
        }
        
        RandomClusterAndBackgroundGenerator generator = new RandomClusterAndBackgroundGenerator();

        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        sr.setSeed( System.currentTimeMillis() );
        
        AxisIndexer indexer = null;
        
        indexer = generator.createIndexerWithRandomPoints(sr, 0, 300, 100, 300,
            3, 30, 60, 0.1f);
                
        String filePath = writeFile(indexer);
        
        IInputFileReader reader = new InputFileReader();
        
        MainRunner runner = new MainRunner();
        TwoPointCorrelation tpc = runner.run(reader, new String[]{"--file", filePath});        
        assertNotNull(tpc);
        assertTrue(tpc.getNumberOfGroups() > 0);
        
        /*
         *  System.out.println("     optional:  --twosigma");
            System.out.println("     optional:  --threesigma");
            System.out.println("     optional:  --background 0.3 (requires backgrounderror guesstimate at least)");
            System.out.println("     optional:  --backgrounderror 0.1");
         */
        
        runner = new MainRunner();
        tpc = runner.run(reader, new String[]{"--file", filePath, "--twosigma"});   
        assertNotNull(tpc);
        assertTrue(tpc.getNumberOfGroups() > 0);
        
        runner = new MainRunner();
        tpc = runner.run(reader, new String[]{"--file", filePath, "--threesigma"});   
        assertNotNull(tpc);
        assertTrue(tpc.getNumberOfGroups() > 0);
        
        runner = new MainRunner();
        tpc = runner.run(reader, new String[]{"--file", filePath, "--background", "10.0", "--backgrounderror", "0.1"});   
        assertNotNull(tpc);
        assertTrue(tpc.getNumberOfGroups() > 0);
        
        // test Main doesn't throw an exception
        algorithms.compGeometry.clustering.twopointcorrelation.Main.main(new String[]{"--file", filePath});
        
    }
    
    /**
     *
     * @param n
     * @return
     * @throws NoSuchAlgorithmException
     * @throws IOException
     */
    protected String writeFile(int n) throws NoSuchAlgorithmException, IOException {

        String filePath = getLocalPath() + "/tmp.txt";

        FileWriter writer = null;
        BufferedWriter out = null;

        SecureRandom srr = SecureRandom.getInstance("SHA1PRNG");
        srr.setSeed( System.currentTimeMillis() );
        long seed = srr.nextLong();

        try {
            File fl = new File(filePath);

            writer = new FileWriter(fl);
            out = new BufferedWriter(writer);
            
            for (int i = 0; i < n; i++) {

                StringBuilder sb = new StringBuilder();
                
                sb.append(i).append("\t").append(srr.nextInt(n)).append("\t").append(0.1).append("\t").append(0.1).append("\n");

                writer.write(sb.toString());
            }

        } finally {
            if (writer != null) {
                writer.close();
            }
            if (out != null) {
                out.close();
            }
        }
        
        return filePath;
    }
    
    /**
     *
     * @param indexer
     * @return
     * @throws NoSuchAlgorithmException
     * @throws IOException
     */
    protected String writeFile(AxisIndexer indexer) throws NoSuchAlgorithmException, IOException {
        
        String filePath = getLocalPath() + "/tmp.txt";
        
        FileWriter writer = null;
        BufferedWriter out = null;

        try {
            File fl = new File(filePath);

            writer = new FileWriter(fl);
            out = new BufferedWriter(writer);
            
            for (int i = 0; i < indexer.getNXY(); i++) {

                StringBuilder sb = new StringBuilder();
                
                sb.append(indexer.getX()[i]).append("\t").append(indexer.getY()[i]).append("\t")
                    .append(indexer.getXErrors()[i]).append("\t").append(indexer.getYErrors()[i]).append("\n");

                writer.write(sb.toString());
            }

        } finally {
            if (writer != null) {
                writer.close();
            }
            if (out != null) {
                out.close();
            }
        }
        
        return filePath;
    }
    
    /**
     *
     * @return
     */
    protected String getLocalPath() {
        
        ClassLoader cls = ResourceFinder.class.getClassLoader();
        URL url = cls.getResource(".");
       
        String filePath = url.getPath();
        
        return filePath;
    }
}
