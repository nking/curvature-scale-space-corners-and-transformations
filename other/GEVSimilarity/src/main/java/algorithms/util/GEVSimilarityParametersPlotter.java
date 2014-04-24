package algorithms.util;

import java.awt.Color;
import java.awt.Paint;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;
import java.util.logging.Logger;
import javax.imageio.ImageIO;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartRenderingInfo;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.NumberTickUnit;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.Range;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;


/**
 * a scatterplot writer to write 3 png files for k, sigma, and mu respectively.
 * 
 * <pre>
 * suggested usage:
 * 
 *     GEVSimilarityParametersPlotter plotter = new GEVSimilarityParametersPlotter(
 *         0, 100, 1E-3, 1E+3, 1E-1, 1E+1, -10, 10);
 *  
 *     plotter.addToPlot(yPointsPt1, kPointsPt1, sigmaPointsPt1, muPointsPt1);
 *     plotter.addToPlot(yPointsPt2, kPointsPt2, sigmaPointsPt2, muPointsPt2);
 * 
 *     plotter.writeToFile(0);
 * 
 *     plotter.clearPlotter();
 * 
 *     plotter.addToPlot(yPointsPt3, kPointsPt3, sigmaPointsPt3, muPointsPt3);
 * 
 *     plotter.writeToFile(1);
 * 
 * This will have written files:
 *     
 * 
 * </pre>
 * 
 * @author nichole
 */
public class GEVSimilarityParametersPlotter {
    
    protected final int yMin;
    protected final int yMax;
    protected final float kMin;
    protected final float kMax;
    protected final float sigmaMin;
    protected final float sigmaMax;
    protected final float muMin;
    protected final float muMax;
    
    protected int plotWidth = 4000;
    protected int plotHeight = 4000;
    
    protected float[] y;
    protected float[] k;
    protected float[] s;
    protected float[] m;
    protected int n = 0;
    
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    private enum PARAM {
        K, SIGMA, MU
    }

    public GEVSimilarityParametersPlotter(int yMinimum, int yMaximum, float kMin, float kMax, 
        float sigmaMin, float sigmaMax, float muMin, float muMax) 
        throws FileNotFoundException, IOException {

        this.yMin = yMinimum;
        this.yMax = yMaximum;
        this.kMin = kMin;
        this.kMax = kMax;
        this.sigmaMin = sigmaMin;
        this.sigmaMax = sigmaMax;
        this.muMin = muMin;
        this.muMax = muMax;
    }
    
    public void clearPlotter() {
        n = 0;
        y = null;
        k = null;
        s = null;
        m = null;
    }

    public void addToPlot(float[] yPoints, float[] kPoints, float[] sigmaPoints, float[] muPoints) {

        int addN = yPoints.length;
        
        if (n == 0) {
            y = new float[addN];
            k = new float[addN];
            s = new float[addN];
            m = new float[addN];
        } else {
            y = Arrays.copyOf(y, n + addN);
            k = Arrays.copyOf(k, n + addN);
            s = Arrays.copyOf(s, n + addN);
            m = Arrays.copyOf(m, n + addN);
        }
        
        for (int i = 0; i < addN; i++) {
            y[i + n] = yPoints[i];
            k[i + n] = kPoints[i];
            s[i + n] = sigmaPoints[i];
            m[i + n] = muPoints[i];
        }
        
        n += addN;
    }

    protected void writeToBinary(String suffix) throws IOException {
        
        BufferedImage img = writeToBinary(PARAM.K);
        String filePath = writeToPng(img, "gev_similarity_k" + suffix + ".png"); 
        log.info("wrote " + filePath);
       
        img = writeToBinary(PARAM.SIGMA);
        filePath = writeToPng(img, "gev_similarity_s" + suffix + ".png"); 
        log.info("wrote " + filePath);
        
        img = writeToBinary(PARAM.MU);
        filePath = writeToPng(img, "gev_similarity_m" + suffix + ".png"); 
        log.info("wrote " + filePath);
    }
    
    protected String writeToPng(BufferedImage img, String fileName) throws IOException {
        
        String dirPath = ResourceFinder.findDirectory("target");
        String filePath = dirPath + System.getProperty("file.separator") + fileName;
        
        File file = new File(filePath);
        
        ImageIO.write(img, "png", file);
        
        return filePath;
    }
   
    private BufferedImage writeToBinary(PARAM param) {
        
        XYSeriesCollection data = new XYSeriesCollection();
        
        boolean useLog = true;
        
        float[] x;
        float xMin, xMax;
        if (param == PARAM.K) {
            x = k;
            xMin = kMin;
            xMax = kMax;
        } else if (param == PARAM.SIGMA) {
            x = s;
            xMin = sigmaMin;
            xMax = sigmaMax;
        } else /*if (param == PARAM.MU)*/ {
            x = m;
            xMin = muMin;
            xMax = muMax;
            useLog = false;
        }
        
        XYSeries series = new XYSeries("");
        for (int i = 0; i < n; i++) {
            series.add(x[i], y[i]);
        }
        data.addSeries(series);
        
        JFreeChart chart = ChartFactory.createScatterPlot(
            null,
            param.toString(), null, 
            (XYDataset)data
        );
        
        Paint bPaint = new Color(255, 255, 255, 0);
        
        chart.setBackgroundPaint(bPaint);
        
        final XYPlot plot = chart.getXYPlot();
        
        Range xRange = new Range(xMin, xMax);
        Range yRange = new Range(yMin, yMax);
        
        NumberAxis yAxis = new NumberAxis();
        NumberAxis xAxis = null;
        if (useLog) {
            xAxis = new LogarithmicAxis(param.toString());
        } else {
            xAxis = new NumberAxis(param.toString());
        }
        xAxis.setRange(xRange, true, false);
        yAxis.setRange(yRange, true, false);
        yAxis.setTickUnit(new NumberTickUnit(1));
        
        
        plot.setDomainAxis(xAxis);
        plot.setRangeAxis(yAxis);
        
        plot.setBackgroundPaint(bPaint);
        plot.setDomainGridlinePaint(Color.BLACK);
        plot.setRangeGridlinePaint(Color.BLACK);
        
        ChartRenderingInfo chartInfo 
            = new ChartRenderingInfo();
        //    = new ChartRenderingInfo(new StandardEntityCollection());
        
        int imageType = BufferedImage.TYPE_INT_ARGB;
        
        // in the absence of an X11 environment, may need to set this:
        System.setProperty("java.awt.headless", "true");
               
        BufferedImage bImg = chart.createBufferedImage(500, 1400, imageType, chartInfo);
        
        return bImg;
    }
    
    public void writeFile() throws IOException {
        if (n == 0) {
            return;
        }
        
        writeToBinary("");
    }

    public void writeFile(int num) throws IOException {
        if (n == 0) {
            return;
        }
        
        writeToBinary("_" + Integer.toString(num));
    }
    
}
