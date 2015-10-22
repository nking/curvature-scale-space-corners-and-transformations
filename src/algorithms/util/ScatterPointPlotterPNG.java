package algorithms.util;

import algorithms.misc.MiscMath;
import java.awt.Color;
import java.awt.Paint;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import javax.imageio.ImageIO;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartRenderingInfo;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.labels.ItemLabelAnchor;
import org.jfree.chart.labels.ItemLabelPosition;
import org.jfree.chart.labels.StandardXYToolTipGenerator;
import org.jfree.chart.labels.XYItemLabelGenerator;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.chart.title.TextTitle;
import org.jfree.chart.title.Title;
import org.jfree.data.Range;
import org.jfree.data.xy.AbstractXYDataset;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.TextAnchor;

/**
 *
 * @author nichole
 */
public class ScatterPointPlotterPNG {

    protected boolean dataMinMaxAreSet = false;
    
    protected float minX, maxX, minY, maxY;

    protected int plotWidth = 4000;
    
    protected int plotHeight = 4000;
    
    BufferedImage imgToWrite = null;
    
    public ScatterPointPlotterPNG(float minX, float maxX, float minY, 
        float maxY) {

        this.minX = minX;
        this.maxX = maxX;
        this.minY = minY;
        this.maxY = maxY;
    }

    public ScatterPointPlotterPNG() {

    }

    public void plot(float minX, float maxX, float minY, float maxY,
        int[] xPoints, int[] yPoints, String plotLabel,
        String xAxisLabel, String yAxisLabel) {

        int n0 = (xPoints != null) ? xPoints.length : 0;
        
        float[] x0 = new float[n0];
    	float[] y0 = new float[n0];
        
        for (int i = 0; i < n0; i++) {
            x0[i] = xPoints[i];
            y0[i] = yPoints[i];
        }
        
        plot(minX, maxX, minY, maxY, x0, y0, null, null, plotLabel,
            xAxisLabel, yAxisLabel);
    }
    
    public void plot(float minX, float maxX, float minY, float maxY,
        PairIntArray xyPoints, String plotLabel,
        String xAxisLabel, String yAxisLabel) {

        int n0 = (xyPoints != null) ? xyPoints.getN() : 0;
        
        float[] x0 = new float[n0];
    	float[] y0 = new float[n0];
        
        for (int i = 0; i < n0; i++) {
            x0[i] = xyPoints.getX(i);
            y0[i] = xyPoints.getY(i);
        }
        
        plot(minX, maxX, minY, maxY, x0, y0, null, null, plotLabel,
            xAxisLabel, yAxisLabel);
    }
    
    public void plot(float minX, float maxX, float minY, float maxY,
        float[] xPoints, float[] yPoints, String plotLabel,
        String xAxisLabel, String yAxisLabel) {

        plot(minX, maxX, minY, maxY, xPoints, yPoints, null, null, plotLabel,
            xAxisLabel, yAxisLabel);
    }

    public void plot(int[] xPoints, float[] yPoints, String plotLabel,
        String xAxisLabel, String yAxisLabel) {
        
        int n0 = (xPoints != null) ? xPoints.length : 0;
      
    	float[] xx = new float[n0];
        for (int i = 0; i < n0; i++) {
    		xx[i] = (float)xPoints[i];
    	}
        
    	plot(xx, yPoints, null, null, plotLabel, xAxisLabel, yAxisLabel);
    }

    public void plot(float[] xPoints, float[] yPoints, String plotLabel,
        String xAxisLabel, String yAxisLabel) {
        
        plot(xPoints, yPoints, null, null, plotLabel, xAxisLabel, yAxisLabel);
    }
    
    public void plot(int[] xPoints, int[] yPoints, String plotLabel,
        String xAxisLabel, String yAxisLabel) {
        
        int n0 = (xPoints != null) ? xPoints.length : 0;
        
        float[] x0 = new float[n0];
    	float[] y0 = new float[n0];
    	for (int i = 0; i < n0; i++) {
    		x0[i] = (float)xPoints[i];
    		y0[i] = (float)yPoints[i];
    	}
        
        plot(x0, y0, null, null, plotLabel, xAxisLabel, yAxisLabel);
    }
    
    public void plot(float[] xPoints, float[] yPoints, float[] xErrPoints, 
        float[] yErrPoints, String plotLabel, String xAxisLabel, String yAxisLabel) {
        
        if (!dataMinMaxAreSet) {

            float minX0 = MiscMath.findMin(xPoints);
            float maxX0 = MiscMath.findMax(xPoints);
            float minY0 = MiscMath.findMin(yPoints);
            float maxY0 = MiscMath.findMax(yPoints);
            
            plot(minX0, maxX0, minY0, maxY0, xPoints, yPoints, 
                xErrPoints, yErrPoints, plotLabel, xAxisLabel, yAxisLabel);
            
        } else {
            
            plot(minX, maxX, minY, maxY, xPoints, yPoints, 
                xErrPoints, yErrPoints, plotLabel, xAxisLabel, yAxisLabel);
        }
    }

    /**
     * 
     * @param xmn
     * @param xmx
     * @param ymn
     * @param ymx
     * @param xPoints
     * @param yPoints
     * @param xErrPoints can be null
     * @param yErrPoints can be null
     * @param plotLabel 
     */
    public void plot(float xmn, float xmx, float ymn, float ymx,
        float[] xPoints, float[] yPoints, float[] xErrPoints, float[] yErrPoints,
        String plotLabel, String xAxisLabel, String yAxisLabel) {

        XYSeriesCollection data = new XYSeriesCollection();
        
        //  ===== add points data =====
        if (xPoints != null) {
            XYSeries series1 = new XYSeries("");
            for (int i = 0; i < xPoints.length; i++) {
                series1.add(xPoints[i], yPoints[i]);
            }
            data.addSeries(series1);
        }

        JFreeChart chart = ChartFactory.createScatterPlot(
            null,
            plotLabel, null, 
            (XYDataset)data
        );
        
        Paint bPaint = Color.WHITE;
        
        chart.setBackgroundPaint(bPaint);
        
        final XYPlot plot = chart.getXYPlot();
        
        Range xRange = new Range(xmn, xmx);
        Range yRange = new Range(ymn, ymx);
        
        NumberAxis yAxis = new NumberAxis(yAxisLabel);
        NumberAxis xAxis = new NumberAxis(xAxisLabel);
        xAxis.setRange(xRange, true, false);
        yAxis.setRange(yRange, true, false);
        //yAxis.setTickUnit(new NumberTickUnit(1));
        
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
        
        imgToWrite = chart.createBufferedImage(500, 500, imageType, chartInfo);
    }
    
    /**
     * 
     * @param xmn
     * @param xmx
     * @param ymn
     * @param ymx
     * @param xPoints
     * @param yPoints
     * @param plotLabel 
     */
    public void plotLabeledPoints(float xmn, float xmx, float ymn, float ymx,
        float[] xPoints, float[] yPoints, String plotLabel, String xAxisLabel, 
        String yAxisLabel) {

        if (xPoints == null) {
            return;
        }
        
        //  ===== add points data =====
        LabeledXYDataset data = new LabeledXYDataset();
        for (int i = 0; i < xPoints.length; i++) {
            data.add(xPoints[i], yPoints[i], Integer.toString(i));
        }

        JFreeChart chart = ChartFactory.createScatterPlot(
            null,
            plotLabel, null, 
            (XYDataset)data
        );
        chart.setTitle(plotLabel);
        
        Paint bPaint = Color.WHITE;
        
        chart.setBackgroundPaint(bPaint);
        
        final XYPlot plot = chart.getXYPlot();
        
        //http://stackoverflow.com/questions/14458402/jfreecart-line-chart-with-text-each-point
        XYItemRenderer renderer = plot.getRenderer();
        renderer.setBaseItemLabelGenerator(new LabelGenerator());
        renderer.setBaseItemLabelPaint(Color.green.darker());
        renderer.setBasePositiveItemLabelPosition(
            new ItemLabelPosition(ItemLabelAnchor.CENTER, TextAnchor.TOP_CENTER));
        renderer.setBaseItemLabelFont(
            renderer.getBaseItemLabelFont().deriveFont(10f));
        renderer.setBaseItemLabelsVisible(true);
        renderer.setBaseToolTipGenerator(new StandardXYToolTipGenerator());
                
        Range xRange = new Range(xmn, xmx);
        Range yRange = new Range(ymn, ymx);
        
        NumberAxis yAxis = new NumberAxis(yAxisLabel);
        NumberAxis xAxis = new NumberAxis(xAxisLabel);
        xAxis.setRange(xRange, true, false);
        yAxis.setRange(yRange, true, false);
        //yAxis.setTickUnit(new NumberTickUnit(1));
        
        plot.setDomainAxis(xAxis);
        plot.setRangeAxis(yAxis);
        
        plot.setBackgroundPaint(bPaint);
        plot.setDomainGridlinePaint(Color.BLACK);
        plot.setRangeGridlinePaint(Color.BLACK);
        
        ChartRenderingInfo chartInfo = new ChartRenderingInfo();
        
        int imageType = BufferedImage.TYPE_INT_ARGB;
        
        // in the absence of an X11 environment, may need to set this:
        System.setProperty("java.awt.headless", "true");
       
        imgToWrite = chart.createBufferedImage(500, 500, imageType, chartInfo);
    }

    public void plot(double[] xPoints, double[] yPoints, String plotLabel) {

        int n0 = (xPoints != null) ? xPoints.length : 0;
        
        float[] xx = new float[n0];
        float[] yy = new float[n0];
    	for (int i = 0; i < n0; i++) {
    		xx[i] = (float)xPoints[i];
            yy[i] = (float)yPoints[i];
    	}
        
        plot(xx, yy, null, null, plotLabel);
    }
    
    private static class LabelGenerator implements XYItemLabelGenerator {
        @Override
        public String generateLabel(XYDataset dataset, int series, int item) {
            LabeledXYDataset labelSource = (LabeledXYDataset) dataset;
            return labelSource.getLabel(series, item);
        }
    }
    
    private static class LabeledXYDataset extends AbstractXYDataset {

        protected static final long serialVersionUID = 1111122222;
        
        private List<Number> x = new ArrayList<Number>();
        private List<Number> y = new ArrayList<Number>();
        private List<String> label = new ArrayList<String>();

        public void add(double x, double y, String label){
            this.x.add(x);
            this.y.add(y);
            this.label.add(label);
        }

        public String getLabel(int series, int item) {
            return label.get(item);
        }

        @Override
        public int getSeriesCount() {
            return 1;
        }

        @SuppressWarnings({"rawtypes"})
        @Override
        public Comparable getSeriesKey(int series) {
            return "labeled";
        }

        @Override
        public int getItemCount(int series) {
            return label.size();
        }

        @Override
        public Number getX(int series, int item) {
            return x.get(item);
        }

        @Override
        public Number getY(int series, int item) {
            return y.get(item);
        }
    }

    public String writeFile() throws IOException {
        return writeToFile("points.png");
    }

    public String writeFile(Integer num) throws IOException {
        return writeToFile( 
            "points" + num.toString() +".png");
    }

    public String writeToFile(String fileName) throws IOException {

        if (imgToWrite == null) {
            throw new IllegalStateException("must add plot data first");
        }
        
        String filePath = ResourceFinder.getAFilePathInCWD(fileName);
        
        File file = new File(filePath);
        
        ImageIO.write(imgToWrite, "png", file);
        
        return filePath;        
    }

}
