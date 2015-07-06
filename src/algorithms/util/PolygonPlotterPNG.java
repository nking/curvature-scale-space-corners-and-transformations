package algorithms.util;

import algorithms.misc.MiscMath;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Paint;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import javax.imageio.ImageIO;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartRenderingInfo;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYPolygonAnnotation;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.NumberTickUnit;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.Range;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/**
 *
 * @author nichole
 */
public class PolygonPlotterPNG {

    protected boolean dataMinMaxAreSet = false;
    
    protected float minX, maxX, minY, maxY;

    protected int plotWidth = 4000;
    
    protected int plotHeight = 4000;
        
    protected JFreeChart chart = null;
    
    protected int lastColor = 0;
    
    public PolygonPlotterPNG(float minX, float maxX, float minY, 
        float maxY, String plotLabel, String xAxisLabel, String yAxisLabel) 
        throws FileNotFoundException, IOException {

        this.minX = minX;
        this.maxX = maxX;
        this.minY = minY;
        this.maxY = maxY;
                
        dataMinMaxAreSet = true;
        
        // ----- create the chart object -----
        chart = ChartFactory.createScatterPlot(
            null,
            plotLabel, null, 
            new XYSeriesCollection()
        );
        
        Paint bPaint = Color.WHITE;
        
        chart.setBackgroundPaint(bPaint);
        
        final XYPlot plot = chart.getXYPlot();
        
        Range xRange = new Range(minX, maxX);
        Range yRange = new Range(minY, maxY);
        
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
        
        // in the absence of an X11 environment, may need to set this:
        System.setProperty("java.awt.headless", "true");
    }

    public void addPolygon(int[] xPoints, int[] yPoints) {

        int n0 = (xPoints != null) ? xPoints.length : 0;
        
        double[] xy = new double[2*n0];
        
        for (int i = 0; i < n0; ++i) {
            xy[2*i] = xPoints[i];
            xy[2*i + 1] = yPoints[i];
        }
        
        Color clr = chooseNextColor();
        
        XYPolygonAnnotation a = new XYPolygonAnnotation(xy, null, null, clr);
        
        XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) 
            chart.getXYPlot().getRenderer();
        
        renderer.addAnnotation(a);
    }
    
    public void addPolygon(float[] xPoints, float[] yPoints) {

        int n0 = (xPoints != null) ? xPoints.length : 0;
        
        double[] xy = new double[2*n0];
        
        for (int i = 0; i < n0; ++i) {
            xy[2*i] = xPoints[i];
            xy[2*i + 1] = yPoints[i];
        }
        
        Color clr = chooseNextColor();
        
        XYPolygonAnnotation a = new XYPolygonAnnotation(xy, null, null, clr);
        
        XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) 
            chart.getXYPlot().getRenderer();
        
        renderer.addAnnotation(a);
    }
    
    public void addPolygon(double[] xPoints, double[] yPoints) {

        int n0 = (xPoints != null) ? xPoints.length : 0;
        
        double[] xy = new double[2*n0];
        
        for (int i = 0; i < n0; ++i) {
            xy[2*i] = xPoints[i];
            xy[2*i + 1] = yPoints[i];
        }
        
        Color clr = chooseNextColor();
        
        XYPolygonAnnotation a = new XYPolygonAnnotation(xy, 
            new BasicStroke(3f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_BEVEL), 
            clr);
        
        XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) 
            chart.getXYPlot().getRenderer();
                
        renderer.addAnnotation(a);
    }
    
    public void addPolygon(double[] xPoints, double[] yPoints, Color clr) {

        int n0 = (xPoints != null) ? xPoints.length : 0;
        
        double[] xy = new double[2*n0];
        
        for (int i = 0; i < n0; ++i) {
            xy[2*i] = xPoints[i];
            xy[2*i + 1] = yPoints[i];
        }
                
        XYPolygonAnnotation a = new XYPolygonAnnotation(xy, 
            new BasicStroke(4f), clr);
        
        XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) 
            chart.getXYPlot().getRenderer();
                
        renderer.addAnnotation(a);
    }
    
    private Color chooseNextColor() {
                
        if (lastColor > 5) {
            lastColor = 0;
        }
        Color c = Color.BLACK;
        switch(lastColor) {
            case 1:
                c = Color.PINK;
                break;
            case 2:
                c = Color.GREEN;
                break;
            case 3:
                c = Color.RED;
                break;
            case 4:
                c = Color.CYAN;
                break;
            case 5:
                c = Color.ORANGE;
                break;
            case 6:
                c = Color.DARK_GRAY;
                break;
            case 7:
                c = Color.BLUE;
                break;
            case 8:
                c = Color.LIGHT_GRAY;
                break;
            case 9:
                c = Color.MAGENTA;
                break;
            default:
                break;
        }            
        lastColor++;
        return c;
    }
    
    public String writeFile() throws IOException {
        return writeToFile("polygon.png");
    }

    public String writeFile(Integer num) throws IOException {
        return writeToFile( 
            "polygon" + num.toString() +".png");
    }

    protected String writeToFile(String fileName) throws IOException {

        ChartRenderingInfo chartInfo 
            = new ChartRenderingInfo();
        //    = new ChartRenderingInfo(new StandardEntityCollection());
        
        int imageType = BufferedImage.TYPE_INT_ARGB;
        
        BufferedImage imgToWrite = chart.createBufferedImage(500, 500, 
            imageType, chartInfo);
        
        String filePath = ResourceFinder.getAFilePathInCWD(fileName);
        
        File file = new File(filePath);
        
        ImageIO.write(imgToWrite, "png", file);
        
        return filePath;        
    }

}
