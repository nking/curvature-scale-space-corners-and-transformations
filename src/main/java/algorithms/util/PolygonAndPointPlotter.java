package algorithms.util;

import algorithms.misc.MiscMath;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

/**
 *
 * @author nichole
 */
public class PolygonAndPointPlotter {

    protected final StringBuffer plotContent;

    protected int plotNumber = 0;

    protected boolean dataMinMaxAreSet = false;
    
    protected float minX, maxX, minY, maxY;

    public PolygonAndPointPlotter(float minX, float maxX, float minY, 
        float maxY) throws FileNotFoundException, IOException {

        plotContent = getTemplateHtmlPlot();

        setDataMinMax(plotContent, minX, maxX, minY, maxY);
    }

    public PolygonAndPointPlotter() throws FileNotFoundException, IOException {

        plotContent = getTemplateHtmlPlot();
    }

    protected void setDataMinMax(StringBuffer plotContent, 
        float minX, float maxX, float minY, float maxY) {

        this.minX = minX;
        this.maxX = maxX;
        this.minY = minY;
        this.maxY = maxY;
        
        StringBuilder dataSB = new StringBuilder();

        dataSB.append("\nvar xmin = ").append(minX).append(";\n");
        dataSB.append("var xmax = ").append(maxX).append(";\n");
        dataSB.append("var ymin = ").append(minY).append(";\n");
        dataSB.append("var ymax = ").append(maxY).append(";\n");

        String srchFor = "/* === DO NOT REMOVE THIS == START DATA */";
        int insertOffset = plotContent.indexOf(srchFor);
        if (insertOffset == -1) {
            throw new IllegalStateException("cannot find END DATA marker");
        }
        insertOffset += srchFor.length();
        plotContent.insert(insertOffset, dataSB.toString());

        dataMinMaxAreSet = true;
    }

    public void addPlot(float minX, float maxX, float minY, float maxY,
        int[] xPoints, int[] yPoints, int[] xPolygon, int[] yPolygon, 
        String plotLabel) {

        int n0 = (xPoints != null) ? xPoints.length : 0;
        int n1 = (xPolygon != null) ? xPolygon.length : 0;
        
        float[] x0 = new float[n0];
    	float[] y0 = new float[n0];
        float[] x1 = new float[n1];
    	float[] y1 = new float[n1];
        
        for (int i = 0; i < n0; i++) {
            x0[i] = xPoints[i];
            y0[i] = yPoints[i];
        }
        
        for (int i = 0; i < n1; i++) {
            x1[i] = xPolygon[i];
            y1[i] = yPolygon[i];
        }
        
        addPlot(minX, maxX, minY, maxY, x0, y0, null, null,
            x1, y1, plotLabel);
    }
    
    public void addPlot(float minX, float maxX, float minY, float maxY,
        float[] xPoints, float[] yPoints, float[] xPolygon, float[] yPolygon, 
        String plotLabel) {

        addPlot(minX, maxX, minY, maxY, xPoints, yPoints, null, null,
            xPolygon, yPolygon, plotLabel);
    }

    public void addPlot(int[] xPoints, float[] yPoints, int[] xPolygon, 
        float[] yPolygon, String plotLabel) {
        
        int n0 = (xPoints != null) ? xPoints.length : 0;
        int n1 = (xPolygon != null) ? xPolygon.length : 0;
      
    	float[] xx = new float[n0];
        for (int i = 0; i < n0; i++) {
    		xx[i] = (float)xPoints[i];
    	}
        
    	float[] xp = new float[n1];
    	for (int i = 0; i < n1; i++) {
    		xp[i] = (float)xPolygon[i];
    	}
        
    	addPlot(xx, yPoints, null, null, xp, yPolygon, plotLabel);
    }

    public void addPlot(float[] xPoints, float[] yPoints, float[] xPolygon, 
        float[] yPolygon, String plotLabel) {
        
        addPlot(xPoints, yPoints, null, null, xPolygon, yPolygon, plotLabel);
    }
    
    public void addPlot(int[] xPoints, int[] yPoints, int[] xPolygon, 
        int[] yPolygon, String plotLabel) {
        
        int n0 = (xPoints != null) ? xPoints.length : 0;
        int n1 = (xPolygon != null) ? xPolygon.length : 0;
        
        float[] x0 = new float[n0];
    	float[] y0 = new float[n0];
    	for (int i = 0; i < n0; i++) {
    		x0[i] = (float)xPoints[i];
    		y0[i] = (float)yPoints[i];
    	}
        float[] x1 = new float[n1];
    	float[] y1 = new float[n1];
    	for (int i = 0; i < n1; i++) {
    		x1[i] = (float)xPolygon[i];
    		y1[i] = (float)yPolygon[i];
    	}
        
        float[] Null = null;
        
        addPlot(x0, y0, Null, Null, x1, y1, plotLabel);
    }
    
    public void addPlot(float[] xPoints, float[] yPoints, int[] xPolygon, 
        int[] yPolygon, String plotLabel) {
        
        int n0 = (xPoints != null) ? xPoints.length : 0;
        int n1 = (xPolygon != null) ? xPolygon.length : 0;
        
        float[] x1 = new float[n1];
    	float[] y1 = new float[n1];
    	for (int i = 0; i < n1; i++) {
    		x1[i] = (float)xPolygon[i];
    		y1[i] = (float)yPolygon[i];
    	}
        
        float[] Null = null;
        
        addPlot(xPoints, yPoints, Null, Null, x1, y1, plotLabel);
    }
    
    public void addPlot(float[] xPoints, float[] yPoints, float[] xErrPoints, 
        float[] yErrPoints, float[] xPolygon, float[] yPolygon, 
        String plotLabel) {
        
        if (!dataMinMaxAreSet) {

            float minX0 = MiscMath.findMin(xPoints);
            float maxX0 = MiscMath.findMax(xPoints);
            float minY0 = MiscMath.findMin(yPoints);
            float maxY0 = MiscMath.findMax(yPoints);
            
            addPlot(minX0, maxX0, minY0, maxY0, xPoints, yPoints, 
                xErrPoints, yErrPoints, xPolygon, yPolygon, plotLabel);
            
        } else {
            
            addPlot(minX, maxX, minY, maxY, xPoints, yPoints, 
                xErrPoints, yErrPoints, xPolygon, yPolygon, plotLabel);
            
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
     * @param xPolygon can be null
     * @param yPolygon can be null
     * @param plotLabel 
     */
    public void addPlot(float xmn, float xmx, float ymn, float ymx,
        float[] xPoints, float[] yPoints, float[] xErrPoints, float[] yErrPoints,
        float[] xPolygon, float[] yPolygon, String plotLabel) {

        StringBuffer dataSB = new StringBuffer("\n");
        
        //  ===== add plotLabel data =====
        dataSB.append("var data_plot_label_").append(plotNumber)
            .append(" = '").append(plotLabel).append("';\n");

        //  ===== add points data =====
        if (xPoints == null) {
            dataSB.append("var data_points_").append(plotNumber)
                .append(" = undefined;\n");
        } else {
            dataSB.append("var data_points_").append(plotNumber).append(" = [\n");
            for (int i = 0; i < xPoints.length; i++) {
                if (i > 0) {
                    dataSB.append(",\n");
                }
                dataSB.append("    {x:").append(xPoints[i]).append(", y:").append(yPoints[i]);
                if (xErrPoints != null) {
                    dataSB.append(", dx:").append(xErrPoints[i]).append(", dy:").append(yErrPoints[i]);
                }
                dataSB.append("}");
            }
            dataSB.append("\n];\n");
        }

        if (xPolygon == null) {
            dataSB.append("var data_polygon_").append(plotNumber)
                .append(" = undefined;\n");
        } else {
            //  ===== add polygon =====
            dataSB.append("var data_polygon_").append(plotNumber).append(" = [\n");
            dataSB.append("    [");
            for (int ii = 0; ii < xPolygon.length; ii++) {
                String xStr = String.format("%.7f", xPolygon[ii]);
                String yStr = String.format("%.7f", yPolygon[ii]);
                if (ii > 0) {
                    dataSB.append(", ");
                }
                dataSB.append("    {x:").append(xStr).append(", y:").append(yStr).append("}");
            }
            dataSB.append("],\n ");
            dataSB.append("];\n");
        }

        dataSB.append("var xmin_").append(plotNumber).append("=")
            .append(xmn).append(";\n");
        dataSB.append("var xmax_").append(plotNumber).append("=")
            .append(xmx).append(";\n");
        dataSB.append("var ymin_").append(plotNumber).append("=")
            .append(ymn).append(";\n");
        dataSB.append("var ymax_").append(plotNumber).append("=")
            .append(ymx).append(";\n");

        // ======= add RENDER statement ==========
        dataSB.append("\nrenderPlot('plot").append(plotNumber)
            .append("', data_points_").append(plotNumber)
            .append(", data_polygon_").append(plotNumber)
            .append(", data_plot_label_").append(plotNumber)
            .append(", ")
            .append(" xmin_").append(plotNumber).append(", ")
            .append(" xmax_").append(plotNumber).append(", ")
            .append(" ymin_").append(plotNumber).append(", ")
            .append(" ymax_").append(plotNumber)
            .append( ");\n\n");

        String srchFor = "/* === DO NOT REMOVE THIS == END DATA */";
        int insertOffset = plotContent.indexOf(srchFor);
        if (insertOffset == -1) {
            throw new IllegalStateException("cannot find END DATA marker");
        }
        plotContent.insert(insertOffset, dataSB.toString());
        dataSB = null;


        // ========== add the PLOT DIVS ==============
        StringBuffer plotDivs = new StringBuffer();
        plotDivs.append("<div id='plot").append(plotNumber).append("' class='plot'></div>\n");


        srchFor = "<!-- === DO NOT REMOVE THIS == END PLOT DIVS -->";
        insertOffset = plotContent.indexOf(srchFor);
        if (insertOffset == -1) {
            throw new IllegalStateException("cannot find END DATA marker");
        }
        plotContent.insert(insertOffset, plotDivs.toString());
        plotDivs = null;

        plotNumber++;
    }

    public void addPlot(double[] xPoints, double[] yPoints, double[] xPolygon, 
        double[] yPolygon, String plotLabel) {

        int n0 = (xPoints != null) ? xPoints.length : 0;
        int n1 = (xPolygon != null) ? xPolygon.length : 0;
        
        float[] xx = new float[n0];
        float[] yy = new float[n0];
    	for (int i = 0; i < n0; i++) {
    		xx[i] = (float)xPoints[i];
            yy[i] = (float)yPoints[i];
    	}
        
        float[] xp = new float[n1];
        float[] yp = new float[n1];
        for (int i = 0; i < n1; i++) {
            xp[i] = (float) xPolygon[i];
            yp[i] = (float) yPolygon[i];
        }
        
        addPlot(xx, yy, null, null, xp, yp, plotLabel);
    }

    protected StringBuffer getTemplateHtmlPlot() throws FileNotFoundException, IOException {
        return getTemplateHtmlPlot("plot_points_and_polygon.html");
    }

    protected StringBuffer getTemplateHtmlPlot(String fileName) throws FileNotFoundException, IOException {

        String path = ResourceFinder.findFileInResources(fileName);

        StringBuffer sb = new StringBuffer();

        BufferedReader in = null;

        try {
            in = new BufferedReader(new FileReader(new File(path)));

            String line = in.readLine();

            while (line != null) {
                sb.append(line).append("\n");
                line = in.readLine();
            }
        } finally {
            if (in != null) {
                in.close();
            }
        }

        return sb;
    }

    public String writeFile() throws IOException {
        return writeToFile(this.plotContent.toString(), "points_and_polygon.html");
    }

    public String writeFile2() throws IOException {
        return writeToFile(this.plotContent.toString(), "points_and_polygon2.html");
    }

    public String writeFile3() throws IOException {
        return writeToFile(this.plotContent.toString(), "points_and_polygon3.html");
    }
    
    public String writeFile(Integer num) throws IOException {
        return writeToFile(this.plotContent.toString(), 
            "points_and_polygon" + num.toString() +".html");
    }
    
    public String writeFile(String fileSuffix) throws IOException {
        return writeToFile(this.plotContent.toString(), 
            "points_and_polygon_" + fileSuffix + ".html");
    }
    
    public String writeFile(long num) throws IOException {
        return writeToFile(this.plotContent.toString(), 
            "points_and_polygon" + Long.toString(num) +".html");
    }

    protected String writeToFile(String fileContent, String fileName) throws IOException {

        return ResourceFinder.writeToCWD(fileContent, fileName);
    }

}
