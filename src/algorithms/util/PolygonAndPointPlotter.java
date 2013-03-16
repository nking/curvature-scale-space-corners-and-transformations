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

    public PolygonAndPointPlotter(float minX, float maxX, float minY, float maxY) throws FileNotFoundException, IOException {

        plotContent = getTemplateHtmlPlot();

        setDataMinMax(plotContent, minX, maxX, minY, maxY);
    }

    public PolygonAndPointPlotter() throws FileNotFoundException, IOException {

        plotContent = getTemplateHtmlPlot();
    }

    protected void setDataMinMax(StringBuffer plotContent, float minX, float maxX, float minY, float maxY) {

        StringBuffer dataSB = new StringBuffer();

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
        float[] xPoints, float[] yPoints, float[] xPolygon, float[] yPolygon, String plotLabel) {

        setDataMinMax(plotContent, minX, maxX, minY, maxY);

        addPlot(xPoints, yPoints, xPolygon, yPolygon, plotLabel);
    }

    public void addPlot(int[] xPoints, float[] yPoints, int[] xPolygon, float[] yPolygon, String plotLabel) {
    	float[] xx = new float[xPoints.length];
    	float[] xp = new float[xPoints.length];
    	for (int i = 0; i < xx.length; i++) {
    		xx[i] = (float)xPoints[i];
    		xp[i] = (float)xPolygon[i];
    	}
    	addPlot(xx, yPoints, xp, yPolygon, plotLabel);
    }

    public void addPlot(float[] xPoints, float[] yPoints, float[] xPolygon, float[] yPolygon, String plotLabel) {
        addPlot(xPoints, yPoints, null, null, xPolygon, yPolygon, plotLabel);
    }

    public void addPlot(float[] xPoints, float[] yPoints, float[] xErrPoints, float[] yErrPoints,
        float[] xPolygon, float[] yPolygon, String plotLabel) {

        StringBuffer dataSB = new StringBuffer();

        //  ===== add plotLabel data =====
        dataSB.append("\n\n").append("var data_plot_label = '").append(plotLabel).append("';\n");

        //  ===== add points data =====
        dataSB.append("\n\n").append("var data_points = [\n");
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


          //  ===== add polygon =====
        dataSB.append("\n").append("var data_polygon = [\n");
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


        // TODO: need the scatter plot for distances from center... for each group


        if (!dataMinMaxAreSet) {

            float minX = MiscMath.findMin(xPoints);
            float maxX = MiscMath.findMax(xPoints);
            float minY = MiscMath.findMin(yPoints);
            float maxY = MiscMath.findMax(yPoints);

            dataSB.append("\n").append("var xmin_").append(plotNumber).append("= 0.0").append(";\n");
            dataSB.append("\n").append("var xmax_").append(plotNumber).append("=")
                .append(maxX).append(";\n");
            dataSB.append("\n").append("var ymin_").append(plotNumber).append("= 0.0").append(";\n");
            dataSB.append("\n").append("var ymax_").append(plotNumber).append("=")
                .append(maxY).append(";\n");

            // ======= add RENDER statement ==========
            dataSB.append("renderPlot('plot").append(plotNumber).append("', data_points, data_polygon, data_plot_label,")
                .append(" xmin_").append(plotNumber).append(", ")
                .append(" xmax_").append(plotNumber).append(", ")
                .append(" ymin_").append(plotNumber).append(", ")
                .append(" ymax_").append(plotNumber)
                .append( ");\n");

        } else {

            // ======= add RENDER statement ==========
            dataSB.append("renderPlot('plot").append(plotNumber).append("', data_points, data_polygon, data_plot_label)\n");
        }

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

    public void addPlot(double[] xPoints, double[] yPoints, double[] xPolygon, double[] yPolygon, String plotLabel) {

        if (!dataMinMaxAreSet) {

            /*float minX = MiscMath.findMin(xPoints);
            float maxX = MiscMath.findMax(xPoints);
            float minY = MiscMath.findMin(yPoints);
            float maxY = MiscMath.findMax(yPoints);

            setDataMinMax(plotContent, minX, maxX, minY, maxY);*/
            throw new IllegalStateException("missing x and y min and maxes for plot");
        }

        StringBuffer dataSB = new StringBuffer();

        //  ===== add plotLabel data =====
        dataSB.append("\n\n").append("var data_plot_label = '").append(plotLabel).append("';\n");


        //  ===== add points data =====
        dataSB.append("\n\n").append("var data_points = [\n");
        for (int i = 0; i < xPoints.length; i++) {
            if (i > 0) {
                dataSB.append(",\n");
            }
            dataSB.append("    {x:").append(xPoints[i]).append(", y:").append(yPoints[i]).append("}");
        }
        dataSB.append("\n];\n");


          //  ===== add polygon =====
        dataSB.append("\n").append("var data_polygon = [\n");
        dataSB.append("    [");
        for (int ii = 0; ii < xPolygon.length; ii++) {
            String xStr = String.format("%.1f", xPolygon[ii]);
            String yStr = String.format("%.1f", yPolygon[ii]);
            if (ii > 0) {
                dataSB.append(", ");
            }
            dataSB.append("    {x:").append(xStr).append(", y:").append(yStr).append("}");
        }
        dataSB.append("],\n ");
        dataSB.append("];\n");


        // TODO: need the scatter plot for distances from center... for each group


        // ======= add RENDER statement ==========
        dataSB.append("renderPlot('plot").append(plotNumber).append("', data_points, data_polygon, data_plot_label)\n");


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

    protected StringBuffer getTemplateHtmlPlot() throws FileNotFoundException, IOException {
        return getTemplateHtmlPlot("plot_points_and_polygon.html");
    }

    protected StringBuffer getTemplateHtmlPlot(String fileName) throws FileNotFoundException, IOException {

        String path = ResourceFinder.findFileInResources(fileName);

        StringBuffer sb = new StringBuffer();

        FileReader reader = null;
        BufferedReader in = null;

        try {
            reader = new FileReader(new File(path));
            in = new BufferedReader(reader);

            String line = in.readLine();

            while (line != null) {
                sb.append(line).append("\n");
                line = in.readLine();
            }
        } finally {
            if (in != null) {
                in.close();
            }
            if (reader != null) {
                reader.close();
            }
        }

        return sb;
    }

    public void writeFile() throws IOException {
        writeToFile(this.plotContent.toString(), "points_and_polygon.html");
    }

    public void writeFile2() throws IOException {
        writeToFile(this.plotContent.toString(), "points_and_polygon2.html");
    }

    public void writeFile3() throws IOException {
        writeToFile(this.plotContent.toString(), "points_and_polygon3.html");
    }

    protected void writeToFile(String fileContent, String fileName) throws IOException {

        String copyFilePath = ResourceFinder.writeToCWD(fileContent, fileName);
    }

}
