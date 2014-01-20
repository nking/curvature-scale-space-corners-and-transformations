package algorithms.compGeometry.clustering.twopointcorrelation;

import algorithms.util.ResourceFinder;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;

/**
  code to produce plots for visualizing the results of and debugging TwoPointVoidStats
   
 * @author nichole
 */
public class TwoPointVoidStatsPlotter {

    protected final StringBuffer plotContent;

    protected int plotNumber = 0;

    public TwoPointVoidStatsPlotter() throws FileNotFoundException, IOException {
        plotContent = getTemplateHtmlPlot();
    }

    public void addTwoPointPlot(float[] xPoints, float[] yPoints,
        int[] point1IndexOfLine, int[] point2IndexOfLine,
        float xmin, float xmax, float ymin, float ymax) {

        StringBuffer dataSB = new StringBuffer();

        dataSB.append("\nvar xmin_").append(plotNumber).append(" = ").append(xmin).append(";\n");
        dataSB.append(  "var xmax_").append(plotNumber).append(" = ").append(xmax).append(";\n");
        dataSB.append(  "var ymin_").append(plotNumber).append(" = ").append(ymin).append(";\n");
        dataSB.append(  "var ymax_").append(plotNumber).append(" = ").append(ymax).append(";\n");

        String srchFor = "/* === DO NOT REMOVE THIS == START DATA */";
        int insertOffset = plotContent.indexOf(srchFor);
        if (insertOffset == -1) {
            throw new IllegalStateException("cannot find END DATA marker");
        }
        insertOffset += srchFor.length();
        plotContent.insert(insertOffset, dataSB.toString());


        dataSB = new StringBuffer();

        //  ===== add points data =====
        dataSB.append("\n\n").append("var data_points_").append(plotNumber).append(" = [\n");
        for (int i = 0; i < xPoints.length; i++) {
            if (i > 0) {
                dataSB.append(",\n");
            }
            dataSB.append("    {x:").append(xPoints[i]).append(", y:").append(yPoints[i]).append("}");
        }
        dataSB.append("\n];\n");


        // ===== add lines data =====
        dataSB.append("\n").append("var data_lines_").append(plotNumber).append(" = [\n");
        for (int i = 0; i < point1IndexOfLine.length; i++) {

            int index1 = point1IndexOfLine[i];
            int index2 = point2IndexOfLine[i];

            String xLeftStr  = String.format("%.1f", xPoints[index1]);
            String yLeftStr  = String.format("%.1f", yPoints[index1]);
            String xRightStr = String.format("%.1f", xPoints[index2]);
            String yRightStr = String.format("%.1f", yPoints[index2]);
            if (i > 0) {
                dataSB.append(",\n");
            }
            dataSB.append("    [{x:").append(xLeftStr).append(", y:").append(yLeftStr).append("}")
                .append(", {x:").append(xRightStr).append(", y:").append(yRightStr).append("}]");
        }
        dataSB.append("\n];\n");


        // ======= add RENDER statement ==========
        dataSB.append("renderPlot('plot").append(plotNumber).append("', ")
            .append("data_points_").append(plotNumber).append(", ")
            .append("data_lines_").append(plotNumber).append(", ")
            .append("xmin_").append(plotNumber).append(", ")
            .append("xmax_").append(plotNumber).append(", ")
            .append("ymin_").append(plotNumber).append(", ")
            .append("ymax_").append(plotNumber)
            .append(")\n");


        srchFor = "/* === DO NOT REMOVE THIS == END DATA */";
        insertOffset = plotContent.indexOf(srchFor);
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

    public void addHistogram(float[] values, float xmax, int nBins) {

        StringBuffer dataSB = new StringBuffer();

        //  ===== add values data =====
        dataSB.append("\n\n").append("var data_values_").append(plotNumber).append(" = [\n");
        for (int i = 0; i < values.length; i++) {
            if (i > 0) {
                dataSB.append(",\n");
            }
            dataSB.append(values[i]);
        }
        dataSB.append("\n];\n");

        dataSB.append("\n var xmax_").append(plotNumber).append(" = ").append(xmax).append(";");

        dataSB.append("\n var nbins_").append(plotNumber).append(" = ").append(nBins).append(";");

        // ======= add RENDER statement ==========
        dataSB.append("renderHistogram('plot").append(plotNumber).append("', ")
            .append("data_values_").append(plotNumber)
            .append(", xmax_").append(plotNumber).append(", nbins_").append(plotNumber)
            .append(")\n");


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
        return getTemplateHtmlPlot("plot_minima_stats.html");
    }

    protected StringBuffer getTemplateHtmlPlot(String fileName) throws FileNotFoundException, IOException {

        StringBuffer sb = new StringBuffer();

        Reader reader = null;
        BufferedReader in = null;

        try {
            String path = ResourceFinder.findFileInResources(fileName);

            reader = new FileReader(new File(path));
            in = new BufferedReader(reader);

            String line = in.readLine();

            while (line != null) {
                sb.append(line).append("\n");
                line = in.readLine();
            }

        } catch(IOException e) {

            ClassLoader cls = ResourceFinder.class.getClassLoader();

            InputStream input = cls.getResourceAsStream(fileName);

            if (input == null) {
                throw new IOException("could not find file " + fileName);
            }

            reader = new InputStreamReader(input);
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
        writeToFile(this.plotContent.toString(), "twoptvoid_stats.html");
    }

    public void writeFile2() throws IOException {
        writeToFile(this.plotContent.toString(), "twoptvoid_stats2.html");
    }

    protected void writeToFile(String fileContent, String fileName) throws IOException {

        String copyFilePath = ResourceFinder.writeToCWD(fileContent, fileName);
    }
}
