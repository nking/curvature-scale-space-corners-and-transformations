package algorithms.compGeometry.clustering.twopointcorrelation;

import algorithms.misc.MiscMath;
import algorithms.util.ArrayPair;
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
<pre>
Convenience class to create a plot to visualize results of the
Two-Point correlation code.

It edits a template made with d3.js and creates an output file that
can be loaded into the browser and viewed without a network connection
(the d3.js script is embedded).
The path to the output file is printed to standard out.

   Usage:
      TwoPointCorrelationPlotter plotter = new TwoPointCorrelationPlotter(xMin, xMax, yMin, yMax);
      plotter.addPlot(instanceOfTwoPointCorrelation);
      plotter.writeFile();

   Note that multiple uses of addPlot will append to existing html file for this instance.
</pre>
@author nichole
*/
public class TwoPointCorrelationPlotter {

    protected final StringBuffer plotContent;

    protected int plotNumber = 0;

    public TwoPointCorrelationPlotter(float minX, float maxX, float minY, float maxY) throws FileNotFoundException, IOException {

        plotContent = getTemplateHtmlPlot();

        setDataMinMax(plotContent, minX, maxX, minY, maxY);
    }

    public TwoPointCorrelationPlotter(float minX, float maxX, float minY, float maxY, String relDirectory) throws FileNotFoundException, IOException {

        plotContent = getTemplateHtmlPlotWRTDir(relDirectory);

        setDataMinMax(plotContent, minX, maxX, minY, maxY);
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
    }

    public void addPlot(TwoPointCorrelation twoPtCorr) {
        addPlot(twoPtCorr, null);
    }

    public void addPlot(TwoPointCorrelation twoPtCorr, String plotLabel2) {

        StringBuffer dataSB = new StringBuffer();

        float[] xPoints = twoPtCorr.getX();
        float[] yPoints = twoPtCorr.getY();

        //  ===== add points data =====
        dataSB.append("\n\n").append("var data_points = [\n");
        for (int i = 0; i < xPoints.length; i++) {
            if (i > 0) {
                dataSB.append(",\n");
            }
            dataSB.append("    {x:").append(xPoints[i]).append(", y:").append(yPoints[i]);
            if (twoPtCorr.getXErrors() != null) {
                dataSB.append(", dx:").append(twoPtCorr.getXErrors()[i]).append(", dy:").append(twoPtCorr.getYErrors()[i]);
            }
            dataSB.append("}");
        }
        dataSB.append("\n];\n");


          //  ===== add group centroids =====
        dataSB.append("\n\n").append("var data_group_centroids = [\n");
        ArrayPair centroids = twoPtCorr.getHullCentroids();
        for (int i = 0; i < twoPtCorr.getNumberOfGroups(); i++) {
            if (i > 0) {
                dataSB.append(",\n");
            }
            dataSB.append("    {x:").append(centroids.getX()[i])
                .append(", y:").append(centroids.getY()[i]).append(", name: ").append(i).append("}");
        }
        dataSB.append("\n];\n");


        dataSB.append("\n").append("var data_polygons = [\n");
        for (int i = 0; i < twoPtCorr.getNumberOfGroups(); i++) {
            ArrayPair hull = twoPtCorr.getGroupHull(i);
            dataSB.append("    [");
            for (int ii = 0; ii < hull.getX().length; ii++) {
                String xStr = String.format("%.1f", hull.getX()[ii]);
                String yStr = String.format("%.1f", hull.getY()[ii]);
                if (ii > 0) {
                    dataSB.append(", ");
                }
                dataSB.append("    {x:").append(xStr).append(", y:").append(yStr).append("}");
            }
            if (i == (twoPtCorr.getNumberOfGroups() - 1)) {
                dataSB.append("]\n");
            } else {
                dataSB.append("],\n ");
            }
        }
        dataSB.append("];\n");

        String sdStr = String.format("'%.6f'", twoPtCorr.getBackgroundDensity());
        dataSB.append("\n").append("var plot_label_").append(plotNumber).append("=").append(sdStr).append("\n;");

        // ======= add RENDER statement ==========
        dataSB.append("renderPlot('plot").append(plotNumber)
            .append("', data_group_centroids, data_points, data_polygons, plot_label_").append(plotNumber)
            .append(");\n");

        if ((twoPtCorr.backgroundStats != null) && (twoPtCorr.backgroundStats instanceof TwoPointVoidStats) ) {

            TwoPointVoidStats stats = (TwoPointVoidStats)twoPtCorr.backgroundStats;

            float[] xh = stats.statsHistogram.getXHist();
            float[] yh = stats.statsHistogram.getYHistFloat();
            float[] yg = stats.bestFit.getOriginalScaleYFit();

            float[] xhe = stats.statsHistogram.getXErrors();
            float[] yhe = stats.statsHistogram.getYErrors();

            float maxx = xh[xh.length - 1] + ((xh[1] - xh[0])/2.f);
            float minx = xh[0] - ((xh[1] - xh[0])/2.f);
            float miny = 0.f;
            float maxy = MiscMath.findMax(yh);

            if (plotLabel2 == null) {

                plotLabel2 = String.format("chisq=%.1f errsq=%.1f chistat=%.1e chisq/err=%.1e  chistat/err=%.1e",
                stats.bestFit.getChiSqSum(), stats.bestFit.getYDataErrSq(),
                stats.bestFit.getChiSqStatistic(),
                (stats.bestFit.getChiSqSum()/stats.bestFit.getYDataErrSq()),
                (stats.bestFit.getChiSqStatistic()/stats.bestFit.getYDataErrSq())
                );
            }

            dataSB.append("\n\n").append("var data_plot_gev_label = '").append(plotLabel2).append("';\n");

            //  ===== add points data =====
            dataSB.append("\n\n").append("var data_gev_points = [\n");
            for (int i = 0; i < xh.length; i++) {
                if (i > 0) {
                    dataSB.append(",\n");
                }
                dataSB.append("    {x:").append(xh[i]).append(", y:").append(yh[i]);
                if (xhe != null) {
                    String dxStr = String.format("%.7f", xhe[i]);
                    String dyStr = String.format("%.7f", yhe[i]);
                    dataSB.append(", dx:").append(dxStr).append(", dy:").append(dyStr);
                }
                dataSB.append("}");
            }
            dataSB.append("\n];\n");

              //  ===== add polygon =====
            dataSB.append("\n").append("var data_gev_polygon = [\n");
            dataSB.append("    [");
            for (int ii = 0; ii < xh.length; ii++) {
                String xStr = String.format("%.7f", xh[ii]);
                String yStr = String.format("%.7f", yg[ii]);
                if (ii > 0) {
                    dataSB.append(", ");
                }
                dataSB.append("    {x:").append(xStr).append(", y:").append(yStr);
                dataSB.append("}");
            }
            dataSB.append("],\n ");
            dataSB.append("];\n");

            // ======= add RENDER statement ==========
            dataSB.append("renderGEVPlot('plot_gev").append(plotNumber)
                .append("', data_gev_points, data_gev_polygon, data_plot_gev_label, ")
                .append(minx).append(",").append(maxx).append(",").append(miny).append(",").append(maxy).append(");\n");

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
        if (twoPtCorr.backgroundStats != null) {
            plotDivs.append("<div id='plot_gev").append(plotNumber).append("' class='plot'></div>\n");
        }


        srchFor = "<!-- === DO NOT REMOVE THIS == END PLOT DIVS -->";
        insertOffset = plotContent.indexOf(srchFor);
        if (insertOffset == -1) {
            throw new IllegalStateException("cannot find END DATA marker");
        }
        plotContent.insert(insertOffset, plotDivs.toString());
        plotDivs = null;

        plotNumber++;
    }

    public void addPlotWithoutHull(TwoPointCorrelation twoPtCorr, String plotLabel2) {

        StringBuffer dataSB = new StringBuffer();

        float[] xPoints = twoPtCorr.getX();
        float[] yPoints = twoPtCorr.getY();

        //  ===== add points data =====
        dataSB.append("\n\n").append("var data_points = [\n");
        for (int i = 0; i < xPoints.length; i++) {
            if (i > 0) {
                dataSB.append(",\n");
            }
            dataSB.append("    {x:").append(xPoints[i]).append(", y:").append(yPoints[i]);
            if (twoPtCorr.getXErrors() != null) {
                dataSB.append(", dx:").append(twoPtCorr.getXErrors()[i]).append(", dy:").append(twoPtCorr.getYErrors()[i]);
            }
            dataSB.append("}");
        }
        dataSB.append("\n];\n");


          //  ===== add group centroids =====
        dataSB.append("\n\n").append("var data_group_centroids = [\n");
        ArrayPair centroids = twoPtCorr.getHullCentroids();
        for (int i = 0; i < twoPtCorr.getNumberOfGroups(); i++) {
            if (i > 0) {
                dataSB.append(",\n");
            }
            dataSB.append("    {x:").append(centroids.getX()[i])
                .append(", y:").append(centroids.getY()[i]).append(", name: ").append(i).append("}");
        }
        dataSB.append("\n];\n");


        dataSB.append("\n").append("var data_groups = [\n");
        for (int i = 0; i < twoPtCorr.getNumberOfGroups(); i++) {
            ArrayPair group = twoPtCorr.getGroup(i);
        
            int nPolygonPoints = group.getX().length;

            dataSB.append("    [");
            for (int ii = 0; ii < nPolygonPoints; ii++) {
                String xStr = String.format("%.1f", group.getX()[ii]);
                String yStr = String.format("%.1f", group.getY()[ii]);
                if (ii > 0) {
                    dataSB.append(",\n");
                }
                dataSB.append("    {x:").append(xStr).append(", y:").append(yStr).append("}");
            }
            if (i == (twoPtCorr.getNumberOfGroups() - 1)) {
                dataSB.append("]\n");
            } else {
                dataSB.append("],\n ");
            }
        }
        dataSB.append("];\n");

        String sdStr = String.format("'%.6f'", twoPtCorr.getBackgroundDensity());
        dataSB.append("\n").append("var plot_label_").append(plotNumber).append("=").append(sdStr).append(";\n");

        // ======= add RENDER statement ==========
        dataSB.append("\nrenderPlotWithoutHull('plot").append(plotNumber)
            .append("', data_group_centroids, data_points, data_groups, plot_label_").append(plotNumber)
            .append(");\n");

        if ((twoPtCorr.backgroundStats != null) && (twoPtCorr.backgroundStats instanceof TwoPointVoidStats) ) {

            TwoPointVoidStats stats = (TwoPointVoidStats)twoPtCorr.backgroundStats;

            float[] xh = stats.statsHistogram.getXHist();
            float[] yh = stats.statsHistogram.getYHistFloat();
            float[] yg = stats.bestFit.getOriginalScaleYFit();

            float[] xhe = stats.statsHistogram.getXErrors();
            float[] yhe = stats.statsHistogram.getYErrors();

            float maxx = xh[xh.length - 1] + ((xh[1] - xh[0])/2.f);
            float minx = xh[0] - ((xh[1] - xh[0])/2.f);
            float miny = 0.f;
            float maxy = MiscMath.findMax(yh);

            if (plotLabel2 == null) {

                plotLabel2 = String.format("chisq=%.1f errsq=%.1f chistat=%.1e chisq/err=%.1e  chistat/err=%.1e",
                stats.bestFit.getChiSqSum(), stats.bestFit.getYDataErrSq(),
                stats.bestFit.getChiSqStatistic(),
                (stats.bestFit.getChiSqSum()/stats.bestFit.getYDataErrSq()),
                (stats.bestFit.getChiSqStatistic()/stats.bestFit.getYDataErrSq())
                );
            }

            dataSB.append("\n\n").append("var data_plot_gev_label = '").append(plotLabel2).append("';\n");

            //  ===== add points data =====
            dataSB.append("\n\n").append("var data_gev_points = [\n");
            for (int i = 0; i < xh.length; i++) {
                if (i > 0) {
                    dataSB.append(",\n");
                }
                dataSB.append("    {x:").append(xh[i]).append(", y:").append(yh[i]);
                if (xhe != null) {
                    String dxStr = String.format("%.7f", xhe[i]);
                    String dyStr = String.format("%.7f", yhe[i]);
                    dataSB.append(", dx:").append(dxStr).append(", dy:").append(dyStr);
                }
                dataSB.append("}");
            }
            dataSB.append("\n];\n");

              //  ===== add polygon =====
            dataSB.append("\n").append("var data_gev_polygon = [\n");
            dataSB.append("    [");
            for (int ii = 0; ii < xh.length; ii++) {
                String xStr = String.format("%.7f", xh[ii]);
                String yStr = String.format("%.7f", yg[ii]);
                if (ii > 0) {
                    dataSB.append(", ");
                }
                dataSB.append("    {x:").append(xStr).append(", y:").append(yStr);
                dataSB.append("}");
            }
            dataSB.append("],\n ");
            dataSB.append("];\n");

            // ======= add RENDER statement ==========
            dataSB.append("renderGEVPlot('plot_gev").append(plotNumber)
                .append("', data_gev_points, data_gev_polygon, data_plot_gev_label, ")
                .append(minx).append(",").append(maxx).append(",").append(miny).append(",").append(maxy).append(");\n");

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
        if (twoPtCorr.backgroundStats != null) {
            plotDivs.append("<div id='plot_gev").append(plotNumber).append("' class='plot'></div>\n");
        }


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
        return getTemplateHtmlPlot("plot_twoptcorrelation.html");
    }

    protected StringBuffer getTemplateHtmlPlotWRTDir(String relDir) throws FileNotFoundException, IOException {
        return getTemplateHtmlPlot(relDir, "plot_twoptcorrelation.html");
    }

    protected StringBuffer getTemplateHtmlPlot(String subDir, String fileName) throws FileNotFoundException, IOException {

        String path = ResourceFinder.findFileInResources(subDir, fileName);

        StringBuffer sb = new StringBuffer();

        Reader reader = null;
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

    public String writeFile() throws IOException {
        return writeToFile(this.plotContent.toString(), "twoptcorrelation.html");
    }
    public String writeFile3() throws IOException {
        return writeToFile(this.plotContent.toString(), "twoptcorrelation3.html");
    }

    protected String writeToFile(String fileContent, String fileName) throws IOException {

        String copyFilePath = ResourceFinder.writeToCWD(fileContent, fileName);
        return copyFilePath;
    }
}
