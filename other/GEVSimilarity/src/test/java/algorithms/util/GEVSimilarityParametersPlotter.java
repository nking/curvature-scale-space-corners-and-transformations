package algorithms.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;

/**
 * a scatterplot writer.  the current output is html.
 * 
 * NOTE:  this will be changed to produce png output instead of html soon
 * because the svg is too large in some of the html files for rendering
 * in a thin client.
 * 
 * @author nichole
 */
public class GEVSimilarityParametersPlotter {

    protected final StringBuffer plotContent;

    protected int plotNumber = 0;
    
    protected final int yMin;
    protected final int yMax;

    public GEVSimilarityParametersPlotter(int yMinimum, int yMaximum, float kMin, float kMax, 
        float sigmaMin, float sigmaMax, float muMin, float muMax) 
        throws FileNotFoundException, IOException {

        this.yMin = yMinimum;
        this.yMax = yMaximum;
        
        plotContent = getTemplateHtmlPlot();
        
        
        StringBuffer dataSB = new StringBuffer();
        dataSB.append("\n\ncreateSVG('plot").append(plotNumber).append("', ")
            .append(yMin).append(", ").append(yMax).append(", ")
            .append(Float.toString(kMin)).append(", ")
            .append(Float.toString(kMax)).append(", ")
            .append(Float.toString(sigmaMin)).append(", ")
            .append(Float.toString(sigmaMax)).append(", ")
            .append(Float.toString(muMin)).append(", ")
            .append(Float.toString(muMax))
            .append( ");\n");
        
        String srchFor = "/* === DO NOT REMOVE THIS == END DATA */";
        int insertOffset = plotContent.indexOf(srchFor);
        if (insertOffset == -1) {
            throw new IllegalStateException("cannot find END DATA marker");
        }
        plotContent.insert(insertOffset, dataSB.toString());
        dataSB = null;
        
        
        // ========== add the PLOT DIVS ==============
        StringBuffer plotDivs = new StringBuffer();
        plotDivs.append("<div id='plot").append(plotNumber).append("' class='plot'>\n");
        plotDivs.append("  <div id='plot").append(plotNumber).append("_k' class='plotInner'></div>\n");
        plotDivs.append("  <div id='plot").append(plotNumber).append("_s' class='plotInner'></div>\n");
        plotDivs.append("  <div id='plot").append(plotNumber).append("_m' class='plotInner'></div>\n");
        plotDivs.append("</div>\n\n");
        
        srchFor = "<!-- === DO NOT REMOVE THIS == END PLOT DIVS -->";
        insertOffset = plotContent.indexOf(srchFor);
        if (insertOffset == -1) {
            throw new IllegalStateException("cannot find END DATA marker");
        }
        plotContent.insert(insertOffset, plotDivs.toString());
        plotDivs = null;

    }

    public void addPlot(float[] y, float[] kPoints, float[] sigmaPoints, float[] muPoints) {

        StringBuffer dataSB = new StringBuffer();

        String kDataStr = "k_points_" + Integer.toString(plotNumber);
        String sDataStr = "s_points_" + Integer.toString(plotNumber);
        String mDataStr = "m_points_" + Integer.toString(plotNumber);

        //  ===== add points data =====
        dataSB.append("\n\n").append("var ").append(kDataStr).append(" = [\n");
        for (int i = 0; i < kPoints.length; i++) {
            if (i > 0) {
                dataSB.append(", ");
                if (i % 10 == 0) {
                    dataSB.append("\n");
                }
            }
            dataSB.append("{x:").append(kPoints[i]).append(", y:").append(y[i]).append("}");            
        }
        dataSB.append("\n];\n");


        dataSB.append("\n\n").append("var ").append(sDataStr).append(" = [\n");
        for (int i = 0; i < sigmaPoints.length; i++) {
            if (i > 0) {
                dataSB.append(", ");
                if (i % 10 == 0) {
                    dataSB.append("\n");
                }
            }
            dataSB.append("{x:").append(sigmaPoints[i]).append(", y:").append(y[i]).append("}");
        }
        dataSB.append("\n];\n");


        dataSB.append("\n").append("var ").append(mDataStr).append(" = [\n");
        for (int i = 0; i < muPoints.length; i++) {
            if (i > 0) {
                dataSB.append(", ");
                if (i % 10 == 0) {
                    dataSB.append("\n");
                }
            }
            dataSB.append("{x:").append(muPoints[i]).append(", y:").append(y[i]).append("}");
        }
        dataSB.append("\n];\n");

        // ======= add RENDER statements ==========
        dataSB.append("\n\nrenderK(").append(kDataStr).append(");\n");
        
        dataSB.append("\n\nrenderS(").append(sDataStr).append(");\n");
        
        dataSB.append("\n\nrenderM(").append(mDataStr).append(");\n");
        
        

        String srchFor = "/* === DO NOT REMOVE THIS == END DATA */";
        int insertOffset = plotContent.indexOf(srchFor);
        if (insertOffset == -1) {
            throw new IllegalStateException("cannot find END DATA marker");
        }
        plotContent.insert(insertOffset, dataSB.toString());
        dataSB = null;

        plotNumber++;
    }

    protected StringBuffer getTemplateHtmlPlot() throws FileNotFoundException, IOException {
        return getTemplateHtmlPlot("plot_gev_similarity_parameters.html");
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
        return writeToFile(this.plotContent.toString(), "gev_similarity_parameters.html");
    }

    public String writeFile(int num) throws IOException {
        return writeToFile(this.plotContent.toString(), "gev_similarity_parameters_" + Integer.toString(num) + ".html");
    }

    protected String writeToFile(String fileContent, String fileName) throws IOException {

        String copyFilePath = ResourceFinder.writeToCWD(fileContent, fileName);
        
        return copyFilePath;
    }

}
