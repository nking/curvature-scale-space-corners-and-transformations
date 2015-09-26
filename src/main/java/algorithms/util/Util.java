package algorithms.util;

import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscMath;
import java.io.IOException;
import java.lang.management.MemoryMXBean;
import java.lang.management.ManagementFactory;

public class Util {
    
    public static String plotHistogram(HistogramHolder hist,
        String label, String outputFileSuffix) throws IOException {
                
        float[] xh = hist.getXHist();
        float[] yh = hist.getYHistFloat();
        
        float yMin = MiscMath.findMin(yh);
        int yMaxIdx = MiscMath.findYMaxIndex(yh);
        if (yMaxIdx == -1) {
            return null;
        }
        float yMax = yh[yMaxIdx];
        
        float xMin = MiscMath.findMin(xh);
        float xMax = MiscMath.findMax(xh);        
                
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();

        plotter.addPlot(
            xMin, xMax, yMin, yMax,
            xh, yh, xh, yh, label);

        return plotter.writeFile(outputFileSuffix);
    }
    
    public String plotLogHistogram(HistogramHolder hist,
        String label, String outputFileSuffix) throws IOException {
                
        float[] xh = hist.getXHist();
        float[] yh = hist.getYHistFloat();
        
        float[] yLogH = new float[yh.length];
        for (int i = 0; i < yh.length; ++i) {
            yLogH[i] = (float)Math.log(yh[i]/Math.log(10));
        }
        
        float yMin = MiscMath.findMin(yLogH);
        int yMaxIdx = MiscMath.findYMaxIndex(yLogH);
        float yMax = yLogH[yMaxIdx];
        
        float xMin = MiscMath.findMin(xh);
        float xMax = MiscMath.findMax(xh);
                        
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();

        plotter.addPlot(
            xMin, xMax, yMin, yMax,
            xh, yLogH, xh, yLogH, label);

        return plotter.writeFile(outputFileSuffix);
    }
     
    /**
    get the memory available to the heap as the total memory available to the JVM
    minus the amount of heap used.  
    Note that it does not know if the amount
    of memory available to the program from the OS is less than the JVM expects
    from command line flags
    */
    public static long getAvailableHeapMemory() {

        MemoryMXBean mbean = ManagementFactory.getMemoryMXBean();
        long heapUsage = mbean.getHeapMemoryUsage().getUsed();
        long totMemory = Runtime.getRuntime().totalMemory();

        long totAvail = totMemory - heapUsage;

        return totAvail;
    }

}
