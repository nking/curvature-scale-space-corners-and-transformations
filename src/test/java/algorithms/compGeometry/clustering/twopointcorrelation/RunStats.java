package algorithms.compGeometry.clustering.twopointcorrelation;

import java.lang.management.ManagementFactory;
import java.lang.management.OperatingSystemMXBean;


/**
 * class to hold basic stats for an algorithm's execution run times
 *
 * @author nichole
 */
public class RunStats {

    protected final double[] loadAverages;
    protected final double[] runtimes;
    protected final long[] freeMemory;
    protected final long totalMemory;
    protected final int[] numProcessors;
    protected int currentMemoryIndex = -1;
    protected int currentRuntimeIndex = -1;

    protected int npoints = -1;

    protected double averageRuntime = -1;
    protected double stdDevRuntime = -1;
    protected double averageFreeMemory = -1;
    protected double stdDevFreeMemory = -1;

    private double averageLoadAverage = -1;
    private double stdDevLoadAverage = -1;

    protected final String name;
    private String runttimeUnit = "seconds";

    public RunStats(int nIterations, String name) {
        this.runtimes = new double[nIterations];
        this.freeMemory = new long[nIterations];
        this.numProcessors = new int[nIterations];
        this.totalMemory = Runtime.getRuntime().totalMemory();
        this.name = name;
        this.loadAverages = new double[nIterations];
    }

    public void recordSystemStats() {
        currentMemoryIndex++;

        Runtime runtime = Runtime.getRuntime();
        long memFree = runtime.freeMemory();
        int nProcessors = runtime.availableProcessors();

        freeMemory[currentMemoryIndex] = memFree;
        numProcessors[currentMemoryIndex] = nProcessors;

        OperatingSystemMXBean osBean = ManagementFactory.getOperatingSystemMXBean();
        loadAverages[currentMemoryIndex] = osBean.getSystemLoadAverage();
    }

    public void addRuntime(double runtime) {
        currentRuntimeIndex++;
        runtimes[currentRuntimeIndex] = runtime;
    }

    public void calculateFirstMomentRuntimeStats() {

        double sum = 0;
        for (int i = 0; i < runtimes.length; i++) {
            sum += runtimes[i];
        }
        this.averageRuntime = sum / runtimes.length;

        sum = 0.;
        for (int i = 0; i < runtimes.length; i++) {
            double diff = runtimes[i] - averageRuntime;
            sum += Math.pow( diff, 2);
        }
        this.stdDevRuntime = Math.sqrt(sum)/runtimes.length;


        sum = 0;
        for (int i = 0; i < freeMemory.length; i++) {
            sum += freeMemory[i];
        }
        this.averageFreeMemory = sum / freeMemory.length;

        sum = 0.;
        for (int i = 0; i < freeMemory.length; i++) {
            double diff = freeMemory[i] - averageFreeMemory;
            sum += Math.pow( diff, 2);
        }
        this.stdDevFreeMemory = Math.sqrt(sum)/freeMemory.length;

        sum = 0;
        for (int i = 0; i < loadAverages.length; i++) {
            sum += loadAverages[i];
        }
        this.averageLoadAverage = sum / loadAverages.length;

        sum = 0.;
        for (int i = 0; i < loadAverages.length; i++) {
            double diff = loadAverages[i] - getAverageLoadAverage();
            sum += Math.pow( diff, 2);
        }
        this.stdDevLoadAverage = Math.sqrt(sum)/loadAverages.length;
    }

    public double getAverageRuntime() {
        return averageRuntime;
    }

    public double getStdDevRuntime() {
        return stdDevRuntime;
    }

    public double getAverageFreeMemory() {
        return averageFreeMemory;
    }

    public double getStdDevFreeMemory() {
        return stdDevFreeMemory;
    }

    public long getTotalMemory() {
        return totalMemory;
    }

    public int getNumberOfProcessors() {
        return numProcessors[0];
    }

    public String getSummary() {

        StringBuffer summary = new StringBuffer();

        double aveRuntime = getAverageRuntime();
        double aveFreeMemory = getAverageFreeMemory();      

        double loadAverage = getAverageLoadAverage();
        summary.append(name).append(": ").append("\n")
            .append(" num processors= ").append(numProcessors[0]).append("\n")
            .append(" ave load avg = ").append(loadAverage).append(" +- ").append(stdDevLoadAverage).append("\n")
            .append(" tot memory = ").append(totalMemory).append("\n")
            .append(" ave free memory = ").append(aveFreeMemory).append(" +- ").append(stdDevFreeMemory).append("\n")
            .append(" ave program run time = ").append(aveRuntime).append(" +- ").append(stdDevRuntime).append(" ").append(getRunttimeUnit()).append("\n")
            .append("\n");

        return summary.toString();
    }

    public String getName() {
        return name;
    }

    /**
     * @return the runttimeUnit
     */
    public String getRunttimeUnit() {
        return runttimeUnit;
    }

    /**
     * @param runttimeUnit the runttimeUnit to set
     */
    public void setRunttimeUnit(String runttimeUnit) {
        this.runttimeUnit = runttimeUnit;
    }

    /**
     * @return the averageLoadAverage
     */
    public double getAverageLoadAverage() {
        return averageLoadAverage;
    }

    /**
     * @return the stdDevLoadAverage
     */
    public double getStdDevLoadAverage() {
        return stdDevLoadAverage;
    }

    /**
     * @return the npoints
     */
    public int getNpoints() {
        return npoints;
    }

    /**
     * @param npoints the npoints to set
     */
    public void setNpoints(int npoints) {
        this.npoints = npoints;
    }
}
