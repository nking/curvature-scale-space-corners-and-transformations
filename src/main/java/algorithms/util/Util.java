package algorithms.util;

import java.lang.management.MemoryMXBean;
import java.lang.management.ManagementFactory;

public class Util {
    
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
