package algorithms.imageProcessing;

/**
 * a wrapper to suggest command line usage options.
 * 
 * @author nichole
 */
public class CMDLineMain {
        
    private static String eol = System.getProperty("line.separator");
        
    public CMDLineMain() {           
    }
    
    public static void main(String[] args) {
        
        if (args == null || args.length == 0) {
            printUsage();
            return;
        }
    }
    
    public static void printUsage() {
                
        System.out.print(eol + "There are two classes to invoke code from the ");
        System.out.println("command line.");
        System.out.println("One for the corner deteector and the other for the image transformer via");
        System.out.println("inflection points.");
        System.out.println(eol + "Corner Detector:");
        System.out.println("     java -classpath 'insert path to scalespace.jar' algorithms.imageProcessing.CMDLineCornerDetector");
        System.out.println(eol + "Image Transformer/Inflection Points Detector:");
        System.out.println("     java -classpath 'insert path to scalespace.jar' algorithms.imageProcessing.CMDLineInflectionUtil");
        System.out.print(eol);
    }
    
}
