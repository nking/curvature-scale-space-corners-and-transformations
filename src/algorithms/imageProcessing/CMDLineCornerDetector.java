package algorithms.imageProcessing;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

/**
 * a wrapper to allow a user to run parts of the application from a command
 * line interface (for corner detection).
 * 
 * @author nichole
 */
public class CMDLineCornerDetector {
    
    private static String sep = System.getProperty("file.separator");
    
    private static String eol = System.getProperty("line.separator");
    
    private static String cwd = System.getProperty("user.dir") + sep;
    
    public CMDLineCornerDetector() {
                        
    }
    
    protected String checkFilePath(String filePath) {
        
        filePath = filePath.trim();
        if (!filePath.startsWith(sep)) {
            filePath = cwd+ filePath;
        }
        
        File file = new File(filePath);
        
        if (!file.exists()) {
            throw new IllegalArgumentException(
                "ERROR could not find file path: " + filePath);
        }
        
        return filePath;
    }    

    protected void writeCornersFile(PairIntArray corners, 
        String inputImageName) throws IOException {
        
        String outFilePath = cwd + "corners_" + inputImageName + ".tsv";
        
        FileWriter fw = null;
        BufferedWriter writer = null;
        
        try {
            File file = new File(outFilePath);
            file.delete();
            file.createNewFile();

            fw = new FileWriter(file);
            writer = new BufferedWriter(fw);
            
            writer.write("#corners in image " + inputImageName);
            writer.write(eol);
            writer.write("#x\ty");
            writer.write(eol);
            
            for (int i = 0; i < corners.getN(); i++) {
                
                String line = String.format("%d\t%d", corners.getX(i), 
                    corners.getY(i));
                writer.write(line);
                writer.write(eol);
                
                if ((i % 10) == 0) {
                    writer.flush();
                }
            }

            writer.flush();

        } finally {

            if (writer != null) {
                writer.close();
            }
            if (fw != null) {
                fw.close();
            }
            
            System.out.println("wrote: " + outFilePath);
        }
    }
    
    public static void main(String[] args) {
        
        if (args == null || args.length == 0) {
            printUsage();
            return;
        }
        
        String filePath = null;
        boolean writeTextOutput = false;
        boolean writeImageOutput = false;
        boolean includeEdgesInImageOutput = false;
        boolean useLineDrawingMode = false;
        
        //TODO:  tidy log statements and logging.properties and add a debug 
        //       flag to print details to stdout
        
        for (int i = 0; i < args.length; i++) {
            String arg = args[i];
            
            if (arg.equalsIgnoreCase("-image1=")) {
                if (i < (args.length - 1)) {
                    i++;
                    filePath = args[i];
                }
            } else if (arg.startsWith("-image1=")) {
                int idx = args[i].indexOf("=");
                if (idx > -1) {
                    filePath =  args[i].substring(idx + 1);
                }
            } else if (arg.equalsIgnoreCase("-text_output")) {
                writeTextOutput = true;
            } else if (arg.equalsIgnoreCase("-mark_image")) {
                writeImageOutput = true;
            } else if (arg.equalsIgnoreCase("-include_edges")) {
                writeImageOutput = true;
                includeEdgesInImageOutput = true;
            } else if (arg.equalsIgnoreCase("-input_is_line_drawing")) {
                useLineDrawingMode = true;
            }
        }
        
        if (filePath == null) {
            System.err.println("-image1= is required");
            System.exit(1);
        }
        
        CMDLineCornerDetector cmdLineInvoker = new CMDLineCornerDetector();
              
        filePath = cmdLineInvoker.checkFilePath(filePath);
       
        
        if (!writeImageOutput) {
            writeTextOutput = true;
        }
        
        String imageFileName = 
            (new File(filePath)).toPath().getFileName().toString();
        imageFileName = imageFileName.substring(0, imageFileName.length() -4);
        
        GreyscaleImage img = null;
        try {
             img = ImageIOHelper.readImageAsGrayScaleG(filePath);
             
             CurvatureScaleSpaceCornerDetector detector = new CurvatureScaleSpaceCornerDetector(img);
               
             if (useLineDrawingMode) {
                 detector.useLineDrawingMode();
             }
             
             detector.findCorners();
        
             PairIntArray corners = 
                 detector.getCornersInOriginalReferenceFrame();
             
             if (writeTextOutput) {                 
                 cmdLineInvoker.writeCornersFile(corners, imageFileName);
             }
             
             if (writeImageOutput) {
                 
                 Image image = ImageIOHelper.readImageAsGrayScale(filePath);
                                  
                 if (includeEdgesInImageOutput) {
                     
                     List<PairIntArray> edges = 
                         detector.getEdgesInOriginalReferenceFrame();
                     
                     ImageIOHelper.addAlternatingColorCurvesToImage(edges, 
                         image);
                 }
                 
                 ImageIOHelper.addCurveToImage(corners, image, 2, 255, 0, 0);
                 
                 String outFilePath = cwd + "corners_" + imageFileName + ".png";
                 
                 ImageIOHelper.writeOutputImage(outFilePath, image);
                 
                 System.out.println(eol + "wrote: " + outFilePath);
             }

             System.out.println(eol);
             
        } catch (Exception e) {
            System.err.println("ERROR: " + e.getMessage());
            System.exit(1);
        }
    }
    
    public static void printUsage() {
                
        System.out.println(eol + "Runner to invoke corner detection within an image.");
        System.out.println(eol + "Required arguments:");
        System.out.println("  -image1=path_to_image");
        System.out.println(eol + "Optional arguments:");
        System.out.println("  -text_output");
        System.out.println("      writes a file called corners.tsv in the");
        System.out.println("      current directory.");
        System.out.println("  -mark_image");
        System.out.println("      writes a file called corners.png in the");
        System.out.println("      current directory. It's the image with");
        System.out.println("      found corners overplotted.");
        System.out.println("  -include_edges");
        System.out.println("      in the written corners.png, includes the");
        System.out.println("      edges in the image.");
        System.out.println("  -input_is_line_drawing");
        System.out.println("      changes internal logic to handle an image");
        System.out.println("      of lines or solid blocks.");
        System.out.println(eol
            + "Note that default is to create corners.tsv in the currect directory");
        System.out.println(eol);
    }
    
}
