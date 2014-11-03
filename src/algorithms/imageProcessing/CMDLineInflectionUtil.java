package algorithms.imageProcessing;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

/**
 * a wrapper to allow a user to run parts of the application from a command
 * line interface (for inflection detection and image transformation calculation).
 * 
 * @author nichole
 */
public class CMDLineInflectionUtil {
    
    private static String sep = System.getProperty("file.separator");
    
    private static String eol = System.getProperty("line.separator");
    
    private static String cwd = System.getProperty("user.dir") + sep;
    
    public CMDLineInflectionUtil() {
    }
    
    protected String checkFilePath(String filePath) {
        
        filePath = filePath.trim();
        if (!filePath.startsWith(sep)) {
            filePath = cwd + filePath;
        }
        
        File file = new File(filePath);
        
        if (!file.exists()) {
            throw new IllegalArgumentException(
                "ERROR could not find file path: " + filePath);
        }
        
        return filePath;
    }
    
    
    private void writeMatchingPoints(PairIntArray xy1, PairIntArray xy1Tr, 
        PairIntArray xy2, TransformationParameters transformationParams, 
        String imageFileName1, String imageFileName2) throws IOException {
        
        String outFilePath = cwd + "matching_points_" + imageFileName1 +
            "_" + imageFileName2 + ".tsv";
        
        FileWriter fw = null;
        BufferedWriter writer = null;
        
        try {
            File file = new File(outFilePath);
            file.delete();
            file.createNewFile();

            fw = new FileWriter(file);
            writer = new BufferedWriter(fw);
            
            writer.write("#image1\t\timage1 transformed\t\timage2");
            writer.write(eol);
            
            for (int i = 0; i < xy1.getN(); i++) {
                
                String line = String.format("%d\t%d\t\t%d\t%d\t\t%d\t%d",
                    xy1.getX(i), xy1.getY(i), 
                    xy1Tr.getX(i), xy1Tr.getY(i), xy2.getX(i), xy2.getY(i));
                
                writer.write(line);
                
                writer.write(eol);
                
                if ((i % 10) == 0) {
                    writer.flush();;
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
            
            System.out.println(eol + "wrote: " + outFilePath);
        }
    }
    
    private void writeEdges(List<PairIntArray> edges1, PairIntArray[] edges1Tr, 
        String imageFileName1, String imageFileName2) throws IOException {
        
        if (edges1.size() != edges1Tr.length) {
            throw new IllegalStateException(
                "each curve in edges1 must have a parallel curve "
                + "with the same number of points in edges1Tr");
        }
        
        String outFilePath = cwd + "transformed_edges_" + imageFileName1 +
            "_" + imageFileName2 + ".tsv";
        
        FileWriter fw = null;
        BufferedWriter writer = null;
        
        try {
            File file = new File(outFilePath);
            file.delete();
            file.createNewFile();

            fw = new FileWriter(file);
            writer = new BufferedWriter(fw);
            
            writer.write("#image1 points in edges\t\ttransformed points to image 2 reference frame");
            writer.write(eol);
            
            for (int i = 0; i < edges1.size(); i++) {
                                    
                PairIntArray edge1 = edges1.get(i);

                PairIntArray edge1Tr = edges1Tr[i];

                if (edge1.getN() != edge1Tr.getN()) {
                throw new IllegalStateException(
                    "each curve in edges1 must have a parallel curve " 
                     + "with the same number of points in edges1Tr");
                }

                for (int j = 0; j < edge1.getN(); j++) {
                    int x1 = edge1.getX(j);
                    int y1 = edge1.getY(j);
                    int x1Tr = edge1Tr.getX(j);
                    int y1Tr = edge1Tr.getY(j);

                    String line = String.format("%d\t%d\t\t%d\t%d",
                        x1, y1, x1Tr, y1Tr);

                    writer.write(line);
                    writer.write(eol);

                    if ((j % 10) == 0) {
                        writer.flush();;
                    }
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
            
            System.out.println(eol + "wrote: " + outFilePath);
        }
    }
    
    private void writeTransformation(TransformationParameters 
        transformationParams, String imageFileName1, String imageFileName2) 
        throws IOException {
        
        String outFilePath = cwd + "transformation_" + imageFileName1 + "_" + 
            imageFileName2 + ".txt";
        
        FileWriter fw = null;
        BufferedWriter writer = null;
        
        try {
            File file = new File(outFilePath);
            file.delete();
            file.createNewFile();

            fw = new FileWriter(file);
            writer = new BufferedWriter(fw);
            
            writer.write("#transformation parameters for " 
                + imageFileName1 + " to " + imageFileName2);
            writer.write(eol);
            
            double scale = transformationParams.getScale();
            double rotRad = transformationParams.getRotationInRadians();
            double rotDeg = transformationParams.getRotationInDegrees();
            double transX = transformationParams.getTranslationX();
            double transY = transformationParams.getTranslationY();
            
            writer.write("#scale\trotRad\trotDeg\ttransX\ttransY");
            
            writer.write(eol);
            
            writer.write(Double.toString(scale));
            writer.write("\t");
            writer.write(Double.toString(rotRad));
            writer.write("\t");
            writer.write(Double.toString(rotDeg));
            writer.write("\t");
            writer.write(Double.toString(transX));
            writer.write("\t");
            writer.write(Double.toString(transY));
            writer.write("\t");
            
            writer.flush();

        } finally {

            if (writer != null) {
                writer.close();
            }
            if (fw != null) {
                fw.close();
            }
            
            System.out.println(eol + "wrote: " + outFilePath);
        }
    }
    
    public static void main(String[] args) {
        System.out.println(Arrays.toString(args));
        if (args == null || args.length == 0) {
            printUsage();
            return;
        }
        
        String filePath1 = null;
        String filePath2 = null;
        boolean writeTextOutput = false;
        boolean writeImageOutput = false;
        boolean includeEdgesInImageOutput = false;
        boolean useLineDrawingMode = false;
        
        //TODO:  tidy log statements and logging.properties and add a debug 
        //       flag to print details to stdout
        
        //TODO:  add a flag to allow a faster less accurate solution that might
        //       be chosen when the camera is fixed so that scale and rotation
        //       don't change.  then use mapper.doNotRefineTransformations()
        
        for (int i = 0; i < args.length; i++) {
            
            String arg = args[i];
            
            if (arg.equalsIgnoreCase("-image1=")) {
                if (i < (args.length - 1)) {
                    i++;
                    filePath1 = args[i];
                }
            } else if (arg.startsWith("-image1=")) {
                int idx = args[i].indexOf("=");
                if (idx > -1) {
                    filePath1 =  args[i].substring(idx + 1);
                }
            } else if (arg.equalsIgnoreCase("-image2=")) {
                if (i < (args.length - 1)) {
                    i++;
                    filePath2 = args[i];
                }
            } else if (arg.startsWith("-image2=")) {
                int idx = args[i].indexOf("=");
                if (idx > -1) {
                    filePath2 =  args[i].substring(idx + 1);
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
        
        if (filePath1 == null) {
            System.err.println("-image1= is required");
            System.exit(1);
        }
        if (filePath2 == null) {
            System.err.println("-image2= is required");
            System.exit(1);
        }
        
        CMDLineInflectionUtil cmdLineInvoker = new CMDLineInflectionUtil();
              
        filePath1 = cmdLineInvoker.checkFilePath(filePath1);
        filePath2 = cmdLineInvoker.checkFilePath(filePath2);
      
        if (!writeImageOutput) {
            writeTextOutput = true;
        }
        
        String imageFileName1 = 
            (new File(filePath1)).toPath().getFileName().toString();
        imageFileName1 = imageFileName1.substring(0, imageFileName1.length() -4);
        
        String imageFileName2 = 
            (new File(filePath2)).toPath().getFileName().toString();
        imageFileName2 = imageFileName2.substring(0, imageFileName2.length() -4);
        
        GreyscaleImage img1 = null;
        GreyscaleImage img2 = null;
        try {
            img1 = ImageIOHelper.readImageAsGrayScaleG(filePath1);
            img2 = ImageIOHelper.readImageAsGrayScaleG(filePath2);

            CurvatureScaleSpaceInflectionMapper mapper = new 
                CurvatureScaleSpaceInflectionMapper(img1, img2);

            if (useLineDrawingMode) {
                mapper.useLineDrawingLineMode();
            }

            TransformationParameters transformationParams
                = mapper.createEuclideanTransformation();
                
            List<PairIntArray> edges1 = 
                mapper.getEdges1InOriginalReferenceFrame();

            Transformer transformer = new Transformer();
            
            PairIntArray xy1 = mapper.getMatchedXY1();
            
            PairIntArray xy1Tr = transformer.applyTransformation(
                transformationParams, new PairIntArray[]{xy1})[0];
            
            PairIntArray xy2 = mapper.getMatchedXY2();
            
            
            PairIntArray[] edges1Tr = transformer.applyTransformation(
                transformationParams, 
                edges1.toArray(new PairIntArray[edges1.size()]));
            
            List<PairIntArray> edges2 = 
                mapper.getEdges2InOriginalReferenceFrame();
            
            if (writeTextOutput) {

                cmdLineInvoker.writeMatchingPoints(xy1, xy1Tr, xy2,
                    transformationParams, imageFileName1, imageFileName2);

                cmdLineInvoker.writeTransformation(transformationParams,
                    imageFileName1, imageFileName2);

                cmdLineInvoker.writeEdges(edges1, edges1Tr, 
                    imageFileName1, 
                    imageFileName2);

                // write matching_points.tsv, transformation.txt, 
                //     transformed_edges.tsv
            }

            if (writeImageOutput) {
                 
                 Image image1 = ImageIOHelper.readImageAsGrayScale(filePath1);
                 
                 Image image2 = ImageIOHelper.readImageAsGrayScale(filePath2);
                 
                 // improve this by using rotation
                 int maxWidth = (image2.getWidth() > image1.getWidth()) ?
                     image2.getWidth() : image1.getWidth();
                 int maxHeight = (image2.getHeight() > image1.getHeight()) ?
                     image2.getHeight() : image1.getHeight();
                 Image image1Tr = transformer.applyTransformation(image1, 
                     transformationParams, maxWidth, maxHeight);
                 
                 if (includeEdgesInImageOutput) {
                     
                     ImageIOHelper.addAlternatingColorCurvesToImage(edges1, 
                         image1);
                     
                     ImageIOHelper.addAlternatingColorCurvesToImage(edges1Tr, 
                         image1Tr);
                     
                     ImageIOHelper.addAlternatingColorCurvesToImage(edges2, 
                         image2);
                 }
                 
                 ImageIOHelper.addCurveToImage(xy1, image1, 2, 255, 0, 0);
                 
                 ImageIOHelper.addCurveToImage(xy2, image2, 2, 255, 0, 0);
                 
                 ImageIOHelper.addCurveToImage(xy1Tr, image1Tr, 2, 255, 0, 0);
                 
                 String outFilePath1 = cwd + "inflection_points_" + 
                     imageFileName1 + ".png";
                 
                 ImageIOHelper.writeOutputImage(outFilePath1, image1);
                 
                 String outFilePath2 = cwd + "inflection_points_" + 
                     imageFileName2 + ".png";
                 
                 ImageIOHelper.writeOutputImage(outFilePath2, image2);
                 
                 String outFilePath1Tr = cwd + "transformed_" + 
                     imageFileName1 + ".png";
                 
                 ImageIOHelper.writeOutputImage(outFilePath1Tr, image1Tr);
                 
                 System.out.println(eol + "wrote: " + outFilePath1);
                 System.out.println(eol + "wrote: " + outFilePath2);
                 System.out.println(eol + "wrote: " + outFilePath1Tr);
                 
             }

        } catch (Exception e) {
            e.printStackTrace();
            System.err.println("ERROR: " + e.getMessage());
            System.exit(1);
        }
    }
    
    public static void printUsage() {
                
        System.out.println(eol + "Runner to invoke transformation of image1 to");
        System.out.println("reference frame of image 2 via inflection points.");
        System.out.println(eol + "Required arguments:");
        System.out.println("  -image1=path_to_image_1");
        System.out.println("  -image2=path_to_image_2");
        System.out.println(eol + "Optional arguments:");
        System.out.println("  -text_output");
        System.out.println("      writes a file called matching_points.tsv in the");
        System.out.println("      current directory.  Also writes transformation.txt");
        System.out.println("      and transformed_edges.tsv");
        System.out.println("  -mark_image");
        System.out.println("      writes 2 files called matching_points.png in the");
        System.out.println("      current directory. It's the images with ");
        System.out.println("      inflection points overplotted.");
        System.out.println("  -include_edges");
        System.out.println("      in the written matching_points.png, includes the");
        System.out.println("      edges in the image.");
        System.out.println("  -input_is_line_drawing");
        System.out.println("      changes internal logic to handle an image");
        System.out.println("      of lines or solid blocks.");
        System.out.println(eol
            + "Note that default is to create matching_points.tsv in the currect directory");
        System.out.println(eol);
    }

}
