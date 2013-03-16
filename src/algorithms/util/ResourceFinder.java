package algorithms.util;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.net.URL;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class ResourceFinder {

    protected final static String sep = System.getProperty("file.separator");

    public static String findFileInResources(String fileName) throws IOException {

        String dirPath = findResourcesDirectory();

        String filePath = dirPath + sep + fileName;

        File f = new File(filePath);
        if (!f.exists()) {
            throw new IOException("could not find file at " + filePath);
        }
        return filePath;
    }

    public static String findResourcesDirectory() throws IOException {

        return findDirectory("resources");
    }

    public static String findDirectory(String dirName) throws IOException {

        ClassLoader cls = ResourceFinder.class.getClassLoader();

        URL url = cls.getResource(".");
        if (url == null) {
            throw new IOException("base path not found");
        }

        String filePath = url.getPath() + ".." + sep + dirName;

        File f = new File(filePath);
        if (!f.exists()) {
            filePath = url.getPath() + ".." + sep + ".." + sep + dirName;
            f = new File(filePath);
            if (!f.exists()) {
                filePath = url.getPath() + dirName;
                f = new File(filePath);
                if (!f.exists()) {
                    throw new IOException("could not find directory named " + dirName);
                }
            }
        }
        return filePath;
    }

    public static String findFileInTestResources(String fileName) throws IOException {

        try {

            String dirPath = findDirectory("testresources");
            String filePath = dirPath + sep + fileName;

            File f = new File(filePath);
            if (!f.exists()) {
                throw new IOException("could not find file at " + filePath);
            }
            return filePath;

        } catch (IOException e) {

            String dirPath = findDirectory("tests");
            String filePath = dirPath + sep + fileName;

            File f = new File(filePath);
            if (!f.exists()) {
                throw new IOException("could not find file at " + filePath);
            }
            return filePath;
        }
    }

    public static String findFileInCWD(String serializationFileName) throws IOException {

        ClassLoader cls = ResourceFinder.class.getClassLoader();

        URL basedir = cls.getResource(".");
        if (basedir == null) {
            throw new IOException("base path not found");
        }

        String filePath = basedir.getPath() + sep + serializationFileName;

        return filePath;
    }

    public static String getAFilePathInCWD(String fileName) throws IOException {

        ClassLoader cls = ResourceFinder.class.getClassLoader();

        URL basedir = cls.getResource(".");
        if (basedir == null) {
            throw new IOException("base path not found");
        }

        String filePath = basedir.getPath() + sep + fileName;

        return filePath;
    }

    public static String writeToCWD(String fileContent, String fileName) throws IOException {

        String filePath = getAFilePathInCWD(fileName);

        FileWriter fw = null;
        BufferedWriter writer = null;
        try {
            File file = new File(filePath);
            file.delete();
            file.createNewFile();

            fw = new FileWriter(file);
            writer = new BufferedWriter(fw);
            writer.write(fileContent);

            writer.flush();

        } finally {

            if (writer != null) {
                writer.close();
            }
            if (fw != null) {
                fw.close();
            }
            System.out.println( "wrote to file: " + filePath );
        }
        return filePath;
    }

    public static String findTmpDataDirectory() throws IOException {

        return findDirectory("tmpdata");
    }

    public static File findFileInTmpData(String fileName) throws IOException {

        String filePath = getAFilePathInTmpData(fileName);

        File fl = new File(filePath);
        if (!fl.exists()) {
            throw new IOException("could not find file at " + filePath);
        }
        return fl;
    }

    public static String getAFilePathInTmpData(String fileName) throws IOException {

        String baseDir = findTmpDataDirectory();

        String filePath = baseDir + sep + fileName;

        return filePath;
    }

}
