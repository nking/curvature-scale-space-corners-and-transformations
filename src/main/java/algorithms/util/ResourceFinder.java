package algorithms.util;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileWriter;
import java.io.IOException;
import java.net.URL;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class ResourceFinder {

    /**
     *
     */
    protected final static String sep = System.getProperty("file.separator");

    /**
     *
     */
    protected final static Logger log = Logger.getLogger(ResourceFinder.class.getName());

    /**
     *
     * @param fileName
     * @return
     * @throws IOException
     */
    public static String findFileInResources(String fileName) throws IOException {

        String dirPath = findResourcesDirectory();

        String filePath = dirPath + sep + fileName;

        File f = new File(filePath);
        if (!f.exists()) {
            throw new IOException("could not find file at " + filePath);
        }
        return filePath;
    }

    /**
     *
     * @param subDir
     * @param fileName
     * @return
     * @throws IOException
     */
    public static String findFileInResources(String subDir, String fileName) throws IOException {

        String dirPath = findResourcesDirectory(subDir);

        String filePath = dirPath + sep + fileName;

        File f = new File(filePath);
        if (!f.exists()) {
            throw new IOException("could not find file at " + filePath);
        }
        return filePath;
    }

    /**
     *
     * @param subDir
     * @return
     * @throws IOException
     */
    public static String findResourcesDirectory(String subDir) throws IOException {

        String dirPath = findResourcesDirectory();
        
        String filePath = dirPath + sep + subDir;
        
        File f = new File(filePath);
        if (!f.exists()) {
            throw new IOException("could not find resources sub-directory named " + subDir);
        }
        
        return filePath;
    }

    /**
     *
     * @return
     * @throws IOException
     */
    public static String findResourcesDirectory() throws IOException {

        String srcPath = findDirectory("src");

        String path = srcPath + sep + "main" + sep + "resources";

        File f = new File(path);
        if (!f.exists()) {
            throw new IOException("could not find directory named src/main/resources");
        }
        return path;
    }

    /**
     *
     * @param dirName
     * @return
     * @throws IOException
     */
    public static String findDirectory(String dirName) throws IOException {

        ClassLoader cls = ResourceFinder.class.getClassLoader();

        URL url = cls.getResource(".");
        if (url == null) {
            throw new IOException("base path not found for " + dirName);
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
                    // try test/resources
                    filePath = url.getPath() + ".." + sep + ".." + sep + "src"
                        + sep + "test" + sep + "resources" + sep + dirName;
                    f = new File(filePath);
                    if (!f.exists()) {
                        throw new IOException("could not find directory named " + dirName);
                    }
                }
            }
        }
        return filePath;
    }

    /**
     *
     * @return
     * @throws IOException
     */
    public static String findTestDirectory() throws IOException {

        String srcPath = findDirectory("src");

        String filePath = srcPath + sep + "test";

        File f = new File(filePath);
        if (!f.exists()) {
            throw new IOException("could not find directory named src/test");
        }
        return filePath;
    }

    /**
     *
     * @param fileName
     * @return
     * @throws IOException
     */
    public static String getAFilePathInTestResources(String fileName) throws IOException {

        String dirPath = findTestDirectory();
        String filePath = dirPath + sep + "resources" + sep + fileName;

        return filePath;
    }

    /**
     *
     * @param fileName
     * @return
     * @throws IOException
     */
    public static String findFileInTestResources(String fileName) throws IOException {

        String filePath = getAFilePathInTestResources(fileName);

        File f = new File(filePath);
        if (!f.exists()) {
            throw new IOException("could not find file at " + filePath);
        }
        return filePath;
    }
    
    /**
     *
     * @param fileNamePrefix
     * @return
     * @throws IOException
     */
    public static File[] findFilesInTestResources(String fileNamePrefix) throws IOException {
        
        String dirPath = findTestDirectory() + sep + "resources";
        
        File resDir = new File(dirPath);

        File[] files = resDir.listFiles(new PrefixFileFilter(fileNamePrefix));
        
        return files;
    }
    
    /**
     *
     */
    protected static class PrefixFileFilter implements FileFilter {

        /**
         *
         */
        protected final String prefix;

        /**
         *
         * @param prefixForFileName
         */
        public PrefixFileFilter(String prefixForFileName) {
            this.prefix = prefixForFileName;
        }
        @Override
        public boolean accept(File pathname) {
            if (pathname.getName().startsWith(prefix)) {
                return true;
            }
            return false;
        }
    }

    /**
     *
     * @param fileName
     * @return
     * @throws IOException
     */
    public static String getAFilePathInCWD(String fileName) throws IOException {

        ClassLoader cls = ResourceFinder.class.getClassLoader();

        URL basedir = cls.getResource(".");

        if (basedir != null) {
            return basedir.getPath() + sep + fileName;
        }

        String cwd = System.getProperty("user.dir");

        if (cwd == null) {
            throw new IOException("base path not found for " + fileName);
        }

        String filePath = cwd + sep + fileName;

        return filePath;
    }

    /**
     *
     * @param fileContent
     * @param fileName
     * @return
     * @throws IOException
     */
    public static String writeToTestResources(String fileContent, String fileName) throws IOException {

        String filePath = getAFilePathInTestResources(fileName);

        return writeDataToDirectory(fileContent, filePath);
    }

    /**
     *
     * @param fileContent
     * @param fileName
     * @return
     * @throws IOException
     */
    public static String writeToCWD(String fileContent, String fileName) throws IOException {

        String filePath = getAFilePathInCWD(fileName);

        return writeDataToDirectory(fileContent, filePath);
    }

    /**
     *
     * @param fileContent
     * @param fileName
     * @return
     * @throws IOException
     */
    public static String writeToTmpData(String fileContent, String fileName) throws IOException {

        String filePath = getAFilePathInTmpData(fileName);

        return writeDataToDirectory(fileContent, filePath);
    }

    /**
     *
     * @param fileContent
     * @param filePath
     * @return
     * @throws IOException
     */
    protected static String writeDataToDirectory(String fileContent, String filePath) throws IOException {

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
            log.info( "wrote to file: " + filePath );
        }
        return filePath;
    }

    /**
     *
     * @return
     * @throws IOException
     */
    public static String findTmpDataDirectory() throws IOException {

        return findDirectory("tmpdata");
    }

    /**
     *
     * @param fileName
     * @return
     * @throws IOException
     */
    public static File findFileInTmpData(String fileName) throws IOException {

        String filePath = getAFilePathInTmpData(fileName);

        File fl = new File(filePath);
        if (!fl.exists()) {
            throw new IOException("could not find file at " + filePath);
        }
        return fl;
    }

    /**
     *
     * @param fileName
     * @return
     * @throws IOException
     */
    public static String getAFilePathInTmpData(String fileName) throws IOException {

        String baseDir = findTmpDataDirectory();

        String filePath = baseDir + sep + fileName;

        return filePath;
    }

}
