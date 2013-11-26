package algorithms.compGeometry.clustering.twopointcorrelation;

import algorithms.util.ResourceFinder;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.OutputStream;
import java.util.Calendar;
import java.util.logging.Logger;

/**
 *
 * a utility class for serializing and de-serializing the indexer contents to the file system
 *
 * @author nichole
 */
public class SerializerUtil {

    public static String findTmpDataDirectory() throws IOException {

        String path = ResourceFinder.findTmpDataDirectory();

        return path;
    }

    public static String serializeIndexer(DoubleAxisIndexer dain) throws IOException {

        String resDirectory = SerializerUtil.findTmpDataDirectory();

        String sep = System.getProperty("file.separator");

        boolean includeErrors = (dain.getXErrors() != null);

        String number = Calendar.getInstance().getTime().toString();
        number = number.replaceAll("\\s+", "");

        String fileName = sep + "indexer_" + number + ".bin";

        if (includeErrors) {
            fileName = sep + "indexer_" + number + "_with_errors.bin";
        }

        String filePath = resDirectory + fileName;

        return serializeIndexer(dain, filePath);
    }

    public static String serializeIndexer(DoubleAxisIndexer dain, String filePath) throws IOException {

        File file = new File(filePath);

        boolean includeErrors = (dain.getXErrors() != null);

        OutputStream fos = null;
        ObjectOutputStream oos = null;

        float[] x = dain.getX();
        float[] y = dain.getY();

        try {
            fos = new FileOutputStream(file);

            oos = new ObjectOutputStream(fos);

            oos.writeInt(dain.nXY);

            if (includeErrors) {
                for (int i = 0; i < dain.getNXY(); i++) {
                    oos.writeFloat( x[i] );
                    oos.writeFloat( y[i] );
                    oos.writeFloat( dain.getXErrors()[i] );
                    oos.writeFloat( dain.getYErrors()[i] );
                    oos.flush();
                }
            } else {
                for (int i = 0; i < dain.getNXY(); i++) {
                    oos.writeFloat( x[i] );
                    oos.writeFloat( y[i] );
                    oos.flush();
                }
            }

            oos.flush();

            Logger.getLogger(SerializerUtil.class.getName()).info("wrote indexer to " + filePath);

            return filePath;

        } finally {
            try {
                if (oos != null) {
                    oos.close();
                }
                if (fos != null) {
                    fos.close();
                }
            } catch (IOException e1) {
                Logger.getLogger(SerializerUtil.class.getName()).severe(e1.getMessage());
            }
        }
    }

    public static DoubleAxisIndexer readPersistedPoints(String persistedPointsFilePath) throws IOException {
        return readPersistedPoints(persistedPointsFilePath, persistedPointsFilePath.contains("with_errors"));
    }

    public static DoubleAxisIndexer readPersistedPoints(String persistedPointsFilePath,
        boolean containsErrorArrays) throws IOException {

        File file = new File(persistedPointsFilePath);
        if (!file.exists()) {
            throw new IOException("Could not find file");
        }

        InputStream fis = null;
        ObjectInputStream ois = null;

        try {
            DoubleAxisIndexer indexer = new DoubleAxisIndexer();

            fis = new FileInputStream(file);

            ois = new ObjectInputStream(fis);

            int nXY = ois.readInt();

            float[] x = new float[nXY];
            float[] y = new float[nXY];

            if (containsErrorArrays) {

                float[] ex = new float[nXY];
                float[] ey = new float[nXY];

                for (int i = 0; i < nXY; i++) {
                    x[i] = ois.readFloat();
                    y[i] = ois.readFloat();
                    ex[i] = ois.readFloat();
                    ey[i] = ois.readFloat();
                    //System.out.println("**i=" + i + " x=" + x[i] + " y=" + y[i] + " errx=" + ex[i] + " erry=" + ey[i]);
                }

                indexer.sortAndIndexXThenY(x, y, ex, ey, nXY);

            } else {

                for (int i = 0; i < nXY; i++) {
                    x[i] = ois.readFloat();
                    y[i] = ois.readFloat();
                    //System.out.println("**i=" + i + " x=" + x[i] + " y=" + y[i]);
                }

                indexer.sortAndIndexXThenY(x, y, nXY);
            }

            return indexer;

        } finally {
            if (ois != null) {
                ois.close();
            }
            if (fis != null) {
                fis.close();
            }
        }
    }

}
