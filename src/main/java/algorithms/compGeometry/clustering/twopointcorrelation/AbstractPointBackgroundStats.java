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
 * abstract class for implementations of IPointBackgroundStats to do some of the
 * boiler plate code such as file persistence
 *
 * @author nichole
 */
public abstract class AbstractPointBackgroundStats implements IPointBackgroundStats {

    /**
     * an estimate of the background density which is a linear density
     */
    protected float backgroundDensity;
    
    /**
     * a minimum error for the backgroundDensity
     */
    protected float backgroundDensityError;

    /**
     *
     */
    protected final AxisIndexer indexer;

    /**
     *
     */
    protected boolean didReadPerisistedIndexer = false;

    /**
     *
     */
    protected boolean debug = false;

    /**
     *
     * @param indexedSortedPoints
     */
    public AbstractPointBackgroundStats(AxisIndexer indexedSortedPoints) {
        if (indexedSortedPoints == null) {
            throw new IllegalArgumentException("indexedSortedPoints cannot be null");
        }
        if (indexedSortedPoints.getX() == null) {
            throw new IllegalArgumentException("indexedSortedPoints has to contain sorted data");
        }
        this.indexer = indexedSortedPoints;
    }
    /**
     * load the indexer from the file and note that point errors are expected to be in the file
     * @param persistedIndexerFilePath
     * @throws IOException
     */
    public AbstractPointBackgroundStats(String persistedIndexerFilePath) throws IOException {
        indexer = SerializerUtil.readPersistedPoints(persistedIndexerFilePath, true);
        didReadPerisistedIndexer = true;
    }

    /**
     *
     * @param turnDebugOn
     */
    public void setDebug(boolean turnDebugOn) {
        this.debug = turnDebugOn;
    }

    /**
     *
     * @throws TwoPointVoidStatsException
     */
    public abstract void calc() throws TwoPointVoidStatsException;

    /**
     * use a previously persisted set of two-point background data to run the code.
     *
     * This is most useful if you've used complete sampling before on a large
     * dataset.
     *
     * It has to be used in conjunction with reading the indexer files
     * from persistence.
     *
     * @param persistedBackgroundFileName
     * @throws TwoPointVoidStatsException
     * @throws IOException
     */
    public void calc(String persistedBackgroundFileName) throws TwoPointVoidStatsException, IOException {

        if (!didReadPerisistedIndexer) {
            throw new IllegalStateException("In order to use a persisted set of point densities, need to read the complementary indexer files too.");
        }
        readTwoPointBackground(persistedBackgroundFileName);

        calculateStats();
    }

    /**
     *
     * @throws TwoPointVoidStatsException
     */
    protected abstract void calculateStats() throws TwoPointVoidStatsException;

    /**
     *
     * @return
     * @throws IOException
     */
    public abstract String persistTwoPointBackground() throws IOException;

    /**
     *
     * @param persistedFilePath
     * @return
     * @throws IOException
     */
    public abstract boolean readTwoPointBackground(String persistedFilePath) throws IOException;

    /**
     *
     * @return
     * @throws IOException
     */
    public String persistIndexer() throws IOException {
        return SerializerUtil.serializeIndexer(indexer);
    }

    /**
     *
     * @param oos
     * @throws IOException
     */
    protected abstract void serializeTwoPointBackground(ObjectOutputStream oos) throws IOException;

    /**
     *
     * @param fileRootName
     * @return
     * @throws IOException
     */
    protected String serializeTwoPointBackground(String fileRootName) throws IOException {

        String number = Calendar.getInstance().getTime().toString();
        number = number.replaceAll("\\s+", "");

        String fileName = fileRootName + number + ".bin";
        String filePath = ResourceFinder.getAFilePathInTmpData(fileName);

        File file = new File(filePath);

        OutputStream fos = null;
        ObjectOutputStream oos = null;

        try {
            fos = new FileOutputStream(file);

            oos = new ObjectOutputStream(fos);

            serializeTwoPointBackground(oos);

            return file.getPath();

        } finally {
            try {
                if (fos != null) {
                    fos.close();
                }
                if (oos != null) {
                    oos.close();
                }
            } catch (IOException e1) {
                Logger.getLogger(SerializerUtil.class.getName()).severe(e1.getMessage());
            }
        }
    }

    /**
     *
     * @param oos
     * @throws IOException
     */
    protected abstract void deserializeTwoPointBackground(ObjectInputStream oos) throws IOException;

    /**
     *
     * @param persistedFilePath
     * @return
     * @throws IOException
     */
    protected boolean deserializeTwoPointBackground(String persistedFilePath) throws IOException {

        File file = new File(persistedFilePath);
        if (!file.exists()) {
            throw new IOException("Could not find file");
        }

        InputStream fis = null;
        ObjectInputStream ois = null;

        try {
            fis = new FileInputStream(file);

            ois = new ObjectInputStream(fis);

            deserializeTwoPointBackground(ois);

            return true;

        } catch (IOException e1) {
            Logger.getLogger(SerializerUtil.class.getName()).severe(e1.getMessage());
        } finally {
            try {
                if (fis != null) {
                    fis.close();
                }
                if (ois != null) {
                    ois.close();
                }
            } catch (IOException e1) {
                Logger.getLogger(SerializerUtil.class.getName()).severe(e1.getMessage());
            }
        }

        return false;
    }

    /**
     * @return the backgroundDensity
     */
    public float getBackgroundDensity() {
        return this.backgroundDensity;
    }

    /**
     * return the error on the background density.  NOTE that most implementations
     * will probably return a minimum error.  the formal error would include other terms.
     * 
     * @return the backgroundDensityError
     */
    public float getBackgroundDensityError() {
        return this.backgroundDensityError;
    }
}