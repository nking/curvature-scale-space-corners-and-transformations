package algorithms.util;

import java.io.IOException;

/**
 *
 * @author nichole
 */
public interface IInputFileReader {

    /**
     *
     * @return
     */
    public float[] getX();

    /**
     *
     * @return
     */
    public float[] getXErrors();

    /**
     *
     * @return
     */
    public float[] getY();

    /**
     *
     * @return
     */
    public float[] getYErrors();

    /**
     *
     * @param pathToFile
     * @throws IOException
     */
    public void read(String pathToFile) throws IOException;

}
