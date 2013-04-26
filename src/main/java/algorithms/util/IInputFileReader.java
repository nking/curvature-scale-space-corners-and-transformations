package algorithms.util;

import java.io.IOException;

/**
 *
 * @author nichole
 */
public interface IInputFileReader {

    public float[] getX();

    public float[] getXErrors();

    public float[] getY();

    public float[] getYErrors();

    public void read(String pathToFile) throws IOException;

}
