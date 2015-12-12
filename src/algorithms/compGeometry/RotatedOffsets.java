package algorithms.compGeometry;

/**
 * a class to dynamically populate upon need and cache offset arrays for rotated
 * offsets used in the FeatureMatcher. Note that the methods are not thread safe
 * - they reuse an internal array to avoid constructing a new one for each fetch
 * for a rotation. The user should not allow more than one thread to use the
 * getXOffsets or getYOffsets simultaneously. If that is needed, the class can
 * be altered to be thread safe, but incur the cost of an array construction for
 * each get and a block to check existence of internal data.
 *
 * @author nichole
 */
public class RotatedOffsets {

    private static RotatedOffsets instance = null;

    private static int[][] offsets360X;
    private static int[][] offsets360Y;

    private static long[][] offsets360XL;
    private static long[][] offsets360YL;

    private static int cellDimension = 2;
    private static int range0 = 6;

    private final boolean is64Bit;

    /**
     * <pre>
     * defaults in IntensityFeatures: 
     * cellDim=2
     * nCellsAcross=6 
     * range0 = 2 * (nCellsAcross/2) = 6
     * length of array = (2*nCellsAcross) * (2*nCellsAcross) = 144
     * </pre>
     */
    private int[] xOffsets = new int[range0 * range0];
    private int[] yOffsets = new int[xOffsets.length];

    /**
     * range0, that is 6, can be stored in 3 bits, and the negative portion
     * requires one more bit, so 4 bits total to store an offset.
     */
    private static int itemBits = 4;

    private RotatedOffsets() {

        String arch = System.getProperty("sun.arch.data.model");

        is64Bit = ((arch != null) && arch.equals("64")) ? true : false;

        if (is64Bit) {
            offsets360XL = new long[360][];
            offsets360YL = new long[360][];
            offsets360X = null;
            offsets360Y = null;
        } else {
            offsets360X = new int[360][];
            offsets360Y = new int[360][];
            offsets360XL = null;
            offsets360YL = null;
        }
    }

    public static synchronized RotatedOffsets getInstance() {
        if (instance == null) {
            instance = new RotatedOffsets();
        }
        return instance;
    }

    /**
     * get the offset array for rotation by rotationInDegrees. Note that this
     * method is not thread safe - it reuses an internal array to avoid
     * constructing a new one for each invocation. The user should not allow
     * more than one thread to use getXOffsets simultaneously or more than one
     * return from the method to be retained simultaneously because both returns
     * will be the same instance.
     *
     * @param rotationInDegrees
     * @return
     */
    public int[] getXOffsets(int rotationInDegrees) {

        checkBounds(rotationInDegrees);

        populateIfMissingX(rotationInDegrees);

        if (is64Bit) {

            long[] compresedOffsets = offsets360XL[rotationInDegrees];

            populateOffsets(compresedOffsets, xOffsets);

        } else {

            int[] compresedOffsets = offsets360X[rotationInDegrees];

            populateOffsets(compresedOffsets, xOffsets);
        }

        return xOffsets;
    }

    private void populateIfMissingX(int rotationInDegrees) {

        if (is64Bit) {
            if (offsets360XL[rotationInDegrees] == null) {
                populateXOffsetsLong(rotationInDegrees);
            }
        } else {
            if (offsets360X[rotationInDegrees] == null) {
                populateXOffsetsInt(rotationInDegrees);
            }
        }

    }

    private void populateIfMissingY(int rotationInDegrees) {

        if (is64Bit) {
            if (offsets360YL[rotationInDegrees] == null) {
                populateYOffsetsLong(rotationInDegrees);
            }
        } else {
            if (offsets360Y[rotationInDegrees] == null) {
                populateYOffsetsInt(rotationInDegrees);
            }
        }

    }

    protected void populateXOffsetsLong(int rotationInDegrees) {

        double rotationInRadians = rotationInDegrees * Math.PI / 180.;
        double mc = Math.cos(rotationInRadians);
        double ms = Math.sin(rotationInRadians);

        //64 bits / itemBits = 16 values per item
        long mask = (1L << itemBits) - 1L;

        // 144 / 16 = 9
        long[] compresedOffsets = new long[9];

        int count = 0;
        for (int dx = -range0; dx < range0; dx += cellDimension) {
            for (int dy = -range0; dy < range0; dy += cellDimension) {

                int count16 = 0;
                long datum16 = 0;

                for (int dxc = 0; dxc < cellDimension; ++dxc) {
                    for (int dyc = 0; dyc < cellDimension; ++dyc) {

                        double xt = (((dx + dxc) * mc) + ((dy + dyc) * ms));

                        int v = (int) Math.round(xt);

                        // if it's a negative number, make the number a positive 
                        // and times two then make the last bit to 1 to make it odd
                        int t = (v < 0) ? ((-1 * v << 1) | 1) : (v << 1);

                        long shift = (long) (itemBits * count16);

                        long shifted = (t & mask) << shift;

                        datum16 += shifted;

                        count16++;

                        if (count16 == 16) {
                            compresedOffsets[count] = datum16;
                            count++;
                            count16 = 0;
                            datum16 = 0;
                        }
                    }
                }
            }
        }

        offsets360XL[rotationInDegrees] = compresedOffsets;
    }

    protected void populateXOffsetsInt(int rotationInDegrees) {

        double rotationInRadians = rotationInDegrees * Math.PI / 180.;
        double mc = Math.cos(rotationInRadians);
        double ms = Math.sin(rotationInRadians);

        //32 bits / itemBits = 8 values per item
        int mask = (1 << itemBits) - 1;

        // 144 / 8 = 18
        int[] compresedOffsets = new int[18];

        int count = 0;
        for (int dx = -range0; dx < range0; dx += cellDimension) {
            for (int dy = -range0; dy < range0; dy += cellDimension) {

                int count8 = 0;
                int datum8 = 0;

                for (int dxc = 0; dxc < cellDimension; ++dxc) {
                    for (int dyc = 0; dyc < cellDimension; ++dyc) {

                        double xt = (((dx + dxc) * mc) + ((dy + dyc) * ms));

                        int v = (int) Math.round(xt);

                        // if it's a negative number, make the number a positive 
                        // and times two then make the last bit to 1 to make it odd
                        int t = (v < 0) ? ((-1 * v << 1) | 1) : (v << 1);

                        int shift = (itemBits * count8);

                        int shifted = (t & mask) << shift;

                        datum8 += shifted;

                        count8++;

                        if (count8 == 8) {
                            compresedOffsets[count] = datum8;
                            count++;
                            count8 = 0;
                            datum8 = 0;
                        }
                    }
                }
            }
        }

        offsets360X[rotationInDegrees] = compresedOffsets;
    }

    protected void populateYOffsetsLong(int rotationInDegrees) {

        double rotationInRadians = rotationInDegrees * Math.PI / 180.;
        double mc = Math.cos(rotationInRadians);
        double ms = Math.sin(rotationInRadians);

        //64 bits / itemBits = 16 values per item
        long mask = (1L << itemBits) - 1L;

        // 144 / 16 = 9
        long[] compresedOffsets = new long[9];

        int count = 0;
        for (int dx = -range0; dx < range0; dx += cellDimension) {
            for (int dy = -range0; dy < range0; dy += cellDimension) {

                int count16 = 0;
                long datum16 = 0;

                for (int dxc = 0; dxc < cellDimension; ++dxc) {
                    for (int dyc = 0; dyc < cellDimension; ++dyc) {

                        double yt = (-((dx + dxc) * ms)) + ((dy + dyc) * mc);

                        int v = (int) Math.round(yt);

                        // if it's a negative number, make the number a positive 
                        // and times two then make the last bit to 1 to make it odd
                        int t = (v < 0) ? ((-1 * v << 1) | 1) : (v << 1);

                        long shift = (long) (itemBits * count16);

                        long shifted = (t & mask) << shift;

                        datum16 += shifted;

                        count16++;

                        if (count16 == 16) {
                            compresedOffsets[count] = datum16;
                            count++;
                            count16 = 0;
                            datum16 = 0;
                        }
                    }
                }
            }
        }

        offsets360YL[rotationInDegrees] = compresedOffsets;
    }

    protected void populateYOffsetsInt(int rotationInDegrees) {

        double rotationInRadians = rotationInDegrees * Math.PI / 180.;
        double mc = Math.cos(rotationInRadians);
        double ms = Math.sin(rotationInRadians);

        //32 bits / itemBits = 8 values per item
        int mask = (1 << itemBits) - 1;

        // 144 / 8 = 18
        int[] compresedOffsets = new int[18];

        int count = 0;
        for (int dx = -range0; dx < range0; dx += cellDimension) {
            for (int dy = -range0; dy < range0; dy += cellDimension) {

                int count8 = 0;
                int datum8 = 0;

                for (int dxc = 0; dxc < cellDimension; ++dxc) {
                    for (int dyc = 0; dyc < cellDimension; ++dyc) {

                        double yt = (-((dx + dxc) * ms)) + ((dy + dyc) * mc);

                        int v = (int) Math.round(yt);

                        // if it's a negative number, make the number a positive 
                        // and times two then make the last bit to 1 to make it odd
                        int t = (v < 0) ? ((-1 * v << 1) | 1) : (v << 1);

                        int shift = (itemBits * count8);

                        int shifted = (t & mask) << shift;

                        datum8 += shifted;

                        count8++;

                        if (count8 == 8) {
                            compresedOffsets[count] = datum8;
                            count++;
                            count8 = 0;
                            datum8 = 0;
                        }
                    }
                }
            }
        }

        offsets360Y[rotationInDegrees] = compresedOffsets;
    }
    
    private void populateOffsets(long[] compresedOffsets, int[] outputOffsets) {

        //64 bits / itemBits = 16 values per item
        long mask = (1L << itemBits) - 1L;

        int count = 0;

        for (int i = 0; i < compresedOffsets.length; ++i) {

            long item = compresedOffsets[i];

            // read each of the 16 datum in the item
            for (int j = 0; j < 16; ++j) {

                long v = (item >> (long) (j * itemBits)) & mask;

                // decode to negative
                long t = ((v & 1) == 1) ? -1 * (v >> 1) : v >> 1;

                outputOffsets[count] = (int) t;

                count++;
            }
        }
    }

    private void populateOffsets(int[] compresedOffsets, int[] outputOffsets) {

        //32 bits / itemBits = 8 values per item
        int mask = (1 << itemBits) - 1;

        int count = 0;

        for (int i = 0; i < compresedOffsets.length; ++i) {

            int item = compresedOffsets[i];

            // read each of the 8 datum in the item
            for (int j = 0; j < 8; ++j) {

                int v = (item >> (j * itemBits)) & mask;

                // decode to negative
                int t = ((v & 1) == 1) ? -1 * (v >> 1) : v >> 1;

                outputOffsets[count] = t;

                count++;
            }
        }
    }

    /**
     * get the offset array for rotation by rotationInDegrees. Note that this
     * method is not thread safe - it reuses an internal array to avoid
     * constructing a new one for each invocation. The user should not allow
     * more than one thread to use getYOffsets simultaneously or more than one
     * return from the method to be retained simultaneously because both returns
     * will be the same instance.
     *
     * @param rotationInDegrees
     * @return
     */
    public int[] getYOffsets(int rotationInDegrees) {

        checkBounds(rotationInDegrees);

        populateIfMissingY(rotationInDegrees);

        if (is64Bit) {

            long[] compresedOffsets = offsets360YL[rotationInDegrees];

            populateOffsets(compresedOffsets, yOffsets);

        } else {

            int[] compresedOffsets = offsets360Y[rotationInDegrees];

            populateOffsets(compresedOffsets, yOffsets);
        }

        return xOffsets;
    }

    private void checkBounds(int rotationInDegrees) {

        if (rotationInDegrees < 0 || rotationInDegrees > 359) {
            throw new IllegalArgumentException(
                "rotationInDegrees must be >=0 and < 360");
        }
    }

}
