package algorithms.compGeometry.clustering.twopointcorrelation;

import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;

public class RandomRoughRangeSearchVoidFinder extends RoughRangeSearchVoidFinder {
    
    @Override
    protected void findVoidsImpl() {
        
        int nSamples = 10;
        
        int nDivisionsPerSide = 10;
        
        try {
            int n = indexer.getNumberOfPoints();

            SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
            sr.setSeed(System.nanoTime());

            // find nSamples non-overlapping regions to sample the area

            /*
             * divide into a grid of nDivisionsPerSide by nDivisionsPerSide and choose nSamples randomly from them
             *
             * Note that the division is in index space rather than value space to
             * use the indexer.
             *
             */
            int binSize = indexer.getNXY()/nDivisionsPerSide;

            /*    col 0
             *     ||
             *     \/
             *      0 |  1 | 2     <=== row 0
             *    ---------------
             *        |    |
             *      3 |  4 | 5
             *    ---------------
             *        |    |
             *      6 |  7 | 8
             */

            int nSSq = nDivisionsPerSide*nDivisionsPerSide;

            // choices are 0 through nSamples-1
            boolean[] selected = new boolean[nSSq];

            // q=0,1,2,3
            int quadNumber = 0;

            for (int i = 0; i < nSamples; i++) {

                quadNumber++;
                if (quadNumber > 4) {
                    quadNumber = 0;
                }

                // to create more evenly distributed random sampling,
                //     will try for each quadrant
                //   3*nDiv:4*nDiv
                //
                //   2*nDiv:3*nDiv
                //
                //   nDiv:2*nDiv
                //
                //   0:nDiv         nDiv:2*nDiv       2*nDiv:3*nDiv       3*nDiv:4*nDiv
                //
                //   0,1   1,1
                //   0,0   1,0

                /*int bin = sr.nextInt(nSSq);
                while (selected[bin]) {
                    bin = sr.nextInt(nSSq);
                }
                selected[bin] = true;
                int row = (bin/nDivisionsPerSide);
                int col = (bin % nDivisionsPerSide);
                */

                int row = 0;
                int col = 0;
                int bin = 0;

                boolean draw = true;
                while (draw) {
                    int dCol = sr.nextInt(nDivisionsPerSide/2);
                    int dRow = sr.nextInt(nDivisionsPerSide/2);
                    switch(quadNumber) {
                        case 0:
                            col = dCol;
                            row = dRow;
                            break;
                        case 1:
                            col = 2*dCol;
                            row = dRow;
                            break;
                        case 2:
                            col = dCol;
                            row = 2*dRow;
                            break;
                        default:
                            col = 2*dCol;
                            row = 2*dRow;
                            break;
                    }
                    bin = col + (row*nDivisionsPerSide);
                    draw = (selected[bin]);
                }
                selected[bin] = true;

                int startX = col*binSize;
                int endX = startX + binSize;
                int yLo = row*binSize;
                int yHi = yLo + binSize;

                if (debug) {
                    /*log.info("[" + quadNumber + "] " + " bin =" + bin + " nDiv=" + nDivisionsPerSide
                        + " " + String.format("  %4d : %4d", col, row)
                        + " " + String.format("  [X %.4f : %.4f] [Y %.4f : %.4f]",
                        indexer.x[indexer.sortedXIndexes[startX]], indexer.x[indexer.sortedXIndexes[endX]],
                        indexer.y[yLo], indexer.y[yHi])  );*/

                    log.info("processIndexedRegion: " + startX + ":" + endX + ":" + yLo + ":" + yHi);
                }

                
                findVoidsRoughRangeSearch(startX, endX, yLo, yHi, 2, 2);
            }

        } catch (NoSuchAlgorithmException e) {
            throw new IllegalStateException(e);
        }
        
    }
}
