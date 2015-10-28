package algorithms.compGeometry;

import algorithms.compGeometry.ButterflySectionFinder.Routes;
import algorithms.imageProcessing.Image;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class ButterflySectionFinderTest extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public void testFindButterflySections() throws Exception {
        
        int w = 512;
        int h = 512;
        
        String fileName = "blob_butterfly_01.dat";

        String filePath = ResourceFinder.findFileInTestResources(fileName);

        PairIntArray closedCurve = Misc.deserializePairIntArray(filePath);
        
        assertNotNull(closedCurve);
        
        assertTrue(closedCurve.getN() > 0);
        
        for (int i = 0; i < 3; ++i) {
        
            closedCurve = Misc.deserializePairIntArray(filePath);
            
            PairInt[] expectedR0 = new PairInt[]{
                new PairInt(378,264), new PairInt(377,263), new PairInt(376,263),
                new PairInt(375,263), new PairInt(374,263)
            };
        
            PairInt[] expectedR1 = new PairInt[]{
                new PairInt(374,261), new PairInt(375,262), new PairInt(376,262),
                new PairInt(377,262), new PairInt(378,262)
            };
            
            if (i == 1) {
                // reverse the x points and the order of expected points
                for (int j = 0; j < closedCurve.getN(); ++j) {
                    int x = closedCurve.getX(j);
                    int y = closedCurve.getY(j);
                    x = w - x;
                    closedCurve.set(j, x, y);
                }
                int n = expectedR0.length;
                PairInt[] tmp = new PairInt[n];
                for (int j = 0; j < expectedR0.length; ++j) {
                    int x = expectedR0[j].getX();
                    int y = expectedR0[j].getY();
                    x = w - x;
                    tmp[n - j - 1] = new PairInt(x, y);
                }
                expectedR0 = tmp;
                n = expectedR1.length;
                tmp = new PairInt[n];
                for (int j = 0; j < expectedR1.length; ++j) {
                    int x = expectedR1[j].getX();
                    int y = expectedR1[j].getY();
                    x = w - x;
                    tmp[n - j - 1] = new PairInt(x, y);
                }
                expectedR1 = tmp;
            } else if (i == 2) {
                // rotate by -90
                double theta = -0.5 * Math.PI;
                double cosine = Math.cos(theta);
                double sine = Math.sin(theta);
                int xc = w >> 1;
                int yc = h >> 1;
                // reverse the x points and the order of expected points
                for (int j = 0; j < closedCurve.getN(); ++j) {
                    int x = closedCurve.getX(j);
                    int y = closedCurve.getY(j);
                    double xt = xc + (((x - xc)*cosine) + ((y - yc)*sine));
                    double yt = yc + (-((x - xc)*sine) + ((y - yc)*cosine));
                    closedCurve.set(j, (int)Math.round(xt), (int)Math.round(yt));
                }
                int n = expectedR0.length;
                for (int j = 0; j < expectedR0.length; ++j) {
                    int x = expectedR0[j].getX();
                    int y = expectedR0[j].getY();
                    double xt = xc + (((x - xc)*cosine) + ((y - yc)*sine));
                    double yt = yc + (-((x - xc)*sine) + ((y - yc)*cosine));
                    expectedR0[j] = new PairInt((int)Math.round(xt), (int)Math.round(yt));
                }
                n = expectedR1.length;
                for (int j = 0; j < expectedR1.length; ++j) {
                    int x = expectedR1[j].getX();
                    int y = expectedR1[j].getY();
                    double xt = xc + (((x - xc)*cosine) + ((y - yc)*sine));
                    double yt = yc + (-((x - xc)*sine) + ((y - yc)*cosine));
                    expectedR1[j] = new PairInt((int)Math.round(xt), (int)Math.round(yt));
                }
                PairInt[] swap = expectedR0;
                expectedR0 = expectedR1;
                expectedR1 = swap;
            }
            
            Image img = new Image(w, h);
            MiscDebug.writeImage(closedCurve, img, 0, "_butterfly3_" + i);
            
            ButterflySectionFinder finder = new ButterflySectionFinder();
        
            List<Routes> sections = finder.findButterflySections(closedCurve);
        
            assertTrue(sections.size() == 1);
        
            Routes routes = sections.get(0);
            
            log.fine("i=" + i);
            
            Iterator<PairInt> iter = routes.getRoute0().iterator();
            StringBuilder sb = new StringBuilder("r0: ");
            while (iter.hasNext()) {
                sb.append(iter.next()).append(" ");
            }
            iter = routes.getRoute1().iterator();
            sb.append("\nr1: ");
            while (iter.hasNext()) {
                sb.append(iter.next()).append(" ");
            }
            sb.append("\nexpected r0:").append(Arrays.toString(expectedR0));
            sb.append("\nexpected r1:").append(Arrays.toString(expectedR1));
            log.fine(sb.toString());
        
            Iterator<PairInt> r = routes.getRoute0().iterator();
            int nIter = 0;
            while (r.hasNext()) {
                PairInt p = r.next();
                PairInt pExpected = expectedR0[nIter];
                assertEquals(pExpected, p);
                nIter++;
            }

            r = routes.getRoute1().iterator();
            nIter = 0;
            while (r.hasNext()) {
                PairInt p = r.next();
                PairInt pExpected = expectedR1[nIter];
                assertEquals(pExpected, p);
                nIter++;
            }
        }
    }
    
    public void testFindButterflySections3() throws Exception {
        
        int w = 258;
        int h = 187;
        String fileName = "blob_butterfly_03.dat";

        String filePath = ResourceFinder.findFileInTestResources(fileName);

        PairIntArray closedCurve = Misc.deserializePairIntArray(filePath);
        
        assertNotNull(closedCurve);
        
        assertTrue(closedCurve.getN() > 0);
        
        for (int i = 0; i < 3; ++i) {
        
            closedCurve = Misc.deserializePairIntArray(filePath);
            
            PairInt[] expectedR0 = new PairInt[]{
                new PairInt(246,126), new PairInt(247,127), new PairInt(248,128),
                new PairInt(249,128)
            };
        
            PairInt[] expectedR1 = new PairInt[]{
                new PairInt(249,130), new PairInt(248,130), new PairInt(247,129),
                new PairInt(246,129)
            };
            
            if (i == 1) {
                // reverse the x points and the order of expected points
                for (int j = 0; j < closedCurve.getN(); ++j) {
                    int x = closedCurve.getX(j);
                    int y = closedCurve.getY(j);
                    x = w - x;
                    closedCurve.set(j, x, y);
                }
                int n = expectedR0.length;
                PairInt[] tmp = new PairInt[n];
                for (int j = 0; j < expectedR0.length; ++j) {
                    int x = expectedR0[j].getX();
                    int y = expectedR0[j].getY();
                    x = w - x;
                    tmp[n - j - 1] = new PairInt(x, y);
                }
                expectedR0 = tmp;
                n = expectedR1.length;
                tmp = new PairInt[n];
                for (int j = 0; j < expectedR1.length; ++j) {
                    int x = expectedR1[j].getX();
                    int y = expectedR1[j].getY();
                    x = w - x;
                    tmp[n - j - 1] = new PairInt(x, y);
                }
                expectedR1 = tmp;
            } else if (i == 2) {
                // rotate by -90
                double theta = -0.5 * Math.PI;
                double cosine = Math.cos(theta);
                double sine = Math.sin(theta);
                int xc = w >> 1;
                int yc = h >> 1;
                // reverse the x points and the order of expected points
                for (int j = 0; j < closedCurve.getN(); ++j) {
                    int x = closedCurve.getX(j);
                    int y = closedCurve.getY(j);
                    double xt = xc + (((x - xc)*cosine) + ((y - yc)*sine));
                    double yt = yc + (-((x - xc)*sine) + ((y - yc)*cosine));
                    closedCurve.set(j, (int)Math.round(xt), (int)Math.round(yt));
                }
                int n = expectedR0.length;
                for (int j = 0; j < expectedR0.length; ++j) {
                    int x = expectedR0[j].getX();
                    int y = expectedR0[j].getY();
                    double xt = xc + (((x - xc)*cosine) + ((y - yc)*sine));
                    double yt = yc + (-((x - xc)*sine) + ((y - yc)*cosine));
                    expectedR0[j] = new PairInt((int)Math.round(xt), (int)Math.round(yt));
                }
                n = expectedR1.length;
                for (int j = 0; j < expectedR1.length; ++j) {
                    int x = expectedR1[j].getX();
                    int y = expectedR1[j].getY();
                    double xt = xc + (((x - xc)*cosine) + ((y - yc)*sine));
                    double yt = yc + (-((x - xc)*sine) + ((y - yc)*cosine));
                    expectedR1[j] = new PairInt((int)Math.round(xt), (int)Math.round(yt));
                }
            }
            
            Image img = new Image(w, h);
            MiscDebug.writeImage(closedCurve, img, 0, "_butterfly3_" + i);
            
            ButterflySectionFinder finder = new ButterflySectionFinder();
        
            List<Routes> sections = finder.findButterflySections(closedCurve);
        
            assertTrue(sections.size() == 1);
        
            Routes routes = sections.get(0);
            
            Iterator<PairInt> iter = routes.getRoute0().iterator();
            StringBuilder sb = new StringBuilder("r0: ");
            while (iter.hasNext()) {
                sb.append(iter.next()).append(" ");
            }
            iter = routes.getRoute1().iterator();
            sb.append("\nr1: ");
            while (iter.hasNext()) {
                sb.append(iter.next()).append(" ");
            }
            sb.append("\nexpected r0:").append(Arrays.toString(expectedR0));
            sb.append("\nexpected r1:").append(Arrays.toString(expectedR1));
            log.fine(sb.toString());
        
            Iterator<PairInt> r = routes.getRoute0().iterator();
            int nIter = 0;
            while (r.hasNext()) {
                PairInt p = r.next();
                PairInt pExpected = expectedR0[nIter];
                assertEquals(pExpected, p);
                nIter++;
            }

            r = routes.getRoute1().iterator();
            nIter = 0;
            while (r.hasNext()) {
                PairInt p = r.next();
                PairInt pExpected = expectedR1[nIter];
                assertEquals(pExpected, p);
                nIter++;
            }
        }
    }
   
    public void testFindButterflySections2() throws Exception {
        
        int w = 525;
        int h = 525;
        
        String fileName = "tmp_blob_777095183_374_517.dat";

        String filePath = ResourceFinder.findFileInTestResources(fileName);

        PairIntArray closedCurve = Misc.deserializePairIntArray(filePath);
        
        assertNotNull(closedCurve);
        
        assertTrue(closedCurve.getN() > 0);
        
        for (int i = 0; i < 3; ++i) {
        
            closedCurve = Misc.deserializePairIntArray(filePath);
            
            PairInt[] expectedR0 = new PairInt[]{
                new PairInt(256,8), new PairInt(255,9), new PairInt(256,10),
                new PairInt(257,10)
            };
        
            PairInt[] expectedR1 = new PairInt[]{
                new PairInt(255,11), new PairInt(255,10), new PairInt(254,9),
                new PairInt(253,8)
            };
            
            if (i == 1) {
                // reverse the x points and the order of expected points
                for (int j = 0; j < closedCurve.getN(); ++j) {
                    int x = closedCurve.getX(j);
                    int y = closedCurve.getY(j);
                    x = w - x;
                    closedCurve.set(j, x, y);
                }
                for (int j = 0; j < expectedR0.length; ++j) {
                    int x = expectedR0[j].getX();
                    int y = expectedR0[j].getY();
                    x = w - x;
                    expectedR0[j] = new PairInt(x, y);
                }
                for (int j = 0; j < expectedR1.length; ++j) {
                    int x = expectedR1[j].getX();
                    int y = expectedR1[j].getY();
                    x = w - x;
                    expectedR1[j] = new PairInt(x, y);
                }
            } else if (i == 2) {
                // rotate by -90
                double theta = -0.5 * Math.PI;
                double cosine = Math.cos(theta);
                double sine = Math.sin(theta);
                int xc = w >> 1;
                int yc = h >> 1;
                // reverse the x points and the order of expected points
                for (int j = 0; j < closedCurve.getN(); ++j) {
                    int x = closedCurve.getX(j);
                    int y = closedCurve.getY(j);
                    double xt = xc + (((x - xc)*cosine) + ((y - yc)*sine));
                    double yt = yc + (-((x - xc)*sine) + ((y - yc)*cosine));
                    closedCurve.set(j, (int)Math.round(xt), (int)Math.round(yt));
                }
                int n = expectedR0.length;
                for (int j = 0; j < expectedR0.length; ++j) {
                    int x = expectedR0[j].getX();
                    int y = expectedR0[j].getY();
                    double xt = xc + (((x - xc)*cosine) + ((y - yc)*sine));
                    double yt = yc + (-((x - xc)*sine) + ((y - yc)*cosine));
                    expectedR0[j] = new PairInt((int)Math.round(xt), (int)Math.round(yt));
                }
                n = expectedR1.length;
                for (int j = 0; j < expectedR1.length; ++j) {
                    int x = expectedR1[j].getX();
                    int y = expectedR1[j].getY();
                    double xt = xc + (((x - xc)*cosine) + ((y - yc)*sine));
                    double yt = yc + (-((x - xc)*sine) + ((y - yc)*cosine));
                    expectedR1[j] = new PairInt((int)Math.round(xt), (int)Math.round(yt));
                }
            }
            
            Image img = new Image(w, h);
            MiscDebug.writeImage(closedCurve, img, 0, "_butterfly2_" + i);
            
            ButterflySectionFinder finder = new ButterflySectionFinder();
        
            List<Routes> sections = finder.findButterflySections(closedCurve);
        
            assertTrue(sections.size() == 1);
        
            Routes routes = sections.get(0);
            
            log.fine("i=" + i);
            
            Iterator<PairInt> iter = routes.getRoute0().iterator();
            StringBuilder sb = new StringBuilder("r0: ");
            while (iter.hasNext()) {
                sb.append(iter.next()).append(" ");
            }
            iter = routes.getRoute1().iterator();
            sb.append("\nr1: ");
            while (iter.hasNext()) {
                sb.append(iter.next()).append(" ");
            }
            sb.append("\nexpected r0:").append(Arrays.toString(expectedR0));
            sb.append("\nexpected r1:").append(Arrays.toString(expectedR1));
            log.fine(sb.toString());
        
            Iterator<PairInt> r = routes.getRoute0().iterator();
            int nIter = 0;
            while (r.hasNext()) {
                PairInt p = r.next();
                PairInt pExpected = expectedR0[nIter];
                assertEquals(pExpected, p);
                nIter++;
            }

            r = routes.getRoute1().iterator();
            nIter = 0;
            while (r.hasNext()) {
                PairInt p = r.next();
                PairInt pExpected = expectedR1[nIter];
                assertEquals(pExpected, p);
                nIter++;
            }
        }
    }
    
    public void testFindButterflySectionsLargeHoriz() throws Exception {
        
        /*
                                      #      7
                # #  #              #   #    6  <--R0 starts
              #         #  #  # # #     #    5
                # #  #  #  #  # # # # #      4  -->R1 ends
        
            4 5 6 7  8  9  0  1 2 3 4 5 6
        */
        
        for (int i = 0; i < 5; i++) {
            
            PairIntArray closedCurve = null;
            
            if (i == 0) {
                closedCurve = getPattern10();
            } else if (i == 1) {
                closedCurve = getPattern10SwapY();
            } else if (i == 2) {
                closedCurve = getPattern10Crossed();
            } else if (i == 3) {
                closedCurve = getPattern10WrapAround();
            } else {
                closedCurve = getPattern10CrossedSwapY();
            }

            ButterflySectionFinder finder = new ButterflySectionFinder();

            List<Routes> sections = finder.findButterflySections(closedCurve);

            assertTrue(sections.size() == 1);

            Routes routes = sections.get(0);

            assertTrue(routes.route0.size() == 7);
            assertTrue(routes.route1.size() == 7);
            assertNotNull(routes.ep0);
            assertNotNull(routes.ep0End);
            assertNotNull(routes.ep1);
            assertNotNull(routes.ep1End);

            PairInt[] expectedR1 = new PairInt[]{
                new PairInt(8, 4), new PairInt(9, 4), new PairInt(10, 4),
                new PairInt(11, 4), new PairInt(12, 4), new PairInt(13, 4),
                new PairInt(14, 4)
            };
            PairInt[] expectedR0 = new PairInt[]{
                new PairInt(14, 6), new PairInt(13, 5), new PairInt(12, 5),
                new PairInt(11, 5), new PairInt(10, 5), new PairInt(9, 5),
                new PairInt(8, 6)
            };

            if ((i == 1) || (i == 4)) {
                /*
                    # #  #  #  #  # # # # #      4  <-- r0
                  #         #  #  # # #     #    3
                    # #  #              #   #    2  --> r1
                                          #      1
                4 5 6 7  8  9  0  1 2 3 4 5 6
                */
                // expectedR1 is reverse of r0 then swap y: y=8-y
                PairInt[] tmp0 = new PairInt[expectedR1.length];
                for (int ii = 0; ii < tmp0.length; ++ii) {
                    PairInt t = expectedR1[expectedR1.length - 1 - ii];
                    tmp0[ii] = new PairInt(t.getX(), (2*4) - t.getY());
                }
                PairInt[] tmp1 = new PairInt[expectedR0.length];
                for (int ii = 0; ii < tmp1.length; ++ii) {
                    PairInt t = expectedR0[expectedR0.length - 1 - ii];
                    tmp1[ii] = new PairInt(t.getX(), (2*4) - t.getY());
                }
                expectedR0 = tmp0;
                expectedR1 = tmp1;
            }
            
            assertEquals(routes.ep0, expectedR0[0]);
            assertEquals(routes.ep0End, expectedR0[expectedR0.length - 1]);
            assertEquals(routes.ep1, expectedR1[0]);
            assertEquals(routes.ep1End, expectedR1[expectedR1.length - 1]);

            Iterator<PairInt> r = routes.getRoute0().iterator();
            int nIter = 0;
            while (r.hasNext()) {
                PairInt p = r.next();
                PairInt pExpected = expectedR0[nIter];
                assertEquals(pExpected, p);
                nIter++;
            }

            r = routes.getRoute1().iterator();
            nIter = 0;
            while (r.hasNext()) {
                PairInt p = r.next();
                PairInt pExpected = expectedR1[nIter];
                assertEquals(pExpected, p);
                nIter++;
            }
        }
        
    }
    
    public void testFindButterflySectionsLargeVert() throws Exception {
        
        /*
              4 5 6 7  8 
                # #         41
              #     #       40
                #   #       39
                  # #       38
                  # #       37
                  # #       36
                  # #       35
                  # #       34
                #   #       33
                #   #       32
                #   #       31           
                  #         30
              4 5 6 7  8 
        */
        
        for (int i = 3; i < 4; i++) {
            
            PairIntArray closedCurve = null;
            
            if (i == 0) {
                closedCurve = getPattern11();
            } else if (i == 1) {
                closedCurve = getPattern11SwapX();
            } else if (i == 2) {
                closedCurve = getPattern11Crossed();
            } else {
                closedCurve = getPattern11CrossedSwapX();
            }

            ButterflySectionFinder finder = new ButterflySectionFinder();

            List<Routes> sections = finder.findButterflySections(closedCurve);

            assertTrue(sections.size() == 1);

            Routes routes = sections.get(0);

            assertTrue(routes.route0.size() == 7);
            assertTrue(routes.route1.size() == 7);
            assertNotNull(routes.ep0);
            assertNotNull(routes.ep0End);
            assertNotNull(routes.ep1);
            assertNotNull(routes.ep1End);

            PairInt[] expectedR1 = new PairInt[]{
                new PairInt(5, 39), new PairInt(6, 38), new PairInt(6, 37),
                new PairInt(6, 36), new PairInt(6, 35), new PairInt(6, 34),
                new PairInt(5, 33)
            };
            PairInt[] expectedR0 = new PairInt[]{
                new PairInt(7, 33), new PairInt(7, 34), new PairInt(7, 35),
                new PairInt(7, 36), new PairInt(7, 37), new PairInt(7, 38),
                new PairInt(7, 39)
            };

            if ((i == 1) || (i == 3)) {
                
                /*
                      4 5 6 7 8 9 10
                              # #     41
                            #     #   40
                            #   #     39
                            # #       38
                            # #       37
                            # #       36
                            # #       35
                            # #       34
                            #   #     33
                            #   #     32
                            #   #     31           
                              #       30
                      4 5 6 7 8 9 10
                */
                // expectedR1 is reverse of r0 then swap y: y=8-y
                PairInt[] tmp0 = new PairInt[expectedR1.length];
                for (int ii = 0; ii < tmp0.length; ++ii) {
                    PairInt t = expectedR1[expectedR1.length - 1 - ii];
                    tmp0[ii] = new PairInt(14 - t.getX(), t.getY());
                }
                PairInt[] tmp1 = new PairInt[expectedR0.length];
                for (int ii = 0; ii < tmp1.length; ++ii) {
                    PairInt t = expectedR0[expectedR0.length - 1 - ii];
                    tmp1[ii] = new PairInt(14 - t.getX(), t.getY());
                }
                expectedR0 = tmp0;
                expectedR1 = tmp1;
            }
            
            assertEquals(routes.ep0, expectedR0[0]);
            assertEquals(routes.ep0End, expectedR0[expectedR0.length - 1]);
            assertEquals(routes.ep1, expectedR1[0]);
            assertEquals(routes.ep1End, expectedR1[expectedR1.length - 1]);

            Iterator<PairInt> r = routes.getRoute0().iterator();
            int nIter = 0;
            while (r.hasNext()) {
                PairInt p = r.next();
                PairInt pExpected = expectedR0[nIter];
                assertEquals(pExpected, p);
                nIter++;
            }

            r = routes.getRoute1().iterator();
            nIter = 0;
            while (r.hasNext()) {
                PairInt p = r.next();
                PairInt pExpected = expectedR1[nIter];
                assertEquals(pExpected, p);
                nIter++;
            }
        }
        
    }
    
    public void testFindButterflySectionsDiag() throws Exception {        
        
        for (int i = 0; i < 4; i++) {
            
            PairIntArray closedCurve = null;
            
            if (i == 0) {
                closedCurve = getPattern12();
            } else if (i == 1) {
                closedCurve = getPattern12SwapX();
            } else if (i == 2) {
                closedCurve = getPattern12Crossed();
            } else if (i == 3) {
                closedCurve = getPattern12CrossedSwapX();
            }

            ButterflySectionFinder finder = new ButterflySectionFinder();

            List<Routes> sections = finder.findButterflySections(closedCurve);

            assertTrue(sections.size() == 1);

            Routes routes = sections.get(0);

            assertTrue(routes.route0.size() == 7);
            assertTrue(routes.route1.size() == 7);
            assertNotNull(routes.ep0);
            assertNotNull(routes.ep0End);
            assertNotNull(routes.ep1);
            assertNotNull(routes.ep1End);

            PairInt[] expectedR1 = new PairInt[]{
                new PairInt(6, 17), new PairInt(7, 16), new PairInt(8, 15),
                new PairInt(9, 14), new PairInt(10, 13),
                new PairInt(11, 12), new PairInt(12, 11)
            };
            PairInt[] expectedR0 = new PairInt[]{
                new PairInt(12, 13), new PairInt(11, 13), new PairInt(10, 14),
                new PairInt(9, 15), new PairInt(8, 16), new PairInt(7, 17),
                new PairInt(8, 18)
            };
            
            if ((i == 1) || (i == 3)) {
                expectedR1 = new PairInt[]{
                    new PairInt(22, 18), new PairInt(23, 17), new PairInt(22, 16),
                    new PairInt(21, 15), new PairInt(20, 14),
                    new PairInt(19, 13), new PairInt(18, 13)
                };
                expectedR0 = new PairInt[]{
                    new PairInt(19, 12), new PairInt(20, 13), new PairInt(21, 14),
                    new PairInt(22, 15), new PairInt(23, 16), new PairInt(24, 17),
                    new PairInt(25, 18)
                };
            }

            assertEquals(routes.ep0, expectedR0[0]);
            assertEquals(routes.ep0End, expectedR0[expectedR0.length - 1]);
            assertEquals(routes.ep1, expectedR1[0]);
            assertEquals(routes.ep1End, expectedR1[expectedR1.length - 1]);

            Iterator<PairInt> r = routes.getRoute0().iterator();
            int nIter = 0;
            while (r.hasNext()) {
                PairInt p = r.next();
                PairInt pExpected = expectedR0[nIter];
                assertEquals(pExpected, p);
                nIter++;
            }

            r = routes.getRoute1().iterator();
            nIter = 0;
            while (r.hasNext()) {
                PairInt p = r.next();
                PairInt pExpected = expectedR1[nIter];
                assertEquals(pExpected, p);
                nIter++;
            }
        }
        
    }
    
    protected PairIntArray getPattern10() {
        
        /*
                                      #      7
                # #  #              #   #    6
              #         #  #  # # #     #    5
                # #  #  #  #  # # # # #      4     
        
            4 5 6 7  8  9  0  1 2 3 4 5 6
         */
        
        PairIntArray points = new PairIntArray(20);
        points.add(6, 4); points.add(7, 4); points.add(8, 4);
        points.add(9, 4);  points.add(10, 4); points.add(11, 4);
        points.add(12, 4); points.add(13, 4); points.add(14, 4);
        points.add(15, 4);
        
        points.add(5, 5);
        points.add(9, 5); points.add(10, 5); points.add(11, 5);
        points.add(12, 5); points.add(13, 5); points.add(16, 5);
        
        points.add(6, 6);
        points.add(7, 6);
        points.add(8, 6);
        points.add(14, 6);
        points.add(16, 6);
        
        points.add(15, 7);
        
        return points;
    }
    
    protected PairIntArray getPattern10SwapY() {
        
        /*
                # #  #  #  #  # # # # #      4
              #         #  #  # # #     #    3
                # #  #              #   #    2
                                      #      1
            4 5 6 7  8  9  0  1 2 3 4 5 6
         */
        
        PairIntArray p = getPattern10();
        for (int i = 0; i < p.getN(); ++i) {
            int x = p.getX(i);
            int y = (2*4) - p.getY(i);
            p.set(i, x, y);
        }
        
        return p;
    }
    
    protected PairIntArray getPattern10Crossed() {
        
        /*
                                      #      7
                # #  #              #   #    6
              #         #  +  + + #     #    5
                # #  +  +  #  # # + + #      4     
        
            4 5 6 7  8  9  0  1 2 3 4 5 6
         */
        
        PairIntArray points = new PairIntArray(20);
        points.add(6, 4); points.add(7, 4); points.add(8, 4); points.add(9, 4);  
        
        points.add(10, 5); points.add(11, 5); points.add(12, 5); 
        
        points.add(13, 4); points.add(14, 4);
        points.add(15, 4); points.add(16, 5); points.add(16, 6);
        points.add(15, 7); points.add(14, 6);
        
        points.add(13, 5); points.add(12, 4); points.add(11, 4);
        points.add(10, 4); points.add(9, 5);  points.add(8, 6); points.add(7, 6);
        
        points.add(6, 6); points.add(5, 5);
        
        return points;
    }
    
    protected PairIntArray getPattern10WrapAround() {
        
        /*
                                      #      7
                # #  #              #   #    6
              #         #  +  + + #     #    5
                # #  +  +  #  # # + + #      4     
                              ^
                              |
            4 5 6 7  8  9  0  1 2 3 4 5 6
         */
        
        PairIntArray points = new PairIntArray(20);
        points.add(11, 4); points.add(12, 4); points.add(13, 5); 
        points.add(14, 6); points.add(15, 7); points.add(16, 6);
        points.add(16, 5); points.add(15, 4);
        points.add(14, 4); points.add(13, 4);
        points.add(12, 5); points.add(11, 5); points.add(10, 5); 
        
        points.add(9, 4); points.add(8, 4); points.add(7, 4); points.add(6, 4);    
          
        points.add(5, 5); points.add(6, 6); points.add(7, 6);
        points.add(8, 6); points.add(9, 5); points.add(10, 4);        
        
        return points;
    }
    
    protected PairIntArray getPattern10CrossedSwapY() {
        
        PairIntArray p = getPattern10Crossed();
        for (int i = 0; i < p.getN(); ++i) {
            int x = p.getX(i);
            int y = (2*4) - p.getY(i);
            p.set(i, x, y);
        }
        
        return p;
    }
    
    protected PairIntArray getPattern11() {
       
        /*
          4 5 6 7  8 
            # #         41
          #     #       40
            #   #       39
              # #       38
              # #       37
              # #       36
              # #       35
              # #       34
            #   #       33
            #   #       32
            #   #       31           
              #         30
          4 5 6 7  8 
         */
        
        PairIntArray points = new PairIntArray(20);
        points.add(6, 30); points.add(7, 31); points.add(7, 32);
        points.add(7, 33); points.add(7, 34); points.add(7, 35);
        points.add(7, 36); points.add(7, 37); points.add(7, 38); 
        points.add(7, 39); points.add(7, 40);
        
        points.add(6, 41); points.add(5, 41); points.add(4, 40);
        points.add(5, 39);
        points.add(6, 38); points.add(6, 37); points.add(6, 36);
        points.add(6, 35); points.add(6, 34);
        points.add(5, 33); points.add(5, 32); points.add(5, 31);
        return points;
    }

    protected PairIntArray getPattern11Crossed() {
       
        /*
          4 5 6 7  8 
            # #         41
          #     #       40
            #   +       39
              # +       38
              + #       37
              + #       36
              + #       35
              # +       34
            #   #       33
            #   #       32
            #   #       31           
              #         30
          4 5 6 7  8 
        */
        
        PairIntArray points = new PairIntArray(20);
        points.add(6, 30); points.add(7, 31); points.add(7, 32);
        points.add(7, 33); points.add(7, 34); 
        
        points.add(6, 35); points.add(6, 36); points.add(6, 37);
        
        points.add(7, 38); points.add(7, 39); points.add(7, 40);
        points.add(6, 41); points.add(5, 41); points.add(4, 40);
        points.add(5, 39); points.add(6, 38);
        
        points.add(7, 37); points.add(7, 36); points.add(7, 35);
          
        points.add(6, 34); points.add(5, 33); points.add(5, 32); points.add(5, 31);
        
        return points;
    }
    
    private PairIntArray getPattern11SwapX() {
        
        /*
              4 5 6 7 8 9 10
                      # #     41
                    #     #   40
                    #   #     39
                    # #       38
                    # #       37
                    # #       36
                    # #       35
                    # #       34
                    #   #     33
                    #   #     32
                    #   #     31           
                      #       30
              4 5 6 7 8 9 10
        */

        PairIntArray p = getPattern11();
        for (int i = 0; i < p.getN(); ++i) {
            int x = 14 - p.getX(i);
            int y = p.getY(i);
            p.set(i, x, y);
        }
        
        return p;
    }
    
    private PairIntArray getPattern11CrossedSwapX() {
        
        PairIntArray p = getPattern11Crossed();
        for (int i = 0; i < p.getN(); ++i) {
            int x = 14 - p.getX(i);
            int y = p.getY(i);
            p.set(i, x, y);
        }
        
        return p;
    }
    
    protected PairIntArray getPattern12() {
        
        /*
              # #                          20
            #     #                        19
              #     #                      18
                # #                        17
                  # #                      16
                    # #                    15
                      # #                  14
                        # # # #            13
                          #     #          12
                            #     #        11
                              # # #        10
        
            4 5 6 7 8 9 0 1 2 3 4 5 6
         */
        PairIntArray points = new PairIntArray(20);
        points.add(4, 19); points.add(5, 18); points.add(6, 17); points.add(7, 16);
        points.add(8, 15); points.add(9, 14); points.add(10, 13);
        points.add(11, 12); points.add(12, 11);
        points.add(13, 10); points.add(14, 10); points.add(15, 10);
        points.add(15, 11); points.add(14, 12); points.add(13, 13);
        points.add(12, 13); points.add(11, 13);
        
        points.add(10, 14); points.add(9, 15); points.add(8, 16);
        points.add(7, 17); points.add(8, 18); points.add(7, 19);
        points.add(6, 20); points.add(5, 20);
        
        return points;
    }
    
    protected PairIntArray getPattern12Crossed() {
        
        /*
              # #                          20
            #     #                        19
              #     #                      18
                # #                        17
                  # #                      16
                    # #                    15
                      # #                  14
                        + + # #            13
                          #     #          12
                            #     #        11
                              # # #        10
        
            4 5 6 7 8 9 0 1 2 3 4 5 6
         */
        PairIntArray points = new PairIntArray(20);
        
        points.add(11, 12); points.add(12, 11);
        points.add(13, 10); points.add(14, 10); points.add(15, 10);
        points.add(15, 11); points.add(14, 12); points.add(13, 13);
        points.add(12, 13); points.add(11, 13);
        
        points.add(10, 13); points.add(9, 14); points.add(8, 15); 
        points.add(7, 16); points.add(6, 17); points.add(5, 18);         
        points.add(4, 19); points.add(5, 20); points.add(6, 20); 
        points.add(7, 19); points.add(8, 18); points.add(7, 17); 
        points.add(8, 16); points.add(9, 15); points.add(10, 14); 
        
        return points;
    }
    
    private PairIntArray getPattern12SwapX() {
        
        /*
                                                    @ @          20
                                                  @     @        19
                                                @     @          18
                                                  @ @            17
                                                @ @              16
                                              @ @                15
                                            @ @                  14
                                      @ @ @ @                    13
                                    @     @                      12
                                  #     @                        11
                                  # @ @                          10
        
            4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 
                        1                   2                   3
         */

        PairIntArray p = getPattern12();
        for (int i = 0; i < p.getN(); ++i) {
            int x = 30 - p.getX(i);
            int y = p.getY(i);
            p.set(i, x, y);
        }
        
        return p;
    }
    
    private PairIntArray getPattern12CrossedSwapX() {
       
        PairIntArray p = getPattern12Crossed();
        for (int i = 0; i < p.getN(); ++i) {
            int x = 30 - p.getX(i);
            int y = p.getY(i);
            p.set(i, x, y);
        }
        
        return p;
    }
    
    protected PairIntArray getPattern13() {
        
        /*
              # #                          20
            #     #                        19
              #     #                      18
                #   #                      17
                  #   #                    16
                    #   #                  15
                      #   #                14
                        #   # # #          13
                          #       #        12
                            #       #      11
                              # # #        10
        
            4 5 6 7 8 9 0 1 2 3 4 5 6
         */
        
        PairIntArray points = new PairIntArray(20);
        points.add(4, 19); points.add(5, 18); points.add(6, 17); points.add(7, 16);
        points.add(8, 15); points.add(9, 14); points.add(10, 13);
        points.add(11, 12); points.add(12, 11);
        points.add(13, 10); points.add(14, 10); points.add(15, 10);
        points.add(16, 11); points.add(15, 12); points.add(14, 13);
        points.add(13, 13); points.add(12, 13);
        
        points.add(11, 14); points.add(10, 15); points.add(9, 16);
        points.add(8, 17); points.add(8, 18); points.add(7, 19);
        points.add(6, 20); points.add(5, 20);
        
        return points;
    }
}
