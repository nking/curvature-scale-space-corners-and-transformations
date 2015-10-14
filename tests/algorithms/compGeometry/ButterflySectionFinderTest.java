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
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
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
        
        int w = 374;
        int h = 517;
        
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
}
