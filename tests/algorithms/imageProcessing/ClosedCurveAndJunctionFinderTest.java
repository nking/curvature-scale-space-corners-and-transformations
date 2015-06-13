package algorithms.imageProcessing;

import algorithms.util.PairIntArray;
import algorithms.util.PairIntArrayWithColor;
import java.util.ArrayList;
import java.util.List;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class ClosedCurveAndJunctionFinderTest {
    
    public ClosedCurveAndJunctionFinderTest() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }

    @Test
    public void testFindClosedCurves() {
        
        List<PairIntArray> edges = new ArrayList<PairIntArray>();
        edges.add(getTestEdge1());
        edges.add(getTestEdge2());
        
        ClosedCurveAndJunctionFinder instance = new ClosedCurveAndJunctionFinder();
        
        instance.findClosedCurves(edges);

        assertTrue(edges.get(0) instanceof PairIntArrayWithColor);
        assertTrue(((PairIntArrayWithColor)edges.get(0)).getColor() == 1);
        
        assertTrue(edges.get(1) instanceof PairIntArrayWithColor);
        assertTrue(((PairIntArrayWithColor)edges.get(1)).getColor() == 1);
        
        edges.clear();
        PairIntArray e1 = getTestEdge1();
        PairIntArray e2 = getTestEdge2();
        e1.reverse();
        e2.reverse();
        edges.add(e1);
        edges.add(e2);
        
        instance.findClosedCurves(edges);

        assertTrue(edges.get(0) instanceof PairIntArrayWithColor);
        assertTrue(((PairIntArrayWithColor)edges.get(0)).getColor() == 1);
        
        assertTrue(edges.get(1) instanceof PairIntArrayWithColor);
        assertTrue(((PairIntArrayWithColor)edges.get(1)).getColor() == 1);
        
        edges.clear();
        PairIntArray notClosed = new PairIntArray();
        edges.add(notClosed);
        instance = new ClosedCurveAndJunctionFinder();        
        instance.findClosedCurves(edges);        
        assertFalse(edges.get(0) instanceof PairIntArrayWithColor);
        
        //   clockwise circle
        PairIntArray xy = getCWSquareishCircle();
        edges.clear();
        edges.add(xy);
        instance = new ClosedCurveAndJunctionFinder();        
        instance.findClosedCurves(edges);        
        assertTrue(edges.get(0) instanceof PairIntArrayWithColor);
        
        //    counter clockwise circle
        xy = getCWSquareishCircle();
        xy.reverse();
        edges.clear();
        edges.add(xy);
        instance = new ClosedCurveAndJunctionFinder();        
        instance.findClosedCurves(edges);        
        assertTrue(edges.get(0) instanceof PairIntArrayWithColor);
    }
    
    private PairIntArray getCWSquareishCircle() {
        
        PairIntArray xy = new PairIntArray();
        
        xy.add(17, 24-13);
        xy.add(18, 24-13);
        xy.add(19, 24-14);
        xy.add(20, 24-14);
        xy.add(21, 24-15);
        xy.add(22, 24-16);
        xy.add(22, 24-17);
        xy.add(23, 24-18);
        xy.add(23, 24-19);
        xy.add(22, 24-20);
        xy.add(22, 24-21);
        xy.add(21, 24-22);
        xy.add(20, 24-23);
        xy.add(19, 24-23);
        xy.add(18, 24-24);
        xy.add(17, 24-24);
        xy.add(16, 24-23);
        xy.add(15, 24-23);
        xy.add(14, 24-22);
        xy.add(13, 24-21);
        xy.add(13, 24-20);
        xy.add(12, 24-19);
        xy.add(12, 24-18);
        xy.add(13, 24-17);
        xy.add(13, 24-16);
        xy.add(14, 24-15);
        xy.add(15, 24-14);
        xy.add(16, 24-14);
        
        return xy;
    }
    
    private PairIntArray getTestEdge1() {
        
        PairIntArray curve = new PairIntArray();
        curve.add(216, 246);
        curve.add(217, 246);
        curve.add(218, 245);
        curve.add(219, 244);
        curve.add(220, 243);
        curve.add(221, 242);
        curve.add(222, 241);
        curve.add(223, 240);
        curve.add(224, 239);
        curve.add(225, 238);
        curve.add(226, 237);
        curve.add(227, 236);
        curve.add(228, 235);
        curve.add(229, 234);
        curve.add(230, 233);
        curve.add(231, 232);
        curve.add(232, 231);
        curve.add(233, 230);
        curve.add(234, 229);
        curve.add(235, 228);
        curve.add(236, 227);
        curve.add(237, 226);
        curve.add(238, 225);
        curve.add(239, 224);
        curve.add(240, 223);
        curve.add(241, 222);
        curve.add(242, 221);
        curve.add(243, 220);
        curve.add(244, 219);
        curve.add(245, 218);
        curve.add(246, 217);
        curve.add(247, 216);
        curve.add(248, 215);
        curve.add(249, 214);
        curve.add(250, 215);
        curve.add(251, 215);
        curve.add(252, 216);
        curve.add(253, 217);
        curve.add(254, 218);
        curve.add(255, 219);
        curve.add(256, 220);
        curve.add(257, 221);
        curve.add(258, 222);
        curve.add(259, 223);
        curve.add(260, 223);
        curve.add(261, 224);
        curve.add(262, 225);
        curve.add(263, 226);
        curve.add(264, 227);
        curve.add(265, 228);
        curve.add(266, 229);
        curve.add(267, 230);
        curve.add(268, 231);
        curve.add(269, 232);
        curve.add(270, 233);
        curve.add(271, 234);
        curve.add(272, 235);
        curve.add(272, 236);
        curve.add(272, 237);
        curve.add(272, 238);
        curve.add(272, 239);
        curve.add(271, 240);
        curve.add(271, 241);
        curve.add(270, 242);
        curve.add(269, 243);
        curve.add(268, 244);
        curve.add(267, 245);
        curve.add(266, 246);
        curve.add(265, 247);
        curve.add(264, 248);
        curve.add(263, 249);
        curve.add(262, 250);
        curve.add(261, 251);
        curve.add(260, 252);
        curve.add(259, 253);
        curve.add(258, 254);
        curve.add(257, 255);
        curve.add(256, 256);
        curve.add(255, 257);
        curve.add(254, 258);
        curve.add(253, 259);
        curve.add(252, 260);
        curve.add(251, 261);
        curve.add(250, 262);
        curve.add(249, 263);
        curve.add(248, 264);
        curve.add(247, 265);
        curve.add(246, 266);
        curve.add(245, 267);
        curve.add(244, 268);
        curve.add(243, 269);
        curve.add(242, 270);
        curve.add(241, 270);
        curve.add(240, 270);
        curve.add(239, 270);
        curve.add(238, 270);
        curve.add(237, 269);
        curve.add(236, 268);
        curve.add(235, 267);
        curve.add(234, 266);
        curve.add(233, 265);
        curve.add(232, 264);
        curve.add(231, 263);
        curve.add(230, 262);
        curve.add(229, 261);
        curve.add(228, 261);
        curve.add(227, 260);
        curve.add(226, 259);
        curve.add(225, 258);
        curve.add(224, 257);
        curve.add(223, 256);
        curve.add(222, 255);
        curve.add(221, 254);
        curve.add(220, 253);
        curve.add(219, 252);
        curve.add(218, 251);
        curve.add(217, 250);
        curve.add(216, 249);
        curve.add(216, 248);
        
        return curve;
    }
    
    private PairIntArray getTestEdge2() {
        
        PairIntArray curve = new PairIntArray();
        curve.add(257, 32);
        curve.add(258, 33);
        curve.add(259, 34);
        curve.add(260, 35);
        curve.add(261, 36);
        curve.add(262, 37);
        curve.add(263, 38);
        curve.add(264, 39);
        curve.add(265, 40);
        curve.add(266, 41);
        curve.add(267, 42);
        curve.add(268, 43);
        curve.add(269, 44);
        curve.add(270, 45);
        curve.add(271, 46);
        curve.add(272, 47);
        curve.add(273, 48);
        curve.add(274, 48);
        curve.add(275, 49);
        curve.add(276, 50);
        curve.add(277, 51);
        curve.add(278, 52);
        curve.add(279, 53);
        curve.add(280, 54);
        curve.add(281, 55);
        curve.add(282, 56);
        curve.add(283, 57);
        curve.add(284, 58);
        curve.add(285, 59);
        curve.add(286, 60);
        curve.add(287, 61);
        curve.add(288, 62);
        curve.add(289, 63);
        curve.add(290, 64);
        curve.add(291, 65);
        curve.add(292, 65);
        curve.add(293, 66);
        curve.add(294, 67);
        curve.add(295, 68);
        curve.add(296, 69);
        curve.add(297, 70);
        curve.add(298, 71);
        curve.add(299, 72);
        curve.add(300, 73);
        curve.add(301, 74);
        curve.add(302, 75);
        curve.add(303, 76);
        curve.add(304, 77);
        curve.add(305, 78);
        curve.add(306, 79);
        curve.add(307, 80);
        curve.add(308, 81);
        curve.add(309, 82);
        curve.add(310, 83);
        curve.add(311, 84);
        curve.add(312, 85);
        curve.add(313, 86);
        curve.add(314, 87);
        curve.add(314, 88);
        curve.add(314, 89);
        curve.add(314, 90);
        curve.add(313, 91);
        curve.add(313, 92);
        curve.add(312, 93);
        curve.add(311, 94);
        curve.add(311, 95);
        curve.add(310, 96);
        curve.add(309, 97);
        curve.add(308, 98);
        curve.add(308, 99);
        curve.add(307, 100);
        curve.add(306, 101);
        curve.add(305, 102);
        curve.add(305, 103);
        curve.add(304, 104);
        curve.add(303, 105);
        curve.add(302, 106);
        curve.add(302, 107);
        curve.add(301, 108);
        curve.add(300, 109);
        curve.add(299, 110);
        curve.add(299, 111);
        curve.add(298, 112);
        curve.add(297, 113);
        curve.add(297, 114);
        curve.add(296, 115);
        curve.add(295, 115);
        curve.add(294, 115);
        curve.add(293, 115);
        curve.add(292, 115);
        curve.add(291, 115);
        curve.add(290, 114);
        curve.add(289, 113);
        curve.add(288, 112);
        curve.add(287, 111);
        curve.add(286, 110);
        curve.add(285, 109);
        curve.add(284, 108);
        curve.add(283, 107);
        curve.add(282, 106);
        curve.add(281, 106);
        curve.add(280, 105);
        curve.add(279, 104);
        curve.add(278, 103);
        curve.add(277, 102);
        curve.add(276, 101);
        curve.add(275, 100);
        curve.add(274, 99);
        curve.add(273, 98);
        curve.add(272, 97);
        curve.add(271, 96);
        curve.add(270, 95);
        curve.add(269, 94);
        curve.add(268, 93);
        curve.add(267, 93);
        curve.add(266, 92);
        curve.add(265, 91);
        curve.add(264, 90);
        curve.add(263, 89);
        curve.add(262, 88);
        curve.add(261, 87);
        curve.add(260, 86);
        curve.add(259, 85);
        curve.add(258, 84);
        curve.add(257, 83);
        curve.add(256, 82);
        curve.add(255, 81);
        curve.add(254, 80);
        curve.add(253, 79);
        curve.add(252, 78);
        curve.add(251, 77);
        curve.add(250, 76);
        curve.add(249, 75);
        curve.add(248, 74);
        curve.add(247, 73);
        curve.add(246, 72);
        curve.add(245, 71);
        curve.add(244, 70);
        curve.add(243, 69);
        curve.add(242, 68);
        curve.add(241, 67);
        curve.add(240, 66);
        curve.add(240, 65);
        curve.add(239, 64);
        curve.add(239, 63);
        curve.add(239, 62);
        curve.add(240, 61);
        curve.add(240, 60);
        curve.add(241, 59);
        curve.add(241, 58);
        curve.add(242, 57);
        curve.add(242, 56);
        curve.add(243, 55);
        curve.add(244, 54);
        curve.add(244, 53);
        curve.add(245, 52);
        curve.add(245, 51);
        curve.add(246, 50);
        curve.add(246, 49);
        curve.add(247, 48);
        curve.add(247, 47);
        curve.add(248, 46);
        curve.add(248, 45);
        curve.add(249, 44);
        curve.add(249, 43);
        curve.add(250, 42);
        curve.add(251, 41);
        curve.add(251, 40);
        curve.add(252, 39);
        curve.add(252, 38);
        curve.add(253, 37);
        curve.add(253, 36);
        curve.add(254, 35);
        curve.add(254, 34);
        curve.add(254, 33);
        curve.add(254, 32);
        
        return curve;
    }
    
    public static void main(String[] args) {
        
        try {
            
            ClosedCurveAndJunctionFinderTest test = 
                new ClosedCurveAndJunctionFinderTest();
            
            test.testFindClosedCurves();
                        
        } catch (Exception e) {
            e.printStackTrace();
            System.err.println("ERROR: " + e.getMessage());
        }
    }
}
