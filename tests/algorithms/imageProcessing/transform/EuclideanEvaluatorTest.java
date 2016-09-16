package algorithms.imageProcessing.transform;

import algorithms.util.PairInt;
import java.util.HashSet;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class EuclideanEvaluatorTest extends TestCase {
    
    public EuclideanEvaluatorTest() {
    }

    public void testCalculateF1Score() {
        
        Set<PairInt> templateSet = new HashSet<PairInt>();
        templateSet.add(new PairInt(4, 4));
        templateSet.add(new PairInt(7, 4));
        templateSet.add(new PairInt(10, 20));
        
        Set<PairInt> set2 = new HashSet<PairInt>();
        set2.add(new PairInt(4, 4));
        set2.add(new PairInt(8, 4));
        set2.add(new PairInt(15, 25));
        
        int tolerance = 3;
        
        //tp = 2
        //fp = 1
        //fn = 1
        //precision = 2/(2 + 1);
        //recall = 2/(2 + 1);
        //f = 2 * precision * recall/(precision + recall) = 0.3333

        EuclideanEvaluator eval = new EuclideanEvaluator();
        
        float expResult = 0.6667F;
        float result = eval.calculateF1Score(templateSet, set2, tolerance);
        
        assertTrue(Math.abs(result - expResult) < 0.01);
    }
    
}
