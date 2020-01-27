package algorithms.connected;

import algorithms.util.SimpleLinkedListNode;
import java.util.List;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class StronglyConnectedComponentsTest extends TestCase {
    
    public StronglyConnectedComponentsTest() {
    }
    
    @Override
    protected void setUp() throws Exception {
       
    }

    @Override
    protected void tearDown() throws Exception {
        
    }
    
    public void test0() {
        // test from Cormen et al. Fig 22.9
     
        SimpleLinkedListNode[] g = new SimpleLinkedListNode[8];
        for (int i = 0; i < g.length; ++i) {
            g[i] = new SimpleLinkedListNode();
        }
        // a  b  c  d  e  f  g  h
        // 0  1  2  3  4  5  6  7
        g[0].insert(1);
        
        g[1].insert(5);
        g[1].insert(4);
        g[1].insert(2);
        
        g[2].insert(6);
        g[2].insert(3);
        
        g[3].insert(7);
        g[3].insert(2);
        
        g[4].insert(5);
        g[4].insert(0);
        
        g[5].insert(6);
        
        g[6].insert(7);
        g[6].insert(5);
        
        g[7].insert(7);
        
        StronglyConnectedComponents scc = new StronglyConnectedComponents();
        List<SimpleLinkedListNode> components = scc.find(g);
        
        assertEquals(4, components.size());
        
        // a  b  c  d  e  f  g  h
        // 0  1  2  3  4  5  6  7
        SimpleLinkedListNode c3 = components.get(0);
        assertEquals(7, c3.getKey());
        
        SimpleLinkedListNode c4 = components.get(1);
        assertEquals(6, c4.getKey());
        assertEquals(5, c4.getNext().getKey());
        
        SimpleLinkedListNode c2 = components.get(2);
        assertEquals(2, c2.getKey());
        assertEquals(3, c2.getNext().getKey());
        
        SimpleLinkedListNode c1 = components.get(3);
        assertEquals(0, c1.getKey());
        assertEquals(1, c1.getNext().getKey());
        assertEquals(4, c1.getNext().getNext().getKey());
       
    }
    
}
