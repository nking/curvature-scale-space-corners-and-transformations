package thirdparty.ods;

/**
 *
 * @author nichole
 */
public class XFastTrieNode<T> extends BinaryTrieNode<T> {
    
    int prefix;

    @Override
    public boolean equals(Object u) {
        
        boolean t0 = (u instanceof XFastTrieNode<?>)
             && this.prefix == ((XFastTrieNode<T>) u).prefix;
        
        return t0;
    }

    @Override
    public int hashCode() {
        return prefix;
    }
        
}
