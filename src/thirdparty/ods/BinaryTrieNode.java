package thirdparty.ods;

import java.lang.reflect.Array;

/**
 *
 * @author nichole
 */
public class BinaryTrieNode<T> {
    T x;
    BinaryTrieNode<T> parent = null;
    BinaryTrieNode<T>[] child = null;
    BinaryTrieNode<T> jump = null;

    @SuppressWarnings({"unchecked", "rawtypes"})
    public BinaryTrieNode() {
        Class cls = this.getClass();
        child = (BinaryTrieNode<T>[]) Array.newInstance(cls, 2);
    }
  
    @Override
    public String toString() {
        return "{" + String.valueOf(x) + "}";
    }
    
    public String toString2() {
        
        StringBuilder sb = new StringBuilder();
        sb.append("[hash=").append(Integer.toString(hashCode()));
        sb.append(" x=");
        if (x != null) {
            sb.append(x);
        }
        sb.append(" p=");
        if (parent != null) {
            sb.append(Integer.toString(parent.hashCode()));
        }
        sb.append(" c[0]=");
        if (child[0] != null) {
            sb.append(Integer.toString(child[0].hashCode()));
        }
        sb.append(" c[1]=");
        if (child[1] != null) {
            sb.append(Integer.toString(child[1].hashCode()));
        }
        sb.append(" jmp=");
        if (jump != null) {
            sb.append(Integer.toString(jump.hashCode()));
        }
        sb.append("]");
        
        return sb.toString();
    }

}
