package algorithms;

/**
 * a minimum priority queue that has 
 * amortized runtime complexity of O(lg2(lg2(n)))
 * for operation extractMin which is O(lg2)n)) in
 * standard priority queue and min heap implementations.
 * 
 * This follows the paper of 
 * ""On RAM Priority Queues
 * by Thorup
 * 
 * An improvement to this can be implemented from
 * the description of multi-level buckets in
 * "Simple Shortest Path Algorithm with Linear Average Time"
 * by Goldberg.
 * 
 * @author nichole
 */
public class PriorityQueueThorup {
    
    /*
    u = 2^w is the size of the universe and eps is any psitive constant
    
    pg 4
    
    priority qieie for small integers
    
    has n keys
    
    the insert and extractMin operations depend upon the bit size of
    the numbers
        ( w / log (n) ) bit integers 
    the runtime complexity is O(log log n)
    
    O(log k)
    two sorted lists of words that each have k keys in them at the most
    
    the two sorted lists can be merged into a single sorted list stored
        in two words
        
    for n >= k, given two lists of n integers that are (w/k) bit
    spread over n/k words,
    they can be merged into 2 * n / k words in time O(n/k * log k)
    
    <
    
    */
}
