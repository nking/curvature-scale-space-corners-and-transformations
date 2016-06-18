package thirdparty.ods;

/*
The class is from the open datastructures source code
http://opendatastructures.org/ods-java.pdf

"The source code available there is released under a Creative Commons
Attribution license, meaning that anyone is free to share: to copy, distribute
and transmit the work; and to remix: to adapt the work, including
the right to make commercial use of the work. The only condition on
these rights is attribution: you must acknowledge that the derived work
contains code and/or text from opendatastructures.org.
https://github.com/patmorin/ods
*/

/**
 * The interface is adapted from the open datastructures source code
http://opendatastructures.org/ods-java.pdf

"The source code available there is released under a Creative Commons
Attribution license, meaning that anyone is free to share: to copy, distribute
and transmit the work; and to remix: to adapt the work, including
the right to make commercial use of the work. The only condition on
these rights is attribution: you must acknowledge that the derived work
contains code and/or text from opendatastructures.org.
http://github.com/patmorin/ods
 * @author morin
 * 
 * @param <T>
 * @see SortedSSet<T>
 */
public interface SSet<T> {
	
	/**
	 * @return the number of elements in this SSet
	 */
	public int size();

	/**
	 * Find the smallest element in the SSet that is greater than or equal to x.
	 * 
	 * @param x
	 * @return the smallest element in the SSet that is greater than or equal to
	 *         x or null if no such element exists
	 */
	public T find(T x);

	/**
	 * Add the element x to the SSet
	 * 
	 * @param x
	 * @return true if the element was added or false if x was already in the
	 *         set
	 */
	public boolean add(T x);

	/**
	 * Remove the element x from the SSet
	 * 
	 * @param x
	 * @return true if x was removed and false if x was not removed (because x
	 *         was not present)
	 */
	public boolean remove(T x);

	/**
	 * Clear the SSet, removing all elements from the set
	 */
	public void clear();

    /**
     * find the key of the node before the value x.
     * @param x
     * @return 
     */
    public T predecessor(T x);
    
    /**
     * find the key of the node after the value x.
     * @param x
     * @return 
     */
    public T successor(T x);
	
    /**
     * find the minimum key within the nodes.
     * @return 
     */
    public T minimum();
    
    /**
     * find the maximum key within the nodes.
     * @return 
     */
    public T maximum();
    
}
