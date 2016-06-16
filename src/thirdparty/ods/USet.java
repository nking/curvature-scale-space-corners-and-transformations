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
public interface USet<T> extends Iterable<T> {
	public int size();
	public boolean add(T x);
	public T remove(T x);
	public T find(T x);
	public void clear();
}
