package tools;

import java.io.Serializable;
import java.util.Collection;
import java.util.Iterator;


/**
 * An implementation of <code>Collection</code> that stores exactly
 * 3 non-null objects and is not mutable.  The 3 values are treated as 
 * indistinguishable, i.e. <code>equals</code> and <code>hashCode</code> don't distinguish 
 * the order of the 3 values as long as they are the same.
 * They respect <code>equals</code> and may be used as indices or map keys.<p>
 * Note that they do not protect from malevolent behavior: if one 
 * object in the tuple is mutable, then it can be changed with the usual bad
 * effects.
 * The implementation is a copy of edu.uci.ics.jung.graph.util.Pair except for
 * the (important) fact that the order of the values here doesn't matter (in 
 * Pair the order does matter)
 * @author duarte
 */
public final class Triplet<T> implements Collection<T>, Serializable {

	private static final long serialVersionUID = -3293312441342412805L;
	private T first;
	private T second;
	private T third;
	private int hash; // the cached hash value

	/**
	 * Creates a <code>Triplet</code> from the specified elements.
	 * @param value1 the first value in the new <code>Triplet</code>
	 * @param value2 the second value in the new <code>Triplet</code>
	 * @param value3 the third value in the new <code>Triplet</code>
	 * @throws IllegalArgumentException if an argument is null
	 */
	public Triplet(T value1, T value2, T value3) {
		if(value1 == null || value2 == null || value3 == null) 
			throw new IllegalArgumentException("Triplet cannot contain null values");
		first = value1;
		second = value2;
		third = value3;
	}

	/**
	 * Creates a Triplet from the passed Collection.
	 * The size of the Collection must be 3.
	 * @param values the elements of the new <code>Triplet</code>
	 * @throws IllegalArgumentException if the input collection is null,
	 * contains null values, or has != 3 elements.
	 */
	public Triplet(Collection<? extends T> values) 
	{
		if (values == null)
			throw new IllegalArgumentException("Input collection cannot be null");
		if (values.size() == 3)
		{
			if(values.contains(null)) 
				throw new IllegalArgumentException("Triplet cannot contain null values");
			Iterator<? extends T> iter = values.iterator();
			first = iter.next();
			second = iter.next();
			third = iter.next();
		}
		else
			throw new IllegalArgumentException("Triplet may only be created from a Collection of exactly 3 elements");

	}

	/**
	 * Creates a <code>Triplet</code> from the passed array.
	 * The size of the array must be 3.
	 * @throws IllegalArgumentException if the input array is null,
	 * contains null values, or has != 3 elements.
	 */
	public Triplet(T[] values)
	{
		if (values == null)
			throw new IllegalArgumentException("Input array cannot be null");
		if (values.length == 3)
		{
			if(values[0] == null || values[1] == null || values[2] == null) 
				throw new IllegalArgumentException("Triplet cannot contain null values");
			first = values[0];
			second = values[1];
			third = values[2];
		}
		else
			throw new IllegalArgumentException("Triplet may only be created from an " +
			"array of 3 elements");
	}

	/**
	 * Returns the first element.
	 */
	public T getFirst() 
	{
		return first;
	}

	/**
	 * Returns the second element.
	 */
	public T getSecond() 
	{
		return second;
	}
	
	/**
	 * Returns the third element.
	 * @return
	 */
	public T getThird() {
		return third;
	}

	@SuppressWarnings("unchecked")
	@Override
	public boolean equals( Object o ) {
		if (o == this)
			return true;

		if (! (o instanceof Triplet)) {
			return false;
		}
		Triplet otherTriplet = (Triplet) o;
		Object otherFirst = otherTriplet.getFirst();
		Object otherSecond = otherTriplet.getSecond();
		Object otherThird = otherTriplet.getThird();

		if (
		(this.first  == otherFirst  || 
				(this.first != null  && this.first.equals(otherFirst)))   
		&&
		(this.second == otherSecond || 
				(this.second != null && this.second.equals(otherSecond)))
		&&
		(this.third == otherThird ||
				(this.third != null && this.third.equals(otherThird)))
		) return true;
		
		if (
		(this.first  == otherSecond  || 
				(this.first != null  && this.first.equals(otherSecond)))   
		&&
		(this.second == otherFirst || 
				(this.second != null && this.second.equals(otherFirst)))
		&&
		(this.third == otherThird ||
				(this.third != null && this.third.equals(otherThird)))
		) return true;
		
		if (
		(this.first  == otherThird  || 
				(this.first != null  && this.first.equals(otherThird)))   
		&&
		(this.second == otherSecond || 
				(this.second != null && this.second.equals(otherSecond)))
		&&
		(this.third == otherFirst ||
				(this.third != null && this.third.equals(otherFirst)))
		) return true;

		if (
		(this.first  == otherFirst  || 
				(this.first != null  && this.first.equals(otherFirst)))   
		&&
		(this.second == otherThird || 
				(this.second != null && this.second.equals(otherThird)))
		&&
		(this.third == otherSecond ||
				(this.third != null && this.third.equals(otherSecond)))
		) return true;

		return false;
	}

	@Override
	public int hashCode() {
		if (hash==0) { // we use the cached value if already calculated
			hash = (first==null ? 0 : first.hashCode())
			+ (second==null ? 0 : second.hashCode())
			+ (third==null ? 0 : third.hashCode());
		}
		return hash;
	}

	@Override
	public String toString()
	{
		return "<" + first.toString() + ", " + second.toString() + ", " + third.toString() + ">";
	}

	public boolean add(T o) {
		throw new UnsupportedOperationException("Triplets cannot be mutated");
	}

	public boolean addAll(Collection<? extends T> c) {
		throw new UnsupportedOperationException("Triplets cannot be mutated");
	}

	public void clear() {
		throw new UnsupportedOperationException("Triplets cannot be mutated");
	}

	public boolean contains(Object o) {
		return (first == o || first.equals(o) || second == o || second.equals(o) || third == o || third.equals(o));
	}

	public boolean containsAll(Collection<?> c) {
		if (c.size() > 3)
			return false;
		Iterator<?> iter = c.iterator();
		Object c_first = iter.next();
		Object c_second = iter.next();
		Object c_third = iter.next();
		return this.contains(c_first) && this.contains(c_second) && this.contains(c_third);
	}

	public boolean isEmpty() {
		return false;
	}

	public Iterator<T> iterator() {
		return new TripletIterator();
	}

	public boolean remove(Object o) {
		throw new UnsupportedOperationException("Triplets cannot be mutated");
	}

	public boolean removeAll(Collection<?> c) {
		throw new UnsupportedOperationException("Triplets cannot be mutated");
	}

	public boolean retainAll(Collection<?> c) {
		throw new UnsupportedOperationException("Triplets cannot be mutated");
	}

	public int size() {
		return 3;
	}

	public Object[] toArray() {
		Object[] to_return = new Object[3];
		to_return[0] = first;
		to_return[1] = second;
		to_return[2] = third;
		return to_return;
	}

	@SuppressWarnings("unchecked")
	public <S> S[] toArray(S[] a) {
		S[] to_return = a;
		Class<?> type = a.getClass().getComponentType();
		if (a.length < 3)
			to_return = (S[])java.lang.reflect.Array.newInstance(type, 3);
		to_return[0] = (S)first;
		to_return[1] = (S)second;
		to_return[2] = (S)third;

		if (to_return.length > 3)
			to_return[3] = null;
		return to_return;
	}

	private class TripletIterator implements Iterator<T>
	{
		int position;

		private TripletIterator()
		{
			position = 0;
		}

		public boolean hasNext()
		{
			return position < 3;
		}

		public T next()
		{
			position++;
			if (position == 1)
				return first;
			else if (position == 2)
				return second;
			else if (position == 3)
				return third;
			else
				return null;
		}

		public void remove()
		{
			throw new UnsupportedOperationException("Triplets cannot be mutated");
		}
	}
}

