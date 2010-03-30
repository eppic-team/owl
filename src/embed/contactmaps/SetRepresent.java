package embed.contactmaps;

import java.util.*;
//import Analysis.Numerical;
public class SetRepresent<T> {
	
	private HashSet<SetElement<T>> element_set;
	
	public SetRepresent (){
		setFields();
	}
	
	public void setFields (){
		element_set = new HashSet<SetElement<T>> ();
	}
	public boolean add (T element){
		SetElement<T> el = new SetElement<T>(element);
		if(contains(element)) return false;
		else {element_set.add(el); return true;}
	}
	public boolean addAll (T[] array){
		int length = array.length;
		if(length > 0){
			boolean tester = false;
			for(int i = 0; i < length; i++){
				if(add(array[i])) tester = true;
			}
			return tester;
		}
		else return false;
	}
	public boolean addAll (Collection<T> coll){
		if(coll != null && coll.size() > 0){
			Iterator<T> it = coll.iterator();
			boolean tester = false;
			while(it.hasNext()){
				T element = it.next();
				if(add(element)) tester = true;
			}
			return tester;
		}
		else return false;
	}
	public boolean contains(T element){
		if(element_set.contains(new SetElement<T>(element))) return true;
		else return false;
	}
	public boolean containsAll(Collection<T> coll){
		if(coll != null && coll.size()>0){
			boolean tester = true;
			Iterator<T> it = coll.iterator();
			while(it.hasNext()){
				T element = it.next();
				if(contains(element) && tester);
				else {tester = false; break;}
			}
			return tester;
		}
		else return false;
	}
	public boolean remove (T element){
		if(contains(element)){
			return element_set.remove(element);
		}
		else return false;
	}
	public boolean removeAll (Collection<T> coll){
		if(containsAll(coll)) return element_set.removeAll(coll);
		else return false;
	}
	public Iterator<SetElement<T>> iterator(){
		if(element_set != null){
			return element_set.iterator();
		}
		else return null;
	}
	public String toString (){
		if(element_set != null){
			Iterator<SetElement<T>> it = iterator();
			StringBuffer strbf = new StringBuffer();
			while(it.hasNext()){
				SetElement<T> element = it.next();
				if(it.hasNext()) strbf.append(element).append("\n");
				else strbf.append(element);
			}
			if(strbf.length() > 0) return strbf.toString();
			else return "";
		}
		return "";
	}
	
	public static void main(String[] args){
		Numerical<?>[] ar = {new Numerical<Integer>(new Integer(1)),new Numerical<Double>(new Double(1.0)),new Numerical<ComplexDouble>(new ComplexDouble(1.0,0.0))};
		SetRepresent<Object> set = new SetRepresent<Object>();
		set.addAll(ar);
		System.out.println(set.toString());
	}
}