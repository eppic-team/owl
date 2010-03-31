package embed.contactmaps;

public class SetElement<T> {
	
	private T element;
	
	private Class<T> type;
	
	public SetElement (){
		element = null; type = null;
	}
	
	public SetElement(T element){
		setFields(element);
	}
	
	public SetElement(SetElement<T> el){
		this((T) el.element);
	}
	
	@SuppressWarnings("unchecked")
	public void setFields (T element){
		this.element = element;
		type = (Class<T>) element.getClass();
	}
	
	public boolean equals (SetElement<T> el){
		if(type.getName().matches(el.type.getName())){
			T el1 = (T) getElement(); T el2 = (T) el.getElement();
			if(el1.equals(el2)) return true;
			else return false;
		}
		else return false;
	}
	
	public T getElement (){
		T nelement = null; nelement = element;
		return nelement;
	}
	
	public Class<T> getType (){
		return type;
	}
	
	public String toString(){
		return ("Type: "+getType().getName()+" & ref. code: "+getElement().toString());
	}

	public static void main(String[] args){
		Integer val = new Integer (1);
		Double val1 = new Double(1.0);
		SetElement<Integer> set1 = new SetElement<Integer>(val);
		SetElement<Double> set2 = new SetElement<Double>(val1);
		if(set1.equals(set2)) System.out.println("Both instances are equal!");
		else System.out.println("Both instances are not equal!");
	}
}
