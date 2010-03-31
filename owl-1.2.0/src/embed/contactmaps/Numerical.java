package embed.contactmaps;

import java.util.*;

/**
 * Generalized wrapper class: wraps all kinds of numerical value classes, provided with an
 * internal coding, i.e.
 * <p><tt>0 == java.lang.Integer</tt></p>
 * <p><tt>1 == java.lang.Double</tt></p>
 * <p><tt>2 == Analysis.ComplexDouble</tt></p>
 * <b>Note, that if the object does not fit one of these classes</b>
 * a <code> {@link NumericalException}</code>
 * is thrown. The two types <tt>Integer </tt> and <tt>Double</tt> are the only upcastable
 * types, wrapped by this class. Any attempt to upcast a <tt>ComplexDouble</tt>
 * typed instance causes a {@link NumericalException}.
 * @author hendrik1
 */
public class Numerical<T> {
    /**
     * All canonical numerical classes, that will be wrapped by this class.
     */
    private static final String canonicals = "java.lang.Integer,java.lang.Double,embed.contactmaps.ComplexDouble";
    /** The standard threshold, defining the minimal distinction distance of two
     * double floating points. If <tt>Math.abs(value1 - value2)</tt> is
     * equal to or less than <tt>threshold = 1.0e-14</tt>, both <tt>value1</tt> and <tt>value2</tt>
     * are considered equalt.
     */
    private static final double threshold  = 1.0e-14;
    /** The type of the class*/
    private Class<T> type;
    /** The element*/
    private T element;
    /** Upcasting checker: only Integers and Doubles are considered 'upcastable'*/
    private boolean upcastable;
    /** Downcasting checker: only Double representing Integers, with a tolerance of <tt>1.0e-14</tt>
     * or ComplexDouble representing Doubles or Integers are downcastable.
     */
    private boolean downcastable;
    /** A set of Integers representing types, to which an instance can be downcasted.*/
    private HashSet<Integer> down_caster;
    private HashSet<Integer> up_caster;
    /** Type coding: 0 = Integer, 1 = Double and 2 = ComplexDouble*/
    private int type_code;
    /**
     * One parameter constructor: creates an instance of this class. Note, that if
     * the type of this object <tt>obj</tt> is neither <tt>Integer</tt>, <tt>Double</tt>
     * nor <tt>ComplexDouble</tt> a <tt>NumericalException</tt> is thrown.
     * @param obj the object representing a canonical numerical value
     * @throws Analysis.NumericalException if the type of the object is neither <tt>Integer</tt>,
     * <tt>Double</tt> nor <tt>ComplexDouble</tt>
     */
    public Numerical (T obj) throws NumericalException {
        setFields(obj);
    }
    /**
     * One parameter constructor: creates an instance of this class by copying <tt>num</tt>.
     * Since the parameter has already been checked for type match, no further type checking is
     * done.
     * @param num a numerical value representation
     */
    public Numerical (Numerical<T> num) {
        this((T) num.getElement());
    }
    /**
     * Initial setter: wraps the object <tt>obj</tt> into an instance of this
     * class, by first checking whether <tt>obj</tt> is of type
     * <p><tt>java.lang.Integer</tt>, <tt>java.lang.Double</tt> or <tt>Analysis.ComplexDouble</tt>.</p>
     * If the method <code>{@link Numerical#isNumerical() }</code> returns <tt>true</tt>,
     * the initialization is completed and the instance can be further processed.
     * If <tt>false</tt> is returned, a <tt>NumericalException</tt> thrown.
     * @param obj an object representing a numerical value
     * @throws Analysis.NumericalException if the type of the object is neither <tt>Integer</tt>,
     * <tt>Double</tt> nor <tt>ComplexDouble</tt>
     */
    @SuppressWarnings("unchecked")
	public void setFields (T obj) throws NumericalException {
        type = (Class<T>) obj.getClass();
        if(isNumerical()){
            copyElement(obj);
            isCastable();
        }
        else throw new NumericalException("Any parameter must be of type 'Integer'," +
                " 'Double' or 'ComplexDouble'!");
    }
    /**
     * Internal setter: copies the object <tt>obj</tt>.
     * @param obj an object representing a numerical value
     */
    private void copyElement(T obj){
        element = obj;
    }
    /**
     * Checks, whether this instance can be upcasted or downcasted to one of the canonical
     * numerical types. That is:
     * <p><tt>Integer</tt>, <tt>Double</tt> are considered 'upcastable'
     * <p><tt>ComplexDouble</tt> and <tt>Double</tt> are considered 'downcastable',
     * only if the residual <tt>resi := value - Math.floor(value)</tt> is less or equal
     * to <tt>1.0e-14</tt>, where <p><tt>value = ((ComplexDouble) getElement()).getRealPart()</tt> and
     * <p><tt>value = ((Double) getElement()).doubleValue()</tt>, respectively.
     * <p>Additionally, an instance of type <tt>ComplexDouble</tt>
     * is only downcastable, if the method call <tt>getElement().getImpart()</tt> returns
     * an <tt>double</tt> element, less than or equal to <tt>1.0e-14</tt>.
     */
    private void isCastable (){
        String[] types = canonicals.split(",");
        for(int i = 0; i < 3; i++){
            if(types[i].matches(getType().getName())) type_code = i;
        }
        down_caster = new HashSet<Integer> (4);
        up_caster   = new HashSet<Integer> (4);
        switch(type_code){
            case 0: {
                upcastable   = true;
                downcastable = false;
                up_caster.add(new Integer(1));up_caster.add(new Integer(2));
                break;
            }
            case 1: {
                upcastable   = true;
                double value = ((Double) getElement()).doubleValue();
                double resi  = value - Math.floor(value);
                up_caster.add(new Integer(2));
                if(resi<=threshold) {downcastable = true; down_caster.add(new Integer(0));}
                break;
            }
            case 2: {
                upcastable   = false;
                ComplexDouble value = (ComplexDouble) getElement();
                if(Math.abs(value.getImPart())<=threshold){
                    down_caster.add(new Integer (1));
                    if(Math.abs(value.getRealPart() - Math.floor(value.getRealPart()))<=threshold){
                        downcastable = true;
                        down_caster.add(new Integer(0));
                    }
                }
                break;
            }
        }
    }
    public boolean isInvertible(){
        if(type_code == 0) return false;
        else return true;
    }
    /*------------------------------------getters-------------------------------*/
    /**
     * Addition method: takes <tt>this</tt> object and <tt>another_num</tt> and returns
     * the sum of <tt>this + another_num</tt>. Note, that due to internal implementation
     * upcasting is not required to avoid operation connected errors/exceptions.
     * @param another_num an other numerical value, to add to <tt>this</tt>
     * @return the sum of <tt>this + another_num</tt>
     */
    @SuppressWarnings("unchecked")
	public Numerical<?> add(Numerical<T> another_num){
        int type1 = getClassType(), type2 = another_num.getClassType();
        Numerical<?> sum = null;
        if(type1==type2){
            switch(type1){
                case 0 : {
                    sum = new Numerical<Integer>(((Integer) getElement()) + ((Integer) another_num.getElement()));
                    break;
                }
                case 1 : {
                    sum = new Numerical<Double>(((Double) getElement()) + ((Double) another_num.getElement()));
                    break;
                }
                case 2 : {
                    sum = new Numerical<ComplexDouble>(((ComplexDouble) getElement()).add((ComplexDouble) another_num.getElement()));
                    break;
                }
            }
            return sum;
        }
        else{
            int type3 = 0;
            boolean latter = false;
            if(type1<type2) {type3 = type2; latter = true;}
            else type3 = type1;
            switch(type3){
                case 1 : {
                    Numerical<Double> n_num = null; Numerical<Double> m_num = null;
                    if(latter) {n_num = (Numerical<Double>)upcast((Class<T>) getClassType(1)); m_num = new Numerical<Double> ((Numerical<Double>) another_num);}
                    else {n_num = (Numerical<Double>)another_num.upcast((Class<T>)getClassType(1));m_num = new Numerical<Double>((Numerical<Double>)this);}
                    sum = new Numerical<Double> (((Double) n_num.getElement()) + ((Double) m_num.getElement()));
                    break;
                }
                case 2 : {
                    Numerical<ComplexDouble> n_num = null; Numerical<ComplexDouble> m_num = null;
                    if(latter) {n_num = (Numerical<ComplexDouble>)upcast((Class<T>)getClassType(2)); m_num = new Numerical<ComplexDouble> ((Numerical<ComplexDouble>)another_num);}
                    else {n_num = (Numerical<ComplexDouble>)another_num.upcast((Class<T>)getClassType(2));m_num = new Numerical<ComplexDouble>((Numerical<ComplexDouble>)this);}
                    sum = new Numerical<ComplexDouble> (((ComplexDouble) n_num.getElement()).add((ComplexDouble) m_num.getElement()));
                    break;
                }
            }
            return sum;
        }
    }
    @SuppressWarnings("unchecked")
	public boolean equals (Numerical<T> numerical){
        if(type_code == numerical.type_code){
            int ntype = numerical.type_code;
            boolean tester = false;
            switch (ntype){
                case 0 : {
                    Integer val1 = (Integer) getElement();
                    Integer val2 = (Integer) numerical.getElement();
                    if(val1.intValue() == val2.intValue()) tester = true;
                    else tester = false;
                    break;
                }
                case 1 : {
                    Double val1 = (Double) getElement();
                    Double val2 = (Double) numerical.getElement();
                    if(Math.abs(val1.doubleValue()-val2.doubleValue())<=threshold) tester = true;
                    else tester = false;
                    break;
                }
                case 2 : {
                    ComplexDouble val1 = (ComplexDouble) getElement();
                    ComplexDouble val2 = (ComplexDouble) numerical.getElement();
                    if(val1.equals(val2)) tester = true;
                    else tester = false;
                    break;
                }
            }
            return tester;
        }
        else{
            boolean tester = false;
            int type1 = type_code, type2 = numerical.type_code, upcastid = 0;
            int upcst = 0;
            Numerical<?> upnum = null; Numerical<?> kept = null;
            if(type1 < type2) upcastid = 1;
            else upcastid = 2;
            switch(upcastid){
                case 1 : {
                    upnum = upcast((Class<T>)getClassType(type2));
                    kept  = new Numerical<Double> ((Numerical<Double>)numerical);
                    upcst = type2;
                    break;
                }
                case 2 : {
                    upnum = numerical.upcast((Class<T>)getClassType(type1));
                    kept  = new Numerical<ComplexDouble> ((Numerical<ComplexDouble>)this);
                    upcst = type1;
                    break;
                }
            }
            switch(upcst){
                case 1 : {
                    Double val1 = (Double) upnum.getElement();
                    Double val2 = (Double) kept.getElement();
                    if(Math.abs(val1.doubleValue()-val2.doubleValue())<=threshold) tester = true;
                    else tester = false;
                    break;
                }
                case 2 : {
                        ComplexDouble val1 = (ComplexDouble) upnum.getElement();
                        ComplexDouble val2 = (ComplexDouble) kept.getElement();
                        if(val1.equals(val2)) tester = true;
                        else tester = false;
                        break;
                    }
            }
            return tester;
        }
    }
    /**
     * Checks, whether the instance of this class represents a numerical value type, i.e.
     * <p><tt>java.lang.Integer</tt>, <tt>java.lang.Double</tt> or <tt>Analysis.ComplexDouble</tt>.</p>
     * Returns true, if <tt>type.getName()</tt> returns a String instance of the
     * abovementioned types. More precisely:
     * <p><tt>Object obj = new Object(...); Class type = obj.getClass();</tt></p>
     * <p><tt>String str = "java.lang.Integer,java.lang.Double,Analysis.ComplexDouble"</tt>;//<i>the canonical numerical types</i></p>
     * <p><tt>str.contains(type.getName())</tt>//<i>returns true, if type equals one of these classes</i></p>
     * operates similarly to this method.
     * @return true, iff the type of this object is canonical numerical
     */
    public boolean isNumerical(){
        String type_nm = type.getName();
        if(canonicals.contains(type_nm)) return true;
        else return false;
    }
    /**
     * Returns an integer value, coding for the three canonical numerical classes,
     * wrapped by this class. The coding is as follows:
     * <p><tt>0 == Integer</tt>
     * <p><tt>1 == Double</tt>
     * <p><tt>2 == ComplexDouble</tt>
     * @return an int value representing the canonical numerical standard classes
     */
    public int getClassType (){
        return type_code;
    }
    /**
     * Returns the type code of <tt>clazz</tt>. That is:
     * <p><tt>0 == java.lang.Integer</tt>, <tt>0 == java.lang.Double</tt> and
     * <tt>2 == Analysis.ComplexDouble</tt>
     * Except for 2, all types are considered as 'upcastable'. If the type is
     * not compatible, a <tt>NumericalException</tt> is thrown.
     * @param clazz the type
     * @return an integer value between 0 and 2 (inclusively)
     * @throws NumericalException if the type of <tt>clazz</tt> is incompatible
     */
    public int getClassType (Class<T> clazz){
        String[] types = canonicals.split(",");
        int type_val = -1;
        for(int i = 0; i < 3; i++){
            if(types[i].matches(clazz.getName())) type_val = i;
        }
        if(0<=type_val && type_val<=2) return type_val;
        else throw new NumericalException("Incompatible type!");
    }
    /**
     * Returns the class type coded by <tt>type_code</tt> parameter. <b>Note, that</b>
     * the parameter must be integer value of <tt>0,1</tt> or <tt>2</tt>, where the
     * parameter value stands for the following types:
     * <p><tt>0 == Integer</tt>
     * <p><tt>1 == Double</tt> and
     * <p><tt>2 == ComplexDouble</tt>.
     * <p>If the parameter is beyond the defined range, a <tt>NumericalException</tt> is
     * thrown.
     * @param type_code the type code parameter
     * @return the class type of a standard numerical value representation
     * @throws Analysis.NumericalException if the parameter <tt>type_code</tt> is out of the
     * range <tt>Z and [0,2]</tt>
     */
    public Class<?> getClassType (int type_code) throws NumericalException {
        if(0 <= type_code && type_code <= 2){
            Class<?> clazz = null;
            switch(type_code){
                case 0 : {
                    clazz = (new Integer(0)).getClass();
                    break;
                }
                case 1 : {
                    clazz = (new Double (0.0)).getClass();
                    break;
                }
                case 2 : {
                    clazz = (new ComplexDouble(0.0,0.0)).getClass();
                    break;
                }
            }
            return clazz;
        }
        else throw new NumericalException ("Class type code only integer values " +
                "0, 1 or 2, but was set to "+type_code+"!");
    }
    /**
     * Returns a copy of this numerical value representing object. <b>Note, that
     * the return value is of type</b> <tt>Object</tt>, so any rereferencing to
     * the former type requires casting. More precisely:
     * <p><tt>Integer  value = new Integer (1);</tt></p>
     * <p><tt>Numerical  num = new Numerical (value)</tt></p>
     * <p><tt>Integer nvalue = (Integer) num.getElement();</tt>//<i>returns a new object of type </i><tt>Integer</tt></p>
     * @return the element of this numerical value representation
     */
    public T getElement (){
        T el = null; el = element;
        return el;
    }
    /**
     * Returns the type of this instance.
     * @return the type
     */
    public Class<T> getType (){
        return type;
    }
    /**
     * Method to perform downcast operation on this instance, to the specified type of <tt>clazz</tt>.
     * <b>Note</b>, that this method call fails in case a floating point representation is attempted to downcast,
     * but causes a loss of precision (like a <tt>Double val = new Double (1.000001)</tt> value
     * exceeds the tolerance of <tt>1.0e-14</tt> to the value <tt>Integer ival = new Integer(1)</tt>).
     * Additionally, only <tt>ComplexDouble</tt> with imaginary part equal to zero (or below 
     * tolerance) are considered downcastable. If an instance violates those restrictions, a
     * {@link NumericalException} is thrown.
     * @param clazz the class representing the type for the downcast
     * @return a new Numerical object of the same type as <tt>clazz</tt>
     */
    public Numerical<?> downcast (Class<T> clazz){
        if(down_caster.size() > 0 || downcastable){
            int former = getClassType(), latter = getClassType(clazz);
            if(down_caster.contains(new Integer (latter))){
                Numerical<?> num = null;
                switch(latter){
                    case 0 : {
                        if(former == 1){
                            double value = ((Double) getElement()).doubleValue();
                            num = new Numerical<Integer> (new Integer((int) value));
                        }
                        else{
                            double value = ((ComplexDouble) getElement()).getRealPart();
                            num = new Numerical<Integer> (new Integer((int) value));
                        }
                        break;
                    }
                    case 1 : {
                        double value = ((ComplexDouble) getElement()).getRealPart();
                        num = new Numerical<Double> (new Double(value));
                        break;
                    }
                }
                return num;
            }
            else throw new NumericalException ("Downcast operation impossible due loss of precision!");
        }
        else throw new NumericalException ("No such operation possible!");
    }
    /**
     * Multiplication method: takes <tt>this</tt> object and <tt>another_num</tt> and returns
     * the product of <tt>this * another_num</tt>. Note, that due to internal implementation
     * upcasting is not required to avoid operation connected errors/exceptions.
     * @param another_num an other numerical value, to add to <tt>this</tt>
     * @return the sum of <tt>this + another_num</tt>
     */
    @SuppressWarnings("unchecked")
	public Numerical<?> multiply(Numerical<T> another_num){
        int type1 = getClassType(), type2 = another_num.getClassType();
        Numerical<?> sum = null;
        if(type1==type2){
            switch(type1){
                case 0 : {
                    sum = new Numerical<Integer>(((Integer) getElement()) * ((Integer) another_num.getElement()));
                    break;
                }
                case 1 : {
                    sum = new Numerical<Double>(((Double) getElement()) * ((Double) another_num.getElement()));
                    break;
                }
                case 2 : {
                    sum = new Numerical<ComplexDouble>(((ComplexDouble) getElement()).multiply((ComplexDouble) another_num.getElement()));
                    break;
                }
            }
            return sum;
        }
        else{
            int type3 = 0;
            boolean latter = false;
            if(type1<type2) {type3 = type2; latter = true;}
            else type3 = type1;
            switch(type3){
                case 1 : {
                    Numerical<Double> n_num = null; Numerical<Double> m_num = null;
                    if(latter) {n_num = (Numerical<Double>)upcast((Class<T>)getClassType(1)); m_num = new Numerical<Double> ((Numerical<Double>)another_num);}
                    else {n_num = (Numerical<Double>) another_num.upcast((Class<T>)getClassType(1));m_num = new Numerical<Double>((Numerical<Double>)this);}
                    sum = new Numerical<Double> (((Double) n_num.getElement()) * ((Double) m_num.getElement()));
                    break;
                }
                case 2 : {
                	Numerical<ComplexDouble> n_num = null; Numerical<ComplexDouble> m_num = null;
                    if(latter) {n_num = (Numerical<ComplexDouble>)upcast((Class<T>)getClassType(2)); m_num = new Numerical<ComplexDouble> ((Numerical<ComplexDouble>)another_num);}
                    else {n_num = (Numerical<ComplexDouble>)another_num.upcast((Class<T>)getClassType(2));m_num = new Numerical<ComplexDouble>((Numerical<ComplexDouble>)this);}
                    sum = new Numerical<ComplexDouble> (((ComplexDouble) n_num.getElement()).multiply((ComplexDouble) m_num.getElement()));
                    break;
                }
            }
            return sum;
        }
    }
    /**
     * Returns a String representation of this object of class <tt>Numerical</tt>.
     * The String representation is done, by calling the element return method {@link
     * #getElement()} casting the return object to its type and call its canonical
     * <tt>toString()</tt> method.
     * @return a String represention the element value
     */
    @Override
    public String toString(){
        return ((T) getElement()).toString();
    }
   
    /**
     * Returns an upcasted version of this instance as a new <tt>Numerical</tt>
     * instance. <b>Note</b>, that the only upcastable types are of type
     * <p><tt>java.lang.Integer</tt> or <tt>java.lang.Double</tt>.</p>
     * If this instance is not of this type, an <tt>NumericalException</tt> is thrown.
     * @param clazz the type to which this instance is upcasted
     * @return a new upcasted version of this instance
     * @throws Analysis.NumericalException if the instance is not upcastable
     */
    public Numerical<?> upcast (Class<T> clazz) throws NumericalException {
        if(upcastable || up_caster.size()>0){
            int up_cast = getClassType(clazz), former = getClassType();
            Numerical<?> nvalue = null;
            if(up_caster.contains(new Integer(up_cast))){
                switch(up_cast){
                    case 1 : {
                        double value = (double)((Integer) getElement()).intValue();
                        nvalue = new Numerical<Double>(new Double(value));
                        break;
                    }
                    case 2 : {
                        if(former == 0){
                            double value = (double)((Integer) getElement()).intValue();
                            nvalue = new Numerical<ComplexDouble>(new ComplexDouble(value,0.0));
                        }
                        else{
                            double value = ((Double) getElement()).doubleValue();
                            nvalue = new Numerical<ComplexDouble>(new ComplexDouble(value,0.0));
                        }
                        break;
                    }
                }
                return nvalue;
            }
            else throw new NumericalException ("No upcast operation possible!");
        }
        else throw new NumericalException ("Element is not upcastable!");
    }
        
    @SuppressWarnings("unchecked")
	public static void main(String[] args){
        Integer val = new Integer (1);
        Numerical<Integer> num1 = new Numerical<Integer> (val);
        Numerical<ComplexDouble> num2 = new Numerical<ComplexDouble> (new ComplexDouble(1.0,-0.0));
        Integer value = (Integer) num1.getElement();
        //Double nvalue = (Double) num1.getElement();
        for(int i = 0; i < 3; i++){
            value = value+value;
            System.out.println(value.toString());
        }
        val = val+val;
        Numerical<ComplexDouble> num3 = (Numerical<ComplexDouble>) num1.upcast((Class)num2.getType());
        Numerical<ComplexDouble> num45 = (Numerical<ComplexDouble>)num2.add(num3);
        Numerical<ComplexDouble> num54 = (Numerical<ComplexDouble>)num45.multiply(num45);
        System.out.println("Initial value: "+num1.toString()
                +"\nTesting: "+val.toString()+"\nupcast: "+num3.toString()
                +"\nsum: "+num45.toString()+"\nproduct: "+num54.toString());
        try{
            Class<Integer> ctype = (Class<Integer>) (new Integer(0)).getClass();
            Numerical<Integer> num4 = (Numerical<Integer>)num2.downcast((Class)ctype);
            System.out.println(num4.toString());
        }
        catch(NumericalException e){
            e.printStackTrace();
            System.err.println("object cant be downcasted!");
        }
    }
}
/**
 * Thrown to indicate, that the attempt to initialize an instance of the class {@link Numerical}
 * on an object failed, due to incompatiblity of the types.
 * @author hendrik1
 */
class NumericalException extends IllegalArgumentException {
    /**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	/**One parameter constructor: an instance of this class is initialized with
     * a detailed <tt>message</tt>*/
        public NumericalException (String message){
            super(message);
        }
    }

