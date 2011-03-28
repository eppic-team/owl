package owl.core.util.actionTools;

/**
 * @author Lars Petzold
 * 
 * The purpose of this class is to provide an abstract strategy to retrieve
 * data from a source and store it in target object. Therefore, you have to 
 * implement function <code>retrieve(Object obj)</code> as it suits best your 
 * requirements. The data to be stored might violate the storable data type,
 * hence, this function throws <code>java.lang.ClassCastExceptions</code>.
 */

public abstract class Retriever extends Action {
    public Retriever(Object obj) {
	super(obj);
    }
    public abstract void retrieve(Object obj) throws ClassCastException;
}
