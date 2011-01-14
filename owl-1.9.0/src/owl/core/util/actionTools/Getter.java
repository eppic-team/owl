package owl.core.util.actionTools;

/**
 * Performs action and returns a value. The object being set is supposed to be 
 * connected to the action but not necessarily.
 * 
 * @author Lars Petzold
 */
public abstract class Getter extends Action {
    public Getter(Object obj) {
    	super(obj);
    }
    
    public abstract Object get() throws GetterError;
}
