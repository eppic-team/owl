package owl.core.util.actionTools;

/**
 * Abstract class designated to check conditions connected to the set object.
 * 
 * @author Lars Petzold
 */
public abstract class Checker extends Action {
    public Checker(Object obj) {
	super(obj);
    }
    public abstract boolean check();
}
