package actionTools;

/**
 * Abstract class designated to perform any kind of task connected to the set 
 * object.
 * 
 * @author Lars Petzold
 *
 */
public abstract class Doer extends Action {
    public Doer(Object obj) {
	super(obj);
    }
    public abstract void doit();
}