package owl.core.util.actionTools;

/**
 * Basic action class. Provides the storage and retrievel of an object.
 * @author Lars Petzold
 */

public class Action {
    Object obj;

    public Action(Object obj) {
	this.obj = obj;
    }

    public synchronized Object getObject() {
	return obj;
    }

    public synchronized void setObject(Object obj) {
	this.obj = obj; 
    }
}
