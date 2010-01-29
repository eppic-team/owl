package actionTools;


public abstract class TinkerStatusNotifier extends Action {

	public TinkerStatusNotifier(Object obj) {
		super(obj);
	}
	
	public abstract void sendStatus(tinker.TinkerRunner.STATE s);
	
	public abstract void filesDone(int i);
}
