package owl.core.connections;

import owl.core.structure.features.SecondaryStructure;

/**
 * Abstract class which retrieves progress updates from JPredConnection and forwards them by means of
 * its setStatus, setError and setResult methods.
 * @author stehr
 *
 */
public abstract class JPredProgressRetriever {

	public static enum Status {
		
		START("Idle"), 
		SENDING("Connecting to server"), 
		JOBID("Job ID retrieved"), 
		CHECKING("Checking job status"), 
		WAITING("Job waiting"), 
		RUNNING("Job running"), 
		RESULT("Job finished"), 
		DONE("Results retrieved");
		
		String msg;
		
		Status(String msg) {
			this.msg = msg;
		}
		
		public String toString() {
			return msg;
		}
	};
	
	protected Object parent;	// parent which will respond to status changes or errors in the way specified
								// by the implemented methods setStatus and setError
	
	public JPredProgressRetriever(Object parent) {
		this.parent = parent;
	}
	
	public abstract void setStatus(Status stat);
	
	public abstract void setError(Throwable e);
	
	public abstract void setResult(SecondaryStructure ss);
	
}
