package owl.core.util;

public class GEExitStatus {

	/**
	 * failed: GE's internal indicator that a job runs or fails to run (for instance before submission)
	 */
	private int failed;
	
	/**
	 * exit_status: the actual exit status of the job itself
	 */
	private int exitStatus;
	
	public GEExitStatus(int failed, int exitStatus) {
		this.failed = failed;
		this.exitStatus = exitStatus;
	}
	
	/**
	 * Gets the failed GE value.
	 * This is GE's internal indicator that a job runs or fails to run (for instance before submission).
	 * It's 0 if successful, >0 otherwise
	 */
	public int getFailed() {
		return failed;
	}
	
	/**
	 * Gets the exit status of the job executed by GE.
	 * @return
	 */
	public int getExitStatus() {
		return exitStatus;
	}

}
