package actionTools;

/**
 * @author Lars Petzold
 * 
 * Abstract class to implement "inline runnables" for the thread executer by 
 * implementing function {@link #implRun()}.
 */
public abstract class Runner implements Runnable {
	
	/**
	 * Function to be implemented. This function is invoked in function 
	 * {@link #run()} and, therefore, is designated to provide everything to be 
	 * done in the thread.
	 */
	public abstract void implRun();
	
	/**
	 * The actual runner. Invokes function {@link #implRun()}.
	 */
	public void run() {
		implRun();
	}

}
