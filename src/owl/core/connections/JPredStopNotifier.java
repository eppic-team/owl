package owl.core.connections;

/**
 * Holds a flag which a JPredConnection can query to stop execution.
 * Should be registered using JPredConnection.setStopNotifier. Then
 * another thread can directly set the flag pleaseStop to true which
 * will cause the JPredConnection to stop at the next possible instance.
 * @author stehr
 */
public class JPredStopNotifier {
	
	public volatile boolean pleaseStop;
	
	/**
	 * Returns the value of the stop flag.
	 * @return the value of the stop flag.
	 */
	public boolean getStop() {
		return this.pleaseStop;
	}
	
	/**
	 * Sets the value of the stop flag.
	 * @param stop the new value
	 */
	public void setStop(boolean stop) {
		this.pleaseStop = stop;
	}
}
