package graphAveraging;

/**
 * 
 * A square in phi/psi angle space. Think of a square region in a Ramachandran plot.
 * We use 1st, 2nd dimension to refer to phi and psi dimensions 
 *
 */
public class ConsensusSquare {
 
	private ConsensusInterval consInterval1stDim; // the interval containing the phi angles
	private ConsensusInterval consInterval2ndDim; // the interval containing the psi angles

	public ConsensusSquare (ConsensusInterval consInterval1stDim, ConsensusInterval consInterval2ndDim) {
		this.consInterval1stDim = consInterval1stDim;
		this.consInterval2ndDim = consInterval2ndDim;
		checkInput();
	}
	
	private void checkInput() {
		if (consInterval1stDim.getVoteCount()!=consInterval2ndDim.getVoteCount()) {
			throw new IllegalArgumentException("Different voter count for the 2 ConsensusInterval passed to constructor of this ConsensusSquare");
		}
		for (String voter:consInterval1stDim.getVoters()) {
			if (!consInterval2ndDim.contains(voter)) {
				throw new IllegalArgumentException("Different voters found in the 2 ConsensusInterval passed to constructor of this ConsensusSquare");
			}
		}
	}
 
	public boolean equals(Object o) {
		if (! (o instanceof ConsensusSquare)) {
			return false;
		}	
		ConsensusSquare other = (ConsensusSquare) o;
		if (! other.getConsInterval1stDim().equals(this.consInterval1stDim)) {
			return false;
		}
		if (! other.getConsInterval2ndDim().equals(this.consInterval2ndDim)) {
			return false;
		}
		return true;
	}
	
	public int getVoteCount() {
		// we've already checked with checkInput than the 2 dimensions have the same vote count, so we can safely return the count from the first
		return this.consInterval1stDim.getVoteCount();
	}
	
	public ConsensusInterval getConsInterval1stDim () {
		return consInterval1stDim;
	}
	
	public ConsensusInterval getConsInterval2ndDim () {
		return consInterval2ndDim;
	}
	
}
