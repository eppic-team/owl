package graphAveraging;

import java.util.ArrayList;
import java.util.Collections;

import proteinstructure.Interval;

/**
 * A class to represent a consensus phi/psi angle interval.
 * Keeps information of the interval, the voteCount, voters and angle values.
 * Equality is defined (equals method) by equal voteCounts, voters and values 
 */
public class ConsensusInterval extends Interval {
	
	private int voteCount;
	private ArrayList<String> voters;
	private ArrayList<Integer> votersIndices;
	private ArrayList<Double> values; // the angle values within this interval
	
	public ConsensusInterval () {
		super(0,0);
		voters = new ArrayList<String>();
		values = new ArrayList<Double>();
		votersIndices = new ArrayList<Integer>();
		voteCount = 0;
	}
	
	public ConsensusInterval (int beg, int end) {
		super(beg,end);
		voters = new ArrayList<String>();
		values = new ArrayList<Double>();
		votersIndices = new ArrayList<Integer>();
		voteCount = 0;
	}
	
	public ConsensusInterval (int beg, int end, ArrayList<String> voters, ArrayList<Double> values, ArrayList<Integer> votersIndices, int voteCount) {
		super(beg, end);
		this.voters = voters;
		this.values = values;
		this.votersIndices = votersIndices;
		this.voteCount = voteCount;
	}		
	
	public void addVoter(int voterIndex, String voter, double value) {
		voteCount++;
		voters.add(voter);
		values.add(value);
		votersIndices.add(voterIndex);
	}
	
	public int getVoteCount() {
		return voteCount;
	}
	
	public ArrayList<String> getVoters() {
		return voters;
	}
	
	public ArrayList<Double> getValues() {
		return values;
	}
	
	public ArrayList<Integer> getVotersIndices() {
		return votersIndices;
	}
	
	public boolean contains(String voter) {
		return this.voters.contains(voter);
	}
	
	/**
	 * Returns true if this ConsensusInterval has same voteCount, same voters and same values than other
	 */
	public boolean equals(Object o) {
		if (! (o instanceof ConsensusInterval)) return false;
		ConsensusInterval other = (ConsensusInterval) o;
		if (other.voteCount!=this.voteCount)
			return false;
		for (String voter: voters) {
			if (!other.voters.contains(voter))
				return false;
		}
		for (double value:values) {
			if (!other.values.contains(value)) 
				return false;
		}
		return true;
	}
	
	/**
	 * Re-centers the interval of this IntervalCandidate to (maxValue-minValue)/2
	 * @param angleInterval
	 */
	public void reCenterInterval(int angleInterval) {
		
		// we know that all of our values are within the angleInterval.
		// now to solve the wrapping problem we must be sure that we take them in the right direction, i.e.
		// if it's a normal interval (non-wrapping) then max-min will be smaller than angleInterval. But if the interval wraps around
		// then max-min is bigger than the interval (because it's the complementary angle) 
		double max = Collections.max(values);
		double min = Collections.min(values);
		double intervWidth = angleDistance(max, min);

		if ((max-min)<=angleInterval) { // normal non-wrapping case
			// A 'graphical' view of this (x is center, '' the min/max(b/e) values and || the limits of our interval (B/E)): 
			//   |   '  x  '   |
			//   B   b     e   E
			//   <------w------>
			//       <--d-->
			// so we want to get B, E for our final interval. Being d=interWidth and w=angleInterval then: B=b+d/2-w/2, E=e-d/2+w/2
			int B = (int) (min+intervWidth/2-angleInterval/2);
			int E = (int) (max-intervWidth/2+angleInterval/2);
			this.beg = B;
			this.end = E;
		} else { 					// wrapping interval
			// we have to swap max/min because beginning and end are swapped
			int B = (int) (max+intervWidth/2-angleInterval/2);
			int E = (int) (min-intervWidth/2+angleInterval/2);
			this.beg = B;
			this.end = E;			
		}
	}
	
	protected static double wrapAngle(double angle) {
		return angle+360.0;
	}
	
	protected static double unwrapAngle(double angle){
		if (angle>180.0) return angle-360.0;
		else return angle;
	} 
	
	/**
	 * Converts any angle (positive, negative, > or < than 360) to an angle in the first cycle (0 to 360)
	 * @param angle
	 * @return
	 */
	private static double toFirst(double angle) {
		return (angle%360+360)%360;
	}

	/**
	 * From the 2 possible angle distances between 2 given angles (one <180, the other >180)
	 * returns the one below 180. Input angles can be in any order and can be positive/negative or in any 
	 * cycle (i.e. above 360 or below -360)
	 * @param angle1
	 * @param angle2
	 * @return
	 */
	protected static double angleDistance(double angle1, double angle2) {	
		return Math.min(toFirst(angle2-angle1),toFirst(angle1-angle2));
	}
}