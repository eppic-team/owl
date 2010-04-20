package owl.gmbp;

import java.util.Vector;

import edu.uci.ics.jung.graph.util.Pair;

public class Gmbp {
	
	private Vector<Pair<Float>> phiRanges;			// permanent list of currently selected phi ranges
	private Vector<Pair<Float>> thetaRanges;       // permanent list of currently selected theta ranges
	private Vector<Pair<Integer>> selContacts;         // permanent list of currently selected and referred contacts
	
	
	public void setPhiRanges(Vector<Pair<Float>> phiRanges) {
		this.phiRanges = phiRanges;
	}
	public Vector<Pair<Float>> getPhiRanges() {
		return phiRanges;
	}
	public void setThetaRanges(Vector<Pair<Float>> thetaRanges) {
		this.thetaRanges = thetaRanges;
	}
	public Vector<Pair<Float>> getThetaRanges() {
		return thetaRanges;
	}
	public void setSelContacts(Vector<Pair<Integer>> selContacts) {
		this.selContacts = selContacts;
	}
	public Vector<Pair<Integer>> getSelContacts() {
		return selContacts;
	}
	
	

}