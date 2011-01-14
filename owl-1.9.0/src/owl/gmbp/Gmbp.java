package owl.gmbp;

import java.util.Vector;

import owl.core.runners.tinker.TinkerConstraint;

import edu.uci.ics.jung.graph.util.Pair;

public class Gmbp {
	
	private Vector<Pair<Double>> phiRanges;			// permanent list of currently selected phi ranges
	private Vector<Pair<Double>> thetaRanges;       // permanent list of currently selected theta ranges
	private Vector<Pair<Integer>> selContacts;         // permanent list of currently selected and referred contacts
	
	public boolean hasConstraints() {
		if (phiRanges == null || thetaRanges == null) {
			return false;
		}
		return this.phiRanges.size() > 0 && this.thetaRanges.size() > 0;
	}
	
	public void setLambdaRanges(Vector<Pair<Double>> phiRanges) {
		this.phiRanges = phiRanges;
	}
	public Vector<Pair<Double>> getLambdaRanges() {
		return phiRanges;
	}
	public void setPhiRanges(Vector<Pair<Double>> thetaRanges) {
		this.thetaRanges = thetaRanges;
	}
	public Vector<Pair<Double>> getPhiRanges() {
		return thetaRanges;
	}
	public void setSelContacts(Vector<Pair<Integer>> selContacts) {
		this.selContacts = selContacts;
	}
	public Vector<Pair<Integer>> getSelContacts() {
		return selContacts;
	}

	public TinkerConstraint[] getConstraints() {
		TinkerConstraint[] constraints = new TinkerConstraint[selContacts.size()*2];
		for (int i = 0; i < selContacts.size(); i++) {
			constraints[i*2] = new TinkerConstraint(selContacts.elementAt(i).getFirst(), selContacts.elementAt(i).getSecond() ,thetaRanges.elementAt(i).getFirst(),thetaRanges.elementAt(i).getSecond(),
					500.0, TinkerConstraint.CONSTRAINT.GMBPTHETA);
			constraints[i*2+1] = new TinkerConstraint(selContacts.elementAt(i).getFirst(), selContacts.elementAt(i).getSecond() ,phiRanges.elementAt(i).getFirst(),phiRanges.elementAt(i).getSecond(),
					500.0, TinkerConstraint.CONSTRAINT.GMBPPHI);
		}
		return constraints;
		
	}
	
	

}