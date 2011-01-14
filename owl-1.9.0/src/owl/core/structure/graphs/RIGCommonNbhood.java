package owl.core.structure.graphs;

import java.util.TreeMap;
import java.util.Collections;

import owl.core.structure.AAinfo;

public class RIGCommonNbhood extends TreeMap<Integer,RIGNode> {
	

	private static final long serialVersionUID = 1L;

	private RIGNode iNode;
	private RIGNode jNode;
	private boolean connected;
	
	public RIGCommonNbhood(RIGNode iNode, RIGNode jNode, boolean connected){
		super();
		this.iNode = iNode;
		this.jNode = jNode;
		this.connected = connected;
	}

	public RIGNode getFirstNode() {
		return iNode;
	}
	
	public RIGNode getSecondNode() {
		return jNode;
	}
	
	public boolean areNodesConnected() {
		return connected;
	}
	
	public String getMotif(){
		String motif="";
		int min=Math.min(Math.min(iNode.getResidueSerial(),jNode.getResidueSerial()), Collections.min(this.keySet()));
		int max=Math.max(Math.max(iNode.getResidueSerial(),jNode.getResidueSerial()), Collections.max(this.keySet()));
		int gapSize=0;
		String gap="";
		for (int i=min;i<=max;i++){
			if (this.containsKey(i)){
				motif+=gap;
				motif+=AAinfo.threeletter2oneletter(this.get(i).getResidueType());
				gapSize=0;
				gap="";
			} else if (i==iNode.getResidueSerial()){
				motif+=gap;
				motif+="x";
				gapSize=0;
				gap="";
			} else if (i==jNode.getResidueSerial()){
				motif+=gap;
				motif+="y";
				gapSize=0;
				gap="";
			} else {
				gapSize++;
				gap="_{"+gapSize+"}";
			}
		}
		return motif;
	}

	public String getMotifFullGaps(){
		String motif="";
		int min=Math.min(Math.min(iNode.getResidueSerial(),jNode.getResidueSerial()), Collections.min(this.keySet()));
		int max=Math.max(Math.max(iNode.getResidueSerial(),jNode.getResidueSerial()), Collections.max(this.keySet()));
		for (int i=min;i<=max;i++){
			if (this.containsKey(i)){
				motif+=AAinfo.threeletter2oneletter(this.get(i).getResidueType());	
			} else if (i==iNode.getResidueSerial()){
				motif+="x";
			} else if (i==jNode.getResidueSerial()){
				motif+="y";				
			} else {
				motif+="_";
			}
		}
		return motif;
	}

	public String getMotifNoGaps(){
		String motif="";
		int min=Math.min(Math.min(iNode.getResidueSerial(),jNode.getResidueSerial()), Collections.min(this.keySet()));
		int max=Math.max(Math.max(iNode.getResidueSerial(),jNode.getResidueSerial()), Collections.max(this.keySet()));
		for (int i=min;i<=max;i++){
			if (this.containsKey(i)){
				motif+=AAinfo.threeletter2oneletter(this.get(i).getResidueType());	
			} else if (i==iNode.getResidueSerial()){
				motif+="x";
			} else if (i==jNode.getResidueSerial()){
				motif+="y";				
			} 
		}
		return motif;
	}
	
	public String getCommaSeparatedResSerials(){
		String ressers="";
		for (int resser:this.keySet()){
			ressers += resser+","; 
		}
		// we chop off the last comma
		if (ressers.length()!=0) { 
			ressers = ressers.substring(0, ressers.length()-1);
		}
		return ressers;
	}
	
	public String toString(){
		if (this.isEmpty()) return "";
		else return this.getMotif();
	}
}
