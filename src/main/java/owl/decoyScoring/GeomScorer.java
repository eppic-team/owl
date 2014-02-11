package owl.decoyScoring;

import java.io.IOException;
import java.sql.SQLException;

import javax.vecmath.Point3d;

import edu.uci.ics.jung.graph.util.Pair;
import owl.core.structure.PdbChain;
import owl.core.structure.AaResidue;
import owl.core.structure.graphs.RIGEdge;
import owl.core.structure.graphs.RIGGeometry;
import owl.core.structure.graphs.RIGNode;
import owl.core.structure.graphs.RIGraph;
import owl.gmbp.CMPdb_sphoxel;

public class GeomScorer extends Scorer {
	
	private static final double DEFAULT_CUTOFF = 8.0;
	private static final int DEFAULT_MIN_SEQ_SEP = 3;
	private static final String DEFAULT_CT = "Cb";
	public static final int FCT_EC = 0;
	public static final int FCT_SUM_LOS = 1;
	public static final int FCT_EC_LOSg0 = 2;
	public static final int FCT_EC_LOSs0 = 3;
	public static final int FCT_SUM_LOSnorth = 4;
	public static final int FCT_SUM_LOSsouth = 5;
	private static final int INVALID_LOS = -666666;
	
//	protected RIGraph graph;
//	protected RIGGeometry graphGeom; 
	protected String archiveFN=null;
//	private boolean scoreGeom = false;
	private int scoringFctType;
	

	private boolean normaliseScore = false;
	private int directNbThres = 1;

	public GeomScorer(String archiveFN, String ct, double cutoff, int minSeqSep){
		this.archiveFN = archiveFN;
		this.ct = ct;
		this.cutoff = cutoff;
		this.minSeqSep = minSeqSep;
	}
	
	public GeomScorer(String ct, double cutoff, int minSeqSep){
		this.archiveFN = null;
		this.ct = ct;
		this.cutoff = cutoff;
		this.minSeqSep = minSeqSep;
	}
	
	public GeomScorer(String archiveFN){
		this.archiveFN = archiveFN;
		this.ct = DEFAULT_CT;
		this.cutoff = DEFAULT_CUTOFF;
		this.minSeqSep = DEFAULT_MIN_SEQ_SEP;
	}
	
	public double scoreIt(PdbChain pdb) {
		RIGraph graph = pdb.getRIGraph(this.ct, this.cutoff);
		graph.restrictContactsToMinRange(minSeqSep);
//		switch (scoringFctType){
//		case FCT_EC:
//		case FCT_SUM_LOS:
//		case FCT_EC_LOSg0:
//		case FCT_EC_LOSs0:
//		case FCT_SUM_LOSnorth:
//		case FCT_SUM_LOSsouth:			
//		}
		
		if (scoringFctType == FCT_EC){
			return scoreIt(pdb, graph);
		}
		else {
			RIGGeometry graphGeom = new RIGGeometry(graph, pdb);
			return scoreIt(pdb, graph, graphGeom);
		}
//		if (scoreGeom){
//			RIGGeometry graphGeom = new RIGGeometry(graph, pdb.getResidues());
//			return scoreIt(pdb, graph, graphGeom);			
//		}
//		else{
//			return scoreIt(pdb, graph);
//		}
	}
	
	public double scoreIt(PdbChain pdb, RIGraph graph){
		int numEdges=0;
				
		for (RIGEdge edge:graph.getEdges()){
			Pair<RIGNode> nodes = graph.getEndpoints(edge);
			RIGNode iNode = nodes.getFirst();
			RIGNode jNode = nodes.getSecond();
			
			int iNum = iNode.getResidueSerial();
			int jNum = jNode.getResidueSerial();
			
			if (Math.abs(iNum-jNum)>this.directNbThres){
				numEdges++;
			}
		}
		
		return numEdges;
	}
	
	public double scoreIt(PdbChain pdb, RIGraph graph, RIGGeometry graphGeom){
		double score=0;
		int numSummedScores=0;
		CMPdb_sphoxel sphoxel = null;
		if (archiveFN==null){
			try {
				sphoxel = new CMPdb_sphoxel();
			} catch (SQLException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
		}
		
		for (RIGEdge edge:graph.getEdges()){
			Pair<RIGNode> nodes = graph.getEndpoints(edge);
			RIGNode iNode = nodes.getFirst();
			RIGNode jNode = nodes.getSecond();
			
			int iNum = iNode.getResidueSerial();
			int jNum = jNode.getResidueSerial();
			
			AaResidue atomsI = (AaResidue)pdb.getResidue(iNum);
			AaResidue atomsJ = (AaResidue)pdb.getResidue(jNum);
			
			if (atomsI.containsAtom("CA") && atomsJ.containsAtom("CA")){
				Point3d iCoord=atomsI.getAtom("CA").getCoords();
				Point3d jCoord=atomsJ.getAtom("CA").getCoords();
				
				double resDist=Math.sqrt((Math.pow((iCoord.x-jCoord.x), 2)+Math.pow((iCoord.y-jCoord.y), 2)+Math.pow((iCoord.z-jCoord.z), 2)));
				
				if (Math.abs(iNum-jNum)>this.directNbThres){
					boolean inOnNorthernHem = graphGeom.isContactOnNorthernHemisphere((AaResidue)pdb.getResidue(iNum), (AaResidue)pdb.getResidue(jNum));
					if (archiveFN!=null){
						try {
							// position of j with respect to i
							double val = graphGeom.getLogOddsScore((AaResidue)pdb.getResidue(iNum), (AaResidue)pdb.getResidue(jNum), resDist, archiveFN);
							val = testLOSforValidity(val, inOnNorthernHem);
							if (val!=INVALID_LOS){
								score += val;
								numSummedScores++;
							}
							// position of i with respect to j
							val = graphGeom.getLogOddsScore((AaResidue)pdb.getResidue(jNum),(AaResidue)pdb.getResidue(iNum), resDist, archiveFN);
							val = testLOSforValidity(val, inOnNorthernHem);
							if (val!=INVALID_LOS){
								score += val;
								numSummedScores++;
							}
						} catch (NumberFormatException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						} catch (IOException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
					} else{
						try {
							double val = sphoxel.getLogOddsScore((AaResidue)pdb.getResidue(iNum), (AaResidue)pdb.getResidue(jNum), resDist, graphGeom);
							val = testLOSforValidity(val, inOnNorthernHem);
							if (val!=INVALID_LOS){
								score += val;
								numSummedScores++;
							}
							val = sphoxel.getLogOddsScore((AaResidue)pdb.getResidue(jNum), (AaResidue)pdb.getResidue(iNum), resDist, graphGeom);
							val = testLOSforValidity(val, inOnNorthernHem);
							if (val!=INVALID_LOS){
								score += val;
								numSummedScores++;
							}
						} catch (SQLException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
					}				
				}
//				System.out.print(iNum+"  "+jNum+"  "+"Intermediate logOddsScores are : "+score+"\n");				
			}			
		}
//		System.out.println("Final logOddsScore for protein "+pdb.getPdbCode()+" is : "+score+" norm:"+(score/i)+"  for "+numSummedScores+" chosen contacts");
		
		if (this.normaliseScore)
			score = score/numSummedScores; // normalise score
		
		if (scoringFctType==FCT_EC || scoringFctType==FCT_EC_LOSg0 || scoringFctType==FCT_EC_LOSs0)
			return numSummedScores;
		else
			return score;
	}
	
	private double testLOSforValidity(double val, boolean contactIsNorth){
		switch (scoringFctType){
		case FCT_EC:
		case FCT_SUM_LOS:
		case FCT_EC_LOSg0:
			if (val<0)
				val = INVALID_LOS;
			break;
		case FCT_EC_LOSs0:
			if (val>0)
				val = INVALID_LOS;
			break;
		case FCT_SUM_LOSnorth:
			if (!contactIsNorth)
				val = INVALID_LOS;
			break;
		case FCT_SUM_LOSsouth:
			if (contactIsNorth)
				val = INVALID_LOS;
			break;
		}
		
		return val;
	}

	public void setNormaliseScore(boolean normaliseScore) {
		this.normaliseScore = normaliseScore;
	}

	public boolean isNormaliseScore() {
		return normaliseScore;
	}

	public void setDirectNbThres(int directNbThres) {
		this.directNbThres = directNbThres;
	}

	public int getDirectNbThres() {
		return directNbThres;
	}

	public int getScoringFctType() {
		return scoringFctType;
	}

	public void setScoringFctType(int scoringFctType) {
		this.scoringFctType = scoringFctType;
	}

}
