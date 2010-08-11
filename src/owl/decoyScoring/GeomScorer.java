package owl.decoyScoring;

import java.io.IOException;
import java.sql.SQLException;

import javax.vecmath.Point3d;

import edu.uci.ics.jung.graph.util.Pair;
import owl.core.structure.Pdb;
import owl.core.structure.graphs.RIGEdge;
import owl.core.structure.graphs.RIGGeometry;
import owl.core.structure.graphs.RIGNode;
import owl.core.structure.graphs.RIGraph;
import owl.gmbp.CMPdb_sphoxel;

public class GeomScorer extends Scorer {
	
	private static final double DEFAULT_CUTOFF = 8.0;
	private static final int DEFAULT_MIN_SEQ_SEP = 3;
	private static final String DEFAULT_CT = "Cb";
	
//	protected RIGraph graph;
//	protected RIGGeometry graphGeom; 
	protected String archiveFN=null;
	
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
	
	public double scoreIt(Pdb pdb) {
		RIGraph graph = pdb.getRIGraph(this.ct, this.cutoff);
		graph.restrictContactsToMinRange(minSeqSep);
		RIGGeometry graphGeom = new RIGGeometry(graph, pdb.getResidues());
		return scoreIt(pdb, graph, graphGeom);
	}
	
	public double scoreIt(Pdb pdb, RIGraph graph, RIGGeometry graphGeom){
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
			
			Point3d iCoord=pdb.getResidue(iNum).getAtomsMap().get("CA").getCoords();
			Point3d jCoord=pdb.getResidue(jNum).getAtomsMap().get("CA").getCoords();
			
			double resDist=Math.sqrt((Math.pow((iCoord.x-jCoord.x), 2)+Math.pow((iCoord.y-jCoord.y), 2)+Math.pow((iCoord.z-jCoord.z), 2)));
			
			if (Math.abs(iNum-jNum)>this.directNbThres){
				if (archiveFN!=null){
					try {
						// position of j with respect to i
						double val = graphGeom.getLogOddsScore(pdb.getResidue(iNum), pdb.getResidue(jNum), resDist, archiveFN);
						if (val!=-666666){
							score += val;
							numSummedScores++;
						}
						// position of i with respect to j
						val = graphGeom.getLogOddsScore(pdb.getResidue(jNum), pdb.getResidue(iNum), resDist, archiveFN);
						if (val!=-666666){
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
						double val = sphoxel.getLogOddsScore(pdb.getResidue(iNum), pdb.getResidue(jNum), resDist, graphGeom);
						if (val!=-666666){
							score += val;
							numSummedScores++;
						}
						val = sphoxel.getLogOddsScore(pdb.getResidue(jNum), pdb.getResidue(iNum), resDist, graphGeom);
						if (val!=-666666){
							score += val;
							numSummedScores++;
						}
					} catch (SQLException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}				
			}
//			System.out.print(iNum+"  "+jNum+"  "+"Intermediate logOddsScores are : "+score+"\n");
		}
//		System.out.println("Final logOddsScore for protein "+pdb.getPdbCode()+" is : "+score+" norm:"+(score/i)+"  for "+numSummedScores+" chosen contacts");
		
		if (this.normaliseScore)
			score = score/numSummedScores; // normalise score
		
		return score;
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

}
