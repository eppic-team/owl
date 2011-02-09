package owl.core.structure.graphs;

//import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.SQLException;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.Map.Entry;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import edu.uci.ics.jung.graph.util.Pair;

import owl.core.structure.AminoAcid;
import owl.core.structure.Residue;
import owl.core.structure.features.SecStrucElement;
import owl.gmbp.CMPdb_sphoxel;
import owl.gmbp.CSVhandler;
import owl.gmbp.GmbpGeometry;

//RIGGeometry: class that computes and stores geometry information for graph based on rotation and translation invariant framework.

public class RIGGeometry {
	
	private RIGraph graph;
	private TreeMap<Integer, Residue> residues;
//	private HashMap<String,Vector3d> coord_sph_rotated;  // holds geometry (translated and rotated contact coordinates of CA-position
//						//	of j-Residue with respect to central iResidue of each edge of the graph, defined by residue serials of edge)
	private HashMap<Pair<Integer>, Vector3d> coord_sph_rotated; // holds geometry (translated and rotated contact coordinates of CA-position
							//of j-Residue with respect to central iResidue of each edge of the graph, defined by residue serials of edge)
	
	private String atomType = "CA";
	
	
	/**
	 * Constructs a RIGraph with a sequence but no edges
	 * @param sequence
	 */
	public RIGGeometry(RIGraph graph, TreeMap<Integer, Residue> residues) {
		this.graph = graph;
		this.residues = residues;
		this.atomType = "CA";
		
		initialiseGeometry();
	}
	
	/**
	 * Constructs a RIGraph with a sequence but no edges
	 * @param sequence
	 */
	public RIGGeometry(RIGraph graph, TreeMap<Integer, Residue> residues, String atomType) {
		this.graph = graph;
		this.residues = residues;
		this.atomType = atomType;
		
		initialiseGeometry();
	}
	
	private void initialiseGeometry(){
		
//		System.out.println("Geometry of graph with "+this.graph.getEdgeCount()+" edges:");
		int edgeNum = 0;
		GmbpGeometry gmbp = new GmbpGeometry();

		PrintWriter logOut = null;
//		String logFileFN = "/Volumes/StruPPi/CASP8/server_models/Geometry.csv";
//		logFileFN = "/Users/vehlow/GeometryTestTest.csv";
//		try {
//			logOut = new PrintWriter(new FileWriter(logFileFN));
//		} catch (IOException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
		
//		coord_sph_rotated = new HashMap<String,Vector3d>();
		coord_sph_rotated = new HashMap<Pair<Integer>, Vector3d>();
		for (RIGEdge edge:graph.getEdges()) {
			edgeNum++;
			// extract (x,y,z) coordinates for nodes of both end of edge
			Pair<RIGNode> nodes = this.graph.getEndpoints(edge);
			RIGNode iNode = nodes.getFirst();
			RIGNode jNode = nodes.getSecond();
			int iNum = iNode.getResidueSerial();
			int jNum = jNode.getResidueSerial();
			String iResType = iNode.getResidueType();
			String jResType = jNode.getResidueType();
			
			
			Residue iRes = this.residues.get(iNum);
			Residue jRes = this.residues.get(jNum);
				
			
			// translate coordinates with rotation and translation invariant framework
			HashMap<String,Point3d> iCoord = gmbp.getTheResidueCoord(iRes);
			HashMap<String,Point3d> jCoord = gmbp.getTheResidueCoord(jRes);
						
			Vector3d coord_I = new Vector3d(0,0,0);
			Vector3d coord_J = new Vector3d(0,0,0);
			// LEAVE this EDGE and CONTINUE with next one if the above condition is not satisfied.
			HashMap<String,Vector3d> atom_coord_rotated_I = new HashMap<String, Vector3d>();
			HashMap<String,Vector3d> atom_coord_rotated_J = new HashMap<String, Vector3d>();
			if (iCoord!=null && jCoord!=null) {
				// METHOD FOR EXTRACTING THE NEIGHBOR'S TRANSLATED-ROTATED COORDINATES //
				atom_coord_rotated_I = gmbp.getNeighborsTransRotatedCoord(jCoord, iCoord, jResType, iResType, jRes, iRes, true);
				atom_coord_rotated_J = gmbp.getNeighborsTransRotatedCoord(iCoord, jCoord, iResType, jResType, iRes, jRes, true);
				//========= C-ALPHA-coordinates ================// 
				coord_I = atom_coord_rotated_I.get(this.atomType); //("CA");
				coord_J = atom_coord_rotated_J.get(this.atomType); //("CA");
			}
			else {
				continue;
			}
			
			// GET the SPHERICAL COORDINATES for CA, C, CB, and CG using METHOD "getSphericalFromCartesian", 
			// if these coordinates exist (are non-zero).
			Vector3d coord_sph_I = new Vector3d(0,0,0);
			Vector3d coord_sph_J = new Vector3d(0,0,0);
			if (!coord_J.equals(new Vector3d(0.0,0.0,0.0))) {
				coord_sph_J = gmbp.getSphericalFromCartesian(coord_J); // (r,theta,phi) // (r, phi,lambda)
			}
			if (!coord_I.equals(new Vector3d(0.0,0.0,0.0))) {
				coord_sph_I = gmbp.getSphericalFromCartesian(coord_I); // (r,theta,phi) // (r, phi,lambda)
			}
			
			// Save translated and rotated coordinate of contact
			Pair<Integer> key = new Pair<Integer>(jNum, iNum);
			coord_sph_rotated.put(key, coord_sph_I);
			
			key = new Pair<Integer>(iNum, jNum);
			coord_sph_rotated.put(key, coord_sph_J);
			if (!coord_sph_rotated.containsKey(key)){
				System.out.println("Something is going wrong");
			}
			
//			System.out.println("TransRotCoord Cartesian: "+coord.x+","+coord.y+","+coord.z
//					+" SPH: "+coord_sph.x+","+coord_sph.y+","+coord_sph.z);
//			System.out.printf("TransRotCoord i->j Cartesian: %s Spherical: %s   j->i Cartesian: %s Spherical: %s \n", coord_J, coord_sph_J, coord_I, coord_sph_I);
			char iResT = AminoAcid.three2one(iResType);
			char jResT = AminoAcid.three2one(jResType);
			char iSS, jSS;
			SecStrucElement ssElem = iNode.getSecStrucElement();
			if (ssElem == null)
				iSS = 'A';
			else
				iSS = ssElem.getType();
			ssElem = jNode.getSecStrucElement();
			if (ssElem == null)
				jSS = 'A';
			else
				jSS = ssElem.getType();
			if (logOut!=null){
				logOut.println(iNum+","+jNum+","+iResT+","+jResT+","+iSS+","+jSS+","+coord_sph_J.x+","+coord_sph_J.y+","+coord_sph_J.z);
				logOut.println(jNum+","+iNum+","+jResT+","+iResT+","+jSS+","+iSS+","+coord_sph_I.x+","+coord_sph_I.y+","+coord_sph_I.z);
			}
//			System.out.printf("TransRotCoord i->j Spherical: %s   j->i Spherical: %s \n", coord_sph_I, coord_sph_J);
//			System.out.printf("TransRotCoord i->j Spherical: %s   j->i Spherical: %s \n", coord_sph_I, coord_sph_J);
		}
		if (logOut!=null)
			logOut.close();

	}
	
	public void printGeom(){
		for (Entry<Pair<Integer>, Vector3d> entry: this.coord_sph_rotated.entrySet()){
			System.out.println(entry.getKey().getFirst()+"_"+entry.getKey().getSecond()+":"+entry.getValue().x);
		}
	}
	
	public boolean isContactOnNorthernHemisphere(Residue iRes, Residue jRes){
		int iNum=iRes.getSerial();
        int jNum=jRes.getSerial();
        
        double phi = this.coord_sph_rotated.get(new Pair<Integer>(iNum, jNum)).y; //[0:Math.PI]
//        double lambda = this.coord_sph_rotated.get(new Pair<Integer>(iNum, jNum)).z;  //[-Math.PI:Math.PI]
        
        if (phi<(Math.PI/2))   // north
        	return true;
        else
        	return false;
	}
	
	public double getLogOddsScore(Residue iRes, Residue jRes, double resDist, String archiveFN) throws NumberFormatException, IOException{
		ZipFile zipfile = null;
        try {
            zipfile = new ZipFile(archiveFN);
        } catch (IOException e1) {
            // TODO Auto-generated catch block
            e1.printStackTrace();
        }
        
        String radiusPrefix="";
        if (/*2.0<=resDist &&*/ resDist<5.6)
            radiusPrefix="rSR";        
        if (5.6<=resDist && resDist<9.2)
            radiusPrefix="rMR";        
        if (9.2<=resDist /*&& resDist<12.8*/)
            radiusPrefix="rLR";
        
        String fn = "";
        String iSecSType="";

        if (iRes.getSsElem() == null)
            iSecSType = "a";
        else{
        	//
//        	String[] sstypes; // 
        	char[] types = {'H', 'S', 'O'};
        	String ssType = String.valueOf(iRes.getSsElem().getType());
        	if (types.toString().indexOf(ssType)>=0 && types.toString().indexOf(ssType)<types.length)
        		iSecSType = ssType.toLowerCase();
        	else {
//        		System.out.println("unknown ssType: "+ssType);
        		iSecSType = "a";
        	}
        }
        
        int iNum=iRes.getSerial();
        int jNum=jRes.getSerial();

        double score = 0;
        
        fn = fn+"sphoxelBG_"+iRes.getAaType().getOneLetterCode()+"-"+iSecSType+"_"+jRes.getAaType().getOneLetterCode()+"-"
            +"a"+"_"+radiusPrefix+".csv";
        
        ZipEntry zipentry = zipfile.getEntry(fn);
        if (zipentry!=null){
	        CSVhandler csv = new CSVhandler();
	        double bayesRatios [][][];        
	        bayesRatios = csv.readCSVfile3Ddouble(zipfile, zipentry); // dim = 72x144x3
	
	        int i1=(int)Math.floor((this.coord_sph_rotated.get(new Pair<Integer>(iNum, jNum)).y)/(Math.PI/72));
	        int j1=(int)Math.floor((this.coord_sph_rotated.get(new Pair<Integer>(iNum, jNum)).z)/(Math.PI/72))+72;  
	        
	        // sum up scores over tiles in direct surrounding (3x3)
//	        for (int i=i1-1; i<=i1+1; i++){
//	        	for (int j=j1-1; j<=j1+1; j++){
//	        		int iIndex=i, jIndex=j;
//	            	if (i<0)
//	            		iIndex = i1+1;
//	            	if (i>=bayesRatios.length)
//	            		iIndex = i1-1;
//	            	if (j<0)
//	            		jIndex = bayesRatios[0].length-1;
//	            	if (j>=bayesRatios[0].length)
//	            		jIndex = 0;
//	            	score += bayesRatios[iIndex][jIndex][0];            		
//	            }
//	        }	    
	        score = bayesRatios[i1][j1][0]; 
	        
        }
        else{
        	System.out.println("No SphoxelBG for "+fn);
        	return -666666;
        }	
       
        return   score;  
	}
	
	public double getLogOddsScore(Residue iRes, Residue jRes, double resDist) throws NumberFormatException, IOException, SQLException{
		CMPdb_sphoxel sphoxel = new CMPdb_sphoxel();
		
		double score = sphoxel.getLogOddsScore(iRes, jRes, resDist, this);
		
		return score;
	}
	
	/*
     * This is a method to get the LogOddsScore for all the contact pairs in the RIGraph file. This function accesses 
     * the precalculated log odds scores database in sphoxelBG.zip
     * @params: Residues ires, jRes, double resDist
     * @returns: double score for each contact. 
     * @throws: NumberFormatException, IOException
     */
	public double outputLogOddsScore(Residue iRes, Residue jRes, double resDist) throws NumberFormatException, IOException
        {
        
        String radiusPrefix="";
            if (2.0<=resDist && resDist<5.6)
                radiusPrefix="rSR";
            
            if (5.6<=resDist && resDist<9.2)
                radiusPrefix="rMR";
            
            if (9.2<=resDist /*&& resDist<12.8*/)
                radiusPrefix="rLR";
            
            
//            File dir1 = new File (".");        
//            String fn = "/amd/talyn/1/project/StruPPi/Saurabh/workspace/CMView/src/resources/sphoxelBG"; 
            //                fn = dir1.getCanonicalPath() + "/src/resources/sphoxelBG/";
//            String fn = "/project/StruPPi/Saurabh/workspace/CMView/src/resources/sphoxelBG/";
//            fn = "/Volumes/StruPPi/Saurabh/workspace/CMView/src/resources/sphoxelBG/";
            String zipP = "/Users/vehlow/Documents/workspace/CMView/src/resources/sphoxelBG/";
            String archiveFN = zipP + "SphoxelBGs.zip";
//            System.out.println("archiveFN= "+archiveFN);
            ZipFile zipfile = null;
            try {
                zipfile = new ZipFile(archiveFN);
            } catch (IOException e1) {
                // TODO Auto-generated catch block
                e1.printStackTrace();
            }
            String fn = "";
            String iSecSType="";

            if (iRes.getSsElem() == null)
                {
                iSecSType = "a";
                }
            else
                {
                iSecSType = String.valueOf(iRes.getSsElem().getType()).toLowerCase();
                }
            
            int i=iRes.getSerial();
            int j=jRes.getSerial();
            
            fn = fn+"sphoxelBG_"+iRes.getAaType().getOneLetterCode()+"-"+iSecSType+"_"+jRes.getAaType().getOneLetterCode()+"-"
                +"a"+"_"+radiusPrefix+".csv";
            System.out.println(i+"  "+j+"  "+resDist);
            System.out.println("Filename= "+fn);
            
            ZipEntry zipentry = zipfile.getEntry(fn);
            CSVhandler csv = new CSVhandler();
            double bayesRatios [][][];            
            
            bayesRatios = csv.readCSVfile3Ddouble(zipfile, zipentry);

            
//            double logOddsScores [][] = new double[bayesRatios.length][bayesRatios[0].length];
            
//            System.out.println(this.coord_sph_rotated.get(i+"_"+j).y+"  "+this.coord_sph_rotated.get(i+"_"+j).z);
//            int i1=(int)Math.floor((this.coord_sph_rotated.get(i+"_"+j).y)/(Math.PI/72));
//            int j1=(int)Math.floor((this.coord_sph_rotated.get(i+"_"+j).z)/(Math.PI/72))+72;  
            int i1=(int)Math.floor((this.coord_sph_rotated.get(new Pair<Integer>(i, j)).y)/(Math.PI/72));
            int j1=(int)Math.floor((this.coord_sph_rotated.get(new Pair<Integer>(i, j)).z)/(Math.PI/72))+72;  
            
            
//            for (int i2=0; i2<bayesRatios.length; i2++)
//                {
//                for (int j2=0; j2<bayesRatios[i2].length; j2++)
//                    {
//                    logOddsScores[i2][j2] = bayesRatios[i2][j2][0];
//                    }
//                }
            
//            System.out.println(i+"  "+j+"  "+logOddsScores[i1][j1]);
            return   bayesRatios[i1][j1][0];  //logOddsScores[i1][j1];            
        }
	
	// --------- getters ----------
//	public HashMap<String,Vector3d> getRotatedCoordOfContacts(){
//		return this.coord_sph_rotated;
//	}
	
	public HashMap<Pair<Integer>,Vector3d> getRotCoordOfContacts(){
		return this.coord_sph_rotated;
	}
	
	public String getAtomType(){
		return this.atomType;
	}

}
