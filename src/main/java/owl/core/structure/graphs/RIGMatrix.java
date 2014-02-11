
package owl.core.structure.graphs;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.Random;

import javax.vecmath.GMatrix; 

import edu.uci.ics.jung.graph.util.Pair;
import owl.core.structure.PdbChain;
import owl.core.structure.PdbAsymUnit;
import owl.core.structure.PdbLoadException;
import owl.core.structure.graphs.RIGNode;
import owl.core.structure.graphs.RIGraph;
import owl.core.util.FileFormatException;
import owl.core.util.MySQLConnection;

/**
 * Representation of a contact map in matrix form
 * 
 * @author lappe
 */

public class RIGMatrix {
	
	// Object definition by GMatrix M=(n x n), CharArra[n] for the sequence 
	private GMatrix M;
	private char[] S;
	private String pdbId;
	private String chainCode;
	private String ct;
	private double cutoff;
	static int maxRank=21;
	private static final String CIFREPODIR = "/path/to/mmCIF/gz/all/repo/dir";
	private static final String DEF_CT     = "Ca";//default contact type
	private static final double DEF_CUTOFF = 9.0;//default cutoff distance
	
	public RIGMatrix (){};//Zero-constructor
	
	public RIGMatrix (String pdbCode, String chainCode) throws FileFormatException, IOException, PdbLoadException{//
		this();
		setPdbId(pdbCode);
		setChainCode(chainCode);
		setContactType(DEF_CT);
		setCutoff(DEF_CUTOFF);
		loadFromFile();
	}

	public RIGMatrix( int n) { // Constructor for an empty RIGMatrix size n 
		M = new GMatrix( n, n); 
		M.setIdentity();
		S = new char[n];
		for (int i=0; i<n; i++) S[i]='-'; 
	} // end constructing empty RIGMatrix 
	
	public RIGMatrix( RIGraph G)  { // Constructing RIGMatrix from a RIGRaph object G 
		// creating new matrix size <nrofnodes x nrofnodes>
		int n = G.getObsLength(); 
		// System.out.println("n="+n); 
		M = new GMatrix( n, n); 
		G.convert2GMatrix(M); 
		// attach the sequence
		setSequence(G.getSequence());
		setPdbId(G.getPdbCode());
		setChainCode(G.getChainCode());
		setContactType(G.getContactType());
		setCutoff(G.getCutoff());
		 // next node serial i 
	} // end converting RIGraph RIGMatrix 	

	public RIGMatrix( RIGMatrix Original) { // Constructing RIGMatrix as copy from Oroginal 
		// creating new matrix size <nrofnodes x nrofnodes>
		int n = Original.M.getNumCol(); 
		// System.out.println("n="+n); 
		M = new GMatrix( n, n); 
		M.setZero(); 
		M.add( Original.M); 
		S = new char[n];
		setSequence(Original.getSequence());
		setPdbId(Original.pdbId);
		setChainCode(Original.chainCode);
		setContactType(Original.ct);
		setCutoff(Original.cutoff);
		 // next node serial i 
	} // end copy RIGMatrix 	
	
	public RIGMatrix( String fFile, int VL) {
		String si="", sj="", sw="", seq="", pdb="", pdbcc="", chain=".", model=".", ct="", cutoff="";
		int i=0, n=0;
		boolean doneHeader=false;
		if (VL>0) System.out.print("\nLoading from cmap file <"+fFile+">"); 
		try {
			BufferedReader in = new BufferedReader(new FileReader(fFile));
			String str;
			while ((str = in.readLine()) != null) {
				if (VL>0) System.out.print("\n"+(++i)+">"+str);
				if (str.startsWith("#")) {
					// Info field parsing 
					if (str.startsWith("#SEQUENCE:")) seq=str.substring(10).trim();
					if (str.startsWith("#PDB:")) pdb=str.substring(6);
					if (str.startsWith("#PDB CHAIN CODE:")) pdbcc=str.substring(17);
					if (str.startsWith("#CHAIN:")) chain=str.substring(8);
					if (str.startsWith("#MODEL:")) model=str.substring(8);
					if (str.startsWith("#CT:")) ct=str.substring(5);
					if (str.startsWith("#CUTOFF:")) cutoff=str.substring(9);
				} else {					
					if (VL>0 && doneHeader) { // verbosity level : print Header 
						doneHeader = true; 
						System.out.println("\nSaving Header :"); 
						System.out.println("#SEQUENCE:"+seq);
						System.out.println("#PDB:"+pdb);
						System.out.println("#PDB CHAIN CODE:"+pdbcc);
						System.out.println("#CHAIN:"+chain);
						System.out.println("#MODEL:"+model);
						System.out.println("#CT:"+ct);
						System.out.println("#CUTOFF:"+cutoff);	
						n=seq.length(); 
						M = new GMatrix( n, n); 
						M.setIdentity();
						setSequence(seq);
						setPdbId(pdb);
						setChainCode(pdbcc);
						setContactType(ct);
						setCutoff(Double.parseDouble(cutoff));
					} // end if 
					// parsing the contact lines  
					// i, j, w : edge line extraction 
					si=""; sj=""; sw="";
					si=str.substring( 0, str.indexOf("\t")); 
					if (VL>1) System.out.print("\nsi:"+si ); 
					sj=str.substring(str.indexOf("\t")+1,str.lastIndexOf("\t")); 
					if (VL>1) System.out.print("\tsj:"+sj);
					sw=str.substring(str.lastIndexOf("\t")+1); 
					if (VL>1) System.out.print("\tsw:"+sw);
					// store i, j, w 
					M.setElement( Integer.parseInt(si), Integer.parseInt(sj), Double.parseDouble(sw)); 
					M.setElement( Integer.parseInt(sj), Integer.parseInt(si), Double.parseDouble(sw)); 
				} // end if startswith # 
			} // end while 
			in.close();
		} catch (IOException e) {
			System.out.println("ERROR while loading "+fFile+" >> "+e);
		} // end try/catch 
	} // end of constructor loadfromCmapFile

	
	public void loadFromFile() throws IOException, FileFormatException, PdbLoadException{ 
		// pdbID, chainCode, ct & cutoff have to be pre-set for this to work properly
		// pass them as parameters and convert back as constructor method ? 
		// throw errors ?
		File cifFile = new File(System.getProperty("java.io.tmpdir"),pdbId+".cif");
		cifFile.deleteOnExit();
		PdbAsymUnit.grabCifFile(CIFREPODIR, null, pdbId, cifFile, false);				
		PdbAsymUnit fullpdb = new PdbAsymUnit(cifFile);
		
		PdbChain pdb = fullpdb.getChain(chainCode);
		RIGraph G = pdb.getRIGraph(ct, cutoff);
		if(S==null) setSequence(G.getSequence());
		M = new GMatrix(S.length,S.length);
		G.convert2GMatrix(M);
		
	}
	
	public void setSequence (String sequence){
		if(sequence==null) throw new NullPointerException ("No null String accepted!");
		else{
			int length = sequence.length();
			S = new char[length];
			for(int i = 0; i < length; i++){
				S[i] = sequence.charAt(i);
			}
		}
	}
	
	public void randomizeSequence () {
		int n = M.getNumCol(), j=0;
		char t; 
		Random RG = new Random(); 
		for(int i = 0; i < n; i++){
			j = RG.nextInt( n);
			t = S[i]; 
			S[i] = S[j];
			S[j] = t; 
		} // next aa 
	} // end randomizeSequence
	
	// simple setters and getters 
	public void setPdbId (String pdbId) { this.pdbId = new String (pdbId); }  
	public String getPdbId (){ return new String(pdbId); }
	public void setChainCode (String chainCode) { this.chainCode = new String (chainCode); }
	public String getChainCode () { return new String (chainCode); }
	public void setContactType (String ct) { this.ct = new String (ct); }
	public String getContactType () { return new String (ct); }
	public void setCutoff(double cutoff) { this.cutoff = cutoff; }
	public double getCutoff() { return cutoff; }	
	// handing over some basic GMatrix functionality 
	public void setM( GMatrix Mx) { M.set(Mx); }
	public GMatrix getM( ) { return M; }; 
	public void mul( RIGMatrix Mx) { M.mul( Mx.M); }
	public void add( RIGMatrix Mx) { M.add( Mx.M); }
	public void setZero() { M.setZero(); }
	public void setIdentity() { M.setIdentity(); }
	public void setElement( int i, int j, double v) { M.setElement( i, j, v); }
	public double getElement( int i, int j) { return M.getElement( i, j); }
	public int getNumCol() { return M.getNumCol(); }
	
	public String getSequence() {
		if(S==null) throw new NullPointerException ("Sequence field null!");
		else{
			String seq = new String ();
			int length = S.length;
			for(int i = 0; i < length; i++){
				seq+= new Character(S[i]).toString();
			}
			return seq;
		}			
	}

	public void printMatrixFormat( String form, String seqSep, RIGMatrix trueC) {
		int n=M.getNumCol(); 
		double v = 0.0; 
		System.out.print("  ");
		for (int i=0; i<n; i++) System.out.print(seqSep+S[i]);
		System.out.println();
		for (int i=0; i<n; i++) {
			System.out.print(S[i]+" ");
			for (int j=0; j<n; j++) {
				v = M.getElement(i, j);
				System.out.format( form, v);
				if (trueC.M.getElement(i,j) > 0.0) System.out.print("*");
				else System.out.print(" ");
			} // next j 
			System.out.println();
		} // next i
	} // end printMatrixFormat
	
	public void printMatrixFormat( String form, String seqSep) {
		int n=M.getNumCol(); 
		double v = 0.0; 
		System.out.print("  ");
		for (int i=0; i<n; i++) System.out.print(seqSep+S[i]);
		System.out.println();
		for (int i=0; i<n; i++) {
			System.out.print(S[i]+" ");
			for (int j=0; j<n; j++) {
				v = M.getElement(i, j); 
				System.out.format( form, v);
			} // next j 
			System.out.println();
		} // next i
	} // end printMatrixFormat

	public String getNbString( int i)  { // returns the nbstring according to entries of M[i]>0.0 
		String nbs="";
		for (int j=0; j<M.getNumCol(); j++)  { // j is in GMatrix coords: starts from 0
			if ( j==i ) nbs+="x"; 
			else if (M.getElement(i, j)>0.0) nbs+=S[j]; 
		} // next j 
		return nbs; 
	} // end getNbString
	
	public void listNbStrings( )  { // returns the nbstring according to entries of M[i]>0.0 
		for (int i=0; i<M.getNumCol(); i++)  { 
			System.out.println(i+" "+S[i]+" "+this.getNbString( i)); 
		} // next i 
	} // end listNbStrings
	
	public void randomize( double fp, double fn)  { // false positive & negative rates 
		int n=M.getNumCol(); 
		double r=0.0; 
		Random rand = new Random();
		for (int i=0; i<n; i++) {
			for (int j=i; j<n; j++) {
				// random number 0..1
				r = rand.nextDouble();
				if (M.getElement( i, j)>0)  { // i,j is a positive nonzero entry : apply FN rate 
					if (r<fn) { 
						M.setElement( i, j, 0.0);
						M.setElement( j, i, 0.0);
					} // end if we turn this contact as a FN 
				} else { // empty: apply FP rate 
					if (r<fp) { 
						M.setElement( i, j, 1.0);
						M.setElement( j, i, 1.0);
					} // end if we turn this non-contact into a FP entry  
				} // end if contact 
			} // next j 
		} // next i  
	} // end randomizeGMatrix
	
	public double getMax() { // maximal value in the Matrix 
		double max = M.getElement(0, 0);
		int n=M.getNumCol(); 
		double v = 0.0; 
		for (int i=0; i<n; i++) {
			for (int j=0; j<n; j++) {
				v = M.getElement(i, j); 
				if (v>max) max=v; 
			} // next j 
		} // next i
		return max; 
	} // get Max
	

	public Pair<Integer> locateMaxDiag(double v) { // locates first entry v with max.seq.separation (maxdiagonal)
		int d=0, i=0, j=0, n=M.getNumCol();
		Pair<Integer> ij = new Pair<Integer>(-1, -1);
		for (d=n; d>0; d--) { // go through diagonal d 
			j=d;  
			for (i=0; i<n-d; i++)  {
				if (M.getElement( i,j)==v) {
					// return i, j as a pair;
					ij = new Pair<Integer>(i, j);
					break; 
				}
				j++; 
			} // next i 
		} // next diagonal 
		return ij; 
	} // end locateMax 
	
	
	public void setMax( double maxval) { // sets all entries to given maximal value where value > cutoff maxval 
		int n=M.getNumCol(); 
		for (int i=0; i<n; i++) {
			for (int j=0; j<n; j++) {
				if ( M.getElement(i, j) > maxval) M.setElement(i, j, maxval); 
			} // next j 
		} // next i
	} // set Max
	
	public double getMin() { // minimal value 
		double min = M.getElement(0, 0);
		int n=M.getNumCol(); 
		double v = 0.0; 
		for (int i=0; i<n; i++) {
			for (int j=0; j<n; j++) {
				v = M.getElement(i, j); 
				if (v<min) min=v; 
			} // next j 
		} // next i
		return min; 
	} // get Min

	public void setMin( double minval) { // sets all entries to minimal value where value < cutoff minval 
		int n=M.getNumCol(); 
		for (int i=0; i<n; i++) {
			for (int j=0; j<n; j++) {
				if ( M.getElement(i, j) < minval) M.setElement(i, j, minval); 
			} // next j 
		} // next i
	} // set Min
	
	public double getSum() { // sums over all entries in M 
		double sum = 0.0;
		int n=M.getNumCol(); 
		for (int i=0; i<n; i++) {
			for (int j=i; j<n; j++) {
				double entry = M.getElement(i, j);
				if(i==j){ if(entry!=0) sum+=entry;}
				else{ if(entry!=0) sum+=2d*entry;}
			} // next j 
		} // next i
		return sum; 
	} // get Sum
	
	public double getSuMul( RIGMatrix Mask) { // sums over all (this * Mask) ==> i.e. sums the entries for all native contacts given its' binary contact map
		double sum = 0.0;
		int n=M.getNumCol(); 
		for (int i=0; i<n; i++) {
			for (int j=0; j<n; j++) {
				sum+= M.getElement(i, j)*Mask.M.getElement(i, j); 
			} // next j 
		} // next i
		return sum; 
	} // get SuMul

	public RIGMatrix overlap( RIGMatrix cmp) { // sums over all (this * Mask) ==> i.e. sums the entries for all native contacts given its' binary contact map
		int n=M.getNumCol(); 
		RIGMatrix inBoth = new RIGMatrix( cmp); 
		inBoth.setZero(); 
		for (int i=0; i<n; i++) {
			for (int j=0; j<n; j++) {
				if (M.getElement(i, j)>0 && cmp.M.getElement(i, j)>0) cmp.M.setElement(i, j, 1.0); 
			} // next j 
		} // next i
		return inBoth; 
	} // get overlap
	
	public int cmo( RIGMatrix cmp) { // sums over all (this * Mask) ==> i.e. sums the entries for all native contacts given its' binary contact map
		int count = 0;
		int n=M.getNumCol(); 
		for (int i=0; i<n; i++) {
			for (int j=0; j<n; j++) {
				if (M.getElement(i, j)>0 && cmp.M.getElement(i, j)>0) count++; 
			} // next j 
		} // next i
		return count; 
	} // get cmo
	
	public int cmo_weighted( RIGMatrix cmp) { // same as cmo but weighted by (abs(i-j)) hence long range contacts contribute more to the cmo
		int count = 0;
		int n=M.getNumCol(); 
		for (int i=0; i<n; i++) {
			for (int j=0; j<n; j++) {
				if (M.getElement(i, j)>0 && cmp.M.getElement(i, j)>0) count+=Math.abs(i-j); 
			} // next j 
		} // next i
		return count; 
	} // get cmo
	
	public void multiply( RIGMatrix muli) { // sums over all (this * Mask) ==> i.e. sums the entries for all native contacts given its' binary contact map
		int n=M.getNumCol(); 
		for (int i=0; i<n; i++) {
			for (int j=0; j<n; j++) {
				M.setElement(i, j, M.getElement(i, j)*muli.M.getElement(i, j)); 
			} // next j 
		} // next i 
	} // get Sum
	
	public double getAvg() { // average value of all entries in the Matrix 
		double sum = getSum();
		double n=(double)M.getNumCol();
		return sum*Math.pow(n,-2.0);
	} // get Avg

	public double getAvgN0() { // average value of all non-zero entries in the Matrix 
		double sum = 0.0, c=0.0, v=0.0;
		int n=M.getNumCol(); 
		for (int i=0; i<n; i++) {
			for (int j=0; j<n; j++) {
				v = M.getElement(i, j);
				sum +=v; 
				if(v!=0) c++; 	
			} // next j 
		} // next i
		v=0.0; 
		if (c!=0) v=(double)((double)sum/(double)c);
		return v; 	
	} // get AvgN0 

	public double getStdDevN0(){
		double avg = getAvgN0(), sumStdDev = 0D; // what does 0D do? 
		GMatrix dif = new GMatrix (M);
		double n = (double) dif.getNumRow();
		//double counter = 0.0;
		for(int i = 0; i < n-1;i++){
			for(int j = i; j < n; j++){
				if(i==j) sumStdDev += Math.pow(dif.getElement(i, j) - avg,2d);//square of difference average vs. actual entry
				else sumStdDev += 2.0*Math.pow(dif.getElement(i, j) - avg,2d);//square of difference average vs. actual entry
				//since dif is symmetric, only the upper half needs to be considered 
				//counter++;
			}
		}
		double sqrMin = Math.pow(n, 2d)-1;
		return Math.pow(sumStdDev/sqrMin,0.5);
	} // end getStdDevN0
	
	/**
	 * Computes the standard deviation of this instance.
	 * @return
	 */
	public double getStdDev (){
		double avg = getAvg(), sumStdDev = 0D;
		GMatrix dif = new GMatrix (M);
		double n = (double) dif.getNumRow();
		//double counter = 0.0;
		for(int i = 0; i < n-1;i++){
			for(int j = i; j < n; j++){
				if(i==j) sumStdDev += Math.pow(dif.getElement(i, j) - avg,2d);//square of difference average vs. actual entry
				else sumStdDev += 2.0*Math.pow(dif.getElement(i, j) - avg,2d);//square of difference average vs. actual entry
				//since dif is symmetric, only the upper half needs to be considered 
				//counter++;
			}
		}
		double sqrMin = Math.pow(n, 2d)-1;
		return Math.pow(sumStdDev/sqrMin,0.5);
	}
	
	public boolean isSymmetric() { // sanity-check, returns true if the Matrix is really symmetric 
		boolean sym = true;
		int n=M.getNumCol(); 
		for (int i=0; i<n; i++) {
			for (int j=i; j<n; j++) {
				if( M.getElement(i, j)!=M.getElement(j, i)) sym=false; 
			} // next j 
		} // next i
		return sym; 
	} // end isSymmetric() 
	
	public void setDiagonal( double val) { // sets all (i,i) to value val
		int n=M.getNumCol(); 
		for (int i=0; i<n; i++) {
			M.setElement(i, i, val); 
		} // next i
	} // end setDiagonal() 
	
	public void reScale( double newmin, double newmax)  { // rescales the entries in Matrix myM such that min..max maps to newmin..nexmax 
		int n=M.getNumCol(); 
		double v=0.0;
		// determine max and min entry 
		double max = getMax(); 
		double min = getMin(); 
		// rescaling min..max to 0..1
		if (max>min && newmax>newmin) {
			for (int i=0; i<n; i++) {
				for (int j=0; j<n; j++) {
					v = M.getElement(i, j);
					if (v==0) {
						v  = -1;
					} else {
						v = newmin+( ((v-min)/(max-min))*(newmax-newmin)) ;
					}
					M.setElement( i, j, v); 
				} // next j 
			} // next i 
		} // end if max> min 	
	} // end rescale 

	public int listNodeScores( MySQLConnection conx, int VL ) throws SQLException  { // returns the nbstring according to entries of M[i]>0.0 
		String nbs="x", centres="";
		int nullrank=0, rank=0, sumdelta=0;
		for (int i=0; i<M.getNumCol(); i++)  { 
			nbs = this.getNbString( i);
			centres= ""+S[i]; 
			nullrank=getRank( conx, "x", centres); 
			rank=getRank( conx, nbs, centres); 
			if (VL>0) System.out.println(i+"\t"+centres+"\t"+nullrank+"\t"+nbs+"\t"+rank+" \t>d= "+(nullrank-rank));
			sumdelta+=(nullrank-rank); 
		} // next i 
		if (VL>0) System.out.println("Sum node delta = "+sumdelta); 
		return sumdelta; 
	} // end listNbStrings

	public int getRank( MySQLConnection conn, String nbs, String centRes) throws SQLException {
		int rank=maxRank; 
		Statement stmt = conn.createStatement();
		ResultSet rsst = stmt.executeQuery("select locate('"+centRes+"', rvector) from mw.rvecs10 where nbstring='"+nbs+"';"); 
		if (rsst.next()) {	
			rank = rsst.getInt( 1); 
		} // end if 
		// if (rank==0) rank = maxRank;
		if (rank<=0) rank=maxRank; // rank = getRank( conn, "x", centRes); // as an alternative: deliver nullrank if we don't know; 
		rsst.close(); 
		stmt.close(); 
		return rank; 
	} // end of getRank 
		
	public RIGMatrix scoreDeltaMul( MySQLConnection conx, RIGMatrix S) throws SQLException { // calculates deltaRank[i,j]*S[i,j] for (i!=j) && S[i,j]!=0.0 
		String i_priorNbS, j_priorNbS, i_postNbS, j_postNbS, i_res, j_res;
		int i_priorank, i_postrank, j_priorank, j_postrank, i_delta, j_delta;
		int i=0, j=0, n=M.getNumCol(); 
		double current=0.0, Svalue=0.0, sumdelta=0.0; 
		RIGMatrix C = new RIGMatrix( this); // copy C to play out changes in cmap 
		RIGMatrix D = new RIGMatrix( this); // D contains the resulting deltaRank * S entries for all [i,j]>0.0 
		D.M.setZero();  
		for ( i=0; i<n; i++) {
			for ( j=i+1; j<n; j++) {
				Svalue=S.getElement(i, j); // stricter version would be to only score for the non-zero entries in this 
				if( Svalue!=0.0) {	// there is a score S and the contact is in the Originally selected subset  
					// save current status
					current=C.M.getElement( i, j); 
					// get rank with and without (i,j), then restore previous status 
					// Unsetting (i,j) gives prior
					// System.out.print("D["+i+","+j+"]=");
					C.M.setElement(i, j, 0.0); 
					C.M.setElement(j, i, 0.0); 
					i_priorNbS = C.getNbString( i);
					j_priorNbS = C.getNbString( j);
					// Set i,j in s gives posterior 
					C.M.setElement(i, j, 1.0); 
					C.M.setElement(j, i, 1.0); 
					i_postNbS = C.getNbString( i);
					j_postNbS = C.getNbString( j);
					// determine prior&post-ranks for (i,j) 
					i_res= new String( C.S, i, 1); 
					i_priorank= getRank( conx, i_priorNbS, i_res); 
					i_postrank= getRank( conx, i_postNbS, i_res);
					if (i_priorank<maxRank && i_postrank<maxRank) i_delta=i_priorank-i_postrank;
					else i_delta=0; // only calc delta if both prior and post have a valid rank 
					j_res= new String( C.S, j, 1); 
					j_priorank= getRank( conx, j_priorNbS, j_res); 
					j_postrank= getRank( conx, j_postNbS, j_res); 
					if (j_priorank<maxRank && j_postrank<maxRank) j_delta=j_priorank-j_postrank;
					else j_delta=0; // only calc delta if both prior and post have a valid rank 
					j_delta=j_priorank-j_postrank;
					sumdelta=((double)((double)i_delta+(double)j_delta)) * Svalue;
					D.M.setElement( i, j, sumdelta); 
					D.M.setElement( j, i, sumdelta);
					// reset copy C to original 
					C.M.setElement(i, j, current); 
					C.M.setElement(j, i, current);
				} // end if SVlaue > 0 hence calculation necessary 
			} // next j 
		} // next i
		return D; 
	} // end calculateDelta
	

	public RIGMatrix scoreDeltaRank( MySQLConnection conx ) throws SQLException { // calculates deltaRank for each i,j (i!=j)  
		String i_priorNbS, j_priorNbS, i_postNbS, j_postNbS, i_res, j_res;
		int i_priorank, i_postrank, j_priorank, j_postrank, i_delta, j_delta, sumdelta;
		int i=0, j=0, n=M.getNumCol(); 
		double current=0.0; 
		RIGMatrix C = new RIGMatrix( this); // copy C to play out changes in cmap 
		RIGMatrix D = new RIGMatrix( this); // D contains the resulting deltaRank entries 
		D.M.setZero();  
		for ( i=0; i<n; i++) {
			for ( j=i+1; j<n; j++) {
				// save current status
				current=C.M.getElement( i, j); 
				// get rank with and without (i,j), then restore previous status from S 
				// Unsetting (i,j) gives prior
				// System.out.print("D["+i+","+j+"]=");
				C.M.setElement(i, j, 0.0); 
				C.M.setElement(j, i, 0.0); 
				i_priorNbS = C.getNbString( i);
				j_priorNbS = C.getNbString( j);
				// Set i,j in s gives posterior 
				C.M.setElement(i, j, 1.0); 
				C.M.setElement(j, i, 1.0); 
				i_postNbS = C.getNbString( i);
				j_postNbS = C.getNbString( j);
				
				// determine prior&post-ranks for (i,j) 
				i_res= new String( C.S, i, 1); 
				i_priorank= getRank( conx, i_priorNbS, i_res); 
				i_postrank= getRank( conx, i_postNbS, i_res);
				i_delta=i_priorank-i_postrank; 
				j_res= new String( C.S, j, 1); 
				j_priorank= getRank( conx, j_priorNbS, j_res); 
				j_postrank= getRank( conx, j_postNbS, j_res); 
				j_delta=j_priorank-j_postrank;
				sumdelta=i_delta+j_delta;

//				System.out.println("i_priorNbS "+i_priorNbS);
//				System.out.println("j_priorNbS "+j_priorNbS);
//				System.out.println("i_postNbS "+i_postNbS);
//				System.out.println("j_postNbS "+j_postNbS);
//				System.out.println("i_priorank "+i_priorank);
//				System.out.println("i_postrank "+i_postrank);
//				System.out.println("i_delta "+i_delta);
//				System.out.println("j_priorank "+j_priorank);
//				System.out.println("j_postrank "+j_postrank);
//				System.out.println("j_delta "+j_delta);
//				System.out.println("sumdelta "+sumdelta);
//				System.out.println(" "+i+" "+i_res+" ["+i_priorNbS+"-"+i_postNbS+"]=("+i_priorank+"-"+i_postrank+")="+i_delta+" + "+j+" "+j_res+" ["+j_priorNbS+"-"+j_postNbS+"]=("+j_priorank+"-"+j_postrank+")="+j_delta+" ==> "+sumdelta); 
				
				D.M.setElement( i, j, sumdelta); 
				D.M.setElement( j, i, sumdelta);
				// reset copy C to original 
				C.M.setElement(i, j, current); 
				C.M.setElement(j, i, current);
			} // next j 
		} // next i
		return D; 
	} // end calculateDelta
	
	
	/**
	 * Converts an RIGMatrix instance to an RIGraph instance. Note, that
	 * only contacts present in the RIGMatrix instance are used as contacts
	 * in the RIGraph.  
	 * @return an RIGraph
	 */
	public RIGraph convert2RIGraph (){
		RIGraph rig = new RIGraph();
		rig.setChainCode(chainCode);
		rig.setContactType(ct);
		rig.setCutoff(cutoff);
		rig.setPdbCode(pdbId);
		rig.setSequence(getSequence());
		int length = M.getNumCol();
		for(int i = 0; i < length-1;i++){
			for(int j = i+1; j < length; j++){
				if(M.getElement(i,j)!=0){
					String aa1 = (new Character(S[i])).toString();
					String aa2 = (new Character(S[j])).toString();
					RIGNode node1 = new RIGNode (i+1,aa1);
					RIGNode node2 = new RIGNode (j+1,aa2);
					rig.addVertex(node1);
					rig.addVertex(node2);
					rig.addEdgeIJ(i+1, j+1);
				}
			}
		}
		return rig;
	}
	
	public int makeBinaryByThreshold( double t) { // sets every Mij=1 if Mij>=t, Mij=0 otherwise  
		int n=M.getNumCol(), c=0; 
		for (int i=0; i<n; i++) {
			for (int j=0; j<n; j++) {
				if (M.getElement( i, j) >= t)  { // i,j is above t : set to 1 
					M.setElement( i, j, 1.0);
					c++; // increment counter of made contacts 
				} else { // set to 0  
					M.setElement( i, j, 0.0);
				} // end if above t 
				if (i==j) M.setElement( i, j, 1.0); // keep the backbone connectivity intact 
			} // next j 
		} // next i  
		return c; 
	} // end randomizeGMatrix
	
	public int calcTPFN( String what, RIGMatrix tru, double t) { // counts entries in M as positive if >= threshold t 
		// by comparing myM at threshold t with truM 
		int TP=0, FP=0, FN=0, TN=0, r=0;  
		int n=M.getNumCol(); 
		double v=0.0, w=0.0; 

		for (int i=0; i<n; i++) {
			for (int j=0; j<n; j++) {
				v = tru.M.getElement( i, j); 
				w = M.getElement(i, j);
				if (w>=t)  { // we would predict a positive since w>=t
					if (v>0) TP++; 
					else FP++; 
				} else 	 { // we would predict a negative 
					if (v>0) FN++; 
					else TN++; 
				} // end if positive 
			} // next j 
		} // next i 
		// System.out.format("\tTP %4d\tFP %4d\tFN %4d\tTN %4d", TP, FP, FN, TN); 
		if (what=="TP") r= TP; 
		if (what=="FP") r= FP;
		if (what=="FN") r= FN;
		if (what=="TN") r= TN;
		return r; 
	} // end count 
	
	/**
	 * @param args
	 * load PDB -> orig. adj. Matrix 
	 * create subset S (several alternative Methods, seq.separation, per diagonal, uniform noise etc ...
	 * score i.e. sum(deltaRank * CNSize) for entries in S, whole Matrix, CNSize > 3 or > avg.  
	 * @throws FileFormatException 
	 */
	public static void main(String[] args) throws PdbLoadException, IOException, SQLException, FileFormatException {
		MySQLConnection conn;
		conn = new MySQLConnection( "localhost", "lappe", "apple", "mw");
		String cType = "Cb"; // contact type like "Ca"/"Cb"/"ALL"
		double cutoff = 7.5; 
		String myCode = args[0]; 
		String myChain = args[1];
		String infile = "/Users/lappe/_PDBs/"+myCode+".pdb";	
		double fp=0.03, fn=0.87; 
		// Load PDB structure from file 
		System.out.println("Loading PDB "+myCode+" chain "+myChain+" from file "+infile+"."); 
		PdbAsymUnit fullpdb = new PdbAsymUnit(new File(infile));
		System.out.println("original PdbChain object created.");
		PdbChain original = fullpdb.getChain(myChain); 
		System.out.println("chain loaded.");
		RIGraph origraph = original.getRIGraph(cType, cutoff);
		System.out.println("Converted to RIG at contactDef("+cType+","+cutoff+")"); 	
		// convert to Matrix 
		RIGMatrix A = new RIGMatrix( origraph); 
		System.out.println("original contact Matrix A =\n"); 
		A.printMatrixFormat(" %.0f", " ");
		// A.listNbStrings(); 
		// generate random subset
		System.out.println("Subset at fp="+fp+" and fn="+fn+" Sub =\n");
		RIGMatrix Sub = new RIGMatrix( A); 
		Sub.randomize( fp, fn); 		
		Sub.printMatrixFormat(" %.0f", " ");
		// Sub.listNbStrings(); 
		
		// Scoring prototype 
		RIGMatrix Delta = new RIGMatrix( Sub);
		Delta=Sub.scoreDeltaRank( conn); 
		System.out.println("deltaRank Matrix D =\n"); 
		Delta.printMatrixFormat(" %3.0f", "   ");
		System.out.println( "sumDeltaRank ="+Delta.getSum()); 
		System.out.println( "native deltaRankScore "+Delta.getSuMul(A)+"/ #of contacts="+A.getSum()+" = "+Delta.getSuMul(A)/A.getSum());
		System.out.println( "subset deltaRankScore "+Delta.getSuMul(Sub)+"/ #of contacts="+Sub.getSum()+" = "+Delta.getSuMul(Sub)/Sub.getSum());
	
		// calcscore only for selected subset weighted by S^2
		RIGMatrix Square= new RIGMatrix(Sub); // SCube = Sub
		Square.M.mul( Sub.M); // SCube = Sub^2
		System.out.println("SubSquared =\n");
		Square.printMatrixFormat(" %3.0f", "   ");
		// SQuare & SCubed 
		RIGMatrix SCubed = new RIGMatrix( Square); // SCube = Sub^3
		SCubed.M.mul( Sub.M); // SCube = Sub^2
		System.out.println("SubCubed =\n"); 
		SCubed.printMatrixFormat(" %3.0f", "   ");
		
		SCubed.multiply(Delta); 
		System.out.println("SubCubed * Delta =\n"); 
		SCubed.printMatrixFormat(" %3.0f", "   ");
		
		SCubed.multiply( A); 
		System.out.println("Sum SCubed * Delta (native) = "+SCubed.getSum()+"\n"); 
		SCubed.printMatrixFormat(" %3.0f", "   ");
		
		SCubed.multiply( Sub); 
		System.out.println("Sum SCubed * Delta (selected Subset) = "+SCubed.getSum()+"\n"); 
		SCubed.printMatrixFormat(" %3.0f", "   ");
		
		// local search (gradient) to optimize position for each contact in its vincinity given SCubed*delta
		
		RIGraph outgraph = SCubed.convert2RIGraph();
		outgraph.writeToFile(infile+".cmap");  
		
		// still to implement: make binary by threshold 
		// still to implement: average over all nonzero entries makes only sense for strictly positive matrices  
		// implement scoring function(s) outside i.e. optimats
		// sum over all = sum ( this); 
		/* 3. loop to get better 
		 * 3.1. local moves / gradient descent / deduce promising next steps from current setting rather than blind guessing 
		 * 4. more optimizations by GA / breeding / parallelization -> requires functions to read from / write to DB 
		 */ 
		
		conn.close(); // closing the Database connection
		/*
		MySQLConnection con = new MySQLConnection();
		PdbChain pdb = new PdbasePdb("1bkr",DEF_DB, con);
		pdb.load("A");
		RIGraph rig = pdb.getRIGraph("Ca",9.0);
		con.close();
		RIGMatrix mat = new RIGMatrix ("1bkr","A");
		RIGraph cp = mat.convert2RIGraph();
		double av = mat.getAvg();
		double st = mat.getStdDev();
		double er = mat.CMError();		
		System.out.println("average: "+av+"\tstandard deviation: "+st+"\terror: "+er+"\n"+rig.getEdges().containsAll(cp.getEdges())); 
		*/ 
	} // end main 

	public void normalize() {
		// TODO Auto-generated method stub
		
	}
} // end class 
