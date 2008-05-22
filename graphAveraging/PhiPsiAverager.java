package graphAveraging;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.TreeMap;
import java.util.TreeSet;

import proteinstructure.Alignment;
import proteinstructure.Interval;
import proteinstructure.PairwiseSequenceAlignment;
import proteinstructure.Pdb;
import proteinstructure.Template;
import proteinstructure.TemplateList;
import proteinstructure.PairwiseSequenceAlignment.PairwiseSequenceAlignmentException;
import tools.MySQLConnection;

public class PhiPsiAverager {

	private Alignment al;								// alignment
	private ArrayList<TemplateWithPhiPsi> templates;

	public PhiPsiAverager(TemplateList templates, Alignment aln) {

		this.al = aln;
		if (!templates.isPdbDataLoaded()) {
			throw new IllegalArgumentException("TemplateList passed to PhiPsiAverager constructor must have PDB data loaded");
		}
		checkSequences(templates);
		getPhiPsi(templates);

	}
	
	/**
	 * Gets the consensus phi/psi angles mapped onto the sequence of the given targetTag
	 * If the targetTag is not in the Alignment an IllegalArgumentException is thrown
	 * @param threshold
	 * @param angleInterval
	 * @param targetTag
	 * @return
	 * @throws IllegalArgumentException if targetTag is not present in this.al
	 */
	public TreeMap<Integer, IntervalCandidate[]> getConsensusPhiPsiOnTarget(double threshold, int angleInterval, String targetTag) {
		if (!al.hasTag(targetTag)) throw new IllegalArgumentException("Given targetTag is not present in alignment");
		TreeMap<Integer, IntervalCandidate[]> phiPsiConsensus = this.getConsensusPhiPsi(threshold, angleInterval);
		TreeMap<Integer, IntervalCandidate[]> phiPsiConsOnTarget = new TreeMap<Integer, IntervalCandidate[]>();
		for (int i:phiPsiConsensus.keySet()) {
			int resser = al.al2seq(targetTag, i);
			if (resser!=-1) {
				phiPsiConsOnTarget.put(resser, phiPsiConsensus.get(i));
			}
		}
		return phiPsiConsOnTarget;
	}
	
	/**
	 * Gets the consensus phi/psi angles for each alignment column with phi/psi consensus.
	 * The return is a TreeMap with keys alignment positions, values an array of size 2 with the phi/psi IntervalCandidates. 
	 * Columns without phi/psi consensus are not present in the TreeMap
	 * @param threshold a value >= 0.5
	 * @param angleInterval
	 * @return
	 * @throws IllegalArgumentException if threshold < 0.5
	 */
	public TreeMap<Integer, IntervalCandidate[]> getConsensusPhiPsi(double threshold, int angleInterval) {

		if (threshold < 0.5) {
			throw new IllegalArgumentException("Threshold for consensus ph/psi must be above or equals to 0.5");
		}

		// we round DOWN given threshold. So for the 2 cases even/odd and limit case threshold=0.5: 
		// - odd case  : 3 templates => voteThreshold=1, then we apply cutoff (see getInterval) with a '>' so that 2 falls within threshold
		// - even case : 4 templates => voteThreshold=2, then we applu cutoff (see getInterval) with a '>' so that 2 doesn't fall within threshold 
		int voteThreshold = (int) (templates.size()*threshold); // the int casting rounds down
		
		TreeMap<Integer, IntervalCandidate[]> bounds = new TreeMap<Integer, IntervalCandidate[]>();

		// we go through each column i in the alignment and find the consensus per column
		for (int i=1; i<=al.getAlignmentLength(); i++) {
			double[] phiAnglesColumnI = new double[templates.size()];
			double[] psiAnglesColumnI = new double[templates.size()];
			for (int j=0;j<templates.size();j++) {
				TemplateWithPhiPsi template = templates.get(j);
				int resser = al.al2seq(template.getId(), i);

				phiAnglesColumnI[j] = Double.NaN;
				psiAnglesColumnI[j] = Double.NaN;
				if (resser!=-1) { // to skip gaps
					if (template.hasPhiPsiAngles(resser)) { // some columns won't have angle data because of unobserved i, i-1 or i+1 residue. Or for N and C-terminals
						phiAnglesColumnI[j] = template.getPhiAngle(resser);
						psiAnglesColumnI[j] = template.getPsiAngle(resser);
					}
				} 
			}
			
			IntervalCandidate intervPhi = getInterval(phiAnglesColumnI, voteThreshold, angleInterval);
			IntervalCandidate intervPsi = null;
			if (intervPhi!=null) { 
				// so there was consensus in phi angles, now we need to check whether psi angles agree for same voters
				intervPsi = getIntervalPsiGivenPhi(psiAnglesColumnI, intervPhi, angleInterval);
				if (intervPsi!=null) { 
					IntervalCandidate[] phipsiintervs = {intervPhi, intervPsi};
					bounds.put(i, phipsiintervs);
				}
			}
		}
		return bounds;
	}

	/**
	 * Gets the IntervalCandidate for the given a list of psi angles values for a 
	 * column in the alignment and a phi-IntervalCandidate.
	 * We get each of the psi values corresponding to each (phi) voted template, if they are within
	 * a angleInterval window then they are taken as our IntervalCandidate (recentering first 
	 * the interval at max-min/2)
	 * 
	 * @param anglesInColumn
	 * @param intervPhi
	 * @param angleInterval 
	 * @return the IntervalCandidate or null if the psi values don't fit withing an angleInterval window
	 */
	private IntervalCandidate getIntervalPsiGivenPhi(double[] anglesInColumn, IntervalCandidate intervPhi, int angleInterval) {
		// we use a psiValues ArrayList to put the psi values corresponding to the given phi voters 
		ArrayList<Double> psiValues = new ArrayList<Double>();
		for (int voterIndex:intervPhi.getVotersIndices()) {
			psiValues.add(anglesInColumn[voterIndex]);
		}
		// now psiValues contains one psi angle for each of the phi voters. We want to see if the fall within the angleInterval cutoff 
		double intervWidth = Collections.max(psiValues) - Collections.min(psiValues);
		if (intervWidth<=angleInterval) {
			// first we initialise with a dummy 0,0 interval, then we use the recenter method to get the interval
			IntervalCandidate intervPsi = 
				new IntervalCandidate(0, 0, intervPhi.voters, psiValues,intervPhi.votersIndices, intervPhi.voteCount);
			intervPsi.reCenterInterval(angleInterval);
			return intervPsi;
		}
		return null;
	}
	
	/**
	 * Gets the biggest (max voters) above voteThreshold IntervalCandidate given 
	 * a list of angle values for a column of the alignment.
	 * We slide a window of size angleInterval over the whole range of angles (-180,180) 
	 * with an integer step (default 1, see step variable) to do an exhaustive (but inefficient!) search.
	 * The final IntervalCandidate contains the interval (of size angleInterval and centered at (max-min)/2), 
	 * the number of votes, the voters and the indices of the voters (same indices as in the TemplateList)
	 * 
	 * @param anglesInColumn
	 * @param voteThreshold
	 * @param angleInterval
	 * @return
	 */
	private IntervalCandidate getInterval(double[] anglesInColumn, int voteThreshold, int angleInterval) {
		// We use a TreeSet to store all interval candidates, so that duplicates are eliminated.
		// Duplicates are defined according to IntervalCandidate.equals(): 2 Intervals are the same if they have the 
		// same size (voteCount) and same members (voters) 
		// Second: with the TreeSet the IntervalCandidates will be sorted according to IntervalCandidate.compareTo(): the 
		// comparison is based purely on the voteCount. Thus it's easy to get out of the set the one with the maximum count.
		TreeSet<IntervalCandidate> intervalCandidates = new TreeSet<IntervalCandidate>();
		int step = 1;
		for (int begin=-180;begin<=180-angleInterval;begin+=step) {
			int end = begin + angleInterval;
			IntervalCandidate intervalCandidate = new IntervalCandidate(begin, end);
			for (int j=0;j<anglesInColumn.length;j++) {
				if (anglesInColumn[j]<=end && anglesInColumn[j]>=begin) {
					intervalCandidate.addVoter(j, templates.get(j).getId(), anglesInColumn[j]);
				}
			}
			if (intervalCandidate.getVoteCount()>voteThreshold) {  // note we use strictly bigger, see comment in getConsensusPhiPsi
				intervalCandidates.add(intervalCandidate);
			}
		}

		if (intervalCandidates.size()>0) {
			// As we are allowing thresholds only above 50% the max voters candidate will 
			// be unique (i.e. same number of voters, same voters, as defined by IntervalCandidate.equals()), there's no possible duplication.
			
			// Because intervalCandidates is a TreeSet sorted according to IntervalCandidate.compareTo, the last element is the one with the max votes
			IntervalCandidate maxVotersCandidate = intervalCandidates.last();
			// now we change the interval to be centered around the center of min/max values
			maxVotersCandidate.reCenterInterval(angleInterval);
			
			// debugging
			//System.out.print(intervalCandidates.size()+" interv. candidates. Max consensus: "+maxVotersCandidate.getVoteCount()+", voters: "+maxVotersCandidate.getVoters()+", vote counts: ");
			//for (IntervalCandidate ic: intervalCandidates) {
			//	System.out.print(ic.getVoteCount()+" ");
			//}
			//System.out.println();
			return maxVotersCandidate;
		}
		return null;
	}
	
	/**
	 * A class to represent a consensus phi/psi angle interval.
	 * Keeps information of the interval, the voteCount, voters and angle values.
	 * Its comparator is on voteCount and its equal is on voteCount, voters and values 
	 */
	public class IntervalCandidate implements Comparable<IntervalCandidate> {
		private Interval interval;
		private int voteCount;
		private ArrayList<String> voters;
		private ArrayList<Integer> votersIndices;
		private ArrayList<Double> values; // the angle values within this interval
		public IntervalCandidate (int beg, int end) {
			this.interval = new Interval(beg,end);
			voters = new ArrayList<String>();
			values = new ArrayList<Double>();
			votersIndices = new ArrayList<Integer>();
			voteCount = 0;
		}
		public IntervalCandidate (int beg, int end, ArrayList<String> voters, ArrayList<Double> values, ArrayList<Integer> votersIndices, int voteCount) {
			this.interval = new Interval(beg,end);
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
		public Interval getInterval() {
			return interval;
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
		/**
		 * Tells whether this IntervalCandidate has > < == voteCount than other
		 */
		public int compareTo(IntervalCandidate o) {
			return new Integer(voteCount).compareTo(o.getVoteCount());
		}
		/**
		 * Returns true if this IntervalCandidate has same voteCount, same voters and same values than other
		 */
		public boolean equals(Object o) {
			if (! (o instanceof IntervalCandidate)) return false;
			IntervalCandidate other = (IntervalCandidate) o;
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
		 * Recenters the interval of this IntervalCandidate to (maxValue-minValue)/2
		 * @param angleInterval
		 */
		public void reCenterInterval(int angleInterval) {
			double max = Collections.max(values);
			double min = Collections.min(values);
			double intervWidth = max - min;
			// A 'graphical' view of this (x is center, '' the min/max(b/e) values and || the limits of our interval (B/E)): 
			// |   '  x  '   |
			// B   b     e   E
			// <------w------>
			//     <--d-->
			// so we want to get B, E for our final interval. Being d=interWidth and w=angleInterval then: B=b+d/2-w/2, E=e-d/2+w/2
			int B = (int) (min+intervWidth/2-angleInterval/2);
			int E = (int) (max-intervWidth/2+angleInterval/2);
			this.interval = new Interval(B,E);
		}
	}
	
	private class TemplateWithPhiPsi extends Template {
		private TreeMap<Integer,double[]> phipsiAngles;
		public TemplateWithPhiPsi (String id, Pdb pdb, TreeMap<Integer,double[]> phipsiAngles) {
			super(id);
			this.phipsiAngles = phipsiAngles;
			this.setPdb(pdb);
		}
		
		public double getPhiAngle(int resser) {
			return phipsiAngles.get(resser)[0];
		}
		public double getPsiAngle(int resser) {
			return phipsiAngles.get(resser)[1];
		}
		public boolean hasPhiPsiAngles(int resser) {
			return phipsiAngles.containsKey(resser);
		}
	}
	
	
	/**
	 * Gets the phi/psi angles for all templates putting them into this.templates.
	 * Checks first that all template have PDB data removing the ones that don't
	 */
	private void getPhiPsi(TemplateList templates) {
		this.templates = new ArrayList<TemplateWithPhiPsi>();
		for (Template template: templates) {
			if (template.hasPdbData()) {
				this.templates.add(new TemplateWithPhiPsi(template.getId(), template.getPdb(), template.getPdb().getAllPhiPsi()));
			}
		}
		if (this.templates.size()<templates.size()) {
			System.out.println("Using only "+this.templates.size()+" templates for psi/phi averaging, out of given "+templates.size());
		}
 	}
	
	/**
	 * Checks that tags and sequences are consistent between this.al and this.templates 
	 *
	 */
	private void checkSequences(TemplateList templates){
		for (Template template :templates){
			if (!al.hasTag(template.getId())){
				System.err.println("Alignment is missing template sequence "+template.getId()+", check the FASTA tags");
				// TODO throw exception
			}
		}
		//if (templates.size()!=al.getNumberOfSequences()-1){
		//	System.err.println("Number of sequences in alignment is different from number of templates +1 ");
		//	// TODO throw exception
		//}

		for (Template template:templates){
			if(!al.getSequenceNoGaps(template.getId()).equals(template.getPdb().getSequence())) {
				System.err.println("Sequence of template (from pdbase) "+template.getId()+" does not match sequence in alignment");
				System.err.println("Trying to align sequences of alignment vs pdbase sequence: ");
				try {
					PairwiseSequenceAlignment alCheck = new PairwiseSequenceAlignment(template.getPdb().getSequence(),al.getSequenceNoGaps(template.getId()),"pdbase","alignment");
					alCheck.printAlignment();
				} catch (PairwiseSequenceAlignmentException e) {
					System.err.println("Error while creating alignment check, can't display an alignment. The 2 sequences are: ");
					System.err.println("graph:     "+template.getPdb().getSequence());
					System.err.println("alignment: "+al.getSequenceNoGaps(template.getId()));
				}
				// TODO throw exception
			}			
		}
	}
	
	public void printAngles(String type, IntervalCandidate intCandidate) {		
		System.out.print(type+": ");
		
		ArrayList<Double> values = intCandidate.getValues();
		ArrayList<Integer> valsIntegers = new ArrayList<Integer>();
		for (double value:values) valsIntegers.add((int) value);
		for (int i=-360;i<=360;i++) { 
			if (i==1 || i==360) {
				System.out.print("#");
			}else if (valsIntegers.contains(i)) {
				System.out.print("'");
			} else if ((int)intCandidate.getInterval().beg==i || (int)intCandidate.getInterval().end==i) {
				System.out.print("|");
			} else {
				System.out.print(" ");
			}
		}
		System.out.println();
	}
	 
	
	/**
	 * To test the class
	 * @param args
	 */
	public static void main(String[] args) throws Exception {
		MySQLConnection conn = new MySQLConnection("white","duarte","nieve");
		
		File templatesFile = new File("/scratch/local/phipsi/T0332.templates");
		File alnFile = new File("/scratch/local/phipsi/T0332.target2templ.fasta");
		File psipredFile = new File("/scratch/local/phipsi/T0332.horiz");
		
		String pdbaseDb = "pdbase";
		TemplateList templates = new TemplateList(templatesFile);
		//Template temp1 = new Template("2f8aA");
		//Template temp2 = new Template("1gp1A");
		//Template temp3 = new Template("2gs3A");
		
		//TemplateList templates = new TemplateList();
		//templates.add(temp1);
		//templates.add(temp2);
		//templates.add(temp3);
		templates.loadPdbData(conn, pdbaseDb);
		
		Alignment aln = new Alignment(alnFile.getAbsolutePath(),"FASTA");
		System.out.println("Secondary structure matching: ");
		aln.writeSecStructMatching(System.out, "T0332", psipredFile , conn, pdbaseDb, "/project/StruPPi/bin/dssp");
		System.out.println();
		System.out.println("phi/psi angles of each sequence in alignment. Legend: line1 aln positions, line2 sequence positions, line3 phi, line4 psi");
		PhiPsiAverager ppAvrg = new PhiPsiAverager(templates,aln);
		
		TreeMap<Integer, IntervalCandidate[]> phipsibounds = ppAvrg.getConsensusPhiPsi(0.5, 20);
		
		// printing phi/psi angles in each of the templates		
		for (TemplateWithPhiPsi template:ppAvrg.templates) {
			TreeMap<Integer,double[]> phipsi = template.phipsiAngles;
			//String sequence = template.getPdb().getSequence();
			System.out.println(template.getId());
			// printing alignment positions
			for (int i=1; i<=aln.getAlignmentLength(); i++) {
				System.out.printf("%5d",i);
			}
			System.out.println();
			// printing residue serials
			for (int i=1; i<=aln.getAlignmentLength(); i++) {
				int resser = aln.al2seq(template.getId(), i);
				if (resser!=-1 && template.hasPhiPsiAngles(resser))
					System.out.printf("%5s",resser);
				else System.out.printf("%5s","");
			}
			System.out.println();
			// printing phis
			for (int i=1; i<=aln.getAlignmentLength(); i++) {
				int resser = aln.al2seq(template.getId(), i);
				if (resser!=-1 && template.hasPhiPsiAngles(resser))
					System.out.printf("%5.0f", phipsi.get(resser)[0]);
				else System.out.printf("%5s","");
			}
			System.out.println();
			// printing psis
			for (int i=1; i<=aln.getAlignmentLength(); i++) {
				int resser = aln.al2seq(template.getId(), i);
				if (resser!=-1 && template.hasPhiPsiAngles(resser))
					System.out.printf("%5.0f", phipsi.get(resser)[1]);
				else System.out.printf("%5s","");
			}
			System.out.println();
		}
		
		// printing consensus intervals
		System.out.println("Alignment positions that have phi/psi consensus.");
		System.out.println("Legend: col1= aln pos, columns before #: phi interval begin, phi interval end, phi values of members. Same for columns after # but for psi");
		for (int alnPos: phipsibounds.keySet()) {
			System.out.printf("%3d %4d %4d ",
					alnPos,
					phipsibounds.get(alnPos)[0].getInterval().beg,
					phipsibounds.get(alnPos)[0].getInterval().end);
			for (double val:phipsibounds.get(alnPos)[0].getValues()) {
				System.out.printf("%4.0f ",val);	
			}
			System.out.printf("# %4d %4d ",
					phipsibounds.get(alnPos)[1].getInterval().beg, 
					phipsibounds.get(alnPos)[1].getInterval().end);
			for (double val:phipsibounds.get(alnPos)[1].getValues()) {
				System.out.printf("%4.0f ",val);	
			}
			System.out.println();
		}
		
	}
}
