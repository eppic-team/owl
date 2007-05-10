package proteinstructure;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.ArrayList;
import java.lang.Math;

public class ContactMap {
	
	// NOTE: residueNums and residueTypes contain also the unobserved or non standard residues
	// 		 residues and nums2serials contain only observed standard amino acids
	//		 If someone wants to check if i_num or i is unobserved,
	//		 then this can be done by checking whether nums2serials has such key (for i_num) or value (for i)
	// WHERE DO WE USE i or i_num:
	// - i: for contact map iteration and sequence separation
	// - i_num: for residue types and common neighbours
	public Integer[] residueNums;  												// array with all residue nums
	public String[] residueTypes;											// array with all 1-letter code residue type
	public TreeMap<Integer,String> residues = new TreeMap<Integer,String>();// map from nums to types
	public HashMap<Integer,Integer> nums2serials = new HashMap<Integer,Integer>(); // map from residue nums to residue serials (i.e. the indexes of the arrays above)
		
	public int t;															// size of contact map
	public int effT;														// effective size of contact map (do not count the SKIPPED cells)
	public int l;															// length of contact map/protein
	public int numObsStandAA = 0;											// Nr. of standard observed amino acids
	public int numContacts = 0;												// Nr. of contacts
	
	public EdgeState[][] CM;	
	public HashMap<Integer,HashMap<Integer,TreeMap<Integer,String>>> CNs; // hashmap containing neighbourhoods, access a single one through getComNbs
	
	public ContactMap() {		
	}

	/**
	 * Constructs a ContactMap passing a list of contacts, the full residue sequence and a map of residue nums to types (only observed standard aas)
	 * @param contacts ArrayList of contacts (as Contact objects) 
	 * @param residues map from num to types (only observed standard aas)
	 * @param sequence String with full sequence (i.e. also unobserved and non-standard residues, denoted by X)
	 */
	public ContactMap(ArrayList<Contact> contacts, TreeMap<Integer,String> residues, String sequence) {
		this.residues=residues;
		// we create residueTypes and residueNums arrays from full sequence string
		residueTypes = new String[sequence.length()];
		residueNums = new Integer[sequence.length()];
		for (int i=0;i<sequence.length();i++){
			residueTypes[i] = String.valueOf(sequence.charAt(i));
			residueNums[i] = i+1;
		}		
		// we initialise nums2serials using residues and residueNums
		nums2serials = new HashMap<Integer,Integer>();
		for (int i=0;i<residueNums.length;i++){
			// nums2serials must only contain the standard observed aa, so only members of residues TreeMap
			if (residues.containsKey(residueNums[i])){
				nums2serials.put(residueNums[i],i);
			}
		}
		// initialisation of l, numObsStandAA, t, effT from 4 above
		l = this.residueNums.length;
		numObsStandAA = this.residues.size();
		t = (int)((l*(l-1))/2);
		effT = (int)((numObsStandAA*(numObsStandAA-1))/2);
		// initialisation of CM
		CM = new EdgeState[l][l];
		// first we fill with NONCONTACTs (or SKIPPEDs if the i,j is not in nums2serials)
		for (int i=0;i<l;i++){
			for (int j=0;j<l;j++) {
				if ((!this.nums2serials.containsValue(i)) || (!this.nums2serials.containsValue(j))) {
					CM[i][j] = EdgeState.SKIPPED;
				} else {
					CM[i][j] = EdgeState.NONCONTACT;
				}
			}
		}
		// then we fill the CONTACTs
		// the contacts ArrayList we pass as argument must contain the full bidirectional list of contacts
		for (Contact cont:contacts){
			int i_num = cont.i;
			int j_num = cont.j;
			CM[nums2serials.get(i_num)][nums2serials.get(j_num)]=EdgeState.CONTACT;
		}
		// and initialise numContacts
		numContacts = contacts.size();
		// finally we initialise the CNs HashMap
		getAllComNbs();
	}
	
	/**
	 * Gets the number of contacts, non contacts and skipped cells due to unobserved residues 
	 * from the part of the ContactMap >=(above) the given diagonal
	 * @param diagonal
	 * @return
	 */	
	public int[] getCMStats(int diagonal) {
		int[] result = new int[3];
		for (int d=Math.max(1, diagonal);d<l;d++){
			for (int i=0, j=i+d;i<l-d+1 && j<l;i++, j++){
				if (CM[i][j] == EdgeState.CONTACT) result[0]++;
				if (CM[i][j] == EdgeState.NONCONTACT) result[1]++;
				if (CM[i][j] == EdgeState.SKIPPED) result[2]++;
			}
		}
		return result;
	}
	
	/**
	 * Gets the common neighbourhood for residues i_num, j_num from the precomputed HashMap. 
	 * If there is no neighborhood between the two residues it returns an empty TreeMap
	 * @param i_num
	 * @param j_num
	 * @return
	 */
	public TreeMap<Integer,String> getComNbs(int i_num, int j_num) {
		// we initialise to an empty TreeMap, if there is no cn for i_num, j_num we return it empty
		TreeMap<Integer,String> cn = new TreeMap<Integer,String>(); 
		
		if (CNs.containsKey(i_num) && CNs.get(i_num).containsKey(j_num)) {
			cn = CNs.get(i_num).get(j_num);
		} 

		return cn;
	}
	
	/**
	 * Gets the common neighbourhood for residues i_num, j_num from the precomputed HashMap
	 * that are >= (above)/<= (below) the given diagonal
	 * If there is no neighborhood between the two residues it returns an empty TreeMap
	 * @param i_num
	 * @param j_num
	 * @param diagonal
	 * @param above
	 * @return
	 */
	public TreeMap<Integer,String> getComNbs(int i_num, int j_num, int diagonal, boolean above) {
		// we initialise to an empty TreeMap, if there is no cn for i_num, j_num we return it empty
		TreeMap<Integer,String> cn = new TreeMap<Integer,String>(); 
		
		int threshold = diagonal;
		
		if (CNs.containsKey(i_num) && CNs.get(i_num).containsKey(j_num)) {
			cn = CNs.get(i_num).get(j_num);
		}
		
		if (!cn.isEmpty())  {
			ArrayList<Integer> nums2eliminate = new ArrayList<Integer>();
			for (int num:cn.keySet()){
				// here we eliminate all the common neighbours that are below or above the given diagonal
				// 	(first storing the keys in an ArrayList and then removing from TreeMap)
				// num, i_num and j_num are nums (i.e. serial numbers from db)	
				if (above) {
					// just for clarification, the opposite condition is: Math.abs(n-i)>=j-i && Math.abs(n-j)>=j-i
					if ((Math.abs(nums2serials.get(num)-nums2serials.get(i_num)) < threshold || (Math.abs(nums2serials.get(num)-nums2serials.get(j_num)) < threshold))) {
						nums2eliminate.add(num);
					}	
				} else {
					// just for clarification, the opposite condition is: Math.abs(n-i)<=j-i && Math.abs(n-j)<=j-i
					if ((Math.abs(nums2serials.get(num)-nums2serials.get(i_num)) > threshold || (Math.abs(nums2serials.get(num)-nums2serials.get(j_num)) > threshold))) {
						nums2eliminate.add(num);
					}	
				}
			}
			for (int num:nums2eliminate){
				cn.remove(num); //now we remove all the detected non-valid common neighbours
			}
		}

		return cn;
	}
	
	/**
	 * Gets the common neighbourhood for residues i_num, j_num from the precomputed HashMap
	 * that are >= (above)/<= (below) the current diagonal
	 * If there is no neighborhood between the two residues it returns an empty TreeMap
	 * @param i_num
	 * @param j_num
	 * @param above (true if above diagonal)
	 * @return
	 */
	public TreeMap<Integer,String> getComNbs(int i_num, int j_num, boolean above) {
		return getComNbs(i_num, j_num, Math.abs(nums2serials.get(j_num)-nums2serials.get(i_num)), above);
	}
	
	/**
	 * Get all common neighbors for this contact map object. Updates the CNs member and also returns it
	 * @return CNs
	 */
	public HashMap<Integer,HashMap<Integer,TreeMap<Integer,String>>> getAllComNbs() {
		// first we reset the existing CNs
		this.CNs = new HashMap<Integer,HashMap<Integer,TreeMap<Integer,String>>>();
		// initialise a cns4j HashMap which will store all cn TreeMaps for each j 
		HashMap<Integer,TreeMap<Integer,String>> cns4j = new HashMap<Integer,TreeMap<Integer,String>>();
		// initialise a cn TreeMap which will store all cns for a given i,j
		TreeMap<Integer,String> cn = new TreeMap<Integer,String>();
		for (int i=0;i<l;i++){
			for (int j=i+1;j<l;j++) {
				for (int k=0;k<l;k++) { // for each i,j we scan all possible contacts k  
					if ((k != i) && (k != j) && 
						CM[Math.min(i,k)][Math.max(i,k)].contact() && CM[Math.min(j,k)][Math.max(j,k)].contact()) {
						cn.put(residueNums[k],residueTypes[k]);
						// we could have used the following, but would be slower
						//updateCNsHashMap(residueNums[i],residueNums[j],residueNums[k]);
					}
				} // end k
				if (!cn.isEmpty()) {
					cns4j.put(residueNums[j], cn);
				}
				cn = new TreeMap<Integer,String>(); // reset cn for next k
			} // end j
			CNs.put(residueNums[i], cns4j);
			cns4j = new HashMap<Integer,TreeMap<Integer,String>>(); // reset cns4j for next i
		} // end i
		
		return this.CNs;		
	}

	/**
	 * Convenience method to update the monster CNs HashMap without thinking too much
	 * @param i_num
	 * @param j_num
	 * @param k_num
	 */
	public void updateCNsHashMap(int i_num, int j_num, int k_num) {
		if (CNs.containsKey(i_num)){
			HashMap<Integer,TreeMap<Integer,String>> cns4j = CNs.get(i_num);
			if (cns4j.containsKey(j_num)){				
				cns4j.get(j_num).put(k_num, residueTypes[nums2serials.get(k_num)]);
			} else {
				TreeMap<Integer,String> cn = new TreeMap<Integer,String>();
				cn.put(k_num, residueTypes[nums2serials.get(k_num)]);
				cns4j.put(j_num,cn);
			}
		} else {
			HashMap<Integer,TreeMap<Integer,String>> cns4j = new HashMap<Integer,TreeMap<Integer,String>>();
			TreeMap<Integer,String> cn = new TreeMap<Integer,String>();
			cn.put(k_num, residueTypes[nums2serials.get(k_num)]);
			cns4j.put(j_num, cn);
			CNs.put(i_num, cns4j);
		}
	}

	/**
	 * Update all common neighbors given contact (this method should be used only for prediction)
	 * @param i
	 * @param j
	 * @param
	 * @return
	 */
	public void updateComNbsGivenContact(int i, int j, boolean below) {
		
		if (CM[i][j].contact()) { 	// if contact predicted just in case
			
			if (below) {			// if known contacts are only below the diagonal
				int d = Math.abs(i-j);				
				int[] idxs = {i, j};				
				for (int idx:idxs){
					for (int r=Math.max(0,idx-d);r<idx;r++) {
						if ((r != i) && CM[r][idx].contact()) {
							updateCNsHashMap(residueNums[Math.min(i,r)],residueNums[Math.max(i,r)],residueNums[idx]);
						}
					}
					
					for (int c=idx+1;c<=Math.min(l-1,idx+d);c++) {
						if ((c != j) && CM[idx][c].contact()) {
							updateCNsHashMap(residueNums[i],residueNums[c],residueNums[idx]);
						}							
					}
				}
			
				for (int b=i+1; b<j; b++) { // update i-j's common neighbors
					if (CM[i][b].contact() && CM[b][j].contact()) {
						updateCNsHashMap(residueNums[i],residueNums[j],residueNums[b]); 
					}
				}	
			} else { 				// if known contacts can be also above the diagonal			
				for (int k=0;k<l;k++) {
					if ((k != i) && (k != j)) { 
						if (CM[i][j].contact() && CM[Math.min(i,k)][Math.max(i,k)].contact()) {
							updateCNsHashMap(residueNums[Math.min(j,k)],residueNums[Math.max(j,k)],residueNums[i]);
						}
						if (CM[i][j].contact() && CM[Math.min(j,k)][Math.max(j,k)].contact()) {
							updateCNsHashMap(residueNums[Math.min(i,k)],residueNums[Math.max(i,k)],residueNums[j]);
						}
						if (CM[Math.min(i,k)][Math.max(i,k)].contact() && CM[Math.min(j,k)][Math.max(j,k)].contact()) {
							updateCNsHashMap(residueNums[i],residueNums[j],residueNums[k]);
						}
					}
				}			
			}
		}
	}
	
	/**
	 * Returns true if this contact map has no contacts at all
	 * @return
	 */
	public boolean hasNoContacts() {
		boolean noContacts = true;
		for (int i=0;i<l;i++){
			for (int j=i+1;j<l;j++){
				if (this.CM[i][j].contact()) {
					noContacts = false;
					return noContacts;
				}
			}
		}
		return noContacts;
	}
	
	/**
	 * Returns whether the given cell has at least N common neighbors
	 * @param i_num
	 * @param j_num
	 * @param N
	 * @return
	 */
	public boolean cellHasNComNbs(int i_num, int j_num, int N){
		TreeMap<Integer,String> cn = getComNbs(i_num, j_num);
		boolean hasNCNs = false;
		if (((N==0) && cn.isEmpty()) || ((N!=0) && cn.size()>=N)) hasNCNs=true;
		return hasNCNs;
	}
	
	/**
	 * Returns number of contacts in the whole contact map with at least N common neighbors
	 * Useful to find out how many contacts don't have common neighbors at all (N=0)
	 * Will only count the ones above given diagonal value (use diagonal=1 for count in whole contact map)
	 * @param N
	 * @param diagonal
	 * @return
	 */
	public int getNumContactsWithNComNbs(int N, int diagonal){
		int numContactsWithNcns = 0;
		for (int d=Math.max(1, diagonal);d<l;d++){
			for (int i=0, j=i+d;i<l-d+1 && j<l;i++, j++){
				if (CM[i][j].contact()){
					if (cellHasNComNbs(residueNums[i], residueNums[j], N)){
						numContactsWithNcns++;
					}
				}
			}
		}
		return numContactsWithNcns;
	}
	
	/**
	 * Puts into a new ContactMap object the part of this ContactMap that is strictly below the given diagonal (i.e. not including the diagonal itself)
	 * The rest is filled with NONCONTACT
	 * @param diagonal
	 * @return
	 */
	//TODO we fill above range diagonals with NONCONTACT, thus we are "forcing" the return object to be a OrigCM but we want this method for ContactMap's, how to solve this?
	public ContactMap cutCMToBelowRange(int diagonal) {
		ContactMap newcm = this.semiDeepCopy(true); 
		for (int d=1;d<Math.min(diagonal,l);d++){
			for (int i=0, j=i+d;i<l-d+1 && j<l;i++, j++){
				newcm.CM[i][j]=this.CM[i][j];
			}
		}

		for (int d=Math.min(diagonal,l);d<l;d++) {
			for (int i=0, j=i+d;i<l-d+1 && j<l;i++, j++){
				newcm.CM[i][j] = EdgeState.NONCONTACT;
			}
		}
		newcm.getAllComNbs(); // finally we re-calculate the common neighbours for the new contact map
		newcm.updateNumContacts(); // and we recount contacts
		return newcm;
	}
	
	/**
	 * From all contacts in smallerCM, gets a contact map which contains only those contacts that are reachable from this ContactMap
	 * @param smallerCM
	 * @return
	 */
	public ContactMap getReachable(ContactMap smallerCM){
		ContactMap newcm = this.semiDeepCopy(true);
		for (int i=0;i<l;i++) {
			for (int j=i+1;j<l;j++){
				// we loop through all contacts of this contact map
				if (this.CM[i][j].contact()) {
					// if for this contact there's at least 1 cn in smallerCM
					if (smallerCM.cellHasNComNbs(residueNums[i],residueNums[j], 1)){
						//	we assign contact, the rest is kept as it was initialised with semiDeepCopy
						newcm.CM[i][j]=EdgeState.CONTACT;  
					}
				}
			}
		}
		newcm.getAllComNbs(); // finally we re-calculate the CNs for the new ContactMap before returning it
		newcm.updateNumContacts(); // and we recount contacts
		return newcm;
	}
	
	/**
	 * Returns a new ContactMap result of subtracting smallerCM ContactMap from this ContactMap (set subtraction)
	 * @param smallerCM
	 * @return
	 */
	public ContactMap subtract (ContactMap smallerCM){
		ContactMap newcm = this.semiDeepCopy(false); // we copied this cm with all its contacts
		for (int i=0;i<l;i++) {
			for (int j=i+1;j<l;j++){
				// if contact exists in bigger set and also in smaller set then we eliminate it
				if (this.CM[i][j].contact() && smallerCM.CM[i][j].contact()) {
					// we assign NONCONTACT
					newcm.CM[i][j]=EdgeState.NONCONTACT;
				}
			}
		}
		newcm.getAllComNbs(); // finally we re-calculate the CNs for the new ContactMap before returning it
		newcm.updateNumContacts(); // and we recount contacts
		return newcm;
	}

	/**
	 * Returns a new ContactMap result of adding secondCM to this ContactMap (set union) 
	 * @param secondCM
	 * @return
	 */
	public ContactMap add (ContactMap secondCM){
		ContactMap newcm = this.semiDeepCopy(false); // we copied this cm with all its contacts
		for (int i=0;i<l;i++) {
			for (int j=i+1;j<l;j++){
				// if contact doesn't exist in first set, but exists in second set then we add the contact to the resulting set
				if (!this.CM[i][j].contact() && secondCM.CM[i][j].contact()) {
					// we assign CONTACT
					newcm.CM[i][j]=EdgeState.CONTACT;
				}
			}
		}
		newcm.getAllComNbs(); // finally we re-calculate the CNs for the new ContactMap before returning it
		newcm.updateNumContacts(); // and we recount contacts
		return newcm;
	}

	/**
	 * The method performs a semi deep copy of a ContactMap. By semi deep I mean: 
	 * primitives are copied, CM and CNs are copied but residueNums, residueTypes,residues,nums2serials are not copied, but only re-referenced
	 * This is because we don't really need to copy the arrays as they should always stay the same between source and destination
	 * If blanckCM is true then CM and CNs will be just initialised to UNKNOWN in CM and blanks in CNs
	 * @return
	 * @param blankCM
	 */
	public ContactMap semiDeepCopy(boolean blankCM){
		ContactMap destinationCM = new ContactMap();
		// primitives we simply copy
		destinationCM.t = this.t;
		destinationCM.effT = this.effT;
		destinationCM.l = this.l;
		destinationCM.numObsStandAA = this.numObsStandAA;
		destinationCM.numContacts = 0; // we wipe it out initially, it stays so unless blankCM is true where is copied from this.numContacts
		// we don't deep copy here, we can reference the same arrays, they will be always the same for source and destination
		// TODO should we deep copy the arrays too??
		destinationCM.residueNums = this.residueNums;
		destinationCM.residueTypes = this.residueTypes;
		destinationCM.residues = this.residues;
		destinationCM.nums2serials = this.nums2serials;
		// the only thing we need to deep copy is the ContactMap
		destinationCM.CM = new EdgeState[l][l];
		if (!blankCM){
			for (int i=0;i<l;i++) {
				for (int j=i+1;j<l;j++){
					destinationCM.CM[i][j] = this.CM[i][j]; //enum copies should be alright
				}
			}
			destinationCM.numContacts=this.numContacts;
		} else {
			for (int i=0;i<l;i++) {
				for (int j=i+1;j<l;j++){ 
					destinationCM.CM[i][j] = EdgeState.UNKNOWN;
				}
			}			
		}
		// for CNs we create a new object and call the getAllComNbs that gets CNs from the already copied ContactMap array
		destinationCM.CNs = new HashMap<Integer,HashMap<Integer,TreeMap<Integer,String>>>();
		if (!blankCM){
			destinationCM.CNs = destinationCM.getAllComNbs();
		}
		return destinationCM;
	}
	
	/**
	 * Semi deep copy from given otherCM to this ContactMap all data member fields 
	 * @param otherCM
	 * @param blankCM true if we want a blank ContactMap matrix, false if we want to copy from otherCM.CM 
	 */
	public void semiDeepCopyFromOtherCM(ContactMap otherCM, boolean blankCM){
		// primitives we simply copy
		this.t = otherCM.t;
		this.effT = otherCM.effT;
		this.l = otherCM.l;
		this.numObsStandAA = otherCM.l;
		this.numContacts = 0;
		// we don't deep copy here, we can reference the same arrays, they will be always the same for source and destination
		// TODO should we deep copy the arrays too??
		this.residueNums = otherCM.residueNums;
		this.residueTypes = otherCM.residueTypes;
		this.residues = otherCM.residues;
		this.nums2serials = otherCM.nums2serials;
		// the only thing we need to deep copy is the ContactMap
		this.CM = new EdgeState[l][l];
		if (!blankCM){
			for (int i=0;i<l;i++) {
				for (int j=i+1;j<l;j++){
					CM[i][j] = otherCM.CM[i][j]; //enum copies should be alright
				}
			}
			this.numContacts=otherCM.numContacts;
		} else {
			for (int i=0;i<l;i++) {
				for (int j=i+1;j<l;j++){ 
					CM[i][j] = EdgeState.UNKNOWN;
				}
			}			
		}
		// for CNs we create a new object and call the getAllComNbs that gets CNs from the already copied ContactMap array
		this.CNs = new HashMap<Integer,HashMap<Integer,TreeMap<Integer,String>>>();
		if (!blankCM){
			this.CNs = this.getAllComNbs();
		}

	}

	/**
	 * To update the numContacts data member, to be used after a ContactMap matrix has been changed in a ContactMap object 
	 *
	 */
	public void updateNumContacts(){
		this.numContacts=0;
		for (int i=0;i<l;i++){
			for (int j=i+1;j<l;j++) {
				if (this.CM[i][j].contact()) this.numContacts++;
			}
		}	
	}
	
	/**
	 * Returns true if this ContactMap and otherCM have exactly the same set of contacts (true/false, disregarding other EdgeStates)
	 * @param otherCM
	 * @return
	 */
	public boolean hasSameContacts (ContactMap otherCM){
		boolean coincides = true;
		for (int i=0;i<l;i++){
			for (int j=i+1;j<l;j++) {
				if ((this.CM[i][j].contact() && !otherCM.CM[i][j].contact()) 
						|| !this.CM[i][j].contact() && otherCM.CM[i][j].contact()) {
					coincides = false;
					return coincides;
				}
			}
		}	 
		return coincides;
	}
	
	/**
	 * To print constraint equations for a contact
	 * @param i_num
	 * @param j_num
	 */
	public void printConstraint(int i_num,int j_num){
		String xi = "x"+i_num;
		String yi = "y"+i_num;
		String zi = "z"+i_num;
		String xj = "x"+j_num;
		String yj = "y"+j_num;
		String zj = "z"+j_num;
		double d = 4.1;
		double d2 = d*d;
		System.out.println("("+xi+"-"+xj+")^2+("+yi+"-"+yj+")^2+("+zi+"-"+zj+")^2<="+d2);
	}
}
