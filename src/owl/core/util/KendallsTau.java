package owl.core.util;

public class KendallsTau {
	

	private double values[];
	private int ranking[];
	private double freq[];
	private double total = 0; 
	private String keys[];
//	private String VStr = ""; 
	
	// ------ Constructors ------------------------
	public KendallsTau(double[] values) { // Constructor 
		this.values = values;
		initialiseKeys();
		initialiseRank();
		this.total = this.values.length;
		recalc();
	}
	
	public KendallsTau(double[] values, int[] rank) { // Constructor 
		if (values.length == rank.length){
			this.values = values;
			initialiseKeys();
			this.ranking = rank;
			this.total = this.values.length;
			recalc();			
		}
		else{	
			System.out.println("Invalid input arguments: arrays should be of same length!");
		}
	}
	
	public KendallsTau(double[] values, String[] keys) { // Constructor 
		this.values = values;
		this.keys = keys;
		initialiseRank();
		this.total = this.values.length;
		recalc();
	}

	public KendallsTau(double[] values, String[] keys, int[] rank) { // Constructor 
		if (values.length == keys.length && values.length == rank.length){
			this.values = values;
			this.keys = keys;
			this.ranking = rank;
			this.total = this.values.length;
			recalc();		
		}
		else{	
			System.out.println("Invalid input arguments: arrays should be of same length!");
		}
	}
	
	// ------ private methods ------------------------
	
	private void initialiseRank() {
		// TODO Auto-generated method stub
		this.ranking = new int[this.values.length];
		// rank with respect to increasing value
		// rank = 1 --> lowest value
		for (int i=0; i<this.values.length; i++){
			int rank = 0;
			for (int j=0; j<this.values.length; j++){
				if (this.values[j] < this.values[i] && i!=j)
					rank++;
			}			
			rank++;
			this.ranking[i] = rank;
		}		
	}
	
	private void initialiseKeys() {
		this.keys = new String[this.values.length];
		for (int i=0; i<this.values.length; i++){
			this.keys[i] = String.valueOf(i+1);
		}
	}
	
	private void recalc(){		
		this.freq = new double[this.values.length];
		total=0; 
		for (int i=0; i<this.values.length; i++) 
			total+=values[i]; // calculate the new total 
		if (total>0) for (int i=0; i<freq.length; i++) 
			freq[i]=(double)(values[i]/total); // recalc the frequency of each aa w.r.t. the new total
	}
	
	
	// ------ public methods (getters and setters)------------------------
	
	public double[] getValues(){
		return this.values;
	}
	public String[] getKeys(){
		return this.keys;
	}
	public int[] getRanking(){
		return this.ranking;
	}
	
	public void reSetCounts( )  { // recalculation of total and frequencies 
		total=0; 
		for (int i=0; i<values.length; i++) { 
			values[i]=0; 
			freq[i]=(double)(0.0); 
		}
	} 
	
	public void setValue(String key, double val)  {
		int index = getIndexOfKey(key);
		setValue(index, val);
	}
	private void setValue(int index, double val)  {
		total = total-values[index];
		values[index] = val; 
		total = total+val; 
		freq[index] = (double)(val/total); 
	}
	
	private int getIndexOfRank(int r){
		for (int i=0; i<this.ranking.length; i++) {
			if( ranking[i]==r) 
				return i;
		} // next i 
		System.out.println("Index will be out of bounds of array ranking");
		return -1; 
	}
	private int getIndexOfKey(String key){
		for (int i=0; i<this.keys.length; i++) {
			if( this.keys[i].equals(key)) 
				return i;
		} // next i 
		System.out.println("Index will be out of bounds of array keys");
		return -1; 
	}
	
	public String getKeyByRank(int r) {
		for (int i=0; i<this.ranking.length; i++) {
			if( ranking[i]==r) 
				return this.keys[i];
		} // next i 
		return "?"; 
	}
	public int getRankByKey(String key){
		int index = getIndexOfKey(key);
		return this.ranking[index];
	}
	public double getValByRank(int r) {
		for (int i=0; i<this.ranking.length; i++) {
			if( ranking[i]==r) 
				return this.values[i];
		} // next i 
		return -666666; 
	} 
	public double getValueByKey(String key){
		int index = getIndexOfKey(key);
		return getValue(index);
	}
	private double getValue(int index) {
		if (index<0 || index >=this.values.length){
			System.out.println("Index "+index+" out of bounds of array values with length "+values.length);
		}
		return this.values[index];
	}
	
	public void setRankAtKey(String key, int rank){
		int index = getIndexOfKey(key);
		this.ranking[index]=rank;
	}
	private void setRankAtIndex(int index, int rank){
		this.ranking[index]=rank;
	}
		
	// ------ public methods (Kendalls tau calculation)------------------------
	
	// calcrank = bubblesort = pairwise exchange of rank labels under condition val_i>val_j 
	public int bubbleSort( )  {
		int i, nrOfXChanges=0, deltaX=-1; 
		double c, d;
//		@SuppressWarnings("unused")
//		String a, b; 
		int indexA, indexB;
		
		while (deltaX != 0) {
			deltaX = 0; 
			for (i=1; i<this.ranking.length; i++) {
				indexA = getIndexOfRank(i);
				indexB = getIndexOfRank(i+1);
//				a = this.keys[indexA]; //a = getKeyByRank(i);
//				b = this.keys[indexB]; //getKeyByRank(i+1);
				c = this.values[indexA]; //getValueByKey(a);
				d = this.values[indexB]; //getValueByKey(b); 
				if (c!=-666666 && d!=-666666 && indexA!=-1 &&indexB!=-1){
					if( c < d ) { // discordant pair -> exchange 
						setRankAtIndex(indexA, i+1); //setRankAtKey(a, i+1);
						setRankAtIndex(indexB, i); //setRankAtKey(b, i+);
						nrOfXChanges++; 
						deltaX++; 
					} else { } // end if discordant 
				}
				else
					System.out.println("Rank "+i+" is not assigned to any value!");
			} // next count comparison 
		} // end while xchanges happen 
		
//		for (i=0; i<values.length; i++)
//			System.out.println(ranking[i]+" "+values[i]);
		
		return nrOfXChanges; 
	} // end of bubblesort 
	
	// similarly, compare 2 ResultVectors to yield Kendall's tau 
	public double kendallsTau( KendallsTau cmp)  { // by just copying the ranks and resorting
		for (int i=0; i<keys.length; i++) 
			setRankAtKey(keys[i], cmp.getRankByKey(keys[i]));
		int nxc=bubbleSort();
		double fac = this.values.length*(this.values.length-1)/2;
		return (double)((double)nxc/(double)fac); 
	} // end of kendallsTau 

	// kendalls tau B implementation taking care of ties 
	public double kendallsTauB( KendallsTau cmp)  {
		int i, j, C=0, D=0, T=0;
		double a1, a2, b1, b2;
		double tauB=(double)0.0; 
		String keyA, keyB;
		for (i=0; i<keys.length; i++) { // retrieving all aa-pairs <i,j> : (a,b) 
			keyA = keys[i];
			for (j=i+1; j<keys.length; j++) {			
				keyB = keys[j]; 
//				System.out.print("\n<"+i+","+j+"> \t ["+a+","+b+"] ");
				a1 = this.values[i]; // this.getCount( a); 
				b1 = this.values[j]; // this.getCount( b); 
				a2 = cmp.getValueByKey(keyA); 
				b2 = cmp.getValueByKey(keyB); 
				if (a1==b1 || a2==b2) { 
					// we have a tie on one or both sides 
					T++; 
				} else {
					if( (a1 < b1 && a2 < b2) || (a1 > b1 && a2 > b2) ) { // concordance 
						C++;
					} else { // discordant pair  
						D++;
					} // end if (dis-)cordant pairs 
				} // end if no tie 
			} // next j  				
		} // next i 
		
		if ( ( C+D ) > 0) 
			tauB = (double) ( ((double)(C-D)) / ((double)(C+D)) ); 
		else tauB=(double)0.0;
		
		return tauB; 
	} // end of kendallsTauB 

	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		// Test method
		int ranking[] = new int[5]; //[]{1,2,3,4,5,6,7,8,9,10};
		for (int i=0; i<ranking.length; i++)
			ranking[i] = i+1;
		double values1[] = new double[]{1.1, 1.3, 1.5, 1.6, 1.7};
		double values2[] = new double[]{1.8, 1.7, 1.4, 1.2, 1.0};		
//		String keys[] = new String[]{};
		double tB = 0;
		
		if (ranking.length==values1.length && ranking.length==values2.length) // && ranking.length==keys.length)
			System.out.println("Correct Initialization");
		else
			System.out.println("Wrong Initialization");
		
//		for (int i=0; i<values1.length; i++)
//			System.out.println(ranking[i]+" "+values1[i]+" "+values2[i]);		
		
		KendallsTau kt1 = new KendallsTau(values1, ranking);
		KendallsTau kt2 = new KendallsTau(values2, ranking);
		for (int i=0; i<values1.length; i++)
			System.out.println(kt1.getRanking()[i]+" "+kt1.getValues()[i]+"   "+kt2.getRanking()[i]+" "+kt2.getValues()[i]);
		int xChanges = kt1.bubbleSort();
		System.out.println("nr of Exchanges total kt1: "+xChanges);
		xChanges = kt2.bubbleSort();
		System.out.println("nr of Exchanges total kt2: "+xChanges);
		for (int i=0; i<values1.length; i++)
			System.out.println(kt1.getRanking()[i]+" "+kt1.getValues()[i]+"   "+kt2.getRanking()[i]+" "+kt2.getValues()[i]);
		
//		// compare the two RVs 
//		double ktVal = kt1.kendallsTau(kt2);
//		System.out.println("normalized Kendalls tau = (x*2)/((n-1)*n) = "+ktVal); 
//		for (int i=0; i<values1.length; i++)
//			System.out.println(kt1.getRanking()[i]+" "+kt1.getValues()[i]+"   "+kt2.getRanking()[i]+" "+kt2.getValues()[i]);
				
		tB= kt1.kendallsTauB( kt2);
		System.out.println("calculating Kendalls tau B : "+tB); 
		for (int i=0; i<values1.length; i++)
			System.out.println(kt1.getRanking()[i]+" "+kt1.getValues()[i]+"   "+kt2.getRanking()[i]+" "+kt2.getValues()[i]);
//		System.out.println("nr of Exchanges total kt2: "+kt2.bubbleSort());
//		for (int i=0; i<values1.length; i++)
//			System.out.println(kt1.getRanking()[i]+" "+kt1.getValues()[i]+"   "+kt2.getRanking()[i]+" "+kt2.getValues()[i]);
//		tB= kt1.kendallsTauB( kt2);
//		System.out.println("calculating Kendalls tau B : "+tB); 
	}

}
