package proteinstructure;

import java.util.TreeMap;
import java.util.Collections;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class NodeNbh extends TreeMap<Integer,String> {
	

	private static final long serialVersionUID = 1L;
	
	public static final String centralLetter="x";

	// central residue
	public int central_resser;
	public String central_resType;
	
	/**
	 * Specific NodeNbh: is a neighbourhood in a specific structure e.g. ABCxEFG with x=D in position 25
	 * @param resser
	 * @param resType
	 */
	public NodeNbh(int resser, String resType){
		super();
		this.central_resser=resser;
		this.central_resType=resType;
	}

	public String getMotifFullGaps(){
		String motif="";
		int min=Math.min(central_resser, Collections.min(this.keySet()));
		int max=Math.max(central_resser, Collections.max(this.keySet()));
		for (int i=min;i<=max;i++){
			if (this.containsKey(i)){
				motif+=AA.threeletter2oneletter(this.get(i));
			} else if (i==central_resser){
				motif+=centralLetter;
			} else {
				motif+="_";
			}
		}
		return motif;
	}
	
	public String getMotif(){
		String motif="";
		int min=Math.min(central_resser, Collections.min(this.keySet()));
		int max=Math.max(central_resser, Collections.max(this.keySet()));
		int gapSize=0;
		String gap="";
		for (int i=min;i<=max;i++){
			if (this.containsKey(i)){
				motif+=gap;
				motif+=AA.threeletter2oneletter(this.get(i));
				gapSize=0;
				gap="";
			} else if (i==central_resser){
				motif+=gap;
				motif+=centralLetter;
				gapSize=0;
				gap="";
			} else {
				gapSize++;
				gap="_{"+gapSize+"}";
			}
		}
		return motif;
	}

	// **************************** static methods ****************************** //

	public static String motif2motifFullGaps(String motif){
		String newmotif="";
		if (motif.contains("_")){
			String[] tokens = motif.split("_\\{(\\d+)\\}");
			String[] gaps = new String[tokens.length-1];
			Pattern p = Pattern.compile("_\\{(\\d+)\\}");
			Matcher m = p.matcher(motif);
			int index=0;
			while (m.find()){
				int gapLength = Integer.parseInt(m.group(1));
				String gap = "";
				for (int i=1;i<=gapLength;i++){
					gap+="_";
				}
				gaps[index]=gap;
				index++;
			}
			for (int i=0;i<tokens.length;i++){
				newmotif+=tokens[i];
				if (i!=tokens.length-1){ // there are tokens.length-1 gaps, so we don't want the last i
					newmotif+=gaps[i];
				}
			}
		} else {
			newmotif=motif;
		}		
		return newmotif;
	}

	public static String motif2motifNoGapLength(String motif){
		String newmotif="";
		if (motif.contains("_")){
			newmotif=motif.replaceAll("_\\{\\d+\\}", "_");
		} else {
			newmotif=motif;
		}		
		return newmotif;
	}

	public static String motif2motifNoGaps(String motif){
		String newmotif="";
		if (motif.contains("_")){
			newmotif=motif.replaceAll("_\\{\\d+\\}", "");
		} else {
			newmotif=motif;
		}		
		return newmotif;
	}

	public static String motif2regexStr(String motif){
		String motifng = motif2motifNoGaps(motif);
		String regex = "";
		for (char letter:motifng.toCharArray()){
			regex+=letter;
			regex+=".*";
		}
		return regex;
	}
	
	public static boolean matchMotifs(String bgmotif, String querymotif, int compareType){
		// -- types:
		// 1. substring exact match
		// 2. substring match ignoring gap length
		// 3. substring match ignoring gaps
		// 4. wildcard padded substring match
		switch (compareType) {
		case 1: 
			if (bgmotif.contains(querymotif)) return true;
		case 2:
			if (motif2motifNoGapLength(bgmotif).contains(motif2motifNoGapLength(querymotif))) return true;
		case 3:
			if (motif2motifNoGaps(bgmotif).contains(motif2motifNoGaps(querymotif))) return true;
		case 4:
			if (bgmotif.matches(motif2regexStr(querymotif))) return true;
		}
		return false;
	}
	
}
