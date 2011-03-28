package owl.core.sequence;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import owl.core.structure.AminoAcid;

/**
 * Class containing static methods to translate DNA nucleotide sequences
 * to protein sequences.
 * Using the translation tables from NCBI found at:
 * http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
 * 
 * @author duarte_j
 *
 */
public class Translator {

	private static final String RESOURCES_DIR = "/owl/core/sequence/"; // path to data files to be retrieved with getResourceAsStream
	private static final String GENCODE_FILE = "geneticCode.dat";
	
	private static final Map<Integer,Map<Codon,AminoAcid>> gcmaps = parseNCBIGenCodeFile(); // map of ncbi's genetic code identifiers to maps of codons to aminoacids
	
	/**
	 * Parses NCBI genetic codes file
	 * @return
	 */
	private static Map<Integer,Map<Codon,AminoAcid>> parseNCBIGenCodeFile() {
		
		Map<Integer,Map<Codon,AminoAcid>> allmaps= new HashMap<Integer,Map<Codon,AminoAcid>>();
		
		InputStream inp = null;
		inp = Translator.class.getResourceAsStream(RESOURCES_DIR+GENCODE_FILE);
		if(inp == null) {
			System.err.println("Fatal Error! Resource " + RESOURCES_DIR+GENCODE_FILE + " not found. Could not initialize genetic code maps. Exiting.");
			System.exit(1);
		}

		
		
		try {
			BufferedReader br = new BufferedReader(new InputStreamReader(inp));
			String line;
			Pattern p = Pattern.compile("\\s+id\\s+(\\d+).*");
			Pattern p1 = Pattern.compile("\\s+ncbieaa\\s+\"(.*)\".*");
			Pattern p2 = Pattern.compile("\\s+--\\sBase1\\s+([^\\s]+).*");
			Pattern p3 = Pattern.compile("\\s+--\\sBase2\\s+([^\\s]+).*");
			Pattern p4 = Pattern.compile("\\s+--\\sBase3\\s+([^\\s]+).*");

			int id = -1;
			String ncbieaa = null;
			String base1 = null, base2 = null, base3 = null;

			while((line=br.readLine())!=null){
				if (line.startsWith("--")) continue;
				if (line.trim().isEmpty()) continue;
				Matcher m =p.matcher(line);
				Matcher m1 =p1.matcher(line);
				Matcher m2 =p2.matcher(line);
				Matcher m3 =p3.matcher(line);
				Matcher m4 =p4.matcher(line);
				if (m.matches()){
					id = Integer.parseInt(m.group(1));
				} else if (m1.matches()) {
					ncbieaa = m1.group(1);
				} else if (m2.matches()) {
					base1 = m2.group(1);
				} else if (m3.matches()) {
					base2 = m3.group(1);
				} else if (m4.matches()) {
					base3 = m4.group(1);
					allmaps.put(id,initMap(ncbieaa,base1,base2,base3));
					id = -1;
					ncbieaa = null;
					base1 =null;
					base2 = null;
					base3 = null;
				}

			}
			br.close();
		} catch (IOException e) {
			System.err.println("Fatal error! Could not read the genetic code from resource file "+RESOURCES_DIR+GENCODE_FILE);
			System.err.println(e.getMessage());
			System.exit(1);
		}
		
		return allmaps;
	}
	
	private static Map<Codon,AminoAcid> initMap(String ncbieaa, String base1, String base2, String base3) {
		Map<Codon,AminoAcid> map = new HashMap<Codon, AminoAcid>();
		for (int i=0;i<ncbieaa.length();i++) {
			AminoAcid aa = AminoAcid.getByOneLetterCode(ncbieaa.charAt(i));
			char b1 = Character.toLowerCase(base1.charAt(i));
			char b2 = Character.toLowerCase(base2.charAt(i));
			char b3 = Character.toLowerCase(base3.charAt(i));
			map.put(new Codon(b1,b2,b3),aa);
		}
		return map;
	}
	
	/**
	 * Translate given codon to its corresponding AminoAcid based on given GeneticCodeType
	 * @param gcType
	 * @param codon
	 * @return
	 * @throws IllegalArgumentException if codon length is not 3
	 * @throws TranslationException if an ambiguous nucleotide code is found in the codon
	 */
	public static AminoAcid translate(GeneticCodeType gcType, String codon) throws TranslationException {
		return translate(gcType,new Codon(codon.toLowerCase()));
	}

	/**
	 * Translate given codon to its corresponding AminoAcid based on given GeneticCodeType
	 * @param gcType
	 * @param codon
	 * @return
	 * @throws TranslationException if an ambiguous nucleotide code is found in the codon
	 */
	public static AminoAcid translate(GeneticCodeType gcType, Codon codon) throws TranslationException {
		if (codon.isAmbiguous()) {
			throw new TranslationException("Ambiguous nucleotide in codon '"+codon+"'");
		}
		return gcmaps.get(gcType.getId()).get(codon);		
	}
	
	/**
	 * Translate given DNA nucleotide sequence on given reading frame to a protein sequence 
	 * based on given GeneticCodeType 
	 * @param gcType
	 * @param sequence
	 * @param rf
	 * @return
	 */
	public static Sequence translate(GeneticCodeType gcType, Sequence sequence, ReadingFrame rf) throws TranslationException {
		String seq = sequence.getSeq();
		StringBuffer sb = new StringBuffer();
		
		String scanningSeq = seq;
		if (rf.isReverse()) {
			scanningSeq = new StringBuffer(seq).reverse().toString();
		}
		int startIndex = Math.abs(rf.getNumber())-1;
		for (int i=startIndex;i<seq.length() && i+3<=seq.length();i+=3){
			sb.append(translate(gcType,scanningSeq.substring(i,i+3)).getOneLetterCode());
		}
		Sequence translated = new Sequence(sequence.getName(),sb.toString());
		return translated;
	}
}
