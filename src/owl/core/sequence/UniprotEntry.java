package owl.core.sequence;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import owl.core.connections.EmblWSDBfetchConnection;
import owl.core.connections.NoMatchFoundException;
import owl.core.connections.UniProtConnection;

import uk.ac.ebi.kraken.interfaces.uniprot.DatabaseType;
import uk.ac.ebi.kraken.interfaces.uniprot.GeneEncodingType;
import uk.ac.ebi.kraken.interfaces.uniprot.NcbiTaxonomyId;
import uk.ac.ebi.kraken.interfaces.uniprot.Organelle;
import uk.ac.ebi.kraken.interfaces.uniprot.UniProtEntry;
import uk.ac.ebi.kraken.interfaces.uniprot.dbx.embl.Embl;

/**
 * A uniprot entry.
 * The class stores the taxonomy id data and the EMBL CDS ids (coding sequence) and sequences
 *  
 * @author duarte_j
 *
 */
public class UniprotEntry {

	private String uniId;
	private Sequence uniprotSeq;
	private List<String> taxIds; // many taxonomy ids for a uniprot entry should be a very pathological case. Should be safe to use the first one always
	private List<String> emblCdsIds;
	private List<Sequence> emblCdsSeqs; 
	private GeneEncodingType geneEncodingOrganelle; // the organelle where this gene is encoded (important for genetic code): if gene encoded in nucleus this is null
	
	private boolean  repCDScached;
	private Sequence representativeCDS; // cached after running getRepresentativeCDS()

	
	public UniprotEntry(String uniprotId){
		this.uniId = uniprotId;
	}
	
	/**
	 * @return the uniId
	 */
	public String getUniId() {
		return uniId;
	}

	/**
	 * @param uniId the uniId to set
	 */
	public void setUniId(String uniId) {
		this.uniId = uniId;
	}

	/**
	 * @return the uniprotSeq
	 */
	public Sequence getUniprotSeq() {
		return uniprotSeq;
	}

	/**
	 * @param uniprotSeq the uniprotSeq to set
	 */
	public void setUniprotSeq(Sequence uniprotSeq) {
		this.uniprotSeq = uniprotSeq;
	}

	/**
	 * @return the taxIds
	 */
	public List<String> getTaxIds() {
		return taxIds;
	}

	/**
	 * @param taxIds the taxIds to set
	 */
	public void setTaxIds(List<String> taxIds) {
		this.taxIds = taxIds;
	}

	/**
	 * @return the emblCdsIds
	 */
	public List<String> getEmblCdsIds() {
		return emblCdsIds;
	}

	/**
	 * @param emblCdsIds the emblCdsIds to set
	 */
	public void setEmblCdsIds(List<String> emblCdsIds) {
		this.emblCdsIds = emblCdsIds;
	}

	/**
	 * @return the emblCdsSeqs
	 */
	public List<Sequence> getEmblCdsSeqs() {
		return emblCdsSeqs;
	}

	/**
	 * @param emblCdsSeqs the emblCdsSeqs to set
	 */
	public void setEmblCdsSeqs(List<Sequence> emblCdsSeqs) {
		this.emblCdsSeqs = emblCdsSeqs;
	}

	/**
	 * @return the geneEncodingOrganelle
	 */
	public GeneEncodingType getGeneEncodingOrganelle() {
		return geneEncodingOrganelle;
	}

	/**
	 * @param geneEncodingOrganelle the geneEncodingOrganelle to set
	 */
	public void setGeneEncodingOrganelle(GeneEncodingType geneEncodingOrganelle) {
		this.geneEncodingOrganelle = geneEncodingOrganelle;
	}

	
	/**
	 * Retrieves from UniprotKB the sequence, taxonomy and EMBL CDS ids data,
	 * by using the remote Uniprot API
	 */
	public void retrieveUniprotKBData() {
		this.taxIds = new ArrayList<String>();
		this.emblCdsIds = new ArrayList<String>();
		Set<String> tmpEmblCdsIdsSet = new TreeSet<String>();
		
		UniProtConnection uniprotConn = new UniProtConnection();
		UniProtEntry entry = null;
		try {
			entry = uniprotConn.getEntry(uniId);
		} catch (NoMatchFoundException e) {
			System.err.println("Warning: couldn't find uniprot id "+uniId+" through Uniprot JAPI");
			return;
		}
		
		this.setUniprotSeq(new Sequence(this.getUniId(),entry.getSequence().getValue()));
		
		for(NcbiTaxonomyId ncbiTaxId:entry.getNcbiTaxonomyIds()) {
			taxIds.add(ncbiTaxId.getValue());
		}
		
		Collection<Embl> emblrefs = entry.getDatabaseCrossReferences(DatabaseType.EMBL);
		for(Embl ref:emblrefs) {
			String emblCdsIdWithVer = ref.getEmblProteinId().getValue();
			if (!emblCdsIdWithVer.equals("-")) { // for non annotated genomic dna cds sequences the identifier is '-', we ignore them
				String emblCdsId = emblCdsIdWithVer.substring(0, emblCdsIdWithVer.lastIndexOf("."));
				//emblCdsIds.add(emblCdsId);
				// to ensure there are no duplicate ids (it happens sometimes) we use a set
				tmpEmblCdsIdsSet.add(emblCdsId);
			}
		}
		this.emblCdsIds.addAll(tmpEmblCdsIdsSet);
		List<Organelle> orglls = entry.getOrganelles();
		if (orglls.size()>0) {
			this.geneEncodingOrganelle = orglls.get(0).getType();
			if (orglls.size()>1) {
				for (Organelle orgll:orglls){ 
					if (!orgll.getType().equals(this.geneEncodingOrganelle)) {
						System.err.println("Warning! Different gene encoding organelles for Uniprot "+this.uniId);
					}
				}
			}
		}
	}
	
	/**
	 * Retrieves from EMBL DB fetch web service the EMBL CDS sequences
	 * @param cacheFile a FASTA file containing the sequences to retrieve. If present and if
	 * it contains the required sequence then it is read from cacheFile. If null or file
	 * does not exist or file older than {@link EmblWSDBfetchConnection.MAX_CACHE_AGE} then the sequences are 
	 * retrieved from EMBL DB fetch 
	 * @throws IOException
	 */
	public void retrieveEmblCdsSeqs(File cacheFile) throws IOException {
		this.emblCdsSeqs = new ArrayList<Sequence>();
		
		try {
			this.emblCdsSeqs = EmblWSDBfetchConnection.fetchEMBLCDS(this.emblCdsIds, cacheFile);
		} catch (NoMatchFoundException e) {
			// this is unlikely to happen here, that's why we don't write a better error message
			System.err.println("Couldn't retrieve EMBL CDS sequences for EMBL cds ids"); 
		}

	}

	/**
	 * Returns true if there is at least one associated CDS for which we could 
	 * retrieve a sequence
	 * @return
	 */
	public boolean hasCDS() {
		return (this.emblCdsSeqs!=null && emblCdsSeqs.size()>0); 
	}
	
	/**
	 * Gets a single representative translated coding sequence from the set of CDS 
	 * retrieved from EMBL. The sequence found is then cached and subsequent calls
	 * to this method take the sequence from the cached value.
	 * The CDSs are translated on their 6 frames, then aligned to the protein sequence
	 * and if 100% matches exist one of them is returned. If no 100% matches exist the 
	 * best one is returned (TODO must still check what are tolerable mismatches)
	 * @return the translated CDS sequence or null if no correct match to the protein
	 * sequence can be found
	 */
	public Sequence getRepresentativeCDS() {
		if (repCDScached) {
			return representativeCDS;
		}
		
		if (this.geneEncodingOrganelle!=null) {
			System.err.println("Warning! The entry "+this.getUniId()+" is not encoded in nucleus!");
		}

		List<Sequence> fullMatchingCDSs = new ArrayList<Sequence>();
		List<TranslatedSequence> nonFullMatchingCDSs = new ArrayList<TranslatedSequence>();
		for (Sequence cds:emblCdsSeqs) {
			
			ProteinToCDSMatch matching = null;
			try {
				matching = new ProteinToCDSMatch(this.getUniprotSeq(), cds, GeneticCodeType.STANDARD);
			} catch (TranslationException e) {
				System.err.println("Couldn't translate embl CDS "+cds.getName());
				System.err.println(e.getMessage());
				continue; // try the next one
			}
			if (matching.hasFullMatch()) {
				fullMatchingCDSs.add(cds);
			} else {
				nonFullMatchingCDSs.add(matching.getBestTranslation()); 
			}
		}
		if (fullMatchingCDSs.size()>1) {
			HashMap<String,Sequence> map = new HashMap<String, Sequence>();
			for (Sequence cds:fullMatchingCDSs) {
				if (!map.containsKey(cds.getSeq())) {
					map.put(cds.getSeq(), cds); // eliminating identical cdss
				}
			}
			
			Sequence rep = map.values().iterator().next();
			if (map.size()>1) {
				System.err.println("Warning! Multiple fully matching CDSs for uniprot entry "+this.getUniId()+". Using first CDS only ("+rep.getSecondaryAccession()+")"); 
			}
			representativeCDS = rep; // returning the first value (the one value if map.size()==1)
			
		} else if (fullMatchingCDSs.size()==1) {
			representativeCDS = fullMatchingCDSs.get(0);
		// NO full matching CDS, we try with the best matches
		} else if (nonFullMatchingCDSs.size()>=1) {
			TranslatedSequence best = Collections.max(nonFullMatchingCDSs);
			int mismatches = best.getAln().getLength()-best.getAln().getIdentity();
			System.err.println("Warning! No perfect CDS matches for uniprot entry "+this.getUniId()+". Using the best match with "+mismatches+" mismatches.");
			best.getAln().printAlignment();
			//TODO decide what are tolerable mismatching sequences. If tolerable return it, if not return null
			representativeCDS = best.getSequence();
		} else {
			// no good CDS found
			representativeCDS = null;
		}
		repCDScached = true;
		return representativeCDS;
	}
	
	public void checkEmblCDSMatching() {
		if (this.geneEncodingOrganelle!=null) {
			System.err.println("Warning! The entry "+this.getUniId()+" is not encoded in nucleus!");
		}
		int countFullMatches = 0;
		for (Sequence cds:emblCdsSeqs) {
			
			ProteinToCDSMatch matching = null;
			try {
				matching = new ProteinToCDSMatch(this.getUniprotSeq(), cds, GeneticCodeType.STANDARD);
			} catch (TranslationException e) {
				System.err.println("Couldn't translate embl CDS "+cds.getName());
				System.err.println(e.getMessage());
				continue; // try the next one
			}
			matching.printSummary(99.9999f);
			if (matching.hasFullMatch()) {
				countFullMatches++;
			}
		}
		if (countFullMatches>1) {
			System.out.println(this.getUniId()+": "+countFullMatches+" perfect CDS matches");
		}
	}


}
