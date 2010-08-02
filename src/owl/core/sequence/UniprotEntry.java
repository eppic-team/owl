package owl.core.sequence;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import owl.core.connections.EmblWSDBfetchConnection;
import owl.core.connections.NoMatchFoundException;
import owl.core.connections.UniProtConnection;
import owl.core.features.Feature;
import owl.core.features.FeatureType;
import owl.core.features.HasFeatures;
import owl.core.features.InvalidFeatureCoordinatesException;
import owl.core.features.OverlappingFeatureException;
import owl.core.util.Interval;
import owl.core.util.IntervalSet;

import uk.ac.ebi.kraken.interfaces.uniprot.DatabaseType;
import uk.ac.ebi.kraken.interfaces.uniprot.NcbiTaxon;
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
public class UniprotEntry implements HasFeatures {
	
	private static final Logger LOGGER = Logger.getLogger(UniprotEntry.class);

	private static final float MIN_TOLERATED_ID = 0.95f; // the minimum sequence identity for the CDS sequence to be considered as a representative CDS of the uniprot entry
	
	private static final String GENE_ENC_ORG_PLASMID = "Plasmid";
	
	private String uniId;
	private Sequence uniprotSeq;
	private String taxId;
	private List<String> taxons; // the taxonomy as returned by UniprotEntry.getTaxonomy(). First value is the most general (kingdom), last is the most specific (species).
	private List<String> emblCdsIds;
	private List<Sequence> emblCdsSeqs; 
	private String geneEncodingOrganelle; // the organelle where this gene is encoded (important for genetic code): if gene encoded in nucleus this is null
	
	private boolean  repCDScached;
	private ProteinToCDSMatch representativeCDS; // cached after running getRepresentativeCDS()

	private Map<FeatureType, Collection<Feature>> features;
	
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
	 * @return the taxId
	 */
	public String getTaxId() {
		return taxId;
	}

	/**
	 * @param taxId the taxId to set
	 */
	public void setTaxId(String taxId) {
		this.taxId = taxId;
	}
	
	public List<String> getTaxons() {
		return taxons;
	}
	
	public String getFirstTaxon() {
		return this.taxons.get(0);
	}
	
	public String getLastTaxon() {
		return this.taxons.get(this.taxons.size()-1);
	}
	
	public void setTaxons(List<String> taxons) {
		this.taxons = taxons;
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
	public String getGeneEncodingOrganelle() {
		return geneEncodingOrganelle;
	}

	/**
	 * @param geneEncodingOrganelle the geneEncodingOrganelle to set
	 */
	public void setGeneEncodingOrganelle(String geneEncodingOrganelle) {
		this.geneEncodingOrganelle = geneEncodingOrganelle;
	}

	/**
	 * Return the sequence length of this uniprot entry.
	 * @return
	 */
	public int getLength() {
		return this.getUniprotSeq().getLength();
	}
	
	/**
	 * Retrieves from UniprotKB the sequence, taxonomy and EMBL CDS ids data,
	 * by using the remote Uniprot API
	 */
	public void retrieveUniprotKBData() {
		this.taxons = new ArrayList<String>();
		this.emblCdsIds = new ArrayList<String>();
		Set<String> tmpEmblCdsIdsSet = new TreeSet<String>();
		
		UniProtConnection uniprotConn = new UniProtConnection();
		UniProtEntry entry = null;
		try {
			entry = uniprotConn.getEntry(uniId);
		} catch (NoMatchFoundException e) {
			LOGGER.warn("Couldn't find uniprot id "+uniId+" through Uniprot JAPI");
			return;
		}
		
		this.setUniprotSeq(new Sequence(this.getUniId(),entry.getSequence().getValue()));
		
		List<NcbiTaxonomyId> ncbiTaxIds = entry.getNcbiTaxonomyIds();
		if (ncbiTaxIds.size()>1) {
			LOGGER.warn("More than one taxonomy id for uniprot entry "+this.uniId);
		}
		this.taxId = ncbiTaxIds.get(0).getValue();
		for (NcbiTaxon ncbiTax:entry.getTaxonomy()) {
			taxons.add(ncbiTax.getValue());
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
			this.geneEncodingOrganelle = orglls.get(0).getType().getValue();
			if (orglls.size()>1) {
				for (Organelle orgll:orglls){ 
					if (!orgll.getType().equals(this.geneEncodingOrganelle)) {
						LOGGER.warn("Different gene encoding organelles for Uniprot "+this.uniId);
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
	 * best one is returned (as long as it is above the {@value #MIN_TOLERATED_ID}) value)
	 * @return the ProteinsToCDSMatch object corresponding to the chosen translated CDS 
	 * or null if no correct match to the protein sequence can be found
	 */
	public ProteinToCDSMatch getRepresentativeCDS() {
		if (repCDScached) {
			return representativeCDS;
		}

		if (this.geneEncodingOrganelle!=null) {
			LOGGER.info("Non-empty uniprot gene encoding organelle field: entry "+this.getUniId()+" encoded in: "+this.geneEncodingOrganelle);
			if (!this.geneEncodingOrganelle.equals(GENE_ENC_ORG_PLASMID)) { 
				LOGGER.warn("The entry "+this.getUniId()+" is not encoded in nucleus or plasmid! Not using its CDS.");
				//TODO this is a temporary solution, eventually we should let it in and use the appropriate genetic code if needed
				representativeCDS = null;
				repCDScached = true;
				return representativeCDS;
			}
		}
		
		List<ProteinToCDSMatch> allMatchings = new ArrayList<ProteinToCDSMatch>();
		for (Sequence cds:emblCdsSeqs) {
			
			ProteinToCDSMatch matching = null;
			try {
				matching = new ProteinToCDSMatch(this.getUniprotSeq(), cds, GeneticCodeType.STANDARD);
				allMatchings.add(matching);
			} catch (TranslationException e) {
				LOGGER.warn("Couldn't translate embl CDS "+cds.getName());
				LOGGER.warn(e.getMessage());
				continue; // try the next one
			}
		}		

		if (allMatchings.size()>0) {
			List<TranslatedFrame> nonFullMatchingCDSs = new ArrayList<TranslatedFrame>(); 
			for (ProteinToCDSMatch matching:allMatchings) {
				if (matching.hasFullMatch()) {
					representativeCDS = matching;
					repCDScached = true;
					return representativeCDS;
				}
				nonFullMatchingCDSs.add(matching.getBestTranslation());	
			}
			TranslatedFrame bestTranslation = Collections.max(nonFullMatchingCDSs);
			for (ProteinToCDSMatch matching:allMatchings){
				if (matching.getBestTranslation()==bestTranslation) {
					int mismatches = bestTranslation.getNumMismatches();
					int gaps = bestTranslation.getNumGaps();
					if (bestTranslation.getPercentIdentity()/100.0f<MIN_TOLERATED_ID) {
						LOGGER.warn("No fully matching CDSs for uniprot entry "+this.getUniId()+
								". Best match "+matching.getCDSName()+" (reading frame "+bestTranslation.getReadingFrame().getNumber()+
								") not good enough ("+String.format("%5.1f", bestTranslation.getPercentIdentity())+" identity). ");					
						LOGGER.info("Alignment of best translation:\n"+bestTranslation.getAln().getFormattedAlignmentString());
						representativeCDS = null;
					} else {
						representativeCDS = matching;
						LOGGER.warn("No fully matching CDSs for uniprot entry "+this.getUniId()+". Using the best match '"+matching.getCDSName()+
								"' (reading frame "+bestTranslation.getReadingFrame().getNumber()+") with "+mismatches+" mismatches. Identity: "+
								String.format("%4.1f",bestTranslation.getPercentIdentity()));
						if (gaps>0) {
							LOGGER.warn(gaps + " gaps in the CDS to uniprot alignment. Will discard this CDS.");						
							representativeCDS = null;
						}
						if (matching.hasStopCodonsInBestTranslation()) {
							LOGGER.warn("Translation contains STOP codons. Will discard this CDS");
							representativeCDS = null;
						}
						LOGGER.info("Alignment of best translation:\n"+bestTranslation.getAln().getFormattedAlignmentString());
					}
					repCDScached = true;
					return representativeCDS;
				}
			}
		}
		repCDScached = true;
		representativeCDS = null;
		return representativeCDS;
	}
	
	/**
	 * Given a protein sequence position (starting at 0), tells whether the corresponding CDS translation 
	 * properly matches (identical aminoacids) at that position
	 * @param i protein sequence index starting at 0
	 * @return
	 */
	public boolean isReliablePosition(int i) {
		return this.getRepresentativeCDS().getBestTranslation().isMatch(i);
	}
	
	public GeneticCodeType getGeneticCodeType() {
		return GeneticCodeType.getByOrganelleAndOrganism(this.geneEncodingOrganelle,this.getLastTaxon());
	}

	/*------------------------ HasFeature interface implementation -----------------------*/
	
	public boolean addFeature(Feature feature) throws InvalidFeatureCoordinatesException,
														OverlappingFeatureException {

		boolean result = false;
		FeatureType ft = feature.getType();	// will throw a null pointer exception if feature == null
		if(this.features == null) {
			this.features = new HashMap<FeatureType, Collection<Feature>>();
		}
		Collection<Feature> fc = this.features.get(ft);
		if(fc==null) {
			fc = new LinkedList<Feature>();
			features.put(ft, fc);
		}
		IntervalSet intervSet = feature.getIntervalSet();
		for (Interval interv:intervSet){
			if (interv.beg<1 || interv.end>this.getLength())
				throw new InvalidFeatureCoordinatesException("Feature being added "+feature.getDescription()+" of type "+feature.getType()+" contains invalid coordinates for this UniprotEntry.\n"
						+"Interval: "+interv+". Length of this UniprotEntry: "+this.getLength());
		}
		for (Feature f:fc) {
			if (intervSet.overlaps(f.getIntervalSet())) {
				throw new OverlappingFeatureException("Feature being added "+feature.getDescription()+" of type "+feature.getType()+" overlaps existing feature "+f.getDescription()+" of type "+f.getType()+"\n" +
						"New interval set: "+intervSet+". Existing interval set: "+f.getIntervalSet());
			}
		}
		result = fc.add(feature);
		return result;	
	}

	public Collection<FeatureType> getFeatureTypes() {
		return features.keySet();
	}

	public Collection<Feature> getFeatures() {
		Collection<Feature> allfeatures = new LinkedList<Feature>();
		for (Collection<Feature> coll:features.values()) {
			allfeatures.addAll(coll);
		}
		return allfeatures;
	}

	public Collection<Feature> getFeaturesForPositon(int position) {
		Collection<Feature> result = new LinkedList<Feature>(); 
		for(Feature f:this.getFeatures()) {
			if(f.getIntervalSet().getIntegerSet().contains(position)) result.add(f);
		}
		return result;		
	}

	public Collection<Feature> getFeaturesOfType(FeatureType featureType) {
		return features.get(featureType);
	}

	public Collection<Feature> getFeaturesOfTypeForPosition(FeatureType featureType, int position) {
		Collection<Feature> result = new LinkedList<Feature>(); 
		if(this.features.get(featureType) != null) {
			for(Feature f:this.features.get(featureType)) {
				if(f.getIntervalSet().getIntegerSet().contains(position)) result.add(f);
			}
		}
		return result;
	}


}
