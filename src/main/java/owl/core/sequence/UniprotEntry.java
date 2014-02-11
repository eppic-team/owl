package owl.core.sequence;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import owl.core.connections.NoMatchFoundException;
import owl.core.connections.UniProtConnection;
import owl.core.features.Feature;
import owl.core.features.FeatureType;
import owl.core.features.HasFeatures;
import owl.core.features.InvalidFeatureCoordinatesException;
import owl.core.features.OverlappingFeatureException;
import owl.core.util.Interval;
import owl.core.util.IntervalSet;

import uk.ac.ebi.kraken.interfaces.uniprot.NcbiTaxon;
import uk.ac.ebi.kraken.interfaces.uniprot.NcbiTaxonomyId;
import uk.ac.ebi.kraken.interfaces.uniprot.Organelle;
import uk.ac.ebi.kraken.interfaces.uniprot.UniProtEntry;

/**
 * A uniprot entry.
 * The class stores the sequence and taxonomy data 
 *  
 * @author duarte_j
 *
 */
public class UniprotEntry implements HasFeatures, Serializable {
	
	private static final long serialVersionUID = 1L;

	private static final Log LOGGER = LogFactory.getLog(UniprotEntry.class);
		
	private String uniId;
	private Sequence uniprotSeq;
	private int taxId;
	private List<String> taxons; // the taxonomy as returned by UniprotEntry.getTaxonomy(). First value is the most general (kingdom), last is the most specific (species).
	private String geneEncodingOrganelle; // the organelle where this gene is encoded (important for genetic code): if gene encoded in nucleus this is null
	

	private transient Map<FeatureType, Collection<Feature>> features;
	
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
	public int getTaxId() {
		return taxId;
	}

	/**
	 * @param taxId the taxId to set
	 */
	public void setTaxId(int taxId) {
		this.taxId = taxId;
	}
	
	/**
	 * The taxonomy as returned by UniprotEntry.getTaxonomy(). 
	 * First value is the most general (domain of life), last is the most specific (species).
	 * @return
	 */
	public List<String> getTaxons() {
		return taxons;
	}
	
	/**
	 * Returns the domain of life (what used to be kingdom) for this entry 
	 * @return
	 */
	public String getFirstTaxon() {
		return this.taxons.get(0);
	}
	
	/**
	 * Returns the most specific taxonomy annotation (species) for this entry
	 * @return
	 */
	public String getLastTaxon() {
		return this.taxons.get(this.taxons.size()-1);
	}
	
	public void setTaxons(List<String> taxons) {
		this.taxons = taxons;
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
	 * Retrieves from UniprotKB the sequence, taxonomy and gene encoding organelle 
	 * by using the remote Uniprot API
	 */
	public void retrieveUniprotKBData() throws NoMatchFoundException {
		this.taxons = new ArrayList<String>();
		
		UniProtConnection uniprotConn = new UniProtConnection();
		UniProtEntry entry = uniprotConn.getEntry(uniId);
		
		this.setUniprotSeq(new Sequence(this.getUniId(),entry.getSequence().getValue()));
		
		List<NcbiTaxonomyId> ncbiTaxIds = entry.getNcbiTaxonomyIds();
		if (ncbiTaxIds.size()>1) {
			LOGGER.warn("More than one taxonomy id for uniprot entry "+this.uniId);
		}
		this.taxId = Integer.parseInt(ncbiTaxIds.get(0).getValue());
		for (NcbiTaxon ncbiTax:entry.getTaxonomy()) {
			taxons.add(ncbiTax.getValue());
		}

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
	
	public GeneticCodeType getGeneticCodeType() {
		return GeneticCodeType.getByOrganelleAndOrganism(this.geneEncodingOrganelle,this.getLastTaxon());
	}
	
	/**
	 * Returns true if this UniprotEntry belongs to same domain of life (Bacteria, Archaea, Eukaryota) 
	 * as the one given
	 * @param other
	 * @return
	 */
	public boolean isInSameDomainOfLife(UniprotEntry other) {
		return this.getFirstTaxon().equals(other.getFirstTaxon());
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
