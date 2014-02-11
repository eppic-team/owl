package owl.core.connections;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import owl.core.connections.NoMatchFoundException;
import owl.core.sequence.UnirefEntry;

import uk.ac.ebi.kraken.interfaces.uniprot.DatabaseCrossReference;
import uk.ac.ebi.kraken.interfaces.uniprot.DatabaseType;
import uk.ac.ebi.kraken.interfaces.uniprot.NcbiTaxon;
import uk.ac.ebi.kraken.interfaces.uniprot.NcbiTaxonomyId;
import uk.ac.ebi.kraken.interfaces.uniprot.UniProtEntry;
import uk.ac.ebi.kraken.interfaces.uniprot.dbx.embl.Embl;
import uk.ac.ebi.kraken.interfaces.uniprot.dbx.pdb.Pdb;
import uk.ac.ebi.kraken.interfaces.uniprot.dbx.phosphosite.PhosphoSite;
import uk.ac.ebi.kraken.interfaces.uniprot.features.BindingFeature;
import uk.ac.ebi.kraken.interfaces.uniprot.features.Feature;
import uk.ac.ebi.kraken.interfaces.uniprot.features.FeatureType;
import uk.ac.ebi.kraken.model.blast.JobStatus;
import uk.ac.ebi.kraken.model.blast.parameters.DatabaseOptions;
import uk.ac.ebi.kraken.uuw.services.remoting.EntryIterator;
import uk.ac.ebi.kraken.uuw.services.remoting.EntryRetrievalService;
import uk.ac.ebi.kraken.uuw.services.remoting.Query;
import uk.ac.ebi.kraken.uuw.services.remoting.UniProtJAPI;
import uk.ac.ebi.kraken.uuw.services.remoting.UniProtQueryBuilder;
import uk.ac.ebi.kraken.uuw.services.remoting.UniProtQueryService;
import uk.ac.ebi.kraken.uuw.services.remoting.blast.BlastData;
import uk.ac.ebi.kraken.uuw.services.remoting.blast.BlastInput;

/**
 * Our interface to the Uniprot Java API.
 * 
 * @author stehr
 */
public class UniProtConnection {
	
	private static final Log LOGGER = LogFactory.getLog(UniProtConnection.class);
	
	/*--------------------------- member variables --------------------------*/
	private EntryRetrievalService entryRetrievalService;
	private UniProtQueryService uniProtQueryService;
	
	private HashSet<String> nonReturnedIdsLastMultipleRequest;
	
	/*----------------------------- constructors ----------------------------*/
	
	public UniProtConnection() {
		// Create entry retrieval service
		entryRetrievalService = UniProtJAPI.factory.getEntryRetrievalService();
	    // Create UniProt query service
	    uniProtQueryService = UniProtJAPI.factory.getUniProtQueryService();
	}
	
	/*---------------------------- public methods ---------------------------*/
	
	/**
	 * Primary method to retrieve a Uniprot entry by its ID. The entry object
	 * can be used for subsequent method calls or for functions provided by the
	 * Uniprot API for which we do not have our own implementation.
	 * @param uniProtId
	 * @return 
	 * @throws NoMatchFoundException if no match returned by UniProt JAPI
	 */
	public UniProtEntry getEntry(String uniProtId) throws NoMatchFoundException {
		UniProtEntry entry = (UniProtEntry) entryRetrievalService.getUniProtEntry(uniProtId);
		if (entry==null) throw new NoMatchFoundException("No Uniprot entry found for uniprot id "+uniProtId);
		return entry;
	}

	/**
	 * Convenience method to get a owl.core.sequence.UnirefEntry object from UniProt JAPI 
	 * given a uniprot id.
	 * The UnirefEntry object returned contains the uniprot id, tax id, taxons and sequence.
	 * The information returned with this object can be obtained by calling {@link #getEntry(String)} and
	 * then extracting the different parts of the info from the returned JAPI's UniProtEntry
	 * @param uniProtId
	 * @return
	 * @throws NoMatchFoundException if no match returned by UniProt JAPI
	 */
	public UnirefEntry getUnirefEntry(String uniProtId) throws NoMatchFoundException {
		List<String> taxons = new ArrayList<String>();
		
		UniProtEntry entry = getEntry(uniProtId);
		String sequence = entry.getSequence().getValue();
		
		List<NcbiTaxonomyId> ncbiTaxIds = entry.getNcbiTaxonomyIds();
		if (ncbiTaxIds.size()>1) {
			LOGGER.warn("More than one taxonomy id for uniprot entry "+uniProtId);
		}
		int ncbiTaxId = Integer.parseInt(ncbiTaxIds.get(0).getValue());
		for (NcbiTaxon ncbiTax:entry.getTaxonomy()) {
			taxons.add(ncbiTax.getValue());
		}
		UnirefEntry uniref = new UnirefEntry();
		uniref.setUniprotId(uniProtId);
		uniref.setNcbiTaxId(ncbiTaxId);
		uniref.setTaxons(taxons);
		uniref.setSequence(sequence);
		return uniref;
	}
	
	/**
	 * Gets a list of Uniprot entries as an Iterator given a list of Uniprot identifiers.
	 * If any of the input entries can not be retrieved through JAPI then they will be
	 * missing in the returned iterator. The user must check for those.  
	 * @param idsList a list of uniprot ids
	 * @return
	 */
	public EntryIterator<UniProtEntry> getMultipleEntries(List<String> idsList) {
	    Query query = UniProtQueryBuilder.buildIDListQuery(idsList);
	    return uniProtQueryService.getEntryIterator(query);
	}
	
	/**
	 * Convenience method to get a List of owl.core.sequence.UnirefEntry given a List of uniprot ids.
	 * Analogous to {@link #getUnirefEntry(String)} but for multiple entries.
	 * If the JAPI does not return all requested ids a warning is logged and the list of non-returned 
	 * ids can be retrieved through {@link #getNonReturnedIdsLastMultipleRequest()}
	 * @param uniprotIds
	 * @return
	 * @throws IOException
	 */
	public List<UnirefEntry> getMultipleUnirefEntries(List<String> uniprotIds) throws IOException {
		
		List<UnirefEntry> unirefEntries = new ArrayList<UnirefEntry>();
		
		EntryIterator<UniProtEntry> entries = getMultipleEntries(uniprotIds);

		for (UniProtEntry entry:entries) {
			String uniId = entry.getPrimaryUniProtAccession().getValue();
			if (!uniprotIds.contains(uniId)) { // TODO this could be more efficient by using a Map, is it necessary?
				// this happens if the JAPI/server are really broken and return records that we didn't ask for (actually happened on the 09.02.2011!!!)
				throw new IOException("Uniprot JAPI server returned an unexpected record: "+uniId);
			}
			String sequence = entry.getSequence().getValue();

			List<NcbiTaxonomyId> ncbiTaxIds = entry.getNcbiTaxonomyIds();
			if (ncbiTaxIds.size()>1) {
				LOGGER.warn("More than one taxonomy id for uniprot entry "+uniId);
			}
			int ncbiTaxId = Integer.parseInt(ncbiTaxIds.get(0).getValue());
			
			List<String> taxons = new ArrayList<String>();
			for(NcbiTaxon ncbiTaxon:entry.getTaxonomy()) {
				taxons.add(ncbiTaxon.getValue());
			}
			UnirefEntry uniref = new UnirefEntry();
			uniref.setUniprotId(uniId);
			uniref.setNcbiTaxId(ncbiTaxId);
			uniref.setTaxons(taxons);
			uniref.setSequence(sequence);
			unirefEntries.add(uniref);
		}
		
		// now we check if the query to uniprot JAPI did really return all requested uniprot ids
	    HashSet<String> returnedUniIds = new HashSet<String>();
	    for (UnirefEntry uniref:unirefEntries) {
			returnedUniIds.add(uniref.getUniprotId());
	    }
	    nonReturnedIdsLastMultipleRequest = new HashSet<String>();
	    for (String uniprotId:uniprotIds){
	    	if (!returnedUniIds.contains(uniprotId)) {
	    		nonReturnedIdsLastMultipleRequest.add(uniprotId);
	    		LOGGER.warn("Information for uniprot ID "+uniprotId+" could not be retrieved with the Uniprot JAPI.");
	    	}
	    }
		
		return unirefEntries;
	}
	
	public HashSet<String> getNonReturnedIdsLastMultipleRequest() {
		return nonReturnedIdsLastMultipleRequest;
	}
	
	/**
	 * Return the UniProt version this connection is connected to
	 * @return a Uniprot version string
	 */
	public String getVersion() {
		return UniProtJAPI.factory.getVersion();
	}
	
	/**
	 * Use the Uniprot blast interface to find an entry by sequence.
	 * @param sequence the query sequence
	 * @return the uniprot entry with the best match to the query sequence
	 */
	public UniProtEntry getHumanEntryBySequence(String sequence) {
	    //Create a blast input with a Database and sequence
	    BlastInput input = new BlastInput(DatabaseOptions.UNIPROT_HUMAN, sequence);
	    //Submitting the input to the service will return a job id
	    String jobid = uniProtQueryService.submitBlast(input);
	    
	    System.out.print("Waiting for Blast result...");
	    
	    //Use this jobid to check the service to see if the job is complete
	    while (!(uniProtQueryService.checkStatus(jobid) == JobStatus.FINISHED)) {
		    try {
		      //Sleep a bit before the next request
		          System.out.print(".");
		          Thread.sleep(5000);
		    } catch (InterruptedException e) {
		          e.printStackTrace();
		    }
	    }
	    //The blast data contains the job information and the hits with entries
	    BlastData<UniProtEntry> blastResult = uniProtQueryService.getResults(jobid);
	    UniProtEntry bestHit = blastResult.getBlastHits().get(0).getEntry();
		String description = blastResult.getBlastHits().get(0).getHit().getDescription();
		long length = blastResult.getBlastHits().get(0).getHit().getLength();
		float seqId = blastResult.getBlastHits().get(0).getHit().getAlignments().get(0).getIdentity();
	    if(length != sequence.length()) {
	    	System.err.println("Warning: Blast hit is not full length");
	    }
	    if(seqId < 100) {
	    	System.err.println("Warning: Blast hit has sequence identity < 100%");
	    }
	    if(!sequence.equals(bestHit.getSequence().getValue())) {
	    	System.err.println("Warning: Retrieved sequence is not identical to query sequence");
	    }
	    System.out.println();
	    System.out.print("Found: ");
	    System.out.println(description);
	    //System.out.println("length = " + length);
	    //System.out.println("seqid = " + seqId);	    
		return bestHit;
	}
	
	/*---------------------------- static methods ---------------------------*/
	
	/**
	 * Returns a collection of pdb cross references for the given Uniprot entry.
	 * Notes: 
	 * - If a xref contains multiple chains, only the first in returned.
	 * @return a (possible empty) collection of pdb xref objects
	 */
	public static Collection<UniProtPdbRef> getPdbRefs(UniProtEntry entry) {
		ArrayList<UniProtPdbRef> ret = new ArrayList<UniProtPdbRef>();
		String uniprotId = entry.getPrimaryUniProtAccession().getValue();
		// get PDB cross references
		Collection<Pdb> refs = entry.getDatabaseCrossReferences(DatabaseType.PDB);
		for(Pdb ref:refs) {
			String pdbCode = ref.getPdbAccessionNumber().getValue();
			String chainRecordsStr = ref.getPdbChains().getValue();	
			String[] chainRecords = chainRecordsStr.split(",");	// in rare cases, there can be multiple of these
			String chainRecord = chainRecords[0];				// take first, ignore others, TODO: log WARNING
			String[] fields = chainRecord.trim().split("=");	// e.g. A/B/C/D=319-360
			if(fields.length < 2) {
				System.err.println("Error parsing xref from Uniprot for " + uniprotId + ". Sequence range not found.");
				continue;
			}
			String[] chains = fields[0].split("/");
			String[] begEnd = fields[1].split("-");
			int beg = Integer.parseInt(begEnd[0]);
			int end = Integer.parseInt(begEnd[1]);	    		
			String method = ref.getPdbMethod().getValue();
			String resStr = ref.getPdbResolution().getValue();
			double res = resStr.equals("-")?Double.NaN:Double.parseDouble(resStr.split(" ")[0]);
			//DEBUG: System.out.printf("%s:%s:%.2f:%s:%d-%d\n", pdbCode, method, res, chains[0], beg, end);			
			String geneName = entry.getGenes().get(0).getGeneName().getValue();
			UniProtPdbRef xref = new UniProtPdbRef(geneName, uniprotId, pdbCode, chains, method, res, beg, end);
			ret.add(xref);
		}
		return ret;
	}
	
	/**
	 * Returns the start location of an ATP binding site, or 0 if no such site is annotated in the entry.
	 */
	public static int getAtpBindingSite(UniProtEntry entry) {
		for(Feature f:entry.getFeatures(FeatureType.BINDING)) {
			BindingFeature f2 = (BindingFeature) f;
			if(f2.getFeatureDescription().getValue().startsWith("ATP")) {
				System.out.println("Found binding feature: " + f2.getFeatureDescription().getValue());
				System.out.println("Location: " + f2.getFeatureLocation().getStart() + "-" + f2.getFeatureLocation().getEnd());
				return f2.getFeatureLocation().getStart();
			}
		}
		return 0;
	}
	
	/**
	 * Returns the URL for accessing the PhosphoSitePlus record for this entry.
	 * @param entry the Uniprot entry
	 * @return the PhosphoSite URL or null if no entry was found
	 */
	public static String getPhosphoSiteUrl(UniProtEntry entry) {
		
		String baseUrl="http://www.phosphosite.org/uniprotAccAction.do?id=";
		
		// obtain id (most likely just the sp_id but just to make sure...)
		String ac = null;
		List<DatabaseCrossReference> list = entry.getDatabaseCrossReferences(DatabaseType.PHOSPHOSITE);
		if(list != null && list.size() > 0) {
			PhosphoSite phosphoSiteXRef = ((PhosphoSite)list.get(0));
			ac = phosphoSiteXRef.toString();
		} else {
			ac = entry.getUniProtId().getValue();
		}
		return baseUrl + ac;
	}
	
	/*--------------------------------- main --------------------------------*/
	
	public static void main(String[] args) {

		UniProtConnection uniprot = new UniProtConnection();
		
	    //Use the factory to print out the version.
	    System.out.println("UniProt Version = " + uniprot.getVersion());

	    //Retrieve UniProt entry by its accession number
	    try {
	    	String queryId = "P12830";
	    	System.out.println("query = " + queryId);
	    	UniProtEntry entry = uniprot.getEntry(queryId); // CDH1
		    // print PDB cross references
	    	if(entry != null) {
		    	System.out.println("entry = " + entry.getUniProtId().getValue());
		    	System.out.println("sp_id = " + entry.getPrimaryUniProtAccession().getValue());
		    	
	    		Collection<Pdb> refs = entry.getDatabaseCrossReferences(DatabaseType.PDB);
	    		for(Pdb ref:refs) {
	    			System.out.println(ref);
	    		}
	    		Collection<Embl> emblrefs = entry.getDatabaseCrossReferences(DatabaseType.EMBL);
	    		for(Embl ref:emblrefs) {
	    			System.out.println(ref.getDatabase()+" "+ref.getEmblAccessionNumber()+" "+ref.getEmblMoleculeType()+" "+ref.getEmblProteinId());
	    			
	    		}
	    	}
	    } catch (NoMatchFoundException e) {
	    	System.err.println(e.getMessage());
	    }
	    System.exit(0);
	    
	    // get hit by blast
	    String seq = "MELAALCRWGLLLALLPPGAASTQVCTGTDMKLRLPASPETHLDMLRHLYQGCQVVQGNLELTYLPTNASLSFLQDIQEVQGYVLIAHNQVRQVPLQRLRIVRGTQLFEDNYALAVLDNGDPLNNTTPVTGASPGGLRELQLRSLTEILKGGVLIQRNPQLCYQDTILWKDIFHKNNQLALTLIDTNRSRACHPCSPMCKGSRCWGESSEDCQSLTRTVCAGGCARCKGPLPTDCCHEQCAAGCTGPKHSDCLACLHFNHSGICELHCPALVTYNTDTFESMPNPEGRYTFGASCVTACPYNYLSTDVGSCTLVCPLHNQEVTAEDGTQRCEKCSKPCARVCYGLGMEHLREVRAVTSANIQEFAGCKKIFGSLAFLPESFDGDPASNTAPLQPEQLQVFETLEEITGYLYISAWPDSLPDLSVFQNLQVIRGRILHNGAYSLTLQGLGISWLGLRSLRELGSGLALIHHNTHLCFVHTVPWDQLFRNPHQALLHTANRPEDECVGEGLACHQLCARGHCWGPGPTQCVNCSQFLRGQECVEECRVLQGLPREYVNARHCLPCHPECQPQNGSVTCFGPEADQCVACAHYKDPPFCVARCPSGVKPDLSYMPIWKFPDEEGACQPCPINCTHSCVDLDDKGCPAEQRASPLTSIISAVVGILLVVVLGVVFGILIKRRQQKIRKYTMRRLLQETELVEPLTPSGAMPNQAQMRILKETELRKVKVLGSGAFGTVYKGIWIPDGENVKIPVAIKVLRENTSPKANKEILDEAYVMAGVGSPYVSRLLGICLTSTVQLVTQLMPYGCLLDHVRENRGRLGSQDLLNWCMQIAKGMSYLEDVRLVHRDLAARNVLVKSPNHVKITDFGLARLLDIDETEYHADGGKVPIKWMALESILRRRFTHQSDVWSYGVTVWELMTFGAKPYDGIPAREIPDLLEKGERLPQPPICTIDVYMIMVKCWMIDSECRPRFRELVSEFSRMARDPQRFVVIQNEDLGPASPLDSTFYRSLLEDDDMGDLVDAEEYLVPQQGFFCPDPAPGAGGMVHHRHRSSSTRSGGGDLTLGLEPSEEEAPRSPLAPSEGAGSDVFDGDLGMGAAKGLQSLPTHDPSPLQRYSEDPTVPLPSETDGYVAPLTCSPQPEYVNQPDVRPQPPSPREGPLPAARPAGATLERPKTLSPGKNGVVKDVFAFGGAVENPEYLTPQGGAAPQPHPPPAFSPAFDNLYYWDQDPPERGAPPSTFKGTPTAENPEYLGLDVPV";
	    UniProtEntry entry2 = uniprot.getHumanEntryBySequence(seq);
	    
	    if (entry2 != null) {
		      System.out.println("entry2 = " + entry2.getUniProtId().getValue());
		      if(getAtpBindingSite(entry2) > 0) {
		    	  System.out.println("ATP binding site found.");
		      }
		}
	    
	    System.out.println("done.");
	}
}
