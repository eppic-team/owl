package owl.core.connections;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import owl.core.connections.NoMatchFoundException;

import uk.ac.ebi.kraken.interfaces.uniprot.DatabaseCrossReference;
import uk.ac.ebi.kraken.interfaces.uniprot.DatabaseType;
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
 * TODO: Merge with owl.core.connections.UniprotConnection
 * @author stehr
 */
public class UniProtConnection {
	
	/*--------------------------- member variables --------------------------*/
	EntryRetrievalService entryRetrievalService;
	UniProtQueryService uniProtQueryService;
	
	/*----------------------------- constructors ----------------------------*/
	
	public UniProtConnection() {
		// TODO: Create logger
		//Create entry retrieval service
		entryRetrievalService = UniProtJAPI.factory.getEntryRetrievalService();
	    // Create UniProt query service
	    uniProtQueryService = UniProtJAPI.factory.getUniProtQueryService();
	}
	
	/*---------------------------- public methods ---------------------------*/
	
	/**
	 * Primary method to retrieve a Uniprot entry by its ID. The entry object
	 * can be used for subsequent method calls or for functions provided by the
	 * Uniprot API for which we do not have our own implementation.
	 */
	public UniProtEntry getEntry(String uniProtId) throws NoMatchFoundException {
		UniProtEntry entry = (UniProtEntry) entryRetrievalService.getUniProtEntry(uniProtId);
		if (entry==null) throw new NoMatchFoundException("No Uniprot entry found for uniprot id "+uniProtId);
		return entry;
	}

	/**
	 * Gets a list of Uniprot entries as an Iterator given a list of Uniprot identifiers.
	 * @param idsList a list of uniprot ids
	 * @return
	 */
	public EntryIterator<UniProtEntry> getMultipleEntries(List<String> idsList) {

	    Query query = UniProtQueryBuilder.buildIDListQuery(idsList);

	    return uniProtQueryService.getEntryIterator(query);
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
	    while (!(uniProtQueryService.checkStatus(jobid) == JobStatus.DONE)) {
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
	 * @return a (possible empty) collection of pdb xref objects
	 */
	public static Collection<UniprotPdbRef> getPdbRefs(UniProtEntry entry) {
		ArrayList<UniprotPdbRef> ret = new ArrayList<UniprotPdbRef>();		
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
	    	UniProtEntry entry = uniprot.getEntry("P12830"); // CDH1
	    	System.out.println("entry = " + entry.getUniProtId().getValue());
		    // print PDB cross references
	    	if(entry != null) {
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
