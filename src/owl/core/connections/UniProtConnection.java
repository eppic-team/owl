package owl.core.connections;

import java.util.Collection;
import java.util.List;

import uk.ac.ebi.kraken.interfaces.uniprot.DatabaseType;
import uk.ac.ebi.kraken.interfaces.uniprot.UniProtEntry;
import uk.ac.ebi.kraken.interfaces.uniprot.dbx.embl.Embl;
import uk.ac.ebi.kraken.interfaces.uniprot.dbx.pdb.Pdb;
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
	 * Gets a single Uniprot entry given its Uniprot identifier. 
	 * @param uniProtId
	 * @return
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
	 * @return
	 */
	public String getVersion() {
		return UniProtJAPI.factory.getVersion();
	}
	
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
