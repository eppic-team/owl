package features;

/**
 * Enumeration of features implementing the Feature interface.
 * Each subclass implementing the interface should have a corresponding FeatureType.
 * 
 * See also: {@link Feature}
 * 
 * @author stehr
 */
public enum FeatureType {
	GENERAL,		// A general, not further specified feature
	CSA,			// A catalytic site from Catalytic Site Atlas
	UNIPROT,		// A feature based on a Uniprot annotation
	PROSITE,		// A feature based on a matching Prosite motif
	PHOSPHOSITE;	// A feature based on a modification site from PhoshoSitePlus  
	//SCOP,			// A Scop domain annotation
	//EC;			// A catalytic domain as defined by EC
}
