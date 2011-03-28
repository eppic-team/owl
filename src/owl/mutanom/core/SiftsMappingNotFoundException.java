package owl.mutanom.core;

/**
 * Exception being thrown if no SIFTS mapping for a given Uniprot ID has been found.
 * TODO: Merge with owl.core.connection.MappingNotFoundException
 * @author stehr
 */
public class SiftsMappingNotFoundException extends Exception {
	private static final long serialVersionUID = 1L;
	public SiftsMappingNotFoundException(String message) {
		super(message);
	}
}
