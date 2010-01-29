package proteinstructure;

import java.io.IOException;

/**
 * An interface for multiple structural aligners.
 * See implementation PaulStructAligner
 * @author duarte
 *
 */
public interface StructAligner {
	
	/**
	 * Structurally aligns the given templates. 
	 * The templates must have pdb data loaded if the structural alignment method 
	 * is based purely on the structures. Otherwise they must have graph data loaded if 
	 * the structural alignment is based on contact maps.
	 * @param templates
	 * @return
	 * @throws StructAlignmentError if a problem occurs while running the structural alignment
	 * @throws IOException
	 */
	public Alignment alignStructures(TemplateList templates) throws StructAlignmentError, IOException;
	
	/**
	 * Structurally aligns the given pdbs 
	 * @param pdbs
	 * @param tags
	 * @return
	 * @throws StructAlignmentError if a problem occurs while running the structural alignment
	 * @throws IOException
	 */
	public Alignment alignStructures(Pdb[] pdbs, String[] tags) throws StructAlignmentError, IOException;
	
	/**
	 * Performs structural alignment of the given contact maps (graphs)
	 * @param graphs
	 * @param tags
	 * @return
	 * @throws StructAlignmentError if a problem occurs while running the structural alignment
	 * @throws IOException
	 */
	public Alignment alignStructures(RIGraph[] graphs, String[] tags) throws StructAlignmentError, IOException;

}
