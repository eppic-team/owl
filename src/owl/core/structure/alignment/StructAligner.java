package owl.core.structure.alignment;

import java.io.IOException;

import owl.core.sequence.alignment.MultipleSequenceAlignment;
import owl.core.structure.Pdb;
import owl.core.structure.TemplateList;
import owl.core.structure.graphs.RIGraph;

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
	 * @throws StructAlignmentException if a problem occurs while running the structural alignment
	 * @throws IOException
	 */
	public MultipleSequenceAlignment alignStructures(TemplateList templates) throws StructAlignmentException, IOException;
	
	/**
	 * Structurally aligns the given pdbs 
	 * @param pdbs
	 * @param tags
	 * @return
	 * @throws StructAlignmentException if a problem occurs while running the structural alignment
	 * @throws IOException
	 */
	public MultipleSequenceAlignment alignStructures(Pdb[] pdbs, String[] tags) throws StructAlignmentException, IOException;
	
	/**
	 * Performs structural alignment of the given contact maps (graphs)
	 * @param graphs
	 * @param tags
	 * @return
	 * @throws StructAlignmentException if a problem occurs while running the structural alignment
	 * @throws IOException
	 */
	public MultipleSequenceAlignment alignStructures(RIGraph[] graphs, String[] tags) throws StructAlignmentException, IOException;

}
