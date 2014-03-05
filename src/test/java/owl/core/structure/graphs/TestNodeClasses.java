package owl.core.structure.graphs;


import static org.junit.Assert.*;

import javax.vecmath.Point3d;

import org.junit.Test;

import owl.core.structure.AaResidue;
import owl.core.structure.AminoAcid;
import owl.core.structure.Atom;


public class TestNodeClasses {

	@Test
	public void testAtom() {
		
		Atom atom1 = new Atom(1,"CA","C",new Point3d(0,0,0),new AaResidue(AminoAcid.ALA, 1, null), 1.0, 0.0);
		Atom atom2 = new Atom(1,"CA","C",new Point3d(0,0,0),new AaResidue(AminoAcid.ALA, 1, null), 1.0, 0.0);
		Atom atom3 = new Atom(1,"CB","C",new Point3d(0,0,0),new AaResidue(AminoAcid.ALA, 1, null), 1.0, 0.0);
		Atom atom4 = new Atom(2,"CA","C",new Point3d(0,0,0),new AaResidue(AminoAcid.ALA, 1, null), 1.0, 0.0);
		Atom atom5 = new Atom(1,"CA","C",new Point3d(0,0,0),new AaResidue(AminoAcid.SER, 1, null), 1.0, 0.0);
		Atom atom6 = new Atom(1,"CA","C",new Point3d(0,0,0),new AaResidue(AminoAcid.ALA, 2, null), 1.0, 0.0);
		
		assertEquals(atom1, atom2);
		assertEquals(atom1.hashCode(), atom2.hashCode());
		
		assertNotSame(atom1,atom3);
		assertNotSame(atom1.hashCode(),atom3.hashCode());
		
		assertNotSame(atom1, atom4);
		assertNotSame(atom1.hashCode(), atom4.hashCode());

		assertNotSame(atom1, atom5);
		assertNotSame(atom1.hashCode(), atom5.hashCode());
		
		assertNotSame(atom1, atom6);
		assertNotSame(atom1.hashCode(), atom6.hashCode());

	}
	
	@Test
	public void testAIGNode() {
		
		AIGNode node1 = new AIGNode(1, "CA", new RIGNode(1, "ALA", null));
		AIGNode node2 = new AIGNode(1, "CA", new RIGNode(1, "ALA", null));
		
		AIGNode node3 = new AIGNode(2, "CA", new RIGNode(1, "ALA", null));
		AIGNode node4 = new AIGNode(1, "CB", new RIGNode(1, "ALA", null));
		AIGNode node5 = new AIGNode(1, "CA", new RIGNode(2, "ALA", null));
		AIGNode node6 = new AIGNode(1, "CA", new RIGNode(1, "SER", null));
		
		assertEquals(node1, node2);
		// TODO hashCode is not implemented!
		//assertEquals(node1.hashCode(), node2.hashCode());
		
	
		assertNotSame(node1, node3);
		// TODO hashCode is not implemented!
		//assertNotSame(node1.hashCode(), node3.hashCode());
		
		assertNotSame(node1, node4);
		// TODO hashCode is not implemented!
		//assertNotSame(node1.hashCode(), node4.hashCode());

		assertNotSame(node1, node5);
		// TODO hashCode is not implemented!
		//assertNotSame(node1.hashCode(), node5.hashCode());

		assertNotSame(node1, node6);
		// TODO hashCode is not implemented!
		//assertNotSame(node1.hashCode(), node6.hashCode());
		
	}
	
	@Test
	public void testRIGNode() {
		
		RIGNode node1 = new RIGNode(1, "ALA", null);
		RIGNode node2 = new RIGNode(1, "ALA", null);
		RIGNode node3 = new RIGNode(2, "ALA", null);
		RIGNode node4 = new RIGNode(1, "SER", null);
		
		assertEquals(node1, node2);
		// TODO hashCode is not implemented!
		//assertEquals(node1.hashCode(), node2.hashCode());	

		assertNotSame(node1, node3);
		// TODO hashCode is not implemented!
		//assertNotSame(node1.hashCode(), node3.hashCode());

		assertNotSame(node1, node4);
		// TODO hashCode is not implemented!
		//assertNotSame(node1.hashCode(), node4.hashCode());
		
		
	}
	
	@Test
	public void testRICGNode() {
		
		RICGNode node1 = new RICGNode(1, "ALA", "A");
		RICGNode node2 = new RICGNode(1, "ALA", "A");
		
		RICGNode node3 = new RICGNode(2, "ALA", "A");
		RICGNode node4 = new RICGNode(1, "SER", "A");
		RICGNode node5 = new RICGNode(1, "ALA", "B");
		
		assertEquals(node1, node2);
		assertEquals(node1.hashCode(), node2.hashCode());
		
	
		assertNotSame(node1, node3);
		assertNotSame(node1.hashCode(), node3.hashCode());
		
		assertNotSame(node1, node4);
		assertNotSame(node1.hashCode(), node4.hashCode());

		assertNotSame(node1, node5);
		assertNotSame(node1.hashCode(), node5.hashCode());

		
	}
	
}
