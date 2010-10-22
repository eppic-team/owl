package owl.tests.core.structure;


import java.util.Collection;

import javax.vecmath.Matrix3d;
import javax.vecmath.Matrix4d;

import junit.framework.Assert;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import Jama.Matrix;

import owl.core.structure.CrystalCell;
import owl.core.structure.SpaceGroup;
import owl.core.structure.SymoplibParser;

public class SpaceGroupTest {

	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}
	
	@Test
	public void testTransfConversion() {
		Collection<SpaceGroup> allSGs = SymoplibParser.getAllSpaceGroups().values();
		
		int countEn = 0;
		int countNonEn = 0;
		int countSpecial = 0;
		for (SpaceGroup spaceGroup:allSGs) {

			if (spaceGroup.isEnantiomorphic() && spaceGroup.getId()<1000) {
				countEn++;
			}
			if (!spaceGroup.isEnantiomorphic() && spaceGroup.getId()<1000) {
				countNonEn++;
			}
			if (spaceGroup.getId()>1000) {
				countSpecial++;
			}
			for (int i=0;i<spaceGroup.getNumOperators();i++){
				CrystalCell unitCell = spaceGroup.getBravLattice().getExampleUnitCell();
				Matrix4d m = spaceGroup.getTransformation(i);
				double traceBefore = m.m00+m.m11+m.m22+m.m33;

				unitCell.transfToOrthonormal(m);
				double traceAfter = m.m00+m.m11+m.m22+m.m33;
				Assert.assertEquals(traceBefore, traceAfter);
				
				if (spaceGroup.isEnantiomorphic()) {
					Matrix3d rot = new Matrix3d(m.m00,m.m01,m.m02,m.m10,m.m11,m.m12,m.m20,m.m21,m.m22);
					
					// determinant must be 1
					Assert.assertEquals(1.0,rot.determinant());
					
					// at least 1 eigenvalue must be 1
					double[][] ar = {{m.m00,m.m01,m.m02},{m.m10,m.m11,m.m12},{m.m20,m.m21,m.m22}};
					Matrix mat = new Matrix(ar);
					double[] eigenv = mat.svd().getSingularValues();
					Assert.assertTrue(eigenv[0]==1.0 || eigenv[1]==1.0 || eigenv[2]==1.0);
					
					// transpose must be equals to inverse
					Matrix3d rotTransp = new Matrix3d();
					Matrix3d rotInv = new Matrix3d();
					rotTransp.transpose(rot);
					rotInv.invert(rot);
					//if (!rotTransp.equals(rotInv)) {
					//	System.out.println(spaceGroup.getId()+" "+i);
					//}
					Assert.assertEquals(rotTransp,rotInv);
					
				} 		

			}	
		}
		
		Assert.assertEquals(266, allSGs.size()); // the total count must be 266
		Assert.assertEquals(65, countEn);        // enantiomorphic groups (protein crystallography groups)
		Assert.assertEquals(165, countNonEn);    // i.e. 266-65
		Assert.assertEquals(36, countSpecial);   // the rest of the groups present un symop.lib (sometimes used in PDB)
	}

	// to debug the testing code (run as java program so that we can use normal debugger)
	public static void main(String[] args) throws Exception {
		SpaceGroupTest test = new SpaceGroupTest();
		test.testTransfConversion();
	}

}
