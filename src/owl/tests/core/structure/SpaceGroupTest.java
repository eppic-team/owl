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

	private static final double DELTA = 0.000001;
	
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
				Matrix4d mT = unitCell.transfToOrthonormal(m);
				
				// traces before and after transformation must coincide
				Assert.assertEquals(m.m00+m.m11+m.m22+m.m33, mT.m00+mT.m11+mT.m22+mT.m33, DELTA);

				if (spaceGroup.isEnantiomorphic() && spaceGroup.getId()<1000) {
					Matrix3d rot = new Matrix3d(mT.m00,mT.m01,mT.m02,mT.m10,mT.m11,mT.m12,mT.m20,mT.m21,mT.m22);

					// determinant must be 1
					Assert.assertEquals(1.0,rot.determinant(), DELTA);

					// at least 1 eigenvalue must be 1 (there's one direction that remains unchanged under rotation)
					double[][] ar = {{mT.m00,mT.m01,mT.m02},{mT.m10,mT.m11,mT.m12},{mT.m20,mT.m21,mT.m22}};
					Matrix mat = new Matrix(ar);
					double[] eigenv = mat.svd().getSingularValues();
					Assert.assertTrue(eq(eigenv[0],1.0) || eq(eigenv[1],1.0) || eq(eigenv[2],1.0));

					// transpose must be equals to inverse
					Matrix3d rotTransp = new Matrix3d();
					Matrix3d rotInv = new Matrix3d();
					rotTransp.transpose(rot);
					rotInv.invert(rot);
					
					assertMatrixEquals(rotTransp,rotInv);
				} 		


			}
		}
		
		Assert.assertEquals(266, allSGs.size()); // the total count must be 266
		Assert.assertEquals(65, countEn);        // enantiomorphic groups (protein crystallography groups)
		Assert.assertEquals(165, countNonEn);    // i.e. 266-65-36
		Assert.assertEquals(36, countSpecial);   // the rest of the groups present un symop.lib (sometimes used in PDB)
	}
	
	private static void assertMatrixEquals(Matrix3d m1, Matrix3d m2) {
		for (int i=0;i<3;i++) {
			for (int j=0;j<3;j++) {
				 Assert.assertEquals(m1.getElement(i, j),m2.getElement(i, j),DELTA); 
			}
		}
	}
	
	private static boolean eq(double d1, double d2) {
		return (Math.abs(d1-d2)<DELTA);
	}

	// to debug the testing code (run as java program so that we can use normal debugger)
	public static void main(String[] args) throws Exception {
		SpaceGroupTest test = new SpaceGroupTest();
		test.testTransfConversion();
	}

}
