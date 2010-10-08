package owl.tests.core.structure;


import javax.vecmath.Matrix4d;

import junit.framework.Assert;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

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
		
		for (SpaceGroup spaceGroup:SymoplibParser.getAllSpaceGroups().values()) {
			Assert.assertEquals(SpaceGroup.getMatrixFromAlgebraic(spaceGroup.getTransfAlgebraic(0)),spaceGroup.getTransformation(0));
			int i =1;
			for (Matrix4d transf:spaceGroup.getTransformations()){
				Assert.assertEquals(SpaceGroup.getAlgebraicFromMatrix(transf), 
						SpaceGroup.getAlgebraicFromMatrix(SpaceGroup.getMatrixFromAlgebraic(SpaceGroup.getAlgebraicFromMatrix(transf))));
				
				Assert.assertEquals(SpaceGroup.getMatrixFromAlgebraic(spaceGroup.getTransfAlgebraic(i)),transf);
				i++;
			}
			
		}
	}

}
