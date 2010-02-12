package tests.proteinstructure;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.TreeSet;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import tools.Interval;

public class IntervalTest {

	private static final String[] TP = {"1","10","10,1","10-9","10-10,5","1,6-7,8-9"};
	private static final ArrayList<TreeSet<Integer>> TP_SETS = initialiseTpSets();
	private static final String[] TN = {"", "-",",","1,1,","1,-2","10-9-7",",4","4--5","1.5","1,,2"," "};

	private static ArrayList<TreeSet<Integer>> initialiseTpSets() {
		ArrayList<TreeSet<Integer>> sets = new ArrayList<TreeSet<Integer>>();
		sets.add(new TreeSet<Integer>(Arrays.asList(1)));
		sets.add(new TreeSet<Integer>(Arrays.asList(10)));
		sets.add(new TreeSet<Integer>(Arrays.asList(1,10)));
		sets.add(new TreeSet<Integer>());
		sets.add(new TreeSet<Integer>(Arrays.asList(5,10)));
		sets.add(new TreeSet<Integer>(Arrays.asList(1,6,7,8,9)));
		return sets;
	}
	
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
	public void testIsValidSelectionString() {
		
		System.out.println("Positives:");
		for(String selStr:TP) {
			assertTrue(Interval.isValidSelectionString(selStr));
		}
		System.out.println("Negatives:");
		for(String selStr:TN) {
			assertTrue(!Interval.isValidSelectionString(selStr));
		}
		
	}
	
	@Test
	public void testParseSelectionString() {
		System.out.println("Parsing:");
		for(int i=0;i<TP.length;i++) {
			System.out.print(TP[i] + "\t\t");
			TreeSet<Integer> nodeSet = Interval.parseSelectionString(TP[i]);
			assertNotNull(nodeSet);
			System.out.println(nodeSet);
			for (int member:nodeSet) {
				assertTrue(TP_SETS.get(i).contains(member));
			}
			for (int member:TP_SETS.get(i)) {
				assertTrue(nodeSet.contains(member));
			}
		}
		
	}
}
