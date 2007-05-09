package proteinstructure;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.HashMap;
import java.util.ArrayList;
import tools.MySQLConnection;

public class AA {
	
	public final static String MYSQLSERVER="white";
	public final static String MYSQLUSER=getUserName();
	public final static String MYSQLPWD="nieve";
	public final static String DB = "aa_info";

	//public static final ArrayList<String> AAS = aas();
	//public static final HashMap<String,String> THREELETTER2ONELETTER = threeletter2oneletter();

	/** get user name from operating system (for use as database username) */
	private static String getUserName() {
		String user = null;
		user = System.getProperty("user.name");
		if(user == null) {
			System.err.println("Could not get user name from operating system. Exiting");
			System.exit(1);
		}
		return user;
	}
	
	public static HashMap<String,String> getThreeletter2oneletter() {
		HashMap<String,String> three2oneletter = new HashMap<String,String>();
		three2oneletter.put("CYS", "C");
		three2oneletter.put("ASP", "D");
		three2oneletter.put("SER", "S");
		three2oneletter.put("GLN", "Q");
		three2oneletter.put("LYS", "K");
		three2oneletter.put("ILE", "I");
		three2oneletter.put("PRO", "P");
		three2oneletter.put("THR", "T");
		three2oneletter.put("PHE", "F");
		three2oneletter.put("ALA", "A");
		three2oneletter.put("GLY", "G");
		three2oneletter.put("HIS", "H");
		three2oneletter.put("GLU", "E");
		three2oneletter.put("LEU", "L");
		three2oneletter.put("ARG", "R");
		three2oneletter.put("TRP", "W");
		three2oneletter.put("VAL", "V");
		three2oneletter.put("ASN", "N");
		three2oneletter.put("TYR", "Y") ;
		three2oneletter.put("MET", "M");
		return three2oneletter;
	}
	
	public static String threeletter2oneletter(String three) {
		HashMap<String,String> three2oneletter = getThreeletter2oneletter();
		return three2oneletter.get(three);
	}
	
	public static ArrayList<String> aas() {
		HashMap<String,String> three2oneletter = getThreeletter2oneletter();
		ArrayList<String> aas = new ArrayList<String>();
		for (String aa:three2oneletter.keySet()) {
			aas.add(aa);
		}
		return aas;
	}
	
	public static HashMap<String,ArrayList<String>> getaas2atoms() {
		HashMap<String,ArrayList<String>> aas2atoms = new HashMap<String,ArrayList<String>>();
		MySQLConnection conn = new MySQLConnection(MYSQLSERVER,MYSQLUSER,MYSQLPWD,DB);
		String sql="SELECT code_3_letter, group_concat(chem_atom_name) FROM atom_names GROUP BY code_3_letter ORDER BY code_3_letter";
		try {
			Statement stmt = conn.createStatement();
			ResultSet rsst = stmt.executeQuery(sql);
			while (rsst.next()) {
				String res_type=rsst.getString(1);
				String atomsStr = rsst.getString(2);
				String[] atoms = atomsStr.split(",");
				ArrayList<String> atomsAL = new ArrayList<String>();
				for (String atom:atoms){
					if (! atom.equals("OXT")){
						atomsAL.add(atom);
					}
				}
				aas2atoms.put(res_type, atomsAL);	
			} 
			rsst.close();
			stmt.close();
		} catch (SQLException e) {
			e.printStackTrace();
		}
		conn.close();
		return aas2atoms;
		
	}
	
	public static HashMap<String,String[]> ct2atoms(String ct) {
		ArrayList<String> aas = aas();
		HashMap<String,String[]> ct2atoms = new HashMap<String,String[]>();
		if (ct.equals("Ca")){
			String[] atoms = {"CA"};
			for (String aa:aas) {
				ct2atoms.put(aa, atoms);
			}
		} 
		else if (ct.equals("Cb")){
			String[] atoms = {"CB"};
			for (String aa:aas) {
				ct2atoms.put(aa, atoms);
			}
			atoms = new String[1];
			atoms[0]="CA";
			ct2atoms.put("GLY", atoms);
		}
		else if (ct.equals("C")){
			String[] atoms = {"C"};
			for (String aa:aas) {				
				ct2atoms.put(aa, atoms);
			}			
		}
		else if (ct.equals("Cg")){
			String[] atoms = {"SG"};
			ct2atoms.put("CYS", atoms);
			atoms = new String[1];
			atoms[0]= "CG";
			ct2atoms.put("ASP", atoms);
			atoms = new String[1];
			atoms[0]="OG";
			ct2atoms.put("SER", atoms);
			atoms = new String[1];
			atoms[0]="CG";
			ct2atoms.put("GLN", atoms);
			atoms = new String[1];
			atoms[0]="CD";
			ct2atoms.put("LYS", atoms);
			atoms = new String[1];
			atoms[0]="CG1";
			ct2atoms.put("ILE", atoms);
			atoms = new String[1];
			atoms[0]="CG";
			ct2atoms.put("PRO", atoms);
			atoms = new String[1];
			atoms[0]="OG1";
			ct2atoms.put("THR", atoms);
			atoms = new String[1];
			atoms[0]="CG";
			ct2atoms.put("PHE", atoms);
			ct2atoms.put("ALA", new String[0]);
			ct2atoms.put("GLY", new String[0]);
			atoms = new String[1];
			atoms[0]="CG";
			ct2atoms.put("HIS", atoms);
			atoms = new String[1];
			atoms[0]="CG";
			ct2atoms.put("GLU", atoms);
			atoms = new String[1];
			atoms[0]="CG";
			ct2atoms.put("LEU", atoms);
			atoms = new String[1];
			atoms[0]="NE";
			ct2atoms.put("ARG", atoms);
			atoms = new String[1];
			atoms[0]="CD2";
			ct2atoms.put("TRP", atoms);
			atoms = new String[1];
			atoms[0]="CG1";
			ct2atoms.put("VAL", atoms);
			atoms = new String[1];
			atoms[0]="CG";
			ct2atoms.put("ASN", atoms);
			atoms = new String[1];
			atoms[0]="CG";
			ct2atoms.put("TYR", atoms);
			atoms = new String[1];
			atoms[0]="CG";
			ct2atoms.put("MET", atoms);			
		}
		else if (ct.equals("BB")){
			String[] atoms = {"CA", "N", "C", "O"};
			for (String aa:aas) {
				ct2atoms.put(aa, atoms);
			}
		}
		else if (ct.equals("SC")){
			HashMap<String,ArrayList<String>> aas2atoms = getaas2atoms();			
			for (String aa:aas) {
				ArrayList<String> SCatoms =aas2atoms.get(aa);
				SCatoms.remove("CA");
				SCatoms.remove("N");
				SCatoms.remove("C");
				SCatoms.remove("O");
				String[] SCatomsar= new String[SCatoms.size()];
				SCatomsar=SCatoms.toArray(SCatomsar);
				ct2atoms.put(aa, SCatomsar);
			}
		}
		else if (ct.equals("ALL")){
			HashMap<String,ArrayList<String>> aas2atoms = getaas2atoms();
			for (String aa:aas) {
				ArrayList<String> ALLatoms = aas2atoms.get(aa);
				String[] ALLatomsar= new String[ALLatoms.size()];
				ALLatomsar=ALLatoms.toArray(ALLatomsar);				
				ct2atoms.put(aa,ALLatomsar);
			}
		}		
		return ct2atoms;
	}
}
