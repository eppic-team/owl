package proteinstructure;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import sequence.BlastHit;
import sequence.BlastHitList;
import sequence.GTGHit;
import sequence.GTGHitList;
import tools.MySQLConnection;

/**
 * A set of Templates to be used for homology modelling
 * See also RIGEnsemble.
 */
public class TemplateList implements Iterable<Template> {

	public static final String IDS_REGEX1 = "^(\\d\\w\\w\\w)(\\w)";
	public static final String IDS_REGEX2 = "^(\\d\\w\\w\\w)_(\\w)";
	public static final String IDS_REGEX3 = "^(\\d\\w\\w\\w)\\s+(\\w)";
	
	private ArrayList<Template> list;
	private boolean pdbDataLoaded;

	/**
	 * Constructs a new empty TemplateList
	 */
	public TemplateList() {
		list = new ArrayList<Template>();
		pdbDataLoaded = false;
	}
	
	/**
	 * Constructs a new TemplateList given a BlastHitList
	 * @param hits
	 * @param conn
	 */
	public TemplateList(BlastHitList hits) {
		list = new ArrayList<Template>();
		Iterator<BlastHit> it = hits.iterator();
		while (it.hasNext()) {
			this.add(new Template(it.next()));
		}
		pdbDataLoaded = false;
	}

	/**
	 * Constructs a new TemplateList given a GTGHitList 
	 * @param hits
	 */
	public TemplateList(GTGHitList hits) {
		list = new ArrayList<Template>();
		Iterator<GTGHit> it = hits.iterator();
		while (it.hasNext()) {
			this.add(new Template(it.next()));
		}
		pdbDataLoaded = false;
	}
	
	/**
	 * Constructs a TemplateList from a file with ids in the format pdbCode+chain, e.g. 1abcA
	 * @param idsFile
	 * @throws IOException
	 */
	public TemplateList(File idsFile) throws IOException {
		list = new ArrayList<Template>();
		readIdsFromFile(idsFile);
		pdbDataLoaded = false;
	}
	
	/**
	 * Constructs a TemplateList given an array of ids in the format pdbCode+chain, e.g. 1abcA
	 * @param ids
	 */
	public TemplateList(String[] ids) {
		list = new ArrayList<Template>();
		for (String id:ids) {
			this.add(new Template(id));
		}
		pdbDataLoaded = false;
	}
	
	/**
	 * Reads a file with a list of ids, creating Templates for them and adding them to this list
	 * Lines starting with # will be considered as comments
	 * @param idsFile
	 * @throws IOException
	 */
	private void readIdsFromFile(File idsFile) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(idsFile));
		String line;
		while ((line=br.readLine())!=null) {
			if (line.length()==0) continue;
			if (line.startsWith("#")) continue;
			Pattern p = Pattern.compile("\\d\\w\\w\\w\\w");
			Matcher m = p.matcher(line);
			if (m.matches()) {
				this.add(new Template(line));
			}
		}
		br.close();
	}
	
	/**
	 * Adds a Template to this list
	 * @param template
	 */
	public void add(Template template) {
		list.add(template);
	}
	
	/**
	 * Returns the numbers of templates in this list
	 * @return
	 */
	public int size() {
		return list.size();
	}
	
	/**
	 * Prints some info about each of the templates in this list  
	 */
	public void print() {
		for (Template template:list) {
			template.print();
		}
	}
	
	/**
	 * Returns an iterator over this list
	 * @return
	 */
	public Iterator<Template> iterator() {
		return this.list.iterator();
	}
	
	/**
	 * Returns an array with all ids of the templates in this list. IDs are in the format pdbCode+chain, e.g. 1abcA 
	 * @return
	 */
	public String[] getIds() {
		String[] ids = new String[this.size()];
		for (int i=0;i<this.size();i++) {
			ids[i] = list.get(i).getId();
		}
		return ids;
	}
	
	/**
	 * Tells whether this TemplateList contains the Template represented by the given templateId
	 * @param templateId
	 * @return
	 */
	public boolean contains(String templateId) {
		for (Template template:this) {
			if (template.getId().equals(templateId)) return true;
		}
		return false;
	}
	
	/**
	 * Gets a Template given its templateId or null if Template not present in this list
	 * @param id
	 * @return
	 */
	public Template getTemplate(String templateId) {
		for (Template template:this) {
			if (template.getId().equals(templateId)) return template;
		}
		return null;
	}
	
	/**
	 * Gets the rank (from 1 to max) of the given templateId in this TemplateList
	 * or 0 if templateId is not present in list.
	 * @param templateId
	 * @return
	 */
	public int getRank(String templateId) {
		for (int i=0;i<this.size();i++) {
			if (list.get(i).getId().equals(templateId)) return i+1;
		}
		return Integer.MAX_VALUE; // for sorting we want a null rank to be the MAX integer so it sorts at the end
	}
	
	/**
	 * Compare this TemplatesList agains the given one, printing out the ids of the common templates.
	 * Note: also the templates unique to each list are found, but we don't output them yet 
	 * @param otherTemplateList
	 */
	public void compare(TemplateList otherTemplateList) {
		ArrayList<String> common = new ArrayList<String>();
		ArrayList<String> uniqueToThis = new ArrayList<String>();
		//ArrayList<String> uniqueToOther = new ArrayList<String>();
		for (Template template: this) {
			if (otherTemplateList.contains(template.getId())) {
				common.add(template.getId());
			} else {
				uniqueToThis.add(template.getId());
			}
		}
		//for (Template template: otherTemplateList) {
		//	if (!this.contains(template.getId())) {
		//		uniqueToOther.add(template.getId());
		//	}
		//}
		System.out.println(common.size()+" common templates: ");
		for (String id:common) {
			System.out.print(id+ " ");
		}
		System.out.println();
		//System.out.println(common.size()+" unique to templates: ");
		//System.out.println(common.size()+" unique to templates: ");
	}
	
	/**
	 * Writes the ids of all Templates in this list to given File
	 * @param outFile
	 * @throws IOException
	 */
	public void writeIdsToFile(File outFile) throws IOException {
		PrintWriter Out = new PrintWriter(outFile);
		for (Template template:this) {
			Out.println(template.getId());
		}
		Out.close();
	}
	
	/**
	 * Loads the PDB data (structures) for this TemplateList
	 * Tries to load all data and throws exception if something goes wrong. 
	 * If there's no PDB code for a certain Template, then we simply print a warning 
	 * and skip it. We consider subsequently that data were loaded.
	 * @param conn
	 * @param pdbaseDb
	 * @throws PdbLoadError
	 * @throws SQLException
	 */
	public void loadPdbData(MySQLConnection conn, String pdbaseDb) throws PdbLoadError, SQLException {
		for (Template template:this) {
			try {
				template.loadPdbData(conn, pdbaseDb);
			} catch (PdbCodeNotFoundError e) {
				System.err.println("Warning: no PDB data for template "+template.getId()+". Error: "+e.getMessage());
			}
		}
		pdbDataLoaded = true;
	}
	
	/**
	 * Tells whether there is PDB data loaded for Templates in this TemplateList, 
	 * i.e. if loadPdbData has been called. 
	 * @return
	 */
	public boolean isPdbDataLoaded() {
		return pdbDataLoaded;
	}
	
	/*--------------------------------- static methods ----------------------------------*/
	
	/**
	 * Reads a list file containing a list of pdb codes and chain codes in 3 possible formats:
	 * - 1 column pdbCodes+chainCodes, e.g. 1bxyA
	 * - 1 column underscore-separated pdbCodes and chainCodes, e.g. 1bxy_A
	 * - 2 colums tab/spaces-separated pdbCodes and chainCodes, e.g. 1bxy A or 1bxy   A
	 * See the IDS_REGEX constants of this class for the regex that we are using.
	 * A mix of the formats is also tolerated.
	 * Chain codes can only be a 1 letter code (so we must use an "A" for NULL codes)
	 * @param listFile
	 * @return an array of pdbCodes(lower case)+chainCodes(conserving case) in the format 1bxyA
	 * @throws IOException
	 */
	public static String[] readIdsListFile(File listFile) throws IOException {
		ArrayList<String> codesAL = new ArrayList<String>(); 

		BufferedReader fileIn = new BufferedReader(new FileReader(listFile));
		String line;
		int lineCount=0;
		while((line = fileIn.readLine()) != null) {
			lineCount++;
			if (line.length()!=0 && !line.startsWith("#")) {
				Pattern p1 = Pattern.compile(IDS_REGEX1);
				Matcher m1 = p1.matcher(line);
				Pattern p2 = Pattern.compile(IDS_REGEX2);
				Matcher m2 = p2.matcher(line);
				Pattern p3 = Pattern.compile(IDS_REGEX3);
				Matcher m3 = p3.matcher(line);				

				if (m1.matches()) {
					codesAL.add(m1.group(1).toLowerCase()+m1.group(2));
				} 
				else if (m2.matches()) {
					codesAL.add(m2.group(1).toLowerCase()+m2.group(2));
				} 
				else if (m3.matches()){
					codesAL.add(m3.group(1).toLowerCase()+m3.group(2));
				}
				else {
					System.err.println("Line "+lineCount+" in list file "+listFile+" is not in any of the recognised pdbCode+chainCode formats");
				}
			
			}
		}
		String[] codes = new String[codesAL.size()];
		codesAL.toArray(codes);
		return codes;
	}
	
	/**
	 * Filters out from a given list of PDB ids (in format pdbCode+pdbChainCode, e.g. 1abcA)
	 * the ones that were released after the given date (in format yyyymmdd, e.g. 20060425) 
	 * @param conn
	 * @param pdbaseDb
	 * @param ids
	 * @param date date in format yyyymmdd
	 * @return a new array containing the filtered out list of PDB ids
	 * @throws SQLException
	 */
	public static String[] filterIdsByMaxReleaseDate(MySQLConnection conn, String pdbaseDb, String[] ids, String date) throws SQLException {
		Statement stmt = conn.createStatement();
		String tempTable = pdbaseDb +".temp" + System.currentTimeMillis(); // hopefully a unique name (although with temporary tables not really needed)
		String sql = "CREATE TEMPORARY TABLE "+tempTable+" (id char(5))";
		stmt.executeUpdate(sql);
		for (String id:ids) {
			sql = "INSERT INTO "+tempTable+" (id) VALUES ('"+id+"')";
			stmt.executeUpdate(sql);
		}
		sql = "SELECT DISTINCT a.id " +
				" FROM "+tempTable+" AS a " +
				 " JOIN "+pdbaseDb+".struct AS b " +
				 " JOIN "+pdbaseDb+".database_pdb_rev AS c" +
				 " ON (upper(substring(id,1,4))=b.entry_id AND b.entry_key=c.entry_key)" +
				" WHERE CAST(replace(c.date2,'-','') AS UNSIGNED) >= "+date;
		ResultSet rsst = stmt.executeQuery(sql);
		ArrayList<String> filteredIdsAL = new ArrayList<String>();
		while (rsst.next()) {
			filteredIdsAL.add(rsst.getString(1));
		}
		
		rsst.close();
		
		stmt.executeUpdate("DROP TABLE "+tempTable);
		stmt.close();
		
		String[] filteredIds = new String[filteredIdsAL.size()];
		filteredIdsAL.toArray(filteredIds);
		return filteredIds;
	}

}
