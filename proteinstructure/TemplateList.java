package proteinstructure;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import sequence.BlastHit;
import sequence.BlastHitList;

/**
 * A set of Templates to be used for homology modelling
 *
 */
public class TemplateList {

	ArrayList<Template> list = new ArrayList<Template>();

	/**
	 * Constructs a new empty TemplateList
	 */
	public TemplateList() {
		list = new ArrayList<Template>();
	}
	
	/**
	 * Constructs a TemplateList given a BlastHitList
	 * @param hits
	 */
	public TemplateList(BlastHitList hits) {
		list = new ArrayList<Template>();
		Iterator<BlastHit> it = hits.iterator();
		while (it.hasNext()) {
			this.add(new Template(it.next()));
		}
	}
	
	/**
	 * Constructs a TemplateList from a file with ids in the format pdbCode+chain, e.g. 1abcA
	 * @param idsFile
	 * @throws IOException
	 */
	public TemplateList(File idsFile) throws IOException {
		list = new ArrayList<Template>();
		readIdsFromFile(idsFile);
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
	}
	
	/**
	 * Reads a file with a list of ids, creating Templates for them and adding them to this list
	 * @param idsFile
	 * @throws IOException
	 */
	private void readIdsFromFile(File idsFile) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(idsFile));
		String line;
		while ((line=br.readLine())!=null) {
			if (line.length()==0) continue;
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
	
}
