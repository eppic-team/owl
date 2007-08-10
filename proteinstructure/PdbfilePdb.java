package proteinstructure;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.vecmath.Point3d;

public class PdbfilePdb extends Pdb {
	
	private static final String UNKNOWN_STRING ="Unknown";
	private static final String NULL_chainCode = "A";
	
	private String pdbfile;

	/**
	 * Constructs Pdb object reading from pdb file given pdb chain code. Model will be DEFAULT_MODEL
	 * @param pdbfile
	 * @param pdbChainCode
	 * @throws FileNotFoundException
	 * @throws IOException
	 * @throws PdbfileFormatError
	 * @throws PdbChainCodeNotFoundError 
	 */
	public PdbfilePdb (String pdbfile, String pdbChainCode) throws FileNotFoundException, IOException, PdbfileFormatError, PdbChainCodeNotFoundError {
		this(pdbfile,pdbChainCode,DEFAULT_MODEL);
	}
	
	/**
	 * Constructs Pdb object reading from pdb file given pdb chain code and model serial
	 * @param pdbfile
	 * @param pdbChainCode
	 * @param model_serial
	 * @throws FileNotFoundException
	 * @throws IOException
	 * @throws PdbfileFormatError
	 * @throws PdbChainCodeNotFoundError 
	 */
	public PdbfilePdb (String pdbfile, String pdbChainCode, int model_serial) throws FileNotFoundException, IOException, PdbfileFormatError, PdbChainCodeNotFoundError {
		this.pdbfile = pdbfile;
		this.model=model_serial;
		this.pdbCode=UNKNOWN_STRING; // we initialise to unknown in case we don't find it in pdb file 
		this.pdbChainCode=pdbChainCode.toUpperCase();	// our convention: chain codes are upper case
		// we set chainCode to pdbChainCode except for case NULL where we use "A"
		this.chainCode=pdbChainCode;
		if (pdbChainCode.equals("NULL")) this.chainCode=NULL_chainCode;
		
		this.sequence=""; // we initialise it to empty string, then is set in read_pdb_data_from_file 
		
		// we initialise the resser2secstruct and secstruct2resinterval Maps to empty, if no sec structure info is found then it remains empty
		this.resser2secstruct = new HashMap<Integer, String>();
		this.secstruct2resinterval = new TreeMap<String, Interval>();
		
		read_pdb_data_from_file();
		
		// when reading from pdb file we have no information of residue numbers or author's (original) pdb residue number, so we fill the mapping with the residue numbers we know
		//TODO eventually we could assign our own internal residue numbers when reading from pdb and thus this map would be used
		this.resser2pdbresser = new HashMap<Integer, String>();
		this.pdbresser2resser = new HashMap<String, Integer>();
		for (int resser:resser2restype.keySet()){
			resser2pdbresser.put(resser, String.valueOf(resser));
			pdbresser2resser.put(String.valueOf(resser), resser);
		}
		
		// initialising atomser2atom from resser_atom2atomserial
		atomser2atom = new HashMap<Integer, String>();
		for (String resser_atom:resser_atom2atomserial.keySet()){
			int atomserial = resser_atom2atomserial.get(resser_atom);
			String atom = resser_atom.split("_")[1];
			atomser2atom.put(atomserial,atom);
		}
	}

	
	
	/**
	 * To read the pdb data (atom coordinates, residue serials, atom serials) from file.
	 * chainCode gets set to same as pdbChainCode, except if input chain code NULL then chainCode will be 'A'
	 * pdbCode gets set to the one parsed in HEADER or to 'Unknown' if not found
	 * sequence gets set to the sequence read from ATOM lines (i.e. observed residues only)
	 * No insertion codes are parsed or taken into account at the moment. Thus files with 
	 * insertion codes will be incorrectly read
	 * @param pdbfile
	 * @throws FileNotFoundException
	 * @throws IOException
	 * @throws PdbfileFormatError
	 * @throws PdbChainCodeNotFoundError  
	 */
	private void read_pdb_data_from_file() throws FileNotFoundException, IOException, PdbfileFormatError, PdbChainCodeNotFoundError {
		resser_atom2atomserial = new HashMap<String,Integer>();
		resser2restype = new HashMap<Integer,String>();
		atomser2coord = new HashMap<Integer,Point3d>();
		atomser2resser = new HashMap<Integer,Integer>();
		Pattern p;
		Matcher m;
		boolean empty = true; // controls whether we don't find any atom line for given pdbChainCode and model
		// we set chainCodeStr (for regex) to pdbChainCode except for case NULL where we use " " (NULL is a blank chain code in pdb files)
		String chainCodeStr=pdbChainCode;
		if (pdbChainCode.equals("NULL")) chainCodeStr=" ";
		
		int thismodel=DEFAULT_MODEL; // we initialise to DEFAULT_MODEL, in case file doesn't have MODEL lines 
		BufferedReader fpdb = new BufferedReader(new FileReader(new File(pdbfile)));
		int linecount=0;
		String line;
		// read first line
		if((line = fpdb.readLine()) != null ) {
			linecount = 1;
			// HEADER
			p = Pattern.compile("^HEADER");
			m = p.matcher(line);
			if (m.find()){
				Pattern ph = Pattern.compile("^HEADER.{56}(\\d\\w{3})");
				Matcher mh = ph.matcher(line);
				if (mh.find()) {
					pdbCode=mh.group(1).toLowerCase();
				}
			} else { // header not found
				// check whether this is a Casp prediction file
				p = Pattern.compile("^PFRMAT TS");
				m = p.matcher(line);
				if(m.find()) {
					// ok, it is
					pdbCode = "CASP";
				} else {
					// a HEADER is the minimum we ask at the moment for a pdb file to have, if we don't find it in line 1 we throw an exception
					throw new PdbfileFormatError("The pdb file "+pdbfile+" doesn't seem to have the right format");
				}
			}
		} else {
			throw new PdbfileFormatError("The file "+pdbfile+" is empty.");
		}
		// read further lines
		while ((line = fpdb.readLine() ) != null ) {
			linecount++;
			// SEQRES
			//SEQRES   1 A  348  VAL ASN ILE LYS THR ASN PRO PHE LYS ALA VAL SER PHE
			p = Pattern.compile("^SEQRES.{5}"+chainCodeStr);
			m = p.matcher(line);
			if (m.find()){
				for (int i=19;i<=67;i+=4) {
					if (!line.substring(i, i+3).equals("   ")) {
						sequence+= AA.threeletter2oneletter(line.substring(i, i+3));
					}
				}
			}
			// SECONDARY STRUCTURE
			// helix
			//HELIX    1   1 LYS A   17  LEU A   26  1
			//							helix ser				beg res ser					end res ser
			p = Pattern.compile("^HELIX..(...).{9}"+chainCodeStr+".(....).{6}"+chainCodeStr+".(....)");
			m = p.matcher(line);
			if (m.find()){
				int serial = Integer.valueOf(m.group(1).trim());
				int beg = Integer.valueOf(m.group(2).trim());
				int end = Integer.valueOf(m.group(3).trim());
				String ssId = "H"+serial;
				secstruct2resinterval.put(ssId, new Interval(beg,end));
				for (int i=beg;i<=end;i++){
					if (resser2secstruct.containsKey(i)){// if already assigned we print a warning and then assign it
						//System.err.println("Inconsistency in secondary structure assignment. " +
						//		"Residue "+i+" is getting reassigned from "+resser2secstruct.get(i)+" to "+ssId);
					}
					resser2secstruct.put(i,ssId);
				}
			}
			// sheet
			//SHEET    2   A 5 ILE A  96  THR A  99 -1  N  LYS A  98   O  THR A 107
			//                       strand ser sheet id			 beg res ser                 end res ser
			p = Pattern.compile("^SHEET..(...).(...).{7}"+chainCodeStr+"(....).{6}"+chainCodeStr+"(....)");
			m = p.matcher(line);
			if (m.find()){
				int strandSerial = Integer.valueOf(m.group(1).trim());
				String sheetId = m.group(2).trim();
				int beg = Integer.valueOf(m.group(3).trim());
				int end = Integer.valueOf(m.group(4).trim());
				String ssId = "S"+sheetId+strandSerial;
				secstruct2resinterval.put(ssId, new Interval(beg,end));
				for (int i=beg;i<=end;i++){
					if (resser2secstruct.containsKey(i)){// if already assigned we print a warning and then assign it
						//System.err.println("Inconsistency in secondary structure assignment. " +
						//		"Residue "+i+" is getting reassigned from "+resser2secstruct.get(i)+" to "+ssId);
					}
					resser2secstruct.put(i,ssId);
				}
			}
			// we've stored the sec structure info in the strands2begEnd and sheets2strands maps.
			// the assignment to resser2secstruct is done when we reach the ATOM lines, see below
			// turn
			//TURN     1 S1A GLY A  16  GLN A  18     SURFACE
			//							turn ser				beg res ser					end res ser
			p = Pattern.compile("^TURN...(...).{9}"+chainCodeStr+"(....).{6}"+chainCodeStr+"(....)");
			m = p.matcher(line);
			if (m.find()){
				int serial = Integer.valueOf(m.group(1).trim());
				int beg = Integer.valueOf(m.group(2).trim());
				int end = Integer.valueOf(m.group(3).trim());
				String ssId = "T"+serial;
				secstruct2resinterval.put(ssId, new Interval(beg,end));
				for (int i=beg;i<=end;i++){
					if (resser2secstruct.containsKey(i)){// if already assigned we print a warning and then assign it
						//System.err.println("Inconsistency in secondary structure assignment. " +
						//		"Residue "+i+" is getting reassigned from "+resser2secstruct.get(i)+" to "+ssId);
					}
					resser2secstruct.put(i,ssId);
				}
			}			
			// MODEL
			p = Pattern.compile("^MODEL\\s+(\\d+)");
			m = p.matcher(line);
			if (m.find()){
				thismodel=Integer.parseInt(m.group(1));
			}
			if (thismodel!=model) continue; // we skip reading of atom lines if we are not in the desired model
			// ATOM
			p = Pattern.compile("^ATOM");
			m = p.matcher(line);
			if (m.find()){
				//                                 serial    atom   res_type      chain 	   res_ser     x     y     z
				Pattern pl = Pattern.compile(".{6}(.....).{2}(...).{1}(...).{1}"+chainCodeStr+"(.{4}).{4}(.{8})(.{8})(.{8})",Pattern.CASE_INSENSITIVE);
				Matcher ml = pl.matcher(line);
				if (ml.find()) {
					empty=false;
					int atomserial=Integer.parseInt(ml.group(1).trim());
					String atom = ml.group(2).trim();
					String res_type = ml.group(3).trim();
					int res_serial = Integer.parseInt(ml.group(4).trim());
					double x = Double.parseDouble(ml.group(5).trim());
					double y = Double.parseDouble(ml.group(6).trim());
					double z = Double.parseDouble(ml.group(7).trim());
					Point3d coords = new Point3d(x,y,z);
					ArrayList<String> aalist=AA.aas();
					if (aalist.contains(res_type)) {
						atomser2coord.put(atomserial, coords);
						atomser2resser.put(atomserial, res_serial);
						resser2restype.put(res_serial, res_type);
						ArrayList<String> atomlist = aas2atoms.get(res_type);
						if (atomlist.contains(atom)){
							resser_atom2atomserial.put(res_serial+"_"+atom, atomserial);
						}
					}					
				}
			}
		}
		fpdb.close();
		if (empty) {
			//System.err.println("Couldn't find any atom line for given pdbChainCode: "+pdbChainCode+", model: "+model);
			throw new PdbChainCodeNotFoundError("Couldn't find any ATOM line for given pdbChainCode: "+pdbChainCode+", model: "+model);
		}
		if (sequence.equals("")){
			// if we couldn't read anything from SEQRES then we read it from the resser2restype HashMap
			// NOTE: we must make sure elsewhere that there are no unobserved residues, we can't check that here!
			ArrayList<Integer> ressers = new ArrayList<Integer>();
			for (int resser:resser2restype.keySet()) {
				ressers.add(resser);
			}
			Collections.sort(ressers);
			for (int resser:ressers){
				String oneletter = AA.threeletter2oneletter(resser2restype.get(resser));
				sequence += oneletter;
			}
		}
	}
	
}
