package owl.core.structure.io;

import java.util.ArrayList;

public class CifFieldInfo {

	/**
	 * The field name, e.g. _struct
	 */
	private String fieldName;
	
	/**
	 * Each of the subfields, e.g. for _struct subfields are entry_id, title, etc
	 */
	private ArrayList<String> subFields; 
	
	/**
	 * Data present in each subfield (in same line or multiline), 
	 * if loop field then only the last subfield will contain the
	 * whole of the data.
	 * Note we use StringBuilder in order to do appending as fast 
	 * as possible (StringBuffer is thread-safe, but we don't need that
	 * here. Since java 1.5 StringBuilder is available, not thread-safe 
	 * but faster) 
	 */
	private ArrayList<StringBuilder> subFieldsData;

	/**
	 * Flags loop fields
	 */
	private boolean isLoop;
	
	/**
	 * The index of the data stream being tokenised. Points to current character while 
	 * tokenising data
	 */
	private int dataIdx;
	
	
	public CifFieldInfo (String fieldName) {
		this.subFields = new ArrayList<String>();
		this.subFieldsData = new ArrayList<StringBuilder>();
		this.fieldName = fieldName;
	}

	public String getFieldName() {
		return fieldName;
	}

	public void setFieldName(String fieldName) {
		this.fieldName = fieldName;
	}

	public void addSubField(String subField) {
		this.subFields.add(subField);
		this.subFieldsData.add(new StringBuilder());
	}
	
	/**
	 * Adds data to last subfield, to be used in loop cases where
	 * all data is placed together only after the last subfield header
	 * @param data
	 */
	public void addDataToLastSubFieldDataBuffer(String data) {
		this.subFieldsData.get(subFields.size()-1).append(data);
	}
	
	public int getNumSubFields() {
		return this.subFields.size();
	}
	
	public int getIndexForSubField(String subField) {
		for (int i=0;i<getNumSubFields();i++) {
			if (subFields.get(i).equals(subField)) return i;
		}
		return -1;
	}
	
	public boolean isSubFieldPresent(String subField) {
		int index = getIndexForSubField(subField);
		if (index<0) return false;
		return true;
	}
	
	public String getSubFieldData(String subField) {
		return this.subFieldsData.get(getIndexForSubField(subField)).toString();
	}
	
	/**
	 * Gets the data from the last subfield, to be used in loop cases where 
	 * all the data is stored in last subfield only.
	 * @return
	 */
	private StringBuilder getLastSubFieldData() {
		return this.subFieldsData.get(this.subFieldsData.size()-1);
	}
	
	public StringBuilder getSubFieldData(int index) {
		return this.subFieldsData.get(index);
	}
	
	public boolean isLoop() {
		return isLoop;
	}

	public void setLoop(boolean isLoop) {
		// once we know it is a loop we initialise the dataIdx to 1
		// we initalise to 1 because the parsing/tokenising needs to start from the char after the first '\n'
		this.dataIdx = 1;
		this.isLoop = isLoop;
	}

	/**
	 * Tells whether the field is not present at all in the scanned mmCIF file
	 * @return
	 */
	public boolean isEmpty() {
		return this.subFields.isEmpty();
	}
	
	/**
	 * Checks consistency of the data structure read:
	 * - subFields and subFieldsData length coincide
	 * - in loop cases, all data stored in last subField
	 * 
	 * @return true if data structure consistent, false otherwise 
	 */
	public boolean checkConsistency() {
		if (subFields.size()!=subFieldsData.size()) {
			return false;
		}
		if (isLoop) {
			for (int i=0;i<subFieldsData.size();i++) {
				// in loops all subFieldsData are empty (with only a new line char) except last one
				if ((i!=subFieldsData.size()-1) && 
						subFieldsData.get(i).length()>1) {
					return false;
				}
				
			}
		}
		return true;
	}
	
	/**
	 * Returns the next tokens for next record in a loop field. Use together
	 * with {@link #hasMoreData()}
	 * @return
	 */
	public String[] getNextTokens() {
		if (!isLoop()) throw new IllegalArgumentException("Called loop method in a non-loop field");
		return tokeniseData(getLastSubFieldData(), getNumSubFields());
	}
	
	/**
	 * In a loop field it will return true while more records are present.
	 * Use together with {@link #getNextTokens()} to go through all records
	 * of a loop field.
	 * @param data
	 * @return
	 */
	public boolean hasMoreData() {
		if (!isLoop()) throw new IllegalArgumentException("Called loop method in a non-loop field");
		return (dataIdx < getLastSubFieldData().length()-1);
	}
	
	/**
	 * Processes the data present in a subField and returns a trimmed String.
	 * The data can be either one-line or multi-line (;; quoted or non-quoted)
	 * @param subField
	 * @return
	 */
	public String getFinalSubFieldData(String subField) {
		
		if (isLoop) throw new IllegalArgumentException("Called non-loop method in a loop field");
		
		int index = getIndexForSubField(subField);
		
		String finalStr = null;
		
		// data can be single line or multiline, we can catch it by checking first character: if multiline 1st char is a '\n'
		
		if (getSubFieldData(index).charAt(0)!='\n'){ // single line

			finalStr = getSubFieldData(index).toString().trim();

		} else { // multi line

			String[] tokens = tokeniseSubFieldData(index);
			
			//TODO should we trim here too?
			finalStr = tokens[0];

		}
		
		// TODO we should remove quotes here: '' "" and maybe also ;; (check in CiffileParser)
		
		return finalStr;
	}
	
	private String[] tokeniseSubFieldData(int index) {
		dataIdx = 1;
		return tokeniseData(this.subFieldsData.get(index), 1);
	}
	
	/**
	 * Splits a cif data string into its individual tokens returning a String array with all tokens
	 * Takes care of all particularities of the format of data in the ciffiles:
	 *  - fields within records are separated by spaces
	 *  - spaces can be used within quoted strings (with single or double quotes)
	 *  - free style with all characters allowed if something is quoted with \n; ;\n 
	 * 
	 * The given data StringBuilder must start with '\n'
	 * The global variable dataIdx must be initialized to 1 (in order to start at the 
	 * character after the first '\n') before starting a read
	 * 
	 * The java class StreamTokenizer could have done all this, but it was limited to do all that we needed to do
	 * 
	 * @param data
	 * @param numberTokens
	 * @return
	 */
	private String[] tokeniseData(StringBuilder data, int numberTokens) {
		
		//TODO there's no error handling within this tokeniser,
		//     what if something goes wrong? how can we catch errors? can we catch for instance arrays out of bounds?
		
		String[] tokens = new String[numberTokens];
		// initialise tokens to empty strings
		for (int i=0; i<numberTokens;i++){
			tokens[i]="";
		}
		
		int i = 0; // token index
		char lastChar = 0; // ' '
		char quoteChar = 0;
		while (true) {
			
			char currentChar = data.charAt(dataIdx++);
			
			// '' quoting
			if (quoteChar!=';' && currentChar=='\'' && (lastChar==' ' || lastChar=='\n' || lastChar==0)){
				quoteChar = '\'';
			}
			else if (quoteChar!=';' && currentChar==' ' && lastChar=='\''){
				quoteChar = 0;
			}
			// "" quoting
			if (quoteChar!=';' && currentChar=='"' && (lastChar==' ' || lastChar=='\n' || lastChar==0)){
				quoteChar = '"';
			}
			else if (quoteChar!=';' && currentChar==' ' && lastChar=='"'){
				quoteChar = 0;
			}			
			// ;; quoting (multi-line quoting)
			if (quoteChar!=';' && currentChar==';' && (lastChar=='\n' || lastChar==0)){ 
				quoteChar = ';';
			}
			else if (quoteChar==';' && currentChar==';' && lastChar=='\n'){
				quoteChar = 0;
			}
			
			// reading field
			if (quoteChar==0) { // not within quotes
				if (currentChar==' ' || currentChar=='\n') { 
					if (currentChar!=lastChar && !(currentChar=='\n' && lastChar==' ')) i++; // we only increment when we move from a non-space to a space or from non-space to \n
				} else {
					tokens[i]+=currentChar;
					// if we are adding the last ; of a ;;-quoted string then strip the starting ';' and ending "\n;" out 
					if (currentChar==';' && lastChar=='\n' && tokens[i].startsWith(";") && tokens[i].endsWith("\n;")) {
						tokens[i]=tokens[i].replaceFirst("^;", "");
						tokens[i]=tokens[i].replaceFirst("\n;","");
					}									
				} 
			} else {			// within quotes (of type '', ""  or  ;;)
				tokens[i]+=currentChar;
				// if string is surrounded by '' or "" then strip them out (except when string is length 1 and thus beginning and end are quoteChar)
				if (tokens[i].length()!=1 && 
						tokens[i].startsWith(Character.toString(quoteChar)) && 
						tokens[i].endsWith(Character.toString(quoteChar))) 
					tokens[i]=tokens[i].replaceAll(Character.toString(quoteChar), "");

			}
			
			lastChar = currentChar;
			 
			if (i==numberTokens) {
				// for the last record of an element it is important to have read up to the end of the line (including the '\n'), 
				// This is needed because:
				// 1) next iteration needs to start on the character after the '\n'
				// 2) condition "while (dataIdx<data string length)-1" 
				while (true) {
					// first check if we are already at the end, we can't go further: in this case
					// we'll stay at the '\n', that's why the condition "while (dataIdx<data string length)-1" has the -1					
					if (dataIdx>=data.length()-1) break;
					
					currentChar = data.charAt(dataIdx++);
					if (currentChar=='\n'){ 
						break;
					}
					
				}
				return tokens;
			}
			
		}
	}

	
}
 