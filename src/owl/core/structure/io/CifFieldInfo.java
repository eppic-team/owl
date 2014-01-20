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
	public StringBuilder getLastSubFieldData() {
		return this.subFieldsData.get(this.subFieldsData.size()-1);
	}
	
	public StringBuilder getSubFieldData(int index) {
		return this.subFieldsData.get(index);
	}
	
	public boolean isLoop() {
		return isLoop;
	}

	public void setLoop(boolean isLoop) {
		this.isLoop = isLoop;
	}

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
}
 