package proteinstructure;

public enum EdgeState { 
	UNKNOWN, NONCONTACT, CONTACT, SKIPPED, SKIPPED_GIVEN, CN_NOT_IN_BG, NOCN, ORIG_NONCONTACT_GIVEN, ORIG_CONTACT_GIVEN;
    public boolean contact() {
    	boolean contact = false;
    	if (this.equals(EdgeState.CONTACT) || this.equals(EdgeState.ORIG_CONTACT_GIVEN)) {
    		contact = true;
    	}
    	return contact; 
    }	
}	