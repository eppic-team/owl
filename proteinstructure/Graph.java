package proteinstructure;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.io.IOException;
import java.util.ArrayList;


public class Graph {

	ArrayList<Contact> contacts;
	double cutoff;
	String ct;
	
	public Graph (ArrayList<Contact> contacts,double cutoff,String ct) {
		this.contacts=contacts;
		this.cutoff=cutoff;
		this.ct=ct;
	}
	
	public void write_contacts_to_file (String outfile) throws IOException {
		PrintStream Out = new PrintStream(new FileOutputStream(outfile));
		for (Contact pair:contacts){
			int i_resser=pair.i;
			int j_resser=pair.j;
			Out.println(i_resser+"\t"+j_resser);
		}
		Out.close();		
	}
}
