import java.util.regex.Matcher;
import java.util.regex.Pattern;


public class testRegex {

	private static final Pattern REGEX = Pattern.compile("^.*\\|.*\\|([^. ]+)(?:.\\d+)?\\s.*$"); //
	
	public static void main(String[] args) {
		
		String string = ">emblcds|X90765.2|CAA62291.3 Thermus aquaticus ribosomal protein L30 ";
		
		Matcher m = REGEX.matcher(string);
		if (m.matches()) {
			System.out.println("Match!");
			System.out.println(m.group(1));
		}
	
	}

}
