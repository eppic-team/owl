package casp.tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class splitGTGFile {

	public static void main(String[] args) throws Exception {

		File file = new File(args[0]);
		
		BufferedReader br = new BufferedReader(new FileReader(file));

		String line;
		int count=0;
		String target = null;
		String lastTarget = "";
		
		PrintWriter pw = null;
		
		while ((line=br.readLine())!=null) {
			if (line.equals("")) continue;
			if (line.contains("<PRE>")) continue;
			Pattern p = Pattern.compile("# target=([Tt]\\d\\d\\d\\d)");
			Matcher m = p.matcher(line);
			if (m.find()) {
				target = m.group(1);
				if (!lastTarget.equals(target)) {
					if (pw!=null) pw.close();
					count++;
					pw = new PrintWriter(new File(target+".gtg"));					
				}
				
				lastTarget = target;
			}
			pw.println(line);
		}		
		
		br.close();
	}

}
