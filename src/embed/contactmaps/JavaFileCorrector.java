package embed.contactmaps;

import java.io.*;

public class JavaFileCorrector {
	
	private static final String emb = "embed";
	private static final String pro = "proteinstructure";
	private static final String too = "tool";
	
	private String content;
	private String file_name;
	private String dir;
	
	public JavaFileCorrector (String dir_name, String file_name) throws IOException{
		File file = new File (dir_name+file_name);
		if(file.exists() && file.isFile()){
			BufferedReader reader = new BufferedReader(new FileReader(file));
			String linereader = "", content = "";
			while((linereader = reader.readLine())!=null){
				if(linereader.contains("src."+emb)) content += linereader+"\n";
				else{
					if(linereader.contains(emb)) content += "src."+linereader+"\n";
					else content += linereader+"\n";
				}
				if(linereader.contains("src."+pro)) content += linereader+"\n";
				else{
					if(linereader.contains(pro)) content += "src."+linereader+"\n";
					else content += linereader+"\n";
				}
				if(linereader.contains("src."+too)) content += linereader+"\n";
				else{
					if(linereader.contains(too)) content += "src."+linereader+"\n";
					else content += linereader+"\n";
				}
			}
			this.file_name = file_name.replaceAll(".java", "");
			dir = new String(dir_name);
			replaceSubstrings();
		}
	}
	
	public void replaceSubstrings (){
		if(content.contains(emb) && !content.contains("src."+emb)) content.replaceAll(emb,"src."+emb);
		if(content.contains(pro) && !content.contains("src."+pro)) content.replaceAll(pro, "src."+pro);
		if(content.contains(too) && !content.contains("src."+too)) content.replaceAll(too, "src."+too);
	}
	
	public void writeShallowCopy() throws IOException{
		FileOutputStream strm = new FileOutputStream(dir+file_name+"2.java");
		PrintStream print = new PrintStream(strm);
		print.println(content);
		print.close();
		strm.close();
		System.out.println("Content written to shallow copy!");
	}
	
	public static void main (String[] args) throws IOException{
		String pa = "/project/StruPPi/gabriel/workspace_old/aglappe/src/embed/contactmaps/";
		String fl = "Metric2.java";
		JavaFileCorrector jfc = new JavaFileCorrector(pa,fl);
		jfc.writeShallowCopy();
	}

}
