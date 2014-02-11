package owl.core.util;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.OutputStream;
import java.nio.channels.FileChannel;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.*;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

/**
 * Collection of small useful helper methods.
 * @author stehr
 *
 */
public class Goodies {
	
	public static final boolean ASCENDING = true;
	public static final boolean DESCENDING = false;
	
	public static final String MD5_ALGORITHM = "MD5";

	
	/**
	 * Sorts a map by the values and returns a sorted map with the right ordering.
	 * @param map the initial map
	 * @param ascending if true, map will be sorted in ascending order of values, otherwise in descending order
	 * @return the sorted map
	 */
	public static <K,V extends Comparable<V>> LinkedHashMap<K,V> sortMapByValue(Map<K,V> map, boolean ascending) {
		List<Map.Entry<K,V>> list = new LinkedList<Map.Entry<K,V>>(map.entrySet());
		if(ascending) {
		Collections.sort(list, new Comparator<Map.Entry<K,V>>() {
			public int compare(Map.Entry<K,V> o1, Map.Entry<K,V> o2) {
				return (o1.getValue()).compareTo(o2.getValue());				
			}
		});
		} else {
			Collections.sort(list, new Comparator<Map.Entry<K,V>>() {
				public int compare(Map.Entry<K,V> o1, Map.Entry<K,V> o2) {
					return -(o1.getValue()).compareTo(o2.getValue());				
				}
			});			
		}

		LinkedHashMap<K,V> result = new LinkedHashMap<K,V>();
		for (Iterator<Map.Entry<K,V>> it = list.iterator(); it.hasNext();) {
			Map.Entry<K,V> entry = it.next();
			result.put(entry.getKey(), entry.getValue());
		}
		return result;
	}
	
	/**
	 * Copies a file, does not copy the metadata (file permissions etc)
	 * Java does not have a method for this in the File (or any other) class. There is
	 * only a File.renameTo() which is very flaky (does not work across file systems and
	 * so on). Java 7 (finally!) will probably come with a copyTo method in the 
	 * java.nio.file.Path class.
	 * See http://today.java.net/pub/a/today/2008/07/03/jsr-203-new-file-apis.html
	 * @param srcFile
	 * @param destFile
	 * @throws IOException
	 */
	public static void copyFile(File srcFile, File destFile) throws IOException {
		InputStream in = new FileInputStream(srcFile);
		OutputStream out = new FileOutputStream(destFile);

		byte[] buf = new byte[1024];
		int len;
		while ((len = in.read(buf)) > 0){
			out.write(buf, 0, len);
		}
		in.close();
		out.close();

	}
	
	/**
	 * Another implementation of basic file copy. Supposedly faster by using channels (I havn't compared though).
	 * Copies source file to dest file. Does not handle permissions etc.
	 * Found here: http://stackoverflow.com/questions/106770/standard-concise-way-to-copy-a-file-in-java
	 * @param srcFile
	 * @param destFile
	 * @throws IOException
	 */
	public static void copyFileFast(File sourceFile, File destFile) throws IOException {
		 if(!destFile.exists()) {
		  destFile.createNewFile();
		 }

		 FileChannel source = null;
		 FileChannel destination = null;
		 try {
		  source = new FileInputStream(sourceFile).getChannel();
		  destination = new FileOutputStream(destFile).getChannel();
		  destination.transferFrom(source, 0, source.size());
		 }
		 finally {
		  if(source != null) {
		   source.close();
		  }
		  if(destination != null) {
		   destination.close();
		  }
		}
	}
		
	/**
	 * Compute the MD5 sum of a given input String using the java.security library.
	 * @param input
	 * @return
	 */
	public static final String computeMD5Sum(String input) {

		if (input == null) {
			throw new IllegalArgumentException("Input cannot be null!");
		}

		StringBuffer sbuf = new StringBuffer();
		MessageDigest md = null;
		try {
			md = MessageDigest.getInstance(MD5_ALGORITHM);
		} catch (NoSuchAlgorithmException e) {
			System.err.println("Unexpected error while computing md5 hash");
			System.err.println(e.getMessage());
			System.exit(1);
		}
		byte [] raw = md.digest(input.getBytes());

		for (int i = 0; i < raw.length; i++) {
			int c = (int) raw[i];
			if (c < 0) {
				c = (Math.abs(c) - 1) ^ 255;
			}
			String block = toHex(c >>> 4) + toHex(c & 15);
			sbuf.append(block);
		}

		return sbuf.toString();

	}

	private static final String toHex(int s) {
		if (s < 10) {
			return new StringBuffer().
			append((char)('0' + s)).
			toString();
		} else {
			return new StringBuffer().
			append((char)('A' + (s - 10))).
			toString();
		}
	}

	/**
	 * To serialize to file a given Serializable object
	 * @param serializedFile
	 * @param obj
	 * @throws IOException
	 */
	public static void serialize(File serializedFile, Object obj) throws IOException {
		FileOutputStream fileOut = new FileOutputStream(serializedFile);
		ObjectOutputStream out = new ObjectOutputStream(fileOut);
		out.writeObject(obj);
		out.close();
		fileOut.close();
	}

	/**
	 * To deserialize from file a Serializable object. The returned object must be cast to the appropriate class.
	 * @param serialized
	 * @return
	 * @throws IOException
	 * @throws ClassNotFoundException
	 */
	public static Object readFromFile(File serialized) throws IOException, ClassNotFoundException {
		FileInputStream fileIn = new FileInputStream(serialized);
		ObjectInputStream in = new ObjectInputStream(fileIn);
		Object obj = in.readObject();
		in.close();
		fileIn.close();
		return obj;
	}
	
	/**
	 * Gunzips (decompresses gzip) given gzFile into outFile using java.util gzip implementation
	 * @param inFile
	 * @param outFile
	 * @throws IOException
	 */
	public static void gunzipFile(File inFile, File outFile) throws IOException {
		GZIPInputStream zis = new GZIPInputStream(new FileInputStream(inFile));
		FileOutputStream os = new FileOutputStream(outFile);
		int b;
		while ( (b=zis.read())!=-1) {
			os.write(b);
		}
		zis.close();
		os.close();
	}
	
	/**
	 * Gzips given inFile into outFile gzip file using java.util gzip implementation
	 * @param inFile
	 * @param outFile
	 * @throws IOException
	 */
	public static void gzipFile(File inFile, File outFile) throws IOException {
		GZIPOutputStream zos = new GZIPOutputStream(new FileOutputStream(outFile));
		FileInputStream is = new FileInputStream(inFile);

		int b;
		while ( (b=is.read())!=-1) {
			zos.write(b);
		}
		zos.close();
		is.close();
	}
}
