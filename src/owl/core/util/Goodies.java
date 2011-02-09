package owl.core.util;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.*;

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

}