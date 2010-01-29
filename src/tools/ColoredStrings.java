package tools;

/**
 * Utility class for colored text output. This class contains static
 * methods to format strings such that they will be output in different
 * colors. This is done by adding ANSI escape codes to the string which
 * may or may not be supported by the terminal. The behaviour is
 * therefore platform and environment dependent. 
 * 
 * Note:
 * - The returned strings contain control characters which may cause
 *   problems if written to files or other output devices. It is
 *   therefore recommended to use colored strings only for writing
 *   to stdout.
 * - Any attributes that were set in the console before will be
 *   overwritten when printing Strings formatted by this class.
 *   
 * @author Henning Stehr, stehr@molgen.mpg.de
 * @version 0001 2006/Mar/29
 */
/**
 * @author stehr
 *
 */
public class ColoredStrings {
	
	/**
	 * An Ansi control code to specify the text foreground color.
	 * @author stehr
	 * @version 0001 2006/Mar/29
	 */
	public enum AnsiForeground {
		BLACK	(30),
		RED		(31),
		GREEN	(32),
		YELLOW	(33),
		BLUE	(34),
		MAGENTA	(35),
		CYAN	(36),
		WHITE	(37);
		int ansiCode;
		
		AnsiForeground(int ansiCode) {
			this.ansiCode = ansiCode;
		}
	}

	/**
	 * An Ansi control code to specify the text background color.
	 * @author stehr
	 * @version 0001 2006/Mar/29
	 */
	public enum AnsiBackground {
		BLACK	(40),
		RED		(41),
		GREEN	(42),
		YELLOW	(43),
		BLUE	(44),
		MAGENTA	(45),
		CYAN	(46),
		WHITE	(47);
		int ansiCode;
		
		AnsiBackground(int ansiCode) {
			this.ansiCode = ansiCode;
		}
	}	

	/**
	 * An Ansi control code to specify the text intensity.
	 * @author stehr
	 * @version 0001 2006/Mar/29
	 */	
	public enum AnsiIntensity {
		NORMAL(22),
		BOLD(1),
		DIM(2);
		
		int ansiCode;
		
		AnsiIntensity(int mode) {
			this.ansiCode = mode;
		}
	}

	/**
	 * An Ansi control code to specify the text decoration.
	 * @author stehr
	 * @version 0001 2006/Mar/29
	 */		
	public enum AnsiDecoration {
		UNDERLINE(4),
		BLINK(5),
		REVERSE(7),
		BLANK(8),
		OVERSTRIKE(9);
		
		int ansiCode;
		
		AnsiDecoration(int mode) {
			this.ansiCode = mode;
		}
	}
	
	/**
	 * Get a unicode character containing an ansi control sequence. An
	 * ansi control sequence can contain multiple ansi control codes.
	 * The supported codes are defined in the enums AnsiForeground,
	 * AnsiBackground, AnsiIntensity and AnsiDecoration.
	 * @param ansiCodes An array or a list of integers specifiying
	 * ansi control codes.
	 * @return A string containing the specified ANSI control sequence.
	 */
	private static String getAnsiControlSequence(int...ansiCodes) {
		String newStr = (char) 27 + "[";
		if(ansiCodes.length > 0) {
			newStr += ansiCodes[0];
			for(int i = 1; i < ansiCodes.length; i++) {
				newStr += ";" + ansiCodes[i];
			}
		}
		newStr += "m";		
		return newStr;
	}
	
	/**
	 * Get Ansi control sequence to reset all font attributes.
	 * Print the returned string to reset the console to normal text.
	 * @return String containing ANSI escape sequence to reset font attributes.
	 */
	public static String getConsoleResetString() {
		return getAnsiControlSequence(0);
	}
	
	/**
	 * Add ANSI format attributes to string. Supported attributes include
	 * foreground color, background color, bold, dim and some
	 * others such as blinking defined in enum AnsiDecoration.
	 * @param s String to be formatted with ANSI attributes
	 * @return String with ANSI format attributes.
	 */
	public static String getColoredString(String originalString, 
					AnsiForeground fg, AnsiBackground bg,
					AnsiIntensity ai, AnsiDecoration...ads) {
		
		String newStr = "";
		int[] decorationCodes;
		
		// add ansi codes to array
		decorationCodes = new int[3 + ads.length];
		decorationCodes[0] = fg.ansiCode;
		decorationCodes[1] = bg.ansiCode;
		decorationCodes[2] = ai.ansiCode;
		
		// convert decoration objects to integers
		for(int i = 3; i < 3 + ads.length; i++) {
			decorationCodes[i] = ads[i-3].ansiCode;
		}
		
		// add ansi control sequences to string
		newStr += getAnsiControlSequence(decorationCodes);
		newStr += originalString;
		newStr += getConsoleResetString();
		
		return newStr;
	}
	
	/**
	 * Convenience method to add a foreground color attribute to a string.
	 * @param s The string to be formatted.
	 * @return A string with an ANSI foreground color attribute.
	 */
	public static String black(String s) {
		int ansiCode = AnsiForeground.BLACK.ansiCode;
		return getAnsiControlSequence(ansiCode) + s + getConsoleResetString();
	}
	
	/**
	 * Convenience method to add a foreground color attribute to a string.
	 * @param s The string to be formatted.
	 * @return A string with an ANSI foreground color attribute.
	 */
	public static String red(String s) {
		int ansiCode = AnsiForeground.RED.ansiCode;
		return getAnsiControlSequence(ansiCode) + s + getConsoleResetString();
	}
	
	/**
	 * Convenience method to add a foreground color attribute to a string.
	 * @param s The string to be formatted.
	 * @return A string with an ANSI foreground color attribute.
	 */
	public static String green(String s) {
		int ansiCode = AnsiForeground.GREEN.ansiCode;
		return getAnsiControlSequence(ansiCode) + s + getConsoleResetString();
	}
	
	/**
	 * Convenience method to add a foreground color attribute to a string.
	 * @param s The string to be formatted.
	 * @return A string with an ANSI foreground color attribute.
	 */
	public static String yellow(String s) {
		int ansiCode = AnsiForeground.YELLOW.ansiCode;
		return getAnsiControlSequence(ansiCode) + s + getConsoleResetString();
	}
	
	/**
	 * Convenience method to add a foreground color attribute to a string.
	 * @param s The string to be formatted.
	 * @return A string with an ANSI foreground color attribute.
	 */
	public static String blue(String s) {
		int ansiCode = AnsiForeground.BLUE.ansiCode;
		return getAnsiControlSequence(ansiCode) + s + getConsoleResetString();
	}
	
	/**
	 * Convenience method to add a foreground color attribute to a string.
	 * @param s The string to be formatted.
	 * @return A string with an ANSI foreground color attribute.
	 */
	public static String magenta(String s) {
		int ansiCode = AnsiForeground.MAGENTA.ansiCode;
		return getAnsiControlSequence(ansiCode) + s + getConsoleResetString();
	}
	
	/**
	 * Convenience method to add a foreground color attribute to a string.
	 * @param s The string to be formatted.
	 * @return A string with an ANSI foreground color attribute.
	 */
	public static String cyan(String s) {
		int ansiCode = AnsiForeground.CYAN.ansiCode;
		return getAnsiControlSequence(ansiCode) + s + getConsoleResetString();
	}
	
	/**
	 * Convenience method to add a foreground color attribute to a string.
	 * @param s The string to be formatted.
	 * @return A string with an ANSI foreground color attribute.
	 */
	public static String white(String s) {
		int ansiCode = AnsiForeground.WHITE.ansiCode;
		return getAnsiControlSequence(ansiCode) + s + getConsoleResetString();
	}
	
}
