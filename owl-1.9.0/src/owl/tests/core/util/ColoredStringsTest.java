package owl.tests.core.util;

import org.junit.Test;

import owl.core.util.ColoredStrings;


public class ColoredStringsTest {
	
	static final ColoredStrings.AnsiForeground 
		DefaultAnsiForeground = ColoredStrings.AnsiForeground.BLACK;
	static final ColoredStrings.AnsiBackground 
		DefaultAnsiBackground = ColoredStrings.AnsiBackground.WHITE;
	static final ColoredStrings.AnsiIntensity 
		DefaultAnsiIntensity = ColoredStrings.AnsiIntensity.NORMAL;
	
	/**
	 * Test method for class ColoredStrings
	 */
	@Test
	public void testColoredStrings() {
		
		System.out.println("Testing ColoredStrings");
		
		// testing convenience methods
		System.out.println(ColoredStrings.red("Text in red"));
		System.out.println(ColoredStrings.green("Text in green"));
		System.out.println(ColoredStrings.yellow("Text in yellow"));
		System.out.println(ColoredStrings.blue("Text in blue"));
		System.out.println(ColoredStrings.magenta("Text in magenta"));
		System.out.println(ColoredStrings.cyan("Text in cyan"));
		System.out.println(ColoredStrings.white("Text in white"));
		System.out.println(ColoredStrings.black("Text in black"));
		
		// testing no attributes
		System.out.println(ColoredStrings.getConsoleResetString());
		
		// testing foreground attributes
		for(ColoredStrings.AnsiForeground fg: ColoredStrings.AnsiForeground.values()) {			
			String str = fg.name() + " foreground";
			String fStr = ColoredStrings.getColoredString(str, fg, DefaultAnsiBackground, DefaultAnsiIntensity);
			System.out.println(fStr);
		}
		
		// testing background attributes
		for(ColoredStrings.AnsiBackground bg: ColoredStrings.AnsiBackground.values()) {			
			String str = bg.name() + " background";
			String fStr = ColoredStrings.getColoredString(str, DefaultAnsiForeground, bg, DefaultAnsiIntensity);
			System.out.println(fStr);
		}		
		
		// testing text intensities
		for(ColoredStrings.AnsiIntensity ai: ColoredStrings.AnsiIntensity.values()) {			
			String str = ai.name() + " intensity";
			String fStr = ColoredStrings.getColoredString(str, DefaultAnsiForeground, DefaultAnsiBackground, ai);
			System.out.println(fStr);
		}
		
		// testing decorations
		for(ColoredStrings.AnsiDecoration ad: ColoredStrings.AnsiDecoration.values()) {			
			String str = ad.name() + " decoration";
			String fStr = ColoredStrings.getColoredString(str, DefaultAnsiForeground, DefaultAnsiBackground, DefaultAnsiIntensity, ad);
			System.out.println(fStr);
		}
	
		// testing multiple decorations
		String str = "Multiple decorations";
		String fStr = ColoredStrings.getColoredString(str, 
				DefaultAnsiForeground, 
				DefaultAnsiBackground, 
				DefaultAnsiIntensity, 
				ColoredStrings.AnsiDecoration.BLINK,
				ColoredStrings.AnsiDecoration.OVERSTRIKE,
				ColoredStrings.AnsiDecoration.REVERSE,
				ColoredStrings.AnsiDecoration.UNDERLINE);
		System.out.println(fStr);
		
		// done
		System.out.println("Done testing ColoredStrings");
		
	}
	
	public static void main(String args[]) {
		ColoredStringsTest myself = new ColoredStringsTest();
		myself.testColoredStrings();
	}

}
