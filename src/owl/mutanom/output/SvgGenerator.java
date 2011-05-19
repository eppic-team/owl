package owl.mutanom.output;

import java.awt.BasicStroke;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.GradientPaint;
import java.awt.Rectangle;
import java.awt.Graphics2D;
import java.awt.Color;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Line2D;
import java.awt.geom.RoundRectangle2D;
import java.io.File;
import java.io.FileOutputStream;
import java.io.Writer;
import java.io.OutputStreamWriter;
import java.io.IOException;



import org.apache.batik.svggen.SVGGraphics2D;
import org.apache.batik.dom.GenericDOMImplementation;

import org.w3c.dom.Document;
import org.w3c.dom.DOMImplementation;

import owl.core.features.Feature;
import owl.core.features.FeatureType;
import owl.core.util.Interval;
import owl.mutanom.core.Gene;
import owl.mutanom.core.Mutation;
import owl.mutanom.core.Substructure;
import owl.mutanom.core.TargetList;
import owl.mutanom.core.Substructure.SubstructureType;


public class SvgGenerator {
	
	/**
	 * Only for testing
	 * @param g2d
	 */
    public void paint(Graphics2D g2d) {
    	
    	// global parameters
    	int width = 1000;	
    	int minorTickLen = 3;
    	int majorTickLen = 5;
    	//int capLen = 7;
    	int titleXOffset = -5;
    	int titleYOffset = 30;
    	int pdbYOffset = 0;
    	int pdbHeight = 10;
    	int mutLength = 15;
    	int mutYOffset = 12;
    	int baselineX = 10;	// in screen coordinates
    	
    	// input data
    	int geneLength = 766;	// sequence length
    	int baselineY = 50;	// in screen coordinates
    	int pdbStart = 350;	// in sequence coordinates
    	int pdbEnd = 680;	// in sequence coordinates
    	double factor = 1.0 * width / geneLength;
    	
    	// graphics settings
        final BasicStroke myStroke = new BasicStroke(1.0f, 
                                              BasicStroke.CAP_BUTT, 
                                              BasicStroke.JOIN_MITER
                                              );
        g2d.setStroke(myStroke);

    	// raw ruler baseline
    	g2d.setPaint(Color.black);
//    	g2d.drawLine(baselineX, baselineY, baselineX + width, baselineY);
//    	g2d.drawLine(baselineX, baselineY-capLen, baselineX, baselineY+capLen);
//    	g2d.drawLine(baselineX+width, baselineY-capLen, baselineX+width, baselineY+capLen);
    	
    	// draw tickmarks
    	for(int i = 1; i <= geneLength; i++) {
    		if(i % 100 == 0) {
    			int x = baselineX + (int) Math.round(i*factor);
    			//g2d.drawLine(x, baselineY - majorTickLen, x, baselineY + majorTickLen);
    			g2d.drawLine(x, baselineY + 2, x, baselineY + 2 + majorTickLen * 2);
    		} else
    		if(i % 10 == 0) {
    			int x = baselineX + (int) Math.round(i*factor);
    			//g2d.drawLine(x, baselineY - minorTickLen, x, baselineY + minorTickLen);
    			g2d.drawLine(x, baselineY + 2, x, baselineY + 2 + minorTickLen * 2);
    		}
    	}
    	
    	// draw structure region 1
    	int x1 = baselineX + (int) Math.round(factor * pdbStart);
    	int x2 = baselineX + (int) Math.round(factor * pdbEnd);
    	g2d.setPaint(new GradientPaint(0.0f, baselineY - 0.7f * pdbHeight, new Color(0.4f,0.8f,1.0f), 0.0f, baselineY, new Color(0.0f,0.0f,0.8f), true));
    	g2d.fill(new Rectangle(x1, baselineY - pdbHeight, x2-x1, pdbHeight));
    	int ellWidth = 4;
    	g2d.fill(new Ellipse2D.Double(x1 - 0.5 * ellWidth, baselineY-pdbHeight, ellWidth, pdbHeight));
    	g2d.setPaint(new Color(0.2f,0.4f,1.0f));
    	g2d.fill(new Ellipse2D.Double(x2 - 0.5 * ellWidth, baselineY-pdbHeight, ellWidth, pdbHeight));
    	
    	// draw structure region 2
    	pdbStart = 100;
    	pdbEnd = 250;
    	x1 = baselineX + (int) Math.round(factor * pdbStart);
    	x2 = baselineX + (int) Math.round(factor * pdbEnd);
    	g2d.setPaint(new GradientPaint(0.0f, baselineY - 0.7f * pdbHeight, new Color(0.4f,0.8f,1.0f), 0.0f, baselineY, new Color(0.0f,0.0f,0.8f), true));
		g2d.fill(new RoundRectangle2D.Double(x1, baselineY - pdbHeight + pdbYOffset, x2-x1, pdbHeight + pdbYOffset,5,5));
//    	g2d.fill(new Rectangle(x1, baselineY - pdbHeight, x2-x1, pdbHeight));
//    	ellWidth = 4;
//    	g2d.fill(new Ellipse2D.Double(x1 - 0.5 * ellWidth, baselineY-pdbHeight, ellWidth, pdbHeight));
//    	g2d.setPaint(new Color(0.2f,0.4f,1.0f));
//    	g2d.fill(new Ellipse2D.Double(x2 - 0.5 * ellWidth, baselineY-pdbHeight, ellWidth, pdbHeight));
    	
    	// draw mutations
    	g2d.setPaint(Color.red);
    	int[] mut = {50,51,52,53,55,75,80,150,560,563,565};
    	for(int i = 0; i < mut.length; i++) {
    		int x = baselineX + (int) Math.round(factor * mut[i]);
    		g2d.drawLine(x,baselineY - mutYOffset - mutLength, x, baselineY - mutYOffset);
    	}
    	
    	// draw gene name
    	g2d.setPaint(Color.black);
    	g2d.drawString("BRAF", baselineX + titleXOffset, baselineY - titleYOffset);
    	
//        g2d.setPaint(Color.red);
//        g2d.fill(new Rectangle(10, 10, 100, 100));
    }

    public void paintSeqOverview(Graphics2D g2d, TargetList tl) {
    	
    	// global parameters
    	int width = 800;	// width of plot on screen (not including baselineX) - 800 for print, 1200 for presentation
    	int baselineX = 10;	// offset from left side of screen in screen coordinates
		int baselineY = -10;	// in screen coordinates
    	int baselineYStep = 95;	// in screen coordinates
    	// ticks and caps
    	int minorTickLen = 2;	// half length of minor tick marks
    	int majorTickLen = 4;	// half length of major tick marks
    	int capLen = 6;			// half length of terminal caps
    	int ticklabelYOffset = 4; // space between label and baselineY
    	int ticklabelXOffset = 0; // additional offset of tick labels from tick mark
    	// gene names
    	int titleXOffset = -5;		// label offset from baselineX
    	int titleYOffset = 43;		// label offset from baselineY
    	// structured regions
    	int pdbYOffset = -7;			// offset from baselineY
    	int pdbHeight = 11;			// height of pdb rectangle
    	int rectRadius = pdbHeight/2; // edge-radius for rounded rectangles
    	boolean useGradient = false; // whether to paint rounded rectangles with a gradient
    	// mutations
    	int mutLength = 15;			// length of mutation marker line
    	int mutYOffset = 20;		// offset from baselineY
    	// domains
    	boolean showDomains = false;	// query and display domains from pDomains
    	int domYOffset = 25;		// offset from baselineY
    	int domHeight = 5;			// height of box depicting domain

    	// colors & fonts
    	Color rulerColor = Color.black;
    	//Color pdbColor = new Color(0.2f, 0.4f, 1.0f);
    	Color mutColor = new Color(0.8f, 0.0f, 0.0f);
    	Color geneNameColor = Color.black;
    	Color pdbColor = new Color(0.25f, 0.45f, 0.9f);	// used if useGradient = false
    	Color pdbGradStart = new Color(0.4f,0.6f,1.0f);
    	Color pdbGradEnd = new Color(0.1f,0.3f,0.8f);
    	Color pdbLabelColor = Color.white;
    	Color domainColor = Color.green;
    	
    	Font geneNamefont = new Font("Arial", Font.PLAIN, 14);
    	Font tickLabelFont = new Font("Arial", Font.PLAIN, 10);
    	Font pdbLabelFont = new Font("Arial", Font.BOLD, 10);
    	
    	// stroke settings 
    	// for screen
    	// final BasicStroke mutStroke = new BasicStroke(1.0f);	// for screen
    	// final BasicStroke tickStroke = new BasicStroke(1.0f);	// for screen
    	
    	// for print (page width)
        // final BasicStroke mutStroke = new BasicStroke(0.5f);   // for print
        // final BasicStroke tickStroke = new BasicStroke(0.3f);  // for print
        
        // for print (single column)
        final BasicStroke mutStroke = new BasicStroke(0.3f);
        final BasicStroke tickStroke = new BasicStroke(0.3f);        
        
        // draw genes
    	for(Gene g:tl.getTargets()) {
    	
    		// input data
    		int geneLength = g.getLength();	// sequence length
    		baselineY += baselineYStep;
    		double factor = 1.0 * width / geneLength;

    		// draw gene name
	    	g2d.setFont(geneNamefont);
    		g2d.setPaint(geneNameColor);
    		g2d.drawString(g.getGeneName(), baselineX + titleXOffset, baselineY - titleYOffset);
    		
    		// raw ruler baseline
            g2d.setStroke(tickStroke);
    		g2d.setPaint(rulerColor);
    		g2d.draw(new Line2D.Double(baselineX, baselineY, baselineX + width, baselineY));
    		// left cap
    		g2d.draw(new Line2D.Double(baselineX, baselineY, baselineX+5, baselineY-capLen));
    		g2d.draw(new Line2D.Double(baselineX, baselineY, baselineX+5, baselineY+capLen));
    		// right cap
    		g2d.draw(new Line2D.Double(baselineX+width, baselineY, baselineX+width-5, baselineY-capLen));
    		g2d.draw(new Line2D.Double(baselineX+width, baselineY, baselineX+width-5, baselineY+capLen));
    		
    		// draw tickmarks
    		int majorStep;	// this is where numeric labels are written
    		int minorStep;	// this is where minor tick marks are written
    		boolean abbK;	// if true, abbreviate 1000 as 1k etc.
    		if(geneLength > 1000) {
    			majorStep = 100;
    			minorStep = 10;
    			abbK = true;
    		} else
    		if(geneLength > 200) {
    			majorStep = 100;
    			minorStep = 10;
    			abbK = false;
    		} else {
    			majorStep = 10;
    			minorStep = 1;
    			abbK = false;
    		}
    		
    		// set font and font metric for tick labels
        	g2d.setFont(tickLabelFont);
    	    FontMetrics tickLabelmetrics = g2d.getFontMetrics(tickLabelFont);
    		
    		for(int i = 2; i < geneLength; i++) {
    			if(i % majorStep == 0) {
    				double x = baselineX + (i-1)*factor;
    				g2d.draw(new Line2D.Double(x, baselineY - majorTickLen, x, baselineY + majorTickLen));
    				String tickLabel = abbK?String.format("%3.1fk",0.001*i):String.format("%d",i);
    	    	    int hgt = tickLabelmetrics.getHeight();
    	    	    int adv = tickLabelmetrics.stringWidth(tickLabel);
    				g2d.drawString(tickLabel, (float) x - adv / 2 + ticklabelXOffset, (float) 1.0 * baselineY + hgt + ticklabelYOffset);
    			} else
    				if(i % minorStep == 0) {
    					double x = baselineX + (i-1)*factor;
    					g2d.draw(new Line2D.Double(x, baselineY - minorTickLen, x, baselineY + minorTickLen));
    				}
    		}

    		// draw structure regions (see also Gene.writeSubstructureData)
			g2d.setFont(pdbLabelFont);
			FontMetrics pdbLabelmetrics = g2d.getFontMetrics(pdbLabelFont);
    		SubstructureType restrictToType = SubstructureType.XRAY;
			for (Substructure ss : g.getSubstructures()) {
				if(restrictToType==null || ss.getType() == restrictToType) {
					int pdbStart = ss.getBegPos();
					int pdbEnd = ss.getEndPos();
					String pdbLabel = (ss.getType()==SubstructureType.PREDICTION?ss.getTemplatesStr():ss.getPdbCode().toUpperCase()+" "+ss.getChainCode());
					double x1 = baselineX + factor * pdbStart;
					double x2 = baselineX + factor * pdbEnd-1;
					if(useGradient) g2d.setPaint(new GradientPaint(0.0f, baselineY - 0.7f * pdbHeight, pdbGradStart, 0.0f, baselineY, pdbGradEnd, true));
					else g2d.setPaint(pdbColor);
					g2d.fill(new RoundRectangle2D.Double(x1, baselineY - pdbHeight + pdbYOffset, x2-x1, pdbHeight,rectRadius,rectRadius));
					//g2d.fill(new Rectangle.Double(x1, baselineY - pdbHeight + pdbYOffset, x2-x1, pdbHeight + pdbYOffset));
					// draw pdb label
					g2d.setPaint(pdbLabelColor);
    	    	    int adv = pdbLabelmetrics.stringWidth(pdbLabel);
					g2d.drawString(pdbLabel, (int) Math.round( 0.5 * (x1 + x2)) - adv / 2, baselineY + pdbYOffset - 2);
					
//					// draw domains within structural regions
//					if(showDomains && ss.type != SubstructureType.PREDICTION) {
//						PDomainsConnection pDomains = new PDomainsConnection();
//						Map<String, IntervalSet> domainMap = pDomains.getDomains(ss.getPdbCode(), ss.getChainCode(), domainType);
//						if(domainMap.size() > 0) {
//							for(String domName:domainMap.keySet()) {
//								IntervalSet intervals = domainMap.get(domName);
//								//System.out.print(domainType.toString() + " domain " + domName + ": " + intervals + " ");
//								for(Interval domInt:intervals) {
//									int domPdbStart = domInt.beg;
//									int domPdbEnd = domInt.end;
//									int domStart = Substructure.ALIGNMENT_UNDEFINED;
//									while(domStart == Substructure.ALIGNMENT_UNDEFINED && domPdbStart < domInt.end) {		// if starting residue is unobserved, try next
//										domStart = ss.mapPdbResser2Uniprot(domPdbStart++);
//									}
//									int domEnd = Substructure.ALIGNMENT_UNDEFINED;
//									while(domEnd == Substructure.ALIGNMENT_UNDEFINED && domPdbEnd > domInt.beg) {		// // if final residue is unobserved, try previous
//										domEnd = ss.mapPdbResser2Uniprot(domPdbEnd--);
//									}
//									if(domStart == Substructure.ALIGNMENT_UNDEFINED || domEnd == Substructure.ALIGNMENT_UNDEFINED) {
//										System.err.println("Warning: Domain " + domInt.beg + "-" + domInt.end + " could not be mapped to Uniprot sequence.");
//									}
//									//System.out.print(domStart + "-" + domEnd + " ");
//									double dx1 = baselineX + factor * domStart;
//									double dx2 = baselineX + factor * domEnd;
//									g2d.setPaint(domainColor);
//									g2d.fill(new RoundRectangle2D.Double(dx1, baselineY - domHeight + domYOffset, dx2-dx1, domHeight,rectRadius,rectRadius));
//								}
//								//System.out.println();
//							}
//						} else {
//							System.out.println("No domains found for " + ss.getPdbCode() + ss.getChainCode());
//						}
//					}
				}
			}
			
			// draw domains
			if(showDomains) {
				g2d.setPaint(domainColor);
				for(Feature f:g.getFeatures()) {
					if(f.getType() == FeatureType.SDOMAIN) {
						//String domLabel = f.getDescription();
						System.out.print(f.toString() + " ");
						for(Interval intv:f.getIntervalSet()) {
							int domStart = intv.beg;
							int domEnd = intv.end;
							System.out.print(domStart + "-" + domEnd + " ");
							double dx1 = baselineX + factor * domStart;
							double dx2 = baselineX + factor * domEnd;
							g2d.fill(new RoundRectangle2D.Double(dx1, baselineY - domHeight + domYOffset, dx2-dx1, domHeight,rectRadius,rectRadius));						
						}
						System.out.println();
					}
				}
			}

    		// draw mutations (see also Gene.writeMutationData)
            g2d.setStroke(mutStroke);
    		g2d.setPaint(mutColor);
			for (Mutation m : g.getMutations()) {
				Substructure ss = g.getSubstructure(m.getPos());
				if(restrictToType==null || (ss != null && ss.getType() == restrictToType) && ss.getPdb().containsStdAaResidue(ss.mapUniprotResser2Cif(m.getPos()))) {
					int mutPos = m.getPos();
					//String mutLabel = m.getMutStr();
	    			double x = baselineX + factor * mutPos;
	    			g2d.draw(new Line2D.Double(x,baselineY - mutYOffset - mutLength, x, baselineY - mutYOffset));					
				}
			}		
    	}
    }    
    
    
    public static void writeSeqOverviewForPrint(File outFile, TargetList targets) {
        // Get a DOMImplementation.
        DOMImplementation domImpl =
            GenericDOMImplementation.getDOMImplementation();

        // Create an instance of org.w3c.dom.Document.
        String svgNS = "http://www.w3.org/2000/svg";
        Document document = domImpl.createDocument(svgNS, "svg", null);

        // Create an instance of the SVG Generator.
        SVGGraphics2D svgGenerator = new SVGGraphics2D(document);

        // Ask the test to render into the SVG Graphics2D implementation.
        SvgGenerator test = new SvgGenerator();
        test.paintSeqOverview(svgGenerator, targets);

        // Finally, stream out SVG to the standard output using
        // UTF-8 encoding.

        try {
        	FileOutputStream fout = new FileOutputStream(outFile);
        	boolean useCSS = true; // we want to use CSS style attributes
        	Writer out = new OutputStreamWriter(fout, "UTF-8");
        	svgGenerator.stream(out, useCSS);    	
        } catch (IOException e) {
        	System.err.println("Error writing to SVG file " + outFile);
        }
    }
    
    public static void main(String[] args) throws IOException {

        // Get a DOMImplementation.
        DOMImplementation domImpl =
            GenericDOMImplementation.getDOMImplementation();

        // Create an instance of org.w3c.dom.Document.
        String svgNS = "http://www.w3.org/2000/svg";
        Document document = domImpl.createDocument(svgNS, "svg", null);

        // Create an instance of the SVG Generator.
        SVGGraphics2D svgGenerator = new SVGGraphics2D(document);

        // Ask the test to render into the SVG Graphics2D implementation.
        SvgGenerator test = new SvgGenerator();
        test.paint(svgGenerator);

        // Finally, stream out SVG to the standard output using
        // UTF-8 encoding.
        
        FileOutputStream fout = new FileOutputStream("/project/StruPPi/henning/projects/mutanom/analysis/results/saved_results/test/test.svg");
        
        boolean useCSS = true; // we want to use CSS style attributes
        Writer out = new OutputStreamWriter(fout, "UTF-8");
        svgGenerator.stream(out, useCSS);
    }
}