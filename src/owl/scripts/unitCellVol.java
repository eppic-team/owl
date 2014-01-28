package owl.scripts;


import owl.core.structure.CrystalCell;

public class unitCellVol {

	public static void main(String[] args) {
		if (args.length<6) {
			System.err.println("Usage: unitCellVol <a> <b> <c> <alpha> <beta> <gamma>");
			System.exit(1);
		}
		double a = Double.parseDouble(args[0]);
		double b = Double.parseDouble(args[1]);
		double c = Double.parseDouble(args[2]);
		double alpha = Double.parseDouble(args[3]);
		double beta = Double.parseDouble(args[4]);
		double gamma = Double.parseDouble(args[5]);
		CrystalCell cell = new CrystalCell(a, b, c, alpha, beta, gamma);
		System.out.printf("%8.2f\n",cell.getVolume());
	}
}
