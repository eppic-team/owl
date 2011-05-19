package owl.core.structure;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import javax.vecmath.Vector3d;

public class PdbUnitCell implements Iterable<PdbAsymUnit>{

	private List<PdbAsymUnit> units;
	
	public PdbUnitCell() {
		this.units = new ArrayList<PdbAsymUnit>();
	}
	
	public void addUnit(PdbAsymUnit unit) {
		units.add(unit);
	}
	
	public List<PdbAsymUnit> getAllAsymUnits() {
		return units;
	}

	public int getNumAsymUnits() {
		return units.size();
	}
	
	/**
	 * Gets the asym unit by index.
	 * @param i
	 * @return
	 */
	public PdbAsymUnit getAsymUnit(int i) {
		return units.get(i);
	}
	
	/**
	 * Translates this PdbUnitCell to the given unit cell (direction).
	 * e.g. doCrystalTranslation(new Vector3d(1,1,1)) will translate this PdbUnitCell to 
	 * crystal cell (1,1,1), considering always this PdbUnitCell's cell to be (0,0,0)
	 * @param direction
	 */
	public void doCrystalTranslation(Vector3d direction) {
		for (PdbAsymUnit pdb:units) {
			pdb.doCrystalTranslation(direction);
		}		
		
	}
	
	public Iterator<PdbAsymUnit> iterator() {
		return units.iterator();
	}
	
	public PdbUnitCell copy() {
		PdbUnitCell newCell = new PdbUnitCell();
		for (PdbAsymUnit pdb:this) {
			PdbAsymUnit newPdb = pdb.copy();
			newCell.addUnit(newPdb);
		}
		return newCell;
	}
}
