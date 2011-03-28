reconstruct

Protein contact map reconstruction using the TINKER package

http://www.molgen.mpg.de/~lappe/reconstruct

INSTALLATION

	You need Java 1.6 or newer (available from http://java.sun.com).

	1. Get TINKER and PRM files from http://dasher.wustl.edu/tinker/
	
	   Linux TINKER executables are available at: 
	       http://dasher.wustl.edu/tinker/downloads/linux.tar.gz
	   Force Field Parameter files are available at:
	       http://dasher.wustl.edu/tinker/distribution/params/
	       
	   Unfortunately TINKER is written in statically allocated FORTRAN so 
	   memory is allocated at program startup and there must be enough of it 
	   available or the program will fail. The allocation sizes are controlled
	   by static variables in sizes.i file (in the TINKER source distribution).
	   
	   The default binaries provided in the TINKER web site are
	   compiled with too restrictive static values. In practice this means that
	   the "distgeom" program will only run for very small proteins.
	   
	   Thus in order to be able to run the reconstruct program for reasonably 
	   sized proteins you will need  to download the TINKER source code at:
	   http://dasher.wustl.edu/tinker/downloads/tinker-5.1.02.tar.gz
	   modify the sizes.i constants file and recompile it.
	   The parameters to modify are MAXGEO and MAXATM and MAXKEY. 
	   Values MAXATM=100000, MAXGEO=10000 and MAXKEY=20000 should suffice. These 
	   values would already require a memory allocation of >1GB
	   
	   If the programs can not run because they can't allocate enough memory one
	   work-around is to increase the size of your swap file. TINKER only 
	   allocates a lot of memory at startup but most of the time doesn't actually
	   use most of it. 
	
	2. Edit the file reconstruct.cfg and set the parameters TINKER_BIN_DIR and 
	   PRM_FILE. The only type of PRM_FILE supported is AMBER force field 
	   parameter files.
	   A per-user reconstruct.cfg file can be placed in the user's home 
	   directory.
	      
	3. Run it 
	   ./reconstruct 


PARALLEL VERSION, SUN GRID ENGINE

For the parallel version to work (option -A) the path to the SGE root directory 
needs to be set in the reconstruct shell script (sgeroot variable)

CONTACT MAP FILES

The contact map files that the reconstruct program reads are simple text files with 
a few headers and 3 columns: 1st for i residue numbers, 2nd for j residue numbers and 
3rd for weights (currently ignored).
An example file is provided: sample.cm (Cbeta 8A cutoff contact map for PDB 1bxy, chain A)
The format is the same used by our CMView contact map visualization program (see 
http://www.molgen.mpg.de/~lappe/cmview)
The headers are essential for the reconstruct program to work. The parameters of the 
contact map: SEQUENCE, CONTACT TYPE (CT) and CUTOFF are read from the headers.
 