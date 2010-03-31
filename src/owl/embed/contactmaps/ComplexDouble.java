package owl.embed.contactmaps;

/**
 * This class can compute basic algebraic and analytical entities of elements of
 * the field of all complex numbers <tt>z = x + i y in C</tt>, where
 * <tt>i^2 = -1</tt>. It provides methods to convert the <tt>(x, y)</tt> representation
 * to angular representation <tt>(r,phi)</tt>, where <tt>r >= 0 </tt> is the absolute
 * value and <tt>phi \in [-pi/2,pi/2]</tt> is the angle in the complex number plane.
 *
 * @author gmueller
 */
public class ComplexDouble {

    //private static final Double NaN = Double.NaN;

    private Double re_z;		//the real part, i.e. Re(z)

    private Double im_z;		//the imaginary part, i.e. Im(z)

    private Double arg;			//the angular argument, i.e. arctan (Im(z)/Re(z)) in [-pi/2,pi/2]

    private Double abs;			//the absolute value, i.e. Re(z)^2+Im(z)^2

    /**
     * Zero parameter constructor: this constructor should only be used, if right after calling this constructor
     * one of the setters <code>{@link #setCoeffs(ComplexDouble)}</code>,<code>{@link #setCoeffs(Double, Double)}</code> or 
     * <code>{@link #setCoeffs(double, double)}</code> is invoked. Example code:
     * <p>
     * <p><tt>double x = 0.0, y = 0.0</tt></p>
     * <p><tt>ComplexDouble c = new ComplexDouble();</tt></p>
     * <p><tt>c.setCoeffs(x,y);</tt></p>
     * </p>
     */
    public ComplexDouble (){
    	re_z = new Double(0.0);
    	im_z = new Double(0.0);
    	arg  = new Double(0.0);
    	abs  = new Double(0.0);
    };
    /**
     * One parameter constructor: initializes an instance of this class, by copying
     * the original <code>{@link ComplexDouble}</code> instance.
     * @param z the original complex number
     */
    public ComplexDouble (ComplexDouble z){
        this.setCoeffs(z);
    }

    /**
     * Two parameter constructor: initializes an instance of this class, which
     * represents a complex number <tt>z = x + i y</tt>.
     * @param x the real part of <tt>z</tt>
     * @param y the imaginary part of <tt>z</tt>
     */
    public ComplexDouble (double x, double y){
        this.setCoeffs(x,y);
    }
    
    /**
     * Two parameter constructor: initializes an instance of this class, by using
     * the argument and absolute value as parameters.
     * @param ang the argument (angle in the complex number plane)
     * @param abs_val the absolute value
     * @throws IllegalArgumentException if the absolute value is less than zero
     */
    public ComplexDouble (Double ang, Double abs_val) throws IllegalArgumentException {
    	this.setCoeffs(ang,abs_val);
    }

    /**
     * Setter: initializes the fields for real part, imaginary part, angle and
     * absolute value.
     * @param x the real part
     * @param y the imaginary part
     */
    public void setCoeffs (double x, double y){
        this.re_z = new Double (x);
        this.im_z = new Double (y);
        if(x != 0.0){
            this.arg  = new Double (Math.atan(y/x));
        }
        else{
            if(y == 0.0){
                this.arg  = new Double (0.0);
            }
            else{
                if(y > 0) this.arg = new Double (Math.PI/2.0);
                else this.arg = new Double (-Math.PI/2.0);
            }
        }
        this.abs  = new Double (Math.sqrt(x*x + y*y));
    }

    /**
     * Setter: copies the parameter <tt>z</tt>.
     * @param z a complex (double) number
     */
    public void setCoeffs (ComplexDouble z){
        this.re_z = new Double (z.getRealPart());
        this.im_z = new Double (z.getImPart());
        this.arg  = new Double (z.getArg());
        this.abs  = new Double (z.getAbs());
    }
    
    /**
     * Setter: initializes an instance of this class by using the
     * argument and absolute value as parameters.
     * @param ang the argument (angle in the complex number plane)
     * @param abs_val the absolute value of the complex number
     * @throws IllegalArgumentException if the absolute value is less than zero
     */
    public void setCoeffs(Double ang, Double abs_val) throws IllegalArgumentException {
    	double angle = ang.doubleValue(), absolute_val = abs_val.doubleValue();
    	if(abs_val.doubleValue() >= 0.0){
    		this.abs  = new Double(absolute_val);
    		this.arg  = new Double(angle);
    		this.re_z = new Double(absolute_val*Math.cos(angle));
    		this.im_z = new Double(absolute_val*Math.sin(angle));
    	}
    	else throw new IllegalArgumentException ("The absolute value must always be greater than or equal to zero!");
    }

    /**
     * Adds this instance to the parameter <tt>z</tt> and returns a new instance
     * of this class.
     * @param z a complex number
     * @return the sum as a new instance
     */
    public ComplexDouble add(ComplexDouble z){
        double x = this.getRealPart(), y = this.getImPart();
        return new ComplexDouble (x + z.getRealPart(), y + z.getImPart());
    }

    /**
     * Converts this instance to a <tt>Double</tt> array, where the first entry equals
     * the real part and the second equals the imaginary part.
     * @return a <tt>Double</tt> array representation of this object
     */
    public Double[] convertComplexToDoubleArray (){
        Double[] ar = {new Double (this.getRealPart()), new Double (this.getImPart())};
        return ar;
    }
    
    /**
     * Converts this instance to its conjugate: that is <tt>phi : C -> C</tt>
     * the non-trivial automorphism of the field of complex numbers, mapping each complex <tt>z</tt>
     * onto its conjugate <tt>z*</tt>, s.t. <tt>z z* = |z|^2</tt>, <tt>z + z* = 2 Re(z)</tt> and
     * <tt>z - z* = 2 Im(z)</tt>.
     * @return the conjugate
     */
    public ComplexDouble conjugate (){
        return new ComplexDouble (getRealPart(),-getImPart());
    }

    /**
     * Getter: converts a <tt>Double</tt> to a <tt>ComplexDouble</tt> by taking
     * the parameter <tt>real_part</tt> as the real part of the complex number and
     * sets the imaginary part to zero.
     * @param real_part the real part
     * @return a <tt>ComplexDouble</tt> with zero as imaginary part
     */
    public ComplexDouble convertRealDoubleToComplex (Double real_part){
        return new ComplexDouble (real_part.doubleValue(), 0.0);
    }

    /**
     * Getter: converts a <tt>Double</tt> to a <tt>ComplexDouble</tt> by taking
     * the parameter <tt>im_part</tt> as the imaginary part of the complex number and
     * sets the real part to zero.
     * @param im_part the imaginary part
     * @return a <tt>ComplexDouble</tt> with zero as real part
     */
    public ComplexDouble convertImaginaryDoubleToComplex (Double im_part){
        return new ComplexDouble (0.0, im_part.doubleValue());
    }

    /**
     * Multiplies this instance to the parameter <tt>z</tt> and returns a new
     * instance of the class.
     * @param z a complex number
     * @return the product as a new instance
     */
    public ComplexDouble multiply (ComplexDouble z){
        double x = getRealPart(), y = getImPart();
        return new ComplexDouble (x*z.getRealPart() - y * z.getImPart(), x * z.getImPart() + y * z.getRealPart());
    }
    
    /**
     * Multiplies a <tt>double</tt> to ComplexDouble and returns the product as a new instance.
     * @param scalar a <tt>double</tt> scalar
     * @return the scalar multiple of this instance
     */
    public ComplexDouble multiply (double scalar){
    	return new ComplexDouble(getRealPart()*scalar,getImPart()*scalar);
    }

    /**
     * Returns the <tt>k</tt> power of this instance as new instance of this class.
     * @param k an integer
     * @return the power as a new instance
     */
    public ComplexDouble pow (int k){
        double r  = getAbs(), an   = getArg();
        double r1 = Math.pow(r, k), an1 = ((double) k)*an;
        double x = r1 * Math.cos(an1), y = r1 * Math.sin(an1);
        return new ComplexDouble (x,y);
    }

    /**
     * Returns the multiplicative inverse of this instance.
     * @return the inverse, if the absolute value is greater than zero
     */
    public ComplexDouble inverse (){
        if(abs != 0.0){
            double sabs = 1.0/Math.pow(abs.doubleValue(),2.0), an = -arg.doubleValue();
            double x = sabs*Math.cos(an), y = sabs*Math.sin(an);
            return new ComplexDouble (x,y);
        }
        else{
            return new ComplexDouble (Double.NaN, Double.NaN);
        }

    }

    /**
     * Getter: returns the argument field.
     * @return the argument
     */
    public double getArg (){
        return (new Double (arg)).doubleValue();
    }

    /**
     * Getter: returns the absolute value.
     * @return the absolute value
     */
    public double getAbs (){
        return (new Double (abs)).doubleValue();
    }

    /**
     * Getter: returns the real part
     * @return the real part
     */
    public double getRealPart(){
        return (new Double (re_z)).doubleValue();
    }

    /**
     * Getter: returns the imaginary part
     * @return the imaginary part
     */
    public double getImPart (){
        return (new Double (im_z)).doubleValue();
    }

    /**
     * Subtracts the ComplexNumber <tt>z</tt> from this instance and returns the
     * difference as a new instance of this class. So, the result equals:
     * <p>
     * <tt> this - z = result</tt>
     * </p>
     * @param z
     * @return the complex difference
     */
    public ComplexDouble subtract (ComplexDouble z){
        return new ComplexDouble (getRealPart() - z.getRealPart(), getImPart() - z.getImPart());
    }

    /**
     * Getter: overrides canonical <tt>toString()</tt> an
     * returns this instance in real part - imaginary part
     * representation.
     * @return this instance as a String
     */
    @Override
    public String toString (){
        double real = getRealPart(), ima = getImPart();
        if(real != 0.0){
        	if(ima != 0.0 && Math.abs(ima) != 1.0) {
        		if(ima > 0.0) return new String (real + " + "+ima+" i");
        		else return new String (real + " - "+Math.abs(ima)+" i");
        	}
        	else{
        		if(Math.abs(ima) == 1.0){
        			if(ima == 1.0) return new String (real + " + i");
        			else return new String (real + " - i");
        		}
        		else return new String(real+"");
        	}
        }
        else{
        	if(ima != 0.0 && Math.abs(ima) != 1.0) {
        		if(ima > 0.0) return new String (ima+" i");
        		else return new String (" - "+Math.abs(ima)+" i");
        	}
        	else{
        		if(Math.abs(ima) == 1.0){
        			if(ima == 1.0) return new String ("i");
        			else return new String (" - i");
        		}
        		else return new String(real+"");
        	}
        }
        
    }
    
    /**
     * Compares this instance with <tt>z</tt> and returns true if and only if both, the real and
     * imaginary part are equal.
     * @param z
     * @return true, if and only if the real and imaginary part are equal
     */
    public boolean equals (ComplexDouble z){
    	if(getRealPart() == z.getRealPart() && getImPart() == z.getImPart()) return true;
    	else return false;
    }
    
    public static void fillMatrix (ComplexDouble[][] mat){
    	int row = mat.length, col = mat[0].length;
    	for(int i = 0; i < row; i++){
    		for(int j = 0; j < col; j++){
    			if(mat[i][j] == null) mat[i][j] = new ComplexDouble(0.0,0.0);
    		}
    	}
    }
}
