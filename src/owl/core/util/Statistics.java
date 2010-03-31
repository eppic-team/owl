package owl.core.util;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedHashMap;

/**
 * Class to keep some simple statistical methods. Originally created for the pearson function.
 */
public class Statistics {
	
	/**
	 * Returns the mean of the numbers in the given collection.
	 * Note that all values are casted to double before summation.
	 * @param values a collection of numbers (Integer, Double, ...)
	 * @return the mean
	 */
	public static <T extends Number> double mean(Collection<T> values) {
		double sum = 0;
		for(T n:values) {
			sum += (Double) n;
		}
		double result = sum / values.size();
		return result;
	}
	
	/**
	 * Returns the standard deviation of the numbers in the given collection.
	 * Note that all values are casted to double before summation.
	 * @param values a collection of numbers (Integer, Double, ...)
	 * @return the standard deviation
	 */	
	public static <T extends Number> double sd(Collection<T> values) {
		double mean = mean(values);
		double sum = 0;
		for(T n:values) {
			sum += ((Double)n - mean) * ((Double)n - mean);
		}
		double result = Math.sqrt(sum / values.size());
		return result;
	}

	/**
	 * Returns the sample standard deviation (using Bessel's correction) of the numbers in the given collection.
	 * Note that all values are casted to double before summation.
	 * @param values a collection of numbers (Integer, Double, ...)
	 * @return the sample standard deviation
	 */	
	public static <T extends Number> double ssd(Collection<T> values) {
		double mean = mean(values);
		double sum = 0;
		for(T n:values) {
			sum += ((Double)n - mean) * ((Double)n - mean);
		}
		double result = Math.sqrt(sum / (values.size() - 1));
		return result;
	}

	/**
	 * Returns the pearson correlation coeffcient for the input vectors x and y.
	 * If x and y are not the same size, the first n values are used where n=min(len(x),len(y)).
	 * Implementation from http://en.wikipedia.org/wiki/Correlation (as of 5/May/2008).
	 * @param a input vector x
	 * @param a	input vector y
	 * @return the pearson correlation coefficient
	 */
	public static double pearson(double[] x, double[] y) {
		int N = Math.min(x.length, y.length);
		double sum_sq_x = 0;
		double sum_sq_y = 0;
		double sum_coproduct = 0;
		double mean_x = x[0];
		double mean_y = y[0];
		double sweep, delta_x, delta_y;
		for(int i=2; i<=N; i++) {
		    sweep = 1.0 * (i - 1) / i;
		    delta_x = x[i-1] - mean_x;
		    delta_y = y[i-1] - mean_y;
		    sum_sq_x += delta_x * delta_x * sweep;
		    sum_sq_y += delta_y * delta_y * sweep;
		    sum_coproduct += delta_x * delta_y * sweep;
		    mean_x += delta_x / i;
		    mean_y += delta_y / i;
		}
		double pop_sd_x = Math.sqrt( sum_sq_x / N );
		double pop_sd_y = Math.sqrt( sum_sq_y / N );
		double cov_x_y = sum_coproduct / N;
		double correlation = cov_x_y / (pop_sd_x * pop_sd_y);
		return correlation;
	}
	
	/**
	 * Returns the spearman correlation coefficient for the input vectors x and y.
	 * @param x
	 * @param y
	 * @return the spearman correlation coefficient
	 * @throws IllegalArgumentException if arrays of different length
	 */
	public static double spearman(double[] x, double[] y) {
		if (x.length!=y.length) throw new IllegalArgumentException("Arrays given to spearman calculation are of different length");
		return pearson(getRanks(x),getRanks(y));
	}
	
	/**
	 * Gets the (ascending) ranks of the given array of values in a double array 
	 * @param values
	 * @return
	 */
	private static double[] getRanks(double[] values) {
		HashMap<Integer, Double> vals = new HashMap<Integer, Double>();
		for (int i=0;i<values.length;i++) {
			vals.put(i, values[i]);
		}
		LinkedHashMap<Integer, Double> valsSorted = Goodies.sortMapByValue(vals, Goodies.ASCENDING);
		ArrayList<Integer> indices = new ArrayList<Integer>(valsSorted.keySet());
		double[] ranks = new double[indices.size()];
		for (int i=0;i<indices.size();i++) {
			ranks[indices.get(i)] = i+1; 
		}
		return ranks;
		
	}
	
	public static void main(String[] args) {
		// testing method pearson()
		double[] x1 = {1.0,2.0,3.0};
		double[] y1 = {4.0,5.0,6.0};
		System.out.println("x = [1,2,3]");
		System.out.println("y = [4,5,6]");
		System.out.println("corr(x,y) = " + pearson(x1,y1));
		
		double[] x2 = {1.0,2.0,3.0};
		double[] y2 = {6.0,5.0,4.0};
		System.out.println("x = [1,2,3]");
		System.out.println("y = [6,5,4]");
		System.out.println("corr(x,y) = " + pearson(x2,y2));
		
		double[] x3 = {0.9572, 0.4854, 0.8003};
		double[] y3 = {0.1419, 0.4218, 0.9157};
		System.out.println("x = [0.9572, 0.4854, 0.8003]");
		System.out.println("y = [0.1419, 0.4218, 0.9157]");
		System.out.println("pearson corr(x,y)  = " + pearson(x3,y3));	// should be -0.1733
		System.out.println("spearman corr(x,y) = " + spearman(x3, y3)); // should be -0.5
	}
}
