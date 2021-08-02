package clustering;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;

public class Pearsonr {

	public double pearsonp(Double[] x, double[] mean) {
		double p;
		int m = x.length;
	    double[][] aa = new double[m][2];
	    for (int i = 0; i < m; i++) {
	           aa[i][0] = x[i];
	           aa[i][1] = mean[i];
	    }
	    RealMatrix matrix = new Array2DRowRealMatrix(aa, false);
	    PearsonsCorrelation pc = new PearsonsCorrelation(matrix);
	    double r = pc.getCorrelationMatrix().getEntry(0, 1);
	    double pvalue = pc.getCorrelationPValues().getEntry(0, 1); 
	    //两列数相关系数的显著性，即一组随机数与其中一组的相关系数的绝对值大于等于这两组数的概率
	    if(r>0){
			p=1-pvalue/2; //一组随机数与其中一组的相关系数低于这两组数相关系数的概率
			//a random variable takes a value greater than this correlation coefficient
		}else{
			p=pvalue/2;
		}
		return p;
	}
	
	public double pearsonr(Double[] x, double[] mean) {

		int m = x.length;
	    double[][] aa = new double[m][2];
	    for (int i = 0; i < m; i++) {
	           aa[i][0] = x[i];
	           aa[i][1] = mean[i];
	    }
	    RealMatrix matrix = new Array2DRowRealMatrix(aa, false);
	    PearsonsCorrelation pc = new PearsonsCorrelation(matrix);
	    double r = pc.getCorrelationMatrix().getEntry(0, 1);
		return r;
	}
	public double pearsonr(Double[] x, Double[] mean) {

		int m = x.length;
	    double[][] aa = new double[m][2];
	    for (int i = 0; i < m; i++) {
	           aa[i][0] = x[i];
	           aa[i][1] = mean[i];
	    }
	    RealMatrix matrix = new Array2DRowRealMatrix(aa, false);
	    PearsonsCorrelation pc = new PearsonsCorrelation(matrix);
	    double r = pc.getCorrelationMatrix().getEntry(0, 1);
		return r;
	}
	public static double pearsonr(double[] x, double[] mean) {

		int m = x.length;
	    double[][] aa = new double[m][2];
	    for (int i = 0; i < m; i++) {
	           aa[i][0] = x[i];
	           aa[i][1] = mean[i];
	    }
	    RealMatrix matrix = new Array2DRowRealMatrix(aa, false);
	    PearsonsCorrelation pc = new PearsonsCorrelation(matrix);
	    double r = pc.getCorrelationMatrix().getEntry(0, 1);
		return r;
	}
	/*
	public static double[] pearsonr(double[] x, double[] y) {
	    int m = x.length;
	    double[][] aa = new double[m][2];
	    for (int i = 0; i < m; i++) {
	           aa[i][0] = x[i];
	           aa[i][1] = y[i];
	    }
	    RealMatrix matrix = new Array2DRowRealMatrix(aa, false);
	    PearsonsCorrelation pc = new PearsonsCorrelation(matrix);
	    double r = pc.getCorrelationMatrix().getEntry(0, 1);
	    double pvalue = pc.getCorrelationPValues().getEntry(0, 1);
	    return new double[]{r, pvalue};
	}
	
	public static void main(String[] args) {
		double[] x = {2.0, 3.3, 3.0, 5.4, 4.6};
		double[] mean = {4.0, 4.2, 4.0, 5.0, 5.5};
		
		double p;

		double[] r = pearsonr(x,mean);
		if(r[0]>0){
			p=1-r[1]/2;
		}else{
			p=r[1]/2;
		}
		System.out.println(r[0]);
		System.out.println(p);
	}
	*/
}
