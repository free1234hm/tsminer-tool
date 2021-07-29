package bprc.nrem;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map.Entry;

import org.apache.commons.math3.distribution.NormalDistribution;

public class MIM {
	
		
		public Integer[] discretize(double[] X){
			int nrow = X.length;
			int nbins = (int)nrow/2;
			Integer[] res = new Integer[nrow];
			double epsilon = 0.01;
			List<Double> col = new ArrayList<Double>();
			for(int i=0;i<nrow;i++) col.add(X[i]);
			Collections.sort(col);
			
			int freq = nrow/nbins;
			int mod = nrow%nbins;
			int splitpoint=freq-1;
			
			Double[] spl = new Double[nbins];
			for(int i=0;i<nbins-1;i++) {
				if( mod>0 ){
					spl[i] = col.get(splitpoint+1);
					mod--;
				}else{
					spl[i] = col.get(splitpoint);
				}
				splitpoint += freq;
			}
			spl[nbins-1] = col.get(nrow-1)+epsilon;
			
			for(int s=0;s<nrow;s++){
				int bin = -1;
		        for(int k=0;k<nbins;k++){
		        	if(X[s]<=spl[k]){
		        		bin=k;
		        		break;
		        	}else{
		        		bin++;
		        	}
		        }
		        res[s] = bin+1;
			}
			return res;
		}
		
	    public Double Entropy(Integer[] X, Integer[] Y, Integer[] Z) {
	    	int sum = X.length;    	
	    	double e = .0;
	        double H = .0;
	        HashMap<String, Integer> freq = new HashMap<String, Integer>();
	        if(Z!=null){
	        	for(int i=0;i<X.length;i++){
	               	String a = X[i]+","+Y[i]+","+Z[i];	
	               	if(freq.get(a) != null){
	               		int count = freq.get(a)+1;
	               		freq.put(a, count);
	               	}else{
	               		freq.put(a, 1);
	               	}
	               } 
	        }else{
	        	if(Y==null){
	           	 for(int i=0;i<X.length;i++){
	                	String a = X[i]+"";	
	                	if(freq.get(a) != null){
	                		int count = freq.get(a)+1;
	                		freq.put(a, count);
	                	}else{
	                		freq.put(a, 1);
	                	}
	                }
	           }else{
	               for(int i=0;i<X.length;i++){
	               	String a = X[i]+","+Y[i];	
	               	if(freq.get(a) != null){
	               		int count = freq.get(a)+1;
	               		freq.put(a, count);
	               	}else{
	               		freq.put(a, 1);
	               	}
	               }  	
	           }
	        }
	        
	        for(Entry<String, Integer> entry:freq.entrySet())
			  {
			   Integer val = entry.getValue();
			   e += -((double)val * Math.log((double)val));
			  }
	    	H=Math.log(sum)+e/(sum);	
	        	
	        return H;  
	    }
	    
	    public Double condinformation(Integer[] X, Integer[] Y, Integer[] S) {
	        double H = .0;
	    	
	    	 if(S==null){
	    	     double Hyx = Entropy(Y,X,null);
	    		 double Hx = Entropy(X,null,null);
	    		 double Hy = Entropy(Y,null,null);
	    	     H = Hx + Hy - Hyx;
	    	     if (H < 0) H=0;       
	    	 }else{
	    		 double Hsxy = Entropy(S,X,Y);
	    		 double Hsx = Entropy(S,X,null);
	    	     double Hsy = Entropy(S,Y,null);
	    	     double Hs = Entropy(S,null,null);
	    	     H = Hsy - Hs - Hsxy + Hsx;
	    	 }
	    	
	    	 return H;
	    }
	    
	    public Double getPvalue(double[] pw1, double[] pw2, double[] tf, Double dmean, Double dvar) {
	    	 Integer[] y1 = discretize(pw1);
	         Integer[] y2 = discretize(pw2);
	         Integer[] x = discretize(tf);
	         double pvalue = 0.0;
	         
	         double Hy1 = Entropy(y1,null,null);
	         double Hy2 = Entropy(y2,null,null);
	         double Hx = Entropy(x,null,null);
	         double Hxy1 = Entropy(x,y1,null);
	         double Hxy2 = Entropy(x,y2,null);
	         double Hy2y1 = Entropy(y2,y1,null);
	         
	         double Hxy2y1 = Entropy(x,y2,y1);
	         double enty1_y2x = Hxy2y1-Hxy2;
	         double enty2_y1x = Hxy2y1-Hxy1;
	         double entx_y1y2 = Hxy2y1-Hy2y1;
	         double Iy1y2 = Hy2+Hy1-Hy2y1;
	         double Iy1y2_x = Hxy1+Hxy2-Hxy2y1-Hx;
	         double denominator = enty1_y2x+enty2_y1x+entx_y1y2;
	         double Iy1y2x = Iy1y2 - Iy1y2_x;
	        
	         if(denominator==0){
	        	 if(Iy1y2x>0){
	        		 pvalue = 0.0;
	        	 }else{
	        		 pvalue = -1.0;
	        	 }
	         }else{
	        	 double I3 = Iy1y2x/denominator;
	        	 Double zscore = (I3-dmean)/dvar;
	        	 NormalDistribution nor = new NormalDistribution();
	        	 pvalue = 1-nor.cumulativeProbability(zscore);
	         }
	         return pvalue;
	    }
	    
	    public Double[] random(double[] x1, double[] x2, double[] x3) {
	    	double[] pvalue = new double[1000];
	    	int N = x1.length;
	    	
	    	Integer[] y1 = discretize(x1);
	    	Integer[] y2 = discretize(x2);
	    	Integer[] x = discretize(x3);

	    	List<Integer> e1 = new ArrayList<Integer>();
	    	List<Integer> e2 = new ArrayList<Integer>();
	    	List<Integer> e3 = new ArrayList<Integer>();
	    	for(int i=0;i<N;i++) {
	    		e1.add(y1[i]);
	    		e2.add(y2[i]);
	    		e3.add(x[i]);
	    	}
	    	
	    	Integer[] Y1 = new Integer[N];
			Integer[] Y2 = new Integer[N];
			Integer[] X = new Integer[N];
			
	    	for(int i=0;i<pvalue.length;i++){
	    		Collections.shuffle(e1);
	    		Collections.shuffle(e2);
	    		Collections.shuffle(e3);
	    		for(int j=0;j<N;j++){
	    			Y1[j] = e1.get(j);
	    			Y2[j] = e2.get(j);
	    			X[j] = e3.get(j);
	    		}
	    		pvalue[i] = mi3(Y1,Y2,X);
	    	}
	    	
	    	double sum=0;
	    	int count=0;
	    	for(int i=0;i<pvalue.length;i++){
	    		if(pvalue[i] != Double.POSITIVE_INFINITY){
	    			sum += pvalue[i];
	    			count++;
	    		}
	    	}
	    	double dmean = sum/count;
	    	double dVar=0;  
	        for(int i=0;i<pvalue.length;i++){
	        	if(pvalue[i] != Double.POSITIVE_INFINITY){
	        		dVar+=(pvalue[i]-dmean)*(pvalue[i]-dmean); 
	        	}
	        }
	    	double var = Math.sqrt(dVar/(count-1));
	    	Double[] result = {dmean, var};
	    	return result;
	    }
	    public Double mi3(Integer[] y1, Integer[] y2, Integer[] x) {
	        double I3 = 0.0;
	        
	        double Hy1 = Entropy(y1,null,null);
	        double Hy2 = Entropy(y2,null,null);
	        double Hx = Entropy(x,null,null);
	        double Hxy1 = Entropy(x,y1,null);
	        double Hxy2 = Entropy(x,y2,null);
	        double Hy2y1 = Entropy(y2,y1,null);
	        
	        double Hxy2y1 = Entropy(x,y2,y1);
	        double enty1_y2x = Hxy2y1-Hxy2;
	        double enty2_y1x = Hxy2y1-Hxy1;
	        double entx_y1y2 = Hxy2y1-Hy2y1;
	        double Iy1y2 = Hy2+Hy1-Hy2y1;
	        double Iy1y2_x = Hxy1+Hxy2-Hxy2y1-Hx;
	        double denominator = enty1_y2x+enty2_y1x+entx_y1y2;
	        double Iy1y2x = Iy1y2 - Iy1y2_x;
	 
	        if(denominator == 0){
	       	 I3 = Double.POSITIVE_INFINITY;
	        }else{
	       	 I3 = Math.abs(Iy1y2x/denominator);	 
	        }
	        return I3;
	   }
	    
}
