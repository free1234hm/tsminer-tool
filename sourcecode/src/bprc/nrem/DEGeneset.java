package bprc.nrem;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.distribution.NormalDistribution;


public class DEGeneset {
	
	public Double[] permutetest(double[] express){
		if(express.length>2){
			Double[] result = new Double[2];
			int length = express.length/2;
			
			List<Integer> list = new ArrayList<Integer>();
			for(int i=0;i<express.length;i++) list.add(i);
			
			double descore = 0;
			for(int i=0;i<length;i++){
				descore += (express[i+length]-express[i])/length;
			}
			result[1] = descore;
			int rancount = rancount(express.length, length);
			double[] random = new double[rancount];
			double sum=0;
			for(int ran=0;ran<rancount;ran++){
				double ranscore = 0;
				Collections.shuffle(list);
				
				for(int i=0;i<length;i++){
					ranscore += (express[list.get(i+length)]-express[list.get(i)])/length;
				}
				random[ran] = Math.abs(ranscore);
				sum += Math.abs(ranscore);
			}
			double dmean = sum/rancount;
	    	double dVar=0;  
	        for(int i=0;i<random.length;i++){
	           dVar+=(random[i]-dmean)*(random[i]-dmean); 
	        }
	        double var = Math.sqrt(dVar/(rancount-1));
	        
	        Double zscore = (Math.abs(descore)-dmean)/var;
	   	    NormalDistribution nor = new NormalDistribution();
	   	    result[0] = 1-nor.cumulativeProbability(zscore);
	   	    return result;
		}else{
			Double[] result = new Double[2];
			int length = express.length/2;	
			double descore = 0;
			for(int i=0;i<length;i++){
				descore += (express[i+length]-express[i])/length;
			}
			result[1] = descore;
			return result;
		}
	}
	
	public double[] permutetest(List<Double> express, int length1, int length2){
		double[] result = new double[2];
		int length = express.size();
		
		List<Integer> list = new ArrayList<Integer>();
		for(int i=0;i<express.size();i++) list.add(i);
		
		double descore1 = 0;
		for(int i=0;i<length1;i++){
			descore1 += express.get(i);
		}
		double avg1 = descore1/length1;
		
		double descore2 = 0;
		for(int i=length1;i<length;i++){
			descore2 += express.get(i);
		}
		double avg2 = descore2/length2;
		double descore = avg2-avg1;
		result[1] = descore;
		int rancount = rancount(length1+length2,length1);
		double[] random = new double[rancount];
		double sum=0;
		for(int ran=0;ran<rancount;ran++){
			Collections.shuffle(list);
			
			double ranscore1 = 0;
			for(int i=0;i<length1;i++){
				ranscore1 += express.get(list.get(i));
			}
			double ranavg1 = ranscore1/length1;
			
			double ranscore2 = 0;
			for(int i=length1;i<length;i++){
				ranscore2 += express.get(list.get(i));
			}
			double ranavg2 = ranscore2/length2;
		
			double ranscore = ranavg2-ranavg1;

			random[ran] = Math.abs(ranscore);
			sum += Math.abs(ranscore);
		}
		double dmean = sum/rancount;
    	double dVar=0;  
        for(int i=0;i<random.length;i++){
           dVar+=(random[i]-dmean)*(random[i]-dmean); 
        }
        double var = Math.sqrt(dVar/(rancount-1));
        
        Double zscore = (Math.abs(descore)-dmean)/var;
   	    NormalDistribution nor = new NormalDistribution();
   	    result[0] = 1-nor.cumulativeProbability(zscore);
        
		return result;
	}
	
	
	public List<Double> ranES(int[] set, int[][] Random, String[][] FC){
		List<Double> randomES = new ArrayList<Double>();

		   for(int random=0;random<Random[0].length;random++){
			   List<Integer> ranlist = new ArrayList<Integer>();
			   double[][] loc=new double[Random.length][1];

			   double ranES=0;
			   double ranNR=0;
			   double ranNM=0;
			   for(int p=0;p<Random.length;p++){
				   ranlist.add(Random[p][random]);
			   }
			   for(int n=0;n<set.length;n++){
			   for(int m=0;m<ranlist.size();m++){
					   if(ranlist.get(m) == set[n]){
						   ranNR=ranNR+Math.abs(Double.parseDouble(FC[m][1]));
						   ranNM=ranNM+1;
						   loc[m][0]=Math.abs(Double.parseDouble(FC[m][1]));
						   break;
					   }
				   }
			   }
			   double ranMiss=loc.length-ranNM;
			   for(int i=0;i<loc.length;i++){
				   double ranhit=0;
				   double ranmiss=0;
				   for(int j=0;j<=i;j++){
					   if(loc[j][0]==0){
						   ranmiss=ranmiss+1/ranMiss;					   
					   }else{
						   ranhit=ranhit+loc[j][0]/ranNR;
					   }
					   
				   }
						  			   
				   if(ranES<Math.abs(ranhit-ranmiss)){
					   ranES=Math.abs(ranhit-ranmiss);
				   }
			   }
			   randomES.add(ranES);
			   ranlist.clear();
			   //location.clear();
		   }
		
		return randomES;	
	}
	
	public int rancount(int len, int get){
		double aa = 1;
		for(int i=len;i>get;i--){
			aa = aa*i/(i-get);
		}
		if(aa>10000){
			return 10000;
		}else{
			return (int)aa;
		}
	}

}
