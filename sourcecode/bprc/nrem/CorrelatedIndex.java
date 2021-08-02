package bprc.nrem;

import java.util.ArrayList;
import java.util.List;

public class CorrelatedIndex {
	public List<List<Integer>> corr(String[][] express, int length, List<Integer> list1, List<Integer> list2){
		List<List<Integer>> genelist = new ArrayList<List<Integer>>();
		List <Integer> resultlist1 = new ArrayList<Integer>();
		List <Integer> resultlist2 = new ArrayList<Integer>();
		resultlist1.addAll(list1);
		resultlist2.addAll(list2);
		
		int column = express[0].length;
		List <Double> xlist = new ArrayList<Double>();
		List <Double> ylist = new ArrayList<Double>();
		for(int j=column-length;j<column;j++){
			double mean=0;
			for(int i=0;i<list1.size();i++){
				for(int m=0;m<express.length;m++){
					if(Integer.parseInt(express[m][0]) == list1.get(i)){
						mean += Double.parseDouble(express[m][j]);
						break;
					}
				}
			}
			mean = mean/list1.size();	
			xlist.add(mean);
		}
		for(int j=column-length;j<column;j++){
			double mean=0;
			for(int i=0;i<list2.size();i++){
				for(int m=0;m<express.length;m++){
					if(Integer.parseInt(express[m][0]) == list2.get(i)){
						mean += Double.parseDouble(express[m][j]);
						break;
					}
				}
			}
			mean = mean/list2.size();	
			ylist.add(mean);
		}
		
		for(int i=0;i<express.length;i++){
			if(!list1.contains(Integer.parseInt(express[i][0])) && !list2.contains(Integer.parseInt(express[i][0]))){
				List <Double> list = new ArrayList<Double>();
				for(int j=column-length;j<column;j++){
					list.add(Double.parseDouble(express[i][j]));	
				}
				double molecular1 = molecular(list, xlist);
				double denominator1 = denominator(list,xlist);
				double cor1 = molecular1/denominator1;
				double molecular2 = molecular(list, ylist);
				double denominator2 = denominator(list,ylist);
				double cor2 = molecular2/denominator2;
				if(cor1 > cor2){
					resultlist1.add(Integer.parseInt(express[i][0]));
				}else{
					resultlist2.add(Integer.parseInt(express[i][0]));
				}
			}
		}
		
		genelist.add(resultlist1);
		genelist.add(resultlist2);
		
		return genelist;
	}
	
	public double molecular(List<Double> xList,List<Double> yList){  
        double result =0.0;  
        double xAverage = 0.0;  
        double temp = 0.0;  
          
        int xSize = xList.size();  
        for(int x=0;x<xSize;x++){  
            temp += xList.get(x);  
        }  
        xAverage = temp/xSize;  
          
        double yAverage = 0.0;  
        temp = 0.0;  
        int ySize = yList.size();  
        for(int x=0;x<ySize;x++){  
            temp += yList.get(x);  
        }  
        yAverage = temp/ySize;  
          
        //double sum = 0.0;  
        for(int x=0;x<xSize;x++){  
            result+=(xList.get(x)-xAverage)*(yList.get(x)-yAverage);  
        }  
        return result;  
    } 
	        
	    public double denominator(List<Double> xList,List<Double> yList){  
	        double standardDifference = 0.0;  
	        int size = xList.size();  
	        double xAverage = 0.0;  
	        double yAverage = 0.0;  
	        double xException = 0.0;  
	        double yException = 0.0;  
	        double temp = 0.0;  
	        for(int i=0;i<size;i++){  
	            temp += xList.get(i);  
	        }  
	        xAverage = temp/size;  
	         
	        temp=0;
	        for(int i=0;i<size;i++){  
	            temp += yList.get(i);  
	        }  
	        yAverage = temp/size;  
	          
	        for(int i=0;i<size;i++){  
	            xException += Math.pow(xList.get(i)-xAverage,2);  
	            yException += Math.pow(yList.get(i)-yAverage,2);  
	        }  
	        //calculate denominator of   
	        return standardDifference = Math.sqrt(xException*yException);  
	    }   
		
}
