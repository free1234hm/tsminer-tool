package bprc.nrem;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;

public class Biclique {
	public List<List<Integer>> biclique(List<String> list1, List<String> list2, Integer p, Integer q) {
		
		HashSet<String> set = new HashSet<String>();
		 for(int i=0;i<list1.size();i++){
			 set.add(list1.get(i));
		 }
		 for(int i=0;i<list2.size();i++){
			 set.add(list2.get(i));
		 }
		 List<String> Allgene = new ArrayList<String>();
		 Iterator<String> iterator=set.iterator();
	    	while(iterator.hasNext()){
	    		Allgene.add(iterator.next());
			}
	    	
	    //if(p < Allgene.size()/50) p=Allgene.size()/50;
	    if(q < Allgene.size()/10) q=Allgene.size()/10;
	    
		Integer[][] matrix = new Integer[Allgene.size()][Allgene.size()];
		
		for(int i=0;i<Allgene.size();i++){
			for(int j=0;j<Allgene.size();j++){
				matrix[i][j]=0;
			}
		}
		for(int i=0;i<list1.size();i++){
			for(int m=0;m<Allgene.size();m++){
				if(Allgene.get(m).equals(list1.get(i))){
					for(int n=0;n<Allgene.size();n++){
						if(Allgene.get(n).equals(list2.get(i))){
							matrix[m][n]=1;
							break;
						}
					}
					break;
				}
			}
		}

		double LR = Math.pow(Math.E, 1.92);
		
		boolean bagaindelete;
		/****************delete all noise gene***************/
		do{
			bagaindelete = false;
			for(int i=0;i<matrix.length;i++){
				int sumx = sumX(i,matrix);
				int sumy = sumY(i,matrix);
				if(sumx<p && sumy<p && (sumx>0 || sumy>0)){
					bagaindelete = true;
					for(int j=0;j<matrix.length;j++){
						matrix[j][i] = 0;
					}
					for(int j=0;j<matrix[0].length;j++){
						matrix[i][j] = 0;
					}
				}
			}
		}while(bagaindelete);
		/***********************************/
		
		List<Integer> result1 = new ArrayList<Integer>();
		List<Integer> result2 = new ArrayList<Integer>();
		for(int i=0;i<matrix.length;i++){
			int sumx = sumX(i,matrix);
			int sumy = sumY(i,matrix);
			if(sumx > sumy*LR){
				result1.add(i);
			}else{
				if(sumy > sumx*LR){					
				result2.add(i);
				}		
			}
		}
		
		
		if(result1.size() >= q && result2.size() >= q){
			
			List<Integer> genelist1 = new ArrayList<Integer>();
			List<Integer> genelist2 = new ArrayList<Integer>();
			for(int i=0;i<result1.size();i++){
				genelist1.add(Integer.parseInt(Allgene.get(result1.get(i))));
			}
			for(int i=0;i<result2.size();i++){
				genelist2.add(Integer.parseInt(Allgene.get(result2.get(i))));
			}
			List<List<Integer>> AIE = new ArrayList<List<Integer>>();
			AIE.add(genelist1);
			AIE.add(genelist2);			
			return AIE;
		}else{		
			return null;
		}
		
	}
	
	public Integer sumX(Integer X, Integer[][] matrix) {
		int result=0;
		for(int i=0;i<matrix[0].length;i++){
			result+=matrix[X][i];
		}		
		return result;
	}
	
	public Integer sumY(Integer Y, Integer[][] matrix) {
		int result=0;
		for(int i=0;i<matrix.length;i++){
			result+=matrix[i][Y];
		}		
		return result;
	}
	
	public static Integer sup(List<Integer> X, Integer[][] matrix) {
		int result=0;
		for(int i=0;i<matrix[0].length;i++){
			int a=0;
			for(int j=0;j<X.size();j++){
				a+=matrix[X.get(j)][i];
			}
			if(a==X.size()) result++;
		}		
		return result;
	}
	
	public static List<Integer> list2(List<Integer> X, Integer[][] matrix) {
		List<Integer> result = new ArrayList<Integer>();
		for(int i=0;i<matrix[0].length;i++){
			int a=0;
			for(int j=0;j<X.size();j++){
				a+=matrix[X.get(j)][i];
			}
			if(a==X.size()) result.add(i);
		}		
		return result;
	}
	
	public static List<Integer> tail(List<Integer> X, List<Integer> Allgene) {
		List<Integer> tail = new ArrayList<Integer>();
		for(int i=0;i<Allgene.size();i++){
				if(X.size()>0){
					if(Allgene.get(i)>X.get(X.size()-1)) tail.add(Allgene.get(i));
				}else{
					tail.add(Allgene.get(i));
				}	
		}		
		return tail;
	}
	
}
