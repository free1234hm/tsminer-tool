package bprc.nrem;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;



public class Negative_correlation {
	
	
    public HashMap<String, List<List<Integer>>> AIE_pattern(String[][] express, List<List<Integer>> upgenes, 
    	List<List<Integer>> downgenes, Integer Maxcount, Integer p, Integer q){

	if(express.length>Maxcount){
		express = Geneselect(express, Maxcount);	
	}
	
	HashMap<String, List<List<Integer>>> AIE = new HashMap<String, List<List<Integer>>>();
	HashMap<String, List<String>> treegene = new HashMap<String, List<String>>();
	List<Integer> up = new ArrayList<Integer>();
	List<Integer> down = new ArrayList<Integer>();
	List<Integer> unchange = new ArrayList<Integer>();
	
	for(int j=2;j<express[0].length-1;j++){
		up.clear();
		down.clear(); 
		unchange.clear();
		
		for(int i=0;i<express.length;i++){
			if(upgenes.get(j-2).contains(Integer.parseInt(express[i][0]))){
				up.add(i);
			}else{
				if(downgenes.get(j-2).contains(Integer.parseInt(express[i][0]))){
					down.add(i);
				}else{
					unchange.add(i);
				}
			}
		}

		if((up.size()>=Math.max(p, q) && (down.size()+unchange.size())>=Math.max(p, q)) 
				|| (down.size()>=Math.max(p, q) && (unchange.size()+up.size())>=Math.max(p, q)) 
				){
		
			AIEpattern aie = new AIEpattern();
			String[][] dichotomy = aie.AIE(j, express, up, down, unchange);
			
				for(int i=0;i<dichotomy.length;i++){
					StringBuffer node = new StringBuffer("");
	    			String genepair = ""; 
	    			if(dichotomy[i][1].equals(0+"")){
	    				for(int m=1;m<dichotomy[0].length;m++){
	        				node.append(dichotomy[i][m]);
	        			}
	    				genepair = dichotomy[i][0];
					}else{
						for(int m=1;m<dichotomy[0].length;m++){
	        				node.append(1-Integer.parseInt(dichotomy[i][m]));
	        			}
						genepair = Arrays.asList(dichotomy[i][0].split(",")).get(1)+","+Arrays.asList(dichotomy[i][0].split(",")).get(0);
					}
	    						
					List<String> genes = new ArrayList<String>();
					genes = treegene.get(node.toString());
					if(genes != null){
					  genes.add(genepair);
					}else{
						genes = new ArrayList<String>();
						genes.add(genepair);
					}
					treegene.put(node.toString(), genes);
				}
		}
	}

    	Iterator iter = treegene.entrySet().iterator();
    	while (iter.hasNext()) {
    	Map.Entry entry = (Map.Entry) iter.next();
    	Object key = entry.getKey();
    	List<String> genes = treegene.get(key);
    		
    		
	    	List<String> genelist1 = new ArrayList<String>();
	    	List<String> genelist2 = new ArrayList<String>();
	    	for(int i=0;i<genes.size();i++){
				genelist1.add(Arrays.asList(genes.get(i).split(",")).get(0));
				genelist2.add(Arrays.asList(genes.get(i).split(",")).get(1));
			}
				Biclique bic = new Biclique();
				List<List<Integer>> matrix = bic.biclique(genelist1, genelist2, p, q);
				if(matrix != null && matrix.size()>0){
					/****************Assign the rest genes based on correlation*******************/
					CorrelatedIndex corr = new CorrelatedIndex();
					List<List<Integer>> matrix2 = corr.corr(express, key.toString().length(), matrix.get(0), matrix.get(1));
					AIE.put(key.toString(), matrix2);
				}	
    	}	  	
    return AIE;
	}
	
    public HashMap<String, List<List<Integer>>> AIE_pattern(String[][] express, List<List<Integer>> DEgenes, 
        	Integer Maxcount, Integer p, Integer q){

    	if(express.length>Maxcount){
    		express = Geneselect(express, Maxcount);	
    	}
    	
    	HashMap<String, List<List<Integer>>> AIE = new HashMap<String, List<List<Integer>>>();
    	HashMap<String, List<String>> treegene = new HashMap<String, List<String>>();
    	List<Integer> change = new ArrayList<Integer>();
    	List<Integer> unchange = new ArrayList<Integer>();
    	
    	for(int j=2;j<express[0].length-1;j++){
    		change.clear();
    		unchange.clear();
    		
    		for(int i=0;i<express.length;i++){
    			if(DEgenes.get(j-2).contains(Integer.parseInt(express[i][0]))){
    				change.add(i);
    			}else{
    				unchange.add(i);
    			}
    		}

    		if(change.size()>=Math.max(p, q) && unchange.size()>=Math.max(p, q)){
    		
    			AIEpattern aie = new AIEpattern();
    			String[][] dichotomy = aie.AIE(j, express, change, unchange);
    			
    				for(int i=0;i<dichotomy.length;i++){
    					StringBuffer node = new StringBuffer("");
    	    			String genepair = ""; 
    	    			if(dichotomy[i][1].equals(0+"")){
    	    				for(int m=1;m<dichotomy[0].length;m++){
    	        				node.append(dichotomy[i][m]);
    	        			}
    	    				genepair = dichotomy[i][0];
    					}else{
    						for(int m=1;m<dichotomy[0].length;m++){
    	        				node.append(1-Integer.parseInt(dichotomy[i][m]));
    	        			}
    						genepair = Arrays.asList(dichotomy[i][0].split(",")).get(1)+","+Arrays.asList(dichotomy[i][0].split(",")).get(0);
    					}
    	    						
    					List<String> genes = new ArrayList<String>();
    					genes = treegene.get(node.toString());
    					if(genes != null){
    					  genes.add(genepair);
    					}else{
    						genes = new ArrayList<String>();
    						genes.add(genepair);
    					}
    					treegene.put(node.toString(), genes);
    				}
    		}
    	}

        	Iterator iter = treegene.entrySet().iterator();
        	while (iter.hasNext()) {
        	Map.Entry entry = (Map.Entry) iter.next();
        	Object key = entry.getKey();
        		List<String> genes = treegene.get(key);
        		
        		
    	    	List<String> genelist1 = new ArrayList<String>();
    	    	List<String> genelist2 = new ArrayList<String>();
    	    	for(int i=0;i<genes.size();i++){
    				genelist1.add(Arrays.asList(genes.get(i).split(",")).get(0));
    				genelist2.add(Arrays.asList(genes.get(i).split(",")).get(1));
    			}
    				Biclique bic = new Biclique();
    				List<List<Integer>> matrix = bic.biclique(genelist1, genelist2, p, q);
    				if(matrix != null && matrix.size()>0){
    					CorrelatedIndex corr = new CorrelatedIndex();
    					List<List<Integer>> matrix2 = corr.corr(express, key.toString().length(), matrix.get(0), matrix.get(1));
    					AIE.put(key.toString(), matrix2);
    				}	
        	}	  	
        return AIE;
    	}
    
	public static String[][] Geneselect (String[][] express, Integer count) {
		String[][] result = new String[count][express[0].length];
		String[][] distance = new String[express.length][2];
		for(int i=0;i<express.length;i++){
		distance[i][0] = express[i][0];
		double min=0;
		double max=0;
		for(int j=1;j<express[0].length;j++){
			if(Double.parseDouble(express[i][j]) > max) max = Double.parseDouble(express[i][j]);
			if(Double.parseDouble(express[i][j]) < min) min = Double.parseDouble(express[i][j]);
		}
		distance[i][1] = (max-min)+"";
		}
		String distance2[][] = BubbleSort(distance, distance.length); //冒泡排序
		for(int i=0; i<count; i++){
			for(int j=0; j<express.length; j++){
				if(distance2[i][0].equals(express[j][0])){
					result[i] = express[j];
					break;
				}
			}
		}
		return result;
	}
	static String[][] BubbleSort(String[][] r, Integer n) //降序冒泡排序
	{
		 int low = 0;   
		    int high= n -1; //设置变量的初始值  
		    String[] tmp;
		    int j;  
		    while (low < high) {  
		        for (j= low; j< high; ++j) //正向冒泡,找到最大者  
		            if (Double.parseDouble(r[j][1])< Double.parseDouble(r[j+1][1])) {  
		                tmp = r[j]; r[j]=r[j+1];r[j+1]=tmp;  
		            }   
		        --high;                 //修改high值, 前移一位  
		        for ( j=high; j>low; --j) //反向冒泡,找到最小者  
		            if (Double.parseDouble(r[j][1])>Double.parseDouble(r[j-1][1])) {  
		                tmp = r[j]; r[j]=r[j-1];r[j-1]=tmp;  
		            }  
		        ++low; //修改low值,后移一位  
		    }   
	return r;
	}
	
}

