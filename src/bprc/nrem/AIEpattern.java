package bprc.nrem;

import java.util.List;

public class AIEpattern {
	public String[][] AIE(int depth, String[][] value, List<Integer> up, List<Integer> down, List<Integer> unchange) {
	
		String[][] matrix = new String[up.size()*down.size()+up.size()*unchange.size()+down.size()*unchange.size()]
				[value[0].length - depth + 1];
		
		int count=0;
		for(int i=0;i<up.size();i++){
			for(int j=0;j<unchange.size();j++){
				matrix[count][0]=value[up.get(i)][0]+","+value[unchange.get(j)][0];
				for(int m=depth;m<value[0].length;m++){
					if(Double.parseDouble(value[up.get(i)][m])>Double.parseDouble(value[unchange.get(j)][m])){
						matrix[count][m-depth+1]=1+"";
					}else{
						matrix[count][m-depth+1]=0+"";
					}
				}
				count++;
			}
		}
		for(int i=0;i<up.size();i++){
			for(int j=0;j<down.size();j++){
				matrix[count][0]=value[up.get(i)][0]+","+value[down.get(j)][0];
				for(int m=depth;m<value[0].length;m++){
					if(Double.parseDouble(value[up.get(i)][m])>Double.parseDouble(value[down.get(j)][m])){
						matrix[count][m-depth+1]=1+"";
					}else{
						matrix[count][m-depth+1]=0+"";
					}
				}
				count++;
			}
		}
		for(int i=0;i<unchange.size();i++){
			for(int j=0;j<down.size();j++){
				matrix[count][0]=value[unchange.get(i)][0]+","+value[down.get(j)][0];
				for(int m=depth;m<value[0].length;m++){
					if(Double.parseDouble(value[unchange.get(i)][m])>Double.parseDouble(value[down.get(j)][m])){
						matrix[count][m-depth+1]=1+"";
					}else{
						matrix[count][m-depth+1]=0+"";
					}
				}
				count++;
			}
		}
		return matrix;
	}
	
	public String[][] AIE(int depth, String[][] value, List<Integer> change, List<Integer> unchange) {
		
		String[][] matrix = new String[change.size()*(change.size()-1)/2+change.size()*unchange.size()]
				[value[0].length - depth + 1];
		
		int count=0;
		for(int i=0;i<change.size();i++){
			for(int j=i+1;j<change.size();j++){
				matrix[count][0]=value[change.get(i)][0]+","+value[change.get(j)][0];
				for(int m=depth;m<value[0].length;m++){
					if(Double.parseDouble(value[change.get(i)][m])>Double.parseDouble(value[change.get(j)][m])){
						matrix[count][m-depth+1]=1+"";
					}else{
						matrix[count][m-depth+1]=0+"";
					}
				}
				count++;
			}
		}
		for(int i=0;i<change.size();i++){
			for(int j=0;j<unchange.size();j++){
				matrix[count][0]=value[change.get(i)][0]+","+value[unchange.get(j)][0];
				for(int m=depth;m<value[0].length;m++){
					if(Double.parseDouble(value[change.get(i)][m])>Double.parseDouble(value[unchange.get(j)][m])){
						matrix[count][m-depth+1]=1+"";
					}else{
						matrix[count][m-depth+1]=0+"";
					}
				}
				count++;
			}
		}
		return matrix;
	}
}
