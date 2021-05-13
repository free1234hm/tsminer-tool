/*******************************************************************************
 * Copyright 2013 Lars Behnke
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *   http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 ******************************************************************************/

package clustering;

import java.util.ArrayList;
import java.util.List;



public class Hierarchical_clustering {

	public static void Hcluster_clustering(double[][] data, String[] genelist) {
		
		double[][] result = new double[data.length][data[0].length];
		String[] resultlist = new String[data.length];
		double[][] distance = new double[data.length][data.length];
		String[] names = new String[data.length];
		
		for(int i=0;i<data.length;i++){
			names[i] = i+"";
			for(int j=i+1;j<data.length;j++){
				double[] x1 = data[i];
				double[] x2 = data[j];
				double dis = EuclideanDistance.getDistance(x1, x2);
				distance[i][j] = dis;
				distance[j][i] = dis;
			}
		}
		
		ClusteringAlgorithm alg = new DefaultClusteringAlgorithm();
		Cluster cluster = alg.performClustering(distance, names, new AverageLinkageStrategy());
		List<Integer> list = new ArrayList<Integer>();
		getlist(cluster,list);

		for(int i=0;i<list.size();i++){
			result[i] = data[list.get(i)];
			resultlist[i] = genelist[list.get(i)];
		}
		
		for(int i=0;i<list.size();i++){
			data[i] = result[i];
			genelist[i] = resultlist[i];
		}	
	}
	
	public static void getlist(Cluster cluster, List<Integer> list) {

        if (cluster != null) {
        	if(cluster.children == null || cluster.children.size() == 0){
        		list.add(Integer.parseInt(cluster.name));
        	}
            for (Cluster child : cluster.getChildren()) {
                getlist(child, list);
            }
        }
    }
	


	
}
