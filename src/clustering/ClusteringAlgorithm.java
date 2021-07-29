package clustering;
import java.util.List;

public interface ClusteringAlgorithm
{

    public Cluster performClustering(double[][] distances, String[] clusterNames,
                                     LinkageStrategy linkageStrategy);

    public Cluster performWeightedClustering(double[][] distances, String[] clusterNames,
                                             double[] weights, LinkageStrategy linkageStrategy);

    public List<Cluster> performFlatClustering(double[][] distances,
    		String[] clusterNames, LinkageStrategy linkageStrategy, Double threshold);
}
