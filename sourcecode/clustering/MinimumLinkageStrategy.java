package clustering;

import java.util.Collection;

//TODO Not working correctly, fix
public class MinimumLinkageStrategy implements LinkageStrategy {

	@Override
	public Distance calculateDistance(Collection<Distance> distances) {
		double min = Double.MAX_VALUE;
		double result;

		for (Distance dist : distances) {
			if(min>dist.getDistance()) {
				min = dist.getDistance();
			}
		}
		if (distances.size() > 0) {
			result = min;
		} else {
			result = 0.0;
		}
		return  new Distance(result);
	}
}