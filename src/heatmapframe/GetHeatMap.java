package heatmapframe;

import java.awt.*;
import java.util.List;


/**
 * <p>This class is a very simple example of how to use the HeatMap class.</p>
 *
 * <hr />
 * <p><strong>Copyright:</strong> Copyright (c) 2007, 2008</p>
 *
 * <p>HeatMap is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.</p>
 *
 * <p>HeatMap is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.</p>
 *
 * <p>You should have received a copy of the GNU General Public License
 * along with HeatMap; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA</p>
 *
 * @author Matthew Beckler (matthew@mbeckler.org)
 * @author Josh Hayes-Sheen (grey@grevian.org), Converted to use BufferedImage.
 * @author J. Keller (jpaulkeller@gmail.com), Added transparency (alpha) support, data ordering bug fix.
 * @version 1.6
 */

public class GetHeatMap{

    public HeatMap heatmap(double[][] data, String[] dsamplemins, String[] genelist) throws Exception
    {
        //Gradient gra = new Gradient();
        data = transpose(data);
        boolean useGraphicsYAxis = true;
        
        //Color[] gradientColors = new Color[]{Color.BLUE, Color.white, Color.red};
        Color[] gradientColors = new Color[]{Color.BLUE, Color.YELLOW};
        Color[] customGradient = Gradient.createMultiGradient(gradientColors, 40);
        HeatMap panel = new HeatMap(data, dsamplemins, genelist, useGraphicsYAxis, customGradient);
        
        // set miscelaneous settings
        panel.setDrawLegend(true);

        panel.setTitle("TF Heatmap");
        panel.setDrawTitle(false);

        panel.setXAxisTitle("Time Points");
        panel.setDrawXAxisTitle(false);

        panel.setYAxisTitle("Gene List");
        panel.setDrawYAxisTitle(false);

        panel.setCoordinateBounds(0, data.length, 0, data[0].length);
        panel.setDrawXTicks(true);
        panel.setDrawYTicks(false);

        panel.setColorForeground(Color.black);
        panel.setColorBackground(Color.white);

        return panel;
    }
    
    public double[][] transpose(double[][] a){
		double b[][] = new double[a[0].length][a.length];
		for(int i=1;i<=b.length;i++){
			for(int j=1;j<=b[0].length;j++){
				b[i-1][j-1] = a[j-1][i-1];
			}
		}
		return b;
	}
}

