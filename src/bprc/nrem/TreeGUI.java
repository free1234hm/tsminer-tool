package bprc.nrem;

import heatmapframe2.DrawMainHeatMap;
import java.awt.Color;
import java.awt.Container;
import java.awt.Dimension;
import java.util.ArrayList;
import java.util.List;
import javax.swing.JFrame;
import bprc.nrem.TSMiner_Timeiohmm.Treenode;

public class TreeGUI {
	
	List<Integer> generank;
	int[] heatmap2gene;
	int[] gene2heatmap;

	
	public TreeGUI(TSMiner_Timeiohmm theTimeiohmm, Treenode treecopy, int npermutationval) {
		
		JFrame frame = new JFrame("TSMiner - Main Interface");
		frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		frame.setLocation(10, 10);
		frame.setMinimumSize(new Dimension(500, 0));
		Container theDialogContainer = frame.getContentPane();
		theDialogContainer.setBackground(Color.white);
		
		double[][] data = theTimeiohmm.theDataSet.data;
		String[] genename = theTimeiohmm.theDataSet.genenames;
		String[] dsamplemins = new String[theTimeiohmm.theDataSet.dsamplemins.length-1];
		for(int i=0;i<dsamplemins.length;i++) dsamplemins[i] = theTimeiohmm.theDataSet.dsamplemins[i+1];
		
		generank = new ArrayList<Integer>();
		double[][] normdata = new double[data.length][data[0].length-1];
        
        getgenerank(treecopy);
        
        heatmap2gene = new int[generank.size()];
        gene2heatmap = new int[generank.size()];
        
        for(int i=0;i<generank.size();i++){
        	heatmap2gene[i] = generank.get(i);
        	gene2heatmap[generank.get(i)] = i;
        	double[] expr = data[generank.get(i)];
        	double largest = Double.MIN_VALUE;
        	for(int j=0;j<expr.length;j++){
        		largest = Math.max(Math.abs(expr[j]), largest);
        	}
        	
        	for(int j=1;j<expr.length;j++){
        		normdata[i][j-1] = (expr[j]/largest)/2+0.5;
        	}
        	
        }
        getgrid(treecopy, 0);
        
        System.out.println("Begain to draw heatmap!");
        DrawMainHeatMap newContentPane1 = new DrawMainHeatMap(theTimeiohmm, treecopy, npermutationval, 
        		dsamplemins, genename, normdata, heatmap2gene);
        theDialogContainer.add(newContentPane1);
		frame.setContentPane(theDialogContainer);
		frame.pack();
		frame.setVisible(true);
	}
	
	public void getgenerank(Treenode ptr) {
		if (ptr.numchildren == 0) {
			generank.addAll(ptr.genelist);
		}
		for (int nchild = 0; nchild < ptr.numchildren; nchild++) {
			getgenerank(ptr.nextptr[nchild]);
		}
	}

	public void getgrid(Treenode ptr, int child) {
		if(ptr != null){	
			int ndepth = ptr.ndepth;
			
			if(ndepth>0){
				int min = Integer.MAX_VALUE;
				int max = Integer.MIN_VALUE;
				for(int g:ptr.genelist){
					int heatmapindex = gene2heatmap[g];
					min = Math.min(heatmapindex, min);
					max = Math.max(heatmapindex, max);
				}
				ptr.heatmaprange = new int[2];
				ptr.heatmaprange[0] = min;
				ptr.heatmaprange[1] = max;
			}
			
			
			for (int nchild = 0; nchild < ptr.numchildren; nchild++) {
				getgrid(ptr.nextptr[nchild], nchild);
			}
		}
	}
	
	

}
