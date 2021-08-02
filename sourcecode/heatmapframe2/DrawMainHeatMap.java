package heatmapframe2;

import javax.swing.BoxLayout;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;

import java.awt.Color;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.Frame;
import java.util.*;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.awt.event.*;

import bprc.core.*;

import javax.swing.*;
import javax.swing.border.TitledBorder;

import bprc.nrem.TSMinerGui_SaveModel;
import bprc.nrem.TSMinerGui_TFTable;
import bprc.nrem.TSMiner_Timeiohmm;
import bprc.nrem.TSMiner_Timeiohmm.Treenode;


/**
 * Class for a table that shows enrichment of TF targets along a path
 */
public class DrawMainHeatMap extends JPanel implements ActionListener {
	JFrame theframe;
	private JFrame saveImageFrame;
	JPanel tfPanel;
	TSMiner_Timeiohmm theTimeiohmm;
	TSMiner_Timeiohmm.Treenode treecopy;
	
	JButton TFbutton;
	JButton saveimage;
	JButton savemodel;
	JButton repaint;
	HeatMap2 tfheatmap;
	int npermutationval;
	JFrame saveModelFrame;

	final static Color bgColor = Color.white;
	final static Color fgColor = Color.black;
	int numrows;
	int colnum;
	NumberFormat nf;
	NumberFormat nf2;
	boolean bsplit;
	DecimalFormat df1;
	DecimalFormat df2;
	
	/**
	 * Constructor - builds the table
	 */
	public DrawMainHeatMap(TSMiner_Timeiohmm theTimeiohmm, Treenode treecopy, int npermutationval, 
			String[] dsamplemins, String[] genename, double[][] data, int[] heatmap2gene) {
		
		//normtf = Util.Hierarchical(normtf);
		//Util.BubbleSort_dec(normtf, normtf.length, ndepth);
		//Hierarchical_clustering.Hcluster_clustering(expression, genename);
		this.theTimeiohmm = theTimeiohmm;
		this.treecopy = treecopy;
		this.npermutationval = npermutationval;
		
		setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
		setBackground(bgColor);
		setForeground(fgColor);
		
		if(data != null && data.length>0){
			
			try {
				double height = genename.length;
				double width = dsamplemins.length*70;
				tfPanel = new JPanel();
				tfPanel.setBorder(new TitledBorder(null,"Tree structure of " + (int)height +" genes",TitledBorder.LEFT,TitledBorder.TOP));
				BoxLayout layout = new BoxLayout(tfPanel, BoxLayout.Y_AXIS);
				tfPanel.setLayout(layout);
				tfPanel.setPreferredSize(new Dimension(Math.min(1000,Math.max(400, (int)width)), Math.min(600, Math.max(500, (int)height))));

				tfheatmap = getheatmap(theTimeiohmm, treecopy, heatmap2gene, 
						data, dsamplemins, genename, npermutationval);
				tfPanel.add(tfheatmap);
				add(tfPanel);
				
			} catch (Exception e1) {
				e1.printStackTrace();
			}
			
			addBottom();
		}
	}

	private void addBottom() {
		TFbutton = new JButton("TF Summary");
		TFbutton.setActionCommand("tftable");
		TFbutton.setMinimumSize(new Dimension(800, 20));
		TFbutton.addActionListener(this);

		saveimage = new JButton("Save Image", Util.createImageIcon("Save16.gif"));
		saveimage.setActionCommand("saveimage");
		saveimage.setMinimumSize(new Dimension(800, 20));
		saveimage.addActionListener(this);
		
		savemodel = new JButton("Save Model", Util.createImageIcon("Save16.gif"));
		savemodel.setActionCommand("savemodel");
		savemodel.setMinimumSize(new Dimension(800, 20));
		savemodel.addActionListener(this);
		
		repaint = new JButton("Repaint");
		repaint.setActionCommand("repaint");
		repaint.setMinimumSize(new Dimension(800, 20));
		repaint.addActionListener(this);
		
		JPanel buttonPanel = new JPanel();
		buttonPanel.add(TFbutton);
		buttonPanel.add(saveimage);
		buttonPanel.add(savemodel);
		buttonPanel.add(repaint);
		
		buttonPanel.setMaximumSize(new Dimension(Integer.MAX_VALUE, 20));
		add(buttonPanel);
	}

	public HeatMap2 getheatmap(TSMiner_Timeiohmm theTimeiohmm, Treenode treecopy, 
			int[] heatmap2gene, double[][] data, String[] dsamplemins, 
			String[] genelist, int npermutationval) throws Exception{
        //Gradient gra = new Gradient();
        data = transpose(data);
        boolean useGraphicsYAxis = true;
    /*
          天蓝(135,206,235);
          青(0,255,255);
          品红(255,0,135);
    */
        //Color[] gradientColors = new Color[]{Color.BLUE, Color.WHITE, Color.RED};
        //Color[] customGradient = Gradient.createMultiGradient(gradientColors, 40);
        Color[] gradientColors = new Color[]{Color.BLUE, Color.YELLOW};
        Color[] customGradient = Gradient.createMultiGradient(gradientColors, 40);
        HeatMap2 panel = new HeatMap2(theframe, theTimeiohmm, treecopy, heatmap2gene, 
        		data, dsamplemins, genelist, npermutationval, 
        		useGraphicsYAxis, customGradient);
        
        // set miscelaneous settings
        panel.setDrawLegend(true);

        panel.setTitle("Tree Structure");
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
	
	/**
	 * Responds to buttons being pressed on the interface
	 */
	public void actionPerformed(ActionEvent e) {
		
		String szCommand = e.getActionCommand();
		if(szCommand.equals("tftable")){
			javax.swing.SwingUtilities.invokeLater(new Runnable() {
				public void run() {
					JFrame frame1 = new JFrame();
					JDialog frame = new JDialog(frame1, "TF Table", true);
					frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
					frame.setLocation(200, 200);
					Container theDialogContainer = frame.getContentPane();
					theDialogContainer.setBackground(Color.white);
					JTabbedPane tabbedPane = new JTabbedPane();
					
					if(theTimeiohmm.activeTF != null){
						for(int i=0;i<theTimeiohmm.theDataSet.dsamplemins.length;i++){
							Iterator iter1 = theTimeiohmm.activeTF.entrySet().iterator();
					    	while (iter1.hasNext()) {
					    	Map.Entry entry = (Map.Entry) iter1.next();
					    	Object key = entry.getKey();
					    	List<String> tfs = theTimeiohmm.activeTF.get(key);
					    	List<String> depvalues = theTimeiohmm.activeDEpvalue.get(key);
					    	if(tfs != null && tfs.size()>0 && key.equals(theTimeiohmm.theDataSet.dsamplemins[i])){
					    		List<String[]> table1 = new ArrayList<String[]>();
					    		List<String[]> table2 = new ArrayList<String[]>();
					    		for(int p=0;p<tfs.size();p++){
					    			String[] table = new String[4];
					    			table[0] = tfs.get(p);
					    			table[1] = Arrays.asList(depvalues.get(p).split("\\t+")).get(0);
					    			table[2] = Arrays.asList(depvalues.get(p).split("\\t+")).get(1);
					    			table[3] = Arrays.asList(depvalues.get(p).split("\\t+")).get(2);
					    			if(table[0].contains("posi")){
					    				table1.add(table);
					    			}else{
					    				table2.add(table);
					    			}
					    		}
					    		TSMinerGui_TFTable newContentPane1 = new TSMinerGui_TFTable(
										treecopy, i, frame, theTimeiohmm, table1, 
										table2, (String)key, npermutationval);
								newContentPane1.setOpaque(true); // content panes must be opaque
								tabbedPane.addTab(key + " promotion", null, newContentPane1,key + " promotion");
					    	}
					    	}
					    	
					    	Iterator iter2 = theTimeiohmm.repressTF.entrySet().iterator();
					    	while (iter2.hasNext()) {
					    	Map.Entry entry = (Map.Entry) iter2.next();
					    	Object key = entry.getKey();
					    	List<String> tfs = theTimeiohmm.repressTF.get(key);
					    	List<String> depvalues = theTimeiohmm.repressDEpvalue.get(key);
					    	if(tfs != null && tfs.size()>0 && key.equals(theTimeiohmm.theDataSet.dsamplemins[i])){
					    		List<String[]> table1 = new ArrayList<String[]>();
					    		List<String[]> table2 = new ArrayList<String[]>();
					    		for(int p=0;p<tfs.size();p++){
					    			String[] table = new String[4];
					    			table[0] = tfs.get(p);
					    			table[1] = Arrays.asList(depvalues.get(p).split("\\t+")).get(0);
					    			table[2] = Arrays.asList(depvalues.get(p).split("\\t+")).get(1);
					    			table[3] = Arrays.asList(depvalues.get(p).split("\\t+")).get(2);
					    			if(table[0].contains("posi")){
					    				table1.add(table);
					    			}else{
					    				table2.add(table);
					    			}
					    		}
					    		TSMinerGui_TFTable newContentPane1 = new TSMinerGui_TFTable(
										treecopy, i, frame,theTimeiohmm, table1, 
										table2, (String)key, npermutationval);
								newContentPane1.setOpaque(true); // content panes must be opaque
								tabbedPane.addTab(key + " repression ", null, newContentPane1,key + " repression");
					    	}
					    	}
						}								
					}
					theDialogContainer.add(tabbedPane);
					frame1.setContentPane(theDialogContainer);
					frame1.pack();
					frame1.setVisible(true);	
				}
			});
		}else if(szCommand.equals("savemodel")){
			javax.swing.SwingUtilities.invokeLater(new Runnable() {
				public void run() {
					if (saveModelFrame == null) {
						saveModelFrame = new JFrame("Save Model to File");
						saveModelFrame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
						saveModelFrame.setLocation(400,300);
						TSMinerGui_SaveModel newContentPane = new TSMinerGui_SaveModel(
								theTimeiohmm,
								treecopy,
								saveModelFrame);
						newContentPane.setOpaque(true);
						// content panes must be opaque
						saveModelFrame.setContentPane(newContentPane);
						// Display the window.
						saveModelFrame.pack();
					} else {
						saveModelFrame.setExtendedState(Frame.NORMAL);
					}
					saveModelFrame.setVisible(true);
				}
			});
		}else if(szCommand.equals("saveimage")){
			javax.swing.SwingUtilities.invokeLater(new Runnable() {
				public void run() {
					if (saveImageFrame == null) {
						saveImageFrame = new JFrame("Save as Image");
						saveImageFrame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
						saveImageFrame.setLocation(400,300);
						TSMinerGui_SaveImage newContentPane = new TSMinerGui_SaveImage(saveImageFrame, tfheatmap);
						newContentPane.setOpaque(true);
						saveImageFrame.setContentPane(newContentPane);
						saveImageFrame.pack();
					} else {
						saveImageFrame.setExtendedState(Frame.NORMAL);
					}
					saveImageFrame.setVisible(true);
				}
			});
		}else if(szCommand.equals("repaint")){
			tfheatmap.zoom = 1;
			tfheatmap.repaint();
		}
	}

	class Dialog implements Runnable{
	    public void run() {
	    	try {
	    		String text1 = "The expression levels of TFs cannot be found in the time-series data";
		        JOptionPane.showMessageDialog(theframe,text1);
	    	}catch (Exception e) {
				e.printStackTrace();
		}
	    }
	}
	
	/**
	 * Converts the value of dval to a String that is displayed on the table
	 */
	public static String doubleToSz(double dval) {
		String szexp;
		double dtempval = dval;
		int nexp = 0;

		NumberFormat nf2 = NumberFormat.getInstance(Locale.ENGLISH);
		nf2.setMinimumFractionDigits(3);
		nf2.setMaximumFractionDigits(3);

		NumberFormat nf1 = NumberFormat.getInstance(Locale.ENGLISH);
		nf1.setMinimumFractionDigits(2);
		nf1.setMaximumFractionDigits(2);

		if (dval <= 0) {
			szexp = "0.000";
		} else {
			while ((dtempval < 0.9995) && (dtempval > 0)) {
				nexp--;
				dtempval = dtempval * 10;
			}

			if (nexp < -2) {
				dtempval = Math.pow(10, Math.log(dval) / Math.log(10) - nexp);
				szexp = nf1.format(dtempval) + "e" + nexp;

			} else {
				szexp = nf2.format(dval);
			}
		}

		return szexp;
	}
	
	class Pathway_pval {
		String pathwayID;
		double pvalue;
		Pathway_pval(String pathwayID, double pvalue) {
			this.pathwayID = pathwayID;
			this.pvalue = pvalue;
		}
	}


    
    public String[][] BubbleSort_increase(String[][] r, Integer n, Integer col) //升序冒泡排序
 	{
 		 int low = 0;   
 		    int high= n - 1; //设置变量的初始值  
 		    String[] tmp;
 		    int j;  
 		    while (low < high) {  
 		        for (j= low; j< high; ++j) //正向冒泡,找到最大者  
 		            if (Double.parseDouble(r[j][col])> Double.parseDouble(r[j+1][col])) {  
 		                tmp = r[j]; r[j]=r[j+1];r[j+1]=tmp;  
 		            }   
 		        --high;                 //修改high值, 前移一位  
 		        for ( j=high; j>low; --j) //反向冒泡,找到最小者  
 		            if (Double.parseDouble(r[j][col])<Double.parseDouble(r[j-1][col])) {  
 		                tmp = r[j]; r[j]=r[j-1];r[j-1]=tmp; 
 		            }  
 		        ++low; //修改low值,后移一位  
 		    }   
 	return r;
 	}
    
    public String[][] BubbleSort_decrease(String[][] r, Integer n, Integer col) //升序冒泡排序
 	{
 		 int low = 0;   
 		    int high= n - 1; //设置变量的初始值  
 		    String[] tmp;
 		    int j;  
 		    while (low < high) {  
 		        for (j= low; j< high; ++j) //正向冒泡,找到最大者  
 		            if (Double.parseDouble(r[j][col])< Double.parseDouble(r[j+1][col])) {  
 		                tmp = r[j]; r[j]=r[j+1];r[j+1]=tmp;  
 		            }   
 		        --high;                 //修改high值, 前移一位  
 		        for ( j=high; j>low; --j) //反向冒泡,找到最小者  
 		            if (Double.parseDouble(r[j][col])>Double.parseDouble(r[j-1][col])) {  
 		                tmp = r[j]; r[j]=r[j-1];r[j-1]=tmp;  
 		            }  
 		        ++low; //修改low值,后移一位  
 		    }   
 	return r;
 	}

    public Integer getrank(List<Double> pvalue, double value){
    	int cc = 1;
    	for(int i=0;i<pvalue.size();i++){
    		if(pvalue.get(i)<value){
    			cc++;
    		}
    	}
    	return cc;
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