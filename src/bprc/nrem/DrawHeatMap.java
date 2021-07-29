package bprc.nrem;

import heatmapframe.GetHeatMap;
import heatmapframe.HeatMap;

import javax.swing.BoxLayout;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextArea;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Frame;
import java.util.*;
import java.util.List;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.awt.*;
import java.awt.event.*;

import javax.swing.*;
import javax.swing.border.TitledBorder;

import clustering.Hierarchical_clustering;
import bprc.core.Util;


/**
 * Class for a table that shows enrichment of TF targets along a path
 */
public class DrawHeatMap extends JPanel implements ActionListener {
	JFrame theframe;
	Component parentComponent;
	private JFrame saveImageFrame;
	private JFrame saveImageFrame2;
	private JFrame saveImageFrame3;
	JPanel tfPanel;
	JPanel tgPanel;
	JPanel pgPanel;
	String[] dsamplemins;
	String[] tfname;
	String[] tgname;
	String[] pgname;
	
	JButton saveTF;
	JButton saveTG;
	JButton savePG;

	double[][] normtf;
	double[][] normtg;
	double[][] normpg;
	
	HeatMap tfheatmap;
	HeatMap tgheatmap;
	HeatMap pgheatmap;
	List<String> misstf;
	

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
	public DrawHeatMap(JFrame frame1, Component parentComponent, String[] dsamplemins, 
			List<String> tflist, List<String> tglist, List<String> pglist, List<String> misstf, int ndepth, 
			List<double[]> tfexpression, List<double[]> tgexpression, List<double[]> pgexpression) {
		// assuming that if called with root node then only one child so stats
		this.theframe = frame1;
		this.parentComponent = parentComponent;
		this.dsamplemins = dsamplemins;
		this.misstf = misstf;
		this.tfname = new String[tflist.size()];
		for(int i=0;i<tflist.size();i++) tfname[i] = tflist.get(i);
		this.tgname = new String[tglist.size()];
		for(int i=0;i<tglist.size();i++) tgname[i] = tglist.get(i);
		this.pgname = new String[pglist.size()];
		for(int i=0;i<pglist.size();i++) pgname[i] = pglist.get(i);
		
		if(tfexpression != null && tfexpression.size()>0){
			normtf = new double[tfexpression.size()][tfexpression.get(0).length];
			for(int i=0;i<tfexpression.size();i++){
				double[] expression = tfexpression.get(i);
				double largest = Double.MIN_VALUE;
		        double smallest = Double.MAX_VALUE;
		        for(int j=0;j<expression.length;j++) {
					largest = Math.max(expression[j], largest);
					smallest = Math.min(expression[j], smallest);
				}
		        double range = largest - smallest;
		        for(int j=0;j<expression.length;j++){
		        	normtf[i][j] = (expression[j] - smallest) / range;
				}
			}
			//normtf = Util.Hierarchical(normtf);
			//Util.BubbleSort_dec(normtf, normtf.length, ndepth);
			Hierarchical_clustering.Hcluster_clustering(normtf, tfname);
		}
		
		if(tgexpression != null && tgexpression.size()>0){
			normtg = new double[tgexpression.size()][tgexpression.get(0).length];
			for(int i=0;i<tgexpression.size();i++){
				double[] expression = tgexpression.get(i);
				double largest = Double.MIN_VALUE;
		        double smallest = Double.MAX_VALUE;
		        for(int j=0;j<expression.length;j++) {
					largest = Math.max(expression[j], largest);
					smallest = Math.min(expression[j], smallest);
				}
		        double range = largest - smallest;
		        for(int j=0;j<expression.length;j++){
		        	normtg[i][j] = (expression[j] - smallest) / range;
				}
			}
			Hierarchical_clustering.Hcluster_clustering(normtg, tgname);
		}
		
		
		if(pgexpression != null && pgexpression.size()>0){
			normpg = new double[pgexpression.size()][pgexpression.get(0).length];
			for(int i=0;i<pgexpression.size();i++){
				double[] expression = pgexpression.get(i);
				double largest = Double.MIN_VALUE;
		        double smallest = Double.MAX_VALUE;
		        for(int j=0;j<expression.length;j++) {
					largest = Math.max(expression[j], largest);
					smallest = Math.min(expression[j], smallest);
				}
		        double range = largest - smallest;
		        for(int j=0;j<expression.length;j++){
		        	normpg[i][j] = (expression[j] - smallest) / range;
				}
			}
			Hierarchical_clustering.Hcluster_clustering(normpg, pgname);
		}
		
		setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
		setBackground(bgColor);
		setForeground(fgColor);
		
		if(normtf != null && normtf.length>0){
			addtfmap();
		}
		if(normtg != null && normtg.length>0){
			addtgmap();
		}
		if(normpg != null && normpg.length>0){
			addpgmap();
		}
		addBottom();
		
		if(misstf != null && misstf.size()>0){
			//addTop2();
			Dialog dialog = new Dialog();
			Thread t1 = new Thread(dialog);
			t1.start();
		}
		
	}

	
	
    private void addtfmap() {
		try {
			
			double numtf = 0;
			if(normtf!=null){
				numtf = normtf.length;
			}
			double numtg = 0;
			if(normtg!=null){
				numtg = normtg.length;
			}
			double numpg = 0;
			if(normpg!=null){
				numpg = normpg.length;
			}
			double	height = 600*(numtf/(numtf+numtg+numpg));
			
			tfPanel = new JPanel();
			tfPanel.setBorder(new TitledBorder(null,"Expression of TFs",TitledBorder.LEFT,TitledBorder.TOP));
			BoxLayout layout = new BoxLayout(tfPanel, BoxLayout.Y_AXIS);
			tfPanel.setLayout(layout);
			tfPanel.setPreferredSize(new Dimension(this.getWidth(), Math.min(300, Math.max(150, (int)height))));

			GetHeatMap heatmap = new GetHeatMap();
			tfheatmap = heatmap.heatmap(normtf, dsamplemins, tfname);
			tfPanel.add(tfheatmap);
			add(tfPanel);
			
		} catch (Exception e1) {
			e1.printStackTrace();
		}
	}
    private void addtgmap() {
		
		try {
			double numtf = 0;
			if(normtf!=null){
				numtf = normtf.length;
			}
			double numtg = 0;
			if(normtg!=null){
				numtg = normtg.length;
			}
			double numpg = 0;
			if(normpg!=null){
				numpg = normpg.length;
			}
			double height = 600*(numtg/(numtf+numtg+numpg));
			tgPanel = new JPanel();
			tgPanel.setBorder(new TitledBorder(null,"Expression of TGs",TitledBorder.LEFT,TitledBorder.TOP));
			BoxLayout layout = new BoxLayout(tgPanel, BoxLayout.Y_AXIS);
			tgPanel.setLayout(layout);
			tgPanel.setPreferredSize(new Dimension(this.getWidth(), Math.min(300, Math.max(150, (int)height))));
			
			GetHeatMap heatmap = new GetHeatMap();
			tgheatmap = heatmap.heatmap(normtg, dsamplemins, tgname);
			tgPanel.add(tgheatmap);
			add(tgPanel);
		} catch (Exception e1) {
			e1.printStackTrace();
		}
	}
    private void addpgmap() {
		
		try {
			double numtf = 0;
			if(normtf!=null){
				numtf = normtf.length;
			}
			double numtg = 0;
			if(normtg!=null){
				numtg = normtg.length;
			}
			double numpg = 0;
			if(normpg!=null){
				numpg = normpg.length;
			}
			
			double height = 600*(numpg/(numtf+numtg+numpg));
			pgPanel = new JPanel();
			pgPanel.setBorder(new TitledBorder(null,"Expression of PGs",TitledBorder.LEFT,TitledBorder.TOP));
			BoxLayout layout = new BoxLayout(pgPanel, BoxLayout.Y_AXIS);
			pgPanel.setLayout(layout);
			pgPanel.setPreferredSize(new Dimension(this.getWidth(), Math.min(300, Math.max(200, (int)height))));
			
			GetHeatMap heatmap = new GetHeatMap();
			pgheatmap = heatmap.heatmap(normpg, dsamplemins, pgname);
			pgPanel.add(pgheatmap);
			add(pgPanel);
		} catch (Exception e1) {
			e1.printStackTrace();
		}
	}

    private void addTop() {
    	JPanel topPanel = new JPanel();
		JLabel Label1 = new JLabel();
		Label1.setText("   The expression levels of " + misstf.size()+" TFs cannot be found in the time-series data:   ");
		JLabel Label2 = new JLabel();
		Label2.setText(misstf+"");
		
		BoxLayout boxLayout = new BoxLayout(topPanel, BoxLayout.Y_AXIS);
		topPanel.setLayout(boxLayout);
		Label1.setAlignmentX(CENTER_ALIGNMENT);
		Label2.setAlignmentX(CENTER_ALIGNMENT);
		
		topPanel.add(Label1);
		topPanel.add(Label2);
		topPanel.setMaximumSize(new Dimension(Integer.MAX_VALUE, 20));
		//topPanel.setBackground(Color.white);
		topPanel.setBorder(BorderFactory.createEtchedBorder());
		add(topPanel);
    }
    
    private void addTop2() {
    	   String text = "The expression of " + misstf.size()+" TFs cannot be found in the time-series data:\n"
    			   +misstf;
    	   JTextArea textArea = new JTextArea(3,20);
    	    textArea.setText(text);
    	    textArea.setWrapStyleWord(true);
    	    textArea.setLineWrap(true);
    	    textArea.setOpaque(false);
    	    textArea.setEditable(false);
    	    textArea.setFocusable(false);
    	    textArea.setBackground(UIManager.getColor("Label.background"));
    	    textArea.setFont(UIManager.getFont("Label.font"));
    	    textArea.setBorder(UIManager.getBorder("Label.border"));
    	    textArea.setPreferredSize(new Dimension(400, 20));
    	    
    	    JPanel topPanel = new JPanel();
    	    BoxLayout boxLayout = new BoxLayout(topPanel, BoxLayout.Y_AXIS);
    		topPanel.setLayout(boxLayout);
    		textArea.setAlignmentX(CENTER_ALIGNMENT);
    		
    	    topPanel.add(textArea);
    	    topPanel.setMaximumSize(new Dimension(Integer.MAX_VALUE, 20));
    		//topPanel.setBackground(Color.white);
    		topPanel.setBorder(BorderFactory.createEtchedBorder());
    		add(topPanel);
    	   
    }
    
	private void addBottom() {
		saveTF = new JButton("Save TF Map", Util.createImageIcon("Save16.gif"));
		saveTF.setActionCommand("heattf");
		saveTF.setMinimumSize(new Dimension(800, 20));
		saveTF.addActionListener(this);

		saveTG = new JButton("Save TG Map", Util.createImageIcon("Save16.gif"));
		saveTG.setActionCommand("heattg");
		saveTG.setMinimumSize(new Dimension(800, 20));
		saveTG.addActionListener(this);
		
		savePG = new JButton("Save PG Map", Util.createImageIcon("Save16.gif"));
		savePG.setActionCommand("heatpg");
		savePG.setMinimumSize(new Dimension(800, 20));
		savePG.addActionListener(this);
		
		JPanel buttonPanel = new JPanel();
		if(normtf != null && normtf.length>0){
			buttonPanel.add(saveTF);
		}
		if(normtg!= null && normtg.length>0){
			buttonPanel.add(saveTG);
		}
		if(normpg!= null && normpg.length>0){
			buttonPanel.add(savePG);
		}
		buttonPanel.setMaximumSize(new Dimension(Integer.MAX_VALUE, 20));
		add(buttonPanel);
	}

	/**
	 * Responds to buttons being pressed on the interface
	 */
	public void actionPerformed(ActionEvent e) {
		
		String szCommand = e.getActionCommand();
		if(szCommand.equals("heattf")){
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
		} else if(szCommand.equals("heattg")){
			javax.swing.SwingUtilities.invokeLater(new Runnable() {
				public void run() {
					if (saveImageFrame2 == null) {
						saveImageFrame2 = new JFrame("Save as Image");
						saveImageFrame2.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
						saveImageFrame2.setLocation(400,300);
						TSMinerGui_SaveImage newContentPane = new TSMinerGui_SaveImage(saveImageFrame2, tgheatmap);
						newContentPane.setOpaque(true);
						saveImageFrame2.setContentPane(newContentPane);
						saveImageFrame2.pack();
					} else {
						saveImageFrame2.setExtendedState(Frame.NORMAL);
					}
					saveImageFrame2.setVisible(true);
				}
			});
			
		} else if(szCommand.equals("heatpg")){
			javax.swing.SwingUtilities.invokeLater(new Runnable() {
				public void run() {
					if (saveImageFrame3 == null) {
						saveImageFrame3 = new JFrame("Save as Image");
						saveImageFrame3.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
						saveImageFrame3.setLocation(400,300);
						TSMinerGui_SaveImage newContentPane = new TSMinerGui_SaveImage(saveImageFrame3, pgheatmap);
						newContentPane.setOpaque(true);
						saveImageFrame3.setContentPane(newContentPane);
						saveImageFrame3.pack();
					} else {
						saveImageFrame3.setExtendedState(Frame.NORMAL);
					}
					saveImageFrame3.setVisible(true);
				}
			});
			
		}
	}

	class Dialog implements Runnable{
	    public void run() {
	    	try {
	    		String text1 = "The expression levels of " + misstf.size()+
	    				" TFs cannot be found in the time-series data:\n" + misstf;
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
    
    private void rank(double[][] normcontrol, double[][] normcase){
    	int n = normcontrol.length;
    	int col = 0;
		
		int low = 0;   
	    int high= n -1; //设置变量的初始值  
	    double[] tmp;
	    int j;  
	    while (low < high) {
	        for (j= low; j< high; ++j) //正向冒泡,找到最大者  
	            if (normcontrol[j][col]< normcontrol[j+1][col]) {  
	                tmp = normcontrol[j]; normcontrol[j]=normcontrol[j+1];normcontrol[j+1]=tmp;
	                tmp = normcase[j]; normcase[j]=normcase[j+1];normcase[j+1]=tmp;
	            }
	        --high;//修改high值, 前移一位  
	        for ( j=high; j>low; --j) //反向冒泡,找到最小者  
	            if (normcontrol[j][col]>normcontrol[j-1][col]) {  
	                tmp = normcontrol[j]; normcontrol[j]=normcontrol[j-1];normcontrol[j-1]=tmp;
	                tmp = normcase[j]; normcase[j]=normcase[j+1];normcase[j+1]=tmp;
	            }  
	        ++low; //修改low值,后移一位  
	    }
    }
}