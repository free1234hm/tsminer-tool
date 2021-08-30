package bprc.nrem;

import bprc.core.*;
import bprc.nrem.TSMiner_Timeiohmm.Treenode;

import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.JTable;

import java.awt.Color;
import java.awt.Component;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.event.ComponentAdapter;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.*;
import java.util.List;
import java.util.Map.Entry;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.awt.datatransfer.*;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.table.*;

import java.math.*;

/**
 * Class for a table that shows enrichment of TF targets along a path
 */
public class TSMinerGui_PathwayTable extends JPanel implements ActionListener {
	private BufferedWriter outXml;
	String[][] ResultList;
	TSMiner_DataSet rawDataSet;
	int ndepth;
	HashMap<String, List<String>> pathwaygene;
	HashMap<String, List<List<String>>> tf_tg;
	HashMap<String, List<List<String>>> tf_pg;
	JDialog dialog;
	JFrame theframe;
	String[] columnNames;
	String[][] tabledata;
	JButton copyButton;
	JButton saveButton;
	JButton mapButton;
	JButton heatmapButton;
	JScrollPane scrollPane;
	TableSorter sorter;
	JTable table;
	final static Color bgColor = Color.white;
	final static Color fgColor = Color.black;
	TSMinerGui_FilterStaticModel hmst;
	int numrows;
	int colnum;
	int temp;
	NumberFormat nf;
	NumberFormat nf2;
	NumberFormat df1;
	NumberFormat df2;
	
	/**
	 * Constructor - builds the table
	 */
	public TSMinerGui_PathwayTable(String[][] ResultList, JDialog dialog, JFrame frame1, 
			HashMap<String, List<String>> pathwaygene, HashMap<String, List<List<String>>> tf_tg, 
			HashMap<String, List<List<String>>> tf_pg, int ndepth, String time, TSMiner_DataSet rawDataSet, 
			int temp) {
		// assuming that if called with root node then only one child so stats
		// for root will be same as first child
		
		this.dialog = dialog;
		this.theframe = frame1;
		this.ResultList = ResultList;
		this.pathwaygene = pathwaygene;
		this.tf_tg = tf_tg;
		this.tf_pg = tf_pg;
		this.rawDataSet = rawDataSet;
		this.ndepth = ndepth;
		this.temp = temp;
		df1 = NumberFormat.getInstance(Locale.ENGLISH);
		df1.setMinimumFractionDigits(2);
		df1.setMaximumFractionDigits(2);
		df2 = NumberFormat.getInstance(Locale.ENGLISH);
		df2.setMinimumFractionDigits(4);
		df2.setMaximumFractionDigits(4);
		
		setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
		setBackground(bgColor);
		setForeground(fgColor);

		numrows = ResultList.length;
		colnum = ResultList[0].length;
		nf2 = NumberFormat.getInstance(Locale.ENGLISH);
		nf2.setMinimumFractionDigits(3);
		nf2.setMaximumFractionDigits(3);
		
		if(temp == 0){
			columnNames = new String[colnum];
			tabledata = new String[numrows][colnum];
			columnNames[0] = "Pathway";
			columnNames[1] = "Num Total";
			columnNames[2] = "Percent";
			columnNames[3] = time+" DE p-value";
			columnNames[4] = time+" DE q-value";
			columnNames[5] = time+" Avg. FC";
			
			for (int nrow = 0; nrow < numrows; nrow++) {
				for(int ncol = 0; ncol<colnum; ncol++){
					tabledata[nrow][ncol] = ResultList[nrow][ncol];
				}
			}
			
			//Color up = new Color(255, 0, 0, 120);
			//Color down = new Color(0, 0, 255, 120);

			sorter = new TableSorter(new TableModelST(tabledata, columnNames));
			table = new JTable(sorter);
			sorter.setTableHeader(table.getTableHeader());
			
			table.setPreferredScrollableViewportSize(new Dimension(800, Math.min((table.getRowHeight() + table.getRowMargin())
							* table.getRowCount(), 400)));

			table.setAutoResizeMode(JTable.AUTO_RESIZE_LAST_COLUMN);
			TableColumn column;
			column = table.getColumnModel().getColumn(0);
			column.setPreferredWidth(200);
			for (int ncolindex = 1; ncolindex < columnNames.length; ncolindex++) {
				column = table.getColumnModel().getColumn(ncolindex);
				column.setPreferredWidth(120);
			}
			column.setMaxWidth(200);
			setColumnColor(table);
			
		}else{
			columnNames = new String[colnum];
			tabledata = new String[numrows][colnum];
			columnNames[0] = "Pathway";
			columnNames[1] = "Num TG";
			columnNames[2] = "Num PG";
			columnNames[3] = "Num Overlap";
			columnNames[4] = "Enrichment Pvalue";
			columnNames[5] = time+" DE p-value";
			columnNames[6] = time+" DE q-value";
			columnNames[7] = time+" Avg. FC";
			
			for (int nrow = 0; nrow < numrows; nrow++) {
				for(int ncol = 0; ncol<colnum; ncol++){
					tabledata[nrow][ncol] = ResultList[nrow][ncol];
				}
			}
			
			//Color up = new Color(255, 0, 0, 120);
			//Color down = new Color(0, 0, 255, 120);

			sorter = new TableSorter(new TableModelST(tabledata, columnNames));
			table = new JTable(sorter);
			sorter.setTableHeader(table.getTableHeader());
			table.setPreferredScrollableViewportSize(new Dimension(700, Math.min((table.getRowHeight() + table.getRowMargin())
							* table.getRowCount(), 400)));

			table.setAutoResizeMode(JTable.AUTO_RESIZE_ALL_COLUMNS);
			TableColumn column;
			
			column = table.getColumnModel().getColumn(0);
			column.setPreferredWidth(200);
			
			for (int ncolindex = 1; ncolindex < columnNames.length; ncolindex++) {
				column = table.getColumnModel().getColumn(ncolindex);
				column.setPreferredWidth(120);
			}
			column.setMaxWidth(200);
			setColumnColor(table);
		}

		scrollPane = new JScrollPane(table);
		add(scrollPane);

		addBottom();
		
		this.addComponentListener(new ComponentAdapter(){
			public void componentResized(ComponentEvent e){
			if(theframe.getWidth()>(620+(columnNames.length-6)*120)){
				TableColumn column = table.getColumnModel().getColumn(columnNames.length-1);
				column.setPreferredWidth(theframe.getWidth()-470-(columnNames.length-6)*120);
			}
			}
			});
	}

	/**
	 * Helper function that adds information displayed at the bottom of the
	 * table information window
	 */
	private void addBottom() {

		nf = NumberFormat.getInstance(Locale.ENGLISH);
		nf.setMinimumFractionDigits(3);
		nf.setMaximumFractionDigits(3);

		heatmapButton = new JButton("Show heatmap");
		heatmapButton.setActionCommand("heatmap");
		heatmapButton.setMinimumSize(new Dimension(800, 20));
		heatmapButton.addActionListener(this);
		
		mapButton = new JButton("Show network");
		mapButton.setActionCommand("network");
		mapButton.setMinimumSize(new Dimension(800, 20));
		mapButton.addActionListener(this);
		
		copyButton = new JButton("Copy Table", Util.createImageIcon("Copy16.gif"));
		copyButton.setActionCommand("copy");
		copyButton.setMinimumSize(new Dimension(800, 20));
		copyButton.addActionListener(this);

		saveButton = new JButton("Save Table", Util.createImageIcon("Save16.gif"));
		saveButton.setActionCommand("save");
		saveButton.setMinimumSize(new Dimension(800, 20));
		saveButton.addActionListener(this);

		JPanel buttonPanel = new JPanel();
		buttonPanel.setBackground(Color.white);
		
		buttonPanel.add(heatmapButton);
		buttonPanel.add(mapButton);
		buttonPanel.add(copyButton);
		buttonPanel.add(saveButton);

		buttonPanel.setMaximumSize(new Dimension(Integer.MAX_VALUE, 20));
		add(buttonPanel);	
	}

	/**
	 * Writes the content of the table to a file specified through pw
	 */
	public void printFile(PrintWriter pw) {
		for (int ncol = 0; ncol < columnNames.length - 1; ncol++) {
			pw.print(columnNames[ncol] + "\t");
		}
		pw.println(columnNames[columnNames.length - 1]);

		for (int nrow = 0; nrow < tabledata.length; nrow++) {
			for (int ncol = 0; ncol < tabledata[nrow].length - 1; ncol++) {
				pw.print(sorter.getValueAt(nrow, ncol) + "\t");
			}
			pw.println(sorter.getValueAt(nrow, columnNames.length - 1));
		}
	}

	/**
	 * Copies the content of the table to the clipboard
	 */
	public void writeToClipboard() {
		StringBuffer sbuf = new StringBuffer();
		for (int ncol = 0; ncol < columnNames.length - 1; ncol++) {
			sbuf.append(columnNames[ncol] + "\t");
		}
		sbuf.append(columnNames[columnNames.length - 1] + "\n");

		for (int nrow = 0; nrow < tabledata.length; nrow++) {
			for (int ncol = 0; ncol < tabledata[nrow].length - 1; ncol++) {
				sbuf.append(sorter.getValueAt(nrow, ncol) + "\t");
			}
			sbuf.append(sorter.getValueAt(nrow, columnNames.length - 1) + "\n");
		}
		// get the system clipboard
		Clipboard systemClipboard = Toolkit.getDefaultToolkit()
				.getSystemClipboard();
		// set the textual content on the clipboard to our
		// Transferable object

		Transferable transferableText = new StringSelection(sbuf.toString());
		systemClipboard.setContents(transferableText, null);
	}

	/**
	 * Responds to buttons being pressed on the interface
	**/
	public void actionPerformed(ActionEvent e) {
		String szCommand = e.getActionCommand();
		if (szCommand.equals("heatmap")) {
			int[] rows = table.getSelectedRows();
			List<String> pathways = new ArrayList<String>();
			for(int i=0;i<rows.length;i++){
				pathways.add(table.getValueAt(rows[i],0).toString());
			}
			if(pathways.size()>0){
				if(temp == 0){
					Set<String> TF = new HashSet<String>();
					Set<String> TG = new HashSet<String>();
					Set<String> PG = new HashSet<String>();
					
					for(int p=0;p<pathways.size();p++){
						List<List<String>> gs1 = tf_tg.get(pathways.get(p));
						if(gs1 != null && gs1.size()>0){
							for(int i=0;i<gs1.size();i++){
								List<String> a = gs1.get(i);
								TF.add(a.get(0));
								TG.add(a.get(1));
						}
						}
						List<List<String>> gs2 = tf_pg.get(pathways.get(p));
						if(gs2 != null && gs2.size()>0){
							for(int i=0;i<gs2.size();i++){
								List<String> a = gs2.get(i);	
								TF.add(a.get(0));
								PG.add(a.get(1));
							}
						}
					}
					
					if(TF!=null && TF.size()>0){
						javax.swing.SwingUtilities.invokeLater(new Runnable() {
							public void run() {
								List<double[]> tfexpression = new ArrayList<double[]>();
								List<double[]> tgexpression = new ArrayList<double[]>();
								List<double[]> pgexpression = new ArrayList<double[]>();
								List<String> tflist = new ArrayList<String>();
								List<String> tglist = new ArrayList<String>();
								List<String> pglist = new ArrayList<String>();
								List<String> misstf = new ArrayList<String>();
								boolean found = false;
								
								for(String tf:TF){
									for(int j=0;j<rawDataSet.genenames.length;j++){
										if(tf.equalsIgnoreCase(rawDataSet.genenames[j])){
											tfexpression.add(rawDataSet.data[j]);
											tflist.add(tf);
											found = true;
											break;
										}
									}
									if(!found){
										misstf.add(tf);
									}
									found = false;
								}

								for(String tg:TG){
									for(int j=0;j<rawDataSet.genenames.length;j++){
										if(tg.equalsIgnoreCase(rawDataSet.genenames[j])){
											tgexpression.add(rawDataSet.data[j]);
											tglist.add(tg);
											break;
										}
									}
								}

								for(String pg:PG){
									for(int j=0;j<rawDataSet.genenames.length;j++){
										if(pg.equalsIgnoreCase(rawDataSet.genenames[j])){
											pgexpression.add(rawDataSet.data[j]);
											pglist.add(pg);
											break;
										}
									}
								}
								
								JFrame frame = new JFrame("HeatMaps");
								frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
								frame.setLocation(20, 50);
								frame.setMinimumSize(new Dimension(500, 0));
								Container theDialogContainer = frame.getContentPane();
								theDialogContainer.setBackground(Color.white);
								
								DrawHeatMap newContentPane1 = new DrawHeatMap(
											frame, frame, rawDataSet.dsamplemins, tflist, tglist, pglist, 
											misstf, ndepth, tfexpression, tgexpression, pgexpression);

								theDialogContainer.add(newContentPane1);
								frame.setContentPane(theDialogContainer);
								frame.pack();
								frame.setVisible(true);
							}
						});
						
					}else{
						JOptionPane.showMessageDialog(null, "The selected pathways interact with no TFs!", "error", JOptionPane.ERROR_MESSAGE);
					}
				}else{
					Set<String> TF = new HashSet<String>();
					Set<String> TG = new HashSet<String>();
					List<String> tflist = new ArrayList<String>();
					List<String> tglist = new ArrayList<String>();
					
					for(int p=0;p<pathways.size();p++){
						List<List<String>> gs1 = tf_tg.get(pathways.get(p));
						if(gs1 != null && gs1.size()>0){
							for(int i=0;i<gs1.size();i++){
								List<String> a = gs1.get(i);
								TF.add(a.get(0));
								TG.add(a.get(1));
						}
						}
					}
					
					if(TF!=null && TF.size()>0){
						javax.swing.SwingUtilities.invokeLater(new Runnable() {
							public void run() {
								List<double[]> tfexpression = new ArrayList<double[]>();
								List<double[]> tgexpression = new ArrayList<double[]>();
								List<String> misstf = new ArrayList<String>();
								boolean found = false;
								
								for(String tf:TF){
									for(int j=0;j<rawDataSet.genenames.length;j++){
										if(tf.equalsIgnoreCase(rawDataSet.genenames[j])){
											tfexpression.add(rawDataSet.data[j]);
											tflist.add(tf);
											found = true;
											break;
										}
									}
									if(!found){
										misstf.add(tf);
									}
									found = false;
								}
								
								for(String tg:TG){
									for(int j=0;j<rawDataSet.genenames.length;j++){
										if(tg.equalsIgnoreCase(rawDataSet.genenames[j])){
											tgexpression.add(rawDataSet.data[j]);
											tglist.add(tg);
											break;
										}
									}
								}
								
								JFrame frame = new JFrame("HeatMaps");
								//JDialog frame = new JDialog(frame1, "Module Table", true);
								frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
								frame.setLocation(20, 50);
								Container theDialogContainer = frame.getContentPane();
								theDialogContainer.setBackground(Color.white);
								
								DrawHeatMap newContentPane1 = new DrawHeatMap(
											frame, frame, rawDataSet.dsamplemins, tflist, tglist, null, 
											misstf, ndepth, tfexpression, tgexpression, null);

								theDialogContainer.add(newContentPane1);
								frame.setContentPane(theDialogContainer);
								frame.pack();
								frame.setVisible(true);
							}
						});
						
					}else{
						JOptionPane.showMessageDialog(null, "The selected pathways interact with no TFs!", "error", JOptionPane.ERROR_MESSAGE);
					}
				}
			}else{
				JOptionPane.showMessageDialog(null, "No Pathway Selected!", "error", JOptionPane.ERROR_MESSAGE); 
			}
		}else if (szCommand.equals("network")) {
			int[] rows = table.getSelectedRows();
			List<String> pathways = new ArrayList<String>();
			for(int i=0;i<rows.length;i++){
				pathways.add(table.getValueAt(rows[i],0).toString());
			}
			if(pathways.size()>0){
				if(temp == 0){
					List<List<String>> gene1gene2 = new ArrayList<List<String>>();
					List<List<String>> genetype = new ArrayList<List<String>>();
					String[][] selectedtable = new String[pathways.size()][6];
					
					for(int p=0;p<pathways.size();p++){
						for(int i=0;i<ResultList.length;i++){
							if(ResultList[i][0].equals(pathways.get(p))){
								selectedtable[p] = ResultList[i];
								break;
							}
						}
						List<List<String>> gs1 = tf_tg.get(pathways.get(p));
						if(gs1 != null && gs1.size()>0){
							for(int i=0;i<gs1.size();i++){
								List<String> a = gs1.get(i);
								if(!gene1gene2.contains(a)){
									gene1gene2.add(a);
								}
								List<String> tftype = new ArrayList<String>();
								tftype.add(a.get(0));
								tftype.add("TF");
								List<String> tgtype = new ArrayList<String>();
								tgtype.add(a.get(1));
								tgtype.add("TG");
								
								if(!genetype.contains(tftype)){
									genetype.add(tftype);
								}
								if(!genetype.contains(tgtype)){
									genetype.add(tgtype);
								}
						}
						}
						
						List<List<String>> gs2 = tf_pg.get(pathways.get(p));
						if(gs2 != null && gs2.size()>0){
							for(int i=0;i<gs2.size();i++){
								List<String> a = gs2.get(i);
								if(!gene1gene2.contains(a)){
									gene1gene2.add(a);
								}
								
								List<String> tgtype = new ArrayList<String>();
								tgtype.add(a.get(0));
								tgtype.add("TF");
								List<String> pgtype = new ArrayList<String>();
								pgtype.add(a.get(1));
								pgtype.add("PG");
								
								if(!genetype.contains(tgtype)){
									genetype.add(tgtype);
								}
								if(!genetype.contains(pgtype)){
									genetype.add(pgtype);
								}
							}
						}
					}
					
					if(gene1gene2.size()>0 && genetype.size()>0){
						javax.swing.SwingUtilities.invokeLater(new Runnable() {
							public void run() {
								drawNet dn = new drawNet();
								dn.drawPic(genetype, gene1gene2, selectedtable, pathwaygene, ndepth, rawDataSet);
							}
						});
					}else{
						JOptionPane.showMessageDialog(null, "No Target Gene!", "error", JOptionPane.ERROR_MESSAGE);
					}
				}else{
					List<List<String>> gene1gene2 = new ArrayList<List<String>>();
					List<List<String>> genetype = new ArrayList<List<String>>();
					String[][] selectedtable = new String[pathways.size()][6];
					
					for(int p=0;p<pathways.size();p++){
						for(int i=0;i<ResultList.length;i++){
							if(ResultList[i][0].equals(pathways.get(p))){
								selectedtable[p] = ResultList[i];
								break;
							}
						}
						
						List<List<String>> gs1 = tf_tg.get(pathways.get(p));
						if(gs1 != null && gs1.size()>0){
							
							for(int i=0;i<gs1.size();i++){
								List<String> a = gs1.get(i);
								
								if(!gene1gene2.contains(a)){
									gene1gene2.add(a);
								}
								List<String> tftype = new ArrayList<String>();
								tftype.add(a.get(0));
								tftype.add("TF");
								List<String> tgtype = new ArrayList<String>();
								tgtype.add(a.get(1));
								tgtype.add("TG");
								
								if(!genetype.contains(tftype)){
									genetype.add(tftype);
								}
								if(!genetype.contains(tgtype)){
									genetype.add(tgtype);
								}
						}
						}

					}

					try{
						outXml = new BufferedWriter(new FileWriter("D:\\genetype.txt"));
						for(int t=0;t<genetype.size();t++){
	
								outXml.write(genetype.get(t).get(0));
								outXml.write("\t");
								outXml.write(genetype.get(t).get(1));
								outXml.newLine();
						}
						outXml.flush(); 
						outXml.close();
						System.out.println("DONE");
						
						outXml = new BufferedWriter(new FileWriter("D:\\gene1gene2.txt"));
						for(int t=0;t<gene1gene2.size();t++){
	
								outXml.write(gene1gene2.get(t).get(0));
								outXml.write("\t");
								outXml.write(gene1gene2.get(t).get(1));
								outXml.newLine();
						}
						outXml.flush(); 
						outXml.close();
						System.out.println("DONE");
					
						
						
					}catch (Exception e1) { 
						System.out.println("FALSE"); 
					e1.printStackTrace(); 
					 }
					
					if(gene1gene2.size()>0 && genetype.size()>0){
						javax.swing.SwingUtilities.invokeLater(new Runnable() {
							public void run() {
								drawNet2 dn = new drawNet2();
								dn.drawPic(genetype, gene1gene2, selectedtable, pathwaygene, ndepth, rawDataSet);
							}
						});
						
					}else{
						JOptionPane.showMessageDialog(null, "No Target Gene!", "error", JOptionPane.ERROR_MESSAGE);
					}
					
				}
				
				
			}else{
				JOptionPane.showMessageDialog(null, "No Pathway Selected!", "error", JOptionPane.ERROR_MESSAGE); 
			}

		} else if (szCommand.equals("copy")) {
			writeToClipboard();
		} else if (szCommand.equals("save")) {
			try {
				int nreturnVal = TSMiner_IO.theChooser.showSaveDialog(this);
				if (nreturnVal == JFileChooser.APPROVE_OPTION) {
					File f = TSMiner_IO.theChooser.getSelectedFile();
					PrintWriter pw = new PrintWriter(new FileOutputStream(f));
					if (szCommand.equals("save")) {
						printFile(pw);
					}
					pw.close();
				}
			} catch (final FileNotFoundException fex) {
				javax.swing.SwingUtilities.invokeLater(new Runnable() {
					public void run() {
						JOptionPane.showMessageDialog(null, fex.getMessage(),
								"Exception thrown", JOptionPane.ERROR_MESSAGE);
					}
				});
				fex.printStackTrace(System.out);
			}
		} else if (szCommand.equals("help")) {
			String szMessage = "This table gives information about the TFs regulating genes "
					+ "on the selected path.  Consult section 4.13 of the user manual for more details on this table.  ";

			Util.renderDialog(dialog, szMessage);// textArea);
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
	
	
	public Double mi3(Integer[] y1, Integer[] y2, Integer[] x) {
        double I3 = 0.0;
        
        MIM mim = new MIM();
        double Hy1 = mim.Entropy(y1,null,null);
        double Hy2 = mim.Entropy(y2,null,null);
        double Hx = mim.Entropy(x,null,null);
        double Hxy1 = mim.Entropy(x,y1,null);
        double Hxy2 = mim.Entropy(x,y2,null);
        double Hy2y1 = mim.Entropy(y2,y1,null);
        
        double Hxy2y1 = mim.Entropy(x,y2,y1);
        double enty1_y2x = Hxy2y1-Hxy2;
        double enty2_y1x = Hxy2y1-Hxy1;
        double entx_y1y2 = Hxy2y1-Hy2y1;
        double Iy1y2 = Hy2+Hy1-Hy2y1;
        double Iy1y2_x = Hxy1+Hxy2-Hxy2y1-Hx;
        double denominator = enty1_y2x+enty2_y1x+entx_y1y2;
        double Iy1y2x = Iy1y2 - Iy1y2_x;
 
        if(denominator == 0){
       	 I3 = Double.POSITIVE_INFINITY;
        }else{
       	 I3 = Math.abs(Iy1y2x/denominator);	 
        }
        return I3;
   }
    
    public Double[] random(double[] x1, double[] x2, double[] x3) {
    	double[] pvalue = new double[1000];
    	int N = x1.length;
    	
    	MIM mim = new MIM();
    	Integer[] y1 = mim.discretize(x1);
    	Integer[] y2 = mim.discretize(x2);
    	Integer[] x = mim.discretize(x3);

    	List<Integer> e1 = new ArrayList<Integer>();
    	List<Integer> e2 = new ArrayList<Integer>();
    	List<Integer> e3 = new ArrayList<Integer>();
    	for(int i=0;i<N;i++) {
    		e1.add(y1[i]);
    		e2.add(y2[i]);
    		e3.add(x[i]);
    	}
    	
    	Integer[] Y1 = new Integer[N];
		Integer[] Y2 = new Integer[N];
		Integer[] X = new Integer[N];
		
    	for(int i=0;i<pvalue.length;i++){
    		Collections.shuffle(e1);
    		Collections.shuffle(e2);
    		Collections.shuffle(e3);
    		for(int j=0;j<N;j++){
    			Y1[j] = e1.get(j);
    			Y2[j] = e2.get(j);
    			X[j] = e3.get(j);
    		}
    		pvalue[i] = mi3(Y1,Y2,X);
    	}
    	
    	double sum=0;
    	int count=0;
    	for(int i=0;i<pvalue.length;i++){
    		if(pvalue[i] != Double.POSITIVE_INFINITY){
    			sum += pvalue[i];
    			count++;
    		}
    	}
    	double dmean = sum/count;
    	double dVar=0;  
        for(int i=0;i<pvalue.length;i++){
        	if(pvalue[i] != Double.POSITIVE_INFINITY){
        		dVar+=(pvalue[i]-dmean)*(pvalue[i]-dmean); 
        	}
        }
    	double var = Math.sqrt(dVar/(count-1));
    	Double[] result = {dmean, var};
    	return result;
    }

    public static void setColumnColor(JTable table) {
        try
        {
            int cols = table.getColumnCount();
            DefaultTableCellRenderer tcr = new DefaultTableCellRenderer(){
                private static final long serialVersionUID = 1L;
                public Component getTableCellRendererComponent(JTable table,Object value, boolean isSelected, boolean hasFocus,int row, int column){
                	if(table.getValueAt(row, cols-1) != null && 
                			!table.getValueAt(row, cols-1).toString().equals("NA")){
                		if(Double.parseDouble(table.getValueAt(row, cols-1).toString()) > 0)
                            setBackground(Color.YELLOW);
                        else{
                        	setBackground(Color.CYAN);
                        }
                	}else{
                		setBackground(Color.LIGHT_GRAY);
                	}
                	
                    return super.getTableCellRendererComponent(table, value,isSelected, hasFocus, row, column);
                }
            };
            for(int i = 0; i < table.getColumnCount(); i++) {
                table.getColumn(table.getColumnName(i)).setCellRenderer(tcr);
            }
            tcr.setHorizontalAlignment(JLabel.CENTER);
        }
        catch (Exception e){
            e.printStackTrace();
        }
    }
    
    public String[][] BubbleSort(String[][] r, Integer n, Integer col) //升序冒泡排序
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

    public Integer getrank(Double[] pvalue, double value){
    	int cc = 1;
    	for(int i=0;i<pvalue.length;i++){
    		if(pvalue[i]<value){
    			cc++;
    		}
    	}
    	return cc;
    }
    
}