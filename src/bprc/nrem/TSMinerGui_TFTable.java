package bprc.nrem;

import bprc.core.*;
import edu.umd.cs.piccolo.nodes.PPath;
import edu.umd.cs.piccolo.PNode;

import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.SpinnerNumberModel;
import javax.swing.SpringLayout;

import java.awt.Color;
import java.awt.Component;
import java.awt.Container;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.*;
import java.util.List;
import java.util.Map.Entry;
import java.text.NumberFormat;
import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.awt.datatransfer.*;

import javax.swing.*;
import javax.swing.table.*;

import java.lang.reflect.InvocationTargetException;
import java.math.*;

/**
 * Class for a table that shows enrichment of TF targets along a path
 */
public class TSMinerGui_TFTable extends JPanel implements ActionListener {

	TSMiner_Timeiohmm theTimeiohmm;
	TSMiner_Timeiohmm.Treenode treeroot;
	JDialog theframe;
	JDialog theframe2;
	JTable tableposi;
	JTable tablenega;
	String[] columnNames;
	String[][] tabledata1;
	String[][] tabledata2;
	JButton copyButton1;
	JButton saveButton1;
	JButton copyButton2;
	JButton saveButton2;
	JButton pathwayButton;
	JLabel l1 = new JLabel("Minimum number of pathway genes interacting with the selected TF(s): ");
	JLabel l2 = new JLabel("Minimum percentage of pathway genes interacting with the selected TF(s): ");
	JLabel l4 = new JLabel("%");
	JLabel l3 = new JLabel("Find the pathways interacting with the selected TF(s) ");
	JSpinner j1, j2;
	int mininumber = 20;
	int minipersentage = 50;
	int PGnum = 20;
	int PGper = 50;
	//String label_method[] = { "Mutual interaction measure", "Overlapping genes"};
	//JComboBox comboBox_method = new JComboBox(label_method);
	TableSorter sorter1;
	TableSorter sorter2;
	final static Color bgColor = Color.white;
	final static Color fgColor = Color.black;
	int numrowsposi;
	int numrowsnega;
	int ndepth;
	String time;
	NumberFormat nf7;
	NumberFormat nf2;
	NumberFormat nf3;
	boolean bsplit;
	double enrichmentthreshold = 0.01;
	int npermutationval;
	JDialog executeDialognf = null;
	/**
	 * Constructor - builds the table
	 */
	public TSMinerGui_TFTable(TSMiner_Timeiohmm.Treenode treecopy, 
			int ndepth, JDialog theframe, TSMiner_Timeiohmm theTimeiohmm, 
			List<String[]> table1, List<String[]> table2, String time, int npermutationval) {

		this.theframe = theframe;
		this.theTimeiohmm = theTimeiohmm;
		this.ndepth = ndepth;
		this.treeroot = treecopy;
		this.time = time;
		this.enrichmentthreshold = Double.parseDouble(TSMiner_IO.enrichmentthreshold);
		this.npermutationval = npermutationval;
		setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
		setBackground(bgColor);
		setForeground(fgColor);
		numrowsposi = table1.size();
		numrowsnega = table2.size();
		nf2 = NumberFormat.getInstance(Locale.ENGLISH);
		nf2.setMinimumFractionDigits(2);
		nf2.setMaximumFractionDigits(2);
		nf3 = NumberFormat.getInstance(Locale.ENGLISH);
		nf3.setMinimumFractionDigits(4);
		nf3.setMaximumFractionDigits(4);
		nf7 = NumberFormat.getInstance(Locale.ENGLISH);
		nf7.setMinimumFractionDigits(7);
		nf7.setMaximumFractionDigits(7);
		SpinnerNumberModel mininum = new SpinnerNumberModel(new Integer(mininumber), new Integer(0), new Integer(100), new Integer(1));
	    j1 = new JSpinner(mininum);
		SpinnerNumberModel miniper = new SpinnerNumberModel(new Integer(minipersentage), new Integer(0), new Integer(100), new Integer(1));
		j2 = new JSpinner(miniper);
		columnNames = new String[4];
		columnNames[0] = "TF";
		columnNames[1] = "Enrichment q-value";
		columnNames[2] = "DE q-value";
		columnNames[3] = "Avg. Fold-change";
		
		if(table1.size()>0){
			int numrows = table1.size();
			int colnum = table1.get(0).length;
			tabledata1 = new String[numrows][colnum];
			for (int nrow = 0; nrow < numrows; nrow++) {
				String[] value = table1.get(nrow);
				for(int ncol = 0; ncol<colnum; ncol++){
					tabledata1[nrow][ncol] = value[ncol];
				}
			}
			sorter1 = new TableSorter(new TableModelST(tabledata1, columnNames));
			tableposi = new JTable(sorter1);
			TableColumn column;
			for (int nindex = 0; nindex < columnNames.length; nindex++) {
				column = tableposi.getColumnModel().getColumn(nindex);
				column.setPreferredWidth(100);
			}
			sorter1.setTableHeader(tableposi.getTableHeader());
			setColumnColor(tableposi);
			JScrollPane scrollPane = new JScrollPane(tableposi);
			tableposi.setPreferredScrollableViewportSize(new Dimension(700, Math.min(
					(tableposi.getRowHeight() + tableposi.getRowMargin()) * tableposi.getRowCount(), 150)));
			add(scrollPane);	
		}
		
		if(table2.size()>0){
			int numrows = table2.size();
			int colnum = table2.get(0).length;
			tabledata2 = new String[numrows][colnum];
			for (int nrow = 0; nrow < numrows; nrow++) {
				String[] value = table2.get(nrow);
				for(int ncol = 0; ncol<colnum; ncol++){
					tabledata2[nrow][ncol] = value[ncol];
				}
			}
			sorter2 = new TableSorter(new TableModelST(tabledata2, columnNames));
			tablenega = new JTable(sorter2);
			TableColumn column;
			for (int nindex = 0; nindex < columnNames.length; nindex++) {
				column = tablenega.getColumnModel().getColumn(nindex);
				column.setPreferredWidth(100);
			}
			sorter2.setTableHeader(tablenega.getTableHeader());
			setColumnColor(tablenega);
			JScrollPane scrollPane = new JScrollPane(tablenega);
			tablenega.setPreferredScrollableViewportSize(new Dimension(700, Math.min(
					(tablenega.getRowHeight() + tablenega.getRowMargin()) * tablenega.getRowCount(), 150)));
			add(scrollPane);
		}		
		addBottom();
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
                			String tf = table.getValueAt(row, 0).toString();
                			if(tf.contains("posi")){
                				if(Double.parseDouble(table.getValueAt(row, cols-1).toString()) > 0){
                                    setBackground(Color.YELLOW);
                    			}else{
                                	setBackground(Color.CYAN);
                                }
                			}else{
                				if(Double.parseDouble(table.getValueAt(row, cols-1).toString()) > 0){
                                    setBackground(Color.CYAN);
                    			}else{
                               	setBackground(Color.YELLOW);
                                }
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
	
	/**
	 * Helper function that adds information displayed at the bottom of the
	 * table information window
	 */
	private void addBottom() {
		JPanel countPanel = new JPanel();
		
		String szcountLabel = "There are " + numrowsposi+" positive TFs and "+numrowsnega+" negative TFs";
		JLabel countLabel = new JLabel(szcountLabel);

		countPanel.setBackground(Color.white);
		countPanel.add(countLabel);
		countPanel.setMaximumSize(new Dimension(Integer.MAX_VALUE, 20));
		add(countPanel);
		/*
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
		buttonPanel.add(copyButton);
		buttonPanel.add(saveButton);
		buttonPanel.setMaximumSize(new Dimension(Integer.MAX_VALUE, 20));
		add(buttonPanel);
		*/
		JPanel pathwayPanel1 = new JPanel();
		pathwayPanel1.setBackground(Color.white);
		pathwayPanel1.add(l1);
		pathwayPanel1.add(j1);
		pathwayPanel1.setMaximumSize(new Dimension(Integer.MAX_VALUE, 20));
		add(pathwayPanel1);
		
		JPanel pathwayPanel2 = new JPanel();
		pathwayPanel2.setBackground(Color.white);
		pathwayPanel2.add(l2);
		pathwayPanel2.add(j2);
		pathwayPanel2.add(l4);
		pathwayPanel2.setMaximumSize(new Dimension(Integer.MAX_VALUE, 20));
		add(pathwayPanel2);
		
		pathwayButton = new JButton("Run");
		pathwayButton.setActionCommand("pathway");
		pathwayButton.setMinimumSize(new Dimension(800, 20));
		pathwayButton.addActionListener(this);
		
		JPanel pathwayPanel3 = new JPanel();
		pathwayPanel3.setBackground(Color.white);
		pathwayPanel3.add(l3);
		pathwayPanel3.add(pathwayButton);
		pathwayPanel3.setMaximumSize(new Dimension(Integer.MAX_VALUE, 20));
		add(pathwayPanel3);
	}
    /*
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
		Clipboard systemClipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
		Transferable transferableText = new StringSelection(sbuf.toString());
		systemClipboard.setContents(transferableText, null);
	}
	*/
	/**
	 * Responds to buttons being pressed on the interface
	 */
	public void actionPerformed(ActionEvent e) {
		String szCommand = e.getActionCommand();
		if (szCommand.equals("pathway")) {
			String selecteditem = theTimeiohmm.theDataSet.dsamplemins[ndepth];
			if((tableposi == null || tableposi.getSelectedRows()==null || tableposi.getSelectedRows().length==0)
			&& (tablenega == null || tablenega.getSelectedRows()==null || tablenega.getSelectedRows().length==0)){
				JOptionPane.showMessageDialog(theframe2, 
						"Please select one or more TF(s)!", "Information", JOptionPane.INFORMATION_MESSAGE);
			}else{
				if (theTimeiohmm.raw_gobinding == null || 
						theTimeiohmm.raw_gobinding.reg2GeneBindingIndex[ndepth] == null) {
						throw new IllegalArgumentException(     
								"No gene annotation input file given!");
				}
				/*************************************/

					HashMap<String,Integer> TFindex = theTimeiohmm.raw_reg2DataSetIndex; //用rawDataSet
					int[][] Pathway2Geneindex = theTimeiohmm.raw_gobinding.reg2GeneBindingIndex[ndepth];//用rawDataSet
					TSMiner_DataSet rawDataSet = theTimeiohmm.rawDataSet;
					
					int numgenes = rawDataSet.data.length;
					int num1=(int)(Math.random()*numgenes);
					int num2=(int)(Math.random()*numgenes);
					int num3=(int)(Math.random()*numgenes);
					double[] x1 = rawDataSet.data[num1];
					double[] x2 = rawDataSet.data[num2];
					double[] x3 = rawDataSet.data[num3];
					
					MIM mim = new MIM();
					Double[] permute = mim.random(x1,x2,x3);
					//selected TFs
					
					//int[][] pathwaycount = new int[Pathway2Geneindex.length][3];
					HashMap<Integer, List<List<String>>> tf_tg = new HashMap<Integer, List<List<String>>>();//结果网络中画的overlap节点
					HashMap<Integer, List<List<String>>> tf_pg = new HashMap<Integer, List<List<String>>>();//结果网络中画的overlap节点
					HashMap<String, List<List<String>>> tf_tg2 = new HashMap<String, List<List<String>>>();
					HashMap<String, List<List<String>>> tf_pg2 = new HashMap<String, List<List<String>>>();
					HashMap<Integer, List<Integer[]>> twogene = new HashMap<Integer, List<Integer[]>>();//用于构成三联基因的前两个基因
					HashMap<Integer, List<List<Integer>>> threegene = new HashMap<Integer, List<List<Integer>>>();//各通路的三联基因
					HashMap<Integer, List<Integer>> geneoverlap = new HashMap<Integer, List<Integer>>();//某通路overlap的基因编号
					HashMap<Integer, Double> triple_pval = new HashMap<Integer, Double>();
					HashSet<Integer> triple_miss = new HashSet<Integer>();
					boolean upregulated = false;
					boolean downregulated = false;
					
						HashMap<String,Set<Integer>> regulatedGene = new HashMap<String,Set<Integer>>();
						FindTGs(regulatedGene, treeroot, 0);
						List<String> TFs = new ArrayList<String>();
						if(tableposi != null && tableposi.getSelectedRows()!=null){
							int[] rows1 = tableposi.getSelectedRows();
							for(int i=0;i<rows1.length;i++){
								TFs.add(tableposi.getValueAt(rows1[i],0).toString().toUpperCase());
								String fc = (String) tableposi.getValueAt(rows1[i],3);
								if(!fc.equals("NA")){
									if(Double.parseDouble(fc)>0){
										upregulated = true;
									}else{
										downregulated = true;
									}
								}
							}
						}
						if(tablenega != null && tablenega.getSelectedRows()!=null){
							int[] rows2 = tablenega.getSelectedRows();
							for(int i=0;i<rows2.length;i++){
								TFs.add(tablenega.getValueAt(rows2[i],0).toString().toUpperCase());
								String fc = (String) tablenega.getValueAt(rows2[i],3);
								if(!fc.equals("NA")){
									if(Double.parseDouble(fc)>0){
										upregulated = true;
									}else{
										downregulated = true;
									}
								}
							}
						}
						for(int i=0;i<TFs.size();i++){
							String tfname = TFs.get(i);
							if(regulatedGene.get(tfname) != null && regulatedGene.get(tfname).size() > 0){
								Set<Integer> TGs = regulatedGene.get(tfname);
								if(TFindex.get(tfname) != null){
									int TFr = TFindex.get(tfname);
									for(int p=0;p<Pathway2Geneindex.length;p++){
										if(Pathway2Geneindex[p].length>=5){
												for(int bb:Pathway2Geneindex[p]){
													if(bb == TFr){
														if(geneoverlap.get(p)!=null){
															List<Integer> gene = geneoverlap.get(p);
															if(!gene.contains(bb)){
																gene.add(bb);
																geneoverlap.put(p, gene);
																//pathwaycount[p][1]++;
															}
														}else{
															List<Integer> gene = new ArrayList<Integer>();
															gene.add(bb);
															geneoverlap.put(p, gene);
															//pathwaycount[p][1]++;
														}
													}else if(TGs.contains(bb)){
															Integer[] aa = new Integer[2];
															aa[0] = TFr;
															aa[1] = bb;
															if(twogene.get(p) != null){
																List<Integer[]> tg = twogene.get(p);
																tg.add(aa);
																twogene.put(p, tg);
															}else{
																List<Integer[]> tg = new ArrayList<Integer[]>();
																tg.add(aa);
																twogene.put(p, tg);
															}
															
															if(geneoverlap.get(p)!=null){
																List<Integer> gene = geneoverlap.get(p);
																if(!gene.contains(bb)){
																	gene.add(bb);
																	geneoverlap.put(p, gene);
																	//pathwaycount[p][1]++;
																}
															}else{
																List<Integer> gene = new ArrayList<Integer>();
																gene.add(bb);
																geneoverlap.put(p, gene);
																//pathwaycount[p][1]++;
															}
															
															List<String> TF_TG = new ArrayList<String>();
															TF_TG.add(tfname.split("_")[0]);
															TF_TG.add(rawDataSet.genenames[bb].toUpperCase());
															if(tf_tg.get(p)!=null){
																List<List<String>> tfgene = tf_tg.get(p);
																tfgene.add(TF_TG);
																tf_tg.put(p, tfgene);
															}else{
																List<List<String>> tfgene = new ArrayList<List<String>>();
																tfgene.add(TF_TG);
																tf_tg.put(p, tfgene);
															}
													}
												}
										}
									}
									}else{
										for(int p=0;p<Pathway2Geneindex.length;p++){
											if(Pathway2Geneindex[p].length>=5){
												for(int bb:Pathway2Geneindex[p]){
													if(TGs.contains(bb)){
														if(geneoverlap.get(p)!=null){
															List<Integer> gene = geneoverlap.get(p);
															if(!gene.contains(bb)){
																gene.add(bb);
																geneoverlap.put(p, gene);
																//pathwaycount[p][1]++;
															}
														}else{
															List<Integer> gene = new ArrayList<Integer>();
															gene.add(bb);
															geneoverlap.put(p, gene);
															//pathwaycount[p][1]++;
														}
														List<String> TF_TG = new ArrayList<String>();
														TF_TG.add(tfname.split("_")[0]);
														TF_TG.add(rawDataSet.genenames[bb].toUpperCase());
														if(tf_tg.get(p)!=null){
															List<List<String>> tfgene = tf_tg.get(p);
															tfgene.add(TF_TG);
															tf_tg.put(p, tfgene);
														}else{
															List<List<String>> tfgene = new ArrayList<List<String>>();
															tfgene.add(TF_TG);
															tf_tg.put(p, tfgene);
														}
													}
												}
											}	
										}
								  }
						   }
					}
						
					for(Entry<Integer, List<Integer[]>> entry:twogene.entrySet())
					  {
						Integer p = entry.getKey();
						List<Integer[]> gs = entry.getValue();
						List<List<Integer>> three = new ArrayList<List<Integer>>();
						//List<Integer> threecode = new ArrayList<Integer>();
						List<Integer> genelist = geneoverlap.get(p);
						for(int i=0;i<gs.size();i++){
							int TFr = gs.get(i)[0];
							int TGr = gs.get(i)[1];
								for(int PGr:Pathway2Geneindex[p]){
									if(!genelist.contains(PGr)){
										int code = TFr*TGr*PGr+TFr+TGr+PGr;
										List<Integer> triple = new ArrayList<Integer>();
										triple.add(code);
										triple.add(TFr); //TF index
										triple.add(TGr); //TG index
										triple.add(PGr); //PG index
										three.add(triple);
									}
							  }
						}
						threegene.put(p, three);
					  }
					
					for(Entry<Integer, List<List<Integer>>> entry:threegene.entrySet())
					{
						Integer p = entry.getKey();
						List<List<Integer>> three = entry.getValue();
						int threenum = three.size();
						List<List<Integer>> result = new ArrayList<List<Integer>>();
						//List<Double> pvallist = new ArrayList<Double>();
						for(int i=0;i<threenum;i++){
							List<Integer> triple = three.get(i);
							if(triple_pval.get(triple.get(0)) == null){
								if(!triple_miss.contains(triple.get(0))){
									double[] TFval = rawDataSet.data[triple.get(1)]; //TF expression
									double[] TGval = rawDataSet.data[triple.get(2)]; //TG expression
									double[] PGval = rawDataSet.data[triple.get(3)]; //PG expression
									double pvalue = mim.getPvalue(TFval, TGval, PGval, permute[0], permute[1]);
	                                if(pvalue>-1){
	                                	triple_pval.put(triple.get(0), pvalue);
	                                	double qvalue = pvalue*threenum;
	                                	if(qvalue<=0.01){
	    									result.add(triple);
	    								}
	                                	//pvallist.add(pvalue);
	                                }else{
										triple_miss.add(triple.get(0));
									}
								}
							}else{
								double pvalue = triple_pval.get(triple.get(0));
								double qvalue = pvalue*threenum;
                            	if(qvalue<=0.05){
									result.add(triple);
								}
								//pvallist.add(pvalue);
							}
						}
						if(result.size()>0) threegene.put(p, result);
						/*
						for(int i=0;i<three.size();i++){
							List<Integer> triple = three.get(i);
							if(triple_pval.get(triple.get(0)) != null){
								double pvalue = triple_pval.get(triple.get(0));				
								//int index = getrank(pvallist,pvalue);
								//double qvalue = pvalue*pvallist.size()/index;
								double qvalue = pvalue*pvallist.size();
								if(qvalue<=0.01){
									result.add(triple);
								}
							}
						}
						*/
					}
					
					for(Entry<Integer, List<List<Integer>>> entry:threegene.entrySet())
					  {
						Integer p = entry.getKey();
						List<List<Integer>> three = entry.getValue();
						List<Integer> genelist = geneoverlap.get(p);
						for(int i=0;i<three.size();i++){
							List<Integer> triple = three.get(i);
							int tfr = triple.get(1);
							int pgr = triple.get(3);
							if(!genelist.contains(pgr)) genelist.add(pgr);
								List<String> TF_PG = new ArrayList<String>();
								TF_PG.add(rawDataSet.genenames[tfr].toUpperCase());
								TF_PG.add(rawDataSet.genenames[pgr].toUpperCase());
								if(tf_pg.get(p)!=null){
									List<List<String>> tfpg = tf_pg.get(p);
									if(!tfpg.contains(TF_PG)){
										tfpg.add(TF_PG);
										tf_pg.put(p, tfpg);
									}	
								}else{
									List<List<String>> tfpg = new ArrayList<List<String>>();
									tfpg.add(TF_PG);
									tf_pg.put(p, tfpg);
								}
						}
						geneoverlap.put(p, genelist);
					}
					
					DEGeneset permutation = new DEGeneset();
					List<String[]> ResultList = new ArrayList<String[]>();
					PGnum = Integer.parseInt(j1.getValue().toString());
					PGper = Integer.parseInt(j2.getValue().toString());
					for(Entry<Integer, List<Integer>> entry:geneoverlap.entrySet())
					{
						Integer p = entry.getKey();
						List<Integer> genelist = entry.getValue();
						int sum = genelist.size();
						double per = (double)sum/(double)Pathway2Geneindex[p].length*100;
						
						if(sum>PGnum && per>PGper){
							String[] pathway = new String[6];			
							pathway[0] = theTimeiohmm.raw_gobinding.regNames[p];
							pathway[1] = sum+"";
							pathway[2] = nf2.format(per);
							double[] express = new double[genelist.size()*2];
							if(npermutationval == 1){
								for(int r=0;r<genelist.size();r++){
									express[r] = rawDataSet.data[genelist.get(r)][ndepth-1];
									express[r+genelist.size()] = rawDataSet.data[genelist.get(r)][ndepth];
								}
							}else{
								for(int r=0;r<genelist.size();r++){
									express[r] = 0;
									express[r+genelist.size()] = rawDataSet.data[genelist.get(r)][ndepth];
								}
							}
							Double[] result = permutation.permutetest(express);
							if(result[0] != null){
								pathway[3] = result[0]+"";
							}else{
								pathway[3] = "NA";
							}
							pathway[5] = nf3.format(result[1]);
							if((upregulated && result[1]>0) || (downregulated && result[1]<0)){
								ResultList.add(pathway);
							}
						}
					}
					
					String[][] Resulttable = new String[ResultList.size()][6];
					int index = 0;
					int pathcount=0;
					for(int i=0;i<Resulttable.length;i++){
						String[] aa = ResultList.get(i);
						if(!aa[3].equals("NA")){
							Resulttable[index] = aa;
							index++;
							pathcount++;
						}
					}
					for(int i=0;i<Resulttable.length;i++){
						String[] aa = ResultList.get(i);
						if(aa[3].equals("NA")){
							Resulttable[index] = aa;
							index++;
						}
					}
					Util.BubbleSort_inc(Resulttable, pathcount, 3);
					
					for(int i=0;i<pathcount;i++){
						double fdr = Double.parseDouble(Resulttable[i][3])*pathcount/(i+1);
						if(fdr<1){
							Resulttable[i][4] = nf3.format(fdr);
						}else{
							Resulttable[i][4] = 1+"";
						}
					}
					for(int i=pathcount;i<Resulttable.length;i++){
						Resulttable[i][4] = "NA";
					}

					HashMap<String, List<String>> pathwaygene = new HashMap<String, List<String>>();
					for(Entry<Integer, List<Integer>> entry:geneoverlap.entrySet())
					  {
						Integer path = entry.getKey();
						String name = theTimeiohmm.raw_gobinding.regNames[path];
						List<Integer> gs = entry.getValue();
						ArrayList<String> gs2 = new ArrayList<String>();
						for(int i=0;i<gs.size();i++){
							gs2.add(rawDataSet.genenames[gs.get(i)].toUpperCase());
						}
						pathwaygene.put(name, gs2);
					  }
					
					for(Entry<Integer, List<List<String>>> entry:tf_tg.entrySet())
					  {
						Integer path = entry.getKey();
						List<List<String>> gs = entry.getValue();
						String pathname = theTimeiohmm.raw_gobinding.regNames[path];
						tf_tg2.put(pathname, gs);
					  }
					for(Entry<Integer, List<List<String>>> entry:tf_pg.entrySet())
					  {
						Integer path = entry.getKey();
						List<List<String>> gs = entry.getValue();
						String pathname = theTimeiohmm.raw_gobinding.regNames[path];
						tf_pg2.put(pathname, gs);
					  }
					if(tableposi != null) tableposi.getSelectionModel().clearSelection();
					if(tablenega != null) tablenega.getSelectionModel().clearSelection();
					/*
					if (executeDialognf != null) {
	    	    		executeDialognf.dispose();
	    	    		executeDialognf.setVisible(false);
	    	    	}
					*/
					javax.swing.SwingUtilities.invokeLater(new Runnable() {
						public void run() {
							JFrame pathwayframe = new JFrame();
							JDialog frame = new JDialog(pathwayframe, "Pathway Table", true);
							frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
							frame.setLocation(20, 50);
							Container theDialogContainer = frame.getContentPane();
							theDialogContainer.setBackground(Color.white);
							JTabbedPane tabbedPane = new JTabbedPane();
							if(Resulttable!=null &&Resulttable.length>0){
								TSMinerGui_PathwayTable newContentPane1 = new TSMinerGui_PathwayTable(
										Resulttable, frame, pathwayframe,
										pathwaygene,tf_tg2,tf_pg2, ndepth, selecteditem, rawDataSet, 0);
								newContentPane1.setOpaque(true); // content panes must be opaque
								tabbedPane.addTab(time+" Pathway", null, newContentPane1,time+" Pathway");
								theDialogContainer.add(tabbedPane);
								pathwayframe.setContentPane(theDialogContainer);
								pathwayframe.pack();
								pathwayframe.setVisible(true);
							}else{
								JOptionPane.showMessageDialog(theframe2, 
										"No pathway found", "Information", JOptionPane.INFORMATION_MESSAGE); 
							}
						}
					});
			}
			
		} else if (szCommand.equals("help")) {
			String szMessage = "This table gives information about the TFs regulating genes "
					+ "on the selected path.  Consult section 4.13 of the user manual for more details on this table.  ";

			Util.renderDialog(theframe, szMessage);// textArea);
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
	public Integer getrank(List<Double> pvalue, double value){
    	int cc = 1;
    	for(int i=0;i<pvalue.size();i++){
    		if(pvalue.get(i)!=null && pvalue.get(i)<value){
    			cc++;
    		}
    	}
    	return cc;
    }
    
    public void FindTGs(HashMap<String,Set<Integer>> TFTGs, 
    		TSMiner_Timeiohmm.Treenode treeroot, int depth){
    	
    		if(depth>0 && treeroot.parent.numchildren>1 && treeroot.resultTF != null){
    			String[][] tftable = treeroot.resultTF;
    			HashMap<String,List<Integer>> regulatedGene = treeroot.regulatedGeneindex;
    			
    			for(int i=0;i<tftable.length;i++){
    				String tfname = tftable[i][0].toUpperCase();
    				if(Double.parseDouble(tftable[i][4]) <= enrichmentthreshold){
    					List<Integer> tgs = regulatedGene.get(tfname);
    					if(tgs!=null && tgs.size()>0){
    						if(TFTGs.get(tfname) != null){
        						Set<Integer> tgset = TFTGs.get(tfname);
        						tgset.addAll(tgs);
        						TFTGs.put(tfname, tgset);
        					}else{
        						Set<Integer> tgset = new HashSet<Integer>();
        						tgset.addAll(tgs);
        						TFTGs.put(tfname, tgset);
        					}
    					}
    				}
    			}
        	}
    	
    	if(depth<ndepth){
    		for (int nchild = 0; nchild < treeroot.numchildren; nchild++) {
        		FindTGs(TFTGs, treeroot.nextptr[nchild], depth+1);
            }
    	}
	}
    
}