package bprc.nrem;

import bprc.core.*;
import bprc.nrem.TSMiner_Timeiohmm.Treenode;

import javax.swing.JButton;
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
public class TSMinerGui_EdgeTable2 extends JPanel implements ActionListener {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	TSMiner_Timeiohmm theTimeiohmm;
	JDialog dialog;
	JFrame theframe;
	JFrame theframe2;
	String[] columnNames;
	String[][] tabledata;
	JButton copyButton;
	JButton saveButton;
	JLabel l1 = new JLabel("Minimum number of pathway genes interacted with the selected TFs: ");
	JLabel l2 = new JLabel("Minimum percentage of pathway genes interacted with the selected TFs: ");
	JLabel l4 = new JLabel("%");
	JLabel l3 = new JLabel("Find the pathways associated with the selected TFs ");
	JSpinner j1, j2;
	int mininumber = 20;
	int minipersentage = 50;
	int PGnum = 20;
	int PGper = 50;
	
	JButton pathwayButton;
	JScrollPane scrollPane;
	TableSorter sorter;
	JTable table;
	final static Color bgColor = Color.white;
	final static Color fgColor = Color.black;
	TSMinerGui_FilterStaticModel hmst;
	int numrows;
	TSMiner_Timeiohmm.Treenode ptr;
	int nchild;
	int ndepth;
	int fccount;
	NumberFormat nf7;
	NumberFormat nf2;
	NumberFormat nf3;
	boolean bsplit;
	boolean broot;
	int npermutationval;
	
	/**
	 * Constructor - builds the table
	 */
	public TSMinerGui_EdgeTable2(JDialog dialog, JFrame frame1,
			TSMiner_Timeiohmm theTimeiohmm, TSMiner_Timeiohmm.Treenode ptr,
			int ndepth, int nchild, boolean broot, int npermutationval) {
		this.dialog = dialog;
		this.theframe = frame1;
		this.theTimeiohmm = theTimeiohmm;
		this.ptr = ptr;
		this.ndepth = ndepth;
		this.nchild = nchild;
		this.broot = broot;
		this.npermutationval = npermutationval;
		
		setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
		setBackground(bgColor);
		setForeground(fgColor);
		SpinnerNumberModel mininum = new SpinnerNumberModel(new Integer(mininumber), new Integer(0), new Integer(100), new Integer(1));
	    j1 = new JSpinner(mininum);
		SpinnerNumberModel miniper = new SpinnerNumberModel(new Integer(minipersentage), new Integer(0), new Integer(100), new Integer(1));
		j2 = new JSpinner(miniper);

		numrows = theTimeiohmm.bindingData.regNames.length;
		fccount = theTimeiohmm.theDataSet.data[0].length-ndepth;

		nf2 = NumberFormat.getInstance(Locale.ENGLISH);
		nf2.setMinimumFractionDigits(2);
		nf2.setMaximumFractionDigits(2);
		nf3 = NumberFormat.getInstance(Locale.ENGLISH);
		nf3.setMinimumFractionDigits(4);
		nf3.setMaximumFractionDigits(4);
		nf7 = NumberFormat.getInstance(Locale.ENGLISH);
		nf7.setMinimumFractionDigits(7);
		nf7.setMaximumFractionDigits(7);
		
				if(ptr.nextptr[nchild].deTF[0].length > 5){
					columnNames = new String[7];
					tabledata = new String[numrows][columnNames.length];
					columnNames[0] = "TF";
					columnNames[1] = "Num total";
					columnNames[2] = "Num parent";
					columnNames[3] = "Num sub-path";
					columnNames[4] = "Enrichment p-value";
					columnNames[5] = theTimeiohmm.theDataSet.dsamplemins[ndepth]+" DE p-value";
					columnNames[6] = theTimeiohmm.theDataSet.dsamplemins[ndepth]+" Avg. fold-change";
					int nrowindex = 0;
					
						int hangshu=ptr.nextptr[nchild].deTF.length;
						for (int nrow = 0; nrow < hangshu; nrow++) {
							tabledata[nrowindex][0] = ptr.nextptr[nchild].deTF[nrow][0];
							tabledata[nrowindex][1] = ptr.nextptr[nchild].deTF[nrow][1];
							tabledata[nrowindex][2] = ptr.nextptr[nchild].deTF[nrow][2];
							tabledata[nrowindex][3] = ptr.nextptr[nchild].deTF[nrow][3];
							tabledata[nrowindex][4] = ptr.nextptr[nchild].deTF[nrow][4];
							tabledata[nrowindex][5] = ptr.nextptr[nchild].deTF[nrow][5];
							tabledata[nrowindex][6] = ptr.nextptr[nchild].deTF[nrow][6];
							nrowindex++;
					}
				}else{
					columnNames = new String[5];
					tabledata = new String[numrows][columnNames.length];
					columnNames[0] = "TF";
					columnNames[1] = "Num total";
					columnNames[2] = "Num parent";
					columnNames[3] = "Num sub-path";
					columnNames[4] = "Enrichment p-value";
					int nrowindex = 0;
					int hangshu=ptr.nextptr[nchild].deTF.length;
					for (int nrow = 0; nrow < hangshu; nrow++) {
						tabledata[nrowindex][0] = ptr.nextptr[nchild].deTF[nrow][0];
						tabledata[nrowindex][1] = ptr.nextptr[nchild].deTF[nrow][1];
						tabledata[nrowindex][2] = ptr.nextptr[nchild].deTF[nrow][2];
						tabledata[nrowindex][3] = ptr.nextptr[nchild].deTF[nrow][3];
						tabledata[nrowindex][4] = ptr.nextptr[nchild].deTF[nrow][4];
						nrowindex++;
				    }
				}
		//Color up = new Color(255, 0, 0, 120);
		//Color down = new Color(0, 0, 255, 120);

		sorter = new TableSorter(new TableModelST(tabledata, columnNames));
		table = new JTable(sorter);
		sorter.setTableHeader(table.getTableHeader());
		
		table.setPreferredScrollableViewportSize(new Dimension(800, Math.min((table.getRowHeight() + table.getRowMargin())
						* table.getRowCount(), 400)));

		table.setAutoResizeMode(JTable.AUTO_RESIZE_ALL_COLUMNS);
		TableColumn column;
		
		for (int ncolindex = 0; ncolindex < 4; ncolindex++) {
			column = table.getColumnModel().getColumn(ncolindex);
			column.setPreferredWidth(120);
		}
		column = table.getColumnModel().getColumn(4);
		column.setPreferredWidth(240);
		for (int ncolindex = 5; ncolindex < columnNames.length; ncolindex++) {
			column = table.getColumnModel().getColumn(ncolindex);
			column.setPreferredWidth(200);
		}

		if(ptr.nextptr[nchild].deTF[0].length > 5){
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
	
	/**
	 * Helper function that adds information displayed at the bottom of the
	 * table information window
	 */
	private void addBottom() {
		JPanel countPanel = new JPanel();
		JPanel labelPanel = new JPanel();
		String szcountLabel = "Total number of genes most likely going through this path is "
				+ ptr.nextptr[nchild].numPath;
		if (bsplit) {
			szcountLabel += " ("
					+ nf2.format(100 * (double) ptr.nextptr[nchild].numPath
							/ ptr.numPath) + "% of split genes)";
		}
		JLabel countLabel = new JLabel(szcountLabel);
		countPanel.setBackground(Color.white);
		countPanel.add(countLabel);
		countPanel.setMaximumSize(new Dimension(Integer.MAX_VALUE, 20));
		add(countPanel);
		String szInfo;
		if (broot) {
			szInfo = "Path output distribution at "
					+ theTimeiohmm.theDataSet.dsamplemins[0]
					+ " is Normal(mu =" + nf3.format(ptr.dmean) + ",sigma = "
					+ nf3.format(ptr.dsigma) + ")";
		} else {
			szInfo = "Path output distribution at "
					+ theTimeiohmm.theDataSet.dsamplemins[ptr.ndepth + 1]
					+ " is Normal(mu =" + nf3.format(ptr.nextptr[nchild].dmean)
					+ ",sigma = " + nf3.format(ptr.nextptr[nchild].dsigma) + ")";
		}

		JLabel infoLabel = new JLabel(szInfo);
		labelPanel.setBackground(Color.white);
		labelPanel.add(infoLabel);
		labelPanel.setMaximumSize(new Dimension(Integer.MAX_VALUE, 20));
		add(labelPanel);

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
		
		pathwayButton = new JButton("Run");
		pathwayButton.setActionCommand("pathway");
		pathwayButton.setMinimumSize(new Dimension(800, 20));
		pathwayButton.addActionListener(this);
		
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
		
		JPanel pathwayPanel3 = new JPanel();
		pathwayPanel3.setBackground(Color.white);
		pathwayPanel3.add(l3);
		pathwayPanel3.add(pathwayButton);
		pathwayPanel3.setMaximumSize(new Dimension(Integer.MAX_VALUE, 20));
		add(pathwayPanel3);
		
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
	 */
	public void actionPerformed(ActionEvent e) {
		
		String szCommand = e.getActionCommand();
		if (szCommand.equals("pathway")) {
			String selecteditem = theTimeiohmm.theDataSet.dsamplemins[ndepth];
			if(table.getSelectedRows()==null || table.getSelectedRows().length==0){
				JOptionPane.showMessageDialog(theframe, 
						"Please select one or more TFs!", "Information", JOptionPane.INFORMATION_MESSAGE);
			}else{
					if (theTimeiohmm.raw_gobinding == null || 
						theTimeiohmm.raw_gobinding.reg2GeneBindingIndex[ndepth] == null) {
						throw new IllegalArgumentException(       
								"No gene annotation input file given!");
					}
					
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
					
					Double[] permute = random(x1,x2,x3);
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
					//boolean upregulated = false;
					//boolean downregulated = false;
					
						HashMap<String,List<Integer>> regulatedGene = ptr.nextptr[nchild].regulatedGeneindex; //用rawDataSet
						int[] rows = table.getSelectedRows();
						List<String> TFs = new ArrayList<String>();
						for(int i=0;i<rows.length;i++){
							TFs.add(table.getValueAt(rows[i],0).toString().toUpperCase());
								//String fc = (String) table.getValueAt(rows[i],6);
								//if(!fc.equals("NA")){
								//	if(Double.parseDouble(fc)>0){
								//		upregulated = true;
								//	}else{
								//		downregulated = true;
								//	}
								//}
						}
						
						for(int i=0;i<TFs.size();i++){
							String tfname = TFs.get(i);
							if(regulatedGene.get(tfname) != null && regulatedGene.get(tfname).size() > 0){
								List<Integer> TGs = regulatedGene.get(tfname);
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
							List<Double> pvallist = new ArrayList<Double>();
							for(int i=0;i<three.size();i++){
								List<Integer> triple = three.get(i);
								if(triple_pval.get(triple.get(0)) == null){
									if(!triple_miss.contains(triple.get(0))){
										double[] TFval = rawDataSet.data[triple.get(1)]; //TF expression
										double[] TGval = rawDataSet.data[triple.get(2)]; //TG expression
										double[] PGval = rawDataSet.data[triple.get(3)]; //PG expression
										MIM mim = new MIM();
										double pvalue = mim.getPvalue(TFval, TGval, PGval, permute[0], permute[1]);
		                                if(pvalue>-1){
		                                	triple_pval.put(triple.get(0), pvalue);
		                                	pvallist.add(pvalue);
		                                }else{
											triple_miss.add(triple.get(0));
										}
									}
								}else{
									double pvalue = triple_pval.get(triple.get(0));
									pvallist.add(pvalue);
								}
							}
							List<List<Integer>> result = new ArrayList<List<Integer>>();
							for(int i=0;i<three.size();i++){
								List<Integer> triple = three.get(i);
								if(triple_pval.get(triple.get(0)) != null){
									double pvalue = triple_pval.get(triple.get(0));
									int index = getrank(pvallist,pvalue);
									double qvalue = pvalue*pvallist.size()/index;
									if(qvalue<=0.01){
										result.add(triple);
									}
								}
							}
							if(result.size()>0) threegene.put(p, result);
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
								//if((upregulated && result[1]>0) || (downregulated && result[1]<0)){
									ResultList.add(pathway);
								//}
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
								String segment = theTimeiohmm.theDataSet.dsamplemins[ndepth];
								TSMinerGui_PathwayTable newContentPane1 = new TSMinerGui_PathwayTable(
										Resulttable, frame, pathwayframe,
										pathwaygene,tf_tg2,tf_pg2, ndepth, selecteditem, rawDataSet, 0);
								newContentPane1.setOpaque(true); // content panes must be opaque
								tabbedPane.addTab(segment+" Pathway", null, newContentPane1,segment+" Pathway");
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
		}if (szCommand.equals("copy")) {
			writeToClipboard();
		}else if (szCommand.equals("save")) {
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
    		if(pvalue.get(i)!=null && pvalue.get(i)<value){
    			cc++;
    		}
    	}
    	return cc;
    }
    
}