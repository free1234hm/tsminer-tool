package bprc.nrem;

import edu.uci.ics.jung.algorithms.layout.KKLayout;
import edu.uci.ics.jung.graph.SparseMultigraph;
import edu.uci.ics.jung.visualization.VisualizationViewer;
import edu.uci.ics.jung.visualization.control.DefaultModalGraphMouse;
import edu.uci.ics.jung.visualization.control.GraphMouseListener;
import edu.uci.ics.jung.visualization.control.ModalGraphMouse;
import edu.uci.ics.jung.visualization.decorators.ToStringLabeller;
import edu.uci.ics.jung.visualization.renderers.Renderer.VertexLabel.Position;
import edu.uci.ics.jung.visualization.util.VertexShapeFactory;

import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.Frame;
import java.awt.Paint;
import java.awt.Shape;
import java.awt.Stroke;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;

import javax.swing.GroupLayout;
import javax.swing.GroupLayout.Alignment;
import javax.swing.LayoutStyle.ComponentPlacement;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.JTable;
import javax.swing.SpinnerNumberModel;
import javax.swing.border.EmptyBorder;

import org.apache.commons.collections15.Transformer;


public class drawNet2 {
	private JPanel contentPane;
	private JTable table;
	private JTable table2;
	private JFrame frame;
	private JFrame saveImageFrame;
	private JFrame saveTableFrame;
	private JFrame saveNetworkFrame;
	private JScrollPane scrollPane;
	private ArrayList<String> getgeneList = new ArrayList<String>();
	private ArrayList<String> getTFList = new ArrayList<String>();
	private ArrayList<String> listTF = new ArrayList<String>();
	private ArrayList<String> listTG = new ArrayList<String>();
	private ArrayList<String> listPG = new ArrayList<String>();
	private ArrayList<String> listvertex = new ArrayList<String>();
	private ArrayList<String> upgenes = new ArrayList<String>();
	private ArrayList<String> downgenes = new ArrayList<String>();
	private String[] Names = { "Pathway", "Num TG", "Num PG", "Num Overlap", "Pvalue", "Qvalue" };
	private String[] Names2 = { "TF", "Num target"};
	SparseMultigraph<String, String> initialgraph = null;
	SparseMultigraph<String, String> currentgraph = null;
	Collection<String> allVertex = null;
	Collection<String> allEdge = null;
	private ArrayList<String> vertexList = new ArrayList<String>();
	private ArrayList<String> edgeList = new ArrayList<String>();
		
	public void drawPic(List<List<String>> genetype, List<List<String>> networkfile, String[][] pathwaylist, 
			HashMap<String, List<String>> pathwaygene, int ndepth, TSMiner_DataSet rawDataSet) {
		// TODO Auto-generated method stub
		initialgraph = new SparseMultigraph<String, String>();
		currentgraph = new SparseMultigraph<String, String>();
		SparseMultigraph<String, String> graph = new SparseMultigraph<String, String>();
		
		KKLayout<String, String> layout = new KKLayout<String, String>(graph);

		HashMap<String, HashSet<String>> TFgene = new HashMap<String, HashSet<String>>();
		for(int i = 0; i < networkfile.size(); i++){
			String tf = networkfile.get(i).get(0);
			String tg = networkfile.get(i).get(1);
			if(TFgene.get(tf)!=null){
				HashSet<String> tglist = TFgene.get(tf);
				tglist.add(tg);
				TFgene.put(tf, tglist);
			}else{
				HashSet<String> tglist = new HashSet<String>();
				tglist.add(tg);
				TFgene.put(tf, tglist);
			}
		}
		
		String[][] TFlist = new String[TFgene.size()][2];
		int count=0;
		for(Entry<String, HashSet<String>> entry:TFgene.entrySet())
		  {
			String tf = entry.getKey();
			HashSet<String> tg = TFgene.get(tf);
			TFlist[count][0] = tf;
			TFlist[count][1] = tg.size()+"";
			count++;
		  }
		
		
		for(int i = 0; i < genetype.size(); i++){
			listvertex.add(genetype.get(i).get(0));
			if(genetype.get(i).get(1).equals("TF")){
				listTF.add(genetype.get(i).get(0));
			}else if(genetype.get(i).get(1).equals("TG")){
				listTG.add(genetype.get(i).get(0));
			}else if(genetype.get(i).get(1).equals("PG")){
				listPG.add(genetype.get(i).get(0));
			}
		}
		
		for(int i = 0; i < networkfile.size(); i++){
			initialgraph.addVertex(networkfile.get(i).get(0));
			initialgraph.addVertex(networkfile.get(i).get(1));
			initialgraph.addEdge(networkfile.get(i).get(0)+"-"+networkfile.get(i).get(1), networkfile.get(i).get(0), networkfile.get(i).get(1));
			currentgraph.addVertex(networkfile.get(i).get(0));
			currentgraph.addVertex(networkfile.get(i).get(1));
			currentgraph.addEdge(networkfile.get(i).get(0)+"-"+networkfile.get(i).get(1), networkfile.get(i).get(0), networkfile.get(i).get(1));
			graph.addVertex(networkfile.get(i).get(0));
			graph.addVertex(networkfile.get(i).get(1));
			graph.addEdge(networkfile.get(i).get(0)+"-"+networkfile.get(i).get(1), networkfile.get(i).get(0), networkfile.get(i).get(1));
		}
		
		allVertex = graph.getVertices();
		allEdge = graph.getEdges();
		for(String str : allEdge) {	
			edgeList.add(str);
		}
		for(String str : allVertex) {
			vertexList.add(str);
		}
		
		VisualizationViewer<String, String> vv = new VisualizationViewer<String, String>(layout);
		
		//初始化点颜色
		Transformer<String, Paint> vertexPaint = new Transformer<String, Paint>() {
			public Paint transform(String s) {
				if(listTF.contains(s)){
					return Color.GREEN;
				}else if(listTG.contains(s)){
					return Color.YELLOW;
				}else if(listPG.contains(s)){
					return Color.PINK;
				}else {					
					return null;
				}
			}
		};
		
		Transformer<String, Paint> changePaint = new Transformer<String, Paint>() {
			public Paint transform(String s) {
				if(getgeneList.contains(s)){
					return Color.RED;
				}
				else {					
					return Color.BLACK;
				}
			}
		};
		
		Transformer<String, Paint> DEchangePaint = new Transformer<String, Paint>() {
			public Paint transform(String s) {
				if(upgenes.contains(s)){
					return Color.RED;
				}
				else if(downgenes.contains(s)){					
					return Color.BLUE;
				}else{
					return Color.black;
				}
			}
		};
		
		Transformer<String, Stroke> DEchangeStroke = new Transformer<String, Stroke>() {
			public Stroke transform(String s) {
				if(upgenes.contains(s)){
					return new BasicStroke(2f);
				}else if(downgenes.contains(s)){	
					return new BasicStroke(2f);
				}else {					
					return null;
				}
			}
		};
		
		Transformer<String, Stroke> changeStroke1 = new Transformer<String, Stroke>() {
			public Stroke transform(String s) {
				if(getgeneList.contains(s)){
					return new BasicStroke(1f);
				}
				else {					
					return null;
				}
			}
		};
		
		Transformer<String, Stroke> changeStroke3 = new Transformer<String, Stroke>() {
			public Stroke transform(String s) {
				if(getgeneList.contains(s)){
					return new BasicStroke(2f);
				}
				else {					
					return null;
				}
			}
		};
		
		/*
		Transformer<String, Paint> vertexPaintSelected = new Transformer<String, Paint>() {
			public Paint transform(String s) {
				
				for(String str : allEdge) {	
					graph.removeEdge(str);
				}
				for(String str : allVertex) {
					graph.removeVertex(str);
				}
				
				for(int i = 0; i < getList.size(); i++){ //添加节点
					String vertex = getList.get(i);
					if(listTF.contains(vertex)){
						graph.addVertex(vertex);
					}else{
						graph.addVertex(vertex);
						Collection<String> tfs = initialgraph.getNeighbors(vertex);
						for(String str : tfs) {
							graph.addVertex(str);
							graph.addEdge(str+"-"+vertex, str, vertex);
						}
					}
				}

				vv.getRenderContext().setVertexFillPaintTransformer(vertexPaint);
				vv.getRenderContext().setVertexDrawPaintTransformer(changePaint);
				vv.getRenderContext().setVertexStrokeTransformer(changeStroke3);
				
				scrollPane.repaint();
				
				if(listTF.contains(s)){
					if(getList.contains(s))
						return Color.GREEN;
				}else if(listTG.contains(s)){
					if(getList.contains(s))
						return Color.YELLOW;
				}else if(listPG.contains(s)){
					if(getList.contains(s))
						return Color.ORANGE;
				}
				return Color.WHITE;								
			}
		};
		*/
		
		//顶点变成圆角矩形框
				Transformer<String, Integer> vst1 = new Transformer<String, Integer>() {
					public Integer transform(String s) {
						int len = s.length();
						if (len < 5)
							len = 5;
						return new Integer(len * 8 + 7);
					}
				};
				Transformer<String, Integer> vst2 = new Transformer<String, Integer>() {
					public Integer transform(String s) {
						int len = s.length();
						if (len < 5)
							len = 5;
						return new Integer(len * 9 + 7);
					}
				};
				Transformer<String, Float> vart = new Transformer<String, Float>() {
					public Float transform(String s) {
						//int len = s.length();
						//if (len < 3)
							//len = 3;
						return new Float(0.4); //高宽比
					}
				};
				VertexShapeFactory<String> vsf1 = new VertexShapeFactory<String>(vst1, vart);
				VertexShapeFactory<String> vsf2 = new VertexShapeFactory<String>(vst2, vart);
				Transformer<String, Shape> vstr = new Transformer<String, Shape>() {
					public Shape transform(String s) {
						if(listTF.contains(s)){
							return vsf2.getRegularStar(s, 8);
						}else if(listTG.contains(s)){
							return vsf1.getRoundRectangle(s);
						}else if(listPG.contains(s)){
							return vsf1.getRoundRectangle(s);
						}else {					
							return null;
						}
						
					}
				};
		
		//图形变换模式
		DefaultModalGraphMouse<Integer, String> gm = new DefaultModalGraphMouse<Integer, String>();		
		gm.setMode(ModalGraphMouse.Mode.PICKING);	
		//FRLayout<String, String> layout = new FRLayout<String, String>(graph);	
		vv.getRenderContext().setVertexFillPaintTransformer(vertexPaint); //fill color
		//vv.getRenderContext().setVertexDrawPaintTransformer(vertexPaint);
		vv.getRenderContext().setVertexLabelTransformer(new ToStringLabeller<String>()); 
		vv.getRenderContext().setVertexShapeTransformer(vstr); //shape
		vv.getRenderer().getVertexLabelRenderer().setPosition(Position.CNTR);
		vv.setGraphMouse(gm);
		vv.setBackground(Color.WHITE);
		
		
		//设置对话框属性
		frame = new JFrame("Regulatory Network");
		frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		frame.setBounds(100, 100, 1350, 700);
		
		contentPane = new JPanel();
		contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
		frame.setContentPane(contentPane);
		
		//添加table单击事件
		
		table = new JTable(pathwaylist, Names);			
		table.addMouseListener(new MouseAdapter() {
		     public void mouseClicked(MouseEvent e) { 
		    	 if(table.getValueAt(table.getSelectedRow(),0)!=null && currentgraph.getVertexCount()>0){		
		    		 getgeneList.clear();

		    		 
		    		 
		    		 for(int i = 0; i < pathwaygene.get(table.getValueAt(table.getSelectedRow(), 0)).size(); i++){
		    		    getgeneList.add(pathwaygene.get(table.getValueAt(table.getSelectedRow(), 0)).get(i));
		    		 }
		    		 
		    			for(String str : vertexList) {
		    				graph.removeVertex(str);
		    			}
		    			for(String str : edgeList) {
			                graph.removeEdge(str);
		    			}

						for(int i = 0; i < getgeneList.size(); i++){ //添加节点
							String vertex = getgeneList.get(i);
							if(listvertex.contains(vertex)){
								if(getTFList!=null && getTFList.size()>0){
									if(getTFList.contains(vertex)){
										if(listTF.contains(vertex)){
											graph.addVertex(vertex);
										}else{
											graph.addVertex(vertex);
											Collection<String> tfs = currentgraph.getNeighbors(vertex);
												if(tfs != null && tfs.size()>0){
													for(String str : tfs) {
														if(getTFList.contains(str)){
															graph.addVertex(str);
															graph.addEdge(str+"-"+vertex, str, vertex);
														}
													}
												}
										}
									}
								}else{
									if(listTF.contains(vertex)){
										graph.addVertex(vertex);
									}else{
										graph.addVertex(vertex);
										Collection<String> tfs = currentgraph.getNeighbors(vertex);
											if(tfs != null && tfs.size()>0){
												for(String str : tfs) {
													graph.addVertex(str);
													graph.addEdge(str+"-"+vertex, str, vertex);
												}
											}
									}
								}
							}
						}

						vv.getRenderContext().setVertexFillPaintTransformer(vertexPaint);
						vv.getRenderContext().setVertexDrawPaintTransformer(changePaint);
						vv.getRenderContext().setVertexStrokeTransformer(changeStroke3);
						scrollPane.repaint();
		    		    //vv.getRenderContext().setVertexFillPaintTransformer(vertexPaintSelected);
		    		    scrollPane.updateUI(); 
		    	 }
		     }
		});
		
		
		table2 = new JTable(TFlist, Names2);	
		table2.addMouseListener(new MouseAdapter() {
		     public void mouseClicked(MouseEvent e) { 
		    	 if(table2.getValueAt(table2.getSelectedRow(),0)!=null && currentgraph.getVertexCount()>0){		
		    		 getTFList.clear();
		    		 
		    		 for(String str : vertexList){
		    				graph.removeVertex(str);
		    		 }
		    		 for(String str : edgeList) {
			                graph.removeEdge(str);
		    		 }
		    		 
		    		 //List<String> tflist = new ArrayList<String>();
		    		 //tflist.add("NFATC1");
		    		 //tflist.add("STAT1");
		    		 //tflist.add("NFKB1");
		    		 //tflist.add("IRF9");
		    		 
		    			 //String TFID = tflist.get(i);
		    			 String TFID = (String)table2.getValueAt(table2.getSelectedRow(), 0);
			    		 getTFList.add(TFID);
			    		 for (String str : TFgene.get(TFID)) { 
			    			 getTFList.add(str);
			    		 }

			    			Collection<String> genes = currentgraph.getNeighbors(TFID);
			    			graph.addVertex(TFID);
			    			if(genes != null && genes.size()>0){
			    				for(String str : genes) {
			    					if(getgeneList != null && getgeneList.size()>0){
			    						if(getgeneList.contains(str)){
			    							graph.addVertex(str);
											graph.addEdge(TFID+"-"+str, TFID, str);
			    						}
			    					}else{
			    						graph.addVertex(str);
										graph.addEdge(TFID+"-"+str, TFID, str);
			    					}
			    				}
			    			}
		    		 
		    		

						vv.getRenderContext().setVertexFillPaintTransformer(vertexPaint);
						vv.getRenderContext().setVertexDrawPaintTransformer(changePaint);
						vv.getRenderContext().setVertexStrokeTransformer(changeStroke3);
						scrollPane.repaint();
		    		    //vv.getRenderContext().setVertexFillPaintTransformer(vertexPaintSelected);
		    		    scrollPane.updateUI(); 
		    	 }
		     }
		});
		
		
		Container c = frame.getContentPane();
		GroupLayout mainLayout = new GroupLayout(c);
		
		JScrollPane scrollPane1 = new JScrollPane(table);
		c.add(scrollPane1);
		//scrollPane1.setViewportView(table);
		
		JScrollPane scrollPane12 = new JScrollPane(table2);
		c.add(scrollPane12);
		//scrollPane12.setViewportView(table2)
		
		
		scrollPane = new JScrollPane(vv);
		scrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
		scrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED);
		//panel.add(scrollPane2); 
		c.add(scrollPane);
		
		
		
		JLabel fc1 = new JLabel("Differential expression at ");
		
		String[] labels = new String[rawDataSet.dsamplemins.length-ndepth];
		for(int i=ndepth;i<rawDataSet.dsamplemins.length;i++){
			labels[i-ndepth] = rawDataSet.dsamplemins[i];
		}
		JComboBox comboBox = new JComboBox(labels);
		comboBox.setMaximumSize(new Dimension(100,30));
		
		JLabel fc2 = new JLabel(" with threshold :");
		SpinnerNumberModel FCvalue = new SpinnerNumberModel(new Double(0.0), new Double(0.0), null, new Double(0.1));
		JSpinner FC2 = new JSpinner(FCvalue);
		FC2.setMaximumSize(new Dimension(120,30));
		JButton fcbutton = new JButton();
		fcbutton.setText("DE Genes");
		c.add(fcbutton);
		
		JLabel degree1 = new JLabel("Degree :");
		SpinnerNumberModel degreevalue = new SpinnerNumberModel(new Integer(1), new Integer(0), null, new Integer(1));
		JSpinner degree2 = new JSpinner(degreevalue);
		degree2.setMaximumSize(new Dimension(120,30));
		JButton filter = new JButton();
		filter.setText("Filter TFs");
		c.add(filter);
		
		JButton restore = new JButton();
		restore.setText("Restore");
		c.add(restore);
		
		JButton savemap = new JButton();
		savemap.setText("Save Image");
		c.add(savemap);
		JButton savepathwaytable = new JButton();
		savepathwaytable.setText("Save Pathway Table");
		c.add(savepathwaytable);
		JButton savetftable = new JButton();
		savetftable.setText("Save TF Table");
		c.add(savetftable);
		JButton savenetwork = new JButton();
		savenetwork.setText("Save Network File");
		c.add(savenetwork);
		
		//frame.pack();
		fcbutton.addActionListener(e -> {
			int time = comboBox.getSelectedIndex();
			double threshold = (Double) FC2.getValue();
			upgenes.clear();
			downgenes.clear();
			for(int i = 0; i < genetype.size(); i++){
				for(int j=0;j<rawDataSet.genenames.length;j++){
					if(genetype.get(i).get(0).equalsIgnoreCase(rawDataSet.genenames[j])){
						double fc = rawDataSet.data[j][ndepth+time]-rawDataSet.data[j][ndepth+time-1];
						if(fc>=threshold){ 
							upgenes.add(genetype.get(i).get(0));
						}else if(fc<=-threshold){  
							downgenes.add(genetype.get(i).get(0));
						}
						break;
					}
				}
			}
			//vv.getRenderContext().setVertexFillPaintTransformer(vertexPaint); //fill color
			vv.getRenderContext().setVertexDrawPaintTransformer(DEchangePaint);
			vv.getRenderContext().setVertexStrokeTransformer(DEchangeStroke);
			scrollPane.repaint();
        });
		
		filter.addActionListener(e -> {
			int threshold = (int) degree2.getValue();
			getgeneList.clear();
			getTFList.clear();
			
			for(String str : vertexList) {
				graph.removeVertex(str);
				currentgraph.removeVertex(str);
			}
			for(String str : edgeList) {
                graph.removeEdge(str);
                currentgraph.removeVertex(str);
			}
			
			for(String str : vertexList) {
				if(listTF.contains(str) && initialgraph.getNeighborCount(str)>=threshold){
					Collection<String> neighbors = initialgraph.getNeighbors(str);
					graph.addVertex(str);
					currentgraph.addVertex(str);
					for(String str2 : neighbors) {
						graph.addVertex(str2);
						graph.addEdge(str+"-"+str2, str, str2);
						currentgraph.addVertex(str2);
						currentgraph.addEdge(str+"-"+str2, str, str2);
					}
				}
			}
			vv.getRenderContext().setVertexFillPaintTransformer(vertexPaint); //fill color
			vv.getRenderContext().setVertexDrawPaintTransformer(changePaint);
			vv.getRenderContext().setVertexStrokeTransformer(changeStroke1);
			scrollPane.repaint();
        });
		
		restore.addActionListener(e -> {
			getgeneList.clear();
			getTFList.clear();
			for(int i = 0; i < networkfile.size(); i++){
				graph.addVertex(networkfile.get(i).get(0));
				graph.addVertex(networkfile.get(i).get(1));
				graph.addEdge(networkfile.get(i).get(0)+"-"+networkfile.get(i).get(1), networkfile.get(i).get(0), networkfile.get(i).get(1));
				currentgraph.addVertex(networkfile.get(i).get(0));
				currentgraph.addVertex(networkfile.get(i).get(1));
				currentgraph.addEdge(networkfile.get(i).get(0)+"-"+networkfile.get(i).get(1), networkfile.get(i).get(0), networkfile.get(i).get(1));
			}
			vv.getRenderContext().setVertexFillPaintTransformer(vertexPaint); //fill color
			vv.getRenderContext().setVertexDrawPaintTransformer(changePaint);
			vv.getRenderContext().setVertexStrokeTransformer(changeStroke1);
			scrollPane.repaint();
        });
		
		savemap.addActionListener(e -> {
			javax.swing.SwingUtilities.invokeLater(new Runnable() {
				public void run() {
					if (saveImageFrame == null) {
						saveImageFrame = new JFrame("Save as Image");
						saveImageFrame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
						saveImageFrame.setLocation(400,300);
						TSMinerGui_SaveImage newContentPane = new TSMinerGui_SaveImage(saveImageFrame,scrollPane);
						newContentPane.setOpaque(true);
						// content panes must be opaque
						saveImageFrame.setContentPane(newContentPane);
						// Display the window.
						saveImageFrame.pack();
					} else {
						saveImageFrame.setExtendedState(Frame.NORMAL);
					}
					saveImageFrame.setVisible(true);
				}
			});
			
        });
		
		savepathwaytable.addActionListener(e -> {
			javax.swing.SwingUtilities.invokeLater(new Runnable() {
				public void run() {
					if (saveTableFrame == null) {
						saveTableFrame = new JFrame("Save Table");
						saveTableFrame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
						saveTableFrame.setLocation(400,300);
						TSMinerGui_SaveTable newContentPane = new TSMinerGui_SaveTable(saveTableFrame,table);
						newContentPane.setOpaque(true);
						// content panes must be opaque
						saveTableFrame.setContentPane(newContentPane);
						// Display the window.
						saveTableFrame.pack();
					} else {
						saveTableFrame.setExtendedState(Frame.NORMAL);
					}
					saveTableFrame.setVisible(true);
				}
			});
        });
		
		savetftable.addActionListener(e -> {
			javax.swing.SwingUtilities.invokeLater(new Runnable() {
				public void run() {
					if (saveTableFrame == null) {
						saveTableFrame = new JFrame("Save Table");
						saveTableFrame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
						saveTableFrame.setLocation(400,300);
						TSMinerGui_SaveTable newContentPane = new TSMinerGui_SaveTable(saveTableFrame,table2);
						newContentPane.setOpaque(true);
						// content panes must be opaque
						saveTableFrame.setContentPane(newContentPane);
						// Display the window.
						saveTableFrame.pack();
					} else {
						saveTableFrame.setExtendedState(Frame.NORMAL);
					}
					saveTableFrame.setVisible(true);
				}
			});
        });
		
		savenetwork.addActionListener(e -> {
			javax.swing.SwingUtilities.invokeLater(new Runnable() {
				public void run() {
					if (saveNetworkFrame == null) {
						saveNetworkFrame = new JFrame("Save Network File");
						saveNetworkFrame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
						saveNetworkFrame.setLocation(400,300);
						TSMinerGui_SaveNetwork newContentPane = new TSMinerGui_SaveNetwork(saveNetworkFrame,networkfile);
						newContentPane.setOpaque(true);
						// content panes must be opaque
						saveNetworkFrame.setContentPane(newContentPane);
						// Display the window.
						saveNetworkFrame.pack();
					} else {
						saveNetworkFrame.setExtendedState(Frame.NORMAL);
					}
					saveNetworkFrame.setVisible(true);
				}
			});
        });
		
		GroupLayout.SequentialGroup  DEanalysis1 = mainLayout.createSequentialGroup().addGap(30).addComponent(fc1).addGap(3).addComponent(comboBox).addGap(3).addComponent(fc2).addGap(3).addComponent(FC2).addComponent(fcbutton);
		GroupLayout.SequentialGroup  TFfilter1 = mainLayout.createSequentialGroup().addComponent(degree1).addGap(3).addComponent(degree2).addComponent(filter);
		GroupLayout.SequentialGroup  hParalGroup01 = mainLayout.createSequentialGroup().addGroup(DEanalysis1).addGap(20).addGroup(TFfilter1).addGap(20).addComponent(restore).addGap(5);
		GroupLayout.SequentialGroup hParalGroup02 = mainLayout.createSequentialGroup().addComponent(savemap).addGap(10).addComponent(savepathwaytable).addGap(10).addComponent(savetftable).addGap(10).addGap(10).addComponent(savenetwork).addGap(5);
		GroupLayout.ParallelGroup hParalGroup03 = mainLayout.createParallelGroup(Alignment.TRAILING).addComponent(scrollPane, GroupLayout.DEFAULT_SIZE, 646, Short.MAX_VALUE)
				.addGroup(hParalGroup01);
		GroupLayout.ParallelGroup hParalGroup04 = mainLayout.createParallelGroup(Alignment.TRAILING).addComponent(scrollPane1, GroupLayout.DEFAULT_SIZE, 400, Short.MAX_VALUE)
				.addComponent(scrollPane12, GroupLayout.DEFAULT_SIZE, 400, Short.MAX_VALUE).addGroup(hParalGroup02);
		GroupLayout.SequentialGroup  hParalGroup = mainLayout.createSequentialGroup().addContainerGap().addGroup(hParalGroup03).addGap(15).addGroup(hParalGroup04).addContainerGap();
		mainLayout.setHorizontalGroup(hParalGroup);
		
		GroupLayout.ParallelGroup  DEanalysis2 = mainLayout.createParallelGroup().addComponent(fc1).addComponent(comboBox).addComponent(fc2).addComponent(FC2).addComponent(fcbutton);
		GroupLayout.ParallelGroup  TFfilter2 = mainLayout.createParallelGroup().addComponent(degree1).addComponent(degree2).addComponent(filter);
		GroupLayout.ParallelGroup hParalGroup05 = mainLayout.createParallelGroup().addGroup(DEanalysis2).addGroup(TFfilter2).addComponent(restore);
		GroupLayout.ParallelGroup hParalGroup06 = mainLayout.createParallelGroup().addComponent(savemap).addComponent(savepathwaytable).addComponent(savetftable).addComponent(savenetwork);
		GroupLayout.SequentialGroup hParalGroup07 = mainLayout.createSequentialGroup().addContainerGap().addComponent(scrollPane, GroupLayout.DEFAULT_SIZE, 646, Short.MAX_VALUE)
				.addGap(10).addGroup(hParalGroup05).addGap(10);
		GroupLayout.SequentialGroup hParalGroup08 = mainLayout.createSequentialGroup().addContainerGap().addComponent(scrollPane1, GroupLayout.DEFAULT_SIZE, 400, Short.MAX_VALUE)
				.addGap(10).addComponent(scrollPane12, GroupLayout.DEFAULT_SIZE, 400, Short.MAX_VALUE).addGap(10).addGroup(hParalGroup06).addGap(10);
		GroupLayout.ParallelGroup  hParalGroup2 = mainLayout.createParallelGroup().addGroup(hParalGroup07).addGap(18).addGroup(hParalGroup08);
		mainLayout.setVerticalGroup(hParalGroup2);

		//panel.setLayout(new BorderLayout(0, 0));
		contentPane.setLayout(mainLayout);
		frame.setVisible(true);
	}

}
