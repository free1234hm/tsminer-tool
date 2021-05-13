package bprc.nrem;

import bprc.core.*;

import java.awt.*;
import java.awt.event.*;
import java.io.*;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import javax.swing.event.*;

import java.util.*;
import java.text.*;
import java.net.*;

/**
 * Class implementing the main input interface
 */
public class TSMiner_IO extends JFrame implements ActionListener {
	/** IS TSMiner being run in a non-interactive, batch mode */
	boolean bbatchMode;
	/**
	 * The filename prefix used when saving the model and activites in batch
	 * mode
	 */
	String saveFile;

	String szorganismsourceval;
	String szxrefsourceval;
	boolean btraintest;
	public static final int GOANN = 0;
	public static final int XREF = 1;
	public static final int OBO = 2;
	public static final int MIRNAEXP = 3;
	public static final int FILETYPES = 4;
	boolean[] bdownloading = new boolean[FILETYPES];
	Object lockpd = new Object();
	Object lockxref = new Object();
	int[] npercentdone = new int[FILETYPES];
	boolean[] bexception = new boolean[FILETYPES];
	int nexceptions;
	static boolean bdisplaycurrent = false;
	static String EBIURL = "ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/";
	static String szDefaultFile = "";

	// Main
	static boolean ballowmergeDEF = false;
	static String szStaticFileDEF = "";
	static String szgoFileDEF = "";
	static String szDataFileDEF = "";
	static String szInitFileDEF = "";
	static boolean bspotcheckDEF = false;
	static int nnormalizeDEF = 1;
	static String szGeneAnnotationFileDEF = "";
	static String szCrossRefFileDEF = "";
	static int ndbDEF = 0;
	static int nstaticsourceDEF = 0;
	static int numchildDEF = 3;

	// Repeat
	static Vector<String> vRepeatFilesDEF = new Vector<String>();
	static boolean balltimeDEF = true;

	// Search Options
	static double dCONVERGENCEDEF = 0.01;
	static double dMinScoreDEF = 0.0;
	static double dDELAYPATHDEF = .15;
	static double dDMERGEPATHDEF = .15;
	static double dPRUNEPATHDEF = .15;
	static double dMINSTDDEVALDEF = 0.0;
	static double dPvalueshold = 0.01;
	static double dPvalueshold2 = 0.01;

	static double dMinScoreDIFFDEF = 0;
	static double dDELAYPATHDIFFDEF = 0;
	static double dDMERGEPATHDIFFDEF = 0;
	static double dPRUNEPATHDIFFDEF = 0;
	static int ninitsearchDEF = 0;
	static int npermutationDEF = 0;
	static int nSEEDDEF = 1260;
	static double dNODEPENALTYDEF = 40;
	static boolean bPENALIZEDDEF = true;
	static boolean bstaticsearchDEF = true;

	static String miRNAInteractionDataFile = "";
	static String miRNAExpressionDataFile = "";
	static boolean checkStatusTF = false;

	static boolean checkStatusmiRNA = false;
	static boolean miRNATakeLog = true;
	static boolean miRNAAddZero = false;
	static double miRNAWeight = 1.0;
	static double tfWeight = 0.5;
	static Vector<String> miRNARepeatFilesDEF = new Vector<String>();
	static boolean miRNAalltimeDEF = true;
	static boolean bgetmirnaexp = false;


	// Regulator Scoring
	static String regScoreFile = "";

	// Filtering
	static String szPrefilteredDEF = "";
	static int nMaxMissingDEF = 0;
	static double dMinExpressionDEF = 1;
	static double dMinCorrelationRepeatsDEF = 0;
	static boolean bfilterstaticDEF = false;

	// Gene Annotations
	static String sztaxonDEF = "";
	static String szevidenceDEF = "";
	static boolean bpontoDEF = true;
	static boolean bcontoDEF = true;
	static boolean bfontoDEF = true;
	static String szcategoryIDDEF = "";

	// GO Analysis
	static int nSamplesMultipleDEF = 500;
	static int nMinGoGenesDEF = 5;
	static int nMinGOLevelDEF = 3;
	static boolean brandomgoDEF = true;
	static boolean bmaxminDEF = false;
	static boolean bcompare0DEF = false;
	static boolean bcomparepreDEF = false;
	static double dpercentDEF = 0;

	// GUI
	static boolean brealXaxisDEF = false;
	static double dYaxisDEF = 1;
	static double dXaxisDEF = 1;
	static int nKeyInputTypeDEF = 1;
	static double dKeyInputXDEF = 3;
	static String SZSTATICDIR = "TF-gene file";
	static String SZGODIR = "Gene annotation file";
	static String szGeneOntologyFileDEF = "gene_ontology.obo";
	static double dnodekDEF = 1;
	
	// STSMiner
	static double dProbBindingFunctional = 0.8;

	long s1;
	JRadioButton bfButton;
	JRadioButton fdrButton;
	JRadioButton noneButton;
	JRadioButton bfgoButton;
	JRadioButton randomgoButton;
	JRadioButton initopenButton;
	JRadioButton initsearchButton;
	JRadioButton initnoButton;
	
	static boolean bendsearch = false;

	String szorig1val;
	String szstaticFileval;
	String szgoFileval;
	String szxrefval = szCrossRefFileDEF;
	String szorig2val;
	String szgoval = szGeneAnnotationFileDEF;
	String szgocategoryval = szGeneOntologyFileDEF;
	String szextraval;
	String szcategoryIDval;
	String szinitfileval;
	String sznumchildval;

	String szevidenceval;
	String sztaxonval;
	boolean ballowmergeval;
	boolean bpontoval;
	boolean bcontoval;
	boolean bfontoval;
	String szmaxmissingval;
	String szexpressval;
	String szfilterthresholdval;
	String szlbval;
	String szalphaval;
	String szpercentileval;
	String sznumberprofilesval;
	String szsamplepvalval;
	String szmingoval;
	String szmingolevelval;
	boolean bstaticcheckval;
	int ninitsearchval = 1;
	int npermutationval = 0;
	boolean brandomgoval = true;
	boolean btakelog = false;
	boolean bstaticsearchval;
	int ndb;
	int nxrefcb;
	int nstaticsourcecb;
	int ngosourcecb;
	String szusergann;
	String szuserFileField;
	String szuserFileField2;
	String szuserxref;
	String szepsilonval;
	String sznodepenaltyval;
	String szconvergenceval;
	String szminstddeval;
	static String enrichmentthreshold;
	String dethreshold;
	String szprunepathval;
	String szdelaypathval;
	String szmergepathval;
	String szprunepathvaldiff;
	String szdelaypathvaldiff;
	String szmergepathvaldiff;
	String szseedval=1280+"";
	boolean bmaxminval;
	boolean bfcto0;
	boolean bfctopre;

	boolean balltimemiRNA = false;
	boolean balltime = false;
	boolean bspotincluded;
	boolean badd0 = false;

	static JFileChooser theChooser = new JFileChooser();
	JLabel orig1Label;
	JLabel staticFileLabel;
	JLabel xrefLabel;
	JLabel xrefsourceLabel;
	JLabel categoryIDLabel;
	JLabel extraLabel;
	JLabel compare1Label;
	JLabel compare2Label;
	JLabel alphaLabel;
	JLabel corrmodelLabel;
	JLabel maxmissingLabel;
	JLabel filterthresholdLabel;
	JLabel expressLabel;
	JLabel maxchangeLabel;
	JLabel percentileLabel;
	JLabel lbLabel;
	JLabel samplepvalLabel;
	JLabel goLabel;
	JLabel mingoLabel;
	JLabel mingolevelLabel;
	JLabel randomgoLabel;
	JLabel epsilonLabel;
	JLabel prunepathLabel;
	JLabel delaypathLabel;
	JLabel mergeLabel;
	JLabel epsilonLabeldiff;
	JLabel prunepathLabeldiff;
	JLabel delaypathLabeldiff;
	JLabel mergeLabeldiff;
	JLabel numchildLabel;
	JLabel seedLabel;
	JLabel nodepenaltyLabel;
	JLabel initFileLabel;
	JLabel savedModelOptionsLabel;
	JLabel filterchoiceLabel;
	JLabel modelLabel;
	JLabel convergenceJLabel;
	JLabel minstddevJLabel;

	JComboBox staticsourcecb;
	JComboBox gosourcecb;

	JCheckBox anncheck = new JCheckBox("Annotations", false);
	JCheckBox xrefcheck = new JCheckBox("Cross References", false);
	JDialog executeDialognf = null;
	

	ListDialog miRNARepeatList;
	JButton miRNARepeatHButton = new JButton(Util.createImageIcon("Help16.gif"));

	ButtonGroup miRNAnormGroup = new ButtonGroup();
	JRadioButton miRNAlognormButton;
	JRadioButton miRNAnormButton;
	JRadioButton miRNAnonormButton;
	JButton miRNANormHButton = new JButton(Util.createImageIcon("Help16.gif"));

	static JCheckBox TfCheckBox;
	static JCheckBox miRNACheckBox;
	JButton miRNAScoringHButton = new JButton(Util
			.createImageIcon("Help16.gif"));

	JCheckBox filtermiRNABox;
	JButton miRNAFilterHButton = new JButton(Util.createImageIcon("Help16.gif"));

	JSpinner spinnerMIRNAWeight;
	JButton miRNAWeightHButton = new JButton(Util.createImageIcon("Help16.gif"));
	JSpinner spinnerTFWeight;
	JButton miRNATFWeightHButton = new JButton(Util
			.createImageIcon("Help16.gif"));
	JSpinner spinnerdProbBind;
	
	
	// DECOD gui
	JButton fastaDataFileButton = new JButton("Browse...", Util
			.createImageIcon("Open16.gif"));
	JButton fastaDataHButton = new JButton(Util.createImageIcon("Help16.gif"));
	JTextField fastaDataField;
	JButton decodPathButton = new JButton("Browse...", Util
			.createImageIcon("Open16.gif"));
	JButton decodPathHButton = new JButton(Util.createImageIcon("Help16.gif"));
	JTextField decodPathField;

	// Regulator Scoring GUI
	JButton regScoreFileButton = new JButton("Browse...", Util
			.createImageIcon("Open16.gif"));
	JButton regScoreHButton = new JButton(Util.createImageIcon("Help16.gif"));
	JTextField regScoreField;

	static int NUMCOLS = 42;

	// Strings for the labels
	static Color gray = new Color(235,235,235);
	static Color defaultColor;
	static String[] staticsourceArray = { "User provided" };
	static String[] gosourceArray = { "User provided" };


	JTextField goField;
	JTextField extraField;
	JTextField xrefField;
	JTextField categoryIDField;
	JTextField taxonField;
	JTextField evidenceField;


	JButton orig2Button = new JButton("Browse...", Util.createImageIcon("Open16.gif"));
	JButton miRNARepeatButton = new JButton("miRNA Repeat Data...", Util.createImageIcon("Open16.gif"));
	JButton infoButton = new JButton(Util.createImageIcon("About16.gif"));


	static JFileChooser fc = new JFileChooser();

	JSpinner thespinnermingo;
	JSpinner thespinnermingolevel;
	JSpinner thespinnermaxmissing;
	JSpinner thespinnerfilterthreshold;
	JSpinner thespinnerexpress;
	JSpinner thespinnersamplepval;

	JSpinner thespinnernumchild;
	JSpinner thespinnerconvergence;
	JSpinner thespinnerepsilondiff;

	JCheckBox thelogcheck;
	JDialog theOptions;
	ListDialog theRepeatList;
	String szClusterA = "Execute";
	
	JButton tfgenebutton= new JButton();
	JButton savedmodelbutton= new JButton();
	JButton timeseriesbutton= new JButton();
	JButton repeatButton= new JButton();
	JButton gobutton= new JButton();

	JButton help1= new JButton(Util.createImageIcon("Help16.gif"));
	JButton help2= new JButton(Util.createImageIcon("Help16.gif"));
	JButton help3= new JButton(Util.createImageIcon("Help16.gif"));

	JTextArea runtext;
	JButton Run = new JButton();
	JButton currentButton = new JButton();
	JButton endSearchButton = new JButton();
	
	JTextField staticFileField;
	JTextField initfileField;
	JTextField orig1Field;
	JTextField goFileField;
	
	String labels11[] = { "Log normalize data", "Normalize data", "No Normalization/add 0"};
	JComboBox comboBox3 = new JComboBox(labels11);
	JCheckBox jc = new JCheckBox("Filter gene if it has no binding TF");
	final JSpinner j13;
	final JSpinner j14;
	final JSpinner j15;
	
	JRadioButton mmm1=new JRadioButton("Maximum-minimum",false); 
    JRadioButton mmm2=new JRadioButton("Compare to 0",true);
    JRadioButton mmm3=new JRadioButton("Compare to previous",false);
	
    final JSpinner j19;
	final JSpinner j21;
	final JSpinner j22;
	final JSpinner j23;
	final JCheckBox jc3 = new JCheckBox("Use TF-TG interaction data to train the model",true);
	JRadioButton sss1;  
    JRadioButton sss2;  
    JRadioButton sss3;
    JRadioButton permute1;  
    JRadioButton permute2;  
	
	/**
	 * Class constructor - builds the input interface calls parseDefaults to get
	 * the initial settings from a default settings file if specified
	 */
	public TSMiner_IO() throws FileNotFoundException, IOException {
		super("TSMiner");

		File dir = new File(SZSTATICDIR);
		{
			String[] children = dir.list();
			if (children == null) {	
				/*
				final JFrame fframe = this;
				javax.swing.SwingUtilities.invokeLater(new Runnable() {
					public void run() {
						JOptionPane.showMessageDialog(fframe, "The directory "
								+ SZSTATICDIR + " was not found.",
								"Directory not found",
								JOptionPane.WARNING_MESSAGE);
					}
				});
				*/
			} else {
				staticsourceArray = new String[children.length + 1];
				staticsourceArray[0] = "User Provided";
				for (int i = 0; i < children.length; i++) {
					// Get filename of file or directory
					staticsourceArray[i + 1] = children[i];
				}
			}
		}
		staticsourcecb = new JComboBox(staticsourceArray);
		staticsourcecb.addActionListener(this);
		
		File dir2 = new File(SZGODIR);
		{
			String[] children = dir2.list();
			if (children != null) {
				gosourceArray = new String[children.length + 1];
				gosourceArray[0] = "User Provided";
				for (int i = 0; i < children.length; i++) {
					gosourceArray[i + 1] = children[i];
				}
			}
		}
		gosourcecb = new JComboBox(gosourceArray);
		gosourcecb.addActionListener(this);

		bbatchMode = false;
		saveFile = "";

		try {
			parseDefaults(); //Init Parameters
		} catch (IllegalArgumentException iex) {
			final IllegalArgumentException fiex = iex;
			final JFrame fframe = this;
			javax.swing.SwingUtilities.invokeLater(new Runnable() {
				public void run() {
					JOptionPane.showMessageDialog(fframe, fiex.getMessage(),
							"Exception thrown", JOptionPane.ERROR_MESSAGE);
				}
			});
		} catch (Exception ex) {
			final Exception fex = ex;
			final JFrame fframe = this;
			javax.swing.SwingUtilities.invokeLater(new Runnable() {
				public void run() {
					JOptionPane.showMessageDialog(fframe, fex.toString(),
							"Exception thrown", JOptionPane.ERROR_MESSAGE);
					fex.printStackTrace(System.out);
				}
			});
		}

		Container contentPane = getContentPane();
		BoxLayout layout = new BoxLayout(contentPane, BoxLayout.Y_AXIS);
		contentPane.setLayout(layout);

		JPanel panel1=new JPanel();
		panel1.setBorder(new TitledBorder(null,"Load Data",TitledBorder.LEFT,TitledBorder.TOP));
		panel1.setLayout(new BorderLayout());

		JLabel j1=new JLabel("TF-TG Interaction Source:");
		j1.setBounds(15,30,300,35);
		JLabel j2=new JLabel("Load TF-TG Interaction File:");
		j2.setBounds(15,62,300,35);
		JLabel j7=new JLabel("Gene Annotation Source:");
		j7.setBounds(15,94,200,35);
		JLabel j77=new JLabel("Load Gene Annotation File:");
		j77.setBounds(15,126,200,35);
		JLabel j3=new JLabel("Load Time-series File:");
		j3.setBounds(15,190,300,35);
		JLabel j33=new JLabel("Saved Model File (Optional):");
		j33.setBounds(15,158,300,35);
		JLabel j4=new JLabel("Normalization Method:");
		j4.setBounds(15,222,300,35);
		
		staticsourcecb.setBounds(200, 32, 250, 25);
        gosourcecb.setBounds(200, 96, 250, 25);
		comboBox3.setBounds(200,224,200,25);
        if (nnormalizeDEF == 0) {
        	comboBox3.setSelectedIndex(0);
		} else if (nnormalizeDEF == 1) {
			comboBox3.setSelectedIndex(1);
		} else {
			comboBox3.setSelectedIndex(2);
		}
		
		staticFileField=new JTextField(szStaticFileDEF, JLabel.TRAILING);
		staticFileField.setBounds(200,64,500,25);
		staticFileField.setBorder(BorderFactory.createLineBorder(Color.black));
		staticFileField.setOpaque(true);
		staticFileField.setBackground(Color.white);
		
		goFileField=new JTextField(szgoFileDEF, JLabel.TRAILING);
		goFileField.setBounds(200,128,500,25);
		goFileField.setBorder(BorderFactory.createLineBorder(Color.black));
		goFileField.setOpaque(true);
		goFileField.setBackground(Color.white);
		
		initfileField=new JTextField(szInitFileDEF, JLabel.TRAILING);
		initfileField.setBounds(200,160,500,25);
		initfileField.setBorder(BorderFactory.createLineBorder(Color.black));
		initfileField.setOpaque(true);
		initfileField.setBackground(Color.white);
		
		orig1Field=new JTextField(szDataFileDEF, JLabel.TRAILING);
		orig1Field.setBounds(200,192,500,25);
		orig1Field.setBorder(BorderFactory.createLineBorder(Color.black));
		orig1Field.setOpaque(true);
		orig1Field.setBackground(Color.white);

		tfgenebutton.setText("Load File");
		tfgenebutton.setHideActionText(false);
		tfgenebutton.setBounds(705,64,100,25);
		tfgenebutton.addActionListener(this);
		
		gobutton.setText("Load File");
		gobutton.setHideActionText(true);
		gobutton.addActionListener(this);
		gobutton.setBounds(705,128,100,25);
		
		savedmodelbutton.setText("Load File");
		savedmodelbutton.setHideActionText(true);
		savedmodelbutton.setBounds(705,160,100,25);
		savedmodelbutton.addActionListener(this);
		
		timeseriesbutton.setText("Load File");
		timeseriesbutton.setHideActionText(true);
		timeseriesbutton.addActionListener(this);
		timeseriesbutton.setBounds(705,192,100,25);

		repeatButton.setText("Repeat Data");
		repeatButton.setHideActionText(true);
		repeatButton.addActionListener(this);
		defaultColor = repeatButton.getBackground();
		if (TSMiner_IO.vRepeatFilesDEF.size() >= 1) {
			repeatButton.setBackground(ListDialog.buttonColor);
		}
		repeatButton.setBounds(810,192,120,25);
		
		help1.addActionListener(this);
		help1.setBounds(920,224,40,25);

		panel1.setLayout(null);
		panel1.add(j1); 
		panel1.add(j2); 
		panel1.add(j3); 
		panel1.add(j33);
		panel1.add(j4);
		panel1.add(staticFileField);
		panel1.add(initfileField);
		panel1.add(orig1Field);
		panel1.add(staticsourcecb);
		staticsourcecb.setSelectedIndex(nstaticsourceDEF);
		panel1.add(comboBox3);
		panel1.add(tfgenebutton);
		panel1.add(savedmodelbutton);
		panel1.add(timeseriesbutton);
		panel1.add(repeatButton);
		
		panel1.add(j7);
		panel1.add(gosourcecb);
		gosourcecb.setSelectedIndex(ndbDEF);
		panel1.add(j77);
		panel1.add(goFileField);
		panel1.add(gobutton);
		panel1.add(help1);

        JScrollPane   scrollpanel1   =   new   JScrollPane(panel1, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        scrollpanel1.setBounds(100, 100, 1040, 140);
		panel1.setPreferredSize(new Dimension(scrollpanel1.getWidth() - 50, scrollpanel1.getHeight()*2));
		/*******************************************************/
		JPanel panel3=new JPanel();
		JScrollPane   scrollpanel2   =   new   JScrollPane(panel3, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        scrollpanel2.setBounds(100, 100, 1040, 133);
        panel3.setPreferredSize(new Dimension(scrollpanel2.getWidth() - 50, scrollpanel2.getHeight()*2));
		panel3.setBorder(new TitledBorder(null,"Set Parameters",TitledBorder.LEFT,TitledBorder.TOP));
		panel3.setLayout(new BorderLayout());
		
		
		JPanel panel31=new JPanel();
		panel31.setBorder(new TitledBorder(null,"Filtering Options",TitledBorder.LEFT,TitledBorder.TOP));
		panel31.setBounds(5,25,480,235);
		
		jc.setBounds(10,15,300,35);
		if(bfilterstaticDEF) jc.setSelected(true);
		JLabel j9=new JLabel("Maximum number of missing values:");
		JLabel j10=new JLabel("Minimum correlation between repeats:");
		JLabel j11=new JLabel("Minimum absolute expression change:");
		JLabel j12=new JLabel("Change can be based on:");
		j9.setBounds(10,45,300,35);
		j10.setBounds(10,75,300,35);
		j11.setBounds(10,105,300,35);
		j12.setBounds(10,135,300,35);

		SpinnerNumberModel missingvalue = new SpinnerNumberModel(new Integer(nMaxMissingDEF), new Integer(0), null, new Integer(1));
		j13 = new JSpinner(missingvalue);	
		j13.setBounds(270,50,45,23);		
		SpinnerNumberModel corre = new SpinnerNumberModel(new Double(dMinCorrelationRepeatsDEF), new Double(0), null, new Double(0.05));
		j14 = new JSpinner(corre);	
		j14.setBounds(270,80,45,23);
		SpinnerNumberModel snmodelExpress = new SpinnerNumberModel(new Double(dMinExpressionDEF), new Double(0.0), null, new Double(0.1));
		j15=new JSpinner(snmodelExpress);
		j15.setBounds(270,110,45,23);
		mmm2.setBounds(10,165,120,32);
	    mmm3.setBounds(135,165,170,32);
	    mmm1.setBounds(305,165,155,32); 
	    if(szDefaultFile != null && szDefaultFile.length()>0){
	    	if(bmaxminDEF) {
		    	mmm1.setSelected(true);
		    }else{
		    	mmm1.setSelected(false);
		    }
		    if(bcompare0DEF) {
		    	mmm2.setSelected(true);
		    }else{
		    	mmm2.setSelected(false);
		    }
		    if(bcomparepreDEF) {
		    	mmm3.setSelected(true);
		    }else{
		    	mmm3.setSelected(false);
		    }
	    }
	    
		help2.addActionListener(this);
		help2.setBounds(430,200,40,25);
		
		panel31.setLayout(null);
		panel31.add(jc);
		panel31.add(j9);
		panel31.add(j10);
		panel31.add(j11);
		panel31.add(j12);
		panel31.add(j13);
		panel31.add(j14);
		panel31.add(j15);
		panel31.add(mmm1);
		panel31.add(mmm2);
		panel31.add(mmm3);
		panel31.add(help2);
		

		JPanel panel32=new JPanel();
		panel32.setBorder(new TitledBorder(null,"Model Options",TitledBorder.LEFT,TitledBorder.TOP));
		panel32.setBounds(490,25,480,235);
		
		JLabel maxchild=new JLabel("Maximum number of sub-paths:");
		maxchild.setBounds(10,15,300,35);
		JLabel minSD=new JLabel("Minimum standard deviation:");
		minSD.setBounds(10,45,300,35);
		
        jc3.setBounds(10,75,350,35);
        if(bstaticsearchDEF){
        	jc3.setSelected(true);
        }else{
        	jc3.setSelected(false);
        }
		
		JLabel ja1= new JLabel();
		ja1.setText("Saved model:");
		ja1.setBounds(10,105,110,35);
		
		sss1=new JRadioButton("As final");  
		sss1.setBounds(100,105,90,35);
	    sss2=new JRadioButton("Start from");  
	    sss2.setBounds(190,105,100,35);
	    sss3=new JRadioButton("Do not use");
	    sss3.setBounds(290,105,110,35);    
	    ButtonGroup g1=new ButtonGroup();
	    g1.add(sss1);
	    g1.add(sss2);
	    g1.add(sss3);
	
		if (ninitsearchDEF == 0) {
			sss1.setSelected(true);
		} else if (ninitsearchDEF == 1) {
			sss2.setSelected(true);
		} else {
			sss3.setSelected(true);
		}
		
		
		
		SpinnerNumberModel snnumchild = new SpinnerNumberModel(new Integer(numchildDEF), new Integer(2), null, new Integer(1));
	    j19 = new JSpinner(snnumchild);
		j19.setBounds(205,20,45,23);

		SpinnerNumberModel snminstddev = new SpinnerNumberModel(new Double(dMINSTDDEVALDEF), new Double(0), null, new Double(0.01));
		j21 = new JSpinner(snminstddev);
		j21.setBounds(205,50,45,23);
		
		JLabel threshold=new JLabel("Significance level of enrichment q-value:");
		threshold.setBounds(10,135,300,32);
		SpinnerNumberModel pthreshold = new SpinnerNumberModel(new Double(dPvalueshold), new Double(0), null, new Double(0.01));
		j22 = new JSpinner(pthreshold);
		j22.setBounds(245,140,45,23);
		
		JLabel threshold2=new JLabel("and DE q-value:");
		threshold2.setBounds(295,135,300,32);
		SpinnerNumberModel pthreshold2 = new SpinnerNumberModel(new Double(dPvalueshold2), new Double(0), null, new Double(0.01));
		j23 = new JSpinner(pthreshold2);
		j23.setBounds(385,140,45,23);
		
		 JLabel per1= new JLabel();
		 per1.setText("Permutation test:");
		 per1.setBounds(10,165,140,35);
		 
		 permute1=new JRadioButton("Compare to 0",true);  
		 permute1.setBounds(130,165,140,35);
		 permute2=new JRadioButton("Compare to previous",false);   
		 permute2.setBounds(270,165,160,35);
		 ButtonGroup g2=new ButtonGroup();
		 g2.add(permute1);
		 g2.add(permute2);
		 
		 if (npermutationDEF == 0) {
			 permute1.setSelected(true);
		 } else if (npermutationDEF == 1) {
			 permute2.setSelected(true);
		 }
		
		help3.addActionListener(this);
		help3.setBounds(430,200,40,25);
		
		panel32.setLayout(null);
		panel32.add(maxchild);
		panel32.add(minSD);
		panel32.add(threshold);
		panel32.add(threshold2);
		panel32.add(j19);
		panel32.add(j21);
		panel32.add(j22);
		panel32.add(j23);
		panel32.add(jc3);
		panel32.add(ja1);
		panel32.add(sss1);
		panel32.add(sss2);
		panel32.add(sss3);
		panel32.add(help3);
		
		panel32.add(per1);
		panel32.add(permute1);
		panel32.add(permute2);
		
		panel3.setLayout(null);
		panel3.add(panel31); 
		panel3.add(panel32); 

		/********************************************************/
		
		JPanel panel4=new JPanel();
		panel4.setPreferredSize(new Dimension(900, 180));
		panel4.setBorder(new TitledBorder(null,"Search Model",TitledBorder.LEFT,TitledBorder.TOP));
		panel4.setLayout(new BorderLayout());
		runtext=new JTextArea();
		runtext.setLineWrap(true);
		JScrollPane sp=new JScrollPane(runtext);
		sp.setBounds(10, 20, 500, 100);
		panel4.add(sp);
		
		Run.setText("Run");
		Run.setHideActionText(true);
		Run.addActionListener(this);	
		currentButton.setText("Show Current Model");
		currentButton.setEnabled(false);
		currentButton.addActionListener(this);
		endSearchButton.setText("End Searching");
		endSearchButton.setEnabled(false);
		endSearchButton.addActionListener(this);
		infoButton.addActionListener(this);
		infoButton.setBackground(gray);
		
		JPanel panel5=new JPanel();
		panel5.setPreferredSize(new Dimension(900, 50));
		panel5.add(Run);
		panel5.add(currentButton);
		panel5.add(endSearchButton);
		panel5.add(infoButton);
		
		contentPane.add(scrollpanel1);
		contentPane.add(scrollpanel2);
		contentPane.add(panel4);
		contentPane.add(panel5);

		theRepeatList = new ListDialog(this, TSMiner_IO.vRepeatFilesDEF,
				TSMiner_IO.balltimeDEF, repeatButton, TSMiner_IO.defaultColor,
				TSMiner_IO.defaultColor, TSMiner_IO.fc);
		miRNARepeatList = new ListDialog(this, TSMiner_IO.miRNARepeatFilesDEF,
				TSMiner_IO.miRNAalltimeDEF, miRNARepeatButton, TSMiner_IO.defaultColor,
				TSMiner_IO.defaultColor, TSMiner_IO.fc);
	}

	/**
	 * Assigns the initial settings of the parameters based on the contents of
	 * szDefaultFile
	 */
	public static void parseDefaults() throws FileNotFoundException,
			IOException {
		String szLine;
		BufferedReader br;
		try {
			String szError = "";
			br = new BufferedReader(new FileReader(szDefaultFile));
			while ((szLine = br.readLine()) != null) {
				StringTokenizer st = new StringTokenizer(szLine, "\t");
				if (st.hasMoreTokens()) {
					String sztype = st.nextToken().trim();
					String szvalue = "";
					if (st.hasMoreTokens()) {
						szvalue = st.nextToken().trim();
					}

					if (!szvalue.equals("")) {
						if ((sztype.equalsIgnoreCase("Use_static_input_to_build_model"))
								|| (sztype.equalsIgnoreCase("Use_TF_gene_interaction_data_to_build"))) {
							if (szvalue.equalsIgnoreCase("true")) {
								bstaticsearchDEF = true;
							} else if (szvalue.equalsIgnoreCase("false")) {
								bstaticsearchDEF = false;
							} else {
								szError += "Warning: " + szvalue
										+ " is an unrecognized " + "value for "
										+ sztype + " "
										+ "(expecting true or false)";
							}
						} else if (sztype.equalsIgnoreCase("Static_Input_Data_File")
								|| sztype.equalsIgnoreCase("TF-gene_Interaction_File")) {
							szStaticFileDEF = szvalue;
						} else if (sztype
								.equalsIgnoreCase("Gene_Annotation_File")) {
							szgoFileDEF = szvalue;
						} else if (sztype.equalsIgnoreCase("Data_File")
								|| sztype
										.equalsIgnoreCase("Expression_Data_File")) {
							szDataFileDEF = szvalue;
						} else if (sztype
								.equalsIgnoreCase("Minimum_Standard_Deviation")) {
							dMINSTDDEVALDEF = Double.parseDouble(szvalue);
						} else if (sztype.equalsIgnoreCase("Threshold_of_hypergeometric_test")) {
							dPvalueshold = Double.parseDouble(szvalue);
						} else if (sztype.equalsIgnoreCase("Threshold_of_permutation_test")) {
							dPvalueshold2 = Double.parseDouble(szvalue);
						} else if (sztype.equalsIgnoreCase("Saved_Model_File")) {
							szInitFileDEF = szvalue;
						} else if ((sztype
								.equalsIgnoreCase("Spot_IDs_included_in_the_data_file"))
								|| (sztype
										.equalsIgnoreCase("Spot_IDs_included_in_the_the_data_file"))
								|| (sztype
										.equalsIgnoreCase("Spot_IDs_in_the_data_file"))) {
							bspotcheckDEF = (szvalue.equalsIgnoreCase("true"));
						} else if ((sztype.equalsIgnoreCase("Normalize_Data"))
								|| (sztype.equalsIgnoreCase("Transform_Data"))
								|| (sztype.equalsIgnoreCase("Transform_Data[Log transform data,Linear transform data,Add 0]"))
								|| (sztype.equalsIgnoreCase("Normalize_Data[Log normalize data,Normalize data,No normalization/add 0]"))) {
							try {
								nnormalizeDEF = Integer.parseInt(szvalue);
								if ((nnormalizeDEF < 0) || (nnormalizeDEF > 2)) {
									throw new IllegalArgumentException(
											szvalue+ " is an invalid argument for Normalize_Data");
								}
							} catch (NumberFormatException ex) {
								if (szvalue.equalsIgnoreCase("Log normalize data")) {
									nnormalizeDEF = 0;
								} else if (szvalue.equalsIgnoreCase("Normalize data")) {
									nnormalizeDEF = 1;
								} else if (szvalue.equalsIgnoreCase("No normalization/add 0")) {
									nnormalizeDEF = 2;
								} else {
									throw new IllegalArgumentException(
											szvalue + " is an invalid argument for Normalize_Data");
								}
							}
						} else if ((sztype.equalsIgnoreCase("Change_can_be_based_on[Compare to 0,Compare to Previous,Maximum-Minimum]"))
								|| (sztype.equalsIgnoreCase("Change_can_be_based_on"))) {
							if (szvalue.contains("Maximum-Minimum")) {
								bmaxminDEF = true;
							}
							if (szvalue.contains("Compare to 0")) {
								bcompare0DEF = true;
							} 
							if (szvalue.contains("Compare to Previous")) {
								bcomparepreDEF = true;
							}
							if(!(szvalue.contains("Maximum-Minimum")||szvalue.contains("Compare to 0")||
									szvalue.contains("Compare to Previous"))){
								szError += szvalue
										+ " is an invalid value of "
										+ "Change_should_be_based_on[Maximum-Minimum,Difference From 0]\n";
							}
						} else if (sztype.equalsIgnoreCase("Gene_Annotation_Source")) {
							int numitems = gosourceArray.length;
							try {
								ndbDEF = Integer.parseInt(szvalue);
								if ((ndbDEF < 0) || (ndbDEF >= numitems)) {
									ndbDEF = 0;
								}
							} catch (NumberFormatException ex) {
								boolean bfound = false;
								int nsource = 0;
								while ((nsource < numitems)
										&& (!bfound)) {
									if (gosourceArray[nsource]
											.equalsIgnoreCase(szvalue)) {
										bfound = true;
										ndbDEF = nsource;
									} else {
										nsource++;
									}
								}

								if (!bfound) {
									szError += "Warning: "
											+ szvalue
											+ " is an unrecognized "
											+ "type for Gene Annotation Source\n";
								}
							}
						} else if (sztype.equalsIgnoreCase("TF-gene_Interaction_Source")) {
							int numitems = staticsourceArray.length;
							try {
								nstaticsourceDEF = Integer.parseInt(szvalue);
								if ((nstaticsourceDEF < 0) || (nstaticsourceDEF >= numitems)) {
									nstaticsourceDEF = 0;
								}
							} catch (NumberFormatException ex) {
								boolean bfound = false;
								int nsource = 0;
								while ((nsource < numitems) && (!bfound)) {
									if (((String) staticsourceArray[nsource])
											.equalsIgnoreCase(szvalue)) {
										bfound = true;
										nstaticsourceDEF = nsource;
									} else {
										nsource++;
									}
								}

								if (!bfound) {
									szError += "Warning: "
											+ szvalue
											+ " is an unrecognized "
											+ "type for TF-gene_Interaction_Source";
								}
							}
						} else if (sztype
								.equalsIgnoreCase("Gene_Annotation_File")) {
							szGeneAnnotationFileDEF = szvalue;
						} else if (sztype
								.equalsIgnoreCase("Cross_Reference_File")) {
							szCrossRefFileDEF = szvalue;
						} else if ((sztype
								.equalsIgnoreCase("Repeat_Data_Files(comma delimited list)"))
								|| (sztype
										.equalsIgnoreCase("Repeat_Data_Files"))) {
							vRepeatFilesDEF = new Vector<String>();
							StringTokenizer stRepeatList = new StringTokenizer(
									szvalue, ",");
							while (stRepeatList.hasMoreTokens()) {
								vRepeatFilesDEF.add(stRepeatList.nextToken());
							}
						} else if ((sztype
								.equalsIgnoreCase("Repeat_Data_is_from"))
								|| (sztype
										.equalsIgnoreCase("Repeat_Data_is_from[Different time periods,The same time period]"))) {
							if (szvalue
									.equalsIgnoreCase("Different time periods")) {
								balltimeDEF = true;
							} else if (szvalue
									.equalsIgnoreCase("The same time period")) {
								balltimeDEF = false;
							} else if (!szvalue.equals("")) {
								szError += "WARNING: '"
										+ szvalue
										+ "' is an invalid value for "
										+ "Repeat Data is from it must be either "
										+ "'Different time periods' or 'The same time period'\n";
							}
						} else if (sztype
								.equalsIgnoreCase("Y-axis_Scale_Factor")) {
							dYaxisDEF = Double.parseDouble(szvalue);
						} else if (sztype
								.equalsIgnoreCase("Scale_Node_Areas_By_The_Factor")) {
							dnodekDEF = Double.parseDouble(szvalue);
						} else if (sztype
								.equalsIgnoreCase("X-axis_Scale_Factor")) {
							dXaxisDEF = Double.parseDouble(szvalue);
						} else if ((sztype.equalsIgnoreCase("X-axis_scale"))
								|| (sztype
										.equalsIgnoreCase("X-axis_scale_should_be"))
								|| (sztype
										.equalsIgnoreCase("X-axis_scale[Uniform,Based on Real Time]"))
								|| (sztype
										.equalsIgnoreCase("X-axis_scale_should_be[Uniform,Based on Real Time]"))) {
							if (szvalue.equalsIgnoreCase("Uniform")) {
								brealXaxisDEF = false;
							} else if (szvalue
									.equalsIgnoreCase("Based on Real Time")) {
								brealXaxisDEF = true;
							} else if (!szvalue.equals("")) {
								szError += "WARNING: '"
										+ szvalue
										+ "' is an invalid value for "
										+ "X-axis_scale it must be either 'Uniform'"
										+ "or 'Based on Real Time'.\n";
							}
						} else if (sztype
								.equalsIgnoreCase("Key_Input_X_p-val_10^-X")) {
							dKeyInputXDEF = Double.parseDouble(szvalue);
						} else if ((sztype
								.equalsIgnoreCase("Key_Input_Significance_Based_On["
										+ "Path Significance Conditional on Split,Path Significance Overall,Split Significance]"))
								|| (sztype
										.equalsIgnoreCase("Key_Input_Significance_Based_On["
												+ "Split Significance,Path Significance Conditional on Split,Path Significance Overall]"))) {
							try {
								nKeyInputTypeDEF = Integer.parseInt(szvalue);
								if ((nKeyInputTypeDEF < 0)
										|| (nKeyInputTypeDEF > 2)) {
									throw new IllegalArgumentException(
											szvalue
													+ " is an invalid argument for Key Input Significance Based On");
								} else {
									// so code maps to input order
									if (nKeyInputTypeDEF == 0)
										nKeyInputTypeDEF = 1;
									else if (nKeyInputTypeDEF == 1)
										nKeyInputTypeDEF = 2;
									else if (nKeyInputTypeDEF == 2)
										nKeyInputTypeDEF = 0;

								}
							} catch (NumberFormatException ex) {
								if (szvalue
										.equalsIgnoreCase("Split Significance")) {
									nKeyInputTypeDEF = 0;
								} else if (szvalue
										.equalsIgnoreCase("Path Significance Conditional on Split")) {
									nKeyInputTypeDEF = 1;
								} else if (szvalue
										.equalsIgnoreCase("Path Significance Overall")) {
									nKeyInputTypeDEF = 2;
								} else {
									throw new IllegalArgumentException(
											szvalue
													+ " is an invalid argument for Key_Input_Significance_Based_On");
								}
							}
						} else if (sztype
								.equalsIgnoreCase("Maximum_number_of_paths_out_of_split")) {
							numchildDEF = Integer.parseInt(szvalue);
						} else if ((sztype.equalsIgnoreCase("Split_Seed"))||(sztype.equalsIgnoreCase("Random_Seed"))) {
							nSEEDDEF = Integer.parseInt(szvalue);
						} else if (sztype
								.equalsIgnoreCase("Penalized_likelihood_node_penalty")) {
							dNODEPENALTYDEF = Double.parseDouble(szvalue);
						} else if ((sztype.equalsIgnoreCase("Model_selection_framework"))
								|| (sztype.equalsIgnoreCase("Model_selection_framework[Penalized Likelihood,Train-Test]"))) {
							try {
								int ntempval;
								ntempval = Integer.parseInt(szvalue);
								if ((ntempval < 0) || (ntempval > 1)) {
									bPENALIZEDDEF = (ninitsearchDEF == 0);
									throw new IllegalArgumentException(
											szvalue
													+ " is an invalid argument for Model_selection_framework");
								}
							} catch (NumberFormatException ex) {
								if (szvalue
										.equalsIgnoreCase("Penalized Likelihood")) {
									bPENALIZEDDEF = true;
								} else if (szvalue
										.equalsIgnoreCase("Train-Test")) {
									bPENALIZEDDEF = false;
								} else if (!szvalue.equals("")) {
									szError += "WARNING: '"
											+ szvalue
											+ "' is an invalid value for "
											+ "Model_selection_framework "
											+ "it must be either "
											+ "'Use As Is', 'Start Search From', or 'Do Not Use'\n";
								}
							}
						}  else if ((sztype.equalsIgnoreCase("Prune_path_improvement"))||(sztype.equalsIgnoreCase("Delete_path_score_%"))) {
							dPRUNEPATHDEF = Double.parseDouble(szvalue);
							if (dPRUNEPATHDEF < 0) {
								throw new IllegalArgumentException(szvalue
										+ " is an invalid value for " + sztype
										+ " must be >= 0");
							}
						} else if ((sztype.equalsIgnoreCase("Minimum_score_improvement"))|| (sztype.equalsIgnoreCase("Main_search_score_%"))) {
							dMinScoreDEF = Double.parseDouble(szvalue);
							if (dMinScoreDEF < 0) {
								throw new IllegalArgumentException(szvalue
										+ " is an invalid value for " + sztype
										+ " must be >= 0");
							}
						} else if ((sztype.equalsIgnoreCase("Saved_Model"))|| (sztype.equalsIgnoreCase("Saved_Model[As Final/Start From/Do Not Use]"))) {
							try {
								ninitsearchDEF = Integer.parseInt(szvalue);
								if ((ninitsearchDEF < 0)|| (ninitsearchDEF > 2)) {
									throw new IllegalArgumentException(
											szvalue+ " is an invalid argument for Saved_Model");
								}
							} catch (NumberFormatException ex) {
								if (szvalue.equalsIgnoreCase("As Final")) {
									ninitsearchDEF = 0;
								} else if (szvalue
										.equalsIgnoreCase("Start From")) {
									ninitsearchDEF = 1;
								} else if (szvalue
										.equalsIgnoreCase("Do Not Use")) {
									ninitsearchDEF = 2;
								} else if (!szvalue.equals("")) {
									szError += "WARNING: '"
											+ szvalue
											+ "' is an invalid value for "
											+ "Saved_Model "
											+ "it must be either "
											+ "'Use As Is', 'Start Search From', or 'Do Not Use'\n";
								}
							}
						} else if ((sztype.equalsIgnoreCase("Permutation_Test"))|| (sztype.equalsIgnoreCase("Permutation_Test[Compare to 0/Compare to Previous]"))) {
							try {
								npermutationDEF = Integer.parseInt(szvalue);
								if ((npermutationDEF < 0)|| (npermutationDEF > 1)) {
									throw new IllegalArgumentException(
											szvalue+ " is an invalid argument for Permutation_Test");
								}
							} catch (NumberFormatException ex) {
								if (szvalue.equalsIgnoreCase("Compare to 0")) {
									npermutationDEF = 0;
								} else if (szvalue
										.equalsIgnoreCase("Compare to Previous")) {
									npermutationDEF = 1;
								} else if (!szvalue.equals("")) {
									szError += "WARNING: '"
											+ szvalue
											+ "' is an invalid value for "
											+ "Compare to 0 "
											+ "it must be either "
											+ "'Compare to 0', 'Compare to Previous'\n";
								}
							}
						} else if (sztype
								.equalsIgnoreCase("Filter_Gene_If_It_Has_No_Static_Input_Data")) {
							bfilterstaticDEF = (szvalue.equalsIgnoreCase("true"));
						} else if (sztype.equalsIgnoreCase("Maximum_Number_of_Missing_Values")) {
							nMaxMissingDEF = Integer.parseInt(szvalue);
						} else if (sztype.equalsIgnoreCase("Minimum_Absolute_Log_Ratio_Expression")) {
							dMinExpressionDEF = Double.parseDouble(szvalue);
						} else if (sztype
								.equalsIgnoreCase("Minimum_Correlation_between_Repeats")) {
							dMinCorrelationRepeatsDEF = Double
									.parseDouble(szvalue);
						} else if (sztype
								.equalsIgnoreCase("Pre-filtered_Gene_File")) {
							szPrefilteredDEF = szvalue;
						} else if (sztype
								.equalsIgnoreCase("Include_Biological_Process")) {
							bpontoDEF = (szvalue.equalsIgnoreCase("true"));
						} else if (sztype
								.equalsIgnoreCase("Include_Molecular_Function")) {
							bfontoDEF = (szvalue.equalsIgnoreCase("true"));
						} else if (sztype
								.equalsIgnoreCase("Include_Cellular_Process")) {
							bcontoDEF = (szvalue.equalsIgnoreCase("true"));
						} else if (sztype
								.equalsIgnoreCase("Only_include_annotations_with_these_evidence_codes")) {
							szevidenceDEF = szvalue;
						} else if (sztype
								.equalsIgnoreCase("Only_include_annotations_with_these_taxon_IDs")) {
							sztaxonDEF = szvalue;
						} else if ((sztype.equalsIgnoreCase("Category_ID_File"))
								|| (sztype
										.equalsIgnoreCase("Category_ID_Mapping_File"))) {
							szcategoryIDDEF = szvalue;
						} else if ((sztype
								.equalsIgnoreCase("GO_Minimum_number_of_genes"))
								|| (sztype
										.equalsIgnoreCase("Minimum_number_of_genes"))) {
							nMinGoGenesDEF = Integer.parseInt(szvalue);
						} else if (sztype.equalsIgnoreCase("Minimum_GO_level")) {
							nMinGOLevelDEF = Integer.parseInt(szvalue);
						} else if (sztype
								.equalsIgnoreCase("Minimum_Split_Percent")) {
							dpercentDEF = Double.parseDouble(szvalue);
						} else if (sztype
								.equalsIgnoreCase("Number_of_samples_for_randomized_multiple_hypothesis_correction")) {
							nSamplesMultipleDEF = Integer.parseInt(szvalue);
						} else if ((sztype.equalsIgnoreCase("Multiple_hypothesis_correction_method_enrichment[Bonferroni,Randomization]"))
								|| (sztype.equalsIgnoreCase("Multiple_hypothesis_correction_method[Bonferroni,Randomization]"))) {
							if (szvalue.equalsIgnoreCase("Bonferroni")) {
								brandomgoDEF = false;
							} else if (szvalue
									.equalsIgnoreCase("Randomization")) {
								brandomgoDEF = true;
							} else if (!szvalue.equals("")) {
								szError += "WARNING: '"
										+ szvalue
										+ "' is an invalid value for "
										+ "Correction_Method it must be either 'Bonferroni'"
										+ "or 'Randomization'.\n";
							}
						}
						else if (sztype.equalsIgnoreCase("miRNA-gene_Interaction_Source")) {
							miRNAInteractionDataFile = szvalue;
						} else if (sztype.equalsIgnoreCase("miRNA_Expression_Data_File")) {
							miRNAExpressionDataFile = szvalue;
						} else if (sztype.equalsIgnoreCase("Regulator_Types_Used_For_Activity_Scoring")) {
							if(szvalue.equalsIgnoreCase("None")){
								checkStatusTF = false;
								checkStatusmiRNA = false;
							} else if(szvalue.equalsIgnoreCase("TF")){
								checkStatusTF = true;
								checkStatusmiRNA = false;
							} else if (szvalue.equalsIgnoreCase("miRNA")) {
								checkStatusTF = false;
								checkStatusmiRNA = true;
							} else if (szvalue.equalsIgnoreCase("Both")) {
								checkStatusTF = true;
								checkStatusmiRNA = true;
							} else if (!szvalue.equals("")){
								szError = "WARNING: '"
									+ szvalue
									+ "' is an invalid value for "
									+ "Regulator_Types_Used_For_Activity_Scoring it must be either 'None', 'TF'"
									+ ", 'miRNA' or 'Both'.\n";
							}
						} else if (sztype.equalsIgnoreCase("Normalize_miRNA_Data[Log normalize data,Normalize data,No normalization/add 0]")) {
							if (szvalue.equalsIgnoreCase("Log normalize data")) {
								miRNATakeLog = true;
								miRNAAddZero = false;
							} else if (szvalue.equalsIgnoreCase("Normalize data")) {
								miRNATakeLog = false;
								miRNAAddZero = false;
							} else if (szvalue.equalsIgnoreCase("No normalization/add 0")) {
								miRNATakeLog = false;
								miRNAAddZero = true;
							} 
						} else if ((sztype.equalsIgnoreCase("Repeat_miRNA_Data_is_from"))
								|| (sztype
										.equalsIgnoreCase("Repeat_miRNA_Data_is_from[Different time periods,The same time period]"))) {
							if (szvalue.equalsIgnoreCase("Different time periods")) {
								miRNAalltimeDEF = true;
							} else if (szvalue.equalsIgnoreCase("The same time period")) {
								miRNAalltimeDEF = false;
							} else if (!szvalue.equals("")) {
								szError += "WARNING: '"
										+ szvalue
										+ "' is an invalid value for "
										+ "Repeat miRNA Data is from it must be either "
										+ "'Different time periods' or 'The same time period'\n";
							}
						} else if (sztype.equalsIgnoreCase("Expression_Scaling_Weight")) {
							miRNAWeight = Double.parseDouble(szvalue);
						} else if (sztype.equalsIgnoreCase("Minimum_TF_Expression_After_Scaling")) {
							tfWeight = Double.parseDouble(szvalue);
						} else if (sztype.equalsIgnoreCase("Regulator_Score_File")) {
							regScoreFile = szvalue;
						}  else if (sztype.equalsIgnoreCase("Active_TF_influence")) {
							dProbBindingFunctional = Double
									.parseDouble(szvalue);
							System.out
									.println("Setting active TF influence to "
											+ dProbBindingFunctional);
						} else if ((sztype.charAt(0) != '#')) {
							szError += "WARNING: '" + sztype
									+ "' is an unrecognized variable.\n";
						}
					}
				}
			}
			br.close();
			if (!szError.equals("")) {
				throw new IllegalArgumentException(szError);
			}
		} catch (FileNotFoundException ex) {
		}
	}

	/**
	 * Checks if the two data sets have the same number of rows, time points,
	 * and the gene name matches.
	 */
	public static void errorcheck(TSMiner_DataSet theDataSet1,
			TSMiner_DataSet theOtherSet) {

		if (theDataSet1.numcols != theOtherSet.numcols) {
			throw new IllegalArgumentException(
					"Repeat data set must have same "
							+ "number of columns as original, expecting "
							+ theDataSet1.numcols + " found "
							+ theOtherSet.numcols + " in the repeat");
		} else if (theDataSet1.numrows != theOtherSet.numrows) {
			throw new IllegalArgumentException(
					"Repeat data set must have same "
							+ "number of spots as the original, expecting "
							+ theDataSet1.numrows + " found "
							+ theOtherSet.numrows + " in the repeat");
		} else {
			for (int nrow = 0; nrow < theDataSet1.numrows; nrow++) {
				if (!theDataSet1.genenames[nrow]
						.equals(theOtherSet.genenames[nrow])) {
					throw new IllegalArgumentException("In row " + nrow
							+ " of the repeat set " + "expecting gene symbol "
							+ theDataSet1.genenames[nrow] + " found "
							+ theOtherSet.genenames[nrow]);
				} else if (!theDataSet1.probenames[nrow]
						.equals(theOtherSet.probenames[nrow])) {
					throw new IllegalArgumentException("In row " + nrow
							+ " of the repeat set " + "expecting gene symbol "
							+ theDataSet1.probenames[nrow] + " found "
							+ theOtherSet.probenames[nrow]);
				}
			}
		}
	}

	/**
	 * Checks if origcols and nrepeat cols are the same value, the length of
	 * origgenes and repeatgenes is the same, and the gene names are the same
	 */
	public static void errorcheck(String[] origgenes, String[] repeatgenes,
			int norigcols, int nrepeatcols) {
		if (norigcols != nrepeatcols) {
			throw new IllegalArgumentException(
					"Repeat data set must have same "
							+ "number of columns as original, expecting "
							+ norigcols + " found " + nrepeatcols
							+ " in the repeat");
		} else if (origgenes.length != repeatgenes.length) {
			throw new IllegalArgumentException(
					"Repeat data set must have same "
							+ "number of spots as the original, expecting "
							+ origgenes.length + " found " + repeatgenes.length
							+ " in the repeat");
		} else {
			for (int nrow = 0; nrow < origgenes.length; nrow++) {
				if (!origgenes[nrow].equals(repeatgenes[nrow])) {
					throw new IllegalArgumentException("In row " + nrow
							+ " of the repeat set " + "expecting gene symbol "
							+ origgenes[nrow] + " found " + repeatgenes[nrow]);
				} else if (!origgenes[nrow].equals(repeatgenes[nrow])) {
					throw new IllegalArgumentException("In row " + nrow
							+ " of the repeat set " + "expecting gene symbol "
							+ origgenes[nrow] + " found " + repeatgenes[nrow]);
				}
			}
		}
	}

	/**
	 * 利用输入文件和参数返回一个TSMiner_DataSet类型的数据集
	 */
	synchronized public static TSMiner_DataSet buildset(
			String szorganismsourceval, String szxrefsourceval,
			String szxrefval, String szexp1val, String szgoval,
			String szgocategoryval, int nmaxmissing, double dexpressedval,
			double dmincorrelation, int nsamplespval, int nmingo,
			int nmingolevel, String szextraval, boolean balltime,
			Vector<String> repeatnames, boolean btakelog,
			boolean bspotincluded, boolean badd0, String szcategoryIDval,
			String szevidenceval, String sztaxonval, boolean bpontoval,
			boolean bcontoval, boolean bfontoval, boolean brandomgoval,
			boolean bmaxminval, boolean bfcto0, boolean bfctopre) throws Exception {
		TSMiner_DataSet theDataSetsMerged = null;
		
		if (balltime) {
			TSMiner_DataSet theDataSet1 = new TSMiner_DataSet(szexp1val, nmaxmissing,
					dexpressedval, dmincorrelation, btakelog, bspotincluded,
					false, badd0, bmaxminval, bfcto0, bfctopre, balltime);

			if (theDataSet1.numcols <= 1) {
				theDataSet1 = new TSMiner_DataSet(theDataSet1.filterDuplicates(),
						new TSMiner_GoAnnotations(szorganismsourceval,
								szxrefsourceval, szxrefval, szgoval,
								szgocategoryval, theDataSet1.genenames,
								theDataSet1.probenames, nsamplespval, nmingo,
								nmingolevel, szextraval, szcategoryIDval,
								bspotincluded, szevidenceval, sztaxonval,
								bpontoval, bcontoval, bfontoval, brandomgoval));

				TSMiner_DataSet theDataSet1fm;
				if (theDataSet1.numcols == 1) {
					theDataSet1fm = new TSMiner_DataSet(theDataSet1.filterMissing1point(), theDataSet1.tga);
					theDataSet1fm = new TSMiner_DataSet(theDataSet1fm.filtergenesthreshold1point(), theDataSet1fm.tga);
				} else {theDataSet1fm = theDataSet1;
				}
				return theDataSet1fm;
			} else {
				String[] origgenes = theDataSet1.genenames;
				theDataSet1 = new TSMiner_DataSet(theDataSet1.logratio2(),theDataSet1.tga);
				theDataSet1 = new TSMiner_DataSet(theDataSet1.averageAndFilterDuplicates(), theDataSet1.tga);

				// genevalues in log ratio before averaging stored
				// need for each gene duplicated
				// a mutlidimensional array of time series for each occurence

				int numrepeats = repeatnames.size();

				if (numrepeats > 0) {
					TSMiner_DataSet[] repeatSets = new TSMiner_DataSet[numrepeats];
					for (int nset = 0; nset < numrepeats; nset++) {
						String szfile = (String) repeatnames.get(nset);

						TSMiner_DataSet theOtherSet = new TSMiner_DataSet(szfile,
								nmaxmissing, dexpressedval, dmincorrelation,
								btakelog, bspotincluded, true, badd0,
								bmaxminval, bfcto0, bfctopre, balltime);
						errorcheck(origgenes, theOtherSet.genenames,
								theDataSet1.numcols, theOtherSet.numcols);
						// compute log ratio of each time series first then
						// merge
						// normalize the data
						theOtherSet = new TSMiner_DataSet(theOtherSet.logratio2(),
								theOtherSet.tga);
						theOtherSet = new TSMiner_DataSet(theOtherSet
								.averageAndFilterDuplicates(), theOtherSet.tga);
						// gene values in log ratio before averaging stored
						repeatSets[nset] = theOtherSet;
					}
					theDataSetsMerged = new TSMiner_DataSet(theDataSet1.mergeDataSets(repeatSets), theDataSet1.tga);
					theDataSetsMerged = new TSMiner_DataSet(theDataSetsMerged.filterdistprofiles(theDataSet1, repeatSets),
							theDataSetsMerged.tga);
				} else {
					theDataSetsMerged = theDataSet1;
				}

				theDataSetsMerged = new TSMiner_DataSet(theDataSetsMerged.filterMissing(), theDataSetsMerged.tga);
				theDataSetsMerged = new TSMiner_DataSet(theDataSetsMerged.filtergenesthreshold2(), theDataSetsMerged.tga);

				theDataSetsMerged.tga = new TSMiner_GoAnnotations(
						szorganismsourceval, szxrefsourceval, szxrefval,
						szgoval, szgocategoryval, theDataSet1.genenames,
						theDataSet1.probenames, nsamplespval, nmingo,
						nmingolevel, szextraval, szcategoryIDval,
						bspotincluded, szevidenceval, sztaxonval, bpontoval,
						bcontoval, bfontoval, brandomgoval);

				theDataSetsMerged.addExtraToFilter(theDataSetsMerged.tga);
				return theDataSetsMerged;
			}
		} else {
			
			TSMiner_DataSet theDataSet1 = new TSMiner_DataSet(szexp1val, nmaxmissing,
					dexpressedval, dmincorrelation, btakelog, bspotincluded,
					false, badd0, bmaxminval, bfcto0, bfctopre, balltime);
			
			if (theDataSet1.numcols <= 1) {
				theDataSet1 = new TSMiner_DataSet(theDataSet1.filterDuplicates(),
						new TSMiner_GoAnnotations(szorganismsourceval,
								szxrefsourceval, szxrefval, szgoval,
								szgocategoryval, theDataSet1.genenames,
								theDataSet1.probenames, nsamplespval, nmingo,
								nmingolevel, szextraval, szcategoryIDval,
								bspotincluded, szevidenceval, sztaxonval,
								bpontoval, bcontoval, bfontoval, brandomgoval));
				TSMiner_DataSet theDataSet1fm;
				if (theDataSet1.numcols == 1) {
					theDataSet1fm = new TSMiner_DataSet(theDataSet1.filterMissing1point(), theDataSet1.tga);
					theDataSet1fm = new TSMiner_DataSet(theDataSet1fm.filtergenesthreshold1point(), theDataSet1fm.tga);
				} else {
					theDataSet1fm = theDataSet1;
				}
				return theDataSet1fm;
			} else {
				int numrepeats = repeatnames.size();

				if (numrepeats > 0) {
					TSMiner_DataSet[] repeatSets = new TSMiner_DataSet[numrepeats];
					for (int nset = 0; nset < numrepeats; nset++) {
						String szfile = (String) repeatnames.get(nset);

						TSMiner_DataSet theOtherSet = new TSMiner_DataSet(szfile,
								nmaxmissing, dexpressedval, dmincorrelation,
								btakelog, bspotincluded, true, badd0,
								bmaxminval, bfcto0, bfctopre, balltime);

						errorcheck(theDataSet1, theOtherSet);

						repeatSets[nset] = theOtherSet;
					}
					theDataSetsMerged = new TSMiner_DataSet(theDataSet1.mergeDataSets(repeatSets), theDataSet1.tga);
				} else {
					theDataSetsMerged = theDataSet1;
				}

				theDataSetsMerged = new TSMiner_DataSet(theDataSetsMerged.logratio2(), theDataSetsMerged.tga);
				theDataSetsMerged = new TSMiner_DataSet(theDataSetsMerged.averageAndFilterDuplicates(), theDataSetsMerged.tga);
				// gene values before averaging stored
				theDataSetsMerged = new TSMiner_DataSet(theDataSetsMerged.filterMissing(), theDataSetsMerged.tga);
				theDataSetsMerged = new TSMiner_DataSet(theDataSetsMerged.filtergenesthreshold2(), theDataSetsMerged.tga);

				theDataSetsMerged.tga = new TSMiner_GoAnnotations(
						szorganismsourceval, szxrefsourceval, szxrefval,
						szgoval, szgocategoryval, theDataSet1.genenames,
						theDataSet1.probenames, nsamplespval, nmingo,
						nmingolevel, szextraval, szcategoryIDval,
						bspotincluded, szevidenceval, sztaxonval, bpontoval,
						bcontoval, bfontoval, brandomgoval);
			}
			theDataSetsMerged.addExtraToFilter(theDataSetsMerged.tga);

			return theDataSetsMerged;
		}
	}

	/**
	 * Returns a TSMiner_DataSet based on the provided input parameters for
	 * miRNA.
	 */
	synchronized public static DataSetCore buildMIRNAset(
			String szorganismsourceval, String szxrefsourceval,
			String szxrefval, String szexp1val, String szgoval,
			String szgocategoryval, int nmaxmissing, double dexpressedval,
			double dmincorrelation, int nsamplespval, int nmingo,
			int nmingolevel, String szextraval, boolean balltime,
			Vector<String> repeatnames, boolean btakelog,
			boolean bspotincluded, boolean badd0, String szcategoryIDval,
			String szevidenceval, String sztaxonval, boolean bpontoval,
			boolean bcontoval, boolean bfontoval, boolean brandomgoval,
			boolean bmaxminval, boolean bfcto0, boolean bfctopre) throws Exception {
		DataSetCore theDataSetCoresMerged = null;

		TSMiner_DataSet theDataSet1 = new TSMiner_DataSet(szexp1val, nmaxmissing,
				dexpressedval, dmincorrelation, btakelog, bspotincluded,
				false, badd0, bmaxminval, bfcto0, bfctopre, balltime);
		DataSetCore theDataSet1Core = (DataSetCore) theDataSet1;
		if (theDataSet1.numcols <= 1) {
			theDataSet1Core = theDataSet1Core.filterDuplicates(); //去重

			DataSetCore theDataSet1fmCore;  
			if (theDataSet1Core.numcols == 1) {
				theDataSet1fmCore = theDataSet1Core.filterMissing1point();
				theDataSet1fmCore = theDataSet1fmCore.filtergenesthreshold1point();
			} else {
				theDataSet1fmCore = theDataSet1Core;
			}
			return theDataSet1fmCore;
		}
		
		if (balltime) {  //如果重复实验文件夹下全部文件被选中
			String[] origgenes = theDataSet1.genenames;
			theDataSet1Core = theDataSet1Core.logratio2();
			theDataSet1Core = theDataSet1Core.averageAndFilterDuplicates();
			// genevalues in log ratio before averaging stored
			// need for each gene duplicated
			// a mutlidimensional array of time series for each occurence
			int numrepeats = repeatnames.size();

			if (numrepeats > 0) {
				DataSetCore[] repeatSets = new DataSetCore[numrepeats];
				for (int nset = 0; nset < numrepeats; nset++) {
					String szfile = (String) repeatnames.get(nset);
					TSMiner_DataSet theOtherSet = new TSMiner_DataSet(szfile,
							nmaxmissing, dexpressedval, dmincorrelation,
							btakelog, bspotincluded, true, badd0,
							bmaxminval, bfcto0, bfctopre, balltime);
					errorcheck(origgenes, theOtherSet.genenames,
							theDataSet1.numcols, theOtherSet.numcols);
					// compute log ratio of each time series first then
					// merge
					// normalize the data
					DataSetCore theOtherSetCore = theOtherSet.logratio2();
					theOtherSetCore = theOtherSetCore.averageAndFilterDuplicates();
					// gene values in log ratio before averaging stored
					repeatSets[nset] = theOtherSetCore;
				}
				theDataSetCoresMerged = theDataSet1Core.mergeDataSets(repeatSets);
				theDataSetCoresMerged = theDataSetCoresMerged.filterdistprofiles(theDataSet1, repeatSets);
			} else {
				theDataSetCoresMerged = theDataSet1Core;
			}
			return theDataSetCoresMerged;
		} else {
			int numrepeats = repeatnames.size();

			if (numrepeats > 0) {
				DataSetCore[] repeatSets = new DataSetCore[numrepeats];
				for (int nset = 0; nset < numrepeats; nset++) {
					String szfile = (String) repeatnames.get(nset);
						TSMiner_DataSet theOtherSet = new TSMiner_DataSet(szfile,
							nmaxmissing, dexpressedval, dmincorrelation,
							btakelog, bspotincluded, true, badd0,
							bmaxminval, bfcto0, bfctopre, balltime);
					errorcheck(theDataSet1, theOtherSet);
					repeatSets[nset] = (DataSetCore)theOtherSet;
				}
				theDataSetCoresMerged = theDataSet1.mergeDataSets(repeatSets);
			} else {
				theDataSetCoresMerged = theDataSet1Core;
			}

			theDataSetCoresMerged = theDataSetCoresMerged.logratio2();
			theDataSetCoresMerged = theDataSetCoresMerged.averageAndFilterDuplicates();
			return theDataSetCoresMerged;
		}
	}
	
	/**
	 * A control method that handles the response for when the execute button on
	 * the interface is pressed including building the data set, running the
	 * TSMiner modeling procedure, and displaying the results
	 */
	public void clusterscript(
			String szstaticFileval, String szgoFileval, String szxrefval, String szexp1val, String szgoval,
			String szgocategoryval, String szmaxmissingval,
			String szexpressedval, String szfilterthresholdval,
			String szsamplepval, String szmingoval, String szmingolevelval,
			String szextraval, boolean balltime, Vector<String> repeatnames,
			Vector<String> mirnarepeatnames, boolean btakelog,
			boolean bgetxref, boolean bgetgoann, boolean bspotincluded,
			boolean badd0, String szcategoryIDval, String szinitfileval,
			String szevidenceval, String sztaxonval, boolean bpontoval,
			boolean bcontoval, boolean bfontoval, boolean brandomgoval,
			boolean bmaxminval, boolean bfcto0, boolean bfctopre) throws Exception {
		
		if (nexceptions == 0) {
			if (szstaticFileval.trim().equals("")) {
				throw new IllegalArgumentException(       
						"No transcription factor gene interaction input file given!");
			}
			if ((!szstaticFileval.trim().equals(""))
					&& (!(new File(szstaticFileval)).exists())) {
				throw new IllegalArgumentException("The transcription factor gene interaction input file '"
								+ szstaticFileval + "' cannot be found.");
			}

			if (szexp1val.trim().equals("")) {
				throw new IllegalArgumentException(
						"No time series input data file given!");
			} else if (!(new File(szexp1val)).exists()) {
				throw new IllegalArgumentException(
						"The time series input data file '" + szexp1val+ "' cannot be found.");
			}

			if (szinitfileval.trim().equals("")) {
				szinitfileval = "";
			} else if (!(new File(szinitfileval)).exists()) {
				throw new IllegalArgumentException("The initial model file '"+ szinitfileval + "' cannot be found.");
			}

			if (szcategoryIDval.trim().equals("")) {
				szcategoryIDval = "";
			} else if (!(new File(szcategoryIDval)).exists()) {
				throw new IllegalArgumentException("The category ID file '"+ szcategoryIDval + "' cannot be found.");
			}

			if (szxrefval.trim().equals("")) {
				szxrefval = "";
			} else if ((!bgetxref) && !(new File(szxrefval)).exists()) {
				throw new IllegalArgumentException("The cross reference file '"+ szxrefval + "' cannot be found.");
			}

			if (szgoval.trim().equals("")) {
				szgoval = "";
			} else if ((!bgetgoann) && (!(new File(szgoval)).exists())) {
				throw new IllegalArgumentException("The GO annotation file '"+ szgoval + "' cannot be found.");
			}

			if (szextraval.trim().equals("")) {
				szextraval = "";
			} else if (!(new File(szextraval)).exists()) {
				throw new IllegalArgumentException(
						"The pre-filtered gene list file '" + szextraval+ "' cannot be found.");
			}

			if (szgocategoryval.trim().equals("")) {
				szgocategoryval = "";
			}

			int nmaxmissing;
			try {
				nmaxmissing = Integer.parseInt(szmaxmissingval);
				if (nmaxmissing < 0) {
					throw new IllegalArgumentException(
							"Maximum missing values must be positive");
				}
			} catch (NumberFormatException ex) {
				throw new IllegalArgumentException(
						"Maximum missing values must be an integer");
			}

			for (int nrepeat = 0; nrepeat < repeatnames.size(); nrepeat++) {
				if (!(new File((String) repeatnames.get(nrepeat))).exists()) {
					throw new IllegalArgumentException("The repeat data file '"
							+ repeatnames.get(nrepeat) + "' cannot be found");
				}
			}

			for (int nrepeat = 0; nrepeat < mirnarepeatnames.size(); nrepeat++) {
				if (!(new File((String) mirnarepeatnames.get(nrepeat)))
						.exists()) {
					throw new IllegalArgumentException("The repeat data file '"
							+ mirnarepeatnames.get(nrepeat)
							+ "' cannot be found");
				}
			}

			double dmincorrelation = Double.parseDouble(szfilterthresholdval);
			if ((dmincorrelation < -1.1) || (dmincorrelation > 1.1)) {
				throw new IllegalArgumentException(
						"Correlation Lower Bound for Filtering must be in [-1.1,1.1]");
			}

			double dexpressedval = Double.parseDouble(szexpressedval);
			if (dexpressedval < -0.05) {
				throw new IllegalArgumentException(
						"Expression Value for filter must be >= -0.05");
			}

			int nmingo = Integer.parseInt(szmingoval);
			if (nmingo < 1) {
				throw new IllegalArgumentException(
						"Minimum number of GO genes must be at least 1");
			}

			int nmingolevel = Integer.parseInt(szmingolevelval);
			if (nmingolevel < 1) {
				throw new IllegalArgumentException(
						"Minimum number of GO level must be at least 1");
			}

			int nsamplespval;

			try {
				nsamplespval = Integer.parseInt(szsamplepval);
			} catch (NumberFormatException ex) {
				throw new IllegalArgumentException(
						"Number of samples for p-value correction must be an integer");
			}

			if (nsamplespval < 1) {
				throw new IllegalArgumentException(
						"Number of samples for p-value correction must be positive");
			}
			bendsearch = false;  //结束状态设置为不结束搜索
			
			
			
			
			final TSMiner_DataSet rawDataSet = TSMiner_IO.buildset(
					szorganismsourceval, szxrefsourceval, szxrefval, szexp1val,
					szgoval, szgocategoryval, 0, 0,
					0, nsamplespval, nmingo, nmingolevel,
					szextraval, balltime, repeatnames, btakelog, bspotincluded,
					badd0, szcategoryIDval, szevidenceval, sztaxonval,
					bpontoval, bcontoval, bfontoval, brandomgoval, bmaxminval, bfcto0, bfctopre);
			
			final TSMiner_DataSet thefDataSetfmnel = TSMiner_IO.buildset(  //利用buildset方法处理timeseries得到数据集
					szorganismsourceval, szxrefsourceval, szxrefval, szexp1val,
					szgoval, szgocategoryval, nmaxmissing, dexpressedval,
					dmincorrelation, nsamplespval, nmingo, nmingolevel,
					szextraval, balltime, repeatnames, btakelog, bspotincluded,
					badd0, szcategoryIDval, szevidenceval, sztaxonval,
					bpontoval, bcontoval, bfontoval, brandomgoval, bmaxminval, bfcto0, bfctopre);

			DataSetCore theMIRNADataSet = null;
			//---------------------Search the model-----------------------//			
			
			TSMiner_Timeiohmm thetimehmm = null;

			thetimehmm = new TSMiner_Timeiohmm(rawDataSet, thefDataSetfmnel, szstaticFileval, szgoFileval,
					sznumchildval, szseedval, bstaticcheckval, ballowmergeval, ninitsearchval,
					npermutationval, szinitfileval, runtext, endSearchButton,
					bstaticsearchval, brealXaxisDEF, dYaxisDEF, dXaxisDEF,
					dnodekDEF, nKeyInputTypeDEF, dKeyInputXDEF, dpercentDEF,
					currentButton, "", szminstddeval, enrichmentthreshold, dethreshold,
					staticsourceArray[nstaticsourcecb], checkStatusTF,
					checkStatusmiRNA, miRNAInteractionDataFile,
					theMIRNADataSet, regScoreFile, dProbBindingFunctional, miRNAWeight, tfWeight);
			
			thetimehmm.traverse(thetimehmm.treeptr, 0, true);
			thetimehmm.traverse(thetimehmm.treeptr, 0, false);

			final TSMiner_Timeiohmm fthetimehmm = thetimehmm;

			runtext.append("Display Interface..."+"\n");
			runtext.paintImmediately(runtext.getBounds());
			
			try {
				javax.swing.SwingUtilities.invokeAndWait(new Runnable() {
					public void run() {
						// not copying on last load
						final TSMinerGui theTSMinerGui = new TSMinerGui(fthetimehmm,
								fthetimehmm.treeptr, brealXaxisDEF, dYaxisDEF,
								dXaxisDEF, nKeyInputTypeDEF, dKeyInputXDEF,
								dpercentDEF, "(Final Model)", dnodekDEF, npermutationval);
						edu.umd.cs.piccolo.PCanvas.CURRENT_ZCANVAS = null;
						theTSMinerGui.setLocation(15, 40);
						theTSMinerGui.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

						theTSMinerGui.setVisible(true);
						theTSMinerGui.addWindowListener(new WindowAdapter() {
							public void windowClosing(WindowEvent we) {
								theTSMinerGui.closeWindows();
								if (!theTSMinerGui.bsavedchange) {
									Object[] options = { "Yes", "No" };
									int noption = JOptionPane.showOptionDialog(
													theTSMinerGui,
													"Would you like to save the model?",
													"Question",
													JOptionPane.YES_NO_OPTION,
													JOptionPane.QUESTION_MESSAGE,
													null, options, options[1]);

									if (noption == 0) {
										if (theTSMinerGui.saveModelFrame == null) {
											theTSMinerGui.saveModelFrame = new JFrame(
													"Save Model to File");
											theTSMinerGui.saveModelFrame
													.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
											theTSMinerGui.saveModelFrame
													.setLocation(400, 300);
											TSMinerGui_SaveModel newContentPane = new TSMinerGui_SaveModel(
													theTSMinerGui.theTimeiohmm,
													theTSMinerGui.theTimeiohmm.treeptr,
													theTSMinerGui.saveModelFrame,
													theTSMinerGui);
											newContentPane.setOpaque(true);
											// content panes must be opaque
											theTSMinerGui.saveModelFrame
													.setContentPane(newContentPane);
											// Display the window.
											theTSMinerGui.saveModelFrame.pack();
										} else {
											theTSMinerGui.saveModelFrame
													.setExtendedState(Frame.NORMAL);
										}
										theTSMinerGui.saveModelFrame
												.setVisible(true);
									}
								}
							}
						});

						if (bbatchMode && theTSMinerGui.theTimeiohmm.bindingData.regPriors != null) {
							theTSMinerGui.batchSave(saveFile);
							theTSMinerGui.dispose();
						}
					}
				});
			} catch (InterruptedException iex) {
				runtext.append("visulization error"+"\n");
				runtext.paintImmediately(runtext.getBounds());
			} catch (java.lang.reflect.InvocationTargetException itex) {
				runtext.append("visulization error"+"\n");
				runtext.paintImmediately(runtext.getBounds());
			}

			if (executeDialognf != null) {
				executeDialognf.dispose();
				executeDialognf.setVisible(false);
			}
		}
		long e1 = System.currentTimeMillis();
		System.out.println("Time: " + (e1 - s1) + "ms");

		if (bbatchMode) {
			this.dispose();
		}
	}

	/**
	 * update the input static data
	 */
	public void handlestaticsource() {
		if (staticFileField.isEditable()) {
			szuserFileField = staticFileField.getText(); 
		}

		nstaticsourcecb = staticsourcecb.getSelectedIndex(); 

		if (nstaticsourcecb >= 1) {
			staticFileField.setText(SZSTATICDIR
					+ System.getProperty("file.separator")
					+ staticsourceArray[nstaticsourcecb]);
			staticFileField.setEditable(false);
			tfgenebutton.setEnabled(false);
		} else {
			staticFileField.setText(szuserFileField);
			staticFileField.setEditable(true);
			tfgenebutton.setEnabled(true);
		}
	}
	public void handlestaticsource2() {

		if (goFileField.isEditable()) {
			szuserFileField2 = goFileField.getText(); 
		}
 
		ngosourcecb = gosourcecb.getSelectedIndex(); 

		if (ngosourcecb >= 1) {
			goFileField.setText(SZGODIR
					+ System.getProperty("file.separator")
					+ gosourceArray[ngosourcecb]);
			goFileField.setEditable(false);
			gobutton.setEnabled(false);
		} else {
			goFileField.setText(szuserFileField2);
			goFileField.setEditable(true);
			gobutton.setEnabled(true);
		}
	}


	/**
	 * define the button methods
	 */
	public void actionPerformed(ActionEvent e) {
		Object esource = e.getSource();

		if (esource == endSearchButton) {
			bendsearch = true;
			runtext.append("End Search Requested. Search Will End Soon..."+"\n");
			runtext.paintImmediately(runtext.getBounds());
			endSearchButton.setEnabled(false);
		} else if (esource == currentButton) {
			bdisplaycurrent = true;
			currentButton.setEnabled(false);
		} else if (esource == staticsourcecb) {
			handlestaticsource();
		} else if (esource == gosourcecb) {
			handlestaticsource2();
		} else if (esource == tfgenebutton) {
			int returnVal = fc.showOpenDialog(this);
			if (returnVal == JFileChooser.APPROVE_OPTION) {
				File file = fc.getSelectedFile();
				staticFileField.setText(file.getAbsolutePath());
			}
		} else if (esource == gobutton) {
			int returnVal = fc.showOpenDialog(this);
			if (returnVal == JFileChooser.APPROVE_OPTION) {
				File file = fc.getSelectedFile();
				goFileField.setText(file.getAbsolutePath());
			}
		} else if (esource == savedmodelbutton) {
			int returnVal = fc.showOpenDialog(this);

			if (returnVal == JFileChooser.APPROVE_OPTION) {
				File file = fc.getSelectedFile();
				initfileField.setText(file.getAbsolutePath());
			}
		} else if (esource == timeseriesbutton) {
			int returnVal = fc.showOpenDialog(this);

			if (returnVal == JFileChooser.APPROVE_OPTION) {
				File file = fc.getSelectedFile();
				orig1Field.setText(file.getAbsolutePath());
			}
		} else if (esource == repeatButton) {
			theRepeatList.setLocation(this.getX() + 75, this.getY() + 100);
			theRepeatList.setVisible(true);
		}  else if (esource == miRNARepeatButton) {
			miRNARepeatList.setLocation(this.getX() + 75, this.getY() + 100);
			// setModal forces the miRNA repeat dialog to appear as the top
			// window
			// and blocks focus from other windows
			miRNARepeatList.setModal(true);
			miRNARepeatList.setVisible(true);
		}  else if (esource == fastaDataFileButton) {
			int returnVal = fc.showOpenDialog(this);

			if (returnVal == JFileChooser.APPROVE_OPTION) {
				File file = fc.getSelectedFile();
				fastaDataField.setText(file.getAbsolutePath());
			}
		} else if (esource == decodPathButton) {
			int returnVal = fc.showOpenDialog(this);

			if (returnVal == JFileChooser.APPROVE_OPTION) {
				File file = fc.getSelectedFile();
				decodPathField.setText(file.getAbsolutePath());
			}
		} else if (esource == regScoreFileButton) {
			int returnVal = fc.showOpenDialog(this);

			if (returnVal == JFileChooser.APPROVE_OPTION) {
				File file = fc.getSelectedFile();
				regScoreField.setText(file.getAbsolutePath());
			}
		}
		
		else if (esource == Run) {
			s1 = System.currentTimeMillis();
			szorig1val = orig1Field.getText(); 
			szstaticFileval = staticFileField.getText();  
			szgoFileval = goFileField.getText();
			szinitfileval = initfileField.getText();
			
			int temp=comboBox3.getSelectedIndex();
			if(0==temp){
				btakelog=true;
			}else{
				btakelog=false;
			}
			if(2==temp){
				badd0=true;
			}else{
				badd0=false;
			}
			
			
			bstaticcheckval=false;
			if(jc.isSelected())
			{
				bstaticcheckval=true;
			}
			
			szmaxmissingval=j13.getValue().toString();
			szfilterthresholdval=j14.getValue().toString();
			szexpressval=j15.getValue().toString();

			if (mmm1.isSelected()) {
				bmaxminval=true;            
			}else{
				bmaxminval=false;
			}
			if (mmm2.isSelected()) { 
				bfcto0 = true;              
			}else{
				bfcto0 = false;
			}
			if(mmm3.isSelected()){
				bfctopre = true;   
			}else{
				bfctopre = false;
			}
			
			sznumchildval=j19.getValue().toString();
			szminstddeval=j21.getValue().toString();
			enrichmentthreshold = j22.getValue().toString();
			dethreshold = j23.getValue().toString();
			
			if(jc3.isSelected())
			{
				bstaticsearchval=true;
			}
			
			if (sss1.isSelected()) {
				ninitsearchval = 0; 
			} else if (sss2.isSelected()) {
				ninitsearchval = 1; 
			} else if(sss3.isSelected()){
				ninitsearchval = 2; 
			}
			if (permute1.isSelected()) {
				npermutationval = 0;        
			} else if (permute2.isSelected()) { 
				npermutationval = 1;              
			}

			szxrefval = ""; 
			szorig2val = "";
			szgoval = ""; 
			szgocategoryval = "";
			szextraval = "";
			szcategoryIDval = ""; 
			szsamplepvalval = 1+"";
			szmingoval = 1+""; 
			szmingolevelval = 1+""; 
			bspotincluded = false;
			
			balltime = theRepeatList.allButton.isSelected();
			balltimemiRNA = miRNARepeatList.allButton.isSelected();
			
			this.setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));

			nexceptions = 0;
			final boolean bgetxref = xrefcheck.isSelected();
			final boolean bgetgoann = anncheck.isSelected();

			final JFrame fframe = this;

			Runnable clusterrun = new Runnable() {
				public void run() {
					Run.setEnabled(false);
					try {
						clusterscript(szstaticFileval, szgoFileval, szxrefval, szorig1val,
								szgoval, szgocategoryval, szmaxmissingval,
								szexpressval, szfilterthresholdval,
								szsamplepvalval, szmingoval, szmingolevelval,
								szextraval, balltime, theRepeatList.data,
								miRNARepeatList.data, btakelog, bgetxref,
								bgetgoann, bspotincluded, badd0,
								szcategoryIDval, szinitfileval, szevidenceval,
								sztaxonval, bpontoval, bcontoval, bfontoval,
								brandomgoval, bmaxminval, bfcto0, bfctopre);
					} catch (IllegalArgumentException iex) {
						final IllegalArgumentException fiex = iex;
						iex.printStackTrace(System.out);

						if (executeDialognf != null) {
							executeDialognf.setVisible(false);
							executeDialognf.dispose();
						}

						javax.swing.SwingUtilities.invokeLater(new Runnable() {
							public void run() {
								JOptionPane.showMessageDialog(fframe, fiex.getMessage(), "Error",JOptionPane.ERROR_MESSAGE);
							}
						});
					} catch (Exception ex) {
						final Exception fex = ex;

						if (executeDialognf != null) {
							executeDialognf.setVisible(false);
							executeDialognf.dispose();
						}

						javax.swing.SwingUtilities.invokeLater(new Runnable() {
							public void run() {
								JOptionPane.showMessageDialog(fframe, fex.toString(), "Exception thrown",JOptionPane.ERROR_MESSAGE);
								fex.printStackTrace(System.out);
							}
						});
					}
					Run.setEnabled(true); 
					fframe.setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
				}
			};
			(new Thread(clusterrun)).start();
			
		} else if ((esource == help1) || (esource == help2) || (esource == help3)){
			makeHelpDialog((Component) esource);
		} else if (esource == infoButton) {//The Information and Help button
			String szMessage = "This is version 1.0.0 of TSMiner.\n\n"
					+ "The TSMiner was developed by Mingfei Han.\n"
					+ "The TSMiner is available under a GPL v3.0 license.\n"
					+ "Any questions or bugs found should "
					+ "be emailed to free1234hm@gmail.com, or to free1234hm@163.com.";

			Util.renderDialog(this, szMessage, 50, 100, "Information");
		}
	}

	/**
	 * Help button
	 */
	public void makeHelpDialog(Component esource) {
		String szMessage = "";
		if (esource == help1) {
			szMessage = "#  TF-gene Interaction Source:\n"
		            + "Selecting a file from the menu enters that file name in the 'Load TF-TG Interaction File' field.\n"
					+ "Selecting 'User Provided' leaves the 'Load TF-TG Interaction File' field empty so "
					+ "that a user can specify the desired interaction file. "
					+ "Consult Appendix A of the user manual for the listed TF-gene files."+"\n"
					+"=============================================================================================="+"\n"
					+ "#  Load TF-TG Interaction File:\n"
					+ "File encodes the static TF-gene interaction predictions as "
					+ "input to TSMiner. The file is a tab delmited text file, that can be in one of "
					+ "two formats. Either in a grid format or in a three column format.\n"
					+ "The grid format file is as follows:\n"
					+ "The first column of the file contains the gene identifier symbols. "
					+ "The first row contains the TFs identifiers. "
					+ "The first entry of the first row is a label for the gene symbol column. "
					+ "Each remaining entry corresponds to the relationship between a TF and a "
					+ "gene. Under a binary encoding an entry is 1 if the TF is predicted "
					+ "to regulate the gene and 0 otherwise. If a three way encoding is used an entry is 1 if "
					+ "the TF is predicted to activate the gene, -1 if it is predicted to repress "
					+ "the gene and 0 otherwise. Below is portion of a sample TF-gene interaction input file\n\n"
					+ "ID	ADR1	ARG80	ARG81	ARO80	BAS1	CAD1	CBF1\n"
					+ "YAL053W	0	0	0	0	0	0	1\n"
					+ "YAL054C	0	0	0	0	0	0	1\n"
					+ "YAL055W	0	0	0	0	1	0	0\n\n"
					+ "In the three column format the first column contains "
					+ "the transcription factors, the second column "
					+ "the regulated gene, and the third column input value. The first row is a header row "
					+ "where the header of the first column must be 'TF' column, "
					+ "and the second column must have the header 'Gene'. "
					+ "If a TF-gene pair is not present "
					+ "the input value is assumed to be 0.\n"
					+ "TF	Gene	Input\n"
					+ "BAS1	YAL055W	1\n"
					+ "CBF1	YAL053W	1\n" + "CBF1	YAL054C	1\n"
					+"=============================================================================================="+"\n"
					+ "#  'Gene Annotation Source' and 'Load Gene Annotation File':\n"
					+ "Used to input the files containing pathway gene sets.\n"
					+ "These two options and the file formats are same as the TF-gene inteaction files.\n"
					+ "Consult Appendix B of the user manual for the listed gene annotation files.\n"
                    +"=============================================================================================="+"\n"
                    + "#  Saved Model File(Optional):\n"
                    + "This file is optional that inputs a saved model produced before. "
					+ "A resulting model can be saved by pressing the button 'Save Model' in the 'Display Interface'.\n"
					+ "Depending on the setting of 'Saved Model' in "
					+ "'Model Options', this model is either opened as final model, a search is started from it, or is "
					+ "not used."+"\n"
					+"=============================================================================================="+"\n"
                    + "#  Load Time-series File:\n"
					+ "This entry specifies the gene expression file that includes gene symbols and time-series values.\n"
					+ "The file has the following formatting restrictions:\n"
					+ "* The first line contains the column headers delimited by tabs.\n"
					+ "* The remaining lines contain the gene symbols and the data in sequential order of time.\n"
					+ "* In either the spot or gene field there can be multiple symbols listed delimited by "
					+ "either a pipe ('|'), comma (','), or a semi-colon (';').\n"
					+ "* If a value is missing between two time points then the field should be left empty giving "
					+ "two tabs between the non-missing values.\n"
					+ "SAMPLE FILE included:\n"
					+ "Gene	0h	1h	3h	6h	12h\n"
					+ "YAL053W	-0.027	0.158	0.169	0.193	-0.165\n"
					+ "YAL054C	0.183	-0.068	-0.134	-0.252	0.177\n"
					+ "YAL055W	-0.923	-0.51	-0.718	-0.512	-0.668"+"\n"
					+"=============================================================================================="+"\n"
					+"#  Normalization Method:\n"
					+ "All time-series will be transformed to starts at 0. "
					+ "This can be done in one of three ways based on the option selected to the left.  "
					+ "Given a time series vector of values for a gene (v_0,v_1,...,v_n) the options are:\n"
					+ "1.  'Log normalize data'  \u2212  the vector "
					+ "will be transformed to (0,log\u2082(v_1)\u2212log\u2082(v_0),...,log\u2082(v_n)\u2212log\u2082(v_0)).  "
					+ "Note that any values which are 0 or negative will be treated as missing.\n"
					+ "2.  'Normalize data'  \u2212  the vector "
					+ "will be transformed to (0,v_1\u2212v_0,...,v_n\u2212v_0)\n"
					+ "3.  'No normalization/add 0'  \u2212  a 0 will be inserted transforming the vector to "
					+ "(0,v_0,v_1,...,v_n)\n\n"
					+ "*If the data is not already in log space, then the "
					+ "'Log normalize data' should be selected.\n"
					+ "*If the data is already in log space and "
					+ "a time point 0 experiment was "
					+ "conducted, then the 'Normalize data' option should be selected.\n"
					+ "*If the data is already in log space and no time point 0 experiment was "
					+ "conducted, then the 'No normalization/add 0' option should be selected.";
			Util.renderDialog(this, szMessage, 50, 100);
		}else if(esource == help2) {
			szMessage = "#  Filter gene if it has no TF:\n"
					+ "* If this box is checked then genes are filtered if they are not included in the "
					+ "TF-gene interaction file.\n"
					+ "* If this box is unchecked then genes not included in the file "
					+ "are not filtered and are assumed "
					+ "to have a '0' for every entry.\n"
					+"=============================================================================================="+"\n"
					+"#  Maximum number of missing values:\n"
					+ "* A gene will be filtered if the number of time points for which there are no readings for the gene is "
					+ "greater than this parameter.\n"
					+ "* A gene will also be filtered if its expression value at the first time "
					+ "point is missing and 'log normalize data' or 'normalize data' was selected as the data transformation.\n"
					+"=============================================================================================="+"\n"
					+"#  Minimum correlation between repeats:\n"
					+ "This parameter only applies if there is repeat data over different time periods.\n"
					+ "* If there are two repeat data sets, the correlation value of the gene's expression level "
					+ "between the original and repeat must have "
					+ "a correlation value above this parameter, otherwise the gene will be filtered.\n"
					+ "* If there are three or more repeat sets, the mean pairwise correlation over all data sets "
					+ "must have a correlation value above this parameter, otherwise the gene will be filtered.\n"
					+"=============================================================================================="+"\n"
					+"#  Minimum absolute expression change:\n"
					+"After transformation (Log normalize data, Normalize data, or No Normalization/add 0), "
					+ "if the absolute value of a gene's expression value at every time point "
					+ "is below this threshold, then the gene will be filtered.\n"
					+"=============================================================================================="+"\n"
					+ "#  Change can be based on:\n"
					+ "This option defines how change is defined in the context of gene filter.\n"
					+ "* If 'Maximum\u2212minimum' option is selected, a gene will be filtered if the maximum "
					+ "absolute difference between "
					+ "the values of any two time points "
					+ "after transformation is less than the value of "
					+ "the 'Minimum absolute expression change' parameter.\n"
					+ "* If 'Compare to 0' is selected, a gene "
					+ "will be filtered if the absolute expression change from time point 0 at all time points is less than the "
					+ "value of the 'Minimum absolute expression change' parameter.\n"
					+ "* If 'Compare to previous' is selected, a gene "
					+ "will be filtered if the absolute fold-change between all adjacent time points is less than the "
					+ "value of the 'Minimum absolute expression change' parameter.\n"
					+ "* IF these options are multi-selected, genes lower than all the selected standards will be filtered.";
			Util.renderDialog(this, szMessage, 50, 100);
		}else if(esource == help3) {
			szMessage = "#  Maximum number of sub-paths:\n"
					+ "Determines the maximum number of sub-paths out of any split.\n"
					+"=============================================================================================="+"\n"
					+ "#  Minimum standard deviation:\n"
					+ "This parameter controls the minimum standard deviation on the Gaussian distributions. "
					+ "Increasing this parameter is recommended if applying TSMiner to RNA-seq data to avoid potential overfitting of "
					+ "low variance in expression due to small discrete counts.\n"
					+"=============================================================================================="+"\n"
					+"#  Use TF-TG interaction data to train the model\n"
					+ "* If this box is checked then the transcription TF-gene interaction data is used jointly with the time-series data to infer the model and then assign genes "
					+ "to paths of the model.\n"
					+ "* If this box is unchecked then the time-series data alone is used to infer a model, and "
					+ "the TF-gene interaction predictions are only used as a post-processing step which calculates enrichment scores "
					+ "based on the gene assignments.\n"
					+ "* Using the TF-gene interaction data to infer the model generally gives a more biologically coherent model.\n"
					+ "* Model learning is faster when do not use the TF-gene interaction data.\n"
					+"=============================================================================================="+"\n"
					+ "#  Saved model:\n"
					+ "This option is only relevant if a file is input in 'Saved Model File'.  "
					+ "* If it is set to 'As final' the model in the 'Saved Model File' is opened exactly as the resulting model\n"
					+ "* If it is set to 'Start from' TSMiner will start its search from the model saved.\n"
					+ "* If it is set to 'Do not use' then TSMiner will start model learning from a single chain.\n"
					+"=============================================================================================="+"\n"
					+ "#  Significance level of enrichment q-value and DE q-value:\n"
					+ "The 'Significance level of enrichment q-value' specifies the threshold below which a TF significantly regulates the genes involved in a sub-path.\n"
					+ "The 'Significance level of DE q-value' specifies the threshold below which a TF with significant enrichment q-value significantly up or down-regulates its targets.\n"
					+"=============================================================================================="+"\n"
					+ "#  Permutation test:\n"
					+ "This option specifies the compared time points when calculating the differential expression scores in permutation test.\n"
					+ "* If it is set to 'Compare to previous', the activation score for a TF is calculated based on the differential expression of the gene set it regulates between specific adjacent time points.\n"
					+ "* If it is set to 'Compare to 0', the activation score is based on the differential expression of gene set between specific time point and time 0.\n";
			Util.renderDialog(this, szMessage, 50, 100);
		}
	}

	/**
	 * Places szMessage in thedialog window with the title 'Help'
	 */
	public static void renderDialog(JDialog thedialog, String szMessage,
			int noffsetx, int noffsety) {
		renderDialog(thedialog, szMessage, noffsetx, noffsety, "Help");
	}

	/**
	 * Places szMessage in thedialog window with the title szTitle
	 */
	public static void renderDialog(JDialog thedialog, String szMessage,
			int noffsetx, int noffsety, String szTitle) {
		final JDialog thedialogf = thedialog;
		final JTextArea textAreaf = new JTextArea(szMessage);
		final int noffsetxf = noffsetx;
		final int noffsetyf = noffsety;
		final String szTitlef = szTitle;
		final int nlengthf = szMessage.length();
		javax.swing.SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				JDialog helpDialog = new JDialog(thedialogf, szTitlef, false);
				Container theHelpDialogPane = helpDialog.getContentPane();

				helpDialog.setBackground(Color.white);
				theHelpDialogPane.setBackground(Color.white);

				textAreaf.setLineWrap(true);
				textAreaf.setWrapStyleWord(true);

				textAreaf.setBackground(Color.white);
				textAreaf.setEditable(false);
				JScrollPane jsp = new JScrollPane(textAreaf);
				theHelpDialogPane.add(jsp);

				helpDialog.setLocation(thedialogf.getX() + noffsetxf,
						thedialogf.getY() + noffsetyf);
				if (nlengthf < 600) {
					helpDialog.setSize(700, 150);
				} else if (nlengthf < 1000) {
					helpDialog.setSize(700, 250);
				} else {
					helpDialog.setSize(700, 350);
				}

				helpDialog.setDefaultCloseOperation(DISPOSE_ON_CLOSE);
				helpDialog.setVisible(true);
			}
		});
	}

	/**
	 * Toggles which model selection variables are enabled based on the value of
	 * btraintest
	 */
	private void toggleEnabled(boolean btraintest) {
		mergeLabeldiff.setEnabled(btraintest);
		mergeLabel.setEnabled(btraintest);
		delaypathLabeldiff.setEnabled(btraintest);
		delaypathLabel.setEnabled(btraintest);
		prunepathLabel.setEnabled(btraintest);
		prunepathLabeldiff.setEnabled(btraintest);
		epsilonLabel.setEnabled(btraintest);
		epsilonLabeldiff.setEnabled(btraintest);
		thespinnerepsilondiff.setEnabled(btraintest);
		seedLabel.setEnabled(btraintest);
		nodepenaltyLabel.setEnabled(!btraintest);
	}


	/**
	 * Create the GUI and show it. For thread safety, this method should be
	 * invoked from the event-dispatching thread.
	 */
	private static void createAndShowGUI() throws FileNotFoundException,IOException {
		// Make sure we have nice window decorations.
		// JFrame.setDefaultLookAndFeelDecorated(true);
		// Create and set up the window.
		JFrame frame = new TSMiner_IO();
		frame.setLocation(10, 25);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

		// Display the window.
		frame.pack();
		frame.setVisible(true);
	}


	/**
	 * The main method which when executed will have the input interface created
	 */
	public static void main(String[] args) throws Exception {
		
		boolean bshowusage = true;

		if (args.length == 0) {
			bshowusage = false;
		} else if (args.length == 2) {
			if (args[0].equals("-d")) {
				szDefaultFile = args[1];
				bshowusage = false;
			}
		}

		if (bshowusage) {
			System.out.println("USAGE: java -jar nrem.jar [-d defaultfilename.txt|-b settingsfile.txt outmodelfile.txt]");
			return;
		}
		
			javax.swing.SwingUtilities.invokeLater(new Runnable() {
				public void run() {
					try {
						createAndShowGUI();
					} catch (FileNotFoundException ex) {
						ex.printStackTrace(System.out);
					} catch (IOException ex) {
						ex.printStackTrace(System.out);
					}
				}
			});
	}
}
