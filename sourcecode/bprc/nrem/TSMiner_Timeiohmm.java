package bprc.nrem;

import bprc.core.*;

import javax.swing.*;

import java.util.*;
import java.util.List;
import java.util.Map.Entry;
import java.text.*;
import java.io.*;
import java.awt.*;
import java.math.*;
import java.awt.event.*;

import edu.umd.cs.piccolo.nodes.PText;

/**
 * This class implements the core methods for learning the TSMiner maps
 */
public class TSMiner_Timeiohmm {
	/****************New added Parameters***************/
	List<List<Integer>> upgenes = new ArrayList<List<Integer>>();
	List<List<Integer>> downgenes = new ArrayList<List<Integer>>();
	int Maxcount = 1000;
	int minnc = 10;
	int minbic = 10;
	int countpath;
	/****************New added Parameters***************/
	int nglobaltime;
	double nodepenalty = 10;
	static final boolean BDEBUG = false;
	static final boolean BDEBUGMODEL = false;
	boolean ballowmergeval = false;
	boolean bsavedchange = false;
	boolean BEQUALSTD = false;
	boolean BVITERBI = false;
	boolean BREGTSMiner = false;
	boolean bhasmerge = false;
	double MINPROB = .00000000000000000001;
	double DEFAULTSIGMA = .5;
	double MAXFUTUREITR = 30;
	double MERGEMIN = -.0015;
	double DELAYMIN = -.0015;

	int nrandomseed;
	int nglobaliteration;
	int ninitsearchval;
	int npermutationval;
	int Countsubpath;
	int Countnode;

	double RESPLITMIN = -.0015;

	boolean brealXaxisDEF;
	double dYaxisDEF;
	double dXaxisDEF;
	double dnodekDEF;
	int nKeyInputTypeDEF;
	double dKeyInputXDEF;
	double dpercentDEF;
	JButton currentButton;

	List <Integer> pathh = new ArrayList<Integer>(); 

	static final String SCORESDIR = "TOPAALoc";
	static final String SZDELIM = "|;,";
	TSMiner_NaiveBayes filteredClassifier;  //调用贝叶斯分类器，判断基因是否被删除

	double dprevouterbestlog;
	double dtrainlike;
	boolean bfilterbinding = false;
	double EPSILON = .01;
	double BEPSILON = 0.01;
	int nmaxchild;
	int ntotalcombined;// num filtered and non-filtered

	double[][] CONSTANTA;

	int MINPATH = 5;
	static final double RIDGE = 1;
	static final double MINSIGMA = .0001;
	double dminstddev = 0.5;
	double enrichmentthreshold = 0.01;
	double dethreshold = 0.01;

	String[] testgenenames;
	/**
	 * The subset of theDataSet.data that belongs to the training set??? The
	 * expression data: row are genes, columns are time points in experiments
	 */
	double[][] traindata;
	String[] traingenenames;
	double[][] testdata;

	List<Double> Sigma = new ArrayList<Double>();

	/**
	 * Object containing all of the Regulator-Gene binding data.
	 */
	RegulatorBindingData bindingData;
	GOBindingData raw_gobinding;
	GOBindingData filtered_gobinding;
	RegulatorBindingData trainBinding;

	DataSetCore miRNADataSet;
	boolean scaleTFExp; // should TF expression be scaled?
	boolean scaleMIRNAExp; // should miRNA be used as a regulator
	HashMap<String, Integer> reg2DataSetIndex;
	HashMap<String, Integer> geneToDataSetIndex;
	HashMap<String, Integer> raw_reg2DataSetIndex;
	HashMap<String, Integer> raw_geneToDataSetIndex;

	double expressionScalingFactor = 1.0;
	double minTFExpAfterScaling = 0.5;
	Treenode treeptr;
	NumberFormat nf2;
	NumberFormat nf3;
	NumberFormat nf4;
	double dbestlog;
	Treenode bestTree = null;
	public TSMiner_DataSet theDataSet;
	public TSMiner_DataSet rawDataSet;

	JTextArea runtext;
	int numtotalPath = 1;
	double dbesttrainlike;

	/**
	 * Number of boolean items for the static vector. Equal to the number of TFs
	 */
	int numbits;
	double TRAINRATIO = .75;

	Random theRandom;

	int ntrain;
	int ntest;
	int numcols;
	int nholdout;

	double HOLDOUTRATIO = 0;
	Double[][] holdoutdata;
	double[][][] holdoutpval;
	int[][] holdoutpma;
	int[][][] holdoutpvalIndex;

	boolean[] bholdout;
	BigInteger[] storedbestpath;
	String szinitfileval;
	String szexpname;
	double convergenceval = 0.01;
	String szbinding;
	String szgobinding;
	ArrayList<Color> savedColors = new ArrayList<Color>();
	double dgloballike;
	private int ntotalID;

	Hashtable<Treenode, BackPtrRec> htBackPts = new Hashtable<Treenode, BackPtrRec>();

	String szstaticsourceval;

	HashMap<String,String> gene2FastaMap;
	HashMap<String, List<Double>> time_TF;

	/**
	 * The probability that binding is functional, i.e. that the genes a TF
	 * binds are regulated by the TF. Used to calculated TF activity scores. Can
	 * be set in the TSMiner defaults file.
	 */
	double dProbBindingFunctional = 0.8;
	/**
	 * All TF activiy priors are adjusted by ACTIVITY_EPSILON to avoid divide by
	 * 0 errors when the TF activity prior is 1.
	 */
	static final double ACTIVITY_EPSILON = 1E-8;
	/**
	 * An array of activity scores for all TFs based on how the genes the TF
	 * regulates transition to the next states compared to how all genes into
	 * the split transistion. The max is the max score over this node and its
	 * descendent nodes.
	 */
	double[] dMaxTFActivityScore;
	String regPriorsFile = "";
	
	public HashMap<String, List<String>> activeTF;
	public HashMap<String, List<String>> activeDEpvalue;
	public HashMap<String, List<String>> repressTF;
	public HashMap<String, List<String>> repressDEpvalue;
	HashMap<String, List<Integer>> activeindex;
	HashMap<String, List<Integer>> repressindex;
	HashMap<String, List<Integer>> activegenes;
	HashMap<String, List<Integer>> repressgenes;
	HashMap<Integer, Double[]> minposipvalue;
	HashMap<Integer, Double[]> minnegapvalue;
	HashMap<Integer, Integer> filter2rawIndex;
	HashMap<Integer, Integer> raw2filterIndex;
	
	
	/**
	 * Class constructor - provides the execution control on the method
	 */
	public TSMiner_Timeiohmm(TSMiner_DataSet rawDS, TSMiner_DataSet theDS, String szbinding, String szgobinding,
			String sznumchildval, String szseedval, boolean bstaticcheckval, boolean ballowmergeval,
			int ninitsearchval, int npermutationval, String szinitfileval, JTextArea runtext,
			JButton endSearchButton, boolean bstaticsearchval, final boolean brealXaxisDEF,
			final double dYaxisDEF, final double dXaxisDEF,
			final double dnodekDEF, final int nKeyInputTypeDEF,
			final double dKeyInputXDEF, final double dpercentDEF,
			final JButton currentButton, String szexpname,
			String szminstddeval, String enrichmentthreshold, String dethreshold, 
			String szstaticsourceval, boolean scaleTFExp,
			boolean scaleMIRNAExp, DataSetCore miRNAExpressionDataCore, 
			String regPriors, double dProbBindingFunctional,
			double expressionScalingFactor, double minTFExpAfterScaling)
			throws Exception {
		this.szstaticsourceval = szstaticsourceval;
		this.szbinding = szbinding;
		this.szgobinding = szgobinding;
		this.dminstddev = Double.parseDouble(szminstddeval);
		this.enrichmentthreshold = Double.parseDouble(enrichmentthreshold);
		this.dethreshold = Double.parseDouble(dethreshold);
		this.runtext = runtext;
		this.szexpname = szexpname;
		this.ballowmergeval = ballowmergeval;
		this.brealXaxisDEF = brealXaxisDEF;
		this.dYaxisDEF = dYaxisDEF;
		this.dXaxisDEF = dXaxisDEF;
		this.dnodekDEF = dnodekDEF;
		this.nKeyInputTypeDEF = nKeyInputTypeDEF;
		this.dKeyInputXDEF = dKeyInputXDEF;
		this.dpercentDEF = dpercentDEF;
		this.currentButton = currentButton;
		this.ninitsearchval = ninitsearchval;
		this.npermutationval = npermutationval;
		this.scaleTFExp = scaleTFExp;
		this.scaleMIRNAExp = scaleMIRNAExp;
		this.dProbBindingFunctional = dProbBindingFunctional;
		this.expressionScalingFactor = expressionScalingFactor;
		this.minTFExpAfterScaling = minTFExpAfterScaling;
		this.regPriorsFile = regPriors;

		BREGTSMiner = (!bstaticsearchval);
		this.szinitfileval = szinitfileval;
		this.bfilterbinding = bstaticcheckval;

		nf2 = NumberFormat.getInstance(Locale.ENGLISH);
		nf2.setMinimumFractionDigits(2);
		nf2.setMaximumFractionDigits(2);

		nf3 = NumberFormat.getInstance(Locale.ENGLISH);
		nf3.setMinimumFractionDigits(3);
		nf3.setMaximumFractionDigits(3);
		
		nf4 = NumberFormat.getInstance(Locale.ENGLISH);
		nf4.setMinimumFractionDigits(8);
		nf4.setMaximumFractionDigits(8);
		

		BEPSILON = convergenceval / 100; //>0
		nrandomseed = Integer.parseInt(szseedval);
		theRandom = new Random(nrandomseed);
		nmaxchild = Integer.parseInt(sznumchildval);
		theDataSet = theDS;
		rawDataSet = rawDS; 
		filter2rawIndex = new HashMap<Integer, Integer>();
		for(int i=0;i<theDataSet.genenames.length;i++){
			String gene = theDataSet.genenames[i];
			for(int j=0;j<rawDataSet.genenames.length;j++){
				if(gene.equals(rawDataSet.genenames[j])){
					filter2rawIndex.put(i, j);
					break;
				}
			}
		}
		raw2filterIndex = new HashMap<Integer, Integer>();
		for(int i=0;i<rawDataSet.genenames.length;i++){
			String gene = rawDataSet.genenames[i];
			for(int j=0;j<theDataSet.genenames.length;j++){
				if(gene.equals(theDataSet.genenames[j])){
					raw2filterIndex.put(i, j);
					break;
				}
			}
		}
		
		HashMap<String,Integer> timepointName2Int = new HashMap<String,Integer>(); //键值为样本名称，值为字符串。将各列列头对应到整数，从0开始
		minposipvalue = new HashMap<Integer, Double[]>();
		minnegapvalue = new HashMap<Integer, Double[]>();
		for(int i = 0; i < theDataSet.dsamplemins.length; i++) timepointName2Int.put(theDataSet.dsamplemins[i].toUpperCase(), i);
		miRNADataSet = miRNAExpressionDataCore;
		CONSTANTA = new double[nmaxchild + 1][];
		for (int nindex = 1; nindex <= nmaxchild; nindex++) {
			CONSTANTA[nindex] = new double[nindex];
			for (int njindex = 0; njindex < nindex; njindex++) CONSTANTA[nindex][njindex] = 1.0 / nindex;
		}
		
		if (szbinding.equals("")) {
			BREGTSMiner = true;
		} else if(!bfilterbinding){

			bindingData = new RegulatorBindingData(szbinding, theDataSet, 
					SZDELIM, BDEBUG, regPriors, null, timepointName2Int);
			
			regname2index();
			
			bindingData.adjustBindingData(reg2DataSetIndex, theDataSet, miRNADataSet, 
					scaleMIRNAExp, scaleTFExp, expressionScalingFactor, minTFExpAfterScaling);
			
			numbits = bindingData.numberRegs; //TFs number
			dMaxTFActivityScore = new double[numbits]; //A matrix saving the Maximum score of each TF
			for (int ntf = 0; ntf < bindingData.numberRegs; ntf++) {
				dMaxTFActivityScore[ntf] = 0;
			}
		}else{
			bindingData = new RegulatorBindingData(szbinding, theDataSet, 
					SZDELIM, BDEBUG, regPriors, null, timepointName2Int);
			regname2index();
			bindingData.adjustBindingData(reg2DataSetIndex, theDataSet, miRNADataSet, 
					scaleMIRNAExp, scaleTFExp, expressionScalingFactor, minTFExpAfterScaling);
			theDataSet = filternotf(theDataSet, bindingData);
			
			bindingData = new RegulatorBindingData(szbinding, theDataSet, 
					SZDELIM, BDEBUG, regPriors, null, timepointName2Int);
			regname2index();
			bindingData.adjustBindingData(reg2DataSetIndex, theDataSet, miRNADataSet, 
					scaleMIRNAExp, scaleTFExp, expressionScalingFactor, minTFExpAfterScaling);
			
			numbits = bindingData.numberRegs; //TFs number
			dMaxTFActivityScore = new double[numbits]; //A matrix saving the Maximum score of each TF
			for (int ntf = 0; ntf < bindingData.numberRegs; ntf++) {
				dMaxTFActivityScore[ntf] = 0;
			}
		}
		
		
		if(szgobinding != null && szgobinding.length() > 0){
			raw_gobinding = new GOBindingData(szgobinding, "",
					false, rawDataSet, SZDELIM, BDEBUG, regPriors,
					null, timepointName2Int);
			filtered_gobinding = new GOBindingData(szgobinding, "",
					false, theDataSet, SZDELIM, BDEBUG, regPriors,
					null, timepointName2Int);
		}

		runtext.append("Gene number: "+theDataSet.data.length+"  "+"TF number: "+numbits+"\n");
		runtext.paintImmediately(runtext.getBounds());
		
		splitdata();

		if(minbic < traindata.length/200) minbic = traindata.length/200;
		treeptr = new Treenode();
		int numintervals = traindata[0].length - 1;
		int[] path = new int[numintervals];

		boolean binitfile = (!szinitfileval.equals(""));
		if ((binitfile) && (ninitsearchval == 0)) {
			readInitTree(treeptr);
			computeOrders(treeptr);
			combineTrainAndTest();
			viterbi(treeptr, theDataSet.data);
			runtext.append("Use Saved Model: Done"+"\n");
			runtext.paintImmediately(runtext.getBounds());
			
			if (bindingData.gene2RegBinding != null) {
				runtext.append("================================"+"\n");
				runtext.append("TF enrichment analysis..."+"\n");
				runtext.paintImmediately(runtext.getBounds());	
				computeStats(treeptr, treeptr);
				runtext.append("\n"+"TF enrichment analysis Done."+"\n");
				runtext.paintImmediately(runtext.getBounds());
				computeminparentlevel(treeptr);
			}
			bsavedchange = true;
		} else {
			bsavedchange = false;
			if ((binitfile) && (ninitsearchval == 1)) {
				readInitTree(treeptr);
			} else {
				double[] dsigma = new double[traindata[0].length];
				double[] dmeans = new double[traindata[0].length];
				computeDataStats(traindata, dsigma, dmeans);
				buildEmptyTree(0, treeptr, dsigma, dmeans);
			}

			nodepenalty = Math.log(traindata.length*traindata[0].length);
			runtext.append("================================"+"\n");
			runtext.append("Penalty for each node: "+nodepenalty+"\n");
			runtext.paintImmediately(runtext.getBounds());
			
			currentButton.setEnabled(true);
			endSearchButton.setEnabled(true);
			searchstage1();

			int nremoved = 0;

			if (endSearchButton != null) {
				endSearchButton.setEnabled(false);
			}

			boolean bagaindelay;
			boolean bagaindelayouter;
			int numdelay = 0;
			
			runtext.append("================================"+"\n");
			runtext.append("Delay split..."+"\n");
			runtext.paintImmediately(runtext.getBounds());
			do {
				bagaindelayouter = false;
				for (int ndesiredlevel = 1; ndesiredlevel < numintervals; ndesiredlevel++) {
					do {

						if (BDEBUG) {
							System.out.println("trying to delay "
									+ ndesiredlevel + " " + dbestlog);
						}

						dprevouterbestlog = dbestlog;

						dbestlog = Double.NEGATIVE_INFINITY;
						bagaindelay = traverseanddelay(path, treeptr, ndesiredlevel, treeptr);

						if (dbestlog == Double.NEGATIVE_INFINITY) {
							dbestlog = dprevouterbestlog;
						} else {
							numdelay++;
							treeptr = bestTree;
							String szimprovedelayDIFF = nf3.format(dbestlog
									- dprevouterbestlog);

							runtext.append("Delayed path: "+ numdelay+"  ||  Likelihood improvement: "+szimprovedelayDIFF+"\n");
							runtext.paintImmediately(runtext.getBounds());
						}

						treeptr = bestTree;

						if (bagaindelay) {
							bagaindelayouter = true;
						}
					} while (bagaindelay);
				}
			} while (bagaindelayouter);

			if (ballowmergeval) {
				traverseandmerge(path);
			}
			treeptr = bestTree;

			traverse(bestTree, 0, false);
			if (BDEBUG) {
				System.out.println("Test set " + dbestlog);
				System.out.println("calling combine tranandtest");
			}

			combineTrainAndTest();

			double dtemplog = trainhmm(treeptr, true);
			nglobaltime++;
			int numnodes = countNodes(treeptr, nglobaltime);
			if (BDEBUGMODEL) {
				System.out.println(nodepenalty + "\t final train is\t"
						+ dtemplog + "\t" + numnodes + "\t"
						+ (dtemplog + numnodes * nodepenalty) + "\t"
						+ dgloballike);
			}
			traverse(bestTree, 0, false);

			boolean bagain; 

			runtext.append("================================"+"\n");
			runtext.append("Removed the Path < 5 genes..."+"\n");
			runtext.paintImmediately(runtext.getBounds());
			
			do { //Delete path less than 5 genes
				bagain = false;
				computeOrders(treeptr);
				viterbi(treeptr, traindata);			
				
				int[] bestpath = new int[path.length];
				for (int nindex = 0; nindex < bestpath.length; nindex++) {
					bestpath[nindex] = -1;
				}

				MinPathRec theMinPathRec = traverseanddeleteMinPath(path, bestpath, treeptr);
				
				if (BDEBUG) {
					System.out.println("min is " + theMinPathRec.nval + "\t"
							+ "level is " + theMinPathRec.nlevel);
					for (int nindex = 0; nindex < bestpath.length; nindex++) {
						System.out.println(theMinPathRec.bestpath[nindex]);
					}
				}
				
				if (theMinPathRec.nval < MINPATH && theMinPathRec.nval!=0) {//删除5个基因以下的路径

					deleteMinPath(theMinPathRec.bestpath, theMinPathRec.nlevel, treeptr); //删除包含最少基因的路径
					bagain = true;
					traverse(treeptr, 0, true);

					trainhmm(treeptr, true);

					if (BDEBUG) {
						System.out.println("after retrain");
					}
					traverse(treeptr, 0, true);

					nremoved++;
					numtotalPath--;
					
					runtext.append("Removed Path: "+nremoved+"\n");
					runtext.paintImmediately(runtext.getBounds());
				}
			} while (bagain);
			
			if (bindingData.gene2RegBinding != null) {
				runtext.append("================================"+"\n");
				runtext.append("TF enrichment analysis..."+"\n");
				runtext.paintImmediately(runtext.getBounds());	
				
				computeStats(treeptr, treeptr);
				
				runtext.append("\n"+"TF enrichment analysis Done."+"\n");
				runtext.paintImmediately(runtext.getBounds());
				
				computeminparentlevel(treeptr);
			}
		}
		
		/*
		String[][] TreeModle_DREM = viterbi2(treeptr);
		String[][] result = new String[TreeModle_DREM.length][TreeModle_DREM[0].length+1];
		for(int i=0;i<TreeModle_DREM.length;i++){
			result[i][0] = theDataSet.genenames[i];
			result[i][1] = TreeModle_DREM[i][(TreeModle_DREM[0].length-1)/2];
			for(int j=0;j<(TreeModle_DREM[0].length-1)/2-1;j++){
				result[i][j+2] = TreeModle_DREM[i][j];
			}
			for(int j=(TreeModle_DREM[0].length-1)/2;j<TreeModle_DREM[0].length;j++){
				result[i][j+1] = TreeModle_DREM[i][j];
			}
		}
		try{
			BufferedWriter outXml = new BufferedWriter(new FileWriter("E:/ClusteringResult.txt"));
			for(int p=0;p<result.length;p++){
				for(int q=0;q<result[0].length;q++){
					if(result[p][q] != null && result[p][q].length() > 0){	
					outXml.write(result[p][q]);
					outXml.write("\t");
					}else{
						outXml.write("");
						outXml.write("\t");
					}
				}
					outXml.newLine();
			}
			outXml.flush(); 
			outXml.close();
			System.out.println("DONE"); 
			
			
		}catch (Exception e) { 
			System.out.println("FALSE"); 
		e.printStackTrace(); 
		 }
		*/
		
		path = new int[numintervals+1];
		if(numtotalPath >= 2 && numbits > 0){
			
			activeTF = new HashMap<String, List<String>>();
			activeDEpvalue = new HashMap<String, List<String>>();
			repressTF = new HashMap<String, List<String>>();
			repressDEpvalue = new HashMap<String, List<String>>();
			activeindex = new HashMap<String, List<Integer>>();
			repressindex = new HashMap<String, List<Integer>>();
			activegenes = new HashMap<String, List<Integer>>();
			repressgenes = new HashMap<String, List<Integer>>();
			
			time_TF = new HashMap<String, List<Double>>();
			
			runtext.append("================================"+"\n");
			runtext.append("Permutation test..."+"\n");
			runtext.paintImmediately(runtext.getBounds());
			ESpermute(treeptr, 0, path);
			runtext.append("\n"+"Permutation test Done."+"\n");
			runtext.paintImmediately(runtext.getBounds());
			outputTF(treeptr, 0, path);
		}
		
	}
	
	public TSMiner_DataSet filternotf(TSMiner_DataSet theDataSet, 
			RegulatorBindingData bindingData) {
		int numrows = theDataSet.numrows;
		boolean[] goodrow = new boolean[numrows];
		int ngoodrows = 0;

		for (int nrow = 0; nrow < numrows; nrow++) {
			int[] tf = bindingData.gene2RegBindingIndex[0][nrow];
			if (tf != null && tf.length > 0) {
				goodrow[nrow] = true;
				ngoodrows++;
			} else {
				goodrow[nrow] = false;
			}
		}
		TSMiner_DataSet theDataSetsMerged = theDataSet;
		theDataSetsMerged = new TSMiner_DataSet(theDataSetsMerged.filtergenesgeneral(goodrow, ngoodrows, false),
				theDataSetsMerged.tga);
		return theDataSetsMerged;
	}
	
	private void regname2index() throws Exception{
		//Read the TFs data and Create a hash map that takes regs and maps them to their row in
		// the exp data
		reg2DataSetIndex = new HashMap<String, Integer>();   //TF name to its index in theDataSet  
		geneToDataSetIndex = new HashMap<String, Integer>(); //Gene name to its index in theDataSet
		for (int i = 0; i < theDataSet.genenames.length; i++) geneToDataSetIndex.put(theDataSet.genenames[i].toUpperCase(),new Integer(i));
		for (int i = 0; i < bindingData.regNames.length; i++) {
			String upCaseReg = bindingData.regNames[i].toUpperCase();
			if (geneToDataSetIndex.containsKey(upCaseReg.split("_")[0])) {
				reg2DataSetIndex.put(upCaseReg, geneToDataSetIndex.get(upCaseReg.split("_")[0]));
			} else {
				HashSet<String> syns = bindingData.regSyns.get(upCaseReg);
				boolean found = false;
				if (syns != null) {
					for (String syn : syns) {
						if (geneToDataSetIndex.containsKey(syn.split("_")[0])) {
							if (found) {
								throw new Exception(
										"Multiple syns for "
												+ bindingData.regNames[i]
												+ ",a regulator in the"
												+ " binding data, are in the expression data.");
							}
							reg2DataSetIndex.put(upCaseReg,geneToDataSetIndex.get(syn.split("_")[0]));
							found = true;
						}
					}
				}
			}
		}
		
		raw_reg2DataSetIndex = new HashMap<String, Integer>();   //TF name to its index in the raw DataSet  
		raw_geneToDataSetIndex = new HashMap<String, Integer>(); //Gene name to its index in the raw DataSet
		for (int i = 0; i < rawDataSet.genenames.length; i++) raw_geneToDataSetIndex.put(rawDataSet.genenames[i].toUpperCase(),new Integer(i));
		for (int i = 0; i < bindingData.regNames.length; i++) {
			String upCaseReg = bindingData.regNames[i].toUpperCase();
			if (raw_geneToDataSetIndex.containsKey(upCaseReg.split("_")[0])) {
				raw_reg2DataSetIndex.put(upCaseReg, raw_geneToDataSetIndex.get(upCaseReg.split("_")[0]));
			} else {
				HashSet<String> syns = bindingData.regSyns.get(upCaseReg);
				boolean found = false;
				if (syns != null) {
					for (String syn : syns) {
						if (raw_geneToDataSetIndex.containsKey(syn.split("_")[0])) {
							if (found) {
								throw new Exception(
										"Multiple syns for "
												+ bindingData.regNames[i]
												+ ",a regulator in the"
												+ " binding data, are in the expression data.");
							}
							raw_reg2DataSetIndex.put(upCaseReg,raw_geneToDataSetIndex.get(syn.split("_")[0]));
							found = true;
						}
					}
				}
			}
		}
	}
	
	public void ESpermute(Treenode treehmm, int ndepth, int[] path){
	    if(numbits>0 && ndepth > 0){
	    	if(treehmm.parent.numchildren>1 && treehmm.parent.dpvalEdgeSplit[path[ndepth-1]] != null){
	    			List<Integer> rtime = new ArrayList<Integer>();
	    			getRtime(treehmm, ndepth, rtime);
	    				treehmm.splitES = new Double[numbits][rtime.size()*2];
		    			int tfcount = 0;
		    			if(npermutationval == 1){
		    				for(int i=0;i<numbits;i++){
			    				if(treehmm.parent.dpvalEdgeSplit[path[ndepth-1]][i] <= enrichmentthreshold){
			    					List<Integer> gs = new ArrayList<Integer>();
		    						for(int j=0;j<bindingData.reg2GeneBindingIndex[ndepth][i].length;j++){
			    						if(bindingData.reg2GeneBindingSigned[ndepth][i][j] != 0){
			    						//&& treehmm.genelist.contains(bindingData.reg2GeneBindingIndex[ndepth][i][j])){
			    							gs.add(bindingData.reg2GeneBindingIndex[ndepth][i][j]);
			    						}
			    					}
			    					double[] express = new double[gs.size()*2];
			    					if(express.length>2){
			    						tfcount++;
			    						for(int j=0;j<rtime.size();j++){
			    							String timetf = rtime.get(j)+"_"+i;
			    							if(time_TF.get(timetf) != null){
			    								List<Double> es = time_TF.get(timetf);
			    								treehmm.splitES[i][j] = es.get(0);
				        						treehmm.splitES[i][j+rtime.size()] = es.get(1);	
			    							}else{
			    								for(int m=0;m<gs.size();m++){
				    	    						express[m] = theDataSet.data[gs.get(m)][rtime.get(j)-1];
				    	    	    				express[m+gs.size()] = theDataSet.data[gs.get(m)][rtime.get(j)];
				    	    	    			}
			    								
				    							DEGeneset degs = new DEGeneset();
				        						Double[] es = degs.permutetest(express);
				        						treehmm.splitES[i][j] = es[0];
				        						treehmm.splitES[i][j+rtime.size()] = es[1];
				        						List<Double> ES = new ArrayList<Double>();
				        						ES.add(es[0]);
				        						ES.add(es[1]);
				        						time_TF.put(timetf, ES);
			    							}
			    						}
			    					}
			    					
			    				}
			    			}
		    			}else{
		    				for(int i=0;i<numbits;i++){
			    				if(treehmm.parent.dpvalEdgeSplit[path[ndepth-1]][i] <= enrichmentthreshold){
			    					//int[] gs = bindingData.reg2GeneBindingIndex[ndepth][i];
			    					List<Integer> gs = new ArrayList<Integer>();
		    						for(int j=0;j<bindingData.reg2GeneBindingIndex[ndepth][i].length;j++){
			    						if(bindingData.reg2GeneBindingSigned[ndepth][i][j] != 0){
			    							//&& treehmm.genelist.contains(bindingData.reg2GeneBindingIndex[ndepth][i][j])){
			    							gs.add(bindingData.reg2GeneBindingIndex[ndepth][i][j]);
			    						}
			    					}
			    					double[] express = new double[gs.size()*2];
			    					if(express.length>2){
			    						tfcount++;
			    						for(int j=0;j<rtime.size();j++){
			    							String timetf = rtime.get(j)+"_"+i;
			    							if(time_TF.get(timetf) != null){
			    								List<Double> es = time_TF.get(timetf);
			    								treehmm.splitES[i][j] = es.get(0);
				        						treehmm.splitES[i][j+rtime.size()] = es.get(1);
				        						
			    							}else{
			    								for(int m=0;m<gs.size();m++){
				    	    						express[m] = 0;
				    	    	    				express[m+gs.size()] = theDataSet.data[gs.get(m)][rtime.get(j)];
				    	    	    			}
			    								
				    							DEGeneset degs = new DEGeneset();
				        						Double[] es = degs.permutetest(express);
				        						treehmm.splitES[i][j] = es[0];
				        						treehmm.splitES[i][j+rtime.size()] = es[1];
				        						List<Double> ES = new ArrayList<Double>();
				        						ES.add(es[0]);
				        						ES.add(es[1]);
				        						time_TF.put(timetf, ES);
			    							}
			    						}
			    					}
			    					
			    				}
			    			}
		    			}
		    			for(int j=0;j<rtime.size();j++){
		    				Double[] pvalue = new Double[numbits];
		    				for(int i=0;i<numbits;i++){
		    					if(treehmm.splitES[i][j] != null){
		    						pvalue[i] = treehmm.splitES[i][j];
		    					}
		    				}
		    				for(int i=0;i<numbits;i++){
		    					if(pvalue[i] != null){
		    						int index = getrank(pvalue, pvalue[i]);
			    						treehmm.splitES[i][j] = treehmm.splitES[i][j]*tfcount/index;
			    						if(treehmm.splitES[i][j]>1){
			    							treehmm.splitES[i][j] = 1.0;
			    						}
		    					}
		    					
		    				}
		    			}

	    			Countnode++;
	    			runtext.append(Countnode+" sub-paths done..."+"  ");
	    			runtext.paintImmediately(runtext.getBounds());
	    	}
		 }
	   for (int nchild = 0; nchild < treehmm.numchildren; nchild++) {
	        path[ndepth] = nchild;
	        ESpermute(treehmm.nextptr[nchild], ndepth+1,path);
	   }
	  }
	
	public Integer getrank(Double[] pvalue, double value){
    	int cc = 1;
    	for(int i=0;i<pvalue.length;i++){
    		if(pvalue[i] != null){
        		if(pvalue[i]<value){
        			cc++;
        		}
    		}
    	}
    	return cc;
    }
	
	public void getRtime(Treenode treehmm, int ndepth, List<Integer> rtime) {
			rtime.add(ndepth);
			if(treehmm.numchildren>=1){
				getRtime(treehmm.nextptr[0], ndepth+1, rtime);
			}
	}
	
	public void outputTF(Treenode treehmm, int ndepth, int[] path){
        if(numbits>0 && ndepth>0){
        	String name = "";
            for(int i=0;i<ndepth;i++){
                name=name+path[i];
            }
        	if(treehmm.parent.numchildren>1 && treehmm.parent.dpvalEdgeSplit != null){
        			int enrich = 0;
					for(int i=0;i<numbits;i++){
						double pvalue = treehmm.parent.dpvalEdgeSplit[path[ndepth-1]][i];
						if(pvalue <= enrichmentthreshold) enrich++;
					}
						if(enrich > 0){
	                		treehmm.resultTF = new String[numbits][5+treehmm.splitES[0].length];
	    					for(int i=0;i<numbits;i++){
	    						treehmm.resultTF[i][0] = bindingData.regNames[i];
	    						treehmm.resultTF[i][1] = bindingData.reg2GeneBindingIndex[ndepth][i].length+"";
	    						for (int j = 0; j < bindingData.signedBindingValuesSorted.length; j++) {
	    							  if(bindingData.signedBindingValuesSorted[j] != 0){
	    								treehmm.resultTF[i][2] = treehmm.parent.ncountTotals[i][j]+"";
	    		    					treehmm.resultTF[i][3] = treehmm.parent.ncountvals[i][path[ndepth-1]][j]+"";
	    							  }
	    						 }
	    						treehmm.resultTF[i][4] = treehmm.parent.dpvalEdgeSplit[path[ndepth-1]][i]+"";

	    						for(int j=0;j<treehmm.splitES[i].length;j++){
	    							if(treehmm.splitES[i][j] != null){
	    								treehmm.resultTF[i][5+j] = nf4.format(treehmm.splitES[i][j]);
	    							}else{
	    								treehmm.resultTF[i][5+j] = "NA";
	    							}
	    						}
	    					}
	                	}else{
	                		treehmm.resultTF = new String[numbits][5];
	    					for(int i=0;i<numbits;i++){
	    						treehmm.resultTF[i][0] = bindingData.regNames[i];
	    						treehmm.resultTF[i][1] = bindingData.reg2GeneBindingIndex[ndepth][i].length+"";
	    						for (int j = 0; j < bindingData.signedBindingValuesSorted.length; j++) {
	    							  if(bindingData.signedBindingValuesSorted[j] != 0){
	    								treehmm.resultTF[i][2] = treehmm.parent.ncountTotals[i][j]+"";
	    		    					treehmm.resultTF[i][3] = treehmm.parent.ncountvals[i][path[ndepth-1]][j]+"";
	    							  }
	    						 }
	    						treehmm.resultTF[i][4] = treehmm.parent.dpvalEdgeSplit[path[ndepth-1]][i]+"";
	    					}
	                	}
  					
    				giveDEscore2(treehmm, treehmm, 0);
    				BubbleSort(treehmm.resultTF, treehmm.resultTF.length, 4);
        	}
        	giveDEscore(treehmm, treehmm, 0);  	
    	 }
       for (int nchild = 0; nchild < treehmm.numchildren; nchild++) {
            path[ndepth] = nchild;
            outputTF(treehmm.nextptr[nchild], ndepth+1,path);
       }
      }
	
     public void giveDEscore(Treenode treeroot, Treenode treehmm, int count) {
    		 if(treehmm.parent.numchildren>1 && treehmm.resultTF!=null){
    					  if(treehmm.resultTF[0].length>5){
        					  treeroot.deTF = new String[treehmm.resultTF.length][7];
        					  for(int i=0;i<treehmm.resultTF.length;i++){
        						  treeroot.deTF[i][0] = treehmm.resultTF[i][0];
        						  treeroot.deTF[i][1] = treehmm.resultTF[i][1];
        						  treeroot.deTF[i][2] = treehmm.resultTF[i][2];
        						  treeroot.deTF[i][3] = treehmm.resultTF[i][3];
        						  treeroot.deTF[i][4] = treehmm.resultTF[i][4];
        						  treeroot.deTF[i][5] = treehmm.resultTF[i][count+5];
        						  treeroot.deTF[i][6] = treehmm.resultTF[i][count+5+(treehmm.splitES[0].length/2)];  
        					  }
        				  }else{
        					  treeroot.deTF = new String[treehmm.resultTF.length][5];
        					  for(int i=0;i<treehmm.resultTF.length;i++){
        						  treeroot.deTF[i][0] = treehmm.resultTF[i][0];
        						  treeroot.deTF[i][1] = treehmm.resultTF[i][1];
        						  treeroot.deTF[i][2] = treehmm.resultTF[i][2];
        						  treeroot.deTF[i][3] = treehmm.resultTF[i][3];
        						  treeroot.deTF[i][4] = treehmm.resultTF[i][4];
        					  }
        				  }
    		}else{
    			giveDEscore(treeroot, treehmm.parent, count+1);
    		}
   }

     public void giveDEscore2(Treenode treeroot, Treenode treehmm, int count) {

    		 if(treeroot.resultTF!=null && treeroot.resultTF[0].length > 6){
   			  List<String> activetf = activeTF.get(theDataSet.dsamplemins[treehmm.ndepth]);
   			  List<String> activeDE = activeDEpvalue.get(theDataSet.dsamplemins[treehmm.ndepth]);
   			  List<String> represstf = repressTF.get(theDataSet.dsamplemins[treehmm.ndepth]);
   			  List<String> repressDE = repressDEpvalue.get(theDataSet.dsamplemins[treehmm.ndepth]);
   			  List<Integer> activeinde = activeindex.get(theDataSet.dsamplemins[treehmm.ndepth]);
   			  List<Integer> repressinde = repressindex.get(theDataSet.dsamplemins[treehmm.ndepth]);
   			  
   				for(int i=0;i<numbits;i++){
   					if(!treeroot.resultTF[i][5+count].equals("NA")){
   						//double enrichpvalue = Double.parseDouble(treeroot.resultTF[i][4]);
   						double enrichpvalue = 1;
   						for(int j=0;j<treehmm.ndepth;j++){
   							Double[] minPvalue = minposipvalue.get(j);
   							if(enrichpvalue>minPvalue[i]){
   								enrichpvalue = minPvalue[i];
   							}
   						}
   						double depvalue = treeroot.splitES[i][count];
   						double defc = treeroot.splitES[i][count+treeroot.splitES[0].length/2];
   						
   						if(depvalue <= dethreshold){
   							if(defc > 0){
   								if(activetf == null){
   									activetf = new ArrayList<String>();
   									activeDE = new ArrayList<String>();  
   									activeinde = new ArrayList<Integer>();  
   									  String tt = enrichpvalue+"\t"+nf4.format(depvalue)+"\t"+nf4.format(defc);
   									  activetf.add(treeroot.resultTF[i][0]);
   									  activeDE.add(tt);
   									  activeinde.add(i);
   								  }else{
   									  if(!activetf.contains(treeroot.resultTF[i][0])){
   									  String tt = enrichpvalue+"\t"+nf4.format(depvalue)+"\t"+nf4.format(defc);
   									  activetf.add(treeroot.resultTF[i][0]);
   									  activeDE.add(tt);
   									  activeinde.add(i);
   									  }
   								  }
   							}else{
   								if(represstf == null){
   									represstf = new ArrayList<String>();
   									repressDE = new ArrayList<String>();  
   									repressinde = new ArrayList<Integer>();  
   									  String tt = enrichpvalue+"\t"+nf4.format(depvalue)+"\t"+nf4.format(defc);
   									  represstf.add(treeroot.resultTF[i][0]);
   									  repressDE.add(tt);
   									  repressinde.add(i);
   								  }else{
   									  if(!represstf.contains(treeroot.resultTF[i][0])){
   									  String tt = enrichpvalue+"\t"+nf4.format(depvalue)+"\t"+nf4.format(defc);
   									  represstf.add(treeroot.resultTF[i][0]);
   									  repressDE.add(tt);
   									  repressinde.add(i);
   									  }
   								  }
   							}
   						}
   					}
   			  }
   			  activeTF.put(theDataSet.dsamplemins[treehmm.ndepth], activetf);
   			  activeDEpvalue.put(theDataSet.dsamplemins[treehmm.ndepth], activeDE);
   			  repressTF.put(theDataSet.dsamplemins[treehmm.ndepth], represstf);
   			  repressDEpvalue.put(theDataSet.dsamplemins[treehmm.ndepth], repressDE);
   			  activeindex.put(theDataSet.dsamplemins[treehmm.ndepth], activeinde);
   			  repressindex.put(theDataSet.dsamplemins[treehmm.ndepth], repressinde);
   		  }

		  if(treehmm.nextptr[0] != null){
			  giveDEscore2(treeroot, treehmm.nextptr[0], count+1);
		  }  
    }
     
     static String[][] BubbleSort(String[][] r, Integer n, Integer col) //升序冒泡排序
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

	/**
	 * Sets up the stored data so that the training data includes all data
	 * except any external held out validation data, meaning the data which
	 * might otherwise used for model selection is just for parameter learning
	 * 
	 * @throws Exception
	 */
	public void combineTrainAndTest() throws Exception {
		nholdout = (int) (theDataSet.data.length * HOLDOUTRATIO); //HOLDOUTRATIO=0
		ntrain = (theDataSet.data.length - nholdout);
		numcols = theDataSet.numcols;

		traindata = new double[ntrain][numcols];

		traingenenames = new String[ntrain];

		int ntrainindex = 0;
		boolean[] toKeep = new boolean[theDataSet.data.length];

		for (int nindex = 0; nindex < theDataSet.data.length; nindex++) {
			if (!bholdout[nindex]) {
				for (int ncol = 0; ncol < numcols; ncol++) {
					traindata[ntrainindex][ncol] = theDataSet.data[nindex][ncol];
				}
				traingenenames[ntrainindex] = theDataSet.genenames[nindex];
				toKeep[nindex] = true;
				ntrainindex++;
			}
		}
		trainBinding = new RegulatorBindingData(bindingData, toKeep, ntrain, true);
		trainBinding.loadInstances(nmaxchild, numbits);
	}

	/**
	 * Splits the non-held out validation data into a training set for learning
	 * the model parameters and a test set for model selection
	 * 
	 * @throws Exception
	 */
	public void splitdata() throws Exception {
		// splits data into training and testsets
		nholdout = (int) (theDataSet.data.length * HOLDOUTRATIO); //HOLDOUTRATIO为剔除基因比例，一般为0，即不剔除基因

		ntrain = (theDataSet.data.length - nholdout);
		ntest = 0;
		numcols = theDataSet.data[0].length;

		// build arrays with time-series training data
		// and binding p-value
		testgenenames = new String[ntest];

		testdata = new double[ntest][numcols];
		traindata = new double[ntrain][numcols];

		double[] foldrandom = new double[theDataSet.data.length];
		double[] foldrandomcopy = new double[theDataSet.data.length];

		// If a set of genes to be held out does not exist yet, create one
		if (bholdout == null) {
			holdoutdata = new Double[nholdout][numcols];
			bholdout = new boolean[theDataSet.data.length];

			holdoutpval = new double[numcols][nholdout][numbits];

			holdoutpvalIndex = new int[numcols][nholdout][];

			// Assign random values (0-1) and store a copy of them
			for (int nindex = 0; nindex < foldrandom.length; nindex++) {
				foldrandom[nindex] = theRandom.nextDouble();
				foldrandomcopy[nindex] = foldrandom[nindex];
			}
			Arrays.sort(foldrandomcopy);
			double dcutoff = foldrandomcopy[holdoutpval[0].length];

			int nholdoutindex = 0;

			for (int nindex = 0; nindex < foldrandom.length; nindex++) {
				// The cutoff was chosen to ensure that the desired
				// number of genes are held out
				if (foldrandom[nindex] < dcutoff) {
					for (int ncol = 0; ncol < numcols; ncol++) {
						holdoutpval[ncol][nholdoutindex] = bindingData.gene2RegBinding[ncol][nindex];
						holdoutpvalIndex[ncol][nholdoutindex] = bindingData.gene2RegBindingIndex[ncol][nindex];
						holdoutdata[nholdoutindex][ncol] = theDataSet.data[nindex][ncol];
					}
					bholdout[nindex] = true;
					nholdoutindex++;
				} else {
					bholdout[nindex] = false;
				}
			}
		}

		/*---------------------------------------------*/
		// random drawing a set of ntrain elements from data.length elements
		/*
		 * The method for drawing the random elements is to maintain a set of
		 * ntrain elements marked as true while iterating through the elements.
		 * The set starts with the initial ntrain elements and adds each element
		 * afterwards with probability (ntrain/nonholdout) swapping a random
		 * element out when this happens.
		 */
		int[] includetrain = new int[ntrain];
		boolean[] btrain = new boolean[theDataSet.data.length];
		if (BDEBUG) {
			System.out.println(ntrain + " " + theDataSet.data.length);
		}

		int nincludeindex = 0;
		int nbtrainindex = 0;
		int nonholdout = 0;
		while (nincludeindex < ntrain) {
			if (!bholdout[nbtrainindex]) {
				includetrain[nincludeindex] = nbtrainindex;
				btrain[nbtrainindex] = true;
				nincludeindex++;
				nonholdout++;
			}
			nbtrainindex++;
		}
		if (BDEBUG) {
			System.out.println("done");
		}

		// drawing nrandom elements from a set of numtotalannotated elements
		// where each element is equally likely
		for (; nbtrainindex < btrain.length; nbtrainindex++) {
			if (!bholdout[nbtrainindex]) {
				if (theRandom.nextDouble() < ((double) ntrain / (double) (nonholdout + 1))) {
					int nreplaceindex = (int) Math.floor(ntrain
							* theRandom.nextDouble());
					btrain[nbtrainindex] = true;
					btrain[includetrain[nreplaceindex]] = false;
					includetrain[nreplaceindex] = nbtrainindex;
				} else {
					btrain[nbtrainindex] = false;
				}
				nonholdout++;
			}
		}
		/*------------------------------------------------------*/
		
		trainBinding = new RegulatorBindingData(bindingData, btrain, ntrain, true);
		trainBinding.loadInstances(nmaxchild, numbits);
		

		int ntrainindex = 0;
		int ntestindex = 0;
		for (int nindex = 0; nindex < btrain.length; nindex++) {
			if (!bholdout[nindex]) {
				if (btrain[nindex]) {
					for (int ncol = 0; ncol < numcols; ncol++) {
						traindata[ntrainindex][ncol] = theDataSet.data[nindex][ncol];
					}
					ntrainindex++;
				} else {
					testgenenames[ntestindex] = theDataSet.genenames[nindex];
					for (int ncol = 0; ncol < numcols; ncol++) {
						testdata[ntestindex][ncol] = theDataSet.data[nindex][ncol];
					}
					ntestindex++;
				}
			}
		}
		// Flip the boolean keep array to create test data
		for (int i = 0; i < btrain.length; i++) {
			if (!bholdout[i])
				btrain[i] = !btrain[i];
		}
		getDEgenes();
	}
	
	public void getDEgenes() throws Exception {
		upgenes.clear();
		downgenes.clear();
		int[][] DE1 = new int[traindata.length][traindata[0].length];
		for(int i=0;i<traindata.length;i++){
			double mean = 0;
			for(int j=1;j<traindata[0].length;j++){
				mean += Math.abs(traindata[i][j]-traindata[i][j-1]);
			}
			mean = mean/(traindata[0].length-1);
			for(int j=1;j<traindata[0].length;j++){
                if(traindata[i][j]-traindata[i][j-1]>mean){
					DE1[i][j] = 1;
				}else{
					if(traindata[i][j]-traindata[i][j-1]<-mean){
						DE1[i][j] = -1;
					}
				}
			}
		}
		for(int j=1; j<traindata[0].length;j++){
			List<Integer> up = new ArrayList<Integer>();
			List<Integer> down = new ArrayList<Integer>();
			for(int i=0;i<traindata.length;i++){
				if(DE1[i][j] == 1){
					up.add(i);
				}else{
					if(DE1[i][j] == -1){
						down.add(i);
					}
				}
			}
			upgenes.add(up);
			downgenes.add(down);
	}
}
	
   public double[][] GetmeanSigma(Treenode node, double[][] sum, int ndepth) throws Exception{
		if((int) node.ndepth == ndepth){
			sum[ndepth][0]=sum[ndepth][0]+node.dsigma;
			sum[ndepth][1]=sum[ndepth][1]+1;
			
		}	
		for (int nchild = 0; nchild < node.numchildren; nchild++) {
		Treenode nextnode = node.nextptr[nchild];
		GetmeanSigma(nextnode, sum, ndepth);
		
	}	
		return sum;
   }
	/**
	 * Parses a file specified from the GUI to obtain the FASTA format of genes
	 * for use by DECOD.
	 * 
	 * @param file
	 * @throws IOException
	 */
	public void parseFastaFile(String file) throws IOException {
		gene2FastaMap = new HashMap<String, String>();
		BufferedReader br = new BufferedReader(new FileReader(file));
		String gene;
		String fastaFormat = "";
		
		String line = br.readLine();
		while(line != null && line.charAt(0)=='#')
			line = br.readLine();
		
		while (line != null) {
			if (line.charAt(0) != '>')
				throw new IOException();
			gene = line.substring(1);
			fastaFormat = line + "\n";
			line = br.readLine();
			while (line != null && line.charAt(0) != '>') {
				fastaFormat += line + "\n";
				
				line = br.readLine();
				while(line != null && line.charAt(0)=='#')
					line = br.readLine();
			}
			gene2FastaMap.put(gene.toUpperCase(), fastaFormat);
		}
	}

	/**
	 * Computes the average and standard deviation expression level at each time
	 * point
	 */
	public void computeDataStats(double[][] data,  double[] dmeans,
			double[] dsigmas) {
		double dsum;
		double dsumsq;
		int npoints;
		double dsigmaavg = 0;
		for (int ncol = 0; ncol < data[0].length; ncol++) {
			dsum = 0;
			dsumsq = 0;
			npoints = 0;
			// need to handle 0 or 1 point error
			for (int nrow = 0; nrow < data.length; nrow++) {
					dsum += data[nrow][ncol];
					dsumsq += Math.pow(data[nrow][ncol], 2);
					npoints++;
			}
			dmeans[ncol] = dsum / npoints;
			dsigmas[ncol] = (dsumsq - Math.pow(dsum, 2) / npoints)/ (npoints - 1);
			dsigmaavg += dsigmas[ncol];
		}
		DEFAULTSIGMA = dsigmaavg / data[0].length;
		if (BDEBUG) {
			System.out.println("[[[]]]" + DEFAULTSIGMA);
			System.out.println();
		}
	}

	// //////////////////////////////////////////////////
	/**
	 * Displays the current temporary TSMiner model, from which a search is
	 * performed for an improved model
	 */
	public void displayTempMap() throws Exception {
		TSMiner_IO.bdisplaycurrent = false;
		final TSMiner_Timeiohmm fthishmm = this;
		
		
		if (BDEBUG) {
			System.out.println("before hmmgui");
		}

		final Treenode treecopy = (Treenode) treeptr.clone();

		computeOrders(treecopy);

		viterbi(treecopy, theDataSet.data);
		if (bindingData.gene2RegBinding != null) {
			runtext.append("================================"+"\n");
			runtext.append("TF enrichment analysis..."+"\n");
			runtext.paintImmediately(runtext.getBounds());	
			
			computeStats(treecopy, treecopy);
			
			runtext.append("\n"+"TF enrichment analysis Done."+"\n");
			runtext.paintImmediately(runtext.getBounds());
			
			computeminparentlevel(treecopy);
		}
		
		traverse(bestTree, 0, false);

		try {
			javax.swing.SwingUtilities.invokeAndWait(new Runnable() {
				public void run() {				
					new TreeGUI(fthishmm, treecopy, npermutationval);
				}
			});
		} catch (InterruptedException iex) {
			iex.printStackTrace(System.out);
		} catch (java.lang.reflect.InvocationTargetException itex) {
			itex.printStackTrace(System.out);
		}
		if (currentButton != null) {
			currentButton.setEnabled(true);
		}
	}

	/**
	 * Executes the first phase of the structure search Considers adding and
	 * deleting paths, does not consider merges or delays
	 */
	public void searchstage1() throws Exception {
		runtext.append("================================"+"\n");
		runtext.append("Add Path..."+"\n");
		runtext.paintImmediately(runtext.getBounds());
		
		String szDIFF = "";
		String szpercent = "";
		int numintervals = traindata[0].length - 1;
		int[] path = new int[numintervals]; //path saving
		traverse(treeptr, 0, false);

		dbestlog = trainhmm (treeptr, true); //train parameters of single chain
		
		traverse(treeptr, 0, false);
		double dprevbestlog;
		dprevouterbestlog = Double.NEGATIVE_INFINITY;
		dprevbestlog = Double.NEGATIVE_INFINITY;
		dbesttrainlike = Double.NEGATIVE_INFINITY;
		boolean bendsearchlocal;
		
		if(runtext != null){
		runtext.append("Current path: 1  ||  Current likelihood: "+dbestlog+"\n");
		runtext.paintImmediately(runtext.getBounds());
		}
		
		do {
			
			dprevouterbestlog = dbestlog;
			traverse(treeptr, 0, false);
			bestTree = treeptr;		
			dprevbestlog = dbestlog;	
				
			String[] TreeModle_TSMiner = viterbi(treeptr,traindata);
			HashMap<String, List<Integer>> Path = new HashMap<String, List<Integer>>();
			for(int i=0;i<TreeModle_TSMiner.length;i++){
					String a = TreeModle_TSMiner[i];
					List<Integer> genes = new ArrayList<Integer>();
					genes = Path.get(a);
					if(genes != null){
						  genes.add(i);
						}else{
							genes = new ArrayList<Integer>();
							genes.add(i);
						}
						Path.put(a, genes);
			}
			
			
                for(Entry<String, List<Integer>> entry:Path.entrySet()){
				
	        	String key = entry.getKey();
	        	HashMap<String, List<List<Integer>>> NCpattern;
	        		
	        		List<Integer> genes = entry.getValue();
	        		
		        	String express[][] = new String[genes.size()][theDataSet.data[0].length+1];
		        	for(int i=0;i<genes.size();i++){
		        		express[i][0] = genes.get(i)+"";
		        		for(int j=0;j<theDataSet.data[0].length;j++){ 
		        			express[i][j+1]=theDataSet.data[genes.get(i)][j]+""; 
		        		}
		        	}
		        	
		        	Negative_correlation NC = new Negative_correlation();
		        	NCpattern = NC.AIE_pattern(express, upgenes, downgenes, Maxcount, minnc, minbic);
		        	
		        	
		        	for(Entry<String, List<List<Integer>>> entry1:NCpattern.entrySet()){
		        		String split = entry1.getKey();
		        		List<List<Integer>> NCgenes = entry1.getValue();
		        		int ndepth = path.length-split.toString().length();
	            		/*
	            		for(int m=0;m<NCgenes.size();m++){
	            			for(int j=0;j<Arrays.asList(key.toString().split(",")).size();j++){
			        			path[j] = Integer.parseInt(Arrays.asList(key.toString().split(",")).get(j));
			        		}
		            		traverseandadd(path, treeptr, ndepth, NCgenes.get(m));
	            		}
	            		*/
	            		
	            		for(int j=0;j<Arrays.asList(key.toString().split(",")).size();j++){
		        			path[j] = Integer.parseInt(Arrays.asList(key.toString().split(",")).get(j));
		        		}
	            		traverseandadd(path, treeptr, ndepth, NCgenes.get(0), NCgenes.get(1));
		        	}
	    	    }
	    	
			treeptr = bestTree; 
			if (dbestlog > dprevbestlog) { //添加路径操作后最佳likehood是否大于之前的
				numtotalPath++;
			}

			while ((dbestlog - dprevbestlog) > 0) {
				if (BDEBUG) {
					System.out.println("trying to delete");
				}
				dprevbestlog = dbestlog;
				traverseanddelete(path, treeptr, treeptr, false, 0, 0);
				
				traverse(bestTree, 0, false);
				treeptr = bestTree;
				
				if (dbestlog > dprevbestlog) {
					if(runtext != null){
						runtext.append("Delete path"+"\n");
						runtext.paintImmediately(runtext.getBounds());
					}
					numtotalPath--;
				}
			}
			
			szpercent = nf3.format(100 * (dbestlog - dprevouterbestlog)/ Math.abs(dbestlog))+ "%";
			szDIFF = nf3.format(dbestlog - dprevouterbestlog);

			bendsearchlocal = TSMiner_IO.bendsearch; //bend searching

			if(runtext != null){
				runtext.append("Current path: "+numtotalPath+"  ||  Current likelihood: "+dbestlog+"\n");
				runtext.paintImmediately(runtext.getBounds());
			}
			
			if (TSMiner_IO.bdisplaycurrent) {
				displayTempMap();
			}
			
		} while ((dbestlog > dprevouterbestlog) && (!bendsearchlocal));	
	}

	// //////////////////////////////////////////////////////////////////////

	/**
	 * Assigns to the nminparentlevel in all nodes the level of the most
	 * immediate ancestor with two or more children
	 */
	public void computeminparentlevel(Treenode ptr) {
		if (ptr != null) {
			if ((ptr.parentptrA == null) || (ptr.parentptrA.length == 1)) {
				if (ptr.parent == null) {
					ptr.nminparentlevel = ptr.ndepth;
					ptr.minparentptr = ptr.parent;
				} else {
					if (ptr.parent.numchildren >= 2) {
						ptr.nminparentlevel = ptr.ndepth;
						ptr.minparentptr = ptr.parent;
					} else {
						ptr.nminparentlevel = ptr.parent.nminparentlevel;
						ptr.minparentptr = ptr.parent.minparentptr;
					}
				}
			} else {
				Treenode[] tempptrA = new Treenode[ptr.parentptrA.length];
				int ndepth = ptr.ndepth;
				for (int nindex = 0; nindex < tempptrA.length; nindex++) {
					tempptrA[nindex] = ptr.parentptrA[nindex];
				}

				boolean bfoundit = false;
				while (!bfoundit) {
					int nj = 1;
					boolean ballsame = true;
					while ((ballsame) && nj < tempptrA.length) {
						if (tempptrA[nj] != tempptrA[nj - 1]) {
							ballsame = false;
						} else {
							nj++;
						}
					}

					ndepth--;
					bfoundit = ballsame;
					for (int nindex = 0; nindex < tempptrA.length; nindex++) {
						tempptrA[nindex] = tempptrA[nindex].parent;
					}
				}
				ptr.nminparentlevel = ndepth;
				ptr.minparentptr = tempptrA[0];
			}

			for (int nindex = 0; nindex < ptr.nextptr.length; nindex++) {
				computeminparentlevel(ptr.nextptr[nindex]);
			}
		}
	}

	/**
	 * Responsible with helper functions for merging paths. Paths must share a
	 * common most recent split and no splits after a merge are allowed
	 */
	public void traverseandmerge(int[] path) throws Exception {

		computeNumLeaves(bestTree);
		if (BDEBUG) {
			System.out.println("in traverse and merge path.length = "
					+ path.length);
		}

		traverse(bestTree, 0, false);

		for (int nparentlevel = path.length - 2; nparentlevel >= 0; nparentlevel--) {
			for (int nmergingdepth = path.length; nmergingdepth > nparentlevel + 1; nmergingdepth--) {
				// nmergingdepth is the depth of the merging
				// parent level is the common ancestor
				Treenode origbestTree;
				do {
					dprevouterbestlog = dbestlog;
					dbestlog = Double.NEGATIVE_INFINITY;

					if (BDEBUG) {
						System.out.println("calling traverseandmergehelp "
								+ nmergingdepth + " " + nparentlevel);
					}
					origbestTree = bestTree;
					traverseandmergehelp(path, 0, nmergingdepth, nparentlevel,
							origbestTree, origbestTree);
					if (dbestlog == Double.NEGATIVE_INFINITY)
						dbestlog = dprevouterbestlog;
				} while (origbestTree != bestTree);
			}
		}

	}

	/**
	 * Helper function for traverseandmerge used to recursively search for
	 * places to merge nodes
	 */
	private void traverseandmergehelp(int[] path, int nlevel,
			int nmergingdepth, int nparentlevel, Treenode ptr, Treenode origroot)
			throws Exception {

		if (nlevel <= nmergingdepth) {
			for (int nindex = 0; nindex < ptr.numchildren; nindex++) {
				path[nlevel] = nindex;
				traverseandmergehelp(path, nlevel + 1, nmergingdepth,
						nparentlevel, ptr.nextptr[nindex], origroot);
			}
		}

		if (nlevel == nmergingdepth) {
			if (BDEBUG) {
				System.out.println("nparentlevel = " + nparentlevel
						+ " nmergindepth = " + nmergingdepth + " nlevel = "
						+ nlevel + " mean = " + ptr.dmean);
			}
			mergenode(path, nmergingdepth, nparentlevel, origroot);
		}
	}

	/**
	 * Computes the number of paths to leaves accessible from each node in the
	 * tree pointed to by ptr
	 */
	public int computeNumLeaves(Treenode ptr) {
		if (ptr.numchildren == 0) {
			return 1;
		} else {
			int nsum = 0;
			for (int nindex = 0; nindex < ptr.numchildren; nindex++) {
				nsum += computeNumLeaves(ptr.nextptr[nindex]);
			}

			if (BDEBUG) {
				System.out.println("ptr.dmean = " + ptr.dmean + " nsum = "
						+ nsum);
			}

			ptr.numdownleaves = nsum;

			return nsum;
		}
	}

	// ////////////////////////////////////////////////////////////////////////////////////////
	/**
	 * Returns true iff all nodes at nleveltogo only have one path to the leaves
	 * Stores in hsLevel all nodes at nleveltogo is 0
	 */
	public boolean findnodesAtLevel(Treenode treeptr, int nleveltogo,
			HashSet<Treenode> hsLevel) {
		if (nleveltogo == 0) {
			hsLevel.add(treeptr);
			if (treeptr.numdownleaves >= 2)
				return false;
			else
				return true;
		} else {
			boolean ballsingles = true;
			for (int nindex = 0; nindex < treeptr.numchildren; nindex++) {
				ballsingles = (ballsingles)
						&& (findnodesAtLevel(treeptr.nextptr[nindex],
								nleveltogo - 1, hsLevel));
			}
			return ballsingles;
		}
	}

	// ////////////////////////////////////////////////////////////////////////////////////////////////
	/**
	 * Helper function of traverseandmergehelp responsible for finding valid
	 * merges and evaluating them that are along path1 and at depth
	 * nmergingdepth
	 */
	private void mergenode(int[] path1, int nmergingdepth, int nparentlevel,
			Treenode origroot) throws Exception {
		double dlog;

		if (BDEBUG) {
			System.out.println("nmerginglevel = " + nmergingdepth);
			System.out.println("nparentlevel = " + nparentlevel);
			System.out.print("Path: ");
			for (int i = 0; i < nparentlevel; i++) {
				System.out.print(path1[i] + " ");
			}
			System.out.println();
		}

		Treenode tempptr = origroot;
		for (int nindex = 0; nindex < nparentlevel; nindex++) {
			tempptr = tempptr.nextptr[path1[nindex]];
		}

		if ((tempptr.numchildren >= 2)) {
			HashSet<Treenode> hsLevel = new HashSet<Treenode>();
			boolean ballsingles = findnodesAtLevel(tempptr,
					(nmergingdepth - nparentlevel), hsLevel);
			// ballsingles - is true if no more splits among nodes being merged

			int numatLevel = hsLevel.size();
			// only binary for now
			if (BDEBUG) {
				System.out.println("numatLevel = " + numatLevel);
			}

			if ((numatLevel == tempptr.numchildren) && (ballsingles)) {
				// not allowing any downstream splits from parent node
				// either before the level being merged, should have been
				// accepted before
				// among merged nodes, not allowing any resplits
				// require ballsingles - to be true
				// also cannot have any splits
				// assume binary for now

				if (BDEBUG) {
					System.out.println("before clone");
				}

				traverse(origroot, 0, false);
				Treenode treeroot = (Treenode) origroot.clone();
				if (BDEBUG) {
					System.out.println("after clone");
				}
				traverse(treeroot, 0, false);

				// traversing to the parent
				tempptr = treeroot;
				for (int nindex = 0; nindex < nparentlevel; nindex++) {
					tempptr = tempptr.nextptr[path1[nindex]];
				}

				Treenode parentptr = tempptr;

				// travering to the nodes for merging
				Treenode splitptr = null;
				for (int nindex = nparentlevel; nindex < nmergingdepth; nindex++) {
					tempptr = tempptr.nextptr[path1[nindex]];
				}

				if (BDEBUG) {
					System.out.println(nmergingdepth + " parent mean "
							+ parentptr.dmean + " tempptr.dmean "
							+ tempptr.dmean + " split mean ");
				}

				boolean bvalidmerge = false;
				ArrayList<Treenode> newparentptrA = null;
				for (int notherchild = 0; notherchild < parentptr.numchildren; notherchild++) {
					if (BDEBUG) {
						System.out.println(notherchild + "\t" + bvalidmerge);
					}

					if (notherchild != path1[nparentlevel]) {
						splitptr = parentptr.nextptr[notherchild];

						for (int nindex = nparentlevel + 1; nindex < nmergingdepth; nindex++) {
							splitptr = splitptr.nextptr[0];
						}

						if (splitptr != tempptr) {
							// its possible two nodes have already been merged
							// on this level
							// but there are two paths two them, we do not need
							// to adjust parent
							// pointers if just taking a different path to the
							// main merging node
							if (!bvalidmerge) {
								// we haven't tried merging yet, just do this
								// once to load the merging node
								// parents into the list
								newparentptrA = new ArrayList<Treenode>();
								if (tempptr.parentptrA == null) {
									// just single parent add to parent list
									newparentptrA.add(tempptr.parent);
									if (BDEBUG) {
										System.out.println("A " + tempptr.dmean
												+ " now links to "
												+ tempptr.parent.dmean);
									}
								} else {
									for (int i = 0; i < tempptr.parentptrA.length; i++) {
										// add all parent links to merging state
										newparentptrA
												.add(tempptr.parentptrA[i]);
										if (BDEBUG) {
											System.out
													.println("B "
															+ tempptr.dmean
															+ " now links to "
															+ tempptr.parentptrA[i].dmean);
										}
									}
								}

								bvalidmerge = true;

								if (BDEBUG) {
									if (tempptr.nextptr[0] != null)
										System.out.println("mean is "
												+ tempptr.nextptr[0].dmean);
									System.out.println("VALID MERGE");
								}
							}

							if (BDEBUG) {
								System.out.println(splitptr.dmean + " ");
							}

							if (splitptr.parentptrA == null) {
								splitptr.parent.nextptr[0] = tempptr;
								newparentptrA.add(splitptr.parent);
								if (BDEBUG) {
									System.out.println("C " + tempptr.dmean
											+ " now links to "
											+ splitptr.parent.dmean);
								}
							} else {
								for (int i = 0; i < splitptr.parentptrA.length; i++) {
									splitptr.parentptrA[i].nextptr[0] = tempptr;
									newparentptrA.add(splitptr.parentptrA[i]);
									if (BDEBUG) {
										System.out.println("D " + tempptr.dmean
												+ " now links to "
												+ splitptr.parentptrA[i].dmean);
									}
								}
							}

							if (tempptr.nextptr[0] != null) {
								tempptr.nextptr[0].parent = tempptr;
								tempptr.nextptr[0].parentptrA = null;
							}

							if (parentptr == splitptr.parent) {
								parentptr.numchildren--;
							}
						}
					}
				}

				if (bvalidmerge) {
					// copying over new parent nodes
					int nsize = newparentptrA.size();
					tempptr.parentptrA = new Treenode[nsize];
					for (int i = 0; i < nsize; i++) {
						tempptr.parentptrA[i] = newparentptrA.get(i);
					}
					if (BDEBUG) {
						System.out.println("MERGED TREE");
					}
					traverse(treeroot, 0, false);

					dlog = trainhmm(treeroot, true);

					if (BDEBUG) {
						System.out.println("after training:");
					}

					traverse(treeroot, 0, false);

					double dimprovemerge = (dlog - dprevouterbestlog)
							/ Math.abs(dlog);
					if (BDEBUG) {
						System.out.println("merge* " + dimprovemerge + "\t"
								+ MERGEMIN + "\t" + dlog + "\t"
								+ dprevouterbestlog + "\t" + dbestlog);
					}

					if ((dlog >= dprevouterbestlog) && (dlog > dbestlog)) {
						bestTree = treeroot;
						computeNumLeaves(bestTree);
						dbestlog = dlog;
						if (BDEBUG) {
							System.out.println("in merge " + dbestlog
									+ " best tree is");
							traverse(treeroot, 0, false);
						}

						double dimprove = (dbestlog - dprevouterbestlog)
								/ Math.abs(dbestlog);
						String szimprove = nf3.format(100 * dimprove) + "%";
						String szimproveDIFF = nf3.format(dbestlog
								- dprevouterbestlog);
					}
				}
			}
		}
	}

	// /////////////////////////////////////////////////////////////////////
	/**
	 * A record containing a node id, mean, and standard deviation used for
	 * ordering the children node
	 */
	static class OrderRec {
		int nid;
		double dmean;
		double dsigma;

		OrderRec(double dmean, double dsigma, int nid) {
			this.nid = nid;
			this.dmean = dmean;
			this.dsigma = dsigma;
		}
	}

	/**
	 * Comparator for OrderRec sorts first on dmean, then dsigma, then id all in
	 * increasing order.
	 */
	static class OrderRecCompare implements Comparator<OrderRec> {
		public int compare(OrderRec or1, OrderRec or2) {
			if (or1.dmean < or2.dmean)
				return -1;
			else if (or1.dmean > or2.dmean)
				return 1;
			else if (or1.dsigma < or2.dsigma)
				return -1;
			else if (or1.dsigma > or2.dsigma)
				return 1;
			else if (or1.nid < or2.nid)
				return -1;
			else if (or1.nid > or2.nid)
				return 1;
			else
				return 0;
		}
	}

	// /////////////////////////////////////////////////////////////////////
	/**
	 * Sorts the children of each node in treeptr and descendants based on
	 * OrderRecCompare
	 */
	public void computeOrders(Treenode treeptr) {
		if (treeptr.numchildren >= 1) {
			if (treeptr.numchildren >= 2) {
				OrderRec[] theOrderRecs = new OrderRec[treeptr.numchildren];
				for (int nindex = 0; nindex < treeptr.numchildren; nindex++) {
					theOrderRecs[nindex] = new OrderRec(
							treeptr.nextptr[nindex].dmean,
							treeptr.nextptr[nindex].dsigma, nindex);
				}
				Arrays.sort(theOrderRecs, new OrderRecCompare());
				treeptr.orderA = new int[theOrderRecs.length];
				
				for (int nindex = 0; nindex < treeptr.orderA.length; nindex++) {
					treeptr.orderA[nindex] = theOrderRecs[nindex].nid;
					//System.out.println(treeptr.ndepth+" "+nindex+" "+treeptr.orderA[nindex]);
				}
			}else{
				treeptr.orderA = new int[1];
				treeptr.orderA[0] = 0;
			}
			for (int nindex = 0; nindex < treeptr.numchildren; nindex++) {
				computeOrders(treeptr.nextptr[nindex]);
			}
		}
	}

	// /////////////////////////////////////////////////////////////////////

	/**
	 * Builds the initial tree which is just a single chain with mean and
	 * standard deviation for each node being the global mean and standard
	 * deviation at the time point
	 */
	public void buildEmptyTree(int ncurrdepth, Treenode currnode,
			double[] dmean, double[] dsigma) {
		currnode.dsigma = Math.max(Math.max(dsigma[ncurrdepth], dminstddev),MINSIGMA);
		if (ncurrdepth + 1 < dsigma.length) {
			currnode.numchildren = 1;
			currnode.nextptr[0] = new Treenode(currnode);
			currnode.nextptr[0].dmean = dmean[ncurrdepth+1];
			currnode.nextptr[0].parent = currnode;
			if (BREGTSMiner) {
				currnode.ptrans[0] = 1;
			}
			buildEmptyTree(ncurrdepth + 1, currnode.nextptr[0], dmean, dsigma);
		}
	}

	/**
	 * Initializes the model parameters to those in szinitfileval
	 */
	public void readInitTree(Treenode currnode) throws IOException {

		BufferedReader brinit = new BufferedReader(
				new FileReader(szinitfileval));
		String szfirstline = brinit.readLine();
		boolean bmergedtree = false;

		if (szfirstline.equalsIgnoreCase("MERGED TREE")) {
			szfirstline = brinit.readLine();
			bmergedtree = true;
		}

		if (bmergedtree) {
			if (ninitsearchval == 1) {
				throw new IllegalArgumentException(
						"Saved model has merges, but option allowed to start "
								+ "search from model is selected.\n"
								+ "Search may not begin from model with merges.  "
								+ "Consider option to use the saved model as is.");
			}

			if (!ballowmergeval) {
				throw new IllegalArgumentException(
						"Saved model has merges, but merges are not selected to be allowed.");
			}
		}

		StringTokenizer stfirstline = new StringTokenizer(szfirstline, "\t");
		if (stfirstline.countTokens() < 2) {
			throw new IllegalArgumentException("Not a valid saved model file!");
		}

		String szfirsttoken = stfirstline.nextToken();
		int numcoeff = Integer.parseInt(stfirstline.nextToken());
		if ((numcoeff != bindingData.regNames.length)
				&& (!szfirsttoken.startsWith("REGTSMiner"))) {
			throw new IllegalArgumentException(
					"Number of coefficients ("
							+ numcoeff
							+ ") in the saved model does not "
							+ "match the number of coefficients in the static input file ("
							+ bindingData.regNames.length + ")");
		}

		if (!(szfirsttoken.startsWith("REGTSMiner")) && (BREGTSMiner)) {
			throw new IllegalArgumentException(
					"'Use Transcription Factor-gene Interaction Data to Build Model' under 'Main Search Options' is unselected, "
							+ "but saved model was built using the data.");
		} else if ((szfirsttoken.startsWith("REGTSMiner")) && (!BREGTSMiner)) {
			throw new IllegalArgumentException(
					"'Use Transcription Factor-gene Interaction Data to Build Model' under 'Main Search Options' is selected, "
							+ "but saved model was built without the data.");
		}

		int ndepth = readInitTreeHelp(currnode, brinit, numcoeff, bmergedtree);

		if (ndepth != numcols) {
			throw new IllegalArgumentException(
					"Number of time points in the initial model ("
							+ ndepth
							+ ") does not match number of time points in the data ("
							+ numcols + ")");
		}

		String szLine = brinit.readLine();
		if ((szLine != null) && (szLine.equals("COLORS"))) {
			readColors(brinit);
		}

		if (bmergedtree) {
			HashMap<Integer, Treenode> createdNodes = new HashMap<Integer, Treenode>();
			mergeDuplicates(currnode, createdNodes);
		}
	}

	/**
	 * Helper function for reading the model parameters from a file
	 */
	private int readInitTreeHelp(Treenode currnode, BufferedReader brinit,
			int numcoeff, boolean bmergedtree) throws IOException {
		String szLine = brinit.readLine();
		StringTokenizer st = new StringTokenizer(szLine, "\t");
		int nID;

		if (bmergedtree) {
			nID = Integer.parseInt(st.nextToken());
			currnode.nID = nID;
		}

		currnode.dmean = Double.parseDouble(st.nextToken());
		currnode.dsigma = Double.parseDouble(st.nextToken());
		currnode.numchildren = Integer.parseInt(st.nextToken());

		int numcoefftoread;

		if (currnode.numchildren <= 1) {
			numcoefftoread = 0;
			currnode.binit = false;
		} else if (currnode.numchildren == 2) {
			numcoefftoread = numcoeff + 1;
			currnode.binit = true;
			numtotalPath++;
		} else {
			
			if (currnode.numchildren > nmaxchild) {
				
				throw new IllegalArgumentException(
						"Number of paths out of a node in the initial model file ("
								+ currnode.numchildren
								+ ") exceeds the maximum allowed " + nmaxchild);

			}
			
			numtotalPath += currnode.numchildren - 1;
			numcoefftoread = (numcoeff + 1) * (currnode.numchildren - 1);
			currnode.binit = true;
		}

		if ((numcoefftoread == 0) && (BREGTSMiner)) {
			currnode.ptrans[0] = 1;
		} else if (numcoefftoread > 0) {
			if (BREGTSMiner) {
				for (int nindex = 0; nindex < currnode.numchildren; nindex++) {
					currnode.ptrans[nindex] = Double.parseDouble(brinit.readLine());
				}
			} else {
				double[] dcoeff = new double[numcoefftoread];
				String szcoeff;
				StringTokenizer stcoeff;
				for (int ncoeffindex = 0; ncoeffindex < numcoefftoread; ncoeffindex++) {
					szcoeff = brinit.readLine();
					stcoeff = new StringTokenizer(szcoeff, "\t");
					String szcoeffname = stcoeff.nextToken();
					int nmodrow = ncoeffindex % (bindingData.regNames.length + 1);
					
					if ((nmodrow != 0) && (!szcoeffname.equals(bindingData.regNames[nmodrow - 1]))) {
						throw new IllegalArgumentException(ncoeffindex + " "
								+ (bindingData.regNames.length + 1) + " "
								+ "Transcription factor " + szcoeffname
								+ " in saved model  "
								+ "does not match transcription factor "
								+ bindingData.regNames[nmodrow - 1] + ".");
					}
					dcoeff[ncoeffindex] = Double.parseDouble(stcoeff.nextToken());
				}

				currnode.tranC = new TSMiner_FastLogistic2(
						trainBinding.theInstancesIndex[currnode.ndepth][currnode.numchildren - 2],
						trainBinding.theInstances[currnode.ndepth][currnode.numchildren - 2],
						trainBinding.theInstancesRegIndex[currnode.ndepth][currnode.numchildren - 2],
						trainBinding.theInstancesReg[currnode.ndepth][currnode.numchildren - 2],
						trainBinding.ylabels[currnode.ndepth][currnode.numchildren - 2],
						currnode.numchildren, dcoeff, numbits,
						trainBinding.regTypes);
			}
		}

		if (currnode.numchildren == 0) {
			return 1;
		} else {
			int ndepth;
			currnode.nextptr[0] = new Treenode(currnode);
			currnode.nextptr[0].parent = currnode;
			ndepth = readInitTreeHelp(currnode.nextptr[0], brinit, numcoeff,
					bmergedtree);
			for (int nchild = 1; nchild < currnode.numchildren; nchild++) {
				currnode.nextptr[nchild] = new Treenode(currnode);
				currnode.nextptr[nchild].parent = currnode;
				if (ndepth != readInitTreeHelp(currnode.nextptr[nchild],
						brinit, numcoeff, bmergedtree)) {
					throw new IllegalArgumentException(
							"Invalid saved model file.  Two paths are not of the same length.");
				}
			}
			return (ndepth + 1);
		}
	}

	/**
	 * When reading a model with merges from a file some of the states can be
	 * listed more than once this procedures makes it so that in the internal
	 * representation the state only exists once
	 */
	private void mergeDuplicates(Treenode currnode,
			HashMap<Integer, Treenode> createdNodes) {
		if (currnode != null) {
			createdNodes.put(new Integer(currnode.nID), currnode);
			for (int nchild = 0; nchild < currnode.numchildren; nchild++) {
				Treenode tnode = createdNodes.get(new Integer(
						currnode.nextptr[nchild].nID));
				if (tnode == null) {
					mergeDuplicates(currnode.nextptr[nchild], createdNodes);
				} else {
					currnode.nextptr[nchild] = tnode;
					if (tnode.parentptrA == null) {
						tnode.parentptrA = new Treenode[2];
						tnode.parentptrA[0] = tnode.parent;
						tnode.parentptrA[1] = currnode;
					} else {
						Treenode[] temptr = tnode.parentptrA;
						tnode.parentptrA = new Treenode[tnode.parentptrA.length];
						for (int nindex = 0; nindex < temptr.length; nindex++) {
							tnode.parentptrA[nindex] = temptr[nindex];
						}
						tnode.parentptrA[temptr.length - 1] = currnode;
					}
				}
			}
		}
	}

	/**
	 * Loads saved colors for a TSMiner map
	 */
	public void readColors(BufferedReader brinit) throws IOException {
		savedColors = new ArrayList<Color>();
		String szLine;
		while ((szLine = brinit.readLine()) != null) {
			if (!szLine.trim().equals("")) {
				StringTokenizer st = new StringTokenizer(szLine, "\t");
				if (st.countTokens() != 4) {
					throw new IllegalArgumentException("A color line has "
							+ st.countTokens() + " values expecting 4");
				} else {
					float f1, f2, f3, f4;
					f1 = Float.parseFloat(st.nextToken());
					f2 = Float.parseFloat(st.nextToken());
					f3 = Float.parseFloat(st.nextToken());
					f4 = Float.parseFloat(st.nextToken());
					savedColors.add(new Color(f1, f2, f3, f4));
				}
			}
		}
	}

	

	// ////////////////////////////////////////////////////////////////////
	/**
	 * Deletes a child from the specified path on a cloned version of the root
	 * node
	 */
	public Treenode deletepath(int[] path, int ndesiredlevel, Treenode root) {
		// used to be call deleteleaf
		Treenode treeroot = (Treenode) root.clone();
		Treenode ptr = treeroot;

		for (int nindex = 0; nindex < ndesiredlevel - 1; nindex++) {
			ptr = ptr.nextptr[path[nindex]];
		}

		if (!ptr.nextptr[path[ndesiredlevel - 1]].bchange) {
			return null;
		} else {
			for (int nj = path[ndesiredlevel - 1] + 1; nj < ptr.numchildren; nj++) {
				ptr.nextptr[nj - 1] = ptr.nextptr[nj];
			}

			ptr.numchildren--;
			ptr.binit = false;

			if (BREGTSMiner) {
				double dsum = 0;
				for (int nindex = 0; nindex < ptr.numchildren; nindex++) {
					dsum += ptr.ptrans[nindex];
				}

				for (int nindex = 0; nindex < ptr.numchildren; nindex++) {
					ptr.ptrans[nindex] /= dsum;
				}
			}
		}
		return treeroot;
	}

	// //////////////////////////////////////////////////////////////////////
	/**
	 * With its helper functions adds the specified path to the model
	 */
	public void traverseandadd(int[] path, Treenode origroot, int ndepth, List<Integer> list1, List<Integer> list2)
			throws Exception {
		traverseandaddhelp(path, ndepth, origroot, list1, list2);
	}

	// ////////////////////////////////////////////////////////////////////

	/**
	 * Helper function for traverseandaddhelp recursively tries to add paths to
	 * the model, and then evaluates the improvement
	 */
	private void traverseandaddhelp(int[] path, int nlevel, Treenode origroot, List<Integer> list1,
			List<Integer> list2) 
			throws Exception {
		Treenode treeroot = (Treenode) origroot.clone(); //clone root node
		Treenode ptr = treeroot;
		for (int nindex = 0; nindex < nlevel; nindex++) {
			ptr = ptr.nextptr[path[nindex]]; //to the split node
		}
		
		if (nlevel < path.length && !TSMiner_IO.bendsearch && (ptr != null) && (ptr.numchildren < nmaxchild)) {
				
			if (TSMiner_IO.bdisplaycurrent) {
				displayTempMap();
			}
			
			if (BDEBUG) {
					for (int ni = 0; ni < nlevel; ni++) {
						System.out.print(path[ni] + "\t");
					}
				}
				
				Treenode splittree = splitnode(path, nlevel, origroot, list1, list2);  //在nlevel层面添加分岔

				if (BDEBUG) {
					System.out.println(path[nlevel]);
				}

				if (splittree != null) {   //IF succeed to add path
					traverse(splittree, 0, false);

					double dlog = trainhmm(splittree, false); 

					traverse(splittree, 0, false);
					
					if (BDEBUG) {
						System.out.println("*test " + dlog);
					}
					
					if (dlog > dbestlog) {
						bestTree = splittree;
						dbestlog = dlog;
						dbesttrainlike = dtrainlike; //best train score after adding path without penalty
					}
				}
		}
	}

	// ///////////////////////////////////////////////////////////////////
	/**
	 * Helper function for adding a path in traverseandaddhelp. Responsible for
	 * splitting one path out of a node into two.
	 */
	private Treenode splitnode(int[] path, int ndesiredlevel, Treenode origroot, List<Integer> list1,
			 List<Integer> list2) {

		Treenode treeroot = (Treenode) origroot.clone();
		Treenode ptr1 = treeroot;
		Treenode ptr2 = treeroot;
		
		for (int nindex = 0; nindex < ndesiredlevel; nindex++) {
			ptr1 = ptr1.nextptr[path[nindex]];
			ptr2 = ptr2.nextptr[path[nindex]];
		}

		if (!ptr1.bchange && !ptr2.bchange) { //状态值不可变则不能添加
			return null;
		} else {
			double[][] mean_sigma = Mean_sigma(path, ndesiredlevel, origroot, list1, list2);

				for (int nlevel = ndesiredlevel; nlevel < path.length; nlevel++) {
					ptr2.nextptr[path[nlevel]].dmean = mean_sigma[0][nlevel];
					ptr2.nextptr[path[nlevel]].dsigma = mean_sigma[1][nlevel];
					ptr2 = ptr2.nextptr[path[nlevel]];
				}

				// increasing the number of children			
				ptr1.numchildren++; //ptr的子节点数目+1
				ptr1.binit = false;

				/*
				if (BREGTSMiner) {
					ptr1.ptrans[path[ndesiredlevel]] = 0.5*(ptr1.ptrans[path[ndesiredlevel]]);
					ptr1.ptrans[ptr1.numchildren - 1] = ptr1.ptrans[path[ndesiredlevel]];
					
					if (BDEBUG) {
						System.out.println("probs A " + ptr1.ptrans[0] + "\t" + ptr1.ptrans[1]);
					}
				}
				*/
				if (BREGTSMiner) { //TRUE
					for (int nindex = 0; nindex < ptr1.numchildren - 1; nindex++) {
						ptr1.ptrans[nindex] = ptr1.ptrans[nindex] * (ptr1.numchildren - 1) / (ptr1.numchildren);
					}
					
					ptr1.ptrans[ptr1.numchildren - 1] = 1.0 / (ptr1.numchildren);

					if (BDEBUG) {
						System.out.println("probs A " + ptr1.ptrans[0] + "\t" + ptr1.ptrans[1]);
					}
				}

				if (ndesiredlevel < path.length) {
					path[ndesiredlevel] = ptr1.numchildren - 1;
				}
				
				///////////////////////////////////////////////////////////
				//more sophisticated initialization schemes possible here//
				///////////////////////////////////////////////////////////
					for (int nlevel = ndesiredlevel; nlevel < path.length; nlevel++) {  
						ptr1.nextptr[ptr1.numchildren - 1] = new Treenode(ptr1); //生成ptr1子节点
						ptr1.nextptr[ptr1.numchildren - 1].dmean = mean_sigma[2][nlevel];
						ptr1.nextptr[ptr1.numchildren - 1].dsigma = mean_sigma[3][nlevel];
						
						ptr1 = ptr1.nextptr[ptr1.numchildren - 1];
						ptr1.numchildren = 1;
						if (BREGTSMiner) {
							ptr1.ptrans[0] = 1;
						}
					}
				ptr1.numchildren = 0; //terminate
			}
		return treeroot;
	}

	private double[][] Mean_sigma(int[] path, int ndesiredlevel, Treenode treeroot, List<Integer> list1
			, List<Integer> list2) {
		double[][] mean_sigma = new double[4][path.length];
		
		Treenode ptr = treeroot;
		for (int nindex = 0; nindex < ndesiredlevel; nindex++) {
			ptr = ptr.nextptr[path[nindex]];
		}
		
		for(int i=ndesiredlevel;i<path.length;i++){
			ptr = ptr.nextptr[path[i]]; 
			List<Integer> genelist = ptr.genelist;
			double npoints1=0;
			double dsum1=0;
			double dsumsq1=0;
			double npoints2=0;
			double dsum2=0;
			double dsumsq2=0;
			
			for(int j=0;j<genelist.size();j++){
				if(!list2.contains(genelist.get(j))){
					dsum1 += traindata[genelist.get(j)][i+1];
					dsumsq1 += Math.pow(traindata[genelist.get(j)][i+1], 2);
					npoints1++;
				}
				if(!list1.contains(genelist.get(j))){
					dsum2 += traindata[genelist.get(j)][i+1];
					dsumsq2 += Math.pow(traindata[genelist.get(j)][i+1], 2);
					npoints2++;
				}
			}
			
			/*
			for(int j=0;j<list1.size();j++){
				dsum1 += traindata[list1.get(j)][i+1];
				dsumsq1 += Math.pow(traindata[list1.get(j)][i+1], 2);
				npoints1++;
			}
			for(int j=0;j<list2.size();j++){
				dsum2 += traindata[list2.get(j)][i+1];
				dsumsq2 += Math.pow(traindata[list2.get(j)][i+1], 2);
				npoints2++;
			}
            */
			mean_sigma[0][i] = dsum1 / npoints1;
			mean_sigma[1][i] = Math.sqrt((dsumsq1 - Math.pow(dsum1, 2) / npoints1) / (npoints1 - 1));
			mean_sigma[2][i] = dsum2 / npoints2;
			mean_sigma[3][i] = Math.sqrt((dsumsq2 - Math.pow(dsum2, 2) / npoints2) / (npoints2 - 1));
		}
		return mean_sigma;	
	}
	// ////////////////////////////////////////////////////////////////////
	/**
	 * With its helper function finds the path with the fewest genes assigned
	 * this will become a candidate for deletion
	 */
	public MinPathRec traverseanddeleteMinPath(int[] path, int[] bestpath,
			Treenode ptr) throws Exception {
		return traverseanddeletehelpMinPath(path, 0, ptr);
	}

	// ////////////////////////////////////////////////////////////////////
	/**
	 * A record used for storing the path with the fewest genes assigned
	 */
	static class MinPathRec {
		int nval;
		int nlevel;
		int[] bestpath;

		MinPathRec(int nval, int nlevel, int[] bestpath) {
			this.nval = nval;
			this.nlevel = nlevel;
			this.bestpath = bestpath;
		}
	}

	/**
	 * Deletes the specificed path from the model starting from ndesiredlevel
	 */
	public void deleteMinPath(int[] path, int ndesiredlevel, Treenode root) {
		Treenode ptr = root;

		for (int nindex = 0; nindex < ndesiredlevel - 1; nindex++) {
			ptr = ptr.nextptr[path[nindex]];
		}

		for (int nj = path[ndesiredlevel - 1] + 1; nj < ptr.numchildren; nj++) {
			ptr.nextptr[nj - 1] = ptr.nextptr[nj];
		}

		ptr.numchildren--;
		ptr.binit = false;

		if (BREGTSMiner) {
			double dsum = 0;
			for (int nindex = 0; nindex < ptr.numchildren; nindex++) {
				dsum += ptr.ptrans[nindex];
			}

			for (int nindex = 0; nindex < ptr.numchildren; nindex++) {
				if (dsum > 0) {
					ptr.ptrans[nindex] /= dsum;
				} else {
					ptr.ptrans[nindex] = 0;
				}
			}
		}
	}

	/**
	 * Recursively searches for the path with the fewest genes assigned this
	 * will become a candidate for deletion
	 */
	private MinPathRec traverseanddeletehelpMinPath(int[] path, int nlevel,
			Treenode ptr) throws Exception {
		int nmin = Integer.MAX_VALUE;
		if (nlevel < path.length) {
			MinPathRec theminMinPathRec = null;
			for (int nindex = 0; nindex < ptr.numchildren; nindex++) {
				path[nlevel] = nindex;
				MinPathRec thetempMinPathRec = traverseanddeletehelpMinPath(
						path, nlevel + 1, ptr.nextptr[nindex]);

				if (thetempMinPathRec.nval < nmin) {
					nmin = thetempMinPathRec.nval;
					theminMinPathRec = thetempMinPathRec;
					theminMinPathRec.bestpath[nlevel] = nindex;
					if (thetempMinPathRec.nval == ptr.numPath) {
						thetempMinPathRec.nlevel = nlevel;
					}
				}
			}
			return theminMinPathRec;
		} else {
			MinPathRec theMinPathRec = new MinPathRec(ptr.numPath, nlevel,
					new int[path.length]);
			return theMinPathRec;
		}
	}

	// //////////////////////////////////////////////////////////////////////////
	public boolean traverseanddelay(int[] path, Treenode ptr,
			int ndesiredlevel, Treenode origroot)
			throws Exception {
		if (BDEBUG) {
			System.out.println("entering traverse and delay");
		}
		return traverseanddelayhelp(path, 0, ndesiredlevel, ptr, origroot, 0);
	}

	// ///////////////////////////////////////////////////////////////////////////
	/**
	 * Recursively searches for splits to delay one time point
	 */
	private boolean traverseanddelayhelp(int[] path, int nlevel,
			int ndesiredlevel, Treenode ptr, Treenode origroot,
			int nchild) throws Exception {

		if (BDEBUG) {
			System.out.println("traverseanddelayhelp " + ndesiredlevel + " "
					+ ndesiredlevel);
		}

		boolean breturnval = false;
		if ((nlevel < path.length - 1) && (nlevel < ndesiredlevel)) {
			for (int nindex = 0; nindex < ptr.numchildren; nindex++) {
				if (TSMiner_IO.bdisplaycurrent) {
					displayTempMap();
				}
				
				path[nlevel] = nindex;
				boolean breturnvaltemp;
				breturnvaltemp = traverseanddelayhelp(path, nlevel + 1,
						ndesiredlevel, ptr.nextptr[nindex], origroot,
						nindex);
				breturnval = breturnval || breturnvaltemp;
			}
		}

		// assuming ptr.parent also not null
		if ((ptr != null) && (nchild >= 1) && (nlevel == ndesiredlevel)) {
			for (int nchildcombine = 0; nchildcombine < nchild; nchildcombine++) {
				if (ptr.numchildren
						+ ptr.parent.nextptr[nchildcombine].numchildren <= nmaxchild) {
					double dabsmeandiff = Math.abs(ptr.dmean
							- ptr.parent.nextptr[nchildcombine].dmean);
					NumberFormat nf4 = NumberFormat.getInstance();
					nf4.setMinimumFractionDigits(4);
					nf4.setMaximumFractionDigits(4);
					if (BDEBUG) {
						System.out
								.println(nlevel
										+ " "
										+ nf4.format(ptr.dmean)
										+ " "
										+ nf4
												.format(ptr.parent.nextptr[nchildcombine].dmean)
										+ " "
										+ nf4.format(dabsmeandiff)
										+ " "
										+ nf4.format(ptr.dsigma)
										+ " "
										+ nf4
												.format(ptr.parent.nextptr[nchildcombine].dsigma));
						System.out.println("-----------");
					}

					// mark to try to delete
					if (BDEBUG) {
						System.out.print("dy ");
						for (int ni = 0; ni < nlevel; ni++) {
							System.out.print(path[ni] + " ");
						}
						System.out.println();
						System.out.println("before");
					}

					traverse(origroot, 0, false);
					Treenode delaytree = delaysplit(path, nlevel, origroot,
							nchildcombine);

					if (delaytree != null) {
						double dlog = trainhmm(delaytree, true);

						if (BDEBUG) {
							System.out.println("after " + dlog);
						}
						traverse(delaytree, 0, false);

						if (BDEBUG) {
							System.out.println("test " + dlog);
						}

						if ((dlog >= dprevouterbestlog) && (dlog > dbestlog)) {
							bestTree = delaytree;
							dbestlog = dlog;
							if (BDEBUG) {
								System.out.println("in delay " + dbestlog);
							}

							dbesttrainlike = dtrainlike;
							breturnval = true;
						}
					}
				}
			}
		}
		return breturnval;
	}
	/**
	 * Helper function to traverseanddelayhelp that manipulates a split so that
	 * it occurs one time point later
	 */
	private Treenode delaysplit(int[] path, int ndesiredlevel, Treenode root,
			int nchildcombine) {
		Treenode treeroot = (Treenode) root.clone();
		Treenode ptr = treeroot;

		for (int nindex = 0; nindex < ndesiredlevel - 1; nindex++) {
			ptr = ptr.nextptr[path[nindex]];
		}
		Treenode removeptr = ptr.nextptr[path[ndesiredlevel - 1]];

		if (!removeptr.bchange) {
			return null;
		} else {
			for (int nj = path[ndesiredlevel - 1] + 1; nj < ptr.numchildren; nj++) {
				ptr.nextptr[nj - 1] = ptr.nextptr[nj];
			}

			ptr.numchildren--;
			ptr.binit = false;

			if (BREGTSMiner) {
				double dsum = 0;
				for (int nindex = 0; nindex < ptr.numchildren; nindex++) {
					dsum += ptr.ptrans[nindex];
				}

				for (int nindex = 0; nindex < ptr.numchildren; nindex++) {
					ptr.ptrans[nindex] /= dsum;
				}
			}

			int ncurrnumchild = ptr.nextptr[nchildcombine].numchildren;

			Treenode combineptr = ptr.nextptr[nchildcombine];
			for (int nindex = 0; nindex < removeptr.numchildren; nindex++) {
				combineptr.nextptr[nindex + ncurrnumchild] = removeptr.nextptr[nindex];
				removeptr.nextptr[nindex].parent = combineptr;
			}

			if (BREGTSMiner) {
				for (int nindex = 0; nindex < combineptr.numchildren; nindex++) {
					combineptr.ptrans[nindex] = combineptr.ptrans[nindex] * 0.5;
				}

				for (int nindex = 0; nindex < removeptr.numchildren; nindex++) {
					combineptr.ptrans[nindex + combineptr.numchildren] = removeptr.ptrans[nindex] * 0.5;
				}
			}
			combineptr.numchildren += removeptr.numchildren;
			combineptr.binit = false;
		}
		return treeroot;
	}
	// ////////////////////////////////////////////////////////////////////
	/**
	 * With its helper function searches for the best path to delete Returns
	 * true if an eligible path to delete was found
	 */
	public boolean traverseanddelete(int[] path, Treenode ptr,
			Treenode origroot, boolean bresplit, double dimprovemin,
			double dimprovemindiff) throws Exception {
		return traverseanddeletehelp(path, 0, ptr, origroot, bresplit,
				dimprovemin, dimprovemindiff);
	}

	/**
	 * Helper function that searches for the best path to delete
	 */
	private boolean traverseanddeletehelp(int[] path, int nlevel, Treenode ptr,
			Treenode origroot, boolean bresplit, double dimprovemin,
			double dimprovemindiff) throws Exception {

		boolean breturnval = false;
		if (nlevel < path.length) {
			for (int nindex = 0; nindex < ptr.numchildren; nindex++) {
				if (TSMiner_IO.bdisplaycurrent) {
					displayTempMap();
				}
				
				path[nlevel] = nindex;
				boolean breturnvaltemp = traverseanddeletehelp(path,
						nlevel + 1, ptr.nextptr[nindex], origroot, bresplit,
						dimprovemin, dimprovemindiff);
				breturnval = (breturnvaltemp || breturnval);  //遍历所有节点，只要一个节点满足删除条件，为true
			}
		}

		if ((ptr != null) && (ptr.parent != null) && (ptr.parent.numchildren > 1)) {
			if (BDEBUG) {
				// mark to try to delete
				for (int ni = 0; ni < nlevel; ni++) {
					System.out.print(path[ni] + " ");
				}
				System.out.println();
			}

			Treenode deletetree = deletepath(path, nlevel, origroot);

			if (deletetree != null) {
				double dlog = trainhmm(deletetree, false);

				if (BDEBUG) {
					System.out.println("test " + dlog);
				}

				double dimprovedelete = (dlog - dprevouterbestlog) / Math.abs(dlog);
				if (BDEBUG) {
					System.out.println("del* " + dimprovedelete + "\t"
							+ dimprovemin + "\t" + (dlog - dprevouterbestlog));
					System.out.println("~~" + dlog + "\t" + dbestlog + "\t"
							+ dprevouterbestlog);

					System.out.println("DELETE\t" + dlog + "\t" + dbestlog
							+ "\t" + dimprovedelete + "\t"
							+ (dlog - dprevouterbestlog) / Math.abs(dlog)
							+ "\t" + dimprovemin);
				}

				if ((dimprovemin >= dprevouterbestlog) && (dlog > dbestlog))
				{
					bestTree = deletetree;
					dbestlog = dlog;
					if (BDEBUG) {
						System.out.println("in delete " + dbestlog);
					}
					dbesttrainlike = dtrainlike;
					breturnval = true;
				}
			}
		}
		return breturnval;
	}

	 public static void getnode(Treenode ptr, Treenode selectednode, int[] path, int ndepth, int currentdepth) {
			if(currentdepth < ndepth){
				int nchild = path[currentdepth+1];
				getnode(ptr.nextptr[nchild], selectednode, path, ndepth, currentdepth+1);
			}else if(currentdepth == ndepth){
				selectednode = ptr;
			}
		}
	// ////////////////////////////////////////////////////////////////////
	/**
	 * A Treenode corresponds to a state in the model
	 */
	public class Treenode {
		int ncurrtime;// use to count visited
		int nID;
		int niteration;
		int ninstanceiteration;
		int nrec;
		int ninstancerec;
		int nminparentlevel;
		Treenode minparentptr;
		public int[] heatmaprange;

		Treenode[] parentptrA;
		/** number of leaves descendants of this node, only relevant if split */
		int numdownleaves;
		double dens;
		/** expected sum of the emissions from this node */
		double dEsum;
		/** expected sum of the square of emissions from this node */
		double dEsumsq;
		double dPsum;
		List <Integer> genelist = new ArrayList<Integer>(); 
		/** emission mean */
		double dmean;
		/** emission std */
		double dsigma;
		int[] path;
		double dpredictweight;
		String szgolabel = "";
		String szgenesetlabel = "";
		String sztfsetlabel = "";
		/** non zero recs */
		int nprime;
		/** number of children */
		public int numchildren;
		/** can change state's values */
		boolean bchange;
		public int ndepth;
		int nparentcolorindex;
		boolean[] bvector;

		/** The logistic regression classifier */
		TSMiner_FastLogistic2 tranC;
		double[] recweight;
		double[] dTotalScore;
		double[][] dpvals;
		int[] orderA;
		/**
		 * Stores the number of genes regulated by a TF on a given path with the
		 * given interaction value. Dimensions: <br>
		 * [num TFs][num paths out of split][num distinct binding interaction
		 * values]
		 */
		int[][][] ncountvals;
		/**
		 * Stores the number of genes regulated by a TF on the paths other than
		 * the given path with the given interaction value. Dimensions: <br>
		 * [num TFs][num paths out of split][num distinct binding interaction
		 * values]
		 */
		int[][][] nothersum;
		double[][] dscoreother;
		double[][] dscore;
		double[][] ddiff;
		/**
		 * A TreeSet used to display important TFs based on activity score,
		 * which isn't a measure of significance in the statistical sense.
		 */
		TreeSet<SigTFRecv2>[] tsSigTFActivity;
		TreeSet<SigTFRecv2> tsSigTF;
		TreeSet<SigTFRecv2>[] tsSigTFEdgeSplit;
		TreeSet<SigTFRecv2>[] tsSigTFFull;
		Hashtable<Double, Integer>[] htScore;
		HashMap<String,List<Integer>> regulatedGeneindex;
		/**
		 * Stores the number of genes regulated by a TF on all paths out of a
		 * split with the given interaction value. Dimensions: <br>
		 * [num TFs][num distinct binding interaction values]
		 */
		int[][] ncountTotals;
		Integer[][] dIDEdgeSplit;
		Integer[][] dscoreEdgeSplit;
		double[][] dpvalEdgeSplit;
		double[][] ddiffEdgeSplit;
		double[][] ddiffEdgeFull;
		double[][] dexpectEdgeFull;
		double[][] dexpectEdgeSplit;
		double[][] dratioSplit;
		
		double[][] dgoEdgeFull;
		double[][] gopvalEdgeSplit;
		double[][] godiffEdgeSplit;
		double[][] goexpectEdgeSplit;
		double[][] goratioSplit;


		/**
		 * An array of activity scores for all TFs based on how the genes the TF
		 * regulates transition to the next states compared to how all genes
		 * into the split transistion.
		 */
		double[] dTFActivityScore;

		boolean[] bInNode;
		/** Number of genes going through this state??? */
		int numPath;
		/** probability of transitioning to each of the next states */
		double[] ptrans;
		/** expected number of transitions given the training sequences */
		double[] aprobs;
		/** pointer to each of the next states */
		public Treenode[] nextptr;
		/** pointer back to the parent */
		public Treenode parent;

		PText thepredictText;
		PText goText;
		PText genesetText;
		PText tfsetText;
		/**
		 * forward variable for current variables <br>
		 * p(x_1,...,x_i,pi_i=k)
		 */
		double df;
		/**
		 * backward variable for current variables <br>
		 * p(x_(i+1),...,x_L|pi_i=k)
		 */
		double db;
		boolean binit;
		Double[][] splitES;
		public String[][] resultTF;
		public String[][] deTF;

		/**
		 * Calls inittreenode with a null parent
		 */
		Treenode() {
			inittreenode(null);
		}

		Treenode(Treenode parent) {
			inittreenode(parent);
		}

		/**
		 * Calls inittreenode with parent
		 */
		void inittreenode(Treenode parent) {
			niteration = 0;
			ninstanceiteration = 0;
			nrec = -1;
			dmean = 0;
			dsigma = .5;
			this.parent = parent;
			dEsum = 0;
			dEsumsq = 0;
			dPsum = 0;
			nprime = 0;
			bchange = true;
			binit = false;
			
			if (parent == null) {
				ndepth = 0;
			} else {
				ndepth = parent.ndepth + 1;
			}
			nextptr = new Treenode[nmaxchild];

			if (BREGTSMiner) {
				ptrans = new double[nmaxchild];
				aprobs = new double[nmaxchild];
			}

			bvector = new boolean[numbits];

			if (bindingData.gene2RegBinding != null) {
				dTotalScore = new double[numbits];
			}

			for (int nindex = 0; nindex < numbits; nindex++) {
				bvector[nindex] = false;
			}

			numchildren = 0;
			for (int nindex = 0; nindex < nextptr.length; nindex++) {
				nextptr[nindex] = null;
			}
		}
		/**
		 * For making copy of nodes
		 */
		public Object clone() {
			if (parent == null) {
				if (BDEBUG) {
					System.out.println("in clone........");
				}
				htBackPts = new Hashtable<Treenode, BackPtrRec>();
			}

			Treenode tnode = new Treenode();
			tnode.nID = nID;
			tnode.dmean = dmean;
			tnode.dsigma = dsigma;
			tnode.ndepth = ndepth;
			tnode.dEsum = dEsum;
			tnode.dEsumsq = dEsumsq;
			tnode.dPsum = dPsum;
			tnode.nprime = nprime;
			tnode.numchildren = numchildren;
			tnode.thepredictText = thepredictText;
			tnode.parent = null;
			tnode.bchange = bchange;
			tnode.df = df;
			tnode.db = db;
			tnode.binit = binit;
			tnode.nextptr = new Treenode[nmaxchild];

			if (BREGTSMiner) {
				tnode.ptrans = new double[nmaxchild];
				tnode.aprobs = new double[nmaxchild];
			} else if (tranC != null) {
				tnode.tranC = (TSMiner_FastLogistic2) tranC.clone(); // do we need
				// a clone
				// yes
			}
			
			if(recweight != null){
				tnode.recweight = new double[recweight.length];
			}

			for (int nindex = 0; nindex < nextptr.length; nindex++) {
				if (nextptr[nindex] == null) {
					tnode.nextptr[nindex] = null;
				} else {
					BackPtrRec childNodeRec = htBackPts.get(nextptr[nindex]);

					if (childNodeRec != null) {
						// we've already cloned the child node, just have it
						// link to it
						tnode.nextptr[nindex] = childNodeRec.childNode;
						// have the child link back to this node
						if (BDEBUG) {
							System.out.println("[[" + childNodeRec + " "
									+ childNodeRec.childNode + " "
									+ childNodeRec.childNode.parentptrA);
						}

						int nfindparentindex = 0;
						while (nextptr[nindex].parentptrA[nfindparentindex] != this) {
							nfindparentindex++;
						}
						if (BDEBUG) {
							System.out.println("nfindparentindex = "
									+ nfindparentindex);
						}
						childNodeRec.childNode.parentptrA[nfindparentindex] = tnode;
					} else {
						// cloning a new node
						tnode.nextptr[nindex] = (Treenode) nextptr[nindex].clone();
						tnode.nextptr[nindex].parent = tnode;

						if (nextptr[nindex].parentptrA != null) {
							tnode.nextptr[nindex].parentptrA = new Treenode[nextptr[nindex].parentptrA.length];
							for (int nparentindex = 0; nparentindex < nextptr[nindex].parentptrA.length; nparentindex++) {
								if (nextptr[nindex].parentptrA[nparentindex] == this) {
									if (BDEBUG) {
										System.out.println("~~``"
												+ tnode.nextptr[nindex].dmean
												+ " " + nparentindex);
									}
									tnode.nextptr[nindex].parentptrA[nparentindex] = tnode;
								}
							}

							if (nextptr[nindex].parentptrA.length >= 2) {
								if (BDEBUG) {
									System.out.println(dmean + " "
											+ nextptr[nindex].dmean);
								}
								htBackPts.put(nextptr[nindex], new BackPtrRec(tnode.nextptr[nindex]));
							}
						}
					}
				}

				if (BREGTSMiner) {
					tnode.ptrans[nindex] = ptrans[nindex];
					tnode.aprobs[nindex] = aprobs[nindex];
				}
			}

			tnode.bvector = new boolean[numbits];
			for (int nindex = 0; nindex < numbits; nindex++) {
				tnode.bvector[nindex] = bvector[nindex];
			}

			return tnode;
		}

		/**
		 * 
		 * @return a map from TF names (excluding the " [#]" used to indicate
		 *         the primary path out of the split for the TF) to the best
		 *         activity score seen at any of this node's ancestors.
		 */
		public HashMap<String, Double> getAncestorActivityScores() {
			HashMap<String, Double> bestSeenScores = new HashMap<String, Double>();
			// TODO need to handle the case where there are multiple parents?
			if (parent != null) {
				return parent.getAncestorActivityScores(bestSeenScores);
			} else {
				return bestSeenScores;
			}
		}

		/**
		 * Helper function for the other version of getAncestorActivityScores.
		 * 
		 * @param bestSeenScores
		 *            The best activity scores seen in ancestors of the original
		 *            node before the search reached the current node.
		 * @return a map from TF names (excluding the " [#]" used to indicate
		 *         the primary path out of the split for the TF) to the best
		 *         activity score seen at any of this nodes ancestors or in the
		 *         bestSeenScores.
		 */
		private HashMap<String, Double> getAncestorActivityScores(
				HashMap<String, Double> bestSeenScores) {
			// Enumerate all TFs in the tree set and add their activty
			// score to the map if there was no score for that TF or if
			// the score at the current node is better (i.e. lower, because
			// the activity scores have been transformed to resemble p-values)
			for(int i = 0; i < numchildren; i++)
			{
				Iterator<TSMiner_Timeiohmm.SigTFRecv2> itr = tsSigTFActivity[i].iterator();
				TSMiner_Timeiohmm.SigTFRecv2 theTFRec;
				while (itr.hasNext()) {
					theTFRec = itr.next();
	
					String tfName = theTFRec.szname;
					// Check if the path out of the split has been added as " [#]"
					// to the TF name, and if so remove it
					// need to be done
//					if (tfName.endsWith("]")) {
//						tfName = tfName.substring(0, tfName.lastIndexOf(" ["));
//					}
	
					// Add the TF and score to the map if they are not already
					// present
					if (!bestSeenScores.containsKey(tfName)) {
						// The "pval" is actually the transformed activity score
						bestSeenScores.put(tfName, Double.valueOf(theTFRec.dpval));
					}
					// Otherwise add them if the current score is lower
					// than the previously seen best lower score
					else {
						double bestScore = ((Double) bestSeenScores.get(tfName))
								.doubleValue();
	
						if (theTFRec.dpval < bestScore) {
							bestSeenScores.put(tfName, Double
									.valueOf(theTFRec.dpval));
						}
					}
				} // end iteration over all TF activity scores at this node
			}
			if (parent != null) {
				bestSeenScores = parent
						.getAncestorActivityScores(bestSeenScores);
			}

			return bestSeenScores;
		}
	}

	/**
	 * A record with information about a state a back pointer points to
	 */
	static class BackPtrRec {
		BackPtrRec(Treenode childNode) {
			this.childNode = childNode;
		}
		Treenode childNode;
	}

	// ///////////////////////////////////////////////////////////////////
	/**
	 * Returns true iff treeptr or a descendent has two or more parents
	 */
	public boolean hasMerge(Treenode treeptr) {
		boolean bhasmerge = ((treeptr.parentptrA != null) && (treeptr.parentptrA.length >= 2));
		for (int nchild = 0; (nchild < treeptr.numchildren) && (!bhasmerge); nchild++) {
			bhasmerge = hasMerge(treeptr.nextptr[nchild]);
		}

		return bhasmerge;
	}

	// ////////////////////////////////////////////////////////////////////
	/**
	 * Returns a string with the model parameters
	 */
	public String saveString(Treenode treecopy) {
		StringBuffer szbuf = new StringBuffer("");
		;

		bhasmerge = hasMerge(treecopy);

		if (bhasmerge) {
			szbuf.append("MERGED TREE\n");
		}

		if (BREGTSMiner) {
			szbuf.append("REGTSMiner Num. Coefficients\t"
					+ bindingData.regNames.length + "\n");
		} else {
			szbuf.append("Num. Coefficients\t" + bindingData.regNames.length
					+ "\n");
		}

		HashMap<Treenode, Integer> ptrMap = new HashMap<Treenode, Integer>();
		ntotalID = -1;
		szbuf.append(saveStringHelp(treecopy, ptrMap));
		return szbuf.toString();
	}

	/**
	 * Helper function for generating a string with model parameters that
	 * traverses all the states
	 */
	private StringBuffer saveStringHelp(Treenode ptr,
			HashMap<Treenode, Integer> ptrMap) {
		StringBuffer szbuf = new StringBuffer();
		if (ptr != null) {
			if (bhasmerge) {
				Integer obj = ptrMap.get(ptr);
				int nID;
				if (obj != null) {
					nID = obj.intValue();
				} else {
					ntotalID++;
					nID = ntotalID;
					ptrMap.put(ptr, new Integer(nID));
				}
				szbuf.append(nID + "\t");
			}

			szbuf.append(ptr.dmean + "\t" + ptr.dsigma + "\t" + ptr.numchildren
					+ "\n");
			if (ptr.numchildren > 1) {
				if (BREGTSMiner) {
					for (int nchild = 0; nchild < ptr.numchildren; nchild++) {
						szbuf.append(ptr.ptrans[nchild] + "\n");
					}
				} else {
					szbuf.append(ptr.tranC.saveClassifier(bindingData.regNames));
				}
			}

			for (int nindex = 0; nindex < ptr.numchildren; nindex++) {
				if (ptr.nextptr[nindex] != null) {
					szbuf.append(saveStringHelp(ptr.nextptr[nindex], ptrMap));
				}
			}
		}

		return szbuf;
	}

	/**
	 * Save the activity scores at this node and its children
	 * 
	 * @return
	 */
	public StringBuffer saveActivityScores(Treenode node) {
		StringBuffer buf = new StringBuffer();

		if (node != null && node.numchildren >= 1) {
			if (node.numchildren >= 2) {
				for (int t = 0; t < bindingData.regNames.length; t++) {
					buf.append(bindingData.regNames).append("\t");
					buf.append(node.dTFActivityScore[t]).append("\t");
					buf.append(node.ndepth).append("\n");
				}
			}

			for (int nindex = 0; nindex < node.numchildren; nindex++) {
				if (node.nextptr[nindex] != null) {
					buf.append(saveActivityScores(node.nextptr[nindex]));
				}
			}
		}

		return buf;
	}

	// ////////////////////////////////////////////////////////////////////
	/**
	 * Record of information about predictions
	 */
	static class PredictRec {
		double[] mean;
		double[] meansq;
		double[] var;
		double[] entropy;
	}

	/**
	 * Used to set probability of each state a gene with a TF-binding signature
	 * of binding will be in Also computes a record of other statistics based on
	 * these predictions
	 */
	public PredictRec predictTime(int[] binding, double dclassprob,
			TSMiner_Timeiohmm.Treenode treeptr) throws Exception {
		PredictRec thePredictRec = new PredictRec();
		thePredictRec.mean = new double[traindata[0].length];
		thePredictRec.var = new double[traindata[0].length];
		thePredictRec.meansq = new double[traindata[0].length];
		thePredictRec.entropy = new double[traindata[0].length];
		predictTimeHelp(treeptr, binding, dclassprob, 0, thePredictRec.mean,
				thePredictRec.var, thePredictRec.meansq, thePredictRec.entropy);

		return thePredictRec;

	}

	/**
	 * A recursive helper function to predictTime
	 */
	private void predictTimeHelp(Treenode node, int[] theInstance,
			double dweight, int nstep, double[] meanpredict,
			double[] varpredict, double[] meansqpredict, double[] entropypredict)
			throws Exception {
		if (nstep < meanpredict.length) {
			NumberFormat nf1 = NumberFormat.getInstance(Locale.ENGLISH);

			if (dweight < 1) {
				nf1.setMaximumIntegerDigits(0);
				nf1.setMinimumFractionDigits(2);
				nf1.setMaximumFractionDigits(2);
			} else {
				nf1.setMinimumFractionDigits(1);
				nf1.setMaximumFractionDigits(1);
			}
			node.dpredictweight = dweight;

			if (node.thepredictText != null) {
				node.thepredictText.setText(nf1.format(dweight));
			}

			meanpredict[nstep] += dweight * node.dmean;
			varpredict[nstep] += Math.pow(dweight * node.dsigma, 2);
			meansqpredict[nstep] += dweight * node.dmean * node.dmean;
			entropypredict[nstep] -= dweight * Math.log(dweight) / Math.log(2);

			if (node.binit) {
				double[] pdist;
				if (BREGTSMiner) {
					pdist = node.ptrans;
				} else {
					pdist = node.tranC.distributionForInstance(theInstance);
				}

				if (BDEBUG) {
					if (pdist.length == 1)
						System.out.println("warning " + pdist[0]);
				}

				for (int nchild = 0; nchild < node.numchildren; nchild++) {
					predictTimeHelp(node.nextptr[nchild], theInstance, dweight
							* pdist[nchild], nstep + 1, meanpredict,
							varpredict, meansqpredict, entropypredict);
				}
			} else {
				predictTimeHelp(node.nextptr[0], theInstance, dweight,
						nstep + 1, meanpredict, varpredict, meansqpredict,
						entropypredict);
			}
		}
	}

	// /////////////////////////////////////////////////////////////////////
	/**
	 * Used for debugging to see the contents of the model
	 */
	public void traverse(Treenode root, int ndepth, boolean bprintC) {
		if (BDEBUG) {
			for (int nindex = 0; nindex < ndepth; nindex++) {
				System.out.print("\t");
			}

			System.out.print(root.numPath + "\t" + root.dmean + "\t"
					+ root.dsigma + "\t" + root.numchildren + "\t" + root.binit
					+ "\t" + root.dpredictweight + "\t" + root.ndepth + "\t#"
					+ root.numdownleaves + "\t");

			if (root.parent != null)
				System.out.print("p " + root.parent.dmean);

			if (root.parentptrA != null) {
				for (int i = 0; i < root.parentptrA.length; i++)
					System.out.print("\t"
							+ nf2.format(root.parentptrA[i].dmean));
			}

			if (root.ptrans != null)
				System.out.print("\t" + root.ptrans[0] + "\t" + root.ptrans[1]);

			for (int nindex = ndepth; nindex < theDataSet.data[0].length - 1; nindex++) {
				System.out.print("\t");
			}

			if (root.nextptr[0] != null) {
				System.out.println();

				if ((bprintC) && (root.numchildren > 1)) {
					System.out.println(root.tranC);
				}
			} else {
				System.out.println();
			}

			for (int nindex = 0; nindex < root.numchildren; nindex++) {
				if (root.nextptr[nindex] != null) {
					traverse(root.nextptr[nindex], ndepth + 1, bprintC);
				}
			}
		}
	}

	// ///////////////////////////////////////////////////////////////////////

	/**
	 * Combines the results of the forward and backward algorithm to obtain the
	 * probability of each gene going through each state when using viterbi
	 * training
	 */

	public void instanceAEV(Treenode root, double[] vals, 
			double[] bpvals, int[] bestpath, double dpj, int nrec)
			throws Exception {
		Treenode tempptr = root;

		for (int ndepth = 0; ndepth < bestpath.length; ndepth++) {
			tempptr.dEsum += vals[ndepth];
			tempptr.dPsum += 1;    
			tempptr.dEsumsq += Math.pow(vals[ndepth], 2);
			tempptr.nprime++;

			if (BREGTSMiner) {
				tempptr.aprobs[bestpath[ndepth]]++;
			}

			if (tempptr.numchildren >= 1) {
				tempptr.recweight[tempptr.numchildren * nrec + bestpath[ndepth]] = 1;
			}

			tempptr = tempptr.nextptr[bestpath[ndepth]];
		}
	}
	// /////////////////////////////////////////////////////////////////////////////////
	/**
	 * Initializes all the fields that will be used by instanceAE
	 */
	public void initAE(Treenode root) {
		if (root != null) {
			for (int nchild = 0; nchild < root.numchildren; nchild++) {
				if (BREGTSMiner) {
					root.aprobs[nchild] = 0;
				}
				initAE(root.nextptr[nchild]);
			}
			
			if ((root.recweight == null) || (root.recweight.length != root.numchildren * traindata.length)) {
				root.recweight = new double[root.numchildren * traindata.length];
			}
			
			root.dEsum = 0;
			root.dEsumsq = 0;
			root.dPsum = 0;
			root.nprime = 0; 
			root.niteration = 0;
			root.ninstancerec = 0;
			root.ninstanceiteration = 0;
			root.nrec = -1;
		}
	}

	/**
	 * Combines the results of the forward and backward algorithm to obtain the
	 * probability of each gene going through each state when using baum-welch
	 * training
	 */
	public void instanceAE(Treenode root, int ndepth, double[] vals, 
			double dpj, int nrec) throws Exception {
		if ((root.ninstancerec != nrec) || (root.ninstanceiteration != nglobaliteration)) {
			root.ninstanceiteration = nglobaliteration;
			root.ninstancerec = nrec; //The index of gene
			
			if (root.nextptr[0] != null) {
				for (int nchild = 0; nchild < root.numchildren; nchild++) {
					double dmultval = root.nextptr[nchild].df * root.nextptr[nchild].db / dpj;
					root.recweight[root.numchildren * nrec + nchild] = dmultval;

					if (BREGTSMiner) {
						root.aprobs[nchild] += dmultval;
					}		
					instanceAE(root.nextptr[nchild], ndepth + 1, vals, dpj, nrec);
				}
			}

			double dweight = root.df * root.db / dpj;
			//System.out.println(root.ndepth+" "+root.dmeans[1]+" "+dweight);
				if ((dweight > 0) && (!Double.isInfinite(dweight))) {
					root.dEsum += vals[ndepth] * dweight;
					root.dPsum += dweight;
					root.dEsumsq += vals[ndepth] * vals[ndepth] * dweight;
					root.nprime++;
				}	
			}
		}	

	// ///////////////////////////////////////////////////////////////////////
	/**
	 * Initializes all the fields that will be used by instanceAEV
	 */
	public void initAEV(Treenode root) {
		if (root != null) {
			if (BREGTSMiner) {
				for (int nchild = 0; nchild < root.aprobs.length; nchild++) {
					root.aprobs[nchild] = 0;
				}
			}

			for (int nchild = 0; nchild < root.numchildren; nchild++) {
				if (BREGTSMiner) {
					root.aprobs[nchild] = 0;
				}
				initAEV(root.nextptr[nchild]);
			}

			if ((root.recweight == null) || (root.recweight.length != root.numchildren * traindata.length)) {
				root.recweight = new double[root.numchildren * traindata.length];
			} else {
				for (int i = 0; i < root.recweight.length; i++) {
					root.recweight[i] = 0;
				}
			}

			root.dEsum = 0;
			root.dEsumsq = 0;
			root.dPsum = 0;
			root.nprime = 0;
		}
	}

	// //////////////////////////////////////////////////////////////////////
	/**
	 * Forces children of a parent to both have the average of their standard
	 * deviations
	 */
	public void averageChildrenSigmas(Treenode root) {
		if (root.numchildren >= 2) {
			double dsigmasum = 0;
			double dpsum = 0;
			for (int nchild = 0; nchild < root.numchildren; nchild++) {
				dsigmasum += root.nextptr[nchild].dPsum
						* root.nextptr[nchild].dsigma;
				dpsum += root.nextptr[nchild].dPsum;
				if (BDEBUG) {
					System.out.println("##" + root.nextptr[nchild].dsigma
							+ "\t" + root.nextptr[nchild].dPsum);
				}
			}
			double dsigmaavg = dsigmasum / dpsum;
			if (BDEBUG) {
				System.out.println("&&&" + dpsum + "\t" + root.dPsum + "\t"
						+ dsigmasum + "\t" + dsigmaavg);
			}

			for (int nchild = 0; nchild < root.numchildren; nchild++) {
				root.nextptr[nchild].dsigma = dsigmaavg;
			}
		}

		for (int nchild = 0; nchild < root.numchildren; nchild++) {
			averageChildrenSigmas(root.nextptr[nchild]);
		}
	}

	// ///////////////////////////////////////////////////////////////////////

	/**
	 * Updates the mean, standard deviation, and classifier parameters based on
	 * the current expectation of which genes will go through which paths
	 */
	public void updateParams(Treenode root, int nchildnum, double[][] data,
			int ndepth) throws Exception {
		if (!root.bchange) {
		    if (root.nextptr[0] != null) {
				for (int nchild = 0; nchild < root.numchildren; nchild++) {
					updateParams(root.nextptr[nchild], nchild, data, ndepth + 1);
				}
			}
		} else {
				if (root.dPsum > 0) {
					root.dmean = root.dEsum / root.dPsum;	
				} else {
					root.dmean = 0;
				}
				
				if (root.nprime <= 1) {//root.nprime为全部基因，一定大于1.
					root.dsigma = Math.max(dminstddev, DEFAULTSIGMA);
				} else {
					double dval = (double) root.nprime
							/ (double) (root.nprime - 1)
							* (root.dEsumsq/ root.dPsum - Math.pow(root.dEsum/ root.dPsum, 2));
					if (dval <= 0) {
						root.dsigma = Math.max(dminstddev, DEFAULTSIGMA);
					} else {
						root.dsigma = Math.max(dminstddev, Math.sqrt(dval));
					}
				}
			
			if (root.nextptr[0] != null) {
				if (root.numchildren > 1) {
					if (BREGTSMiner) { //TRUE
						root.binit = true;
					} else {
						if (root.tranC == null) {
							root.tranC = new TSMiner_FastLogistic2(
									trainBinding.theInstancesIndex[root.ndepth][root.numchildren - 2],
									trainBinding.theInstances[root.ndepth][root.numchildren - 2],
									trainBinding.theInstancesRegIndex[root.ndepth][root.numchildren - 2],
									trainBinding.theInstancesReg[root.ndepth][root.numchildren - 2],
									trainBinding.ylabels[root.ndepth][root.numchildren - 2],
									root.recweight, root.numchildren, numbits,
									trainBinding.regTypes);
							root.tranC.setRidge(RIDGE);
						} else {
							root.tranC.reinit(
											trainBinding.theInstancesIndex[root.ndepth][root.numchildren - 2],
											trainBinding.theInstances[root.ndepth][root.numchildren - 2],
											trainBinding.theInstancesRegIndex[root.ndepth][root.numchildren - 2],
											trainBinding.theInstancesReg[root.ndepth][root.numchildren - 2],
											trainBinding.ylabels[root.ndepth][root.numchildren - 2],
											root.recweight, root.numchildren,
											trainBinding.regTypes);
						}

						if (root.nextptr[0].dmean > root.nextptr[1].dmean) {
							root.tranC.childExpStatus = 1;
						} else {
							root.tranC.childExpStatus = -1;
						}
						root.tranC.train();
						root.binit = true;
					}
				}

				double dsum = 0;
				if (BREGTSMiner) {
					for (int nchild = 0; nchild < root.numchildren; nchild++) {
						if (root.nextptr[nchild].bchange) {
							dsum += root.aprobs[nchild];
						}
					}
				}

				for (int nchild = 0; nchild < root.numchildren; nchild++) {
					if (BREGTSMiner) {
						if (dsum > 0) {
							root.ptrans[nchild] = root.aprobs[nchild] / dsum;
							if (BDEBUG) {
								System.out.println("probs: "
										+ root.aprobs[nchild] + " " + dsum
										+ " " + root.ptrans[nchild]);
							}
						} else {
							root.ptrans[nchild] = 0;
						}
						if (Double.isNaN(root.ptrans[nchild]))
							System.out.println("%%%%%% " + root.ptrans[nchild]
									+ "\t" + root.aprobs[nchild] + "\t" + dsum);
					}
					updateParams(root.nextptr[nchild], nchild, data, ndepth + 1);		
				}
			}
		}
	}

	// ////////////////////////////////////////////////////////////////////////////
	/**
	 * Executes the forward algorithm portion of the baum-welch algorithm
	 */
	public void forwardalg(double[] vals, 
			int[][][] theInstanceIndex, double[][][] theInstance, int ni,
			double da, Treenode node, int nrec, Treenode parentptr,
			int instanceindex) throws Exception {
		if (node != null) {

				double dxmudiff = (vals[ni] - node.dmean);
				node.dens = Math.exp(-dxmudiff * dxmudiff / (2 * node.dsigma * node.dsigma))
						/(node.dsigma * StatUtil.TWOPISQRT);
			
			if (node.parentptrA == null) {
				node.df = node.parent.df * da * node.dens;		
			} else { //no run
				if ((node.nrec != nrec)|| (node.niteration != nglobaliteration)) {
					node.df = parentptr.df * da * node.dens;
					if (node.niteration != nglobaliteration)
						node.niteration = nglobaliteration;
					if (node.nrec != nrec)
						node.nrec = nrec;	
				} else {
					node.df += parentptr.df * da * node.dens;
					node.niteration++;
				}
			}

			double[] ptranlogit; 

			if (BREGTSMiner) {
				ptranlogit = node.ptrans;
			} else if ((node.numchildren > 1) && (node.binit)) {
				ptranlogit = node.tranC.distributionForInstance(
						theInstanceIndex[node.ndepth][instanceindex],
						theInstance[node.ndepth][instanceindex]);
			} else {
				ptranlogit = CONSTANTA[node.numchildren];
			}

			for (int nchild = 0; nchild < node.numchildren; nchild++) {
				try {
					forwardalg(vals, theInstanceIndex, theInstance,
							ni + 1, ptranlogit[nchild], node.nextptr[nchild],
							nrec, node, instanceindex);
				} catch (Exception ex) {
					ex.printStackTrace();
					System.out.println(node.tranC);
				}
			}
		}
	}

	// ///////////////////////////////////////////////////////////////////////////
	/**
	 * Executes the backward portion of the baum-welch algorithm
	 */
	public void backalg(double[] vals, int[][][] theInstanceIndex,
			double[][][] theInstance, int ni, Treenode node, int nrec,
			boolean bforward, int instanceindex) throws Exception {
		
		if (node.nextptr[0] == null) {
			node.db = 1;
		} else {
			node.db = 0;
			double[] ptranlogit;
			if (BREGTSMiner) {
				ptranlogit = node.ptrans;
			} else if ((node.numchildren > 1) && (node.binit)) {
				ptranlogit = node.tranC.distributionForInstance(
						theInstanceIndex[node.ndepth][instanceindex],
						theInstance[node.ndepth][instanceindex]);
			} else {
				ptranlogit = CONSTANTA[node.numchildren];
			}
			
			for (int nchild = 0; nchild < node.numchildren; nchild++) {
				Treenode nextnode = node.nextptr[nchild];
				backalg(vals, theInstanceIndex, theInstance, ni + 1, nextnode, nrec, bforward, instanceindex);

				double dens;

				if (bforward) {
					dens = nextnode.dens;
				} else {
					double dxmudiff = (vals[ni + 1] - nextnode.dmean);
					dens = Math.exp(-dxmudiff * dxmudiff
							/ (2 * nextnode.dsigma * nextnode.dsigma))
							/ (nextnode.dsigma * StatUtil.TWOPISQRT);
				}

				if (BDEBUG) {
					if (nchild >= ptranlogit.length) {
						System.out.println("**********" + nchild + "\t"
								+ ptranlogit.length);
					}
				}

				node.db += ptranlogit[nchild] * nextnode.db * dens;
				
				if (node.db < 0) {
					System.out.println("##" + node.db + "\t" + nextnode.dens);
					for (int i = 0; i < ptranlogit.length; i++)
					System.out.print(ptranlogit[i] + "\t");
					System.out.println();
				}
			}
		}
	}


	// //////////////////////////////////////////////////////////////
	/**
	 * Clears any prior assignments of genes to paths through the model
	 */
	public void clearCounts(Treenode ptr) {
		if (ptr != null) {
			ptr.bInNode = new boolean[theDataSet.data.length];
			ptr.numPath = 0;
			ptr.genelist.clear();
			if (bindingData.gene2RegBinding != null) {
				ptr.htScore = (Hashtable<Double, Integer>[]) new Hashtable[numbits];
				ptr.regulatedGeneindex = new HashMap<String,List<Integer>>();	
				ptr.dTotalScore = new double[numbits];
				for (int nbindingindex = 0; nbindingindex < ptr.htScore.length; nbindingindex++) {
					ptr.htScore[nbindingindex] = new Hashtable<Double, Integer>();
				}
			}
			for (int nchild = 0; nchild < ptr.numchildren; nchild++) {
				clearCounts(ptr.nextptr[nchild]);
			}
		}
	}

	/*
	public 	String[][]  viterbi2(Treenode treeptr) throws Exception {
		String[][] TreeModle_DREM = new String[theDataSet.data.length][theDataSet.data[0].length*2];
		int[] bestpath;
		treeptr.genelist.clear();
		bestpath = new int[theDataSet.data[0].length];
		storedbestpath = new BigInteger[theDataSet.data.length];
		clearCounts(treeptr);
		treeptr.numPath = theDataSet.data.length;
		if (BDEBUG) {
			System.out.print("!!!!!!" + theDataSet.data.length);
			if (bindingData.gene2RegBinding != null)
				System.out.println(" " + bindingData.gene2RegBinding[treeptr.ndepth].length);
		}

		for (int nrow = 0; nrow < theDataSet.data.length; nrow++) {
			if (bindingData.gene2RegBinding != null) {
				computevlogit2(theDataSet.data[nrow],
						theDataSet.pmavalues[nrow],
						bindingData.gene2RegBinding,
						bindingData.gene2RegBindingIndex, 0, treeptr, bestpath,
						nrow);  
			} else {
				computevlogit2(theDataSet.data[nrow],
						theDataSet.pmavalues[nrow], null, null, 0, treeptr,
						bestpath, 0); //输出最优bestpath
			}
			

			BigInteger nsum = new BigInteger("0");
			Treenode tempptr = treeptr;  //建立新节点
			tempptr.bInNode[nrow] = true;
			tempptr.genelist.add(nrow); //0节点的genelist加入所有基因

			for (int nindex = 0; nindex < bestpath.length - 1; nindex++) {
				int currdepth = tempptr.ndepth;
				tempptr = tempptr.nextptr[bestpath[nindex]];
				tempptr.bInNode[nrow] = true;
				tempptr.numPath++;//经过此节点的基因数+1
				tempptr.genelist.add(nrow); //存储经过此节点的基因编号

				if (bindingData.gene2RegBinding != null) {
					int ntfindex = 0;
					for (int nbit = 0; nbit < numbits; nbit++) {
						while ((ntfindex < bindingData.gene2RegBindingIndex[currdepth][nrow].length)
								&& (nbit > bindingData.gene2RegBindingIndex[currdepth][nrow][ntfindex])) {
							ntfindex++;
						}

						Double objMap;
						if ((ntfindex < bindingData.gene2RegBindingIndex[currdepth][nrow].length)
								&& (nbit == bindingData.gene2RegBindingIndex[currdepth][nrow][ntfindex])) {
							tempptr.dTotalScore[nbit] += bindingData.gene2RegBinding[currdepth][nrow][ntfindex];
							objMap = new Double(bindingData.gene2RegBindingSigned[currdepth][nrow][ntfindex]);
						} else {
							objMap = new Double(0.0);
						}
						Integer objCount = tempptr.htScore[nbit].get(objMap);
						if (objCount != null) {
							tempptr.htScore[nbit].put(objMap, new Integer(objCount.intValue() + 1));
						} else {
							tempptr.htScore[nbit].put(objMap, new Integer(1));
						}
					}
				}
				// szsum += ""+bestpath[nindex];
				// nsum += (int)
				// Math.pow(nmaxchild,bestpath.length-nindex-1)*bestpath[nindex];
				
				nsum = nsum.add((new BigInteger("" + nmaxchild)).pow(
						bestpath.length - nindex - 1).multiply(
						new BigInteger("" + bestpath[nindex])));
			}
			
			for(int i=0;i<bestpath.length;i++){
				TreeModle_DREM[nrow][i] = bestpath[i]+"";
				TreeModle_DREM[nrow][i+bestpath.length] = theDataSet.data[nrow][i]+"";
			}
			storedbestpath[nrow] = nsum;
		}
		
		return TreeModle_DREM;
	}
	
	public double computevlogit2(double[] vals, int[] pma,
			double[][][] theInstance, int[][][] theInstanceIndex, int ndepth,
			Treenode node, int[] bestpath, int instanceindex) throws Exception {
		double dmax = Double.NEGATIVE_INFINITY;
		int bestchild = 0;
		double dlogout = 0;
		double dval = 0;

		if (node.dsigma > 0) {
			if (pma[ndepth] != 0) {
				if (ndepth == 0) {
					dlogout = 0;
				} else {
					dlogout = Math.log(StatUtil.normaldensity(vals[ndepth],
							node.dmean, node.dsigma));
				}
			} else {
				dlogout = 0;
			}
		}

		if (node.numchildren == 0) {
			return dlogout;
		}

		int[] currbestpath = new int[bestpath.length];

		// log (a) + c(log(b_1) + log(b_2) + ... + log(b_n))
		// = log(a*(b_1*b_2*...*b_n)^c)

		double[] ptranlogit;

		if (BREGTSMiner) {
			ptranlogit = node.ptrans;
		} else if ((node.numchildren > 1) && (node.binit)) {
			ptranlogit = node.tranC.distributionForInstance(
					theInstanceIndex[node.ndepth][instanceindex],
					theInstance[node.ndepth][instanceindex]);
		} else {
			ptranlogit = CONSTANTA[node.numchildren];
		}

		for (int nchild = 0; nchild < node.numchildren; nchild++) {
			dval = dlogout;
			if (node.nextptr[nchild] != null) {
				double dvi = computevlogit2(vals, pma, theInstance,
						theInstanceIndex, ndepth + 1, node.nextptr[nchild],
						bestpath, instanceindex);
				if (ptranlogit[nchild] == 0) {
					dval += Math.log(MINPROB);
				} else {
					dval += Math.log(ptranlogit[nchild]) + dvi;
				}
			}

			if (dval > dmax) {
				dmax = dval;
				bestchild = nchild;
				for (int nindex = 0; nindex < bestpath.length; nindex++) {
					currbestpath[nindex] = bestpath[nindex];
				}
			}
		}

		for (int nindex = 0; nindex < bestpath.length; nindex++) {
			bestpath[nindex] = currbestpath[nindex];
		}
		bestpath[ndepth] = bestchild;

		return dmax;
	}
	*/
	//////////////////////////////////////////////////////////////////////////////
	/**
	 * Determines the most likely path of each gene through the model
	 */
	public 	String[]  viterbi(Treenode treeptr, double[][] data) throws Exception {
		String[] TreeModle_TSMiner = new String[theDataSet.data.length];
		int[] bestpath;
		treeptr.genelist.clear();
		bestpath = new int[data[0].length];
		storedbestpath = new BigInteger[data.length];
		clearCounts(treeptr);
		treeptr.numPath = data.length;
		if (BDEBUG) {
			if (bindingData.gene2RegBinding != null)
				System.out.println(" " + bindingData.gene2RegBinding[treeptr.ndepth].length);
		}
		
		for (int nrow = 0; nrow < data.length; nrow++) {
			if (bindingData.gene2RegBinding != null) {
				computevlogit(data[nrow],
						bindingData.gene2RegBinding,
						bindingData.gene2RegBindingIndex, 0, treeptr, bestpath,
						nrow);  
			} else {
				computevlogit(data[nrow],
						null, null, 0, treeptr,
						bestpath, 0); //输出最优bestpath
			}
			

			BigInteger nsum = new BigInteger("0");
			Treenode tempptr = treeptr;  //建立新节点
			tempptr.bInNode[nrow] = true;
			tempptr.genelist.add(nrow); //0节点的genelist加入所有基因

			for (int nindex = 0; nindex < bestpath.length - 1; nindex++) {
				int currdepth = tempptr.ndepth;
				tempptr = tempptr.nextptr[bestpath[nindex]];
				tempptr.bInNode[nrow] = true;
				tempptr.numPath++;//经过此节点的基因数+1
				tempptr.genelist.add(nrow); //存储经过此节点的基因编号

				if (bindingData.gene2RegBinding != null) {
					int ntfindex = 0;
					for (int nbit = 0; nbit < numbits; nbit++) {
						while ((ntfindex < bindingData.gene2RegBindingIndex[currdepth][nrow].length)
							&& (nbit > bindingData.gene2RegBindingIndex[currdepth][nrow][ntfindex])) {
							ntfindex++;
						}

						Double objMap;
						if ((ntfindex < bindingData.gene2RegBindingIndex[currdepth][nrow].length)
								&& (nbit == bindingData.gene2RegBindingIndex[currdepth][nrow][ntfindex])) {
							
								tempptr.dTotalScore[nbit] += bindingData.gene2RegBinding[currdepth][nrow][ntfindex];
								objMap = new Double(bindingData.gene2RegBindingSigned[currdepth][nrow][ntfindex]);
								String Upertf = bindingData.regNames[nbit].toUpperCase();
								List<Integer> genelist = tempptr.regulatedGeneindex.get(Upertf);
								if(genelist!=null){
									genelist.add(filter2rawIndex.get(nrow));
									tempptr.regulatedGeneindex.put(Upertf, genelist);
								}else{
									genelist = new ArrayList<Integer>();
									genelist.add(filter2rawIndex.get(nrow));
									tempptr.regulatedGeneindex.put(Upertf, genelist);
								}
						} else {
							objMap = new Double(0.0);
							//System.out.println(objMap);
						}
						if(tempptr.htScore[nbit].get(objMap) != null){
							int objCount = tempptr.htScore[nbit].get(objMap);
							tempptr.htScore[nbit].put(objMap, objCount+1);
						}else{
							tempptr.htScore[nbit].put(objMap, 1);
						}
					}
				}
				
				// szsum += ""+bestpath[nindex];
				// nsum += (int)
				// Math.pow(nmaxchild,bestpath.length-nindex-1)*bestpath[nindex];
				
				nsum = nsum.add((new BigInteger("" + nmaxchild)).pow(
						bestpath.length - nindex - 1).multiply(
						new BigInteger("" + bestpath[nindex])));
			}
			String bp = bestpath[0]+"";
			for(int i=1;i<bestpath.length-1;i++) bp = bp+","+String.valueOf(bestpath[i]);
			TreeModle_TSMiner[nrow] = bp;
			storedbestpath[nrow] = nsum;
		}
		
		return TreeModle_TSMiner;
	}

	//////////////////////////////////////////////////////////////////////////////////
	/**
	 * Computes association scores for the transcription factors overall on
	 * paths and on paths out of splits conditional on the set of genes going
	 * into the split
	 */
	public void computeStats(Treenode treeptr, Treenode rootptr) {
		if (treeptr != null) {
			
			treeptr.dscore = new double[numbits][treeptr.numchildren];
			treeptr.ddiff = new double[numbits][treeptr.numchildren];
			treeptr.dpvals = new double[numbits][treeptr.numchildren];
			treeptr.ncountvals = new int[numbits][treeptr.numchildren][bindingData.signedBindingValuesSorted.length];

			if (treeptr.numchildren >= 3) {
				treeptr.nothersum = new int[numbits][treeptr.numchildren][bindingData.signedBindingValuesSorted.length];
				treeptr.dscoreother = new double[numbits][treeptr.numchildren];
			}

			int[][] ncountgrid = new int[2][bindingData.signedBindingValuesSorted.length];
			treeptr.ncountTotals = new int[numbits][bindingData.signedBindingValuesSorted.length];

			treeptr.dTFActivityScore = new double[numbits];

			treeptr.tsSigTFActivity = (TreeSet<SigTFRecv2>[]) new TreeSet[treeptr.numchildren];
			for(int i = 0; i < treeptr.numchildren; i++)
			treeptr.tsSigTFActivity[i] = new TreeSet<SigTFRecv2>(new SigTFRecv2Compare());
			
			treeptr.tsSigTF = new TreeSet<SigTFRecv2>(new SigTFRecv2Compare());
			for (int ntf = 0; ntf < numbits; ntf++) {
				double dnumpaths = 0;
				double dscoretotal = 0;
				for (int nindex = 0; nindex < treeptr.ncountTotals[ntf].length; nindex++) {
					treeptr.ncountTotals[ntf][nindex] = 0;
				}

				for (int nchild = 0; nchild < treeptr.numchildren; nchild++) {
					if (treeptr.nextptr[nchild].htScore == null) {
						for (int nvalindex = 0; nvalindex < treeptr.ncountvals[ntf][nchild].length; nvalindex++) {
							treeptr.ncountvals[ntf][nchild][nvalindex] = 0;
						}
					} else {
						for (int nel = 0; nel < bindingData.signedBindingValuesSorted.length; nel++) {
							Double key = new Double((double) bindingData.signedBindingValuesSorted[nel]);
							Integer objcount = treeptr.nextptr[nchild].htScore[ntf].get(key);
							
							if (objcount != null) {
								int ntempval = ((Integer) objcount).intValue();
								treeptr.ncountvals[ntf][nchild][nel] = ntempval;
								treeptr.ncountTotals[ntf][nel] += ntempval;
							} else {
								treeptr.ncountvals[ntf][nchild][nel] = 0;
							}
						}
					}

					treeptr.dscore[ntf][nchild] = treeptr.nextptr[nchild].dTotalScore[ntf]
							/ treeptr.nextptr[nchild].numPath;
					dnumpaths += treeptr.nextptr[nchild].numPath;
					dscoretotal += treeptr.nextptr[nchild].dTotalScore[ntf];
				}

				if (treeptr.numchildren == 2) {
					treeptr.ddiff[ntf][0] = treeptr.dscore[ntf][0]
							- treeptr.dscore[ntf][1];
					treeptr.dpvals[ntf][0] = computepval(
							treeptr.ncountvals[ntf], Math
									.abs(treeptr.ddiff[ntf][0]));
					treeptr.tsSigTF.add(new SigTFRecv2(
							bindingData.regNames[ntf], treeptr.dpvals[ntf][0]));

				} else if (treeptr.numchildren >= 3) {
					for (int nchild = 0; nchild < treeptr.numchildren; nchild++) {
						// divide by 0 need to check about that
						double ddiff = (dnumpaths - treeptr.nextptr[nchild].numPath);
						double davgother;

						if (ddiff > 0) {
							davgother = (dscoretotal - treeptr.nextptr[nchild].dTotalScore[ntf])
									/ ddiff;
						} else {
							davgother = 0;
						}
						treeptr.dscoreother[ntf][nchild] = davgother;

						treeptr.ddiff[ntf][nchild] = treeptr.dscore[ntf][nchild]
								- davgother;

						for (int nindex = 0; nindex < treeptr.ncountTotals[ntf].length; nindex++) {
							ncountgrid[0][nindex] = treeptr.ncountvals[ntf][nchild][nindex];
							ncountgrid[1][nindex] = treeptr.ncountTotals[ntf][nindex]
									- ncountgrid[0][nindex];
							treeptr.nothersum[ntf][nchild][nindex] = ncountgrid[1][nindex];
						}

						treeptr.dpvals[ntf][nchild] = computepval(ncountgrid,
								Math.abs(treeptr.ddiff[ntf][nchild]));
					}

					for (int nchild = 0; nchild < treeptr.numchildren; nchild++) {
						treeptr.tsSigTF.add(new SigTFRecv2(
								bindingData.regNames[ntf],
								treeptr.dpvals[ntf][treeptr.orderA[nchild]]));
					}
				}
				// Activity scores are only calculated if there is a split at
				// this node
				//System.out.println(bindingData.regPriors);
				if (treeptr.numchildren > 1 && bindingData.regPriors != null) {
					
					/**
					 * The number of genes on a particular path out of a split
					 * regulated by a TF. Regulated means the interaction value
					 * is not 0. Dimension:<br>
					 * [num paths out of split]
					 */
					int[] nRegThisPath = new int[treeptr.numchildren];
					/**
					 * The number of genes on all paths except the given path
					 * out of a split that are regulated by a TF. Regulated
					 * means the interaction value is not 0. Dimension:<br>
					 * [num paths out of split]
					 */
					int[] nRegOtherPaths = new int[treeptr.numchildren];
					/**
					 * The number of genes on all paths out of a split that are
					 * regulated by a TF. Regulated means the interaction value
					 * is not 0.
					 */
					int nRegAllPaths = 0;

					// Initialize nRegThisPath
					for (int nchild = 0; nchild < treeptr.numchildren; nchild++) {
						nRegThisPath[nchild] = 0;

					}

					// Populate nRegThisPath and nRegAllPaths
					for (int nchild = 0; nchild < treeptr.numchildren; nchild++) {
						for (int nel = 0; nel < bindingData.signedBindingValuesSorted.length; nel++) {
							// Only count genes the TF regulates, which are
							// those genes with an interaction value != 0
							if (bindingData.signedBindingValuesSorted[nel] != 0) {
								nRegThisPath[nchild] += treeptr.ncountvals[ntf][nchild][nel];
								nRegAllPaths += treeptr.ncountvals[ntf][nchild][nel];
							}
						}
					}

					// Populate nRegOtherPaths
					for (int nchild = 0; nchild < treeptr.numchildren; nchild++) {
						nRegOtherPaths[nchild] = nRegAllPaths
								- nRegThisPath[nchild];
					}

					// Now determine which path is the primary path out of the
					// split.
					// The primary path is the path with the most regulated
					// genes. In the case of the tie, the path with the fewest
					// genes
					// overall out of the split is the primary path.
					int nMaxPathVal = -1;
					for (int nchild = 0; nchild < treeptr.numchildren; nchild++) {
						nMaxPathVal = Math.max(nMaxPathVal,
								nRegThisPath[nchild]);
					}

					// See how many paths have this max value
					boolean[] bMaxPaths = new boolean[treeptr.numchildren];
					int nMaxPathCount = 0;
					for (int nchild = 0; nchild < treeptr.numchildren; nchild++) {
						if (nRegThisPath[nchild] == nMaxPathVal) {
							bMaxPaths[nchild] = true;
							nMaxPathCount++;
						} else {
							bMaxPaths[nchild] = false;
						}
					}

					// The index of the primay path
					int nPrimPathInd = -1;
					if (nMaxPathCount == 1) {
						for (int nchild = 0; nchild < treeptr.numchildren; nchild++) {
							if (bMaxPaths[nchild]) {
								nPrimPathInd = nchild;
							}
						}
					} else {
						// There was a tie so look at all paths out of the split
						// there were involved in the tie
						// to find the path with the fewest genes
						// If there is another tie, it won't affect the score
						// calculation so pick the first one found
						int nMinPathVal = Integer.MAX_VALUE;

						for (int nchild = 0; nchild < treeptr.numchildren; nchild++) {
							if (bMaxPaths[nchild]
									&& treeptr.nextptr[nchild].numPath < nMinPathVal) {
								nMinPathVal = treeptr.nextptr[nchild].numPath;
								nPrimPathInd = nchild;
							}
						}
					}

					// Calculate the number of genes on the paths that are not
					// the primary path
					int nNumOtherPaths = 0;
					for (int nchild = 0; nchild < treeptr.numchildren; nchild++) {
						if (nPrimPathInd != nchild) {
							nNumOtherPaths += treeptr.nextptr[nchild].numPath;
						}
					}

					// The formula for a TF's activity score is:
					// score = numer / denom
					// where numer = (prob binding is functional)^(num genes TF
					// regulates on primary path)
					// * (prob binding is not functional)^(num genes TF
					// regulates on other paths)
					// * (TF activity prior)
					// and denom = (ratio of all genes on TFs primary path)^(num
					// genes TF regulates on primary path)
					// * (ratio of all genes on other paths)^(num genes TF
					// regulates on other paths)
					// * (1 - TF activity prior)
					// Calculations need to be done in logspace, otherwise if there are too
					// many bound genes the score can go to 0 or NaN incorrectly
					double scoreNumer = Math.log(dProbBindingFunctional) *
							nRegThisPath[nPrimPathInd];
					scoreNumer += Math.log(1 - dProbBindingFunctional) *
							nRegOtherPaths[nPrimPathInd];
					// TODO: multiply activity prior by magnitude of mirna score???
					scoreNumer += Math.log(bindingData.regPriors[ntf]);
					double scoreDenom = Math.log(
							treeptr.nextptr[nPrimPathInd].numPath
									/ ((double) treeptr.numPath)) *
							nRegThisPath[nPrimPathInd];
					scoreDenom += Math.log(nNumOtherPaths
							/ ((double) treeptr.numPath)) *
							nRegOtherPaths[nPrimPathInd];
					scoreDenom += Math.log(1 - bindingData.regPriors[ntf]);
					// Take out of logspace for display and output
					treeptr.dTFActivityScore[ntf] = Math.exp(scoreNumer - scoreDenom);
					
					// Max activity score is either the new TF activity score at
					// this node
					// or the pre-existing max activity score
					dMaxTFActivityScore[ntf] = Math.max(
							treeptr.dTFActivityScore[ntf],
							dMaxTFActivityScore[ntf]);

					// Add the activity score at this node to the TreeSet
					// Instead of storing the score directly, store it in a form
					// that resembles
					// a pvalue so that the existing signifance score threshold
					// and sorting
					// code work
					double pseudoPvalue = Double.POSITIVE_INFINITY;
					if (treeptr.dTFActivityScore[ntf] != 0) {
						pseudoPvalue = 1 / treeptr.dTFActivityScore[ntf];
					}
					// Also with the TF name, store the index of the primary
					// path for that
					// TF coming out of the split
					// TODO store the path index separately
					treeptr.tsSigTFActivity[nPrimPathInd]
							.add(new SigTFRecv2(bindingData.regNames[ntf], pseudoPvalue));
				} else // treeptr.numchilren is 0 or 1 so there is no split to
				// calculate the score
				{
					treeptr.dTFActivityScore[ntf] = 0;
					// No need to add a SigTFRecv2 object to the tsSigTFActivity
					// TreeSet
					// if we are not at a split
				}

			} // end loop through all TFs

			int numrowstable;
			int nsize = bindingData.signedBindingValuesUnsorted.size() - 1;
			numrowstable = bindingData.regNames.length * nsize;
			
			boolean bsplit;
			if (treeptr.numchildren >= 2) {
				treeptr.dpvalEdgeSplit = new double[treeptr.numchildren][numrowstable];
				treeptr.dIDEdgeSplit = new Integer[treeptr.numchildren][numrowstable];
				treeptr.dscoreEdgeSplit = new Integer[treeptr.numchildren][numrowstable];
				treeptr.ddiffEdgeSplit = new double[treeptr.numchildren][numrowstable];
				treeptr.dexpectEdgeSplit = new double[treeptr.numchildren][numrowstable];
				treeptr.tsSigTFEdgeSplit = (TreeSet<SigTFRecv2>[]) new TreeSet[treeptr.numchildren];
				treeptr.dratioSplit = new double[treeptr.numchildren][numrowstable];

				for (int nindex = 0; nindex < treeptr.numchildren; nindex++) {
					treeptr.tsSigTFEdgeSplit[nindex] = new TreeSet<SigTFRecv2>(
							new SigTFRecv2Compare());
				}
				bsplit = true;
			} else {
				bsplit = false;
			}

			treeptr.dexpectEdgeFull = new double[treeptr.numchildren][0];
			treeptr.tsSigTFFull = (TreeSet<SigTFRecv2>[]) new TreeSet[treeptr.numchildren];
			
			Double[] minposi;
			Double[] minnega;

				if(minposipvalue.get(treeptr.ndepth) == null){
					minposi = new Double[bindingData.regNames.length];
					for(int i=0;i<minposi.length;i++) minposi[i] = 1.0;
				}else{
					minposi = minposipvalue.get(treeptr.ndepth);
				}
				if(minnegapvalue.get(treeptr.ndepth) == null){
					minnega = new Double[bindingData.regNames.length];
					for(int i=0;i<minnega.length;i++) minnega[i] = 1.0;
				}else{
					minnega = minnegapvalue.get(treeptr.ndepth);
				}
			
			
				for (int nchild = 0; nchild < treeptr.numchildren; nchild++) {
					int nrowindex = 0;
					for (int nrow = 0; nrow < bindingData.regNames.length; nrow++) {
						int nsumthispath = treeptr.nextptr[nchild].numPath;
						for (int nel = 0; nel < bindingData.signedBindingValuesSorted.length; nel++) {
							if (bindingData.signedBindingValuesSorted[nel] != 0) {
							    if (bsplit) {
									int nsumthisel = treeptr.ncountTotals[nrow][nel];
									int nsumthisparentpath = treeptr.numPath;
									double dexpect = (nsumthispath * nsumthisel)/(double) nsumthisparentpath;
									treeptr.dexpectEdgeSplit[nchild][nrowindex] = dexpect;
									double ddiff = treeptr.ncountvals[nrow][nchild][nel]-dexpect;
									treeptr.ddiffEdgeSplit[nchild][nrowindex] = ddiff;
									treeptr.dpvalEdgeSplit[nchild][nrowindex] = StatUtil
											.hypergeometrictail(
													treeptr.ncountvals[nrow][nchild][nel] - 1,
													nsumthisel, nsumthisparentpath - nsumthisel,
													nsumthispath);
									treeptr.dIDEdgeSplit[nchild][nrowindex] = nrow;
									treeptr.dscoreEdgeSplit[nchild][nrowindex] = bindingData.signedBindingValuesSorted[nel];
									treeptr.dratioSplit[nchild][nrowindex] = (double) treeptr.ncountvals[nrow][nchild][nel]
											/ treeptr.ncountTotals[nrow][nel];

										if(minposi[nrow]>treeptr.dpvalEdgeSplit[nchild][nrowindex]){
											minposi[nrow] = treeptr.dpvalEdgeSplit[nchild][nrowindex];
										}
									
									if (bindingData.signedBindingValuesUnsorted.size() == 2) {
										//TODO: what does this code do MARCEL
										String name = bindingData.regNames[nrow];
										if (bindingData.regTypes[nrow] == RegulatorBindingData.MIRNA
												&& treeptr.tranC.traindataTF[nrow].length != 0
												&& treeptr.tranC.traindataTF[nrow][0] != 0) {
											if (reg2DataSetIndex.get(name.toUpperCase()) != null) {
												treeptr.tsSigTFEdgeSplit[nchild]
														.add(new SigTFRecv2(
																name,
																treeptr.dpvalEdgeSplit[nchild][nrowindex],
																treeptr.dratioSplit[nchild][nrowindex]));
											}

										} else {
											if (bindingData.regTypes[nrow] == RegulatorBindingData.TF) {
												treeptr.tsSigTFEdgeSplit[nchild]
														.add(new SigTFRecv2(
																name,
																treeptr.dpvalEdgeSplit[nchild][nrowindex],
																treeptr.dratioSplit[nchild][nrowindex]));
											}
										}
									} else {
										treeptr.tsSigTFEdgeSplit[nchild]
												.add(new SigTFRecv2(
														bindingData.regNames[nrow]
																+ " "
																+ bindingData.signedBindingValuesSorted[nel],
														treeptr.dpvalEdgeSplit[nchild][nrowindex],
														treeptr.dratioSplit[nchild][nrowindex]));
									}
									nrowindex++;
							    }
							}
						}
					}
				}
				
			minposipvalue.put(treeptr.ndepth, minposi);
			minnegapvalue.put(treeptr.ndepth, minnega);
			Countsubpath++;
			runtext.append(Countsubpath+" nodes done..."+"  ");
			runtext.paintImmediately(runtext.getBounds());
			for (int nchild = 0; nchild < treeptr.numchildren; nchild++) {
				computeStats(treeptr.nextptr[nchild], rootptr);
			}
		}
	}

	// ///////////////////////////////////////////////////////////////////////////////
	/**
	 * Record used in the ranking of transcription factors
	 */
	static class SigTFRecv2 {
		String szname;
		double dpval;
		double dpercent;

		SigTFRecv2(String szname, double dpval, double dpercent) {
			this.szname = szname;
			this.dpval = dpval;
			this.dpercent = dpercent;
		}

		SigTFRecv2(String szname, double dpval) {
			this.szname = szname;
			this.dpval = dpval;
			this.dpercent = -1;
		}
	}

	/**
	 * Compares TFs first based on the value first lower dpval, then great
	 * dpercent, then name
	 */
	static class SigTFRecv2Compare implements Comparator<SigTFRecv2> {
		public int compare(SigTFRecv2 cr1, SigTFRecv2 cr2) {
			if (cr1.dpval < cr2.dpval)
				return -1;
			else if (cr1.dpval > cr2.dpval)
				return 1;
			else if (cr1.dpercent > cr2.dpercent)
				return -1;
			else if (cr1.dpercent < cr2.dpercent)
				return 1;
			else
				return cr1.szname.compareTo(cr2.szname);

		}
	}

	// ///////////////////////////////////////////////////////////////////////////////////////
	/**
	 * Computes the fraction of combinations for vals that could lead to a
	 * greater difference than ddiff
	 */
	public double computepval(int[][] vals, double diff) {
		int[] rowSum = new int[vals[0].length];
		int nsum0 = 0;
		int nsum1 = 0;

		for (int nindex = 0; nindex < vals[0].length; nindex++) {
			nsum0 += vals[0][nindex];
			nsum1 += vals[1][nindex];
			rowSum[nindex] = vals[0][nindex] + vals[1][nindex];
		}

		double dpvalsum = 0;
		int[] ncol0 = new int[vals[0].length];

		if (ncol0.length == 2) {
			for (ncol0[0] = 0; ncol0[0] <= Math.min(rowSum[0], nsum0); ncol0[0]++) {
				double dinnersum = 0;
				ncol0[1] = nsum0 - ncol0[0];

				if (ncol0[1] <= rowSum[1]) {
					double dweightsum0 = 0;
					double dweightsum1 = 0;
					for (int i = 0; i < ncol0.length; i++) {
						dweightsum0 += bindingData.signedBindingValuesSorted[i]*ncol0[i];
						dweightsum1 += bindingData.signedBindingValuesSorted[i]*(rowSum[i] - ncol0[i]);
					}
					double dcurrdiff = Math.abs(dweightsum0 / nsum0
							- dweightsum1 / nsum1);
					if (dcurrdiff >= diff) {
						dinnersum++;
					}
				}
				dpvalsum += StatUtil.hypergeometric(ncol0[0], nsum0, nsum1,
						rowSum[0])
						* dinnersum;
			}
		} else if (ncol0.length == 3) {
			for (ncol0[0] = 0; ncol0[0] <= Math.min(rowSum[0], nsum0); ncol0[0]++) {
				double dinnersum = 0;
				int nremainder0 = nsum0 - ncol0[0];
				int nremainder1 = nsum1 - (rowSum[0] - ncol0[0]);
				for (ncol0[1] = 0; ncol0[1] <= Math.min(rowSum[1], nremainder0); ncol0[1]++) {
					ncol0[2] = nsum0;
					for (int i = 0; i < ncol0.length - 1; i++) {
						ncol0[2] -= ncol0[i];
					}
					if (ncol0[2] <= rowSum[2]) {
						double dweightsum0 = 0;
						double dweightsum1 = 0;
						for (int i = 0; i < ncol0.length; i++) {
							dweightsum0 += bindingData.signedBindingValuesSorted[i]*ncol0[i];
							dweightsum1 += bindingData.signedBindingValuesSorted[i]*(rowSum[i] - ncol0[i]);
						}
						double dcurrdiff = Math.abs(dweightsum0 / nsum0
								- dweightsum1 / nsum1);
						if (dcurrdiff >= diff) {
							dinnersum += StatUtil.hypergeometric(ncol0[1],
									nremainder0, nremainder1, rowSum[1]);
						}
					}
				}
				dpvalsum += StatUtil.hypergeometric(ncol0[0], nsum0, nsum1,
						rowSum[0])
						* dinnersum;
			}
		} else if (ncol0.length == 1) {
			return 1;
		} else {
			throw new IllegalArgumentException(
					"can only handle two or three unique inputs");
		}

		return dpvalsum;
	}

	// ///////////////////////////////////////////////////////////////////////

	/**
	 * Recursively determines the best path through the model and its likelihood
	 */
	public double computevlogit(double[] vals, 
			double[][][] theInstance, int[][][] theInstanceIndex, int ndepth,
			Treenode node, int[] bestpath, int instanceindex) throws Exception {
		double dmax = Double.NEGATIVE_INFINITY;
		int bestchild = 0;
		double dlogout = 0;
		double dval = 0;

			if (node.dsigma > 0) {
				if (ndepth == 0) {
					dlogout = 0;
				} else {
					dlogout = Math.log(StatUtil.normaldensity(vals[ndepth], 
							node.dmean, node.dsigma));
				}
			}

		if (node.numchildren == 0) { //直到没有子节点情况，输出dlogout
			return dlogout;
		}

		int[] currbestpath = new int[bestpath.length];

		// log (a) + c(log(b_1) + log(b_2) + ... + log(b_n))
		// = log(a*(b_1*b_2*...*b_n)^c)

		double[] ptranlogit;

		if (BREGTSMiner) {
			ptranlogit = node.ptrans;
		} else if ((node.numchildren > 1) && (node.binit)) {
			ptranlogit = node.tranC.distributionForInstance(
					theInstanceIndex[node.ndepth][instanceindex],
					theInstance[node.ndepth][instanceindex]);
		} else {
			ptranlogit = CONSTANTA[node.numchildren];
		}

		for (int nchild = 0; nchild < node.numchildren; nchild++) {
			dval = dlogout;
			if (node.nextptr[nchild] != null) {//从根节点延伸到末节点的概率和赋给dval
				double dvi = computevlogit(vals, theInstance,
						theInstanceIndex, ndepth + 1, node.nextptr[nchild],
						bestpath, instanceindex);
				if (ptranlogit[nchild] == 0) {
					dval += Math.log(MINPROB);
				} else {
					dval += Math.log(ptranlogit[nchild]) + dvi;
				}
			} //从根节点延伸到末节点的概率和赋给dval

			if (dval > dmax) {  //如果找到概率更大路径，用dmax存储当年最大概率，currbestpath存储当前最佳路径
				dmax = dval;
				bestchild = nchild;
				for (int nindex = 0; nindex < bestpath.length; nindex++) {
					currbestpath[nindex] = bestpath[nindex];
				}
			}
		}

		for (int nindex = 0; nindex < bestpath.length; nindex++) {
			bestpath[nindex] = currbestpath[nindex];
		}
		bestpath[ndepth] = bestchild;
		return dmax;
	}


	/**
	 * Returns the number of nodes in the tree pointed by root for which its
	 * ncurrtime field does not equal the ncurrtime parameter
	 * 为所有节点的ncurrenttime参数赋值为目前当次的训练次数，并统计新训练的节点数！！
	 */
	public int countNodes(Treenode root, int ncurrtime) {

		if (root == null) {
			return 0;
		} else {
			int nsum;
			if (root.ncurrtime != ncurrtime) {
				root.ncurrtime = ncurrtime;
				nsum = 1;
			} else {
				nsum = 0;
			}

			for (int nchild = 0; nchild < root.numchildren; nchild++) {
				nsum += countNodes(root.nextptr[nchild], ncurrtime);
			}

			return nsum;
		}
	}

	/**
	 * Implements a baum-welch method for training the parameters with the
	 * current model structure
	 */
	public double trainhmm(Treenode treehmm, boolean bpruneexempt)throws Exception {

		double dlike = 0;
		double doldlike;
		initAE(treehmm);
		nglobaliteration = 1;
		dlike = 0;
		double[] ptranlogit;
		
		for (int nrow = 0; nrow < traindata.length; nrow++) {
			treehmm.df = 1;

			double[] theInstance = trainBinding.gene2RegBinding[treehmm.ndepth][nrow]; //nrow行基因对应的tf序列
			int[] theInstanceIndex = trainBinding.gene2RegBindingIndex[treehmm.ndepth][nrow];

			if (BREGTSMiner) {
				ptranlogit = treehmm.ptrans;	
			} else if ((treehmm.numchildren > 1) && (treehmm.binit)) {
				ptranlogit = treehmm.tranC.distributionForInstance(
						theInstanceIndex, theInstance);
			} else {
				ptranlogit = CONSTANTA[treehmm.numchildren];
			}
			
			for (int nchild = 0; nchild < treehmm.numchildren; nchild++) {
			forwardalg(traindata[nrow],              
						trainBinding.gene2RegBindingIndex,
						trainBinding.gene2RegBinding, 1, ptranlogit[nchild],
						treehmm.nextptr[nchild], nrow, treehmm, nrow);
			}
			
			backalg(traindata[nrow], 
					trainBinding.gene2RegBindingIndex,
					trainBinding.gene2RegBinding, 0, treehmm, nrow, true, nrow);

			double dpj = treehmm.db;
			
			if (dpj == 0) {
				dlike += Math.log(MINPROB);
			} else {
				instanceAE(treehmm, 0, traindata[nrow], dpj, nrow); 
				dlike += Math.log(treehmm.db);
			}
		}
		
		updateParams(treehmm, 0, traindata, 0);
		
		int ncount = 0;			
		double dpredictitr;	
		
		do {                     
			nglobaliteration++;
			if (BDEBUG) {
				System.out.println("Global iteration is now " + nglobaliteration);
			}

			doldlike = dlike;
			initAE(treehmm);
			dlike = 0;
			
			for (int nrow = 0; nrow < traindata.length; nrow++) {
				double[] theInstance = trainBinding.gene2RegBinding[treehmm.ndepth][nrow];
				int[] theInstanceIndex = trainBinding.gene2RegBindingIndex[treehmm.ndepth][nrow];

				if (BREGTSMiner) { //不加入tf信息，BREGTSMiner为TRUE
					ptranlogit = treehmm.ptrans;
				} else if ((treehmm.numchildren > 1) && (treehmm.binit)) {
					ptranlogit = treehmm.tranC.distributionForInstance(
							theInstanceIndex, theInstance);
				} else {
					ptranlogit = CONSTANTA[treehmm.numchildren];
				}

				treehmm.df = 1;
				for (int nchild = 0; nchild < treehmm.numchildren; nchild++) {
					forwardalg(traindata[nrow], 
							trainBinding.gene2RegBindingIndex,
							trainBinding.gene2RegBinding, 1,
							ptranlogit[nchild], treehmm.nextptr[nchild], nrow,
							treehmm, nrow);
				}
				    backalg(traindata[nrow], 
						trainBinding.gene2RegBindingIndex,
						trainBinding.gene2RegBinding, 0, treehmm, nrow, true, nrow); //节点的后向概率

				double dpj = treehmm.db;

				if (dpj == 0) {
					dlike += Math.log(MINPROB);
				} else {
					instanceAE(treehmm, 0, traindata[nrow], dpj, nrow);
					dlike += Math.log(treehmm.db);
				}
				
				if (Double.isNaN(dlike)) {
					traverse(treehmm, 0, true);
					System.out.println(nrow + "\t");
					if (!BREGTSMiner) {
						for (int na = 0; na < treehmm.tranC.dcoeff.length; na++) {
							System.out.print(treehmm.tranC.dcoeff[na] + "\t");
						}
					}
					System.out.println();
					System.out.println("treehmm.db = " + treehmm.db);
					throw new Exception();
				}
			}

			updateParams(treehmm, 0, traindata, 0);

			traverse(treehmm, 0, false);
			ncount++;
			dpredictitr = (dbesttrainlike - dlike) / (dlike - doldlike);
			if (BDEBUGMODEL) {
				System.out.println(ncount + " - Train hmm Likelihood: " + dlike
						+ " best " + dbesttrainlike + " " + "oldlike"
						+ doldlike + " " + dpredictitr);
			}
			
		} while (((dlike - doldlike) / Math.abs(dlike) > BEPSILON) && ((bpruneexempt) || (dpredictitr < MAXFUTUREITR)));
		
		if (BDEBUG) {
			System.out.println("Likelihood: " + dlike);
		}

		dtrainlike = dlike;  //记录训练得到的最优likehood
		nglobaltime++;
		int nparams = countNodes(treehmm, nglobaltime); 
		
		if (BDEBUG) {System.out.print(dlike + "\t" + nparams);}
		
		dlike -= nodepenalty * nparams;

		if (BDEBUG) {
			System.out.println("\t" + dlike);
		}
		
		return (dlike);
	}
}