package bprc.nrem;

import bprc.core.*;

import java.io.*;

/**
 * Class extends DataSetCore which contains the underlying data and parameters
 * of the methods with an instance of TSMiner_GoAnnotations
 */
public class TSMiner_DataSet extends DataSetCore {

	TSMiner_GoAnnotations tga;

	/**
	 * Constructor that does a simple copy the contents of theDataSetCore and
	 * tga
	 */
	public TSMiner_DataSet(DataSetCore theDataSetCore, TSMiner_GoAnnotations tga) {
		super(theDataSetCore);
		this.tga = tga;
	}

	/**
	 * Constructor that takes input parameters and calls dataSetReader to read
	 * in the content of szInputFile
	 */
	public TSMiner_DataSet(String szInputFile, int nmaxmissing,
			double dthresholdvalue, boolean btakelog,
			boolean bspotincluded, boolean brepeatset, boolean badd0,
			boolean bmaxminval, boolean bfcto0, boolean bfctopre)
			throws IOException, FileNotFoundException, IllegalArgumentException {
		
		this.szInputFile = szInputFile;
		this.nmaxmissing = nmaxmissing;
		this.dthresholdvalue = dthresholdvalue;
		this.bmaxminval = bmaxminval;
		this.bfcto0 = bfcto0;
		this.bfctopre = bfctopre;
		this.btakelog = btakelog;
		this.bspotincluded = bspotincluded;
		this.badd0 = badd0;

		//调用datasetreader方法
		dataSetReader(szInputFile, nmaxmissing, dthresholdvalue,
				btakelog, bspotincluded, brepeatset, badd0);
	}
}