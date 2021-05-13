# TSMiner


## What is it?
TSMiner is a software program to construct time-specific regulatory networks from time-series expression profiles using two groups of genes:  
(1) genes encoding transcription factors (TFs) that are activated or repressed at a specific time;  
(2) genes involved in biological pathways showing significant mutual interactions with these TFs.  
TSMiner provides a network map with extensive interactive visualization operations to help users explore the time-specific regulatory networks.The core algorithms and interfaces of TSMiner were written in java language.

## The release version

### TSMiner version v1.0.0 (11/28/2018)
The first release version of TSMiner.

## Software requirements

- Java 1.5 or later must be installed.  [Download link](http://www.java.com)
- 1024 MB Java heap memory minimum.
- Operating systems: Platform independent

## Running

- Once Java is installed, to start TSMiner in Microsoft Windows double click on TSMiner.jar.  
- Otherwise start TSMiner from a command line while in the drem directory, enter the command: java -mx1024M  -jar TSMiner.jar
- The data input and parameter settings are described in TSMiner manual.pdf.
- Download link: 

## Files

- Source code and Javadoc api for the source code are included in the 'sourcecode' directory.
- License files are included in 'license' directory.
- Gene ontology files and KEGG/Reactome pathway files are included in the 'Gene annotation file' directory.
- Example TF-gene interaction data are included in 'TF-gene file' directory.
- Example time-series expression data sets are included in 'Time-series data' directory.
- Example saved models are included in 'Saved model' directory.

## Example

**To run a case analysis of mouse liver regeneration (LR) with default settings.**

- If in Microsoft Windows double click on the file 'Case\_analysis\_of\_mouseLR.cmd' and press the 'Run' button in the main interface.
- Otherwise from a command line while in the TSMiner directory type the command: java -mx1024M -jar TSMiner.jar -d Case\_analysis\_of\_mouseLR.txt, then press the 'Run' button.

**The default parameters are shown in Defaults\_for\_case\_analysis.txt.**

##  Contact

  For any questions please contact Dr. Mingfei Han (Email: [free1234hm@163.com](mailto:free1234hm@163.com)).

##  License

- This software product is developed by the National Center for Protein Sciences (Beijing)-Bioinformatics group. The TSMiner software makes use of several existing tools: 
- the Dynamic Regulatory Events Miner  (DREM) ([Mol Syst Biol. 2007;3:74. Epub 2007 Jan 16.](https://www.ncbi.nlm.nih.gov/pubmed/17224918)) which is distributed GPL v3.0 license 
- the Piccolo toolkit ([www.cs.umd.edu/hcil/piccolo](https://www.cs.umd.edu/hcil/piccolo)) which is distributed under a BSD license
- the Batik software ([http://xmlgraphics.apache.org/batik/](http://xmlgraphics.apache.org/batik/)) which is distributed under an Apache license
- Please review and agree to the DREM, Piccolo, and Batik licenses before using the TSMiner software. The detailed licenses are described in the license folder.
