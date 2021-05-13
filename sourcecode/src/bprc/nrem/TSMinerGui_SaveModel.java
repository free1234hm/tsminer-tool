package bprc.nrem;

import bprc.core.*;
import edu.umd.cs.piccolo.nodes.PPath;
import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;
import java.util.*;
import java.text.NumberFormat;
import java.io.*;

import javax.imageio.*;
import java.awt.image.*;
import edu.umd.cs.piccolo.nodes.PImage;

/**
 * Class to encapsulate window used to specify a file to save a TSMiner model
 */ 
public class TSMinerGui_SaveModel extends JPanel
{
	/** Use of this map needs to be synchronized */
	private static HashMap<String, String> synMap = null;
	
	final static Color bgColor = Color.white;
	final static Color fgColor = Color.black;
	final TSMiner_Timeiohmm theTimeiohmm;
	final TSMiner_Timeiohmm.Treenode treecopy;
	final JFileChooser theChooser;

	/**
	 * Class constructor
	 */
	public TSMinerGui_SaveModel(final TSMiner_Timeiohmm theTimeiohmm,final TSMiner_Timeiohmm.Treenode treecopy,
			final JFrame theFrame,final TSMinerGui theTSMinerGui)
	{
		this.theTimeiohmm = theTimeiohmm;
		this.treecopy = treecopy;

		setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
		setForeground(fgColor);
		theChooser = new JFileChooser();
		add(theChooser);
		theChooser.setDialogType(JFileChooser.SAVE_DIALOG);
		theChooser.addActionListener(new ActionListener() 
		{
			public void actionPerformed(ActionEvent e) 
			{
				// set label's icon to the current image
				String state = (String)e.getActionCommand();

				if (state.equals(JFileChooser.CANCEL_SELECTION))
				{
					theFrame.setVisible(false);
					theFrame.dispose();
				}
				else if (state.equals(JFileChooser.APPROVE_SELECTION))
				{
					File f = theChooser.getSelectedFile();

					try
					{
						PrintWriter pw = new PrintWriter(new FileOutputStream(f));
						pw.print(theTimeiohmm.saveString(treecopy));
						pw.println("COLORS");
						pw.print(theTSMinerGui.saveColors());
						pw.close();
						
						// Only need to save the activity scores when running STSMiner
						if(theTimeiohmm.regPriorsFile != null && !theTimeiohmm.regPriorsFile.equals(""))
						{
							// TODO Temporarily commented out these function calls until
							// the rest of the STSMiner code merge is complete.  Need a
							// more elegant way to save activity scores.
							// Made these static so that TSMiner_IO_Batch can use them as well
							/*
							saveActivityScoresDynamic(f, theTimeiohmm, treecopy);
							saveActivityScores(f, theTimeiohmm);
							*/
						}
					}
					catch (final IOException fex)
					{
						javax.swing.SwingUtilities.invokeLater(new Runnable() 
						{
							public void run() 
							{
								JOptionPane.showMessageDialog(null, fex.getMessage(), 
										"Exception thrown", JOptionPane.ERROR_MESSAGE);
							}
						});
						fex.printStackTrace(System.out);
					}
					theTSMinerGui.bsavedchange = true;
					theFrame.setVisible(false);
					theFrame.dispose();
				}
			}
		});			      
	}
	
	
	// TODO need a more permanent solution and error checking
	/**
	 * Uses the filename to create three files:
	 * [filename].model - the model file
	 * [filename].model.activities - the TF activities
	 * [filename].model.activitiesStd - the TF activities with their standard names
	 * @param filename
	 */
	public void batchSave(String filename)
	{
		theChooser.setSelectedFile(new File(filename + ".model"));
		theChooser.approveSelection();
	}
	
	// TODO when integrating STSMiner code, need to use the gene name synonym 
	// functionality built into TSMiner
	/**
	 * A temporary means for retrieving the max TF activity scores
	 * @param modelFile
	 */
	public synchronized static void saveActivityScores(File modelFile, TSMiner_Timeiohmm theTimeiohmm)
		throws IOException
	{
		if(synMap == null)
		{
			// No longer hard code a filename
			synMap = loadMap("", 1, 0, false);
		}
		
		String filename = modelFile.getAbsolutePath() + ".activities";
		PrintWriter writer = new PrintWriter(new FileWriter(filename));
		
		String filenameStd = modelFile.getAbsolutePath() + ".activitiesStd";
		PrintWriter writerStd = new PrintWriter(new FileWriter(filenameStd));
		
		for(int t = 0; t < theTimeiohmm.bindingData.regNames.length; t++)
		{
			double score = theTimeiohmm.dMaxTFActivityScore[t];

			String tf = theTimeiohmm.bindingData.regNames[t].toUpperCase();
			String stdName = tf;
			if(synMap != null && synMap.containsKey(tf))
			{
				String newName = synMap.get(tf);
				if(newName.contains("|"))
				{
					newName = newName.substring(0, newName.indexOf('|'));
				}
				tf = newName;
			}

			writer.println(tf + "\t" + score);
			writerStd.println(stdName + "\t" + score);
		}
		writer.close();
		writerStd.close();
	}
	
	// TODO when integrating STSMiner code, need to use the gene name synonym 
	// functionality built into TSMiner
	/**
	 * A temporary means for retrieving the max TF activity scores at
	 * each time point.  Must be called before saveActivityScores if at all
	 * because it does write its own flag.
	 * @param modelFile
	 */
	public synchronized static void saveActivityScoresDynamic(File modelFile, TSMiner_Timeiohmm theTimeiohmm
			, TSMiner_Timeiohmm.Treenode treecopy) throws IOException
	{
		if(synMap == null)
		{
			// No longer hard code a filename
			synMap = loadMap("", 1, 0, false);
		}
		
		String filename = modelFile.getAbsolutePath() + ".activitiesDynamic";
		PrintWriter writer = new PrintWriter(new FileWriter(filename));
		
		writer.print(theTimeiohmm.saveActivityScores(treecopy));
		
		writer.close();
	}
	
	// TODO when integrating STSMiner code, need to use the gene name synonym 
	// functionality built into TSMiner
	private synchronized static HashMap<String, String> loadMap(String file, int keyInd, int valInd, boolean header)
	{
		try
		{
			BufferedReader reader = new BufferedReader(new FileReader(file));
			HashMap<String, String> map = new HashMap<String, String>();
	
			if(header)
			{
				reader.readLine();
			}
	
			String line;
			while((line = reader.readLine()) != null)
			{
				String[] parts = line.split("\t");
				if(parts.length > Math.max(keyInd, valInd))
				{
					String key = parts[keyInd].toUpperCase();
					String val = parts[valInd].toUpperCase();
	
					if(!val.equals("-") && !val.equals(""))
					{
						if(map.containsKey(key))
						{
							String existing = (String) map.get(key);
							map.put(key, existing + "|" + val);
						}
						else
						{
							map.put(key, val);
						}
					}
				}
				else
				{
					System.err.println("Not enough columns: " + line);
				}
			}			
	
			reader.close();
	
			return map;
		}
		catch(IOException e)
		{
			e.printStackTrace();
			return null;
		}
	}
}