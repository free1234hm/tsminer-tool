package bprc.nrem;

import bprc.core.*;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;

import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.*;
import java.util.List;
import java.text.NumberFormat;
import java.text.ParseException;
import java.awt.*;
import java.awt.event.*;
import java.io.*;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.table.*;

import java.awt.datatransfer.*;

/**
 * Class encapsulates a gene table which is used to display the selected genes
 * in TSMiner
 */
public class GeneTable extends JPanel implements ActionListener {
	TSMiner_DataSet theDataSet;
	String[] columnNames;
	String[][] tabledata;
	double[] davg;
	double[] dstd;
	JTable theTable;
	JScrollPane scrollPane;
	JButton saveButton;
	JButton copyButton;
	JButton savenamesButton, copynamesButton;
	TableSorter sorter;
	final static Color bgColor = Color.white;
	final static Color fgColor = Color.black;
	JFrame theFrame;
	TSMiner_Timeiohmm.Treenode node;

	/**
	 * Constructor builds the gene table interface
	 */
	public GeneTable(JFrame theFrame, TSMiner_DataSet theDataSet,
			TSMiner_Timeiohmm.Treenode node) {
		this.theFrame = theFrame;
		this.theDataSet = theDataSet;
		this.node = node;
		geneTableHelper();
	}

	/**
	 * Helper method to build the gene table
	 */
	private void geneTableHelper() {
		setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
		setBackground(bgColor);
		setForeground(fgColor);
		int numcols = theDataSet.numcols + 2;
		columnNames = new String[numcols];
		columnNames[0] = theDataSet.szGeneHeader;
		columnNames[1] = theDataSet.szProbeHeader;

		for (int ncolindex = 0; ncolindex < theDataSet.numcols; ncolindex++) {
			columnNames[ncolindex + 2] = "" + theDataSet.dsamplemins[ncolindex];
		}

		NumberFormat nf2 = NumberFormat.getInstance(Locale.ENGLISH);
		nf2.setMinimumFractionDigits(2);
		nf2.setMaximumFractionDigits(2);
		
		List<Integer> list = node.genelist;
		tabledata = new String[list.size()][columnNames.length];

		davg = new double[theDataSet.numcols];
		dstd = new double[theDataSet.numcols];
		double[] dsum = new double[theDataSet.numcols];
		double[] dsumsq = new double[theDataSet.numcols];

		for (int i = 0; i < tabledata.length; i++) {
			int ngene = list.get(i);
			tabledata[i][0] = theDataSet.genenames[ngene];
			tabledata[i][1] = theDataSet.probenames[ngene];
			for (int ncol = 0; ncol < theDataSet.numcols; ncol++) {
				dsum[ncol] += theDataSet.data[ngene][ncol];
				dsumsq[ncol] += theDataSet.data[ngene][ncol] * theDataSet.data[ngene][ncol];
				tabledata[i][ncol + 2] = nf2.format(theDataSet.data[ngene][ncol]);
			}
		}

		for (int ni = 0; ni < davg.length; ni++) {
			davg[ni] = dsum[ni] / list.size();
			dstd[ni] = Math.sqrt((dsumsq[ni] - list.size() * davg[ni] * davg[ni])
							/ (list.size() - 1));
		}
		
		sorter = new TableSorter(new TableModelST(tabledata, columnNames));
		theTable = new JTable(sorter);

		sorter.setTableHeader(theTable.getTableHeader());
		theTable.setPreferredScrollableViewportSize(new Dimension(800, 
				Math.min((theTable.getRowHeight() + theTable.getRowMargin()) * theTable.getRowCount(), 400)));

		theTable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
		
		TableColumn column;
		column = theTable.getColumnModel().getColumn(0);
		column.setPreferredWidth(150);
		column = theTable.getColumnModel().getColumn(1);
		column.setPreferredWidth(100);

		for (int ncolindex = 0; ncolindex < theDataSet.numcols; ncolindex++) {
			column = theTable.getColumnModel().getColumn(ncolindex + 2);
			column.setPreferredWidth(45);
		}

		int preferredwidth = Math.min(280+theDataSet.numcols*45, 800);
		theFrame.setPreferredSize(new Dimension(preferredwidth, 600));
		
		// Create the scroll pane and add the table to it.
		scrollPane = new JScrollPane(theTable);
		// Add the scroll pane to this panel.
		add(scrollPane);

		JPanel countPanel = new JPanel();
		String szcountLabel = "Total number of genes selected is " + tabledata.length;
		JLabel countLabel = new JLabel(szcountLabel);
		countPanel.setBackground(Color.white);
		countPanel.add(countLabel);
		countPanel.setMaximumSize(new Dimension(Integer.MAX_VALUE, 20));
		add(countPanel);

		JPanel avgPanel = new JPanel();
		String szavgLabel = "Average expression (" + nf2.format(davg[0]);
		for (int ni = 1; ni < davg.length; ni++) {
			szavgLabel += ", " + nf2.format(davg[ni]);
		}
		szavgLabel += ")";
		JLabel avgLabel = new JLabel(szavgLabel);
		avgPanel.setBackground(Color.white);
		avgPanel.add(avgLabel);
		avgPanel.setMaximumSize(new Dimension(Integer.MAX_VALUE, 20));
		add(avgPanel);

		JPanel stdPanel = new JPanel();
		String szstdLabel = "Standard Deviation expression (" + nf2.format(davg[0]);
		for (int ni = 1; ni < davg.length; ni++) {
			szstdLabel += ", " + nf2.format(dstd[ni]);
		}
		szstdLabel += ")";
		JLabel stdLabel = new JLabel(szstdLabel);
		stdPanel.setBackground(Color.white);
		stdPanel.add(stdLabel);
		stdPanel.setMaximumSize(new Dimension(Integer.MAX_VALUE, 20));
		add(stdPanel);

		copyButton = new JButton("Copy Table", Util
				.createImageIcon("Copy16.gif"));
		copyButton.setActionCommand("copy");
		copyButton.setMinimumSize(new Dimension(800, 20));
		copyButton.addActionListener(this);

		saveButton = new JButton("Save Table", Util
				.createImageIcon("Save16.gif"));
		saveButton.setActionCommand("save");
		saveButton.setMinimumSize(new Dimension(800, 20));
		saveButton.addActionListener(this);

		JPanel buttonPanel = new JPanel();
		buttonPanel.setBackground(Color.white);
		buttonPanel.add(copyButton);
		buttonPanel.add(saveButton);

		copynamesButton = new JButton("Copy Gene Names", Util
				.createImageIcon("Copy16.gif"));
		copynamesButton.setActionCommand("copynames");
		copynamesButton.setMinimumSize(new Dimension(800, 20));
		copynamesButton.addActionListener(this);

		savenamesButton = new JButton("Save Gene Names", Util.createImageIcon("Save16.gif"));
		savenamesButton.setActionCommand("savenames");
		savenamesButton.setMinimumSize(new Dimension(800, 20));
		savenamesButton.addActionListener(this);
		buttonPanel.add(copynamesButton);
		buttonPanel.add(savenamesButton);
		
		JButton helpButton = new JButton(Util.createImageIcon("Help16.gif"));
		helpButton.addActionListener(this);
		helpButton.setActionCommand("help");
		buttonPanel.add(helpButton);

		buttonPanel.setMaximumSize(new Dimension(Integer.MAX_VALUE, 20));
		add(buttonPanel);
	}

	/**
	 * Outputs the content of the PrintWriter
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
	 * Copies the contents of the gene table to the clipboard
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
	 * Outputs the list of genes in the table to a file
	 */
	public void printGeneList(PrintWriter pw) {
		for (int nrow = 0; nrow < tabledata.length; nrow++) {
			pw.println(sorter.getValueAt(nrow, 0));
		}
	}

	/**
	 * Copies the gene names to the clipboard
	 */
	public void writenamesToClipboard() {
		StringBuffer sbuf = new StringBuffer();
		// get the system clipboard
		for (int nrow = 0; nrow < tabledata.length; nrow++) {
			sbuf.append(sorter.getValueAt(nrow, 0) + "\n");
		}

		Clipboard systemClipboard = Toolkit.getDefaultToolkit()
				.getSystemClipboard();
		Transferable transferableText = new StringSelection(sbuf.toString());
		systemClipboard.setContents(transferableText, null);
	}

	/**
	 * Responds to actions on the interface
	 */
	public void actionPerformed(ActionEvent e) {
		String szCommand = e.getActionCommand();

		if (szCommand.equals("copy")) {
			writeToClipboard();
		} else if (szCommand.equals("copynames")) {
			writenamesToClipboard();
		} else if ((szCommand.equals("save"))
				|| (szCommand.equals("savenames"))) {
			try {
				int nreturnVal = TSMiner_IO.theChooser.showSaveDialog(this);
				if (nreturnVal == JFileChooser.APPROVE_OPTION) {
					File f = TSMiner_IO.theChooser.getSelectedFile();
					PrintWriter pw = new PrintWriter(new FileOutputStream(f));
					if (szCommand.equals("save")) {
						printFile(pw);
					} else {
						printGeneList(pw);
					}
					pw.close();
				}
			} catch (FileNotFoundException ex) {
				final FileNotFoundException fex = ex;
				javax.swing.SwingUtilities.invokeLater(new Runnable() {
					public void run() {
						JOptionPane.showMessageDialog(null, fex.getMessage(),
								"Exception thrown", JOptionPane.ERROR_MESSAGE);
					}
				});
				ex.printStackTrace(System.out);
			}
		} else if (szCommand.equals("help")) {
			String szMessage = "The table contains the currently displayed genes ";

			String szlogword = "";
			if (theDataSet.btakelog) {
				szlogword = "log base two ";
			}
			szMessage += ".\n\n The first two columns of the table are: \n"
					+ "*  "
					+ theDataSet.szGeneHeader
					+ " - The name of the gene.\n"
					+ "*  "
					+ theDataSet.szProbeHeader
					+ " - The spot ID(s) associated with the gene.\n"
					+ "The remaining columns contain the "
					+ szlogword
					+ "expression change relative to the first time point.\n\n"
					+ "Note:\n"
					+ "+The table can be sorted by any of the columns by clicking on the column's header.\n"
					+ "+Using the 'Save Table' button the entire table can be saved, or just the gene names "
					+ "with the 'Save Gene Names' button.";
			Util.renderDialog(theFrame, szMessage);
		}
	}
}
