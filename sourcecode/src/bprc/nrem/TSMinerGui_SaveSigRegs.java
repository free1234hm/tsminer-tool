package bprc.nrem;

import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;

import javax.swing.BoxLayout;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;

public class TSMinerGui_SaveSigRegs extends JPanel{
	public TSMinerGui_SaveSigRegs(final TSMinerGui theTSMinerGui, final JFrame theFrame)
	{
		setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
		setForeground(Color.BLACK);
		final JFileChooser theChooser = new JFileChooser();
		add(theChooser);
		theChooser.setDialogType(JFileChooser.SAVE_DIALOG);
		theChooser.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				// set label's icon to the current image
				String state = (String) e.getActionCommand();

				if (state.equals(JFileChooser.CANCEL_SELECTION)) {
					theFrame.setVisible(false);
					theFrame.dispose();
				} else if (state.equals(JFileChooser.APPROVE_SELECTION)) {
					File f = theChooser.getSelectedFile();

					try {
						PrintWriter pw = new PrintWriter(
								new FileOutputStream(f));
						theTSMinerGui.printSigRegs(pw);
						pw.close();
					} catch (final IOException fex) {
						javax.swing.SwingUtilities.invokeLater(new Runnable() {
							public void run() {
								JOptionPane.showMessageDialog(null,
										fex.getMessage(), "Exception thrown",
										JOptionPane.ERROR_MESSAGE);
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
}
