import ij.plugin.*;
import ij.plugin.filter.*;
import ij.gui.*;

import java.awt.event.*;
import java.awt.geom.*;
import java.util.*;
import java.io.*;
import javax.swing.*;
import javax.swing.event.*;
import java.lang.InterruptedException;
import java.lang.System;
import java.lang.Math.*;

///3B plugin cclass to bring up a basic help window.
///@ingroup gPlugin
public class ThreeBHelp implements PlugIn {

	public void run(String arg)
	{
		System.out.println(arg);

		if(arg.equals("about"))
		{
			JFrame w = new JFrame("About 3B");

			JEditorPane l = new JEditorPane();
			try{
				l.setPage(getClass().getResource("about.txt").toString());
			}
			catch(IOException e)
			{}

			JScrollPane s = new JScrollPane(l, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);

			java.awt.Dimension  d = new java.awt.Dimension(600, 600);

			w.add(s);
			w.pack();
			w.setSize(d);
			w.setVisible(true);


		}
		else
		{
			JFrame w = new JFrame("3B help");

			JEditorPane l = new JEditorPane();
			try{
				l.setPage(getClass().getResource("instructions.html").toString());
			}
			catch(IOException e)
			{}

			JScrollPane s = new JScrollPane(l, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);

			java.awt.Dimension  d = new java.awt.Dimension(600, 600);

			w.add(s);
			w.pack();
			w.setSize(d);
			w.setVisible(true);
		}
		
	}

};
