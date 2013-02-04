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
		String path = getClass().getResource("img/").toString();
		Reader htmlstream = new InputStreamReader(getClass().getClassLoader().getResourceAsStream("instructions.html"));
		String html = Util.read(htmlstream);

		html = html.replaceAll("img/", path);

		JFrame w = new JFrame("3B help");

		JEditorPane l = new JEditorPane();
		try{
			l.setPage(getClass().getResource("instructions.html").toString());
		}
		catch(IOException e)
		{}
		//l.resize(new java.awt.Dimension(800,600));

		JScrollPane s = new JScrollPane(l, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);

		java.awt.Dimension  d = new java.awt.Dimension(600, 600);

		w.add(s);
		w.pack();
		w.setSize(d);
		w.setVisible(true);
		
	}

};
