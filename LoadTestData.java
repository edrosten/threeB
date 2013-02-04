import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.plugin.*;
import ij.plugin.frame.*;
import ij.plugin.filter.PlugInFilter;
import ij.*;
import ij.io.*;
import ij.plugin.*;
import ij.plugin.filter.*;

import java.awt.event.*;
import java.awt.geom.*;
import java.util.*;
import java.io.*;
import javax.swing.*;
import javax.swing.event.*;
import java.lang.InterruptedException;
import java.lang.System;
import java.lang.Math.*;
import java.net.*;


///Plugin class to load the standard test data from the JAR file.
///@ingroup gPlugin
public class LoadTestData implements PlugIn {

	public void run(String arg)
	{
		Opener o = new Opener();
		InputStream s = getClass().getClassLoader().getResourceAsStream("test_data.tif");
		ImagePlus im = o.openTiff(s, "test_data.tiff");


		try{
				s.close();
		}
		catch(IOException close_err){
			Toolkit.getDefaultToolkit().beep();
			ij.IJ.showStatus("Error reading image.");
			return;
		}

		im.show();
		im.updateAndDraw();


		System.out.println(getClass().getResource("").toString() + "\n");
		System.out.println(getClass().getResource("").getPath().toString() + "\n");
		System.out.println(getClass().getResource("").getFile().toString() + "\n");

	}

};
