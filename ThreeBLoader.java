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
import java.util.zip.*;
import java.io.*;
import javax.swing.*;
import javax.swing.event.*;
import java.lang.InterruptedException;
import java.lang.System;
import java.lang.Math.*;

class XYPair
{
	public int x, y;

	@Override
	public boolean equals(Object o)
	{
		XYPair p = (XYPair)o;
		return p.x == x && p.y == y;
	}

	int hash(int x) 
	{     
		x = ((x >> 16) ^ x) * 0x45d9f3b;     
		x = ((x >> 16) ^ x) * 0x45d9f3b;     
		x = ((x >> 16) ^ x);     
		return x; 
	}

	@Override
	public int hashCode()
	{
		int seed = hash(x);
		seed ^= hash(y) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed;
	}
};

///Plugin class to load up an old 3B run.
///
///@ingroup gPlugin
public class  ThreeBLoader implements PlugIn {
	
	private void close(InputStream i)
	{
		//Not sure what to do here...
		if(i != null)
		{
			try{
				i.close();
			}
			catch(IOException ioe)
			{
				Toolkit.getDefaultToolkit().beep();
				ij.IJ.showStatus("Error closing file: " + ioe.getMessage());
			}
		}
	}

	public void run(String arg) {
		InputStream in = null;
		String name;
		double ps;
		double ini_FWHM=100., ini_reconstruction_pixel_size=10.;
		boolean show_control_panel=true;
		boolean do_filtering=true;

		try{

			if(arg.equals("test"))
			{
				name = "test_data.txt";

				try{
					in = new GZIPInputStream(getClass().getClassLoader().getResourceAsStream("test_data.txt.gz"));
				}
				catch(IOException ioe)
				{
					Toolkit.getDefaultToolkit().beep();
					ij.IJ.showStatus("Could not open test data: " + ioe.getMessage());
					return;
				}
				ps=100;
			}
			else
			{
				OpenDialog o = new OpenDialog("Open 3B run...", null);
				name = o.getFileName();
				
				try{
					in = new FileInputStream(o.getDirectory() + o.getFileName());
				}
				catch(java.io.FileNotFoundException ferr)
				{
					Toolkit.getDefaultToolkit().beep();
					ij.IJ.showStatus("Error opening file: " + ferr.getMessage());
					return;
				}



				GenericDialog g = new GenericDialog("Pixel size");
				g.addNumericField("Pixel size", 100., 0, 5, "nm");
				g.addNumericField("FWHM (initial) ", ini_FWHM, 0, 5, "nm");
				g.addNumericField("Reconstruction size (initial)", ini_reconstruction_pixel_size, 0, 5, "nm");
				g.addCheckbox("Do Filtering", do_filtering);
				g.addCheckbox("Show control panel", show_control_panel);
				
				((Checkbox)g.getCheckboxes().get(1)).hide();
				
				g.showDialog();


				ps = g.getNextNumber();

				ini_FWHM = g.getNextNumber();
				ini_reconstruction_pixel_size=g.getNextNumber();
				do_filtering = g.getNextBoolean();
				show_control_panel=g.getNextBoolean();
			}

			//Yay @ cargoculting
			InputStreamReader d = new InputStreamReader(in);
			BufferedReader r = new BufferedReader(d);

			String line;
			final ArrayList<Spot> spots = new ArrayList<Spot>();
			int iterations=0;
			Rectangle roi=null;
			final double pixel_size_in_nm_ = ps;

			HashSet<XYPair> filter_hash = new HashSet<XYPair>();

			while((line = r.readLine()) != null)
			{
				String tokens[] = line.split("\\p{Space}+");

				if(tokens[0].equals("PIXELS"))
				{
					if((tokens.length - 1)%2 != 0)
						throw new IOException("corrupt file (bad pixels line)");

					for(int i=1; i < tokens.length; i+= 2)
					{
						int x, y;
						try
						{
							x = Integer.parseInt(tokens[i+0]);
							y = Integer.parseInt(tokens[i+1]);
						}
						catch(NumberFormatException nerr)
						{
							throw new IOException("corrupt file (bad pixel coordinate)");
						}

						if(roi == null)
							roi = new Rectangle(x, y, 1, 1);
						else
							roi.add(x, y);
					}
				}

				if(tokens[0].equals("FILTER"))
				{
					if((tokens.length - 1)%2 != 0)
						throw new IOException("corrupt file (bad filter line)");
					filter_hash.clear();

					for(int i=1; i < tokens.length; i+= 2)
					{
						XYPair xy = new XYPair();
						try
						{
							xy.x = Integer.parseInt(tokens[i+0]);
							xy.y = Integer.parseInt(tokens[i+1]);
						}
						catch(NumberFormatException nerr)
						{
							throw new IOException("corrupt file (bad filter coordinate)");
						}
						
						filter_hash.add(xy);
					}

					System.out.println(filter_hash);
				}

				if(tokens[0].matches("PASS[0-9]+:"))
				{
					iterations++;

					if((tokens.length - 1)%4 != 0)
						throw new IOException("corrupt file (pad pass line)");

					for(int i=1; i < tokens.length; i+= 4)
					{
						Spot s = new Spot();
						try
						{
							s.x = Double.parseDouble(tokens[i+2]);
							s.y = Double.parseDouble(tokens[i+3]);
						}
						catch(NumberFormatException nerr)
						{
							throw new IOException("corrupt file (bad spot position)");
						}
						
						if(!do_filtering)
							spots.add(s);
						else
						{
							XYPair xy = new XYPair();
							xy.x = (int)Math.floor(s.x);
							xy.y = (int)Math.floor(s.y);
							if(filter_hash.isEmpty() || filter_hash.contains(xy))
								spots.add(s);
						}
					}
				}

			}
			if(roi == null)
				throw new IOException("corrupt file (no ROI)");


			if(show_control_panel)
			{
				final String fname = name;
				final Rectangle roi_ = roi;
				final int its = iterations;
				final double ini_FWHM_ = ini_FWHM;
				final double ini_reconstruction_pixel_size_ = ini_reconstruction_pixel_size;
				SwingUtilities.invokeLater(
				new Runnable() {
						public void run() {
								EControlPanel e = new EControlPanel(roi_, pixel_size_in_nm_, fname, null);
								e.append(spots, its);
								e.send_update_canvas_event();
								e.send_status_text_message("Using " + fname + ": " + Integer.toString(its) + " iterations.");
								e.set_reconstructed_pixel_size(ini_reconstruction_pixel_size_);
								e.set_reconstruction_blur_fwhm(ini_FWHM_);
								e.send_update_canvas_event();
						}

				}
				);
			}
			else
			{
				System.out.println("etf\n");
				double zoom=ps / ini_reconstruction_pixel_size; 
				ImagePlus   linear_reconstruction;
				linear_reconstruction = new ImagePlus();
				linear_reconstruction.setProcessor(Reconstruction.reconstruct(roi, zoom, ini_FWHM, ps, spots));
				linear_reconstruction.show();



			}




		}
		catch(IOException another_ioe){
			Toolkit.getDefaultToolkit().beep();
			ij.IJ.showStatus("Error reading data: " + another_ioe.getMessage());
		}
		finally{
			close(in);
		}

	}



}
