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
public class  QuickPalmLoader implements PlugIn {
	
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
		double ini_FWHM=100., ini_reconstruction_pixel_size=10.;
		boolean show_control_panel=true;
		boolean do_filtering=true;

		final double pixel_size_in_nm_ = 200; //Make a fake pixel size because we need one
		try{

			OpenDialog o = new OpenDialog("Open QuickPalm data...", null);
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



			GenericDialog g = new GenericDialog("Reconstruction Parameters");
			g.addNumericField("FWHM (initial) ", ini_FWHM, 0, 5, "nm");
			g.addNumericField("Reconstruction size (initial)", ini_reconstruction_pixel_size, 0, 5, "nm");
			g.addCheckbox("Show control panel", show_control_panel);
			
			((Checkbox)g.getCheckboxes().get(0)).hide();
			
			g.showDialog();

			ini_FWHM = g.getNextNumber();
			ini_reconstruction_pixel_size=g.getNextNumber();
			show_control_panel=g.getNextBoolean();
		

			//Yay @ cargoculting
			InputStreamReader d = new InputStreamReader(in);
			BufferedReader r = new BufferedReader(d);

			String line;
			final ArrayList<Spot> spots = new ArrayList<Spot>();
			int iterations=0;
			Rectangle roi=null;


			HashSet<XYPair> filter_sdddash = new HashSet<XYPair>();
			
			int line_number=0;
			int frames=0;
			double max_val = 0; //Assume square image: almost all scientific imagers are close enough to square
			for(  ; (line = r.readLine()) != null; line_number++)
			{
				//discard the header

				if(line_number == 0)
					continue;
				String tokens[] = line.split("\\p{Space}+");
				
				//Tokens are:
				//number intensity xpx ypx xnm ynm ...

				if((tokens.length) < 15)
					throw new IOException("corrupt file (short line)");

				Spot s = new Spot();
				try
				{
					s.x = Double.parseDouble(tokens[4]) / pixel_size_in_nm_;
					s.y = Double.parseDouble(tokens[5]) / pixel_size_in_nm_;

					//QuickPalm has the frame number as a double, but is always integral
					int frame = (int)Double.parseDouble(tokens[14]);

					if(frame > frames)
						frames = frame;

					if(s.x > max_val)
						max_val = s.x;
					if(s.y > max_val)
						max_val = s.y;
				}
				catch(NumberFormatException nerr)
				{
					throw new IOException("corrupt file (bad spot position)");
				}

				spots.add(s);
			}

			roi = new Rectangle(0, 0, (int)Math.ceil(max_val), (int)Math.ceil(max_val));


			if(show_control_panel)
			{
				final String fname = name;
				final Rectangle roi_ = roi;
				final int its = iterations;
				final double ini_FWHM_ = ini_FWHM;
				final double ini_reconstruction_pixel_size_ = ini_reconstruction_pixel_size;
				final String msg = "Using " + fname + ": " + Integer.toString(frames) + " frames, " + Integer.toString(spots.size()) + " localisations .";
				SwingUtilities.invokeLater(
				new Runnable() {
						public void run() {
								EControlPanel e = new EControlPanel(roi_, pixel_size_in_nm_, fname, null);
								e.append(spots, 100000); // 100000 == frames, to stop "Export" button from being red as a warning. This qualifies as lolhax
								e.send_update_canvas_event();
								e.send_status_text_message(msg);
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
				double zoom=pixel_size_in_nm_ / ini_reconstruction_pixel_size; 
				ImagePlus   linear_reconstruction;
				linear_reconstruction = new ImagePlus();
				linear_reconstruction.setProcessor(Reconstruction.reconstruct(roi, zoom, ini_FWHM, pixel_size_in_nm_, spots));
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
