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
import java.net.*;
import java.io.*;
import javax.swing.*;
import javax.swing.event.*;
import java.lang.InterruptedException;
import java.lang.System;
import java.lang.Math.*;
import java.lang.reflect.*;

//Please accept my apologies for poor Java code.
//This is my first ever Java program.


class ThreeBGlobalConstants
{
	public static int critical_iterations=200;
};


///Utility class to hold a pair of strings
///
///@ingroup gPlugin
class SPair
{
	public String a, b;
}

///Utility calss to hold a number of handy static functions.
///@ingroup gPlugin
class Util
{
	///Read a file into a string. It simply sidcards errors
	///because if the file is missing from the JAR archive, 
	///then extreme badness has happened and I have no idea how to
	///recover.
	///@param in File to be read
	///@returns file contents in a string
	static String read(Reader in)
	{
		try {
			final char[] buffer = new char[100]; 
			StringBuilder out = new StringBuilder(); 
			int read;
			do 
			{   
				read = in.read(buffer, 0, buffer.length);
				if (read>0)
					out.append(buffer, 0, read);
			} while (read>=0);
			return out.toString();
		}
		catch(java.io.IOException err)
		{
			return "";
			//What do we do here?
		}
	}
	
	///The a file name for saving and complete path using ImageJ's file open dialog.
	///Kep re-querying user if the file will be overwritten. Windows already provides
	///the query built in. 
	///@ingroup gPlugin
	///@param t Initial file name
	///@returns filename and full path
	static SPair getFileName(String t)
	{
		//Get a filename to save as, with appropriate warnings for 
		//overwriting files.
		String fname, fullname;
		while(true)
		{
			SaveDialog save = new SaveDialog("Save 3B output", t, ".txt");
			fname = save.getFileName();

			fullname = save.getDirectory() + File.separator + fname;

			if(fname == null)
				break;

			File test = new File(fullname);
			//Windows' open dialog seems to do overwrite confirmation automatically,
			//so there is no need to do it here.
			if(!ij.IJ.isWindows() && test.exists())
			{
				GenericDialog g = new GenericDialog("Overwrite file?");
				g.addMessage("The file \"" + fname + "\" already exists. Continue and overwrite?");
				g.enableYesNoCancel("Yes", "No");
				g.showDialog();

				if(g.wasOKed())
					break;
				else if(g.wasCanceled())
				{
					fname = null;
					break;
				}
			}
			else
				break;
		}

		SPair r = new SPair();
		r.a = fname;
		r.b = fullname;
		return r;
	}
}





//ImageJ  plugins must have an _ in the name. lolwut?

///ImageJ plugin class
///@ingroup gPlugin
public class three_B implements PlugInFilter {
	
	ImagePlus window;
	ByteProcessor mask;
	String arg;

    public int setup(String arg_, ImagePlus img) {
		window = img;
		arg=arg_;

        return ROI_REQUIRED + SUPPORTS_MASKING + STACK_REQUIRED +NO_CHANGES + DOES_16 + DOES_32 + DOES_8G;
    }

	public void run(ImageProcessor ip) {



		//Load the config file contents
		Reader cfgstream = new InputStreamReader(getClass().getClassLoader().getResourceAsStream("multispot5.cfg"));
		String cfg = Util.read(cfgstream);

		try{
			cfgstream.close();
		}
		catch(IOException close_err){
			Toolkit.getDefaultToolkit().beep();
			ij.IJ.showStatus("Error reading config file.");
			return;
		}

		//Some basic error checking if reading of the config file from
		//the JAR archive fails.
		if(cfg == "")
		{
			Toolkit.getDefaultToolkit().beep();
			ij.IJ.showStatus("Error reading config file.");
			return;
		}


		//The image from getMask() is only the size of the ROI
		//We need it to be congruent with the original image, in order
		//to work with the C++ code.
		mask = new ByteProcessor(ip.getWidth(), ip.getHeight());

		int x = ip.getRoi().x;
		int y = ip.getRoi().y;
		
		//Rectangular selections do not have a mask set, so we get
		//a null pointer exception when we try to copy.
		//So, we have to manually set the pixels, rather than just
		//copy them
		try{
			mask.copyBits(ip.getMask(), x, y, Blitter.COPY);
		}
		catch(NullPointerException e)
		{
			for(int r=0; r < ip.getRoi().height; r++)
				for(int c=0; c < ip.getRoi().width; c++)
					mask.set(c+x, r+y, 255);
		}
		
		//Count the number of set pixels in the mask. This is used to 
		//warn the user if the number is not within a reasonable range.
		int count=0;
		for(int r=0; r < mask.getHeight(); r++)
			for(int c=0; c < mask.getWidth(); c++)
				if(mask.get(c, r) != 0)
					count ++;


		//The non-config file parameters are the range of frames to  operate on
		//and the pixel size in nm (the config works in terms of FWHM in pixels).
		//These have to be sepficied whether the basic or advanced dialog is used.
		int  firstfr;
		int lastfr;
		double pixel_size_in_nm;

		ImageStack s = window.getStack();

		if(arg.equals("advanced"))
		{
			AdvancedDialog ad = new AdvancedDialog(cfg, s.getSize());
			if(!ad.wasOKed())
				return;
			cfg = ad.getTextArea1().getText();
			

			/*pixel_size_in_nm = ad.getPixelSize();
			firstfr = ad.getFirstFrame();
			lastfr = ad.getLastFrame();*/


			pixel_size_in_nm = ad.getNextNumber();
			firstfr = (int)ad.getNextNumber();
			lastfr = (int)ad.getNextNumber();
		}
		else
		{
			ThreeBDialog gd = new ThreeBDialog(count, mask.getWidth()*mask.getHeight(), s.getSize());

			gd.showDialog();
		
			if(!gd.wasOKed())
				return;

			/*final double fwhm = gd.getFWHM();
			pixel_size_in_nm = gd.getPixelSize();
			final int initial_spots = gd.getSpots();
			firstfr = gd.getFirstFrame();
			lastfr = gd.getLastFrame();*/

			//We have to use getNextNumber, otherwise macro recording does not work.
			final double fwhm = gd.getNextNumber();
			pixel_size_in_nm = gd.getNextNumber();
			final int initial_spots = (int)gd.getNextNumber();
			firstfr = (int)gd.getNextNumber();
			lastfr = (int)gd.getNextNumber();


			//Compute the parameters of the log-normal prior such that the mode
			//matches the size of the spots.
			final double sigma = (fwhm / pixel_size_in_nm) / (2*Math.sqrt(2*Math.log(2)));
			final double blur_sigma=0.1;

			//s = exp(mu-sig^2)
			//ln s = mu - sig^2
			//mu = ln s + sig^2

			final double blur_mu = Math.log(sigma) + blur_sigma*blur_sigma;
			
			//Initialized from the current time.
			Random rng = new Random();
			//


			//Append stuff to the config file now. This will be parsed later in C++.
			cfg = cfg + "placement.uniform.num_spots=" + Integer.toString(initial_spots) + "\n"
					  + "blur.mu="    + Double.toString(blur_mu) + "\n" 
					  + "blur.sigma=" + Double.toString(blur_sigma) + "\n"
					  + "seed=" + Integer.toString(rng.nextInt(16777216)) + "\n";

		}
		
		//Acquire a filename to save moderately safely.
		SPair f = Util.getFileName(window.getTitle());
		final String fname = f.a;
		final String fullname = f.b;

		if(fname!= null)
		{
			//Create the 3B runner and the control panel, then execute the control panel in the
			//GUI thread.
			final Rectangle roi = ip.getRoi();
			final double pixel_size_in_nm_ = pixel_size_in_nm;
			final ThreeBRunner tbr = new ThreeBRunner(mask, s, cfg, fullname, firstfr, lastfr);

			SwingUtilities.invokeLater(
					new Runnable() {
							public void run() {
									new EControlPanel(roi, pixel_size_in_nm_, fname, tbr);
							}
					}
			);

		}
	}

};

///Control panel which basically presents the .cfg file in a large edit box
///along with additional necessary things (frame range and pixel size). Also
///make the user acknowledge that they might be getting into trouble.
///This is uglier than I would like.
///@ingroup gPlugin
class AdvancedDialog extends GenericDialog implements DialogListener
{
	boolean ok=true;
	int nframes;
	AdvancedDialog(String cfg, int nframes_)
	{
		super("3B Analysis");
		nframes=nframes_;

		addMessage("Advanced configuration");
		addMessage("                                 ");
		addCheckbox("I understand 3B enough to be editing the text below", false);
		addTextAreas(cfg, null, 40, 100);

		addNumericField("Pixel size (nm / pixel)", 100., 1, 10, "nm");
		addNumericField("First frame", 0., 0, 10, ""); 
		addNumericField("Last frame", nframes-1., 0, 10, ""); 
		addDialogListener(this);

		dialogItemChanged(null, null);

		showDialog();
	}

	public boolean dialogItemChanged(GenericDialog gd, java.awt.AWTEvent e)
	{
		boolean v = ((Checkbox)(getCheckboxes().get(0))).getState();

		if(v != ok)
		{
			ok=v;

			if(ok)
			{
				getTextArea1().setEditable(true);
				((Label)getMessage()).setText("Warning: strange behaviour may result!");
				getPixelSizeField().setEditable(true);
				getFirstFrameField().setEditable(true);
				getLastFrameField().setEditable(true);

			}
			else
			{
				getTextArea1().setEditable(false);
				((Label)getMessage()).setText("                                             ");
				getPixelSizeField().setEditable(false);
				getFirstFrameField().setEditable(false);
				getLastFrameField().setEditable(false);
			}
		}


		//Clamp the frames
		int first = getFirstFrame();
		int last  = getLastFrame();

		int nfirst = Math.max(0, Math.min(nframes-1, first));
		int nlast = Math.max(nfirst, Math.min(nframes-1, last));

		if(first != nfirst)
			getFirstFrameField().setText(Integer.toString(nfirst));
		if(last != nlast)
			getLastFrameField().setText(Integer.toString(nlast));
		return ok;
	}

	int parseInt(String s)
	{
		try
		{
			return Integer.parseInt(s);
		}
		catch(Exception e)
		{
			return 0;
		}
	}
	
	public double parseDouble(String s)
	{
		try
		{
			return Double.parseDouble(s);
		}
		catch(Exception e)
		{
			return 0;
		}
	}

	TextField getPixelSizeField()
	{
		return (TextField)(getNumericFields().get(0));
	}

	TextField getFirstFrameField()
	{
		return (TextField)(getNumericFields().get(1));
	}


	TextField getLastFrameField()
	{
		return (TextField)(getNumericFields().get(2));
	}

	double getPixelSize()
	{
		return parseDouble(getPixelSizeField().getText());
	}
	int getFirstFrame()
	{
		return parseInt(getFirstFrameField().getText());
	}

	int getLastFrame()
	{
		return parseInt(getLastFrameField().getText());
	}
}

///Dialog box for starting 3B
///The dialog highlights bad things in red.
///@ingroup gPlugin
//Not pretty. Should I use dialogItemChanged rather than textValueChanged???
class ThreeBDialog extends GenericDialog
{
	int count_, npix, nframes;
	Color c, bg;

	ThreeBDialog(int count, int npix_, int nframes_)
	{
		super("3B Analysis");
		count_ = count;
		npix   = npix_;
		nframes= nframes_;

		addNumericField("Microscope FWHM", 250.0, 1, 10, "nm");
		addNumericField("Pixel size", 100., 1, 10, "nm");
		addNumericField("Initial number of spots", Math.round(count / 10.), 0, 10, "spots");
		addNumericField("First frame", 0., 0, 10, ""); 
		addNumericField("Last frame", nframes-1., 0, 10, ""); 
		addTextAreas("",null, 8,30);
		getTextArea1().setEditable(false);
		c = getFWHMField().getBackground(); //Get the default background colour
		bg = getBackground(); //Dialog background color
		getTextArea1().removeTextListener(this); //To prevent event thrashing when we write messages

		//Process initial warnings
		textValueChanged(null);
	}
	
	//Try to parse a double and return 0 for an empty string.
	public double parseDouble(String s)
	{
		try
		{
			return Double.parseDouble(s);
		}
		catch(Exception e)
		{
			return 0;
		}
	}
	int parseInt(String s)
	{
		try
		{
			return Integer.parseInt(s);
		}
		catch(Exception e)
		{
			return 0;
		}
	}
	
	TextField getFWHMField()
	{
		return (TextField)(getNumericFields().get(0));
	}
	TextField getPixelSizeField()
	{
		return (TextField)(getNumericFields().get(1));
	}
	TextField getSpotsField()
	{
		return (TextField)(getNumericFields().get(2));
	}
	TextField getFirstFrameField()
	{
		return (TextField)(getNumericFields().get(3));
	}
	TextField getLastFrameField()
	{
		return (TextField)(getNumericFields().get(4));
	}

	double getFWHM()
	{
		return parseDouble(getFWHMField().getText());
	}

	double getPixelSize()
	{
		return parseDouble(getPixelSizeField().getText());
	}

	int getSpots()
	{
		return parseInt(getSpotsField().getText());
	}
	
	int getFirstFrame()
	{
		return parseInt(getFirstFrameField().getText());
	}

	int getLastFrame()
	{
		return parseInt(getLastFrameField().getText());
	}

	int getCount()
	{
		return count_;
	}

	public void textValueChanged(TextEvent e)
	{
		boolean long_run=false;
		String err = "";
		//               012345678901234567890123456789012345678901234567890
		if(getCount() > 1000)
		{
			long_run=true;
			err = err + "Warning: large area selected.\n3B will run very slowly.\n";
		}

		if(npix < 2500)
			err = err + "Warning: image is very small. Fitting may be bad because\nimage noise cannot be accurately estimated.\n";

		if(getSpots() > 500)
		{
			err = err + "Warning: large number of spots.\n3B will run very slowly.\n";
			getSpotsField().setBackground(Color.RED);
		}
		else
			getSpotsField().setBackground(c);


		if(getFWHM() < 200)
		{
			err = err + "Warning: unrealistically small\nmicsoscope resolution.\n";
			getFWHMField().setBackground(Color.RED);
		}
		else if(getFWHM() > 350)
		{
			err = err + "Warning: 3B will not work well with\na poorly focussed microscope.\n";
			getFWHMField().setBackground(Color.RED);
		}
		else
			getFWHMField().setBackground(c);


		if(getPixelSize() < 70)
		{
			getPixelSizeField().setBackground(Color.RED);
			err = err + "Warning: Very small pixels specified.\nAre you sure?\n";
		}
		else if(getPixelSize() > 180)
		{
			getPixelSizeField().setBackground(Color.RED);
			err = err + "Warning: 3B will not work well if the camera\nresolution is too poor.\n";
		}
		else
			getPixelSizeField().setBackground(c);

		//Clamp the frames
		int first = getFirstFrame();
		int last  = getLastFrame();

		int nfirst = Math.max(0, Math.min(nframes-1, first));
		int nlast = Math.max(nfirst, Math.min(nframes-1, last));

		if(first != nfirst)
			getFirstFrameField().setText(Integer.toString(nfirst));
		if(last != nlast)
			getLastFrameField().setText(Integer.toString(nlast));

		if(last -first + 1 > 500)
		{
			getFirstFrameField().setBackground(Color.RED);
			getLastFrameField().setBackground(Color.RED);
			err = err + "Warning: large number of frames specified.\n3B will run very slowly and may be inaccurate.\nFewer than 500 frames is strongly recommended.\n200--300 is generally most suitable.\n";
			long_run = true;
		}
		if(last -first + 1 < 150)
		{
			getFirstFrameField().setBackground(Color.RED);
			getLastFrameField().setBackground(Color.RED);
			err = err + "Warning: small number of frames specified.\n3B may be inaccurate.\nAt least 150 frames is recommended.\n";
		}
		else
		{
			getFirstFrameField().setBackground(c);
			getLastFrameField().setBackground(c);
		}


		if(!long_run && (last -first + 1)*getCount() > 200000)
		{
			err = err + "Warning: large amount of data specified.\n"+ 
			            "3B will run very slowly.\n"+ 
			            "Reduce the area and/or number of frames.\n"+
						"We recommend: \n" + 
						"Number of frames*area in pixels < 200,000.";
			getFirstFrameField().setBackground(Color.RED);
			getLastFrameField().setBackground(Color.RED);
		}


		if(!err.equals(""))
			getTextArea1().setBackground(Color.RED);
		else
			getTextArea1().setBackground(bg);




		getTextArea1().setText(err);
		repaint();

	}
}





///Basic spot class, simply contains coordinates.
///@ingroup gPlugin
class Spot
{
	double x, y;
	Spot()
	{
	}

	Spot(double xx, double yy)
	{
		x=xx;
		y=yy;
	}
}


///Listener class which triggers a complete redraw. Since all redraws are
///pretty much equal, there is no need to distinguish them.
///@ingroup gPlugin
class SomethingChanges implements ChangeListener
{
	private EControlPanel c;
	public SomethingChanges(EControlPanel c_)
	{
		c = c_;
	}
	
	
	public void stateChanged(ChangeEvent e)
	{
		c.send_update_canvas_event();
	}
}

///This class makes a floating point slider bar with an edit box next to 
///it for more precision. Also has a reciprocal option for inverse, reciprocal scaling./
///@ingroup gPlugin
//That's not at all bodged in in an unpleasant way.
class FloatSliderWithBox extends JPanel {

	private JSlider slider;
	private JTextField number;
	private JLabel label;
	private GridBagConstraints completePanelConstraints_;

	private int steps=1000000;
	private double min, max;
	private String text;
	private String units;
	private String format = "%8.3f";

	private double value;
	private boolean reciprocal;

	public FloatSliderWithBox(String text_, double min_, double max_, double value_, int cols, boolean rec_)
	{
		
		super( new GridBagLayout() );

		reciprocal = rec_;
		min=min_;
		max=max_;
		text=text_;
		value = value_;

		if(reciprocal)
		{
			min=1/max_;
			max=1/min_;
		}

		slider = new JSlider(0, steps);

		label = new JLabel();
		
		number = new JTextField(cols);

		//Assemble into a panel
		completePanelConstraints_ = new GridBagConstraints();
		completePanelConstraints_.fill = GridBagConstraints.HORIZONTAL;

		completePanelConstraints_.gridx = 0;
		completePanelConstraints_.gridy = 0;
		completePanelConstraints_.weightx=1;
		this.add( slider, completePanelConstraints_ );

		completePanelConstraints_.gridx = 2;
		completePanelConstraints_.gridy = 0;
		completePanelConstraints_.weightx=0;
		this.add(label, completePanelConstraints_ );

		completePanelConstraints_.gridx = 1;
		completePanelConstraints_.gridy = 0;
		completePanelConstraints_.weightx=0;
		//Add some space to the left to move away from slider slightly
		//And some more to the bottom to help with the way they are displayed
		completePanelConstraints_.insets = new Insets(0,5,10,0);
		this.add( number, completePanelConstraints_ );

		this.setBorder(BorderFactory.createTitledBorder(text));


		slider.addChangeListener(new SliderChanged(this));
		number.addActionListener(new TextChanged(this));
		setValue(value);

	}

	public FloatSliderWithBox setUnits(String s)
	{
		units = s;
		setValue(value);
		return this;
	}

	public FloatSliderWithBox setFormat(String s)
	{
		format = s;
		setValue(value);
		return this;
	}

	public void addChangeListener(ChangeListener changeListener){
		slider.addChangeListener( changeListener );
		return;
	}

	void setValue(double v)
	{
		value = v;
		if(reciprocal)
			slider.setValue((int)Math.round(steps * (1/value-min)/(max-min)));
		else
			slider.setValue((int)Math.round(steps * (value-min)/(max-min)));
			
		number.setText(String.format(format, value));
		label.setText(units);
	}

	double getValue()
	{
		return value;
	}

	public double get_value_from_slider() 
	{ 
		if(reciprocal)
			return 1/((slider.getValue() * 1.0 / steps) * (max - min) + min);
		else
			return (slider.getValue() * 1.0 / steps) * (max - min) + min;
	}
	public double get_value_from_text() 
	{ 
		return Double.parseDouble(number.getText());
	}

	class SliderChanged implements ChangeListener
	{
		FloatSliderWithBox f;
		SliderChanged(FloatSliderWithBox f_)
		{
			f = f_;
		}

		public void stateChanged(ChangeEvent e)
		{
			f.setValue(f.get_value_from_slider());
		}
	}

	class TextChanged implements ActionListener
	{
		FloatSliderWithBox f;
		TextChanged(FloatSliderWithBox f_)
		{
			f = f_;
		}

		public void actionPerformed(ActionEvent e)
		{
			f.setValue(f.get_value_from_text());
		}
	}

}


///Close button issues a window close event. Actual closing logic is then
///done in the close event handler.
///@ingroup gPlugin
class CloseButtonListener implements ActionListener
{
	private JFrame f;
	public CloseButtonListener(JFrame fr)
	{
		f = fr;
	}
	
	
	public void actionPerformed(ActionEvent e)
	{
		//Voodoo from Stack Overflow
		WindowEvent wev = new WindowEvent(f, WindowEvent.WINDOW_CLOSING);
		Toolkit.getDefaultToolkit().getSystemEventQueue().postEvent(wev);
	}
}

///Stop 3B thread
///@ingroup gPlugin
class StopButtonListener implements ActionListener
{
	private EControlPanel t;
	public StopButtonListener(EControlPanel t_)
	{
		t = t_;
	}

	public void actionPerformed(ActionEvent e)
	{
		t.issue_stop();
	}
}

///Listener for the export button.
///@ingroup gPlugin
class ExportButtonListener implements ActionListener
{
	private EControlPanel c;
	public ExportButtonListener(EControlPanel c_)
	{
		c = c_;
	}

	public void actionPerformed(ActionEvent e)
	{
		c.export_reconstruction_as_ij();
	}
}

///Control panel for running 3B plugin and providing interactive update.
///@ingroup gPlugin
class EControlPanel extends JFrame implements WindowListener
{
	private ImagePlus   linear_reconstruction;  //Reconstructed image 
	private ImageCanvas canvas;
	private JButton stopButton, exportButton, closeButton;
	private Color exportButtonColor;
	private JLabel status, time_msg;
	private FloatSliderWithBox blur_fwhm;
	private FloatSliderWithBox reconstructed_pixel_size;

	private Panel complete;
	private JPanel buttons;

	private ThreeBRunner tbr;
	private ArrayList<Spot> pts;


	private Rectangle roi;
	private double zoom=0.01; //Sets initial zoom small,  so the ImageCanvas will always be bigger than its initial size
	                          //otherwise it always puts up the zooming rectangle :(
	private double reconstruction_blur_fwhm=.5;

	private double pixel_size_in_nm;
	private String filename;

	//Number of iterations. Used to colourize export button and provide warnings.
	//Architecture is now getting quite messy.
	private int iterations = 0;

	
	EControlPanel(Rectangle roi_, double ps_, String filename_, ThreeBRunner tbr_)
	{
		//Constract superclass
		super(filename_);

		tbr = tbr_;
		filename=filename_;

		roi = roi_;
		pixel_size_in_nm = ps_;
		pts = new ArrayList<Spot>();

	
		//Now generate the dialog box
		complete = new Panel(new GridBagLayout());

		//Create the image viewer
		linear_reconstruction = new ImagePlus();
		linear_reconstruction.setProcessor(reconstruct());
		canvas = new ImageCanvas(linear_reconstruction);
		
		//Make the image scrollable. This seems to work, if the correct cargo-culting
		//is performed with the gridbaglayout fill constraints. I don't really understand why.
		ScrollPane scroll = new ScrollPane();
		scroll.add(canvas);
		complete.add(scroll, canvas_pos());

		//Create the status message
		status = new JLabel();
		if(tbr != null)
			set_status("Running.");
		else
			set_status("Using loaded data.");
		complete.add(status, status_pos());
		
		//Create the ETA massage
		time_msg = new JLabel();
		if(tbr != null)
			set_time("unknown");
		else
			set_time("not running");
		complete.add(time_msg, time_pos());
		

		
		//The two control sliders
		blur_fwhm = new FloatSliderWithBox("Reconstruction blur FWHM", 0, 150, 100.0, 5, false);
		blur_fwhm.setUnits("nm").setFormat("%5.1f");
		blur_fwhm.addChangeListener(new SomethingChanges(this));
		complete.add(blur_fwhm, blur_pos());

		reconstructed_pixel_size = new FloatSliderWithBox("Reconstructed pixel size", 1, 300, 10.0, 5, true);
		reconstructed_pixel_size.setUnits("nm").setFormat("%5.1f");
		reconstructed_pixel_size.addChangeListener(new SomethingChanges(this));
		complete.add(reconstructed_pixel_size, pixel_size_pos());

		//The button bar
		buttons = new JPanel();
		
		exportButton = new JButton("Export...");
		exportButton.addActionListener(new ExportButtonListener(this));
		buttons.add(exportButton);
		exportButtonColor = exportButton.getBackground();
		
		//No point in having a stop button if this is just a viewer
		if(tbr != null)
		{
			stopButton = new JButton("Stop");
			stopButton.addActionListener(new StopButtonListener(this));
			buttons.add(stopButton);
		}

		closeButton = new JButton("Close");
		closeButton.addActionListener(new CloseButtonListener(this));
		buttons.add(closeButton);


		complete.add(buttons, buttons_pos());

		getContentPane().add(complete);
		
		//This class knows how to handle window close events.
		//Don't close the window by default, so we can try to bring up a 
		//confirmation dialog.
		addWindowListener(this);
		setDefaultCloseOperation(DO_NOTHING_ON_CLOSE);
		

		//Cargo culting...
		pack();
		setVisible(true);
		validate();


		//Start the 3B thread, if we have one to run.
		if(tbr != null)
		{
			tbr.register(this);
			Thread t = new Thread(tbr);
			t.start();
		}
	}

	//Functions for generating constraings for gridbaglayout
	private GridBagConstraints canvas_pos()
	{
		GridBagConstraints g = new GridBagConstraints();
		g.gridx = 0;
		g.gridy = 0;
		g.fill = GridBagConstraints.BOTH;
		g.weightx=100;
		g.weighty=100;
		return g;
	}

	private GridBagConstraints status_pos()
	{
		GridBagConstraints g = new GridBagConstraints();
		g.gridx = 0;
		g.gridy = 1;
		g.anchor = GridBagConstraints.FIRST_LINE_START;
		g.fill = GridBagConstraints.BOTH;
		return g;
	}

	private GridBagConstraints time_pos()
	{
		GridBagConstraints g = new GridBagConstraints();
		g.gridx = 0;
		g.gridy = 2;
		g.anchor = GridBagConstraints.FIRST_LINE_START;
		g.fill = GridBagConstraints.BOTH;
		return g;
	}
	private GridBagConstraints buttons_pos()
	{
		GridBagConstraints g = new GridBagConstraints();
		g.gridx = 0;
		g.gridy = 5;
		g.fill = GridBagConstraints.BOTH;
		return g;
	}

	private GridBagConstraints pixel_size_pos()
	{
		GridBagConstraints g = new GridBagConstraints();
		g.gridx = 0;
		g.gridy = 4;
		g.anchor = GridBagConstraints.FIRST_LINE_START;
		g.fill = GridBagConstraints.BOTH;
		return g;
	}
	
	private GridBagConstraints blur_pos()
	{
		GridBagConstraints g = new GridBagConstraints();
		g.gridx = 0;
		g.gridy = 3;
		g.anchor = GridBagConstraints.FIRST_LINE_START;
		g.fill = GridBagConstraints.BOTH;
		return g;
	}
	
	private void set_status(String s)
	{
		status.setText("Status: " + s);
	}

	private void set_time(String s)
	{
		time_msg.setText("Estimated time: " + s);
	}
	
	///Generic redraw function which recomputes the reconstruction.
	///this is a synchronized method since it makes use of the shared writable
	///array of points.
	private synchronized void update_canvas()
	{
		reconstruction_blur_fwhm = blur_fwhm.getValue();
		zoom = pixel_size_in_nm / reconstructed_pixel_size.getValue();

		linear_reconstruction.setProcessor(reconstruct());
		linear_reconstruction.updateImage();
		linear_reconstruction.updateAndRepaintWindow();
		canvas.repaint();	
		//Update the canvas size and viewport to prevent funny things with zooming.
		canvas.setDrawingSize(linear_reconstruction.getWidth(), linear_reconstruction.getHeight());
		canvas.setSourceRect(new Rectangle(linear_reconstruction.getWidth(), linear_reconstruction.getHeight()));
		//validate();
		//pack();

		if(iterations >= ThreeBGlobalConstants.critical_iterations)
		{
			exportButton.setBackground(exportButtonColor);
		}
		else
		{
			exportButton.setBackground(Color.RED);
		}

	}


	//Window event methods. We have to implement all of these to 
	//be a window listener. We don't care about any of these.
	public void windowOpened(WindowEvent e) { }
	public void windowClosed(WindowEvent e){}
	public void windowDeactivated(WindowEvent e) { }
	public void windowActivated(WindowEvent e) { }
	public void windowDeiconified(WindowEvent e) { }
	public void windowIconified(WindowEvent e) { }

	///Send a stop message to the thread if the window is closed
	public void windowClosing(WindowEvent e)
	{
		if(tbr != null && !tbr.has_stopped())
		{
			int confirmed = JOptionPane.showConfirmDialog(null, "3B still running! Closing will terminate the run. Really close?", "User Confirmation", JOptionPane.YES_NO_OPTION);

			if (confirmed == JOptionPane.YES_OPTION)
			{
				dispose();
				tbr.stop_thread();
			}

		}
		else
		{
			//OK to close
			dispose();
		}
	}

	///Generate a fresh ImageJ image in a standard imageJ window, so that ti can be 
	///manipulated using the usual ImageJ commands. The image is size calibrated, so
	///scalebars are easy to add.
	void export_reconstruction_as_ij()
	{
		ImageProcessor export = linear_reconstruction.getProcessor().duplicate();
		ImagePlus export_win = new ImagePlus(filename + " reconstruction", export);
		export_win.getCalibration().pixelWidth = reconstructed_pixel_size.getValue();
		export_win.getCalibration().pixelHeight = reconstructed_pixel_size.getValue();
		export_win.getCalibration().setXUnit("nm");
		export_win.getCalibration().setYUnit("nm");
		export_win.show();
		export_win.updateAndDraw();
	}
	
	///Compute a brand new reconstruction.
	private synchronized FloatProcessor reconstruct()
	{
		//Reconstruct an image which is based around the ROI in size
		int xoff = (int)roi.x;
		int yoff = (int)roi.y;

		int xsize = (int)Math.ceil(roi.width * zoom);
		int ysize = (int)Math.ceil(roi.height * zoom);
		
		//New blank image set to zero
		FloatProcessor reconstructed = new FloatProcessor(xsize, ysize);

		for(int i=0; i < pts.size(); i++)
		{
			
			//Increment the count 
			int xc = (int)Math.floor((pts.get(i).x-xoff) * zoom + 0.5);
			int yc = (int)Math.floor((pts.get(i).y-yoff) * zoom + 0.5);
			float p = reconstructed.getPixelValue(xc, yc);
			reconstructed.putPixelValue(xc, yc, p+1);
		}

		double blur_sigma = reconstruction_blur_fwhm / (2 * Math.sqrt(2 * Math.log(2))) * zoom / pixel_size_in_nm;

		(new GaussianBlur()).blurGaussian(reconstructed, blur_sigma, blur_sigma, 0.005);
		return reconstructed;
	}


	public void issue_stop()
	{
		if(tbr != null)
			tbr.stop_thread();
		stopButton.setEnabled(false);
	}

	//Below here are callbacks used by the ThreeBRunner thread.

	///Callback issued by the runnre thread which appends the latest set of 
	///points to the array. This is synchronized since it involves the shared
	///writable array of points.
	synchronized public void append(final ArrayList<Spot> p, int its)
	{
		for(int i=0; i < p.size(); i++)
			pts.add(p.get(i));

		iterations = its;
	}
	
	
	///Callback for causing a recomputation of the reconstruction. 
	///This is a safe method for telling the class to recompute and redraw the reconstruction.
	///It can be called from anywhere, but ensures that the redraw happens in the GUI thread.
	public void send_update_canvas_event()
	{
		SwingUtilities.invokeLater(
			new Runnable() {
				final public void run() {
						update_canvas();
					}
				}
		);
	}

	
	///Callback which causes a fatal termination of the 3B control panel. This is generally caused
	///by something changing such that the error condition was not seen in the Java code, but was seen
	///in C++. Ought never to be called.
	public void die(final String s)
	{
		SwingUtilities.invokeLater(
			new Runnable() {
				final public void run() {
						ij.IJ.showStatus("3B run terminated due to an error.\n");
						JOptionPane.showMessageDialog(null, "Error: "+s, "Fatal error", JOptionPane.ERROR_MESSAGE);
						dispose();
					}
				}
		);
	}
	
	///Callback to update the status message safely.
	public void send_status_text_message(final String s)
	{
		SwingUtilities.invokeLater(
			new Runnable(){
				final public void run() {
						set_status(s);
					}
				}
		);
	}


	///Callback to update the status message safely.
	public void send_time_text_message(final String s)
	{
		SwingUtilities.invokeLater(
			new Runnable(){
				final public void run() {
						set_time(s);
					}
				}
		);
	}



};


///This class deals with running the actual 3B code
///It should be run in a separate thread
///It will make calls back to a viewer, and accepts 
///termination calls as well.
///@ingroup gPlugin
class ThreeBRunner implements Runnable
{
	ByteProcessor  mask;
	float[][]      pixels;
	private volatile boolean stop=false;
	private volatile boolean stopped=false;
	private boolean fatal=false;
	EControlPanel  cp;
	String config;
	String filename;
	long start_time;
	int its;

	ThreeBRunner(ByteProcessor mask_, ImageStack s, String cfg_, String filename_, int firstfr, int lastfr)
	{
		config = cfg_;
		filename = filename_;
		
		//Take a copy of the mask to stop it getting trashed by other threads
		mask = (ByteProcessor)mask_.duplicate();
		
		//Convert all the input images into float as a common
		//format. They will be scaled differently from float images 
		//loaded by libcvd, but the later normalization will take care of that
		//
		//Oh, and ImageStack counts from 1, not 0
		pixels = new float[lastfr - firstfr + 1][];
		for(int i=firstfr; i <= lastfr; i++)
			pixels[i-firstfr] = (float[])s.getProcessor(i+1).convertToFloat().getPixels();
	}

	public void register(EControlPanel c)
	{
		cp = c;
	}

	void stop_thread()
	{
		stop=true;
		send_message_string("stopping...");
	}

	boolean has_stopped()
	{
		return stopped;
	}


	//Callbacks for the native code
	void send_message_string(String s)
	{
		cp.send_status_text_message(s);
	}

	void send_new_points(float[] points)
	{
		//New points are sent in the form x, y, x, y
		ArrayList<Spot> tmp = new ArrayList<Spot>();
		for(int i=0; i < points.length; i+=2)
		{
			Spot s = new Spot();
			s.x = points[i];
			s.y = points[i+1];
			tmp.add(s);
		}

		cp.append(tmp, its);
		cp.send_update_canvas_event();
			
		//This happens each "iteration", i.e. each pass.
		long current = (new java.util.Date()).getTime();

		if(its > 0)
		{
			if(its < 8)
				cp.send_time_text_message("computing...");
			else
			{
				long time_per_it = (current -start_time) / its;
				int target = (int)Math.ceil(its * 1.0 / ThreeBGlobalConstants.critical_iterations) * ThreeBGlobalConstants.critical_iterations;
				int its_remaining = target-its;

				System.out.println("Its = " + its);
				System.out.println("rem = " + its_remaining);
				System.out.println("time_per_it = " + time_per_it);

				
				long time_remaining_s = (its_remaining * time_per_it)/1000;
				
				long s = time_remaining_s % 60;
				long m = (time_remaining_s / 60)%60;
				long h = time_remaining_s / 3600;

				cp.send_time_text_message( h + "h" + m + "m (for " + target + " iterations)");
			}
			
		}
		
		//Increment after, since it outputs the initial spots as well
		its++;
	}

	void die(String err)
	{
		cp.die(err);
		fatal=true;
	}

	boolean should_stop()
	{
		return stop;
	}

	native void call(String cfg, float[][] images, byte[] mask, int n_images,int rows, int cols, String file);

	public void run()
	{
		System.out.println("About to call...");
		start_time = (new java.util.Date()).getTime();
		its=0;

		call(config, pixels, (byte[])mask.getPixels(), pixels.length, mask.getHeight(), mask.getWidth(), filename);

		System.out.println("Finished.");
		
		if(!fatal)
		{
			cp.send_status_text_message("Finished\n");
			ij.IJ.showStatus("3B run terminated");
		}

		stopped=true;
	}

	//
	// Decodes % encoding.
	// "a%20space" translates to "a space"
	private static String decodePercent(String str) 
	{
		StringBuffer sb = new StringBuffer();
		for (int i = 0; i < str.length(); i++) 
		{
			char c = str.charAt(i);
			if(c == '+')
				sb.append(' ');
			else if(c == '%')
			{
				sb.append((char) Integer.parseInt(str.substring(i + 1, i + 3), 16));
				i += 2;
			}
			else
				sb.append(c);
		}
		return sb.toString();
	}

	/**
	* Adds the specified path to the java library path
	*
	* @param pathToAdd the path to add
	* @throws Exception
	*/
	public static void addLibraryPath(String pathToAdd) throws Exception{
		final Field usrPathsField = ClassLoader.class.getDeclaredField("usr_paths");
		usrPathsField.setAccessible(true);
	 
		//get array of paths
		final String[] paths = (String[])usrPathsField.get(null);
	 
		//check if the path to add is already present
		for(String path : paths) {
			if(path.equals(pathToAdd)) {
				return;
			}
		}
	 
		//add the new path
		final String[] newPaths = Arrays.copyOf(paths, paths.length + 1);
		newPaths[newPaths.length-1] = pathToAdd;
		usrPathsField.set(null, newPaths);
	}
	
	//OK, so Java is a tricky beast.
	//
	//On Windows, one needs DLLs to be in the PATH
	//On Linux, one needs .so's to be in the LD_LIBRARY_PATH or ld config.
	//
	//System.loadLibrary examines the PATH, then does its own path parsing and
	//examines that too (so it seems). As a result, if you modify Java's .so lookup
	//after it starts, then you can load the required .so, but you haven't modified the 
	//LD_LIBRARY_PATH/PATH, so the loader doesn't know where to look for dependencies.
	//
	//So we have to do everything manually, and load them in the correct order.
	//
	//~Yay.~
	//
	//Addendum:
	//
	//Due to the vast amounts of pain involved in building a DLL linking to other
	//DLLs (specifically the stock LAPACK, BLAS, gFortran) which rely on MingW
	//we now cross compile a DLL with no dependencies at all, which means 
	//only a single DLL needs to be loaded now.
	//
	//(Woo hoo.)


	static String get_plugin_dir()
	{
		try{
			String img_url = (new SPair()).getClass().getResource("img").getFile();
			URI jar_url = null;
			
			jar_url= new URI(img_url.substring(0, img_url.length()-5));

			File jar = new File(jar_url);

			System.out.println("File: " + jar.getCanonicalPath());


			String dir = jar.getParentFile().getCanonicalPath() + File.separator;
			System.out.println("Dir: " + dir);

			return dir;
			
		}
		catch(Exception e){
			return "";
		}
	}
	
	
	static UnsatisfiedLinkError try_to_load_dlls(String [] s)
	{
		String dir = get_plugin_dir();
		
		try{
			for(int i=0; i < s.length; i++)
			{
				System.out.println("Loading " + dir+s[i]);
				System.load(dir + s[i]);
			}
			
			return null;
		}
		catch(UnsatisfiedLinkError e)
		{
			System.out.println("Link error: " + e.getMessage());
			return e;
		}

	}



	static 
	{
		//List of SOs/DLLs to try...
		String[] all_shared_objects = {
			"libthreeB_jni_64.so",
			"libthreeB_jni_32.so",
			"threeB_jni_32.dll",
			"threeB_jni_64.dll",
		}; 

		String errors;

		// From http://blog.cedarsoft.com/tag/java-library-path/
		//
		//Explanation At first the system property is updated with the new
		//value. This might be a relative path – or maybe you want to create
		//that path dynamically. The Classloader has a static field (sys_paths)
		//that contains the paths. If that field is set to null, it is
		//initialized automatically.  Therefore forcing that field to null will
		//result into the reevaluation of the library path as soon as
		//loadLibrary() is called…
		
		UnsatisfiedLinkError err;
		try{
			System.out.println("Trying to load the normal way...");
			System.loadLibrary("threeB_jni");
		}
		catch(UnsatisfiedLinkError e)
		{
			System.out.println("First error: " + e.getMessage());

			errors = e.getMessage();
			err = e;
		
			System.out.println("Trying manual loading...");

			String dir = get_plugin_dir();
			
			boolean success=false;

			for(int i=0; !success && i < all_shared_objects.length; i++)
			{
				System.out.println("Loading " + dir+all_shared_objects[i]);
				try{
					System.load(dir + all_shared_objects[i]);
					success=true;
				}
				catch(UnsatisfiedLinkError e2)
				{
					System.out.println("Link error: " + e2.getMessage());
					errors = errors + "\n" + e2.getMessage();
				}
			}

			if(!success)
			{
				JOptionPane.showMessageDialog(null, "Error loading plugin:\n" + errors, "Error loading plugin", JOptionPane.ERROR_MESSAGE);
				throw err;
			}
		}
	}
}
