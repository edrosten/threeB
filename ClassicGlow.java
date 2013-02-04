import ij.plugin.*;
import ij.process.*;
import ij.gui.*;
import ij.*;

import java.awt.image.*;


///Plugin implementing the classic glow/hot
///LUT which seems to be missing from ImageJ
///@ingroup gPlugin
public class ClassicGlow implements PlugIn
{
	static byte red[];
	static byte green[];
	static byte blue[];

	static{
		red = new byte[256];
		green = new byte[256];
		blue = new byte[256];
		
		for(int i=0; i < 256; i++)
		{
			if(i < 85)
				red[i] = (byte)((i*3) & 0xff);
			else
				red[i] = (byte)(255 & 0xff);

			if(i < 85)	
				green[i] = 0;
			else if(i -85 < 85)
				green[i] = (byte)(((i-85)*3) & 0xff);
			else
				green[i] = (byte)(255 & 0xff);

			if(i < 170)
				blue[i] = 0;
			else
				blue[i] = (byte)(((i-170)*3) & 0xff);
		}
	}


	public void run(String arg)
	{
		ImagePlus imp = WindowManager.getCurrentImage();

		//If the image exists, it must be of the correct type, otherwise, make a ramp.
		if (imp!=null)
		{
			if (imp.getType()==ImagePlus.COLOR_RGB)
			{
				IJ.error("Color tables cannot be assiged to RGB Images.");
				return;
			}
		}
		else
		{
			//Create a ramp
			ByteProcessor r = new ByteProcessor(256, 32);

			for(int y=0; y < 32; y++)
				for(int x=0; x < 256; x++)
					r.set(x, y, x);

			imp = new ImagePlus("Glow", r);
			imp.show();
		}


		ImageProcessor ip = imp.getProcessor();
		ColorModel cm = new IndexColorModel(8, 256, red, green, blue);
		ip.setColorModel(cm);
		if (imp.getStackSize()>1) 
			imp.getStack().setColorModel(cm);
		imp.updateAndDraw();
	}
}
