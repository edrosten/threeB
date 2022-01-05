#include <gvars3/instances.h>
#include <gvars3/GUI_readline.h>
#include <gvars3/GStringUtil.h>


#include <cvd/timer.h>
#include <cvd/image_io.h>
#include <cvd/image_interpolate.h>
#include <cvd/glwindow.h>
#include <cvd/draw.h>
#include <cvd/convolution.h>
#include <cvd/morphology.h>
#include <cvd/vector_image_ref.h>
#include <cvd/gl_helpers.h>
#include <cvd/byte.h>

#include <TooN/TooN.h>

#include <algorithm>
#include <cstdlib>
#include <cerrno>

#include <tag/printf.h>
#undef make_tuple

#include <tr1/tuple>

#include "cxxgplot.h"
#include "debug.h"
#include "storm_imagery.h"

using namespace CVD;
using namespace TooN;
using namespace GVars3;
using namespace std;
using namespace std::tr1;

map<string, map<string, string> > watch;

void next_show(void*, string , string )
{
	GV3::get<int>("show")++;
}
void prev_show(void*, string , string )
{
	GV3::get<int>("show")--;
}

void watch_var(void*, string comm, string d)
{
	vector<string> vs = ChopAndUnquoteString(d);

	string w_class="";

	if(vs.size() != 1 && vs.size() != 2)
	{
		cerr << "Error: " << comm << " takes 1 or 2 arguments: " << comm << " gvar [class]\n";
		return;
	}

	if(vs.size() == 2)
		w_class = vs[1];

	GV3::get<int>(w_class + ".force", 0, 1);

	watch[w_class][vs[0]] = GV3::get_var(vs[0]);
	watch[w_class][w_class + ".force"]="0";
}

vector<vector<Vector<4> > > spots, spots2;
void print_current(void*, string comm, string d)
{
	vector<string> vs = ChopAndUnquoteString(d);

	if(vs.size() != 0)
	{
		cerr << "Error: " << comm << " takes no arguments\n";
		return;
	}

	int n = GV3::get<int>("show", spots.size()-1, 1);

	n = max(0, min(n, (int)spots.size()-1));

	for(unsigned int i=0; i < spots[n].size(); i++)
		cout << spots[n][i] << endl;
}

bool watch_update(const string& w_class="")
{
	bool changes=0;

	//GV3::get<bool>("");

	for(map<string, string>::iterator i=watch[w_class].begin(); i != watch[w_class].end(); i++)
	{
		string s = GV3::get_var(i->first);

		if(s != watch[w_class][i->first])
		{
			changes=1;
			watch[w_class][i->first] = s;
		}
	}

	GV3::get<int>(w_class + ".force", 0, 1) = 0;
	watch[w_class][w_class + ".force"]="0";
	return changes;
}




Image<float> average_image(const vector<Image<float> >& ims)
{
	assert_same_size(ims);
	Image<float> r(ims[0].size(), 0);

	for(unsigned int i=0; i < ims.size(); i++)
		transform(r.begin(), r.end(), ims[i].begin(), r.begin(), plus<float>());

	transform(r.begin(), r.end(), r.begin(), bind2nd(multiplies<float>(), 1./ims.size()));
	return r;
}

Image<byte> scale_to_bytes(const Image<float>& im, const double mult=1)
{
        float lo = *min_element(im.begin(), im.end());
        float hi = *max_element(im.begin(), im.end());
        Image<byte> out(im.size());
        for(int r=0; r < out.size().y-0; r++)
                for(int c=0; c < out.size().x-0; c++)
                        out[r][c] = (int)min(255., floor( mult* ( (im[r][c]-lo)*255/(hi-lo))));

        return out;
}

void glDrawCross(const Vector<2>& p, float size=3.0)
{
	glVertex(p + makeVector(0,  size));
	glVertex(p + makeVector(0, -size));
	glVertex(p + makeVector( size, 0));
	glVertex(p + makeVector(-size, 0));
}

double get_zoom(ImageRef imsize, ImageRef winsize, double scale)
{
	ImageRef newsize = winsize;
	ImageRef size=imsize;

	float old_r = (float)imsize.x / imsize.y;
	float new_r = (float)newsize.x / newsize.y;

	double zoom;

	if(new_r > old_r) //Then use the y axis
		zoom = newsize.y / (float)size.y; 
	else
		zoom = newsize.x / (float)size.x;

	return zoom * pow(10, scale/10);;
}

double set_GL_zoom_size(ImageRef imsize, ImageRef winsize, double scale)
{
	ImageRef newsize = winsize;

	glViewport(0, 0, newsize.x, newsize.y);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	double zoom = get_zoom(imsize, winsize, scale);

	glOrtho(-.5/zoom, (newsize.x-1.5)/zoom, (newsize.y-1.5)/zoom, -.5/zoom, -1 , 1);

	glPixelZoom(zoom,-zoom);
	glRasterPos2f(0, 0);
	return zoom;
}


Rgba<byte> hot2(double d, double gutter)
{
	d = max(0., min(d, 0.99999999));

	if(d < gutter)
		return Rgba<byte>(255, 0, 255, (int)floor(d*255 * (1./gutter)));
	else if(d < 1./3.)
	{
		//d sits between gutter and 1/3
		//Scale to 0--1
		double s = (d - gutter)/(1./3 - gutter);
		return Rgba<byte>(255, 0, (int)floor((1-s)*255), 255);
	}
	else if(d < 2./3.)
		return Rgba<byte>(255, (int)floor((d-1./3.)*255*3), 0, 255);
	else
		return Rgba<byte>(255, 255, (int)floor((d-1./3.)*255*3), 255);
}


Rgba<byte> hot(double d, bool opaque)
{
	d = max(0., min(d, 0.99999999));

	if(d < 1./3.)
		if(opaque)
			return Rgba<byte>((int)floor(d * 255*3), 0, 0, 255);
		else
			return Rgba<byte>(255, 0, 0, (int)floor(d*255*3));
	else if(d < 2./3.)
		return Rgba<byte>(255, (int)floor((d-1./3.)*255*3), 0, 255);
	else
		return Rgba<byte>(255, 255, (int)floor((d-1./3.)*255*3), 255);
}

Rgba<byte> gray(double d)
{
	d = max(0., min(d, 0.99999999));
	int i = floor(d*255);
	return Rgba<byte>(i,i,i,i);
}

Rgba<byte> green(double d)
{
	d = max(0., min(d, 0.99999999));
	int i = floor(d*255);
	return Rgba<byte>(0,i,0,255);
}

Rgba<byte> magenta(double d)
{
	d = max(0., min(d, 0.99999999));
	int r = floor(d*255/1.5);
	int g = floor(d*255/1.3);
	int b = floor(d*255);
	return Rgba<byte>(r,g,b,255);
}

Rgba<byte> lightblue(double d)
{
	d = max(0., min(d, 0.99999999));
	int i = floor(d*255);
	return Rgba<byte>(i,0,i,255);
}

Rgba<byte> red(double d)
{
	d = max(0., min(d, 0.99999999));
	int i = floor(d*255);
	return Rgba<byte>(i,0,0,255);
}

pair<Image<Rgba<byte> >, Image<float> > make_hot_map(ImageRef size, const vector<Vector<4> >& spots, bool use_brightness, bool use_max_intensity, double max_brightness, double sigma, double colorscale, double glow_power, int zoom)
{
	Image<float> density(size*zoom);

	density.fill(0);

	for(unsigned int i=0; i < spots.size(); i++)
	{
		double alpha=1;
		if(use_brightness)
		{	
			double brightness;
			
			if(use_max_intensity)
				brightness = spots[i][0] / (spots[i][1]*sqrt(2*M_PI));
			else
				brightness = spots[i][0] / sqrt(2*M_PI);

			alpha = min(1.0, brightness / max_brightness);
		}
		
		density[ir_rounded((spots[i].slice<2,2>() + Ones*.5)*zoom)] += alpha;
	}
	convolveGaussian(density, zoom * sigma);

	double mul = colorscale/ *max_element(density.begin(), density.end());

	Image<Rgba<byte> > out(density.size());

	for(int r=0; r < out.size().y; r++)
		for(int c=0; c < out.size().x; c++)
			out[r][c] = hot(pow(density[r][c] * mul, glow_power), 0);

	return make_pair(out, density);
}


//Turbo-bodge: ripped shamelessly from make_simulations.cc
vector<Vector<4> > some_simulated_spots(ImageRef size)
{

	vector<Vector<4> > spot_parameters;

	int num_spots = GV3::get<int>("simulations.spots", 0, -1);
	double spot_max_brightness = GV3::get<double>("simulations.max_brightness", 0, -1);
	double spot_fwhm = GV3::get<double>("simulations.fwhm", 0, -1);
	double spot_sigma = spot_fwhm / ( 2 * sqrt(2 * log(2)));
	double spot_brightness = spot_max_brightness * spot_sigma * sqrt(2*M_PI);

	double spot_spacing = GV3::get<double>("simulations.spacing", 0, -1);

	for(int s=0; s < num_spots; s++)
	{
		double xpos = size.x / 2.0 + (s - (num_spots-1)/2.0) * spot_spacing;
		spot_parameters.push_back(makeVector(spot_brightness, spot_sigma, xpos, size.y/2.0));
	}
	
	return spot_parameters;
}

vector<pair<Vector<2>, Vector<2> > > get_normal_points(const vector<Vector<4> >& spots, vector<Vector<2> > linescan_coords, double normal_width)
{
	vector<pair<Vector<2>, Vector<2> > > data;
	if(linescan_coords.size() == 2)
	{

		Vector<2> line = linescan_coords[1] - linescan_coords[0];
		Vector<2> line_vec = line / norm_sq(line);
		Vector<2> norm_vec = unit(Matrix<2>(Data(0, 1, -1, 0)) * line);

		double w = normal_width;

		for(unsigned int i=0; i < spots.size(); i++)
		{
			Vector<2> d = spots[i].slice<2,2>() - linescan_coords[0];

			if(d*line_vec > 0 && d* line_vec < 1 && std::abs(d*norm_vec)  < w)
				data.push_back(make_pair(spots[i].slice<2,2>(),makeVector(d*line_vec, d*norm_vec) ));
		}
	}
	return data;
}


vector<Vector<4> > get_normal_points2(const vector<Vector<4> >& spots, vector<Vector<2> > linescan_coords, double normal_width)
{
	vector<Vector<4> > data;
	if(linescan_coords.size() == 2)
	{

		Vector<2> line = linescan_coords[1] - linescan_coords[0];
		Vector<2> line_vec = line / norm_sq(line);
		Vector<2> norm_vec = unit(Matrix<2>(Data(0, 1, -1, 0)) * line);

		double w = normal_width;

		for(unsigned int i=0; i < spots.size(); i++)
		{
			Vector<2> d = spots[i].slice<2,2>() - linescan_coords[0];

			if(d*line_vec > 0 && d* line_vec < 1 && std::abs(d*norm_vec)  < w)
				data.push_back(spots[i]);
		}
	}
	return data;
}


double scaling(const Vector<4>& spot, int use_brightness, int use_max_intensity, double max_brightness)
{
	double alpha=1;
	if(use_brightness)
	{	
		double brightness;
		
		if(use_max_intensity)
			brightness = spot[0] / (spot[1]*sqrt(2*M_PI));
		else
			brightness = spot[0] / sqrt(2*M_PI);

		alpha = min(1.0, brightness / max_brightness);
	}

	return alpha;
}

void do_linescan(const string& prefix, cplot::Plotter& linescan, const vector<Vector<4> >& spots, vector<Vector<2> > linescan_coords, vector<string>& text)
{
	double max_brightness = GV3::get<double>(prefix + ".max_brightness", 1, 1);
	int use_brightness = GV3::get<int>(prefix + ".use_brightness", 0, 1);
	int use_max_intensity = GV3::get<int>(prefix + ".use_max_intensity", 0, 1);

	double sig = GV3::get<float>(prefix + ".glow_sigma");	
	double xoff  = GV3::get<float>(prefix + ".xoff");	
	double yoff  = GV3::get<float>(prefix + ".yoff");	
	int npoints = GV3::get<int>("linescan.points", 1000, 1);
	double pixel_size_in_nm = GV3::get<double>("pixel_size_in_nm");

	linescan.newline("line title \" " + prefix + " linescan\"");

	Vector<2> off = makeVector(xoff, yoff);

	linescan_coords[0] -= off;
	linescan_coords[1] -= off;

	vector<Vector<4> > n_pts = get_normal_points2(spots, linescan_coords, sig * 4);

	double maxval=-1e99;
	int maxpos = -1;
	vector< Vector<2> > scan;

	for(int n=0; n < npoints; n++)
	{
		double p = n * 1.0 /(npoints-1);

		Vector<2> pos = (linescan_coords[0] * (1-p) + p*linescan_coords[1]);


		double sum=0;

		for(unsigned int i=0; i < n_pts.size(); i++)
			sum += exp(-(norm_sq(pos - n_pts[i].slice<2,2>()) / (2*sig*sig))) * scaling(n_pts[i], use_brightness, use_max_intensity, max_brightness);

		sum /= sqrt(2*M_PI*sig*sig);


		scan.push_back(makeVector(p * norm(linescan_coords[1] - linescan_coords[0])*pixel_size_in_nm, sum));
		//linescan.addpt(scan.back());

		if(scan.back()[1] > maxval)
		{
			maxval = scan.back()[1];
			maxpos = n;
		}

		//if(save)
		//	fo << p * norm(linescan_coords[1] - linescan_coords[0]) << " " <<  sum << endl;
	}

	text.resize(npoints);

	const bool nz = GV3::get<int>(prefix + ".linescan.normalize", 0, 1);
	for(int n=0; n < npoints;  n++)
	{
		Vector<2> p = scan[n];
		if(nz)
			p[1] /= maxval;
			
		ostringstream os;
		if(text[n] == "")
			os << p[0];
		os << " " << scan[n][1]/maxval;
		text[n] += os.str();
		linescan.addpt(p);
	}
}


void dynamic_zoom(const ImageRef win_size, const vector<Vector<4> >& spots, const double scale, const Vector<2> offset, Image<float>& density, Image<Rgba<byte> >& hotmap, string prefix="")
{
	//Oh! The humanity!
	if(prefix != "")
		prefix += ".";

	double max_brightness = GV3::get<double>(prefix + "max_brightness", 1, 1);
	int use_brightness = GV3::get<int>(prefix + "use_brightness", 0, 1);
	int use_max_intensity = GV3::get<int>(prefix + "use_max_intensity", 0, 1);

	//This one is really display related and doesn't need to vary with prefix
	int extra_zoom = GV3::get<int>("dynamic_zoom.zoom", 4, 1);
	density.resize(win_size * extra_zoom);

	float xoff = GV3::get<float>(prefix + "xoff")+.5;
	float yoff = GV3::get<float>(prefix + "yoff")+.5;
	float glow_sigma = GV3::get<float>(prefix + "glow_sigma");
	float glow_power = GV3::get<float>(prefix + "glow_power");
	float glow_transparency = GV3::get<float>(prefix + "glow_transparency", 1, 1);
	float gutter = GV3::get<float>(prefix + "glow.ohot.gutter", .1, 1);
	

	string cmap = GV3::get<string>(prefix + "glow_map", "hot", 1);
	//Find which pixes in the original image are viewable
	//cout << "Top left:\n";
	//cout << -offset << endl;
	//cout << -offset + vec(win.size()) / scale << endl;

	density.fill(0);

	for(unsigned int i=0; i < spots.size(); i++)
	{
		double alpha= scaling(spots[i], use_brightness, use_max_intensity, max_brightness);
		ImageRef pos = ir_rounded((spots[i].slice<2,2>() + makeVector(xoff, yoff) + offset)*scale);

		if(density.in_image(pos))
			density[pos] += alpha;
	}
	
	//This one is really display related and doesn't need to vary with prefix
	double sigma = max(scale*extra_zoom*glow_sigma, GV3::get<double>("dynamic_zoom.min_sigma"));

	convolveGaussian(density, sigma);

	
	double mul;
	if(GV3::get<int>(prefix + "dynamic_zoom.use_fixed_scale"))
		mul = GV3::get<float>(prefix + "dynamic_zoom.fixed_scale");
	else
		mul = GV3::get<float>(prefix + "glow_scale")/ *max_element(density.begin(), density.end());

	if(GV3::get<bool>(prefix + "dynamic_zoom.scale_by_number", 0, 1))
		mul /= spots.size();

	
	hotmap.resize(density.size());
	Image<Rgba<byte> > out=hotmap;
	double gp = glow_power;
	
	if(cmap == "ohot")
		for(int r=0; r < out.size().y; r++)
			for(int c=0; c < out.size().x; c++)
			{
				//out[r][c] = Rgba<byte>(r%255, c%255, ((r+c)%2)*255, 128) ; //hot(pow(density[r][c] * mul, gp), 0);
				out[r][c] = hot2(pow(density[r][c] * mul, gp), gutter);
				out[r][c].alpha *= glow_transparency;
				//out[r][c].alpha = min(out[r][c].alpha*1.0f, glow_transparency*255);
			}
	else if(cmap == "gray" || cmap == "grey")
		for(int r=0; r < out.size().y; r++)
			for(int c=0; c < out.size().x; c++)
			{
				out[r][c] = gray(pow(density[r][c] * mul, gp));
				out[r][c].alpha *= glow_transparency;
			}
	else if(cmap == "green")
		for(int r=0; r < out.size().y; r++)
			for(int c=0; c < out.size().x; c++)
			{
				out[r][c] = green(pow(density[r][c] * mul, gp));
				out[r][c].alpha *= glow_transparency;
			}
	else if(cmap == "red")
		for(int r=0; r < out.size().y; r++)
			for(int c=0; c < out.size().x; c++)
			{
				out[r][c] = red(pow(density[r][c] * mul, gp));
				out[r][c].alpha *= glow_transparency;
			}
	else if(cmap == "magenta")
		for(int r=0; r < out.size().y; r++)
			for(int c=0; c < out.size().x; c++)
			{
				out[r][c] = magenta(pow(density[r][c] * mul, gp));
				out[r][c].alpha *= glow_transparency;
			}
	else if(cmap == "lightblue")
		for(int r=0; r < out.size().y; r++)
			for(int c=0; c < out.size().x; c++)
			{
				out[r][c] = lightblue(pow(density[r][c] * mul, gp));
				out[r][c].alpha *= glow_transparency;
			}
	else
		for(int r=0; r < out.size().y; r++)
			for(int c=0; c < out.size().x; c++)
			{
				//out[r][c] = Rgba<byte>(r%255, c%255, ((r+c)%2)*255, 128) ; //hot(pow(density[r][c] * mul, gp), 0);
				out[r][c] = hot(pow(density[r][c] * mul, gp), 0);
				out[r][c].alpha *= glow_transparency;
				//out[r][c].alpha = min(out[r][c].alpha*1.0f, glow_transparency*255);
			}

	
	//Perform atrificial pixellization
	int px = GV3::get<int>(prefix + "pixelate.size");

	if(px > 1)
		for(int r=0; r < out.size().y-px; r+= px)
			for(int c=0; c < out.size().x-px; c+=px)
			{
				Rgba<int> s(0,0,0,0);
				for(int rr=0; rr < px; rr++)
					for(int cc=0; cc < px; cc++)
						s += out[r+rr][c+cc];
				
				s/=(px*px);
				for(int rr=0; rr < px; rr++)
					for(int cc=0; cc < px; cc++)
						out[r+rr][c+cc] = s;
				
			}
	
	//Now draw a density scale
	if(GV3::get<int>(prefix+"density_scale", 1))
	{
		int width=floor(GV3::get<double>(prefix + "density_scale.width", .7) * win_size.x + .5);
		int height=floor(GV3::get<double>(prefix + "density_scale.height", .7) * win_size.y + .5);
		int xs=floor(GV3::get<double>(prefix + "density_scale.x", .1) * win_size.x + .5);
		int ys=floor(GV3::get<double>(prefix + "density_scale.y", .8) * win_size.y + .5);

		xs = max(0, min(xs, out.size().x-1));
		ys = max(0, min(ys, out.size().y-1));

		width  = max(0, min(width+xs, out.size().x-1)-xs);
		height = max(0, min(height+ys, out.size().y-1)-ys);

		double dmax = *max_element(density.begin(), density.end());

		for(int x=0; x < width; x++)
		{
			float d = dmax * x/(width-1);
			Rgba<byte> col;

			
			if(cmap == "ohot")
				col = hot2(pow(d * mul, gp), gutter);
			else if(cmap == "gray" || cmap == "grey")
				col = gray(pow(d * mul, gp));
			else if(cmap == "green")
				col = green(pow(d * mul, gp));
			else if(cmap == "red")
				col = red(pow(d * mul, gp));
			else if(cmap == "lightblue")
				col = lightblue(pow(d * mul, gp));
			else if(cmap == "magenta")
				col = magenta(pow(d * mul, gp));
			else
				col = hot(pow(d * mul, gp), 0);

			col.red = (col.red*col.alpha)/255;
			col.green = (col.green*col.alpha)/255;
			col.blue = (col.blue*col.alpha)/255;
			col.alpha=255;

			for(int y=0; y < height; y++)
				out[y+ys][x+xs] = col;
		}

		Vector<4> border = GV3::get<Vector<4> >(prefix+"density_scale.border.color", Zeros);
		Rgba<byte> c(border[0], border[1], border[2], border[3]);

		int bw = GV3::get<int>(prefix+"density_scale.border.width");

		bw = min(bw, min(width, height));

		for(int x=0; x < width; x++)
			for(int b=0; b < bw; b++)
			{
				out[ys+b][xs+x] = c;
				out[ys+height-b-1][xs+x] = c;
			}

		for(int y=0; y < height; y++)
			for(int b=0; b < bw; b++)
			{
				out[ys+y][xs+b] = c;
				out[ys+y][xs+width-1-b] = c;
			}

	}

}

Rgb<byte> composite(const Rgba<byte>& fg, const Rgb<byte>& bg)
{
	return Rgb<byte>((bg.red * (255 - fg.alpha) + fg.red * fg.alpha)/255,
	                 (bg.green * (255 - fg.alpha) + fg.green * fg.alpha)/255,
                     (bg.blue * (255 - fg.alpha) + fg.blue * fg.alpha)/255);
}

void dump_video(const ImageRef win_size, const vector<vector<Vector<4> > >& spots, const double scale, const Vector<2> offset, bool transparent, const Image<byte>& mean, double pixel_size, const vector<Image<byte> >& all_images)
{
	int video_number = ++GV3::get<int>("video_number", -1, 0);
	string stub = GV3::get<string>("video_stub");
	string dstub = GV3::get<string>("dump.density_stub");

	Image<float> density;
	Image<Rgba<byte> > hotmap;
	Image<Rgb<byte> > output;

	int scalebar = GV3::get<int>("scalebar", 1, 1);
	double scale_width = GV3::get<float>("scalebar.width", 6, 1);
	Vector<2> scale_pos = GV3::get<Vector<2> >("scalebar.pos");
	double scale_size = GV3::get<double>("scalebar.scale", 100, 1);

	bool also_wf = GV3::get<int>("dump.separate_widefield", 0, 1);
	bool use_mean = GV3::get<int>("use_mean", 1, 1);
	bool save_density = GV3::get<int>("dump.save_density", 0, 1);

	glDisable(GL_BLEND);
	for(unsigned int i=0; i < spots.size(); i++)
	{
		dynamic_zoom(win_size, spots[i], scale, offset, density, hotmap);

		SubImage<Rgb<byte> > hotmap2(0, ImageRef_zero, 0), wf_part(0, ImageRef_zero, 0);
		
		if(also_wf)
		{
			output.resize(hotmap.size().dot_times(ImageRef(2,1)));
			hotmap2 = output.sub_image(ImageRef(hotmap.size().x, 0), hotmap.size());
			wf_part = output.sub_image(ImageRef(0, 0), hotmap.size());
		}
		else
		{
			output.resize(hotmap.size());
			hotmap2 = output;
			wf_part = output;
		}

		unsigned int img = min(max(0, (int)i), (int)(all_images.size()-1));
		
		for(int r=0; r < hotmap.size().y; r++)
			for(int c=0; c < hotmap.size().x; c++)
			{
				//.5 for the pixel offset, then .5 to round properly. Oh yeah
				int xc = (int)floor(-offset[0] + c/scale);
				int yc = (int)floor(-offset[1] + r/scale);

				Rgb<byte> wf(0,0,0), bg(0,0,0);

				if((transparent||also_wf) && mean.in_image(ImageRef(xc, yc)))
				{
					int val;
					if(use_mean)
						val = mean[yc][xc];
					else
						val = all_images[img][yc][xc];

					wf.red   = val;
					wf.green = val;
					wf.blue  = val;
				}

				if(transparent)
					bg = wf;

				//hotmap2[r][c].red = (hotmap[r][c].red * (1+ 2*hotmap[r][c].alpha))/510;
				//hotmap2[r][c].green = (hotmap[r][c].green * (1+ 2*hotmap[r][c].alpha))/510;
				//hotmap2[r][c].blue = (hotmap[r][c].blue * (1+ 2*hotmap[r][c].alpha))/510;

				hotmap2[r][c] = composite(hotmap[r][c], bg);

				if(also_wf)
					wf_part[r][c] = wf;
			}

		//Now draw the scale bar
		if(scalebar)
		{
			int npix = (int)floor(scale_size *scale / pixel_size + .5);
			ImageRef pos = ir_rounded(scale_pos.as_diagonal() * vec(win_size));

			ImageRef tl = pos;
			tl.y -= scale_width/2;

			ImageRef sz = ImageRef(npix, scale_width);

			if(hotmap2.in_image(tl) && hotmap2.in_image(tl + sz))
				hotmap2.sub_image(tl, sz).fill(Rgb<byte>(255,255,255));
			else
				cerr << "Warning: scale bar outside image\n";
		}
		
		

		img_save(output, sPrintf(stub, video_number, i));

		if(save_density)
			img_save(density, sPrintf(dstub, video_number, i));
	}
}

template<class C>
int make_image_list(const BasicImage<C>& mean)
{
	GLuint draw_image = glGenLists(1);

	vector<Vector<2> > linescan_coords;

	glNewList(draw_image, GL_COMPILE);
		glEnable(GL_TEXTURE_RECTANGLE_NV);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
		glTexParameterf( GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_MIN_FILTER, GL_NEAREST );
		glTexParameterf( GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
		glPixelStorei(GL_UNPACK_ALIGNMENT,1);
		glTexImage2D(mean, 0, GL_TEXTURE_RECTANGLE_NV);


		glBegin(GL_QUADS);
		glTexCoord2i(0, 0);
		glVertex2i(0,0);
		glTexCoord2i(mean.size().x, 0);
		glVertex2i(mean.size().x,0);
		glTexCoord2i(mean.size().x,mean.size().y);
		glVertex2i(mean.size().x,mean.size().y);
		glTexCoord2i(0, mean.size().y);
		glVertex2i(0, mean.size().y);
		glEnd ();
		glDisable(GL_TEXTURE_RECTANGLE_NV);
	glEndList();

	return draw_image;
}

Image<Rgba<byte> > mask_to_image(const Image<bool>& mask)
{
	Vector<4> v = GV3::get<Vector<4> >("mask.color");
	Image<Rgba<byte> > ret(mask.size());
	Rgba<byte> c(v[0], v[1], v[2], v[3]);
	
	for(int y=0; y < ret.size().y; y++)
		for(int x=0; x < ret.size().x; x++)
			if(mask[y][x])
				ret[y][x] = c;
			else
				ret[y][x] = Rgba<byte>(0,0,0,0);
	return ret;
}

vector<vector<Vector<4> > > read_spots(string name, const SubImage<bool>& mask, int f = -1)
{
	vector<vector<Vector<4> > > spots;
	string sfn = GV3::get<string>(name, "", f);

	if(sfn == "")
		return spots;

	ifstream sf;
	sf.open(sfn.c_str());

	int err = errno;

	if(!sf.good())
	{
		cerr << "ERROR: failed to open " << sfn << ": " <<strerror(err) << endl;
		exit(1);
	}
	//Read in the spots
	for(int l=0;;l++)
	{
		string line;
		getline(sf, line);
		if(line == "")
			break;


		istringstream ln(line);

		string dummy;
		ln >> dummy;

		if(dummy != "[")
		{
			cerr << "Bad line " << l << ". Expected \"[\" : " << line << endl;
			break;
		}
			
		vector<Vector<4> > vv;
		copy(istream_iterator<Vector<4> >(ln), istream_iterator<Vector<4> >(), back_inserter(vv));

		if(mask.size().area())
		{
			vector<Vector<4> > vv2;
			for(unsigned int i=0; i < vv.size(); i++)
			{
				ImageRef p = ir_rounded(vv[i].slice<2,2>());
				if(mask.in_image(p) && mask[p])
					vv2.push_back(vv[i]);
			}
			
			vv = vv2;
		}

		cout << "Line " << l << " has " << vv.size() << " spots\n";
		spots.push_back(vv);
	}
	return spots;
}


void draw_a_dynamic_density_image(const BasicImage<Rgba<byte> >& hotmap, const Vector<2>& hm, const Vector<2>& offset)
{

	glEnable(GL_BLEND);

	glEnable(GL_TEXTURE_RECTANGLE_NV);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
	glTexParameterf( GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_MIN_FILTER, GL_NEAREST );
	glTexParameterf( GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
	glPixelStorei(GL_UNPACK_ALIGNMENT,1);
	glTexImage2D(hotmap, 0, GL_TEXTURE_RECTANGLE_NV);
	glBegin(GL_QUADS);
	glTexCoord2i(0, 0);
	glVertex(-offset + DiagonalMatrix<2>(makeVector(0,0)) * hm);
	glTexCoord2i(hotmap.size().x, 0);
	glVertex(-offset + DiagonalMatrix<2>(makeVector(1,0)) * hm);
	glTexCoord2i(hotmap.size().x,hotmap.size().y);
	glVertex(-offset + DiagonalMatrix<2>(makeVector(1,1)) * hm);
	glTexCoord2i(0, hotmap.size().y);
	glVertex(-offset + DiagonalMatrix<2>(makeVector(0,1)) * hm);
	glEnd ();
	glDisable(GL_TEXTURE_RECTANGLE_NV);
	glDisable(GL_BLEND);
}


class Macro
{
	public:
		Macro()
		:fn("macro", "", 1)
		{}

		void process()
		{
			//Check ic the file is open
			//If not, then attempt to open it.
			if(!f.is_open())
			{
				if(*fn != "")
				{
					f.open(fn->c_str());

					if(!f.is_open())
						cerr << "Error: could not open file \"" << *fn  << "\" : " << strerror(errno) << endl;
				}
			}
			if(!f.is_open())
				return;

			for(;;)
			{
				string line;
				getline(f, line);

				if(!f.good())
				{
					f.close();
					return;
				}
				else if(line == "BREAK")
					return;
				else if(line == "END")
				{
					f.close();
					*fn="";
					return;
				}
				else
					GUI.ParseLine(line);
			}
		}
	

	private:
		ifstream f;	
		gvar3<string> fn;
};






int main(int argc, char** argv)
{
	cplot::Plotter linescan;

	GUI.RegisterCommand("watch", watch_var);
	GUI.RegisterCommand("next", next_show);
	GUI.RegisterCommand("prev", prev_show);
	GUI.RegisterCommand("dump", print_current);
	GUI.LoadFile("spot_viewer.cfg");
	int lastarg = GUI.parseArguments(argc, argv);
	if(lastarg >= argc)
	{	
		cerr << "Specify the images to load\n";
		exit(1);
	}


	Image<bool> mask;
	if(GV3::get<string>("mask_spots", "", 1) != "")
	{
		mask = img_load(GV3::get<string>("mask_spots", "", 1));
		Image<bool>mask2=morphology(mask, getDisc(GV3::get<float>("mask.dilate", 0, 1)), Morphology::Dilate<bool>());
		mask = mask2;
	}
	
	spots = read_spots("spots", mask);
	spots2 = read_spots("spots2", mask, 1);
	
	//Load the raw data, and then load the spot parameters.
	vector<string> files(argv + lastarg, argv + argc);
	vector<Image<float> > ims = load_and_preprocess_images(files);

	vector<Image<byte> > all_images;
	for(unsigned int i=0; i < ims.size(); i++)
		all_images.push_back(scale_to_bytes(ims[i]));

	Image<byte> mean = scale_to_bytes(average_image(ims));
	Image<float> av = average_image(ims);
	cout << (int)*max_element(av.begin(), av.end()) << endl;

	readline_in_current_thread line("> ");

	GLWindow win(mean.size());
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_LINE_SMOOTH);

	gvar3<double> image_brightness("image_brightness", 1, 1);
	double old_image_brightess = *image_brightness;

	gvar3<float> glow_scale("glow_scale", 1, 1);
	gvar3<float> glow_power("glow_power", 1, 1);
	gvar3<float> glow_sigma("glow_sigma", .1, 1);
	gvar3<int> glow_zoom("glow_zoom", 10, 1);
	gvar3<float> linewidth("linewidth", 1, 1);
	gvar3<float> cross_size("cross_size", 1, 1);
	gvar3<float> xoff("xoff", 0, 1);
	gvar3<float> yoff("yoff", 0, 1);
	gvar3<float> zoom_level("zoom_level", 0, 1);
	gvar3<int> use_hot("hot", 0, 1);
	gvar3<int> hot_transparent("hot_transparent", 0, 1);

	gvar3<int> keep_running("keep_running", 1, 1);
	gvar3<bool> dump_a_video("dump_a_video", 0, 1);

	gvar3<int> show_simulation("show_simulation", 0, 1);

	gvar3<int> dump_density("dump_density", 0, 1);
	gvar3<int> find_fwhm("fwhm.find", 0, 1);
	gvar3<double> fwhm("fwhm", 0, 1);
	gvar3<int> play("play", 0, 1);
	gvar3<float> play_rate("play.rate", .5, 1);
	gvar3<int> save_linescan("save_linescan", 0, 1);
	gvar3<int> save_linescan_num("save_linescan_num", -1, 0);
	gvar3<ImageRef> window_size("window.size", ImageRef_zero, 1);
	gvar3<Vector<2> > gui_offset("offset", Zeros, 1);
	gvar3<double> normal_width("normal_width", 1, 1);
	gvar3<double> pixel_size_in_nm("pixel_size_in_nm", 100., 1);
	gvar3<string> mode("line_mode", "linescan", 1);

	gvar3<int> scalebar("scalebar", 1, 1);
	gvar3<int> show("show", spots.size()-1, 1);
	gvar3<int> show_mask("mask.show", 0, 1);

	bool dragging=0;
	bool dragging_3=0;
	ImageRef drag_start;
	Vector<2> offset=Zeros;

	*gui_offset = offset;

	Vector<2> density_offset = Zeros;
	
	double scale=1;

	Image<Rgba<byte> > hotmap;
	Image<float> density;

	Image<Rgba<byte> > hotmap2;
	Image<float> density2;

	GLuint display_list = glGenLists(1);
	glNewList(display_list, GL_COMPILE);
	glEndList();
	
	GLuint draw_mean = make_image_list(mean);

	int draw_image = -1;

	vector<Vector<2> > linescan_coords;
	vector<Vector<2> > old_linescan_coords;

	bool was_playing=0;
	cvd_timer last_frame_shown_at;

	ImageRef old_size = win.size();
	*window_size = old_size;

	Image<Rgba<byte> > maskim = mask_to_image(mask);
	int mask_list = make_image_list(maskim);

	Macro macro;

	while(*keep_running)
	{
		bool redraw=0;
		bool motion=0;
		bool update_linescan=0;
		*gui_offset = offset;
		line.poll();
		macro.process();
		win.make_current();
		GUI_Widgets.process_in_crnt_thread();
		
		if(*gui_offset != offset)
		{
			motion=1;
			offset = *gui_offset;
		}

		if(*window_size != win.size())
		{
			if(win.size() == old_size) //Window size change has not happened
				win.set_size(*window_size);

			*window_size = win.size();
		}
		old_size = win.size();


		{
			//Treat -1e99 as unset
			gvar3<Vector<2> > linescan1("linescan.1", Ones*-1e99, 1);	
			gvar3<Vector<2> > linescan2("linescan.2", Ones*-1e99, 1);	

			Vector<2> l1 = Ones*-1e99;
			Vector<2> l2 = Ones*-1e99;
			Vector<2> l1old = Ones*-1e99;
			Vector<2> l2old = Ones*-1e99;
			 
			if(linescan_coords.size() > 0)
				l1 = linescan_coords[0];
			if(linescan_coords.size() > 1)
				l2 = linescan_coords[1];

			if(old_linescan_coords.size() > 0)
				l1old = old_linescan_coords[0];
			if(old_linescan_coords.size() > 1)
				l2old = old_linescan_coords[1];
		/*	
			cerr << endl << endl << endl << endl;
			#define E(X) cerr << #X << " = " << X << endl
			E(l1);
			E(l1old);
			E(*linescan1);
			
			E(l2);
			E(l2old);
			E(*linescan2);
*/
			

			if(*linescan1 != l1 || *linescan2 != l2)
			{
//				cerr << "Linescan does not match the gvars.\n";
				//Then the linescan coordinates do not match the GVars.

				if(l1 == l1old && l2==l2old) //The linescan didn't move, so the gvar must have changes
				{
//					cerr << "The GVars have changed\n";
					l1=*linescan1;
					l2=*linescan2;
				}

				//Update the gvar, because the linescan changed (or not)
				*linescan1=l1;
				*linescan2=l2;

				if(l1[0] != -1e99)
				{
					linescan_coords.resize(1);
					linescan_coords[0] = l1;
				}

				if(l2[0] != -1e99)
				{
					linescan_coords.resize(2);
					linescan_coords[1] = l2;
				}

			}

			old_linescan_coords = linescan_coords;
		}


		if(*dump_a_video == 1)
		{
			*dump_a_video = 0;
			dump_video(win.size(), spots, set_GL_zoom_size(mean.size(), win.size(), *zoom_level), offset, *hot_transparent, mean, *pixel_size_in_nm, all_images);
		}


		if(*play)
		{
			if(was_playing==0)
				last_frame_shown_at.reset();

			was_playing=1;

			int prev=*show;

			if(last_frame_shown_at.get_time() > *play_rate)
			{
				last_frame_shown_at.reset();
				*show=max(0, min(*show+1, (int)spots.size()-1));
			}

			if(prev != *show)
				redraw=1;
		}
		else
			was_playing=0;

		if(win.has_events())
		{
			vector<GLWindow::Event> e;
			win.get_events(e);


			for(unsigned int i=0; i < e.size(); i++)
			{
				scale = get_zoom(mean.size(), win.size(), *zoom_level);

				if(e[i].type == GLWindow::Event::RESIZE)
					motion=1;
				else if(e[i].type == GLWindow::Event::EVENT)
					motion=1;
				else if(e[i].type == GLWindow::Event::KEY_DOWN && e[i].which=='n')
				{
					(*show)++;
					*show=min(*show, (int)spots.size()-1);
					motion=1;
				}
				else if(e[i].type == GLWindow::Event::KEY_DOWN && e[i].which=='p')
				{
					(*show)--;
					*show=max(*show, 0);
					motion=1;
				}
				else if(e[i].type == GLWindow::Event::MOUSE_DOWN && (e[i].which & 64))
				{
					(*zoom_level)++;
					motion=1;
				}
				else if(e[i].type == GLWindow::Event::MOUSE_DOWN && (e[i].which & 32))
				{
					motion=1;
					(*zoom_level)--;
				}
				else if(e[i].type == GLWindow::Event::MOUSE_DOWN && (e[i].which & 1))
				{
					dragging=1;
					drag_start = e[i].where;
					motion=1;
				}
				else if(e[i].type == GLWindow::Event::MOUSE_DOWN && (e[i].which & 2))
				{
					cout << vec(e[i].where)/scale - offset - makeVector(*xoff, *yoff) << endl;
				}
				else if(e[i].type == GLWindow::Event::MOUSE_MOVE && dragging==1)
				{
					offset+= vec(e[i].where - drag_start) / scale;
					drag_start = e[i].where;
					motion=1;
				}
				else if(e[i].type == GLWindow::Event::MOUSE_UP && e[i].which == 1 &&dragging==1)
				{	
					dragging=0;
					motion=1;
				}
				else if((e[i].type == GLWindow::Event::MOUSE_DOWN && (e[i].which & 4))|| (e[i].type == GLWindow::Event::MOUSE_MOVE && dragging_3==1 ))
				{
					dragging_3=1;
					Vector<2> pos = vec(e[i].where)/scale - offset - makeVector(*xoff, *yoff);

					if(linescan_coords.size() != 2)
						linescan_coords.push_back(pos);
					else 
					{
						if(norm(pos - linescan_coords[0]) < norm(pos - linescan_coords[1]))
							linescan_coords[0] = pos;
						else
							linescan_coords[1] = pos;

						if(linescan_coords[0][0] > linescan_coords[1][0])
							swap(linescan_coords[0], linescan_coords[1]);
					}

					update_linescan=1;

				}
				else if(e[i].type == GLWindow::Event::MOUSE_UP && e[i].which == 4)
				{
					dragging_3 = 0;
				}
			}
		}

		if(watch_update("mask_update"))
		{
			glDeleteLists(mask_list,1);
			maskim = mask_to_image(mask);
			mask_list = make_image_list(maskim);
		}
		
		if(old_image_brightess != *image_brightness)
		{	
			all_images.clear();
			for(unsigned int i=0; i < ims.size(); i++)
				all_images.push_back(scale_to_bytes(ims[i], *image_brightness));

			mean = scale_to_bytes(average_image(ims), *image_brightness);
			
			old_image_brightess = *image_brightness;
			redraw=1;
			glDeleteLists(draw_mean, 1);
			draw_mean = make_image_list(mean);
		}

		if(watch_update("image_to_show") || draw_image == -1)
		{
			if(draw_image != -1 && draw_image != (int)draw_mean)
				glDeleteLists(draw_image, 1);

			if(GV3::get<int>("use_mean", 1, 1))
				draw_image = draw_mean;
			else
				draw_image = make_image_list(all_images[min(max(0, *show), (int)all_images.size()-1)]);
		}

		if(watch_update("motion"))
			motion=1;
		
		if(watch_update())
		{
			redraw=1;
			motion=1;
		}

		if(redraw)
		{
			int n = GV3::get<int>("show", spots.size()-1, 1);
			double max_brightness = GV3::get<double>("max_brightness", 1, 1);
			double weighting = GV3::get<double>("weighting", .1, 1);
			int use_brightness = GV3::get<int>("use_brightness", 0, 1);
			int use_max_intensity = GV3::get<int>("use_max_intensity", 0, 1);

			n = max(0, min(n, (int)spots.size()-1));

			glDeleteLists(display_list, 1);
			display_list = glGenLists(1);
			glNewList(display_list, GL_COMPILE);

			if(*use_hot)
			{
				tie(hotmap, density) = make_hot_map(mean.size(), spots[n], use_brightness, use_max_intensity, max_brightness, *glow_sigma, *glow_scale, *glow_power, *glow_zoom);

				glEnable(GL_BLEND);

				glEnable(GL_TEXTURE_RECTANGLE_NV);
				glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
				glTexParameterf( GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_MIN_FILTER, GL_NEAREST );
				glTexParameterf( GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
				glPixelStorei(GL_UNPACK_ALIGNMENT,1);
				glTexImage2D(hotmap, 0, GL_TEXTURE_RECTANGLE_NV);
				glBegin(GL_QUADS);
				glTexCoord2i(0, 0);
				glVertex2i(0,0);
				glTexCoord2i(hotmap.size().x, 0);
				glVertex2i(mean.size().x,0);
				glTexCoord2i(hotmap.size().x,hotmap.size().y);
				glVertex2i(mean.size().x,mean.size().y);
				glTexCoord2i(0, hotmap.size().y);
				glVertex2i(0, mean.size().y);
				glEnd ();
				glDisable(GL_TEXTURE_RECTANGLE_NV);
				glDisable(GL_BLEND);
				density_offset = Zeros;

			}
			else
			{

				density_offset =Zeros;
				glEnable(GL_BLEND);
				glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

				glColor(GV3::get<Vector<> >("color", "[1 0 0 1]", 1));
				glLineWidth(*linewidth);
				glBegin(GL_LINES);
				for(unsigned int i=0; i < spots[n].size(); i++)
				{
					if(use_brightness)
					{	
						double brightness;
						
						if(use_max_intensity)
							brightness = spots[n][i][0] / (spots[n][i][1]*sqrt(2*M_PI));
						else
							brightness = spots[n][i][0] / sqrt(2*M_PI);

						double alpha = min(1.0, brightness / max_brightness) * weighting;
						glColor4f(1, 0, 0, alpha);
					}

					glDrawCross(spots[n][i].slice<2,2>() + makeVector(*xoff, *yoff), *cross_size);
				}
				glEnd();

				glDisable(GL_BLEND);

			}

			if(*show_simulation)
			{
				glEnable(GL_BLEND);
				glColor(GV3::get<Vector<> >("simulations.color", "[.3 .3 1 1]", 1));
				glLineWidth(GV3::get<float>("simulations.linewidth", 3, 1));
				float cs = GV3::get<float>("simulations.cross_size", 3, 1);
				glBegin(GL_LINES);
				vector<Vector<4> > sim = some_simulated_spots(mean.size());

				for(unsigned int i=0; i < sim.size(); i++)
					glDrawCross(sim[i].slice<2,2>() + makeVector(*xoff, *yoff), cs);
				glEnd();
				glDisable(GL_BLEND);
			}

			glEndList();
		}

		if(*dump_density)
		{
			img_save(density, "reconstructed_density.tiff");
			*dump_density=0;
		}


		if(motion || update_linescan)
		{
			
			if(GV3::get<int>("dynamic_zoom", 0, 1))
			{
				int n = GV3::get<int>("show", spots.size()-1, 1);
				n = max(0, min(n, (int)spots.size()-1));


				scale=set_GL_zoom_size(mean.size(), win.size(), *zoom_level);
				dynamic_zoom(win.size(), spots[n], scale, offset, density, hotmap, "spots1");

				glTranslate(offset);
				glClear(GL_COLOR_BUFFER_BIT);
				if(*hot_transparent == 1)
					glCallList(draw_image);


				string blend = GV3::get<string>("blend");
				if(blend == "add")
					glBlendFunc(GL_ONE, GL_ONE);
				else
					glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
					

				Vector<2> hm = vec(hotmap.size()) / scale;

				if(!spots2.empty())
				{
					int n2 = GV3::get<int>("show2", spots2.size()-1, 1);
					n2 = max(0, min(n2, (int)spots2.size()-1));
					dynamic_zoom(win.size(), spots2[n2], scale, offset, density2, hotmap2, "spots2");
					draw_a_dynamic_density_image(hotmap2, hm, offset);
				}
				
				draw_a_dynamic_density_image(hotmap, hm, offset);


				glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			}
			else
			{
				scale=set_GL_zoom_size(mean.size(), win.size(), *zoom_level);

				glTranslate(offset);
				glClear(GL_COLOR_BUFFER_BIT);
				if(*use_hot == 0 || *hot_transparent == 1)
					glCallList(draw_image);
				glCallList(display_list);
			}
			
			if(*show_mask)
			{
				glEnable(GL_BLEND);
				glCallList(mask_list);
				glDisable(GL_BLEND);
			}


			if(*scalebar)
			{
				glColor(GV3::get<Vector<> >("scalebar.color", Zeros(3), 1));
				glLineWidth(GV3::get<float>("scalebar.width", 6, 1));
				glEnable(GL_BLEND);
				
				Vector<2> corner = (-offset + DiagonalMatrix<2>(vec(win.size())) * (GV3::get<Vector<2> >("scalebar.pos"))/scale );
				glBegin(GL_LINES);
				glVertex(corner);
				glVertex(corner + makeVector(1, 0) * GV3::get<double>("scalebar.scale", 100, 1)/ *pixel_size_in_nm);
				glEnd();

				glDisable(GL_BLEND);
			}
				
			glTranslatef(*xoff, *yoff,0);

			glEnable(GL_BLEND);
			glLineWidth(GV3::get<double>("linescan.cross.width", 1, 1));
			glColor(GV3::get<Vector<> >("linescan.cross.color", "[0 0 1]", 1));
			glBegin(GL_LINES);
			for(unsigned int i=0; i < linescan_coords.size(); i++)
				glDrawCross(linescan_coords[i], GV3::get<double>("linescan.cross.size", 1, 1)/scale);
			glEnd();
			
	
			if(linescan_coords.size() == 2)
			{
				if(*mode == "linescan" || *mode == "normal points" || *mode == "newlinescan")
				{

					if(*mode == "normal points")
					{
						Matrix<2> r = Data(0, 1, -1, 0);
						Vector<2> normal = unit(r*(linescan_coords[1] - linescan_coords[0])) * *normal_width;

						glColor(GV3::get<Vector<> >("linescan.normalarea.color", "[0 0 1]", 1));
						glBegin(GL_QUADS);


						glVertex(linescan_coords[0]+normal);
						glVertex(linescan_coords[1]+normal);
						glVertex(linescan_coords[1]-normal);
						glVertex(linescan_coords[0]-normal);
						glEnd();

						int n = GV3::get<int>("show", spots.size()-1, 1);
						n = max(0, min(n, (int)spots.size()-1));

						vector<pair<Vector<2>, Vector<2> > > p = get_normal_points(spots[n], linescan_coords, *normal_width);

						glBegin(GL_POINTS);
						glColor3f(0,1,0);
						for(unsigned int i=0; i < p.size(); i++)
							glVertex(p[i].first);
						glEnd();


					}

					glLineWidth(GV3::get<double>("linescan.line.width", 1, 1));
					glColor(GV3::get<Vector<> >("linescan.line.color", "[0 0 1]", 1));
					glBegin(GL_LINES);
					glVertex(linescan_coords[0]);
					glVertex(linescan_coords[1]);
					glEnd();

				}
			}
			glDisable(GL_BLEND);

			win.swap_buffers();
		}

		if((watch_update("linescan") || update_linescan) && linescan_coords.size() == 2 )
		{	
			GV3::get<double>("linescan_length") = norm(linescan_coords[1] - linescan_coords[0])* *pixel_size_in_nm;
			int npoints = GV3::get<int>("linescan.points", 1000, 1);
			image_interpolate<Interpolate::Bicubic, float> di(density);
				
			bool save=0;
			ofstream fo;
			if(*save_linescan)
			{
				save=1;

				*save_linescan=0;
				++*save_linescan_num;

				fo.open(sPrintf("linescan-%02i.txt", *save_linescan_num).c_str());
			}

			if(*mode == "newlinescan")
			{
				int n = GV3::get<int>("show", spots.size()-1, 1);
				n = max(0, min(n, (int)spots.size()-1));
				int n2 = GV3::get<int>("show2", spots2.size()-1, 1);
				n2 = max(0, min(n2, (int)spots2.size()-1));

				vector<string> text;

				do_linescan("spots1", linescan, spots[n], linescan_coords, text);

				if(!spots2.empty())
					do_linescan("spots2", linescan, spots2[n2], linescan_coords, text);

				if(GV3::get<int>("newlinescan.wf", 0, 1))
				{
					linescan.newline("line title \"Widefield\"");
				
					image_interpolate<Interpolate::Bicubic, byte> di(mean);
					vector<Vector<2> >	scan;
					Vector<2> off = makeVector(*xoff, *yoff);
					double hi=-1e99;
					double psnm = *pixel_size_in_nm;
					for(int n=0; n < npoints; n++)
					{
						double p = n * 1.0 /(npoints-1);

						Vector<2> pos = ((linescan_coords[0] * (1-p) + p*linescan_coords[1]));

						if(di.in_image(pos))
						{
							scan.push_back(makeVector(p * norm(linescan_coords[1] - linescan_coords[0])*psnm, di[pos]));
							hi = max(hi, 1.0*di[pos]);
						}
					}

					for(int n=0; n < scan.size();  n++)
					{
						Vector<2> p = scan[n];
						p[1] /= hi;
							
						ostringstream os;
						if(text[n] == "")
							os << p[0];
						os << " " << scan[n][1]/hi;
						text[n] += os.str();
						linescan.addpt(p);
					}

				}

				linescan.draw();


				if(save)
					for(unsigned int i=0; i < text.size(); i++)
						fo << text[i] << endl;
			}
			else if(*mode =="linescan")
			{
				linescan.newline("line title \"linescan\"");
				watch_update("linescan"); 


				vector< Vector<2> > scan;

				double maxval=-1e99;
				int maxpos = -1;
				if(GV3::get<int>("dynamic_zoom") == 1)
				{
					int n = GV3::get<int>("show", spots.size()-1, 1);
					n = max(0, min(n, (int)spots.size()-1));

					vector<pair<Vector<2>, Vector<2> > > n_pts = get_normal_points(spots[n], linescan_coords, *glow_sigma * 4);
					double sig = *glow_sigma;

					for(int n=0; n < npoints; n++)
					{
						double p = n * 1.0 /(npoints-1);

						Vector<2> pos = (linescan_coords[0] * (1-p) + p*linescan_coords[1]);


						double sum=0;

						for(unsigned int i=0; i < n_pts.size(); i++)
							sum += exp(-(norm_sq(pos - n_pts[i].first) / (2*sig*sig)));

						sum /= sqrt(2*M_PI*sig*sig);


						scan.push_back(makeVector(p * norm(linescan_coords[1] - linescan_coords[0])**pixel_size_in_nm, sum));
						linescan.addpt(scan.back());

						if(scan.back()[1] > maxval)
						{
							maxval = scan.back()[1];
							maxpos = n;
						}

						if(save)
							fo << p * norm(linescan_coords[1] - linescan_coords[0]) << " " <<  sum << endl;
					}

					if(*find_fwhm)
					{
					
						int i;
						for(i=maxpos; i >= 0; i--)
							if(scan[i][1] < maxval/2)
								break;
						
						unsigned int j;
						for(j=maxpos; j <scan.size(); j++)
							if(scan[j][1] < maxval/2)
								break;

						if(i >=0 && j < scan.size())
						{
							linescan.newline("points title \"\"");
							linescan.addpt(scan[i]);
							linescan.addpt(scan[j]);
							*fwhm = scan[j][0] - scan[i][0];
						}
						else
						{
							*fwhm = 0;
						}
					}
						
				}
				else
				{
					linescan.newline("line title \"linescan\"");
					double hi=0;
					Vector<2> off = makeVector(*xoff, *yoff);
					for(int n=0; n < npoints; n++)
					{
						double p = n * 1.0 /(npoints-1);

						Vector<2> pos = ( off+ (linescan_coords[0] * (1-p) + p*linescan_coords[1]))**glow_zoom;


						if(di.in_image(pos))
						{
							hi = max(di[pos]*1.0, 1.0*hi);
							linescan.addpt(p * norm(linescan_coords[1] - linescan_coords[0]), 1.0* di[pos]);

							if(save)
								fo << p * norm(linescan_coords[1] - linescan_coords[0]) << " " <<  di[pos] << endl;
						}
					}
					
					if(*show_simulation)
					{
						linescan.newline("line title \"Spot positions\"");

						vector<Vector<4> > sim = some_simulated_spots(mean.size());

						Vector<2> line = unit(linescan_coords[1] - linescan_coords[0]);

						for(unsigned int i=0; i < sim.size(); i++)
						{
							linescan.addpt(line * (sim[i].slice<2,2>() - linescan_coords[0]) + .5, 0.);
							linescan.addpt(line * (sim[i].slice<2,2>() - linescan_coords[0]) + .5, hi);
							linescan.skip();
						}
					}
				}

				linescan.draw();
			}
			else if(*mode == "normal points")
			{
				int n = GV3::get<int>("show", spots.size()-1, 1);
				n = max(0, min(n, (int)spots.size()-1));

				vector<pair<Vector<2>, Vector<2> > > p = get_normal_points(spots[n], linescan_coords, *normal_width);
				linescan.newline("points title \"Spot positions\"");
			

				double mean=0, mean_sq=0;
				for(unsigned int i=0; i < p.size(); i++)
				{
					linescan.addpt(p[i].second[0], p[i].second[1]);
					mean += p[i].second[1];
					mean_sq += p[i].second[1] * p[i].second[1];
				}
				linescan.draw();

				mean/=p.size();
				mean_sq/=p.size();

				GV3::get<double>("normal_std",0,1) = sqrt(mean_sq - mean*mean) ** pixel_size_in_nm;


			}
		}
		
		if(GV3::get<int>("screengrab", 0, 1))
		{
			Image<Rgb<byte> > buf(win.size());
			glReadPixels(buf);
			flipVertical(buf);
			cerr << buf.size() << endl;
			
			try{
				img_save(buf, sPrintf("screenshot-%03i.png", GV3::get<int>("grabno", 0, 1)++));
			}
			catch(Exceptions::All e)
			{
				cerr << "Error: " << e.what << endl;
			}
			GV3::get<int>("screengrab")=0;
		}
		usleep(100000);
	}
}


