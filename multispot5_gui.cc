#include <tag/printf.h>
#undef make_tuple

#include <tr1/tuple>
#include <algorithm>
#include <climits>
#include <iomanip>
#include <map>
#include <cvd/image_io.h>
#include <cvd/image_convert.h>
#include <cvd/glwindow.h>
#include <cvd/morphology.h>
#include <cvd/connected_components.h>
#include <cvd/draw.h>
#include <cvd/gl_helpers.h>
#include <cvd/vector_image_ref.h>
#include <cvd/videodisplay.h>
#include <cvd/byte.h>

#include <gvars3/instances.h>
#include <gvars3/GStringUtil.h>
#include <gvars3/GUI_readline.h>

#include "storm_imagery.h"
#include "multispot5.h"
#include "multispot5_place_choice.h"
#include "utility.h"

using namespace std;
using namespace std::tr1;
using namespace CVD;
using namespace GVars3;
using namespace TooN;

double lim(double x)
{
	return min(max(x, 0.), 1.);
}


Image<byte> scale(const SubImage<double>& i, double ctr, double rng)
{
	Image<byte> s(i.size());
	for(int r=0; r < i.size().y; r++)
		for(int c=0; c < i.size().x; c++)
			Pixel::DefaultConversion<float, byte>::type::convert(lim((i[r][c] - ctr)/rng), s[r][c]);
	return s;
}

void draw_bbox(const BBox& bbox)
{
	glBegin(GL_LINES);
	glVertex(bbox.first);
	glVertex2i(bbox.first.x, bbox.first.y + bbox.second.y);

	glVertex2i(bbox.first.x, bbox.first.y + bbox.second.y);
	glVertex(bbox.first+ bbox.second);

	glVertex(bbox.first+ bbox.second);
	glVertex2i(bbox.first.x + bbox.second.x, bbox.first.y);

	glVertex2i(bbox.first.x + bbox.second.x, bbox.first.y);
	glVertex(bbox.first);

	glEnd();
}

map<string, string> watch;

void watch_var(void*, string comm, string d)
{
	vector<string> vs = ChopAndUnquoteString(d);

	if(vs.size() != 1)
	{
		cerr << "Error: " << comm << " takes 1 argument: " << comm << " gvar\n";
		return;
	}

	watch[vs[0]] = GV3::get_var(vs[0]);
}

bool watch_update()
{
	bool changes=0;

	for(map<string, string>::iterator i=watch.begin(); i != watch.end(); i++)
	{
		string s = GV3::get_var(i->first);

		if(s != watch[i->first])
		{
			changes=1;
			watch[i->first] = s;
		}
	}
	
	return changes;
}

void GUI_Pause(int n=0)
{
	if(!GV3::get<int>("headless", 0, 1))
	{
		glFlush();
		gvar3<int> pause("pause", 1, 1);

		if(n != 0)
			*pause = n;
		(*pause)--;
		while(*pause == 0)
		{
			GUI_Widgets.process_in_crnt_thread();
			usleep(10000);
		}
	}
}

/// Graphics class which draws information to the screen using OpenGL.
class GraphicsGL: public FitSpotsGraphics
{
	private:
		std::auto_ptr<GLWindow>  win;
		GraphicsGL(const GraphicsGL&);
		int debug_window_zoom;

		///Generate circle linesegments as pair of vertices for GL_LINES
		///@param p circle centre
		///@param r circle radius
		void glDrawCircle(const Vector<2>& p, float r)
		{
			float theta=0;
			for(;;)
			{
				glVertex(p + r * makeVector(cos(theta), sin(theta)));
				theta +=0.01;
				if(theta > M_PI*2)
					break;
				glVertex(p + r * makeVector(cos(theta), sin(theta)));
			}
			glVertex(p + r * makeVector(1, 0));
		}
		
		///Set up a GL window so that glDrawPixels and glVertex line up
		///and also so that glDrawPixels is zoomed.
		///@param size Window size
		///@param scale Zoom level
		void set_GL_zoom_size(ImageRef size, double scale)
		{
			double right = size.x;
			double bottom = size.y;
			//double my_scale = scale;


			glMatrixMode(GL_MODELVIEW);
			glLoadIdentity();


			glMatrixMode(GL_PROJECTION);

			glLoadIdentity();
			glOrtho(0-.5, right-.5, bottom-.5, -.5, -1 , 1);    // If the origin is the top left 

			glRasterPos2f(-.5, -.5);

			// video is now the same way as upside down to graphics!
			glPixelZoom(scale, -scale);
		}


	public:
		
		GraphicsGL()
		{
			debug_window_zoom = GV3::get<int>("debug.zoom", 3, 1);
		}
			
		virtual void draw_pixels(const vector<ImageRef>& pix, float r, float g, float b, float a)
		{
			glColor4f(r, g, b, a);
			glPointSize(1.5);
			glBegin(GL_POINTS);
			glVertex(pix);
			glEnd();
		}

		virtual void draw_bbox(const BBox& bbox)
		{
			::draw_bbox(bbox);
		}

	
		virtual void draw_krap(const vector<Vector<4> >& spots, const Image<byte>& im, const BBox& box, int N, Vector<4> s)
		{

			glDrawPixels(im);
			glColor3f(1, 0, 0);
			draw_bbox(box);

			glLineWidth(0.3);
			glBegin(GL_LINES);
			for(unsigned int i=0; i < spots.size(); i++)
			{
				glColor4f(0, 1, 0, .3);

				if((int)i == N && s[0] != 1e99)	
					glDrawCircle(s.slice<2, 2>(), s[1]);
				else
					glDrawCircle(spots[i].slice<2, 2>(), spots[i][1]);
			}
			glEnd();

			glLineWidth(1.0);
			glBegin(GL_LINES);
			for(unsigned int i=0; i < spots.size(); i++)
			{
				glColor3f(1, 0, 0);
				if((int) i == N)
					glColor3f(1, 1, 0);

				if((int)i == N && s[0] != 1e99)	
					glDrawCross(s.slice<2, 2>(), 1);
				else
					glDrawCross(spots[i].slice<2, 2>(), 1);
			}
			glEnd();

			glFlush();
		}

		virtual void glDrawCross(const Vector<2>& p, int size)
		{
			glVertex(p + makeVector(0,  size));
			glVertex(p + makeVector(0, -size));
			glVertex(p + makeVector( size, 0));
			glVertex(p + makeVector(-size, 0));
		}


		virtual void swap()
		{
			win->swap_buffers();
		}


		virtual void init(ImageRef log_ratios_size)
		{
			win = auto_ptr<GLWindow>(new GLWindow(log_ratios_size*debug_window_zoom));
			set_GL_zoom_size(log_ratios_size, debug_window_zoom);
			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			glEnable(GL_LINE_SMOOTH);
		}

		virtual ~GraphicsGL()
		{
		}
};


vector<vector<ImageRef> > get_regions(const SubImage<double>& log_ratios)
{
	gvar3<double> radius("radius", 0, 1);

	//Set the liklihood ratio threshold/spot density prior
	//same thing.
	double threshold = GV3::get<double>("threshold", 0, -1);
	int edge = GV3::get<int>("edge", 0, -1);


	//Threshold image
	Image<byte> thresholded(log_ratios.size(), 0);
	for(int r=0; r < thresholded.size().y; r++)
		for(int c=0; c < min(thresholded.size().x, edge); c++)
			thresholded[r][c] = 255 * (log_ratios[r][c] > threshold);
	
	//Dilate
	Image<byte> dilated = morphology(thresholded, getDisc(*radius), Morphology::BinaryDilate<byte>());

	transform(dilated.begin(), dilated.end(), dilated.begin(), bind1st(multiplies<int>(), 255));
	
	//Connected components of dilated image
	vector<ImageRef> fg;
	for(int r=0; r < thresholded.size().y; r++)
		for(int c=0; c < min(thresholded.size().x, edge); c++)
			if(dilated[r][c])
				fg.push_back(ImageRef(c, r));

	vector<vector<ImageRef> > regions;
	connected_components(fg, regions);

	return regions;
}


void mmain(int argc, char** argv)
{
	GUI.RegisterCommand("watch", watch_var);
	GUI.LoadFile("multispot5.cfg");
	int lastarg = GUI.parseArguments(argc, argv);
	if(lastarg >= argc)
	{	
		cerr << "Specify the images to load\n";
		exit(1);
	}

	
	
	//Load the log_ratios image.
	//We will use this as a starting point for searching for spots.
	Image<double> log_ratios;
	try
	{
		log_ratios = img_load(GV3::get<string>("log_ratios", "", -1));
	}
	catch(LogFileParseError e)
	{
		cerr << "Error loading " << GV3::get<string>("log_ratios", "") << ": " << e.what << endl;
		exit(1);
	}

	Vector<2> log_ratios_zoom = GV3::get<Vector<2> >("log_ratios_zoom", "", -1);

	//Load the raw data, and then load the spot parameters.
	vector<string> files(argv + lastarg, argv + argc);
	vector<Image<float> > ims = load_and_normalize_images(files);

	//How far away from eash spot to look:
	//double spot_sigmas = GV3::get<double>("sigmas", 0., -1);

	
	
	gvar3<int> cluster_to_show("cluster_to_show", 0, -1);
	gvar3<int> use_largest("use_largest", 0, 1);

	vector<vector<ImageRef> > regions;

	GLWindow win(log_ratios.size());
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_LINE_SMOOTH);

	readline_in_current_thread line("> ");

	gvar3<bool> start_processing("process", 0, 1);
	gvar3<int> show_thresholded("show_thresholded", 0, 1);

	
	Image<Rgb<byte> > reg(log_ratios.size(), Rgb<byte>(0,0,0));
	
	bool first = true;
	BBox cbox = make_pair(ImageRef_zero, ImageRef_zero);
	for(;;)
	{
		double centre = GV3::get<double>("centre", 0, -1);
		double range = GV3::get<double>("range", 0, -1);

		line.poll();
		win.make_current();
		GUI_Widgets.process_in_crnt_thread();
		
		if(watch_update() || first)
		{
			first = 0;

			regions = get_regions(log_ratios);
			
			//Colourize regions
			reg.fill(Rgb<byte>(0,0,0));
			for(unsigned int i=0; i < regions.size(); i++)
			{
				Rgb<byte> c;
				do
				{
					c.red = rand()%255;
					c.green = rand()%255;
					c.blue = rand()%255;
				}
				while(min(c.red, min(c.green, c.blue)) < 128 && min(abs(c.red - c.green), min(abs(c.red - c.blue), abs(c.green - c.blue))) < 64);

				for(unsigned int j=0; j < regions[i].size(); j++)
					reg[regions[i][j]] = c;
			}

					
		}
		
		
		if(*use_largest && !regions.empty())
		{
			*cluster_to_show=0;
			for(unsigned int i=1; i < regions.size(); i++)
				if(regions[i].size() > regions[*cluster_to_show].size())
					*cluster_to_show = i;
					
		}
		else
			*cluster_to_show = max(min(*cluster_to_show, (int)regions.size() - 1), 0);

		if(!regions.empty())
			cbox = boundingbox(regions[*cluster_to_show]);

		if(!regions.empty() && *start_processing)
		{
			GraphicsGL graphics;

			place_and_fit_spots(ims, regions[*cluster_to_show], log_ratios, GV3::get<string>("save_spots"), graphics);
			*start_processing=0;
		}

		if(*show_thresholded)
			glDrawPixels(reg);
		else
			glDrawPixels(scale(log_ratios, centre, range));

		if(cbox.second != ImageRef_zero)
			draw_bbox(cbox);

		if(win.has_events())
		{
			vector<GLWindow::Event> e;
			win.get_events(e);

			for(unsigned int i=0; i < e.size(); i++)
			{
				if(e[i].type == GLWindow::Event::RESIZE)
				{
					ImageRef newsize = e[i].size;
					ImageRef imsize = log_ratios.size();
					ImageRef size=imsize;

					float old_r = (float)imsize.x / imsize.y;
					float new_r = (float)newsize.x / newsize.y;


					glViewport(0, 0, newsize.x, newsize.y);
					glMatrixMode(GL_PROJECTION);
					glLoadIdentity();
					double zoom;

					if(new_r > old_r) //Then use the y axis
						zoom = newsize.y / (float)size.y; 
					else
						zoom = newsize.x / (float)size.x;

					glOrtho(-.5/zoom, (newsize.x-1.5)/zoom, (newsize.y-1.5)/zoom, -.5/zoom, -1 , 1);

					glPixelZoom(zoom,-zoom);
					glRasterPos2f(0, 0);
				}
			}
		}
		win.swap_buffers();

		usleep(100000);
	}
}

	
int main(int argc, char** argv)
{
	try{
		mmain(argc, argv);
	}
	catch(LogFileParseError e)
	{
		cerr << "Fatal error: " << e.what << endl;
	}
}
	
