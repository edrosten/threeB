#include <sstream>
#include <algorithm>
#include <cvd/image.h>
#include <cvd/image_convert.h>

#include "ThreeBRunner.h"
#include "storm_imagery.h"
#include "multispot5.h"
#include "multispot5_place_choice.h"
#include "utility.h"
#include <gvars3/instances.h>

#include <tag/printf.h>
#undef make_tuple
#include <tr1/tuple>

#ifdef DEBUG
	#include <cvd/image_io.h>
	#define D(X) do{X;}while(0)
#else
	#define D(X) do{}while(0)
#endif


using namespace std;
using namespace std::tr1;
using namespace CVD;
using namespace GVars3;
using namespace tag;
using namespace TooN;


/**
3B User interface for the Java plugin.

This particular UI ferries various messages between Java and C++
via callbacks and polling.

@ingroup gPlugin
*/
class JNIUserInterface: public UserInterfaceCallback
{
	private:
		JNIEnv *env;
		jobject ThreeBRunner_this;
		jmethodID send_message_string;
		jmethodID die;
		jmethodID should_stop;
		jmethodID send_new_points;
		int passes;

	public:
		JNIUserInterface(JNIEnv* env_, jobject jthis)
		:env(env_),ThreeBRunner_this(jthis)
		{
			jclass cls = env->GetObjectClass(jthis);
			
			send_message_string  = env->GetMethodID(cls, "send_message_string", "(Ljava/lang/String;)V");
			die  = env->GetMethodID(cls, "die", "(Ljava/lang/String;)V");

			should_stop  = env->GetMethodID(cls, "should_stop", "()Z");

			send_new_points = env->GetMethodID(cls, "send_new_points", "([F)V");

			passes = GV3::get<int>("main.passes");
		}


		virtual void per_spot(int iteration, int pass, int spot_num, int total_spots)
		{
			send_message(sPrintf("Iteration %i, optimizing  %4i%%", iteration*passes+pass, 100 *spot_num / total_spots));
		}

		virtual void per_modification(int iteration, int spot_num, int total_spots)
		{
			send_message(sPrintf("Iteration %i, modifying  %4i%%", iteration*passes+passes-1, 100 *spot_num / total_spots));
		}

		virtual void per_pass(int , int , const std::vector<TooN::Vector<4> >& spots)
		{
			//Copy data into the correct format
			vector<jfloat> pts_data;
			for(unsigned int i=0; i < spots.size(); i++)
			{
				pts_data.push_back(spots[i][2]);
				pts_data.push_back(spots[i][3]);
			}
			
			//Allocate a java array and copy data into it
			jfloatArray pts = env->NewFloatArray(pts_data.size());
			env->SetFloatArrayRegion(pts, 0, pts_data.size(), pts_data.data());
			
			//Make the call...
			jvalue pts_obj;
			pts_obj.l = pts;

			env->CallVoidMethod(ThreeBRunner_this, send_new_points, pts_obj);
			
			//Free the object
			env->DeleteLocalRef(pts);
		}

		virtual void perhaps_stop()
		{
			bool stop = env->CallBooleanMethod(ThreeBRunner_this, should_stop);
			if(stop)
				throw UserIssuedStop();
		}


		void send_message(const string& s)
		{
			jvalue message_string;
			message_string.l = env->NewStringUTF(s.c_str());
			env->CallVoidMethod(ThreeBRunner_this, send_message_string, message_string);
			env->DeleteLocalRef(message_string.l);
		}

		void fatal(const string& s)
		{
			jvalue message_string;
			message_string.l = env->NewStringUTF(s.c_str());
			env->CallVoidMethod(ThreeBRunner_this, die, message_string);
			env->DeleteLocalRef(message_string.l);
		}
};

///Get a local C++ copy of a java string.
///@ingroup gPlugin
string get_string(JNIEnv *env, jstring js)
{
	const char* str;

	//Covert the config into a string
	str = env->GetStringUTFChars(js, NULL);

	string stdstring(str);
	env->ReleaseStringUTFChars(js, str);


	return stdstring;
}

///Get a local C++ copy of an image from a jbyteArray coming from the guts of ImageJ
///@ingroup gPlugin
Image<jbyte> get_local_copy_of_image(JNIEnv* env, jbyteArray data, int rows, int cols)
{
	//This takes a copy of the pixels (perhaps)
	jbyte* pix = env->GetByteArrayElements(data, NULL);

	BasicImage<jbyte> pix_im(pix, ImageRef(cols, rows));

	Image<jbyte> im;
	im.copy_from(pix_im);
	
	//This frees the pixels if copied, or releases a reference
	env->ReleaseByteArrayElements(data,pix, JNI_ABORT);
	
	return im;
}

///Get a local C++ copy of an image from a jfloatArray coming from the guts of ImageJ
///@ingroup gPlugin
Image<float> get_local_copy_of_image(JNIEnv* env, jfloatArray data, int rows, int cols)
{
	//This takes a copy of the pixels (perhaps)
	float* pix = env->GetFloatArrayElements(data, NULL);

	BasicImage<float> pix_im(pix, ImageRef(cols, rows));

	Image<float> im;
	im.copy_from(pix_im);
	
	//This frees the pixels if copied, or releases a reference
	env->ReleaseFloatArrayElements(data,pix, JNI_ABORT);
	
	return im;
}


///Run the 3B code.
///@ingroup gPlugin
JNIEXPORT void JNICALL Java_ThreeBRunner_call
  (JNIEnv *env, jobject jthis, jstring cfg, jobjectArray images, jbyteArray mask_data, jint n_images, jint rows, jint cols, jstring file)
{
	istringstream config(get_string(env, cfg));
	GUI.ParseStream(config);

	JNIUserInterface ui(env, jthis);
	ui.send_message("Initializing...");
	
	string filename = get_string(env, file);

	//Attmpt to open the file
	ofstream save_spots;
	save_spots.open(filename.c_str());
	int err = errno;

	if(!save_spots.good())
	{
		ui.fatal("failed to open " + filename + ": " + strerror(err));
		return;
	}

	vector<ImageRef> maskir;
	Image<double> maskd;
	{
		Image<jbyte> mask = get_local_copy_of_image(env, mask_data, rows, cols);

		
		#ifdef DEBUG
			img_save(mask, "/tmp/plugin_debug_mask.png");
		#endif


		maskd = convert_image(mask);
		for(ImageRef p(-1, 0); p.next(mask.size()); )
			if(mask[p])
				maskir.push_back(p);
	}	


	vector<Image<float> > ims;

	for(int i=0; i < n_images; i++)
	{
		jfloatArray f = static_cast<jfloatArray>(env->GetObjectArrayElement(images, i));
		ims.push_back(preprocess_image(get_local_copy_of_image(env, f, rows, cols)));
		env->DeleteLocalRef(f);
		
		#ifdef DEBUG
			img_save(ims.back(), sPrintf("/tmp/plugin_debug_image-%05i.tiff", i));
		#endif

	}

	double mean, variance;
	tie(mean, variance) = mean_and_variance(ims);

	for(unsigned int i=0; i < ims.size(); i++)
			transform(ims[i].begin(), ims[i].end(), ims[i].begin(), bind1st(multiplies<double>(), 1/ sqrt(variance)));
	 
	tie(mean, variance) = mean_and_variance(ims);

	//A sanity check.
	cerr << "Rescaled:\n";
	cerr << "mean = " << mean << endl;
	cerr << "std  = " << sqrt(variance) << endl;
	cerr << "Version 1.1" << endl;


	auto_ptr<FitSpotsGraphics> gr = null_graphics();

	place_and_fit_spots(ims, maskir, maskd, save_spots, *gr, ui);
}

