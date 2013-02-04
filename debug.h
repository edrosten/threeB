#ifndef STORM_INCLUDE_DEBUG_H
#define STORM_INCLUDE_DEBUG_H

#include <cassert>

/** Determines that all images in the incoming container are the 
same size, and that the container is not empty
@param images Container to check
@ingroup gDebug
*/
template<class C> void assert_same_size(const C& images)
{
	assert(!images.empty());
	for(typename C::const_iterator i=images.begin(); i != images.end(); i++)
		assert(i->size() == images.front().size());
}


#endif

