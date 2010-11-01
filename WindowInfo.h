// Stores a hsitory of the various window sizes used for finding matching segments in individual pairs

#ifndef WINDOWINFO_H
#define WINDOWINFO_H

#include "BasicDefinitions.h"

class WindowInfo
{
public:
	WindowInfo(unsigned int start1, unsigned int end1, unsigned int size1);
	~WindowInfo();

	//void setStart(unsigned int);
	//void setEnd(unsigned int);
	//void setSize(unsigned int);
	
	unsigned int getStart();
	unsigned int getEnd();
	unsigned int getSize();

	//static void generateWindowSizes(list<WindowInfo*> w);

private:
	unsigned int start;
	unsigned int end;
	unsigned int size;

	
};

#endif
