// Stores a hsitory of the various window sizes used for finding matching segments in individual pairs

#ifndef WINDOWSLIST_H
#define WINDOWSLIST_H

#include "BasicDefinitions.h"
#include <vector>
#include <cstdlib>
#include <time.h>
#include <math.h>
using namespace std;



class WindowsList;
class WindowInfo
{
public:
	WindowInfo(unsigned int, unsigned int , unsigned int );

	//void setStart(unsigned int);
	void setEnd(unsigned int);
	void setSize(unsigned int);
	
	unsigned int getStart();
	unsigned int getEnd();
	unsigned int getSize();

private:
	unsigned int start;
	unsigned int end;
	unsigned int size;	
};


class WindowsList
{
public:
	unsigned int getWindowStart(unsigned int);
	unsigned int getWindowStart();
	unsigned int getWindowEnd(unsigned int);
	unsigned int getWindowEnd();
	unsigned int getWindowSize(unsigned int);
	unsigned int getWindowSize();

	
	int err_hom(unsigned int);
	int err_het(unsigned int);

	void clear();
	void initialize(unsigned int);
	void getNewWindowSize(unsigned int );
	int updateWindowSize(unsigned int);

	unsigned long long calculateMem();
private:
	vector<WindowInfo> windows_list;
	vector<WindowInfo>::iterator it ;
	int currentWindowPos;
};


#endif

//end WindowsList.h
