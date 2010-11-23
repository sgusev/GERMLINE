// WindowsList.cpp: History of various window sizes used for matching 

#include "WindowsList.h"
using namespace std;


void WindowsList::clear()
{
	windows_list.clear();
}

void WindowsList::initialize(unsigned int snps_size)
{
	windows_list.reserve((snps_size/MIN_WINDOW_SIZE)+1);
	currentWindowPos=-1;
}

void WindowsList::getNewWindowSize(unsigned int snps_size)
{
	unsigned int r;
	if (currentWindowPos == -1)
	{	
		r=MIN_WINDOW_SIZE;					
		windows_list.push_back(new WindowInfo(0,r,r));
	}
	else
	{
		if ((snps_size - getWindowEnd() ) <= MIN_WINDOW_SIZE)  	
		{			
			windows_list.push_back(new WindowInfo(getWindowEnd(),snps_size ,snps_size -getWindowEnd()));
		}
		else 
		{
			r=MIN_WINDOW_SIZE;				
			windows_list.push_back(new WindowInfo(getWindowEnd() ,getWindowEnd() +r,r));
		}
	}
	currentWindowPos++;
	num_sets++;
}

int WindowsList::updateWindowSize(unsigned int snps_size)
{
	if(getWindowEnd() != snps_size)
	{ 
		windows_list[currentWindowPos]->setEnd(getWindowEnd()+ WINDOW_FACTOR);
		windows_list[currentWindowPos]->setSize(getWindowSize()+ WINDOW_FACTOR);
		return WINDOW_FACTOR;
	}
	return 0;
}

unsigned int WindowsList::getWindowStart(unsigned int ms)
{	return ((WindowInfo*)windows_list[ms])->getStart();}


unsigned int WindowsList::getWindowEnd(unsigned int ms)
{	return ((WindowInfo*)windows_list[ms])->getEnd(); }

unsigned int WindowsList::getWindowSize(unsigned int ms)
{	return ((WindowInfo*)windows_list[ms])->getSize();}


unsigned int WindowsList::getWindowStart()
{	return ((WindowInfo*)windows_list[currentWindowPos])->getStart();}

unsigned int WindowsList::getWindowEnd()
{	return ((WindowInfo*)windows_list[currentWindowPos])->getEnd();}

unsigned int WindowsList::getWindowSize()
{	return ((WindowInfo*)windows_list[currentWindowPos])->getSize();}


int WindowsList::err_hom(unsigned int ms)
{	return (int) floor(getWindowSize(ms)*MAX_ERR_HOMp/100);}

int WindowsList::err_het(unsigned int ms)
{	return (int) floor(getWindowSize(ms)*MAX_ERR_HETp/100);}

//end WindowsList.cpp