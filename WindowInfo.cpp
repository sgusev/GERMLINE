#include "WindowsList.h"

WindowInfo::WindowInfo(unsigned int n_start, unsigned int n_end, unsigned int n_size)
{
	start = n_start;
	end = n_end;
	size = n_size;
}

unsigned int WindowInfo::getStart()
{
	return start;
}

unsigned int WindowInfo::getEnd()
{
	return end;
}

unsigned int WindowInfo::getSize()
{
	return size;
}

void WindowInfo::setEnd(unsigned int n_end)
{
	end = n_end;
}

void WindowInfo::setSize(unsigned int n_size)
{
	size = n_size;
}


//end WindowInfo.cpp
