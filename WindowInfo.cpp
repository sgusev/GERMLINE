#include "WindowInfo.h"
#include "iostream"
#include "list"

using namespace std;

WindowInfo::WindowInfo(unsigned int start1, unsigned int end1, unsigned int size1)
{
	start = start1;
	end =end1;
	size = size1;
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

//returns a random number between 100 and 200 
//void WindowInfo::generateWindowSizes()
//{
//	//return ((unsigned int) (rand() % 100 +100));
//	cout<<"\n\t\tGenerate Window Sizes ";
//	unsigned int pos =0, r=0;
//	 while (pos < All_SNPS_CURRENT_SIZE)
//	{
//		if ((All_SNPS_CURRENT_SIZE - pos ) <= 1000) 	
//			{cout<<"\nSET SIZE = "<< All_SNPS_CURRENT_SIZE-pos;  pos+= (All_SNPS_CURRENT_SIZE-pos);  }
//		else 
//			{ r= (unsigned int) (rand() % 100 +1000); cout<<"\nSET SIZE = "<< r; pos+= r; }
//		num_sets++; cout <<"\tEND position = "<<pos ;
//	}
//
//}
