// MarkerSet.cpp: a marker set

#include "MarkerSet.h"
#include <iostream>
using namespace std;

//
// MarkerSet(): default constructor
MarkerSet::MarkerSet()
{
	markers.resize(MARKER_SET_SIZE);
}
// constructor for new implementation for Non Haploid & Variable Window Sizes
MarkerSet::MarkerSet(bool flag)
{
	
}


MarkerSet::MarkerSet(const MarkerSet& copy)
{
	markers = copy.markers; // copy.getMarkerBits();
}

void MarkerSet::clear()
{
	markers.reset();
}


boost::dynamic_bitset<>& MarkerSet::getMarkerBits()
{
	return markers;
}

// getMarker(): gets marker
bool MarkerSet::getMarker(int index) const
{
	return markers.test( index );
}

// have a push back instead of set.... 
void MarkerSet::pushback(int index, bool bit)
{
	markers.push_back(bit);
}

void MarkerSet::set(int index , bool bit )
{ 
	markers.set( index , bit );
}

//update this for variable window size.
void MarkerSet::print(ostream& out, unsigned int start, unsigned int end)
{
	for(unsigned int i=start;i<end && (position_ms*MARKER_SET_SIZE)+i < ALL_SNPS.currentSize();i++)
	{
		out << ALL_SNPS.getSNP(position_ms*MARKER_SET_SIZE + i).getVariant( markers[i] );
	}
}

//update this for variable window
void MarkerSet::print(ostream& out)
{
	print( out , 0 , MARKER_SET_SIZE);
}

// equal(): compares another MarkerSet to see if it has the 
bool MarkerSet::equal(MarkerSet * ms)
{
	return markers == ms->markers;
}


int MarkerSet::diff(MarkerSet * ms) const
{
	return (int) ( markers ^ ms->markers ).count();
}


//update this for variable window
// operator<<:overloaded insertion operator
ostream& operator<<(ostream& fout, const MarkerSet& m1)
{
	for (int i = 0; i < MARKER_SET_SIZE; i++)
		fout << m1.getMarker(i);
	return fout;
}

// end MarkerSet.cpp
