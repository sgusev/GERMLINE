// Chromosome.cpp.haplotyped markers for a chromosome

#include "Chromosome.h"
#include <iostream>
using namespace std;


// Chromosome(): default constructor
Chromosome::Chromosome()
{}

MarkerSet * Chromosome::getMarkerSet()
{
	return chromosome[position_ms];
}

MarkerSet * Chromosome::getMarkerSet(unsigned int pos)
{
	return chromosome[pos];
}

void Chromosome::clear()
{
	for ( size_t i = 0 ; i < chromosome.size() ; i++ ) { delete chromosome[i]; }
	chromosome.clear();
	buffer_chromosome.clear();
}

// addMarkerSet(): adds a MarkerSet
void Chromosome::addMarkerSet(MarkerSet * ms)
{
	chromosome.push_back(ms);
}

//overloaded addMarkerSet() : loads the entire marker data into a buffer chromosome
void Chromosome::addMarkers(list<bool>* markers)
{
	buffer_chromosome = *markers;
}

//TODO: update to use WIndowSize instead of MARKER_SET_SIZE
void Chromosome::print_snps(ostream& out, unsigned int start, unsigned int end)
{
	unsigned int p_ms = position_ms;

	unsigned int ms_start = start / MARKER_SET_SIZE;
	unsigned int ms_end = end / MARKER_SET_SIZE;
	if( start % MARKER_SET_SIZE != 0 ) { position_ms = ms_start; chromosome[ms_start++]->print(out,start % MARKER_SET_SIZE,MARKER_SET_SIZE); out << ' '; }
	print(out,ms_start,ms_end);
	if( end % MARKER_SET_SIZE != 0 ) { out << ' '; chromosome[ms_end]->print(out,0,end % MARKER_SET_SIZE); }

	position_ms = p_ms;
}

//TODO: update to use WIndowSize instead of MARKER_SET_SIZE
void Chromosome::print(ostream& out,unsigned int start,unsigned int end)
{
	for(position_ms=start;position_ms<end;position_ms++) 
	{
		if( position_ms > start ) out << ' ';
		chromosome[position_ms]->print(out);
	}
}

//TODO: update to use WIndowSize instead of MARKER_SET_SIZE
ostream& operator<<(ostream &fout, Chromosome& c)
{
	fout << c.getMarkerSet();
	return fout;
}

/////////////////////////////////////////////////////////////////////////////

void Chromosome::updateMarkerSet(unsigned int start, unsigned int end)
{
	MarkerSet* ms = new MarkerSet(true); 
	for(unsigned int i = start;  i<end; i++)
	{
		ms->pushback(buffer_chromosome.front());
		buffer_chromosome.pop_front();				// future implementation for Overlapping windows needs non-modified buffer
	}
	
	chromosome.push_back(ms);
}
void Chromosome::updateMarkerSet(int num_markers)
{
	while(num_markers>0)
	{
		chromosome.back()->pushback(buffer_chromosome.front());
		buffer_chromosome.pop_front();	
		num_markers--;
	}

}

// end Chromosome.cpp

