// Chromosome.cpp.haplotyped markers for a chromosome

#include "Chromosome.h"
#include <iostream>
using namespace std;


// Chromosome(): default constructor
Chromosome::Chromosome()
{mem_chromosome+= sizeof(Chromosome);}

Chromosome::~Chromosome()
{mem_chromosome-= sizeof(Chromosome);}

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
	for ( size_t i = 0 ; i < chromosome.size() ; i++ ) 
	{ 
		mem_markers -=  (sizeof(MarkerSet) + (chromosome[i]->getMarkerBits().num_blocks() * sizeof(unsigned long)));
		delete chromosome[i]; 		
	}
	mem_bufferchr-= ceil((float)buffer_chromosome.size()/8);

	chromosome.clear();
	buffer_chromosome.clear();

}

// addMarkerSet(): adds a MarkerSet
void Chromosome::addMarkerSet(MarkerSet * ms)
{
	chromosome.push_back(ms);
}

//overloaded addMarkerSet() : loads the entire marker data into a buffer chromosome
void Chromosome::addMarkers(vector<bool>* markers)
{
	buffer_chromosome = *markers;
	mem_bufferchr+=  ceil((float)buffer_chromosome.size()/8);
}

//TODO: update to use WIndowSize instead of MARKER_SET_SIZE
void Chromosome::print_snps(ostream& out, unsigned int start, unsigned int end)
{
	int p_ms = position_ms;
	int ms_start = start / MARKER_SET_SIZE;
	int ms_end = end / MARKER_SET_SIZE;
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
		ms->pushback(buffer_chromosome.at(i));
	}
	chromosome.push_back(ms);
	mem_markers+=  (sizeof(MarkerSet) + (ms->getMarkerBits().num_blocks() * sizeof(unsigned long)));
}
void Chromosome::appendMarkerSet(unsigned int end,int num_markers)
{
	mem_markers-=  (sizeof(MarkerSet) + (chromosome.back()->getMarkerBits().num_blocks() * sizeof(unsigned long)));
	while(num_markers>0)
	{
		chromosome.back()->pushback(buffer_chromosome.at(end-num_markers));
		num_markers--;
	}
	mem_markers+=  (sizeof(MarkerSet) + (chromosome.back()->getMarkerBits().num_blocks() * sizeof(unsigned long)));
}

// end Chromosome.cpp

