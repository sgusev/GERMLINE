// KOSIndividualsExtractor.cpp: Manages input from Plink files

#include "PEDIndividualsExtractor.h"
#include <cctype>
#include <iostream>
using namespace std;

// PEDIndividualsExtractor(): default constructor
PEDIndividualsExtractor::PEDIndividualsExtractor()
{}

void PEDIndividualsExtractor::stripWhitespace()
{
	if (stream.is_open())
	{
		char c;
		while((c=stream.peek())!=EOF && isspace(c))
		stream.get();
	}
	else
	{
		cerr << "WARNING:PolymorphicIndividualsExtractor::stripWhiteSpace():stream is not open" << endl;
		valid_flag = false;
	}
}

void PEDIndividualsExtractor::loadInput()
{
	individualsP = &ALL_SAMPLES;
	ALL_SNPS.processMAPFile();
	ALL_SNPS.beginChromosome();
	numberOfMarkers = ALL_SNPS.size();	
	while (!stream.eof() )
	{
		loadIndividuals();  
		stream.seekg(4*(numberOfMarkers-ALL_SNPS.getROIEnd().getMarkerNumber() -1)+1,ios::cur);
	}
	stream.clear();
}

// getInput(): gets individuals from .ped file
void PEDIndividualsExtractor::getInput(string map, string ped)
{
	if (map == "")
	{
		cout << "Please enter the MAP file name" << endl;
		cin >> map_file;
	}
	else map_file = map;
	if (ped == "")
	{
		cout << "Please enter the PED file name" << endl;
		cin >> ped_file;
	}
	else ped_file = ped;
	if ( !ALL_SNPS.setFile( map_file ) )
	{
		cerr << "WARNING:PEDIndividualsExtractor::getInput():cannot open map file" << endl;
		valid_flag = false;
		return;
	}
	stream.open( ped_file.c_str() );
	if ( !stream )
	{
		cerr << "WARNING:PEDIndividualsExtractor::getInput():cannot open ped file" << endl;
		valid_flag = false;
		return;
	}
}

// getIndividuals(): gets the next nuclear family from stream
void PEDIndividualsExtractor::loadIndividuals()
{
	string discard, ID, famID;
	stream >> famID >> ID >> discard >> discard >> discard >> discard;
	if(!stream.good()) return;
	if ( HAPLOID )
	{
		Individual * new_ind[2];
		new_ind[0] = new Individual();
		new_ind[1] = new Individual();
		new_ind[0]->setOffset( stream.tellg() );
		new_ind[1]->setOffset( stream.tellg() );
		new_ind[0]->setID(famID + " " + ID + ".0" );
		new_ind[1]->setID(famID + " " + ID + ".1" );
		//////////////////////////////////////////////////
		loadCompleteMarkerSet(new_ind);
		//////////////////////////////////////////////////
		individualsP->addIndividual( new_ind[0] );
		individualsP->addIndividual( new_ind[1] );
	} else
	{
		Individual* new_ind = new Individual();
		new_ind->setOffset(stream.tellg());
		new_ind->setID(famID + " " + ID);
		//////////////////////////////////////////////////
		loadCompleteMarkerSet(&new_ind);
		//////////////////////////////////////////////////
		individualsP->addIndividual(new_ind);
	}
}


void PEDIndividualsExtractor::loadCompleteMarkerSet(Individual ** p)
{
	stream.seekg(p[0]->getOffset() + 4*ALL_SNPS.getROIStart().getMarkerNumber());

	unsigned int maxsize = ALL_SNPS.currentSize();
	vector<bool>* buffer[2];
	buffer[0] = new vector<bool>(maxsize,false);
	buffer[1] = new vector<bool>(maxsize,false);

	for (unsigned int position = 0; position <  maxsize; position++)
	{
		for(int al=0;al<2;al++){						
			stripWhitespace();
			char marker = stream.peek();

			if ( ALL_SNPS.mapNucleotideToBinary(marker,  position ) == 1 )
				buffer[al]->at(position)=true;
			else
				buffer[al]->at(position)=false;
			stream.get();
		}
	}
	if(HAPLOID)	
	{
		p[0]->addMarkers(TRANS,buffer[0]);
		p[1]->addMarkers(TRANS,buffer[1]);	
	}
	else
	{
		p[0]->addMarkers(UNTRANS,buffer[0]);
		p[0]->addMarkers(TRANS, buffer[1]);
	}
	delete buffer[0];
	delete buffer[1];
}

void PEDIndividualsExtractor::updateMarkerSet(Individual * p,unsigned int start,unsigned int end)
{
	p->updateMarkerSet(start,end);
}

void PEDIndividualsExtractor::appendMarkerSet(Individual * p,unsigned int end, int num_markers)
{
	p->appendMarkerSet(end,num_markers);
}

/////////////////////////////////////////////////////////////////////////////////////

//:TO BE OMITTED
void PEDIndividualsExtractor::getCompleteMarkerSet(Individual * p)
{
	stream.seekg(p->getOffset() + 4*ALL_SNPS.getROIStart().getMarkerNumber() + 4*position_ms*MARKER_SET_SIZE + 1);
	MarkerSet * ms[2];
	ms[0] = new MarkerSet();
	ms[1] = new MarkerSet();

//Memorytally  
mem_markers += (2*(sizeof(MarkerSet) + (ms[0]->getMarkerBits().num_blocks() * sizeof(unsigned long))));

	readMarkerSet( ms );

	p->addMarkerSet(UNTRANS,ms[0]);
	p->addMarkerSet(TRANS,ms[1]);
}

//:TO BE OMITTED
void PEDIndividualsExtractor::readMarkerSet( MarkerSet ** ms )
{
	unsigned int maxsize = ALL_SNPS.currentSize();
	for (int position = 0; position < MARKER_SET_SIZE; position++)
	{
		if(position_ms*MARKER_SET_SIZE+position >= maxsize) break;
		for(int al=0;al<2;al++){
			stripWhitespace();
			char marker = stream.peek();
			if ( ALL_SNPS.mapNucleotideToBinary(marker,position_ms*MARKER_SET_SIZE+position) == 1 )
				ms[al]->set(position , true );
			stream.get();
		}
	}
}

//:TO BE OMITTED
void PEDIndividualsExtractor::getCompleteMarkerSet(Individual * p0 , Individual * p1 )
{
	stream.seekg(p0->getOffset() + 4*ALL_SNPS.getROIStart().getMarkerNumber() + 4*position_ms*MARKER_SET_SIZE + 1);
	MarkerSet * ms[2];
	ms[0] = new MarkerSet();
	ms[1] = new MarkerSet();

	readMarkerSet( ms );

	p0->addMarkerSet(TRANS,ms[0]);
	p1->addMarkerSet(TRANS,ms[1]);
}


//:TOBE OMITTED
void PEDIndividualsExtractor::loadMarkerSet( MarkerSet ** ms )
{
	unsigned int maxsize = ALL_SNPS.currentSize();
	for (unsigned int position = 0; position <  maxsize; position++)
	{
		for(int al=0;al<2;al++){
			stripWhitespace();
			char marker = stream.peek();

			if ( ALL_SNPS.mapNucleotideToBinary(marker,  position ) == 1 )
				ms[al]->pushback(true );
			else
				ms[al]->pushback(false );

			stream.get();
		}
	}
}



// end PEDIndividualsExtractor.cpp

