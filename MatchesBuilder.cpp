// MatchesBuilder.cpp: builds matches from individuals

#include "MatchesBuilder.h"
#include <math.h>

using namespace std;

unsigned int position_ms;
unsigned int num_sets;
WindowsList WINDOWS_LIST;
int ALLOWED_MASKED;


// MatchesBuilder(): default constructor
MatchesBuilder::MatchesBuilder( PolymorphicIndividualsExtractor * pie )
{
	individualsP = & ALL_SAMPLES;
	pieP = pie;
}

// buildMatches(): builds matches from individuals
void MatchesBuilder::buildMatches()
{
	if(VAR_WINDOW)   
	{
		WINDOWS_LIST.initialize(ALL_SNPS_CURRENT_SIZE);
		position_ms = -1;
		do {
			readMatchMarkerSet(); 
		}while(WINDOWS_LIST.getWindowEnd() != ALL_SNPS_CURRENT_SIZE);
	}
	else
	{
	if ( !SILENT ) cout << "Read Markers" << endl;
	ms_start = 0; ms_end = num_sets;
	readAllMarkers();

	if ( !SILENT ) cout << "Match Markers" << endl;
	matchAllMarkers();
	}
}

void MatchesBuilder::printHaplotypes(string fout_name)
{
	ofstream fout( fout_name.c_str() );

	for(individualsP->begin();individualsP->more();)
	{
		individualsP->next()->print(fout,0,num_sets);
	}
	fout.close();
}

// matchAllMarkers(): builds matches for individuals considering all markers
void MatchesBuilder::matchAllMarkers()
{
	for (position_ms = ms_start; position_ms < ms_end -1 ; position_ms++)
	{
		if ( !SILENT ) cerr << "\rMatching Markers - " << (position_ms*100) / (ms_end - 1) << "%" << flush;
		matchMarkerSet();
	}
	if ( !SILENT ) cerr << '\r' << "Matching Markers Complete" << endl;
}

void MatchesBuilder::readAllMarkers()
{    
	for (position_ms = ms_start; position_ms < ms_end; position_ms++)
	{
		if ( !SILENT ) cerr << "\rReading Markers - " << position_ms*100/ms_end << "%" << flush;
		if ( HAPLOID && ROI) readHaploidMarkerSet();	//HAPLOID && ROI , :TO BE OMITTED
		else 
		{ 
			if (VAR_WINDOW)					//HAPLOID or DIPLOID but not ROI, TODO: include ROI 
				readMarkerSet(WINDOWS_LIST.getWindowStart(position_ms), WINDOWS_LIST.getWindowEnd(position_ms));
			else							//DIPLOID and ROI , :TO BE OMITTED 
				readMarkerSet(0,0);
		}
	}
	if ( !SILENT ) cerr << '\r' << "Reading Markers Complete" << endl;
}

//:TO BE OMITTED
void MatchesBuilder::readHaploidMarkerSet()
{
	// Read the individuals two at a time

	Individual * i[2];
	for(individualsP->begin();individualsP->more();)
	{
		i[0] = individualsP->next();
		i[1] = individualsP->next();
		pieP->getCompleteMarkerSet( i[0] , i[1] );
	}
}

void MatchesBuilder::readMarkerSet(unsigned int start, unsigned int end)
{
	Individual * i;
	for(individualsP->begin();individualsP->more();)
	{
		i = individualsP->next();
		if (VAR_WINDOW) pieP->updateMarkerSet(i,start,end);			
		else pieP->getCompleteMarkerSet(i);
	}
}

// matchMarkerSet(): builds matches for individuals considering markers in marker set.
void MatchesBuilder::matchMarkerSet()
{

	// Match:
	for(individualsP->begin();individualsP->more();)
		matchFactory.hash( individualsP->next() );
	
	// Verify:
	matchFactory.assertShares();

	// Extend:
	for(individualsP->begin();individualsP->more();)
		individualsP->next()->assertShares();
	
	matchFactory.initialize();
}

void MatchesBuilder::readMatchMarkerSet()
{
	WINDOWS_LIST.getNewWindowSize(ALL_SNPS_CURRENT_SIZE);
	
	//Read:
	position_ms++;
	readMarkerSet(WINDOWS_LIST.getWindowStart(), WINDOWS_LIST.getWindowEnd());

	//Match:
	bool flag_updateWindow = true;
	while(flag_updateWindow)
	{
		for(individualsP->begin();individualsP->more();)
			matchFactory.hash( individualsP->next() );
		
		//Check Memory Bound:
		if ( matchFactory.calculateMem() >= MEM_BOUND)			
		{//If not enough memory
			cout<<"--Mem Bound crossed ";
			matchFactory.initialize();
			// Update WindowSize and MarkerSet 
			int num_markers = WINDOWS_LIST.updateWindowSize(ALL_SNPS_CURRENT_SIZE);
			if(num_markers == 0 )	// reached end of total markers - no hashing
			{ position_ms--; return; }
			updateMarkerSet(num_markers);
			
		}
		else flag_updateWindow=false;
	}

	//Verify:

	matchFactory.assertShares();

	// Extend:
	for(individualsP->begin();individualsP->more();)
		individualsP->next()->assertShares();
	
	matchFactory.initialize();

	/*_flushall();
	char ch = getchar();*/
	
}

void MatchesBuilder::updateMarkerSet(int num_markers)
{
	Individual * i;
	for(individualsP->begin();individualsP->more();)
	{
		i = individualsP->next();
		pieP->updateMarkerSet(i,num_markers);			
	}
}


// end MatchesBuilder.cpp
