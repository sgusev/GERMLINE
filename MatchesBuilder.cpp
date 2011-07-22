// MatchesBuilder.cpp: builds matches from individuals

#include "MatchesBuilder.h"
#include <math.h>

using namespace std;
unsigned int position_marker=0;
int position_ms;
unsigned int num_sets=0;
WindowsList WINDOWS_LIST;
int ALLOWED_MASKED;


// MatchesBuilder(): default constructor
MatchesBuilder::MatchesBuilder( PolymorphicIndividualsExtractor * pie )
{
	individualsP = & ALL_SAMPLES;
	pieP = pie;
}

unsigned long long MatchesBuilder::calculateMemData()
{
	return matchFactory.calculateMemData();
}

// buildMatches(): builds matches from individuals
void MatchesBuilder::buildMatches()
{
	unsigned long long init_mem = (unsigned long long)(mem_all_matches+mem_bufferchr+mem_chromosome+mem_ind+mem_inds+mem_matchfactory+mem_markers+mem_snps+mem_window);
	mem_expected_data= matchFactory.calculateMemData();

	if(init_mem+mem_expected_data >= MEM_BOUND)
	{
		cerr<<"\n\n..NOT ENOUGH MEMORY TO BEGIN ANALYSIS" << endl;
        return;
	}

	WINDOWS_LIST.initialize(ALL_SNPS_CURRENT_SIZE);
	position_ms = -1;
	if ( !SILENT ) cout << "Read-Match Markers" << endl;
	do {
		if ( !SILENT ) cerr << "\rReading-Matching Markers - " << (WINDOWS_LIST.getWindowEnd()*100) / ALL_SNPS_CURRENT_SIZE << "%" << flush;
		readMatchMarkerSet(); 
	}while(WINDOWS_LIST.getWindowEnd() != ALL_SNPS_CURRENT_SIZE);
	if ( !SILENT ) cerr << '\r' << "Reading-Matching Markers Complete" << endl;
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

void MatchesBuilder::readMarkerSet(unsigned int start, unsigned int end)
{
	Individual * i;
	for(individualsP->begin();individualsP->more();)
	{
		i = individualsP->next();
		pieP->updateMarkerSet(i,start,end);			
	}
}


// matchMarkerSet(): builds matches for individuals considering markers in marker set.
void MatchesBuilder::matchMarkerSet()
{

	// Match:
	for(individualsP->begin();individualsP->more();)
		matchFactory.hash( individualsP->next() );
	unsigned long long total_mem =  matchFactory.calculateMem();
	
	// Verify:
	matchFactory.assertShares();

	// Extend:
	for(individualsP->begin();individualsP->more();)
		individualsP->next()->assertShares();
	
	matchFactory.initialize();
}



void MatchesBuilder::readMatchMarkerSet()
{
    unsigned long long total_mem=0;
	WINDOWS_LIST.getNewWindowSize(ALL_SNPS_CURRENT_SIZE);
	
	//Read:
	position_ms++; position_marker=WINDOWS_LIST.getWindowEnd();
	readMarkerSet(WINDOWS_LIST.getWindowStart(), WINDOWS_LIST.getWindowEnd());

	//Match:
	bool flag_updateWindow = true;
	while(flag_updateWindow)
	{
		for(individualsP->begin();individualsP->more();)
			matchFactory.hash( individualsP->next() );

		if (WINDOWS_LIST.getWindowEnd() == ALL_SNPS_CURRENT_SIZE) LAST_SET=true;
		
		//Check Memory Bound:
		total_mem =  matchFactory.calculateMem();
		mem_expected_data = matchFactory.calculateMemData();
        if ( (total_mem+ mem_expected_data) >= MEM_BOUND)
        {//If Germline has run out of memory
			matchFactory.initialize();
            if ( LAST_SET ) {
				cerr<<"\nNOT ENOUGH MEMORY TO HASH REMAINING DATA: ";
                position_ms--;
                return; }
			else{			// Update WindowSize
			int num_markers= WINDOWS_LIST.updateWindowSize(ALL_SNPS_CURRENT_SIZE);
			appendMarkerSet(WINDOWS_LIST.getWindowEnd(),num_markers);
			position_marker=WINDOWS_LIST.getWindowEnd();
			}
		}
		else flag_updateWindow=false;
	}

	//Verify:
	matchFactory.assertShares();

	// Extend:
	for(individualsP->begin();individualsP->more();)
		individualsP->next()->assertShares();
	
	matchFactory.initialize();

}

void MatchesBuilder::appendMarkerSet(unsigned int end,int num_markers)
{
	Individual * i;
	for(individualsP->begin();individualsP->more();)
	{
		i = individualsP->next();
		pieP->appendMarkerSet(i,end,num_markers);			
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////
//:TO BE OMITTED
// matchAllMarkers(): builds matches for individuals considering all markers
void MatchesBuilder::matchAllMarkers()
{
	for (position_ms = ms_start; position_ms < ms_end ; position_ms++)
	{
		if ( !SILENT ) cerr << "\rMatching Markers - " << (position_ms*100) / (ms_end ) << "%" << flush;
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
			readMarkerSet(WINDOWS_LIST.getWindowStart(position_ms), WINDOWS_LIST.getWindowEnd(position_ms));
	}
	if ( !SILENT ) cerr << '\r' << "Reading Markers Complete" << endl;
}


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



// end MatchesBuilder.cpp
