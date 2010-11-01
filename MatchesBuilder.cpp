// MatchesBuilder.cpp: builds matches from individuals

#include "MatchesBuilder.h"
#include "list"

using namespace std;

unsigned int position_ms;
unsigned int num_sets;
list<WindowInfo*> WINDOWS_LIST;
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
	if(VAR_WINDOW)   // generate window sizes in the range 100 to 149
	{
//	cout<<"\n\t\tGenerate Window Sizes ";			
	unsigned int pos =0, r=0;
	while (pos < ALL_SNPS_CURRENT_SIZE)
	{   // cout <<"\nSTART position = "<<pos ;
		if ((ALL_SNPS_CURRENT_SIZE - pos ) <= 100)  	
		{
			//cout<<"\tSET SIZE = "<< ALL_SNPS_CURRENT_SIZE-pos;  
			WINDOWS_LIST.push_back(new WindowInfo(pos,ALL_SNPS_CURRENT_SIZE,ALL_SNPS_CURRENT_SIZE-pos));
			pos+= (ALL_SNPS_CURRENT_SIZE-pos);  
		}
		else 
		{
			r= (unsigned int) (rand() % 50 +100); 
			//cout<<"\tSET SIZE = "<< r; 
			WINDOWS_LIST.push_back(new WindowInfo(pos,pos+r,r));
			pos+= r; 			
		}
		num_sets++;// cout <<"\tEND position = "<<pos ;
	}
	cout<<"\n\t\tNum of Sets : "<<num_sets;
	_flushall();
	cin.get();
	}
	if ( !SILENT ) cout << "Read Markers" << endl;
	ms_start = 0; ms_end = num_sets;
	readAllMarkers();

	if ( !SILENT ) cout << "Match Markers" << endl;
	matchAllMarkers();

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
	for (position_ms = ms_start; position_ms < ms_end ; position_ms++)
	{
		if ( !SILENT ) cerr << "\rMatching Markers - " << (position_ms*100) / (ms_end - 1) << "%" << flush;
		matchMarkerSet();
	}
	if ( !SILENT ) cerr << '\r' << "Matching Markers Complete" << endl;
}

void MatchesBuilder::readAllMarkers()
{
	list<WindowInfo*>::iterator wind_it ;
	if (VAR_WINDOW) wind_it = WINDOWS_LIST.begin();
	
	for (position_ms = ms_start; position_ms < ms_end; position_ms++)
	{
		if ( !SILENT ) cerr << "\rReading Markers - " << position_ms*100/ms_end << "%" << flush;
		if ( HAPLOID ) readHaploidMarkerSet(); 
		else 
		{ 
			if (VAR_WINDOW)
			{ 
				readMarkerSet(((WindowInfo*)*wind_it)->getStart(), ((WindowInfo*)*wind_it)->getEnd());
				wind_it++;
			}
			else readMarkerSet(0,0);
		}

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

// end MatchesBuilder.cpp
