// GERMLINE.cpp: GEnetic Relationship Miner for LINear Extension

#include "GERMLINE.h"
#include <iostream>


using namespace std;

ofstream MATCH_FILE;
size_t num_samples;
unsigned long num_matches;
SNPs ALL_SNPS;
Individuals ALL_SAMPLES;


// GERMLINE(): default constructor
GERMLINE::GERMLINE()
{}

// mine(): main function for GERMLINE
void GERMLINE::mine( string params )
{
	PolymorphicIndividualsExtractor * pie = inputManager.getPie();
	if (ROI) VAR_WINDOW = false;
	inputManager.getIndividuals();
	if ( ! pie->valid() ) return;

	//cout<<"MeM_BOUND = "<<MEM_BOUND<<endl;
	
	//string out = "D:\\ColumbiaUniversity\\Project-Germline\\commit3.1\\GERMLINE\\test\\test_fixed_30_tweak";
	string out = inputManager.getOutput();  

	num_samples = 0;						
	num_matches = 0;

	pie->loadInput();
	MatchesBuilder mb( pie );

	

	ofstream fout( ( out + ".log" ).c_str() );

	fout << setw(65) << setfill('-') << ' ' << endl << setfill(' ');
	fout << " Welcome to GERMLINE, a tool for detecting long segments shared" << endl;
	fout << " by descent between pairs of individuals in near-linear time." << endl;
	fout << endl;
	fout << " For more details, please see the paper [ PMID: 18971310 ]" << endl;
	fout << " or the web-site [ http://www.cs.columbia.edu/~gusev/germline/ ]" << endl;
	fout << endl;
	fout << " GERMLINE was coded by Alexander Gusev and collaborators in " << endl;
	fout << " Itsik Pe'er's Computational Biology Lab at Columbia University" << endl;
	fout << setw(65) << setfill('-') << ' ' << endl << setfill(' ');
	
	if ( BINARY_OUT ) MATCH_FILE.open( ( out + ".bmatch" ).c_str() , ios::binary );
	else MATCH_FILE.open( ( out + ".match" ).c_str() );
	
	fout << params << endl;
	fout << setw(65) << setfill('-') << ' ' << endl << setfill(' ');
	fout << setw(50) << left << "Minimum match length: " << MIN_MATCH_LEN << " cM" << endl;
	if(VAR_WINDOW) fout << setw(50) << "Allowed mismatch: " << MAX_ERR_HOMp << "% " << MAX_ERR_HETp << "%"<<endl;
	else fout << setw(50) << "Allowed mismatching bits: " << MAX_ERR_HOM << " " << MAX_ERR_HET << endl;
	if (!VAR_WINDOW) fout << setw(50) << "Word size: " << MARKER_SET_SIZE << endl;
	if ( ROI )
		fout << setw(50) << "Target region: " << ALL_SNPS.getROIStart().getSNPID() << " - " << ALL_SNPS.getROIEnd().getSNPID() << endl;
	else
		fout << setw(50) << "Target region: " << "all" << endl;
	
	time_t timer[2]; time( &timer[0] );
	
	if ( DEBUG ) cout << "DEBUG MODE ON" << endl;

	//TODO: update to use Variable windows
	if ( ROI )
	{
		ALL_SNPS.beginChromosome();
		num_sets = (long)ceil((double)ALL_SNPS.currentSize()/(double)MARKER_SET_SIZE);
		mb.buildMatches();
		ALL_SAMPLES.freeMatches();
		ALL_SAMPLES.freeMarkers();
	}
	else
	{
		for ( ALL_SNPS.beginChromosome() ; ALL_SNPS.moreChromosome() ; ALL_SNPS.nextChromosome() )
		{
			
			ALL_SNPS_CURRENT_SIZE = ALL_SNPS.currentSize();
			if (!VAR_WINDOW)
				num_sets = (long)ceil((double)ALL_SNPS.currentSize()/(double)MARKER_SET_SIZE);
					
			mb.buildMatches();
			//cout<<"Num_Sets = "<<num_sets<<endl;
			if ( !SILENT ) cout << "Matches completed ... freeing memory" << endl;

			ALL_SAMPLES.freeMatches();
			ALL_SAMPLES.freeMarkers();
			if (VAR_WINDOW) WINDOWS_LIST.clear();
		}
		
	}

	time( &timer[1] );

	fout << setw(50) << "Total IBD segments: " << num_matches << endl;
	fout << setw(50) << "Total runtime (sec): " << difftime( timer[1] , timer[0] ) << endl;
	fout.close();
	MATCH_FILE.close();

	//TODO: update to use Variable windows
	if ( BINARY_OUT )		
	{
		ofstream bmid_out( ( out + ".bmid" ).c_str() );
		ALL_SNPS.print( bmid_out );
		bmid_out.close();

		ofstream bsid_out( ( out + ".bsid" ).c_str() );
		ALL_SAMPLES.print( bsid_out );
		bsid_out.close();
	}
}


// end GERMLINE.cpp
