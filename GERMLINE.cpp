// GERMLINE.cpp: GEnetic Relationship Miner for LINear Extension

#include "GERMLINE.h"
#include <iostream>
#include <math.h>

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
void GERMLINE::mine( string params, string map, string ped,string outfile)
{
	PolymorphicIndividualsExtractor * pie = inputManager.getPie();
	inputManager.getIndividuals(map, ped);
	if ( ! pie->valid() ) return;

	string out;
	if (outfile == "")
		out = inputManager.getOutput(); 
	else
		out = outfile;

	num_samples = 0;						
	num_matches = 0;

	pie->loadInput();
	MatchesBuilder mb( pie );
	
	unsigned long long init_mem = (unsigned long long)(mem_all_matches+mem_bufferchr+mem_chromosome+mem_ind+mem_inds+mem_matchfactory+mem_markers+mem_snps+mem_window);
	mem_expected_data= mb.calculateMemData();

	if ( (init_mem+mem_expected_data) < MEM_BOUND)
	{
	cerr<<"\n\n\tBeginning Analysis"<<endl;
	
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
	fout << setw(50) << "Allowed mismatch: " << MAX_ERR_HOMp << "% " << MAX_ERR_HETp << "%"<<endl;
	fout << setw(50) << "Minimum Word size: " << MIN_WINDOW_SIZE << endl;
	if ( ROI )
		fout << setw(50) << "Target region: " << ALL_SNPS.getROIStart().getSNPID() << " - " << ALL_SNPS.getROIEnd().getSNPID() << endl;
	else
		fout << setw(50) << "Target region: " << "all" << endl;
	
	time_t timer[2]; time( &timer[0] );
	
	if ( DEBUG ) cout << "DEBUG MODE ON" << endl;

	if ( ROI )
	{
		ALL_SNPS.beginChromosome();
		ALL_SNPS_CURRENT_SIZE = ALL_SNPS.currentSize();
		mb.buildMatches();
		if ( !SILENT ) cout << "\nMatches completed ... freeing memory" << endl;
		ALL_SAMPLES.freeMatches();
		ALL_SAMPLES.freeMarkers();
		WINDOWS_LIST.clear();
	}
	else
	{	
		for ( ALL_SNPS.beginChromosome() ; ALL_SNPS.moreChromosome() ; ALL_SNPS.nextChromosome() )
		{
			MAX_WINDOW_SIZE=0;
			ALL_SNPS_CURRENT_SIZE = ALL_SNPS.currentSize();
			mb.buildMatches();
			if ( !SILENT ) cout << "\nMatches completed ... freeing memory" << endl;
			ALL_SAMPLES.freeMatches();
			ALL_SAMPLES.freeMarkers();
			WINDOWS_LIST.clear(); LAST_SET=false; 
		}		
	}
	time( &timer[1] );

	fout << setw(50) << "Total IBD segments: " << num_matches << endl;
	fout << setw(50) << "Total runtime (sec): " << difftime( timer[1] , timer[0] ) << endl;
	fout.close();
	MATCH_FILE.close();

	if ( BINARY_OUT )		
	{
		ofstream bmid_out( ( out + ".bmid" ).c_str() );
		ALL_SNPS.print( bmid_out );
		bmid_out.close();

		ofstream bsid_out( ( out + ".bsid" ).c_str() );
		ALL_SAMPLES.print( bsid_out );
		bsid_out.close();
	}}
	else
	{
		cerr<<"\n\n\tNot enough memory to load data...Cannot begin analysis";
		cerr<<"\n\tRequire atleast "<<ceil((float)(init_mem+mem_expected_data)/1048576)<<" Mb space"<<endl;
	}
}


// end GERMLINE.cpp
