#include "GERMLINE.h"
#include <string.h>
#include <math.h>

unsigned long long mem_all_matches=0;
unsigned long long mem_bufferchr=0;
unsigned long long mem_chromosome=0;
unsigned long long mem_ind=0;
unsigned long long mem_inds= sizeof(Individuals);
unsigned long long mem_matchfactory = sizeof(MatchFactory);
unsigned long long mem_markers=0;
unsigned long long mem_snps=sizeof(ALL_SNPS);
unsigned long long mem_window=0;
unsigned long long mem_expected_data=0;
unsigned short mem_one_markerset=0;

unsigned int ALL_SNPS_CURRENT_SIZE=0;
unsigned int MAX_WINDOW_SIZE=0;
unsigned int MIN_WINDOW_SIZE = 128;  // initial default window size = 128 SNPs
unsigned long long  MEM_BOUND = (unsigned long long) 8*1024*1024*1024;				// in bytes, default 8Gb
int WINDOW_FACTOR = 1;			//factor (in markers) by which a window length is increased when mem bound is crossed
bool LAST_SET=false;

double MIN_MATCH_LEN = 3;
int MARKER_SET_SIZE = 128;
bool PRINT_MATCH_HAPS = false;
bool ROI = false;
bool HAP_EXT = false;
bool WIN_EXT = false;
bool ALLOW_HOM = false;
bool HOM_ONLY = false;
bool HAPLOID = false;
bool SILENT = false;
bool DEBUG = false;
bool BINARY_OUT = false;
float MAX_ERR_HOMp = 5.0;				//in percentage: default 5 %
float MAX_ERR_HETp = 5.0;				//in percentage: default 5 %
int MAX_ERR_HOM = 1;					//in SNPS: default 1 SNP
int MAX_ERR_HET = 1;					//in SNPS: default 1 SNP


// main(): runs GERMLINE
int main(int argc, char* argv[])
{
	// parse arguments
	string rs_range[2] , gdmap, map, ped,out; out = gdmap = map = ped = rs_range[0] = rs_range[1] = "";
	string params = argv[0];

	bool bad_param = false;
	for(int i=1;i<argc;i++){
		params += " " + string(argv[i]);
		if( strncmp(argv[i], "-min_m", strlen("-min_m")) == 0 && i < argc-1)				MIN_MATCH_LEN = atof(argv[++i]);
		else if( strncmp(argv[i], "-err_hom", strlen("-err_hom")) == 0 && i < argc-1)		MAX_ERR_HOMp = atoi(argv[++i]); 
		else if( strncmp(argv[i], "-err_het", strlen("-err_het")) == 0 && i < argc-1)		MAX_ERR_HETp = atoi(argv[++i]); 
		else if( strncmp(argv[i], "-from_snp", strlen("-from_snp")) == 0 && i < argc-1 )	rs_range[0] = argv[++i];
		else if( strncmp(argv[i], "-to_snp", strlen("-to_snp")) == 0 && i < argc-1 )		rs_range[1] = argv[++i];
		else if( strncmp(argv[i], "-print", strlen("-print")) == 0 )						PRINT_MATCH_HAPS = true;
		else if( strncmp(argv[i], "-silent", strlen("-silent")) == 0 )						SILENT = true;
		else if( strncmp(argv[i], "-debug", strlen("-debug")) == 0 )						DEBUG = true;
		else if( strncmp(argv[i], "-gdmap", strlen("-gdmap")) == 0 && i < argc-1)			gdmap = argv[++i];
		else if( strncmp(argv[i], "-map", strlen("-map")) == 0 && i < argc-1)				map = argv[++i];
		else if( strncmp(argv[i], "-ped", strlen("-ped")) == 0 && i < argc-1)				ped = argv[++i];
		else if( strncmp(argv[i], "-out", strlen("-out")) == 0 && i < argc-1)				out = argv[++i];
		else if( strncmp(argv[i], "-bits", strlen("-bits")) == 0 && i < argc-1)				{MARKER_SET_SIZE = atoi(argv[++i]);
																							MIN_WINDOW_SIZE=MARKER_SET_SIZE;}
		else if( strncmp(argv[i], "-homoz-only", strlen("-homoz-only")) == 0 )				{ ALLOW_HOM = true; HOM_ONLY = true; }
		else if( strncmp(argv[i], "-homoz", strlen("-homoz")) == 0 )						ALLOW_HOM = true;
		else if( strncmp(argv[i], "-bin_out", strlen("-bin_out")) == 0 )					BINARY_OUT = true;
		else if( strncmp(argv[i], "-haploid", strlen("-haploid")) == 0 )					{ HAPLOID = true; HAP_EXT = true; }
		else if( strncmp(argv[i], "-h_extend", strlen("-h_extend")) == 0 )					HAP_EXT = true;
		else if( strncmp(argv[i], "-w_extend", strlen("-w_extend")) == 0 )					WIN_EXT = true;
		else if( strncmp(argv[i], "-membound", strlen("-membound")) == 0 && i < argc-1)		MEM_BOUND = (unsigned long long) atoi(argv[++i])*1024*1024;
		else bad_param = true;
	}
	
	MAX_ERR_HOM = (int) floor(MARKER_SET_SIZE*MAX_ERR_HOMp/100);
	MAX_ERR_HET = (int) floor(MARKER_SET_SIZE*MAX_ERR_HETp/100);

	if(MIN_MATCH_LEN < 0)
	{
		cerr << "-min_m must be non-negative" << endl << endl;
		bad_param = true;
	} else if(MAX_ERR_HOMp < 0 || MAX_ERR_HETp < 0 || MAX_ERR_HOMp > 100 || MAX_ERR_HETp > 100 )
	{
		cerr << "-err_hom,-err_het must be within 0-100%" << endl << endl;
		bad_param = true;
	}

	if(bad_param)
	{
		cerr << "usage: " << argv[0] << "<flags (optional)>" << endl
		<< "flags:" << endl
		<< '\t' << "-silent" << '\t' << "Suppress all output except for warnings and prompts." << endl
		<< '\t' << "-bin_out" << '\t' << "Output in binary format to save space." << endl
		<< '\t' << "-min_m" << '\t' << "Minimum length for match to be used for imputation (in cM or MB)." << endl
		<< '\t' << "-err_hom" << '\t' << "Maximum percentage (0-100%) of mismatching homozygous markers (per slice)." << endl
		<< '\t' << "-err_het" << '\t' << "Maximum percentage (0-100%) of mismatching heterozygous markers (per slice)." << endl
		<< '\t' << "-from_snp" << '\t' << "Start SNP (rsID)." << endl
		<< '\t' << "-to_snp" << '\t' << "End SNP (rsID)." << endl
		<< '\t' << "-haps" << '\t' << "Print the resolved haplotypes in a seperate HAPS file." << endl
		<< '\t' << "-gdmap" << '\t' << "Genetic distance map." << endl
		<< '\t' << "-map" << '\t' << "MAP file." << endl
		<< '\t' << "-ped" << '\t' << "PED file." << endl
		<< '\t' << "-bits" << '\t' << "Slice size." << endl
		<< '\t' << "-homoz" << '\t' << "Allow self matches (homozygosity)" << endl
		<< '\t' << "-homoz-only" << '\t' << "Look for autozygous/homozygous segments only, does not detect IBD" << endl
		<< '\t' << "-haploid" << '\t' << "Treat input individual as two fully phased chromosomes with no recombination\n\t\toutput IDs with 0/1 suffix for chromosome destinction" << endl
		<< '\t' << "-h_extend" << '\t' << "Extend from seeds if *haplotypes* match" << endl
		<< '\t' << "-w_extend" << '\t' << "Extend, one marker at a time, beyong the boundaries of a found match" << endl
		<< '\t' << "-membound" << '\t' << "Memory bound in Mb" << endl;
		return 0;
	}

	if( rs_range[0] != "" && rs_range[1] != "" )
	{
		ROI = true;
		ALL_SNPS.setROI(rs_range);
	}

	if(gdmap != "")
	{
		ALL_SNPS.loadGeneticDistanceMap( gdmap );
	}

	GERMLINE germline;
    germline.mine(params,map, ped,out );
	
    return 1;
}

// end GERMLINE_0001.cpp
