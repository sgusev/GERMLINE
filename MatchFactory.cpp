// MatchFactory.cpp: Generates matches from individuals

#include "MatchFactory.h"


// MatchFactory(): default constructor
MatchFactory::MatchFactory()
{}

void MatchFactory::assertShares()
{
	for ( iter = segments.begin() ; iter != segments.end() ; iter++ )
	{
		iter->second.assertMatches();
	}
}

int MatchFactory::size()
{
	return (int)segments.size();
}

// initialize(): initializes object
void MatchFactory::initialize()
{
	segments.clear();
}	
void MatchFactory::hash( Individual * i )
{
	int haps , het = i->numHet();
	
	if ( het == 0 || HAPLOID ) haps = 1;
	else haps = 2;

	if(VAR_WINDOW)
	{
		if ( ALLOW_HOM && het <= WINDOWS_LIST.err_hom(position_ms) + WINDOWS_LIST.err_het(position_ms) ) 
		i->assertHomozygous();}
	else
	{
		if ( ALLOW_HOM && het <= MAX_ERR_HOM + MAX_ERR_HET ) 
		i->assertHomozygous();}

	if ( HOM_ONLY ) 
		return;
	
	for ( int c = 0 ; c < haps ; c++ )
	{
		boost::dynamic_bitset<>& ms = i->getChromosome(c)->getMarkerSet()->getMarkerBits();
		if ( (iter = segments.find( ms )) == segments.end() )
			segments.insert( make_pair ( ms , Share( i ) ) );
		else
			iter->second.add( i );
	}
}


unsigned int MatchFactory::calculateMem()
{
  unsigned int mem=0, temp=0;
  map < boost::dynamic_bitset<> , Share >::iterator it;
  cout<<"\nposition_ms= "<<position_ms<<"\t"<<"segments= "<<segments.size()<<"\twindow size= "<<segments.begin()->first.size();
  
  for( it=segments.begin(); it!=segments.end(); ++it)
	  temp += (Nchoose2(it->second.size()) * sizeof(Match)) ;
  
  mem+=temp;
  cout<<"-- "<<mem<<" bytes";
  return mem;
}

unsigned int MatchFactory::Nchoose2(unsigned int num)
{
	if (num<2) return 0;
	else return num* (num-1) / 2;
}



// end MatchFactory.cpp
