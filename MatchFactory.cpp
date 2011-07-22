// MatchFactory.cpp: Generates matches from individuals

#include "MatchFactory.h"

// MatchFactory(): default constructor
MatchFactory::MatchFactory()
{mem_matchfactory = sizeof(MatchFactory);}

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
	mem_matchfactory = sizeof(MatchFactory);
	segments.clear();
}	
void MatchFactory::hash( Individual * i )
{
	int haps , het = i->numHet();
	
	if ( het == 0 || HAPLOID ) haps = 1;
	else haps = 2;

	if ( ALLOW_HOM && het <= WINDOWS_LIST.err_hom(position_ms) + WINDOWS_LIST.err_het(position_ms) ) 
		i->assertHomozygous();
	if ( HOM_ONLY ) 
		return;
	
	for ( int c = 0 ; c < haps ; c++ )
	{
		boost::dynamic_bitset<>& ms = i->getChromosome(c)->getMarkerSet()->getMarkerBits();
		if ( (iter = segments.find( ms )) == segments.end() )
		{
			segments.insert( make_pair ( ms , Share( i ) ) );
			mem_matchfactory+= sizeof(pair<boost::dynamic_bitset<>,Share>)
							+(ms.num_blocks()*sizeof(unsigned long))
							+sizeof(Share)
							+(3*sizeof(Individual*));
		}
		else
		{	
			iter->second.add( i );
			mem_matchfactory+= (3*sizeof(Individual*));
		}
	}
}

unsigned long long MatchFactory::calculateMemData()
{
	unsigned long long mem=0;
	int num_windows=0;
	if (num_sets<1)  
	{
		num_windows = ceil((float)ALL_SNPS.currentSize()/ MIN_WINDOW_SIZE);
 		MarkerSet ms = new MarkerSet(false);
        for(int i=0; i<MIN_WINDOW_SIZE; i++)
        	ms.pushback(true);
		mem_one_markerset = (sizeof(MarkerSet) + (ms.getMarkerBits().num_blocks() * sizeof(boost::dynamic_bitset<>::block_type))) ;
	}
	else
		num_windows = ceil((float)(ALL_SNPS.currentSize()- position_marker)/MIN_WINDOW_SIZE);
	
	mem+= (sizeof(WindowInfo) * num_windows );
	mem+= ( num_samples * 2 * num_windows * mem_one_markerset);
	return mem;
}

unsigned long long MatchFactory::calculateMem()
{
  unsigned long long mem=0;
  map < boost::dynamic_bitset<> , Share >::iterator it;
  
  for( it=segments.begin(); it!=segments.end(); ++it)
	  mem += (unsigned long long )(Nchoose2(it->second.size()) * (sizeof(Match)+sizeof(pair<size_t,Match*>))) ;

  return (unsigned long long) mem+mem_all_matches+mem_bufferchr+mem_chromosome+mem_ind+mem_inds+mem_matchfactory+mem_markers+mem_snps+mem_window;
}

unsigned long long  MatchFactory::Nchoose2(unsigned int num)
{
	if (num<2) return 0;
	else return (unsigned long long) num* (num-1) / 2;
}

// end MatchFactory.cpp
