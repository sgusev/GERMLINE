#include "Match.h"

void Match::extendBack()
{
	// save old position values
	unsigned int SAVE_pms = position_ms;
	
	// iterate backwards through genome
	while(position_ms > 0)
	{
		position_ms--;
		if( !approxEqual() )
		{
			position_ms++;
			break;
		}
	}
	start_ms = position_ms;
	// restore saved values
	position_ms = SAVE_pms;
}

bool Match::approxEqual()
{
	// homozygosity check
	if ( node[0] == node[1] )
	{
		if ( ALLOW_HOM )
		{
			if(VAR_WINDOW)
			{
			if ( (int) ( node[0]->getChromosome( 0 )->getMarkerSet()->getMarkerBits() ^ node[1]->getChromosome( 1 )->getMarkerSet()->getMarkerBits() ).count() 
				<= ( WINDOWS_LIST.err_hom(position_ms) + WINDOWS_LIST.err_het(position_ms) ) ) return true; else return false;
			}
			else
			{
				if ( (int) ( node[0]->getChromosome( 0 )->getMarkerSet()->getMarkerBits() ^ node[1]->getChromosome( 1 )->getMarkerSet()->getMarkerBits() ).count() 
				 <= ( MAX_ERR_HOM + MAX_ERR_HET ) ) return true; else return false;
			}
		}
		else
		{
			return false;
		}
	}
	else
	{
		// 1. Haplotype extension
		if ( HAPLOID )
		{
			if(VAR_WINDOW)
			{if ( (int)(node[0]->getChromosome( 0 )->getMarkerSet()->getMarkerBits() ^ node[1]->getChromosome( 0 )->getMarkerSet()->getMarkerBits()).count() <= WINDOWS_LIST.err_hom(position_ms) ) return true;}
			else
			{if ( (int)(node[0]->getChromosome( 0 )->getMarkerSet()->getMarkerBits() ^ node[1]->getChromosome( 0 )->getMarkerSet()->getMarkerBits()).count() <= MAX_ERR_HOM ) return true;}
		} else
		{
			for ( int a = 0 ; a < 2 ; a++ ) {
				for ( int b = 0 ; b < 2 ; b++ ) { 
					if(VAR_WINDOW){
					if ( (int)(node[0]->getChromosome( a )->getMarkerSet()->getMarkerBits() ^ node[1]->getChromosome( b )->getMarkerSet()->getMarkerBits()).count() <= WINDOWS_LIST.err_hom(position_ms) )
					{
						return true;
					}}
					else{
					if ( (int)(node[0]->getChromosome( a )->getMarkerSet()->getMarkerBits() ^ node[1]->getChromosome( b )->getMarkerSet()->getMarkerBits()).count() <= MAX_ERR_HOM )
					{
						return true;
					}}
				}
			}
		}

		if ( HAPLOID || HAP_EXT ) return false;

		// 2. Genotype extension
		// identify common homozygous SNPs
		boost::dynamic_bitset<> mask
			= ( node[0]->getChromosome( 0 )->getMarkerSet()->getMarkerBits() ^ node[0]->getChromosome( 1 )->getMarkerSet()->getMarkerBits() ).flip()
			& ( node[1]->getChromosome( 0 )->getMarkerSet()->getMarkerBits() ^ node[1]->getChromosome( 1 )->getMarkerSet()->getMarkerBits() ).flip();

		// assert that homozygous SNPs are identical
		if(VAR_WINDOW)
		{
			if ( (int) ((node[0]->getChromosome( 0 )->getMarkerSet()->getMarkerBits() ^ node[1]->getChromosome( 0 )->getMarkerSet()->getMarkerBits()) & mask).count() <= WINDOWS_LIST.err_het(position_ms) )
		{
			return true;
		}else return false;}
		else
		{
		if ( (int) ((node[0]->getChromosome( 0 )->getMarkerSet()->getMarkerBits() ^ node[1]->getChromosome( 0 )->getMarkerSet()->getMarkerBits()) & mask).count() <= MAX_ERR_HET )
		{
			return true;
		}else return false;}
		
	}
}


int Match::scanLeft( unsigned int ms )
{
	bool err = false;
	int marker;

	boost::dynamic_bitset<> mask
		= ( node[0]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits() ^ node[0]->getChromosome( 1 )->getMarkerSet(ms)->getMarkerBits() ).flip()
		& ( node[1]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits() ^ node[1]->getChromosome( 1 )->getMarkerSet(ms)->getMarkerBits() ).flip();
	mask = ( node[0]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits() ^ node[1]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits()) & mask;

	if(VAR_WINDOW){
		for( marker = WINDOWS_LIST.getWindowSize(ms) - 1 ; marker >= 0 && !err ; marker-- )
			if ( mask[marker] ) err = true;
	}
	else{
		
		for( marker = MARKER_SET_SIZE - 1 ; marker >= 0 && !err ; marker-- )
			if ( mask[marker] ) err = true;
		
	}

	return marker;
}

int Match::scanRight( unsigned int ms )
{
	bool err = false;
	int marker;

	boost::dynamic_bitset<> mask
		= ( node[0]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits() ^ node[0]->getChromosome( 1 )->getMarkerSet(ms)->getMarkerBits() ).flip()
		& ( node[1]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits() ^ node[1]->getChromosome( 1 )->getMarkerSet(ms)->getMarkerBits() ).flip();
	mask = ( node[0]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits() ^ node[1]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits()) & mask;

	if (VAR_WINDOW){
		
		for( marker = 0 ; marker < WINDOWS_LIST.getWindowSize(ms) && !err ; marker++ )
			if ( mask[marker] ) err = true;
	}
	else{
		if (ms == num_sets-1)
		{
			int size = ALL_SNPS_CURRENT_SIZE - ((num_sets-1)*MARKER_SET_SIZE);
			for( marker = 0 ; marker < size && !err ; marker++ )
				if ( mask[marker] ) err = true;
		}
		else
		{
		for( marker = 0 ; marker < MARKER_SET_SIZE && !err ; marker++ )
			if ( mask[marker] ) err = true;
		}
	}
	
	return marker;
}

int Match::diff( unsigned int ms )
{
	boost::dynamic_bitset<> mask
		= ( node[0]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits() ^ node[0]->getChromosome( 1 )->getMarkerSet(ms)->getMarkerBits() ).flip()
		& ( node[1]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits() ^ node[1]->getChromosome( 1 )->getMarkerSet(ms)->getMarkerBits() ).flip();
	mask = ( node[0]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits() ^ node[1]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits()) & mask;

	return int(mask.count());
}

bool Match::isHom( int n , unsigned int ms )
{
	if(VAR_WINDOW)
	return (int) ( node[n]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits() ^ node[n]->getChromosome( 1 )->getMarkerSet(ms)->getMarkerBits() ).count() <= ( WINDOWS_LIST.err_hom(ms) + WINDOWS_LIST.err_het(ms));
	else
	return (int) ( node[n]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits() ^ node[n]->getChromosome( 1 )->getMarkerSet(ms)->getMarkerBits() ).count() <= ( MAX_ERR_HOM + MAX_ERR_HET );
}

void Match::print( ostream& fout )
{
	// TODO: If match length + current window size < minimum match threshold, don't bother extending forward
	
	// extend this match from both ends
	unsigned int snp_start=0, snp_end=0;
	if(VAR_WINDOW) 
	{	
		snp_start = WINDOWS_LIST.getWindowStart(start_ms);
		snp_end = WINDOWS_LIST.getWindowEnd(end_ms) -1;		
	}
	else
	{
		snp_start= ALL_SNPS.getROIStart().getMarkerNumber() + start_ms * MARKER_SET_SIZE;
		snp_end = ALL_SNPS.getROIStart().getMarkerNumber() + ( end_ms + 1 ) * MARKER_SET_SIZE - 1;
	}

	int marker;

	if ( WIN_EXT )								//fixed for VAR_WINDOW
	{
		// backwards
		if( start_ms > 0 )
		{
			marker = scanLeft( start_ms - 1 );
			if(VAR_WINDOW)
				snp_start -= (WINDOWS_LIST.getWindowSize(start_ms) - marker - 2);		//Check this ?
			else
				snp_start -= (MARKER_SET_SIZE - marker - 2);
		}
	}
	if (( LAST_SET && end_ms == num_sets - 2 ) || WIN_EXT)	//fixed for VAR_WINDOW
	{
		// forwards
		if( end_ms < num_sets - 1 )
		{
			marker = scanRight( end_ms + 1 );
			snp_end += marker;

		}
	}
	

	bool genetic;
	float distance;
	if ( ( distance = ALL_SNPS.getDistance(snp_start,snp_end,genetic)) < MIN_MATCH_LEN ) return;
	// print

	// get hamming distance & ignored bit count
	int dif = 0;
	for( unsigned int i = start_ms; i <= end_ms ; i++) { dif += diff( i ); }
	
	// calculate if homozygous
	bool hom[2];
	if ( node[0] == node[1] ) fout << '\t' << 1 << '\t' << 1;
	else
	{
		for ( int n = 0 ; n < 2 ; n++ )
		{
			hom[n] = true;
			for ( unsigned int i = start_ms ; i<= end_ms && hom ; i++ )
			{
				hom[n] = isHom( n , i );
			}
		}
	}

	if ( BINARY_OUT )
	{
		unsigned int pid[2];
		pid[0] = node[0]->getNumericID();
		pid[1] = node[1]->getNumericID();
		unsigned int sid[2];
		sid[0] = ALL_SNPS.getSNP(snp_start).getMarkerNumber();
		sid[1] = ALL_SNPS.getSNP(snp_end).getMarkerNumber();
		fout.write( (char*) &pid[0] , sizeof( unsigned int ) );
		fout.write( (char*) &pid[1] , sizeof( unsigned int ) );
		fout.write( (char*) &sid[0] , sizeof( unsigned int ) );
		fout.write( (char*) &sid[1] , sizeof( unsigned int ) );
		fout.write( (char*) &dif , sizeof( int ) );
		fout.write( (char*) &hom[0] , sizeof( bool ) );
		fout.write( (char*) &hom[1] , sizeof( bool ) );
	} else
	{
		fout << node[0]->getID() << '\t';
		fout << node[1]->getID() << '\t';
		fout << ALL_SNPS.getSNP(snp_start).getChr() << '\t';
		fout << ALL_SNPS.getSNP(snp_start).getPhysPos() << ' ';
		fout << ALL_SNPS.getSNP(snp_end).getPhysPos() << '\t';
		fout << ALL_SNPS.getSNP(snp_start).getSNPID() << ' ';
		fout << ALL_SNPS.getSNP(snp_end).getSNPID() << '\t';
		fout << ( snp_end - snp_start + 1) << '\t';
		fout << setiosflags(ios::fixed) << setprecision(2) << distance << '\t';
		if ( genetic ) fout << "cM" << '\t'; else fout << "MB" << '\t';
		fout << dif;
		for ( int n = 0 ; n < 2 ; n++ )
			if ( hom[n] ) fout << '\t' << 1; else fout << '\t' << 0;
		fout << endl;
	}
	num_matches++;
}

