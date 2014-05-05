#ifndef TRIGGER_TOWER_GEO
#define TRIGGER_TOWER_GEO

class TriggerTowerGeometry
{
  public:
	TriggerTowerGeometry(  )
	{

	        // Calculate the TT phi width, 2*PI/PhiTTs
	        mDeltaPhi = 6.28318530718/72.;

		// Store the TT eta widths
		for ( int i = 1; i <= 20; ++i )
			mMappingeta[i] = 0.087;

		mMappingeta[21] = 0.09;
		mMappingeta[22] = 0.1;
		mMappingeta[23] = 0.113;
		mMappingeta[24] = 0.129;
		mMappingeta[25] = 0.15;
		mMappingeta[26] = 0.178;
		mMappingeta[27] = 0.15;
		mMappingeta[28] = 0.35;
	}

	// Returns the eta of the centre of the TT indexed by iEta
	double eta( const int& iEta )
	{
		double eta = 0;
		for ( int i = 1; i <= abs( iEta ); ++i )
		{
			eta += mMappingeta[i];
		}
		eta -= mMappingeta[abs( iEta )] / 2;

		if ( iEta > 0 )
			return eta;
		else
			return -eta;
	}

	double phi( const int& iPhi )
	{
	  // Restrict phi to the range: [ -PI, PI ]
	  if ( iPhi < 37 ){
	    return ( double( iPhi ) - 0.5 ) * mDeltaPhi;
	  }
	  else{
	    return ( double( iPhi ) - 0.5 - 72 ) * mDeltaPhi;
	  }
	}

	// Convert phi (radians) to iPhi
	int iPhi( const double& Phi )
	{ 
	  double PhiDeg = Phi * 57.2957795131;
	  int iPhi = int( PhiDeg/5 ) + 1;
	  return iPhi;
	}
	// Convert eta to iEta
	int iEta( const double& Eta )
	{ 
	  // Return error
	  if ( fabs(Eta) > 3){
	    return 0;
	  }

	  int ieta;
	  double fabsEta = fabs( Eta );
	  if ( fabsEta < 0.087*20){
	    ieta = fabsEta/0.087 + 1;
	  }
	  else{
	    
	    double overEta = fabsEta - 20*0.087;
	    int etaBin = 20;
	    while ( overEta > 0 ){
	      etaBin++;
	      overEta -= mMappingeta[ etaBin ];
	    }
	    ieta = etaBin;

	  }


	  if ( Eta < 0 ){
	    return -ieta;
	  }
	  else{
	    return ieta;
	  }
	}

	// Assumes iPhi is valid, abs(deltaIPhi) < 72                                                                                                                                                   
	int addToiPhi( int iPhi, int deltaIPhi ){
	  int newIPhi = iPhi + deltaIPhi;
	  if ( deltaIPhi > 0){
	    if (newIPhi > 72){ newIPhi -= 72;}
	  }
	  if ( deltaIPhi < 0){
	    if (newIPhi < 1){  newIPhi += 72;}
	  }
	  return newIPhi;
	}
	// Assumes iEta is valid, no limit is set on returned iEta                                                                                                                                      
	int addToiEta( int iEta, int deltaIEta ){
	  int newIEta = iEta + deltaIEta;
	  if ( deltaIEta > 0){
	    if ( (iEta < 0) && (deltaIEta >= 0) ){ deltaIEta++; }
	  }
	  if ( deltaIEta < 0){
	    if ( (iEta > 0) && (deltaIEta <= 0) ){ deltaIEta--; }
	  }
	  return newIEta;
	}



	double towerEtaSize( const int& iEta )
	{
		return mMappingeta[abs( iEta )];
	}

	double towerPhiSize( const int& iPhi )
	{
	  return mDeltaPhi;
	}

	double towerArea( const int& iEta, const int& iPhi )
	{
	  return towerEtaSize( iEta ) * towerPhiSize( iPhi );
	}

  private:
	std::map < int, double >mMappingeta;
	double mDeltaPhi;
};

#endif
