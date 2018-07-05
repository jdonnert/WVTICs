#define N_PEANO_TRIPLETS (sizeof(peanoKey)*CHAR_BIT/3)

//typedef __uint128_t peanoKey;

void Print_Int_Bits128 ( const peanoKey val );
void Print_Int_Bits128r ( const peanoKey val );

void Sort_Particles_By_Peano_Key();
peanoKey Peano_Key ( const double, const double, const double );
peanoKey Reversed_Peano_Key ( const double, const double, const double );
void test_peanokey();
