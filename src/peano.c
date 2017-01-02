#include "globals.h"
#include "peano.h"
#include "sort.h"

peanoKey Peano_Key(const double x, const double y, const double z);
static void reorder_particles();

void Print_Int_Bits128(const peanoKey val)
{
    for (int i = 127; i >= 0; i--) {

        printf("%llu",(long long) ( (val & ((peanoKey)1 << i) ) >> i));

        if (i % 3 == 0 && i != 0)
            printf(".");
    }
    printf("\n");fflush(stdout);

    return ;
}

void Print_Int_Bits128r(const peanoKey val)
{
    for (int i = 127; i >= 0; i--) {

        printf("%llu", (long long) ((val & ((peanoKey) 1 << i) ) >> i));

        if (i % 3-2 == 0 && i != 0)
            printf(".");
    }
    printf("\n");fflush(stdout);

    return ;
}

int compare_peanoKeys(const void * a, const void *b)
{
    const peanoKey *x = (const peanoKey *) a;
    const peanoKey *y = (const peanoKey *) b;

    return (int) (*x > *y) - (*x < *y);
}

static peanoKey *Keys = NULL;
static size_t *Idx = NULL;

void Sort_Particles_By_Peano_Key()
{
    if (Keys == NULL)
        Keys = malloc(Param.Npart * sizeof(*Keys));
    else
        memset(Keys, 0, Param.Npart * sizeof(*Keys));

    if (Idx == NULL)
        Idx = malloc(Param.Npart * sizeof(*Idx));
    else
        memset(Idx, 0, Param.Npart * sizeof(*Idx));

	#pragma omp parallel for
    for (int ipart = 0; ipart < Param.Npart; ipart++) {

        double px = P[ipart].Pos[0] / Problem.Boxsize[0];
        double py = P[ipart].Pos[1] / Problem.Boxsize[1];
        double pz = P[ipart].Pos[2] / Problem.Boxsize[2];

        P[ipart].Key = Keys[ipart] = Peano_Key(px, py, pz);
    }

	#pragma omp parallel
    Qsort_Index(1, Idx, Keys, Param.Npart, sizeof(*Keys),
           	 	&compare_peanoKeys);

    reorder_particles();

    return ;
}

static void reorder_particles()
{
    for (int i = 0; i < Param.Npart; i++) {

        if (Idx[i] == i)
            continue;

        int dest = i;

        struct ParticleData Ptmp = P[i];
        struct GasParticleData Sphtmp = SphP[i];

        int src = Idx[i];

        for (;;) {

            P[dest] = P[src];
            SphP[dest] = SphP[src];

            Idx[dest] = dest;

            dest = src;

            src = Idx[dest];

            if (src == i)
                break;
        }

        P[dest] = Ptmp;
        SphP[dest] = Sphtmp;

        Idx[dest] = dest;

    } // for i

    return ;
}

peanoKey Peano_Key(const double x, const double y, const double z)
{
    Assert(x >= 0 && x <= 1, "X coordinate of out range [0,1] have %g", x);
    Assert(y >= 0 && y <= 1, "Y coordinate of out range [0,1] have %g", y);
    Assert(z >= 0 && z <= 1, "Z coordinate of out range [0,1] have %g", z);

    const uint64_t m = 1UL << 63; // = 2^63;

    uint64_t X[3] = { y*m, z*m, x*m };

    /* Inverse undo */

    for (uint64_t q = m; q > 1; q >>= 1 ) {

        uint64_t P = q - 1;

        if( X[0] & q )
            X[0] ^= P;  // invert

        for (int i = 1; i < 3; i++ ) {

            if( X[i] & q ) {

                X[0] ^= P; // invert

            } else {

                uint64_t t = (X[0] ^ X[i]) & P;

                X[0] ^= t;
                X[i] ^= t;

            } // exchange
        }
    }

    /* Gray encode (inverse of decode) */

    for (int i = 1; i < 3; i++ )
        X[i] ^= X[i-1];

    uint64_t t = X[2];

    for (int i = 1; i < 64; i <<= 1 )
        X[2] ^= X[2] >> i;

    t ^= X[2];

    for (int i = 1; i >= 0; i-- )
        X[i] ^= t;

    /* branch free bit interleave of transpose array X into key */

    peanoKey key = 0;

    X[1] >>= 1; X[2] >>= 2;    // lowest bits not important

    for (int i = 0; i < N_PEANO_TRIPLETS+1; i++) {

        uint64_t col = ((X[0] & 0x8000000000000000)
                	  | (X[1] & 0x4000000000000000)
                	  | (X[2] & 0x2000000000000000)) >> 61;

        key <<= 3;

        X[0] <<= 1;
        X[1] <<= 1;
        X[2] <<= 1;

        key |= col;
    }

    key <<= 2;

    return key;
}

/* This constructs the peano key with reversed triplet order. The order in the
 * triplets however is the same ! Also level zero is carried explicitely
 * to ease tree construction. */

peanoKey Reversed_Peano_Key(const double x, const double y, const double z)
{
    Assert(x >= 0 && x <= 1, "X coordinate of out range [0,1] have %g", x);
    Assert(y >= 0 && y <= 1, "Y coordinate of out range [0,1] have %g", y);
    Assert(z >= 0 && z <= 1, "Z coordinate of out range [0,1] have %g", z);

    const uint64_t m = 1UL << 63; // = 2^63;

    uint64_t X[3] = { y*m, z*m, x*m };

    /* Inverse undo */

    for (uint64_t q = m; q > 1; q >>= 1) {

        uint64_t P = q - 1;

        if(X[0] & q)
            X[0] ^= P;  // invert

        for(int i = 1; i < 3; i++ ) {

            if(X[i] & q) {

                X[0] ^= P; // invert

            } else {

                uint64_t t = (X[0] ^ X[i]) & P;

                X[0] ^= t;
                X[i] ^= t;

            } // exchange
        }
    }

    /* Gray encode (inverse of decode) */

    for(int i = 1; i < 3; i++)
        X[i] ^= X[i-1];

    uint64_t t = X[2];

    for(int i = 1; i < 64; i <<= 1)
        X[2] ^= X[2] >> i;

    t ^= X[2];

    for(int i = 1; i >= 0; i--)
        X[i] ^= t;

    /* branch free reversed (!) bit interleave of transpose array X into key */

    peanoKey key = 0;

    X[0] >>= 18; X[1] >>= 19; X[2] >>= 20;    // lowest bits not important

    for (int i = 0; i < N_PEANO_TRIPLETS+1; i++) {

        uint64_t col = ((X[0] & 0x4) | (X[1] & 0x2) | (X[2] & 0x1));

        key <<= 3;

        key |= col;

        X[0] >>= 1;
        X[1] >>= 1;
        X[2] >>= 1;
    }

    key <<= 3; // include level 0

    return key;
}


void test_peanokey()
{
    const double box[3]  = { 1.0, 1, 1};
    double a[3] = { 0 };
    int order = 1;
    float delta = 1/pow(2.0, order);
    int n = roundf(1/delta);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {

                a[0] = (i + 0.5) * delta / box[0];
                a[1] = (j + 0.5) * delta / box[1];
                a[2] = (k + 0.5) * delta / box[2];

                peanoKey stdkey =  Peano_Key(a[0], a[1], a[2]);

                printf("%g %g %g %llu  \n", a[0], a[1], a[2],
						(long long) stdkey );

                Print_Int_Bits128(stdkey);

                printf("\n");
            }
		}
	}

    return ;
}
