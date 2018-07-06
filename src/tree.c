#include "globals.h"
#include "tree.h"

#define NODES_PER_PARTICLE 4 //2
#ifdef NEW_PEANO
#define N_PEANO_TRIPLETS PEANO_ORDER - 1 //Due to patch in reverse probably requires the -1 ???
#else
#define Print Print_Int_Bits128r
#define Printr Print_Int_Bits128
#endif

struct Tree_Node {
    uint32_t Bitfield;  // bit 0-5:level, 6-8:key, 9:local, 10:top, 11-31:free
    int DNext;          // Distance to the next node; or particle -DNext-1
    float Pos[3];       // Node Center
    int Npart;          // Number of particles in node
    float Size[3];
} *Tree;

int NNodes = 0;
int Max_Nodes = 0;

//int MissingLeaf = 0;
int MissingParticleTrailLength = 0;
int MissingParticleTrail[100];

static int key_fragment ( const int node );
static void add_particle_to_node ( const int ipart, const int node );
static bool particle_is_inside_node ( const peanoKey key, const int lvl,        const int node );
static bool particle_is_inside_node_easy ( const int ipart, const int node );
static void create_node_from_particle ( const int ipart, const int parent,
                                        const peanoKey key, const int lvl, int *NNodes );

void gravity_tree_init();
int Level ( const int node );

int Find_ngb ( const int ipart, const float hsml, int ngblist[NGBMAX] )
{
#ifdef BRUTE_FORCE_NGB
    return Find_ngb_simple ( ipart, hsml, ngblist );
#else
    return Find_ngb_tree ( ipart, hsml, ngblist );
#endif
}

int Find_ngb_tree ( const int ipart, const float hsml, int ngblist[NGBMAX] )
{
    /*int missingNode = P[52975].Tree_Parent;
    if ( ipart == 52931 ) {
        printf ( "Start of ngb tree find, missing particle sits in node %d\n", missingNode );
        printf ( "Testing if this node is a parent of the found leaf.\n" );
        bool found = false;
        for ( int n = missingNode; n < NNodes; ) {
            if ( Tree[n].DNext > 0 ) {
                n += Tree[n].DNext;
            } else {
                ++n;
            }
            if ( n == MissingLeaf ) {
                printf ( "Yes it is\n" );
                found = true;
                break;
            }
        }
        if ( !found ) {
            printf ( "No it is not!!!\n" );
        }
    }*/
    const double boxsize[3] = { Problem.Boxsize[0], Problem.Boxsize[1],
                                Problem.Boxsize[2]
                              };
    const double boxhalf[3] = { boxsize[0] / 2, boxsize[1] / 2, boxsize[2] / 2, };
    const float pos_i[3] = {P[ipart].Pos[0], P[ipart].Pos[1], P[ipart].Pos[2]};

    int node = 1;

    int ngbcnt = 0;

    if ( ipart == 52931 ) {
        for ( int m = 0; m < MissingParticleTrailLength; ++m ) {
            //Copied from tree walk below
            int node = MissingParticleTrail[m];
            double r2 = 0.0;
            for ( int p = 0; p < 3; ++p ) {
                double d = pos_i[p] - Tree[node].Pos[p];
                if ( Problem.Periodic[p] ) {
                    while ( d > boxhalf[p] ) { // find closest image
                        d -= boxsize[p];
                    }
                    while ( d < -boxhalf[p] ) {
                        d += boxsize[p];
                    }
                }
                r2 += d * d;
            }
            float dl = 0.5 * sqrt3 * Tree[node].Size[0];
            printf ( "Node distance check along trail at node %d: (r2 = %g <? p2(dl+hsml) = %g) => %d\tpos_i = (%g, %g, %g) - Tree[node].Pos[p] = (%g, %g, %g) and boxhalf = (%g, %g, %g) and node size = (%g, %g, %g)\n", node, r2,
                     p2 ( dl + hsml ), r2 < p2 ( dl + hsml ), pos_i[0], pos_i[1], pos_i[2], Tree[node].Pos[0], Tree[node].Pos[1], Tree[node].Pos[2], boxhalf[0], boxhalf[1], boxhalf[2], Tree[node].Size[0], Tree[node].Size[1], Tree[node].Size[2] );
        }
    }


    for ( ;; ) {
        /*if ( ipart == 52931 ) {
            printf ( "Start of for loop, node %d\n", node );
            if ( node == missingNode ) {
                printf ( "This is the parent node of the missing particle!\n" );
            }
            if ( node == MissingLeaf ) {
                printf ( "This is the leaf where the missing particle sits\n" );
            }
            for ( int m = 0; m < MissingParticleTrailLength; ++m ) {
                if ( node == MissingParticleTrail[m] ) {
                    printf ( "Encountered a node on the missing particles trail\n" );
                }
            }
        }*/

        double r2 = 0.0;

        for ( int p = 0; p < 3; ++p ) {
            double d = pos_i[p] - Tree[node].Pos[p];

            if ( Problem.Periodic[p] ) {
                while ( d > boxhalf[p] ) { // find closest image
                    d -= boxsize[p];
                }
                while ( d < -boxhalf[p] ) {
                    d += boxsize[p];
                }
            }

            r2 += d * d;
        }

        float dl = 0.5 * sqrt3 * Tree[node].Size[0];

        /*if ( ipart == 52931 ) {
            printf ( "Node distance check: (r2 = %g <? p2(dl+hsml) = %g) => %d\n", r2, p2 ( dl + hsml ), r2 < p2 ( dl + hsml ) );
        }*/

        if ( r2 < p2 ( dl + hsml ) ) {
            /*if ( ipart == 52931 ) {
                printf ( "Look into node with DNext = %d\n", Tree[node].DNext );
            }*/

            if ( Tree[node].DNext < 0 ) {
                /*if ( ipart == 52931 ) {
                    printf ( "Node contains particles, look into it\n" );
                }*/

                int first = - ( Tree[node].DNext + 1 );
                int last = first + Tree[node].Npart;

                /*if ( ipart == 52931 ) {
                    printf ( "Iterate jpart from %d to %d\n", first, last );
                }*/

                for ( int jpart = first; jpart < last; jpart++ ) {

                    r2 = 0.0;

                    for ( int p = 0; p < 3; ++p ) {
                        double d = pos_i[p] - P[jpart].Pos[p];

                        if ( Problem.Periodic[p] ) {
                            while ( d > boxhalf[p] ) { // find closest image
                                d -= boxsize[p];
                            }
                            while ( d < -boxhalf[p] ) {
                                d += boxsize[p];
                            }
                        }

                        r2 += d * d;
                    }

                    /*if ( ipart == 52931 && jpart == 52975 ) {
                        printf ( "Encountered the missing pair (%d -> %d)! Distance^2 = %g; Search radius^2 = %g\n", 52931, 52975, r2, hsml * hsml );
                    }*/

                    if ( r2 < hsml * hsml ) {
                        /*if ( ipart == 52931 && jpart == 52975 ) {
                            printf ( "Accept missing pair!!\n" );
                        }*/
                        ngblist[ngbcnt++] = jpart;
                    }

                    if ( ngbcnt == NGBMAX ) {
                        /*if ( ipart == 52931 ) {
                            printf ( "Maximum ngb reached - aborting\n" );
                        }*/
                        return ngbcnt;
                    }
                }
            } else {
                /*if ( ipart == 52931 ) {
                    printf ( "Is not a leaf\n" );
                }*/
            }
            /*if ( ipart == 52931 ) {
                printf ( "Go to next node\n" );
            }*/

            node++;

            if ( node >= NNodes ) {
                break;
            } else {
                continue;
            }
        } // if dx < dl

        /*if ( ipart == 52931 ) {
            printf ( "End of inner loop, transition to next or further node\n" );
        }*/
        node += max ( 1, Tree[node].DNext );

        if ( node >= NNodes ) {
            /*if ( ipart == 52931 ) {
                printf ( "Reached maximum nodes - break main loop\n" );
            }*/
            break;
        }
    }

    /*if ( ipart == 52931 ) {
        printf ( "ngb tree search done\n" );
    }*/

    return ngbcnt;
}

extern float Guess_hsml ( const size_t ipart, const int DesNumNgb )
{
    int node = P[ipart].Tree_Parent;

    float nodeVol = Tree[node].Size[0] * Tree[node].Size[1] * Tree[node].Size[2];
    float numDens = Tree[node].Npart / nodeVol;

#ifdef TWO_DIM
    return pow ( pi / numDens, 1. / 2. );
#else
    return pow ( fourpithird / numDens, 1. / 3. );
#endif
}

/* We use the Peano Key, whose triplets represents the tree on every level. */

void Build_Tree()
{
    printf ( "Starting tree build\n" );
    MissingParticleTrailLength = 0;
    gravity_tree_init();

    const double boxsize[3] = { Problem.Boxsize[0], Problem.Boxsize[1],
                                Problem.Boxsize[2]
                              };
    const double boxhalf[3] = { boxsize[0] / 2, boxsize[1] / 2, boxsize[2] / 2, };

    create_node_from_particle ( 0, 0, 0, 0, &NNodes ); // first

    Tree[0].Pos[0] = boxhalf[0];
    Tree[0].Pos[1] = boxhalf[1];
    Tree[0].Pos[2] = boxhalf[2];

    float px = P[0].Pos[0] / boxsize[0];
    float py = P[0].Pos[1] / boxsize[1];
    float pz = P[0].Pos[2] / boxsize[2];

    peanoKey last_key = Reversed_Peano_Key ( px, py, pz );

    last_key >>= 3;

    for ( int ipart = 1; ipart < Param.Npart; ipart++ ) {

        double px = P[ipart].Pos[0] / boxsize[0];
        double py = P[ipart].Pos[1] / boxsize[1];
        double pz = P[ipart].Pos[2] / boxsize[2];

        peanoKey key = Reversed_Peano_Key ( px, py, pz );
        if ( ipart == 52975 ) {
            printf ( "Reversed peano key of missing particle:\n" );
            //Print_Int_Bits128(key);
            Printr ( key );
            printf ( "Checking non reversed key:\n" );
            //Print_Int_Bits128(Peano_Key(px, py, pz));
            Print ( Peano_Key ( px, py, pz ) );
        }

        int node = 0;	// current node
        int lvl = 0;	// counts current level
        int parent = 0;	// parent of current node

        while ( lvl < N_PEANO_TRIPLETS ) {

            if ( ipart == 52975 ) {
                printf ( "Checking node %d: Via key: %d; Via coordinate comparison: %d; Node next: %d\n",
                         node, particle_is_inside_node ( key, lvl, node ), particle_is_inside_node_easy ( ipart, node ),
                         Tree[node].DNext );
            }

            // @todo this seems to be buggy!!!
            if ( particle_is_inside_node ( key, lvl, node ) ) { // open
                //if ( particle_is_inside_node_easy ( ipart, node ) ) { // open

                if ( Tree[node].Npart == 1 ) { // refine

                    Tree[node].DNext = 0;

                    create_node_from_particle ( ipart - 1, node, last_key, lvl + 1,
                                                &NNodes ); // new son

                    last_key >>= 3;
                }

                add_particle_to_node ( ipart, node );

                parent = node;

                node++; // decline into node

                lvl++;

                key >>= 3;

            } else { // skip node

                if ( Tree[node].DNext == 0 || node == NNodes - 1 ) {
                    break;    // reached end of branch
                }

                no < de += fmax ( 1, Tree[node].DNext );
            }
        } // while (lvl < 42)

        if ( lvl > N_PEANO_TRIPLETS - 1 ) { // particles closer than PH resolution

            P[ipart].Tree_Parent = parent;

            continue;  // tree cannot be deeper
        }

        if ( Tree[node].DNext == 0 ) { // set DNext for internal node
            if ( node == 0 ) {
                printf ( "Set next for node %d from %d to %d\n", node, Tree[node].DNext, NNodes - node );
            }
            Tree[node].DNext = NNodes - node;    // only delta
        }

        create_node_from_particle ( ipart, parent, key, lvl, &NNodes ); // sibling

        last_key = key >> 3;

    } // for ipart

    Tree[0].DNext = 0;

    int stack[N_PEANO_TRIPLETS + 1] = { 0 };
    int lowest = 0;

    for ( int i = 1; i < NNodes; i++ ) {

        int lvl = Level ( i );

        while ( lvl <= lowest ) { // set pointers

            int node = stack[lowest];

            if ( node > 0 ) {
                Tree[node].DNext = i - node;
            }

            stack[lowest] = 0;

            lowest--;
        }

        if ( Tree[i].DNext == 0 ) { // add node to stack

            stack[lvl] = i;

            lowest = lvl;
        }

    } // for


    //Sanity check
    /*for (int i = 0; i < Param.Npart; ++i) {
        const int node = P[i].Tree_Parent;
        if (Tree[node].DNext >= 0) {
            printf("Particle %d belongs to node %d which is not a leaf but points to node %d\n", i, node, Tree[node].DNext);
        }
    }
    printf ( "End of tree build - looking for leaf in which missing particle 52975 sits\n" );
    for ( int n = 0; n < NNodes; ++n ) {
        if ( Tree[n].DNext < 0 ) {
            const int first = - ( Tree[n].DNext + 1 );
            const int last = first + Tree[n].Npart;

            if ( first <= 52975 && last > 52975 ) {
                printf ( "Found leaf for missing particle: %d\n", n );
                MissingLeaf = n;
            }
        }

    }*/


    return ;
}

/* if a particle falls into a node, its peano triplets down to the level
 * of the node are the same. */

static bool particle_is_inside_node ( const peanoKey key, const int lvl,
                                      const int node )
{
    int part_triplet = key & 0x7;

    int node_triplet = key_fragment ( node );

    return node_triplet == part_triplet;
}

static bool particle_is_inside_node_easy ( const int ipart, const int node )
{
    for ( int d = 0; d < 3; ++d ) {
        const double lower = Tree[node].Pos[d] - 0.5 * Tree[node].Size[d];
        const double upper = Tree[node].Pos[d] + 0.5 * Tree[node].Size[d];
        const double pos = P[ipart].Pos[d];

        if ( lower > pos || pos > upper ) {
            return false;
        }
    }
    return true;
}

/* A new node gets the Peano triplet from the containing particle. */

static void create_node_from_particle ( const int ipart, const int parent,
                                        const peanoKey key, const int lvl, int *NNodes )
{
    const int node = ( *NNodes )++;
    if ( ipart == 52975 ) {
        printf ( "Creating new leaf %d for missing particle\n", node );
    }

    Assert ( *NNodes < Max_Nodes, "Too many nodes \n" );

    Tree[node].DNext = -ipart - 1;

    int keyfragment = ( key & 0x7 ) << 6;

    Tree[node].Bitfield = lvl | keyfragment;

    //! @todo might be problematic in 2D!!!
    const int sign[3] = { -1 + 2 * ( P[ipart].Pos[0] > Tree[parent].Pos[0] ),
                          -1 + 2 * ( P[ipart].Pos[1] > Tree[parent].Pos[1] ),
                          -1 + 2 * ( P[ipart].Pos[2] > Tree[parent].Pos[2] )
                        };

    float size[3] = { 	Problem.Boxsize[0] / ( 1ULL << lvl ),
                        Problem.Boxsize[1] / ( 1ULL << lvl ),
                        Problem.Boxsize[2] / ( 1ULL << lvl )
                    };

    Tree[node].Size[0] = size[0];
    Tree[node].Size[1] = size[1];
    Tree[node].Size[2] = size[2];

    if ( node == 75544 ) {
        printf ( "Node %d Parent %d: Size (%g, %g, %g) sign (%d, %d, %d) Particle pos (%g, %g, %g) Parent pos (%g, %g, %g)\n",
                 node, parent, size[0], size[1], size[2], sign[0], sign[1], sign[2],
                 P[ipart].Pos[0], P[ipart].Pos[1], P[ipart].Pos[2],
                 Tree[parent].Pos[0], Tree[parent].Pos[1], Tree[parent].Pos[2] );
    }

    Tree[node].Pos[0] = Tree[parent].Pos[0] + sign[0] * size[0] * 0.5;
    Tree[node].Pos[1] = Tree[parent].Pos[1] + sign[1] * size[1] * 0.5;
    Tree[node].Pos[2] = Tree[parent].Pos[2] + sign[2] * size[2] * 0.5;

    P[ipart].Tree_Parent = parent;

    if ( node == 75544 ) {
        double x = Tree[node].Pos[0] / Problem.Boxsize[0], y = Tree[node].Pos[1] / Problem.Boxsize[1], z = Tree[node].Pos[2] / Problem.Boxsize[2];
        printf ( "Assign %d%d%d as key fragment to node %d on level %d\nActual peano key at relative coords (%g, %g, %g) is\n",
                 ( ( key & 0x7 ) & 0x4 ) >> 2, ( ( key & 0x7 ) & 0x2 ) >> 1, ( key & 0x7 ) & 0x1, node, lvl, x, y, z );
        peanoKey nodeKey = Peano_Key ( x, y, z );
        Print ( nodeKey );
        printf ( "\n" );
    }

    add_particle_to_node ( ipart, node );

    return ;
}

static  void add_particle_to_node ( const int ipart, const int node )
{
    if ( ipart == 52975 ) {
        printf ( "Adding missing particle to node %d\n", node );
        MissingParticleTrail[MissingParticleTrailLength++] = node;

        for ( int d = 0; d < 3; ++d ) {
            double lower = Tree[node].Pos[d] - 0.5 * Tree[node].Size[d];
            double upper = Tree[node].Pos[d] + 0.5 * Tree[node].Size[d];
            printf ( "Is missing particle in node? Checking dimension %d: (%g <? %g <? %g) => %d\n",
                     d, lower, P[ipart].Pos[d], upper, ( lower <= P[ipart].Pos[d] && P[ipart].Pos[d] <= upper ) );
        }
        double px = P[ipart].Pos[0] / Problem.Boxsize[0];
        double py = P[ipart].Pos[1] / Problem.Boxsize[1];
        double pz = P[ipart].Pos[2] / Problem.Boxsize[2];
        peanoKey key = Reversed_Peano_Key ( px, py, pz );

        int level = Tree[node].Bitfield & 0x3F; //first 6 bits

        int part_triplet = ( key >> 3 * level ) & 0x7;
        int node_triplet = key_fragment ( node );

        printf ( "Checking key segments on level %d: from node %d%d%d and from particle %d%d%d => %d\n",
                 level,
                 ( node_triplet & 0x4 ) >> 2, ( node_triplet & 0x2 ) >> 1, node_triplet & 0x1,
                 ( part_triplet & 0x4 ) >> 2, ( part_triplet & 0x2 ) >> 1, part_triplet & 0x1,
                 node_triplet == part_triplet );

        peanoKey nodeKey = Peano_Key ( Tree[node].Pos[0] / Problem.Boxsize[0], Tree[node].Pos[1] / Problem.Boxsize[1],
                                       Tree[node].Pos[2] / Problem.Boxsize[2] );
        printf ( "Full node peano key:\n" );
        //Print_Int_Bits128(nodeKey);
        Print ( nodeKey );
        printf ( "\n" );
    }
    Tree[node].Npart++;

    return ;
}

static int key_fragment ( const int node )
{
    const uint32_t bitmask = 7UL << 6;

    return ( Tree[node].Bitfield & bitmask ) >> 6; // return bit 6-8
}

int Level ( const int node )
{
    return Tree[node].Bitfield & 0x3FUL; // return but 0-5
}


void gravity_tree_init()
{
    Max_Nodes = Param.Npart * NODES_PER_PARTICLE;

    size_t nBytes = Max_Nodes * sizeof ( *Tree );

    if ( Tree == NULL ) {
        Tree = malloc ( nBytes );
    }

    memset ( Tree, 0, nBytes );

    //    printf("\n   Tree Build istart=%d npart=%d maxnodes=%d ",
    //            Omp.ThreadID, Param.Npart[0], Max_Nodes);
    //    fflush(stdout);

    NNodes = 0;

    return ;
}

int Find_ngb_simple ( const int ipart, const float hsml, int ngblist[NGBMAX] )
{
    const double boxsize[3] = { Problem.Boxsize[0], Problem.Boxsize[1],
                                Problem.Boxsize[2]
                              };
    const double boxhalf[3] = { boxsize[0] / 2, boxsize[1] / 2, boxsize[2] / 2 };

    memset ( ngblist, 0, NGBMAX * sizeof ( *ngblist ) );

    int ngbcnt = 0;

    for ( int jpart = 0; jpart < Param.Npart; jpart++ ) {

        double r2 = 0.0;

        for ( int p = 0; p < 3; ++p ) {
            double d = P[ipart].Pos[p] - P[jpart].Pos[p];

            if ( Problem.Periodic[p] ) {
                while ( d > boxhalf[p] ) { // find closest image
                    d -= boxsize[p];
                }
                while ( d < -boxhalf[p] ) {
                    d += boxsize[p];
                }
            }

            r2 += d * d;
        }
        /*if ( ipart == 52931 && jpart == 52975 ) {
            printf ( "Encountered missing particle pair in brute force manner. (r2 = %g <? hsml^2 = %g) => %d\n", r2,
                     hsml * hsml, r2 < hsml * hsml );
        }*/

        if ( r2 < hsml * hsml ) {
            ngblist[ngbcnt++] = jpart;
        }

        if ( ngbcnt == NGBMAX ) {

            //printf("WARNING, ngbcnt == %d, increase NGBMAX ! ",
            //      ngbcnt);

            break;
        }
    }

    return ngbcnt ;
}
