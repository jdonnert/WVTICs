extern void Build_Tree();
extern int Find_ngb_tree ( const size_t, const float, int * );
extern int *Find_ngb_tree_recursive ( size_t, float, int );
int Find_ngb_simple ( const int ipart,  const float hsml, int *ngblist );
extern float Guess_hsml ( const size_t ipart, const int DesNumNgb );
int Ngbcnt ;
int Ngblist[NGBMAX];
