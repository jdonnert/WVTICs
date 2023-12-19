#ifndef TREE_H
#define TREE_H

void Build_Tree();
int Find_ngb ( const int ipart, const float hsml, int ngblist[NGBMAX] );
int Find_ngb_tree ( const int ipart, const float hsml, int ngblist[NGBMAX] );
int Find_ngb_simple ( const int ipart, const float hsml, int ngblist[NGBMAX] );
float Guess_hsml ( const size_t ipart, const int DesNumNgb );

#endif
