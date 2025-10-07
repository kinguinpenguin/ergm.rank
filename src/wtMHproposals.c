/*  File src/wtMHproposals.c in package ergm.rank, part of the Statnet suite of
 *  packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2008-2025 Statnet Commons
 */

#include "ergm_wtMHproposal.h"
#include "ergm_MHstorage.h"
#include "wtchangestats_rank_aux.h"
#include "ergm_dyadgen.h"

#define ENSURE_CHANGED(gen) {double to; do{ to = gen; }while(to == from); return to;}

static double rdunifj(double *param, double from, double fudge){
  ENSURE_CHANGED(floor(runif(param[0], param[1]+1)));
}
static double lddunifj(double *param, double from, double to, double fudge){
  return(-log(param[1]-param[0]));
}
static double lhdunif(double *param, double to){
  return(0);
}

/*********************
 void MH_AlterSwap

 Default MH algorithm for ERGM over complete orderings
*********************/
WtMH_P_FN(MH_AlterSwap){  
  Vertex tail, head1, head2;
  
  if(MHp->ntoggles == 0) { // Initialize AlterSwap 
    MHp->ntoggles=2;
    return;
  }

  unsigned int trytoggle;

  for(trytoggle=0; trytoggle<MAX_TRIES; trytoggle++){
  
    GetRandDyad(&tail, &head1, nwp);
    
    if(nwp->bipartite){
      head2 = 1 + nwp->bipartite + unif_rand() * (nwp->nnodes - nwp->bipartite - 1);
      if(head2 >= head1) head2++;
    }else{
      head2 = 1 + unif_rand() * (nwp->nnodes - 2);
      if(head2 >= tail) head2++;
      if(head2 >= head1){
	head2++;
	if(head2 == tail) head2++; // Increment if landed on top of tail. Note that head2 must have started out < tail, so there is no possibility of head2 being incremented 3 times.
      }
    }
    
    Mtail[0] = Mtail[1] = tail;
    Mhead[0] = head1;
    Mhead[1] = head2;
    
    Mweight[1] = WtGetEdge(Mtail[0],Mhead[0],nwp);
    Mweight[0] = WtGetEdge(Mtail[1],Mhead[1],nwp);

    if(Mweight[0]!=Mweight[1]) break;
  }
  if(trytoggle>=MAX_TRIES){
    MHp->toggletail[0]=MH_FAILED;
    MHp->togglehead[0]=MH_UNSUCCESSFUL;	
  }
}

/*********************
 void MH_AlterSwapPartial

 Default MH algorithm for ERGM over partial orderings
*********************/

WtMH_I_FN(Mi_Disc){
  MH_STORAGE = DyadGenInitializeR(MHp->R, nwp, FALSE);
  MHp->ntoggles = ((DyadGen *) MH_STORAGE)->ndyads!=0 ? 1 : MH_FAILED;
}

WtMH_P_FN(Mp_Disc){
  const double fudge = 0.5; // Mostly comes in when proposing from 0.
  DyadGen *gen = (DyadGen *) MH_STORAGE;
  double (*rj)(double *param, double from, double fudge),
  (*ldj)(double *param, double from, double to, double fudge),
  (*lh)(double *param, double to);
  rj = rdunifj; 
  ldj = lddunifj; 
  lh = lhdunif;
}

/*********************
 void MH_AdjacentAlterSwap

 MH algorithm for ERGMs over complete orderings that selects an ego
 and an alter and promotes the alter up; assumes that the ordering is
 complete.
*********************/
WtMH_P_FN(MH_AdjacentAlterSwap){  
  GET_AUX_STORAGE(0, double *, sm);
  GET_AUX_STORAGE(1, Pair *, udsm);

  if(MHp->ntoggles == 0) { // Initialize AlterSwap 
    MHp->ntoggles=2;
    return;
  }

  Vertex tail, head1;
  GetRandDyad(&tail, &head1, nwp);

  Vertex head2 = udsm[tail][head1].up;

  if(head2 == 0) {
    MHp->toggletail[0] = MH_FAILED;
    MHp->togglehead[0] = MH_CONSTRAINT;
    return;
  }

  Mtail[0] = Mtail[1] = tail;
  Mhead[0] = head1;
  Mhead[1] = head2;

  Mweight[1] = sm[Mtail[0]][Mhead[0]];
  Mweight[0] = sm[Mtail[1]][Mhead[1]];
}

/*********************
 void MH_AdjacentAlterSwapPartial

 MH algorithm for ERGMs over partial orderings that selects an ego
 and an alter and promotes the alter up
*********************/
