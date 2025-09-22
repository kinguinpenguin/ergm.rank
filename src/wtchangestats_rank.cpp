/*  File src/wtchangestats_rank.c in package ergm.rank, part of the Statnet
 *  suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2008-2025 Statnet Commons
 */

extern "C" {
#include "wtchangestat_rank.h"
}

#include "wtchangestats_rank_aux.h"

extern "C" {

WtC_CHANGESTAT_FN(c_edgecov_rank) {
  GET_AUX_STORAGE(0, double *, sm);
  GET_AUX_STORAGE(1, Pair *, udsm);
  Vertex v1=tail;
  Vertex v2=head;
  double v12_old = sm[tail][head];
  double v12_new = weight;

  if (v12_new > v12_old) { // New is above, so iterate upwards
    for (Vertex v3 : UpDownRange(v1, v2, sm, udsm, v12_old, v12_new)) {
      double v123_covdiff=INPUT_PARAM[(v1-1)*N_NODES + (v2-1)] - INPUT_PARAM[(v1-1)*N_NODES + (v3-1)];
      if(v123_covdiff == 0) continue;
      double v13_old = sm[v1][v3];
      if (v12_old < v13_old) CHANGE_STAT[0] += v123_covdiff;
      if (v12_new > v13_old) CHANGE_STAT[0] += v123_covdiff;
    }
  } else { // New is below, so iterate downwards
    for (Vertex v3 : UpDownRange(v1, v2, sm, udsm, v12_old, v12_new)) {
      double v123_covdiff=INPUT_PARAM[(v1-1)*N_NODES + (v2-1)] - INPUT_PARAM[(v1-1)*N_NODES + (v3-1)];
      if(v123_covdiff == 0) continue;
      double v13_old = sm[v1][v3];
      if (v12_old > v13_old) CHANGE_STAT[0] -= v123_covdiff;
      if (v12_new < v13_old) CHANGE_STAT[0] -= v123_covdiff;
    }
  }
  /*GET_AUX_STORAGE(0, double *, sm);
    GET_AUX_STORAGE(Pair *, udsm);
    Vertex v1=tail;
    Vertex v2=head;
    double v12_old = sm[tail][head];
    double v12_new = weight;
    for (Vertex v3=1; v3 <= N_NODES; v3++){
    if(v3==v2 || v3==v1) continue;
    double v123_covdiff=INPUT_PARAM[(v1-1)*N_NODES + (v2-1)] - INPUT_PARAM[(v1-1)*N_NODES + (v3-1)];
    if(v123_covdiff==0) continue; // If covariate value is 0, don't bother looking up the ranking of v3 by v1.
    double v13_old=sm[v1][v3];
    if(v12_old>v13_old)
    CHANGE_STAT[0] -= v123_covdiff;
    if(v12_old<v13_old)
    CHANGE_STAT[0] += v123_covdiff;
    if(v12_new>v13_old)
    CHANGE_STAT[0] += v123_covdiff;
    if(v12_new<v13_old)
    CHANGE_STAT[0] -= v123_covdiff;
    }*/
}

WtS_CHANGESTAT_FN(s_edgecov_rank){
  GET_AUX_STORAGE(void, sm_raw);
  double **sm = (double **)sm_raw;
  for (Vertex v1=1; v1 <= N_NODES; v1++){
    for (Vertex v2=1; v2 <= N_NODES; v2++){
      if(v2==v1) continue;
      double v12=sm[v1][v2];
      for (Vertex v3=1; v3 <= N_NODES; v3++){
	if(v3==v2 || v3==v1) continue;
	double v123_covdiff = INPUT_PARAM[(v1-1)*N_NODES + (v2-1)] - INPUT_PARAM[(v1-1)*N_NODES + (v3-1)];
	if(v123_covdiff!=0 && v12>sm[v1][v3]) // Short-circuit the lookup of ranking of v3 by v1 if covariate is 0.
	  CHANGE_STAT[0] += v123_covdiff;
      }
    }
  }
}


WtC_CHANGESTAT_FN(c_inconsistency_rank){
  GET_AUX_STORAGE(0, double *, sm);
  GET_AUX_STORAGE(1, Pair *, udsm);
  Vertex v1=tail;
  Vertex v2=head;
  double v12_old = sm[tail][head];
  double v12_ref = INPUT_PARAM[(v1-1)*N_NODES+(v2-1)];
  double v12_new = weight;

  for (Vertex v3 : UpDownRange(v1, v2, sm, udsm, v12_old, v12_new)) {
    double v13 = sm[v1][v3];
    double v13_ref=INPUT_PARAM[(v1-1)*N_NODES+(v3-1)];
    if((v12_old>v13)!=(v12_ref>v13_ref)) CHANGE_STAT[0]--;
    if((v12_new>v13)!=(v12_ref>v13_ref)) CHANGE_STAT[0]++;
    if((v13>v12_old)!=(v13_ref>v12_ref)) CHANGE_STAT[0]--;
    if((v13>v12_new)!=(v13_ref>v12_ref)) CHANGE_STAT[0]++;
  }
}

WtS_CHANGESTAT_FN(s_inconsistency_rank){
  GET_AUX_STORAGE(void, sm_raw);
  double **sm = (double **)sm_raw;
  for(Vertex v1=1; v1 <= N_NODES; v1++){
    for(Vertex v2=1; v2 <= N_NODES; v2++){
      if(v2==v1) continue;
      double v12 = sm[v1][v2], v12_ref = INPUT_PARAM[(v1-1)*N_NODES+(v2-1)];
      for(Vertex v3=1; v3 <= N_NODES; v3++){
	if(v3==v2 || v3==v1) continue;
	unsigned int
	  v123 = v12>sm[v1][v3],
	  v123_ref = v12_ref>INPUT_PARAM[(v1-1)*N_NODES+(v3-1)];
	if(v123!=v123_ref) CHANGE_STAT[0]++;
      }
    }
  }
}

WtC_CHANGESTAT_FN(c_inconsistency_cov_rank){
  GET_AUX_STORAGE(0, double *, sm);
  GET_AUX_STORAGE(1, Pair *, udsm);
  unsigned int cov_start = N_NODES*N_NODES;
  Vertex v1=tail;
  Vertex v2=head;
  double v12_ref = INPUT_PARAM[(v1-1)*N_NODES+(v2-1)];
  double v12_old = sm[tail][head];
  double v12_new = weight;
  for (Vertex v3 : UpDownRange(v1, v2, sm, udsm, v12_old, v12_new)) {
    double v123_cov = INPUT_PARAM[cov_start + (v1-1)*N_NODES*N_NODES + (v2-1)*N_NODES + (v3-1)];
    double v132_cov = INPUT_PARAM[cov_start + (v1-1)*N_NODES*N_NODES + (v3-1)*N_NODES + (v2-1)];
    if(v123_cov!=0 || v123_cov!=0){
      double v13=sm[v1][v3];
      double v13_ref=INPUT_PARAM[(v1-1)*N_NODES+(v3-1)];
      if(v123_cov!=0){
        if((v12_old>v13)!=(v12_ref>v13_ref)) CHANGE_STAT[0]-=v123_cov;
        if((v12_new>v13)!=(v12_ref>v13_ref)) CHANGE_STAT[0]+=v123_cov;
      }
      if(v132_cov!=0){
        if((v13>v12_old)!=(v13_ref>v12_ref)) CHANGE_STAT[0]-=v132_cov;
        if((v13>v12_new)!=(v13_ref>v12_ref)) CHANGE_STAT[0]+=v132_cov;
      }
    }
  }
  /*
      for(Vertex v3=1; v3 <= N_NODES; v3++){
	if(v3==v2 || v3==v1) continue;
	double v123_cov = INPUT_PARAM[cov_start + (v1-1)*N_NODES*N_NODES + (v2-1)*N_NODES + (v3-1)];
	double v132_cov = INPUT_PARAM[cov_start + (v1-1)*N_NODES*N_NODES + (v3-1)*N_NODES + (v2-1)];
	if(v123_cov!=0 || v123_cov!=0){
	  double v13=sm[v1][v3];
	  double v13_ref=INPUT_PARAM[(v1-1)*N_NODES+(v3-1)];

	  if(v123_cov!=0){
	    if((v12_old>v13)!=(v12_ref>v13_ref)) CHANGE_STAT[0]-=v123_cov;
	    if((v12_new>v13)!=(v12_ref>v13_ref)) CHANGE_STAT[0]+=v123_cov;
	  }

	  if(v132_cov!=0){
	    if((v13>v12_old)!=(v13_ref>v12_ref)) CHANGE_STAT[0]-=v132_cov;
	    if((v13>v12_new)!=(v13_ref>v12_ref)) CHANGE_STAT[0]+=v132_cov;
	  }
	}
      }
  */
}

WtS_CHANGESTAT_FN(s_inconsistency_cov_rank){
  GET_AUX_STORAGE(void, sm_raw);
  double **sm = (double **)sm_raw;
  unsigned int cov_start = N_NODES*N_NODES;
  for(Vertex v1=1; v1 <= N_NODES; v1++){
    for(Vertex v2=1; v2 <= N_NODES; v2++){
      if(v2==v1) continue;
      double v12 = sm[v1][v2], v12_ref = INPUT_PARAM[(v1-1)*N_NODES+(v2-1)];
      for(Vertex v3=1; v3 <= N_NODES; v3++){
	if(v3==v2 || v3==v1) continue;
	double v123_cov = INPUT_PARAM[cov_start + (v1-1)*N_NODES*N_NODES + (v2-1)*N_NODES + (v3-1)];
	if(v123_cov==0) continue;
	unsigned int
	  v123 = v12>sm[v1][v3],
	  v123_ref = v12_ref>INPUT_PARAM[(v1-1)*N_NODES+(v3-1)];
	if(v123!=v123_ref)
	  CHANGE_STAT[0]+=v123_cov;
      }
    }
  }
}

// Keep track of which have already been triggered, so that we only print once.
bool up_tail_head_k_old = false,
  up_tail_k_head_old = false,
  up_head_k_tail_old = false,
  up_k_head_tail_old = false,
  down_tail_head_k_old = false,
  down_tail_k_head_old = false,
  down_head_k_tail_old = false,
  down_k_head_tail_old = false,
  up_tail_head_k_new = false,
  up_tail_k_head_new = false,
  up_head_k_tail_new = false,
  up_k_head_tail_new = false,
  down_tail_head_k_new = false,
  down_tail_k_head_new = false,
  down_head_k_tail_new = false,
  down_k_head_tail_new = false;

// See, e.g., https://gcc.gnu.org/onlinedocs/cpp/Stringizing.html and
// https://gcc.gnu.org/onlinedocs/cpp/Concatenation.html for
// documentation: these basically construct the above variable names
// (using ##) and their corresponding strings (using #) and call
// Rprintf() the first time each of these if() statements triggers.
//
// After a sufficiently long MCMC run under a proposal that permits
// ties, we can safely assume that if an if() statement hasn't
// triggered, it never will.
#define defer_perm(i, j, l, dir)                                        \
  if (sm[l][j] > sm[l][i] && sm[i][l] > sm[i][j]) { CHANGE_STAT[0]--; }\
  if (GETNEWWTSM(l, j) > GETNEWWTSM(l, i) && GETNEWWTSM(i, l) > GETNEWWTSM(i, j)) { CHANGE_STAT[0]++;}

WtC_CHANGESTAT_FN(c_deference){
  GET_AUX_STORAGE(0, double *, sm);
  GET_AUX_STORAGE(1, Pair *, udsm);
  Vertex vth_old = sm[tail][head];
  Vertex vth_new = weight;
  if (vth_new > vth_old) { // New is above, so iterate upwards
    for (Vertex k : UpDownRange(tail, head, sm, udsm, vth_old, vth_new)) {
      // i or l can be tail
      defer_perm(tail, head, k, up);
      defer_perm(tail, k, head, up);
      defer_perm(head, k, tail, up);
      defer_perm(k, head, tail, up);
    }
  } else { // New is below, so iterate downwards
    for (Vertex k : UpDownRange(tail, head, sm, udsm, vth_old, vth_new)) {
      defer_perm(tail, head, k, down);
      defer_perm(tail, k, head, down);
      defer_perm(head, k, tail, down);
      defer_perm(k, head, tail, down);
    }
  }
  /*
      for(Vertex v1=1; v1 <= N_NODES; v1++){
	for(Vertex v3=1; v3 <= N_NODES; v3++){
	  if(v3==v1) continue;
	  double v31_old = sm[v3][v1];
	  double v13_old = sm[v1][v3];
	  double v31_new = GETNEWWTOLD_M(v3,v1,v31_old);
	  double v13_new = GETNEWWTOLD_M(v1,v3,v13_old);

	  for(Vertex v2=1; v2 <= N_NODES; v2++){
	    if(v2==v1 || v3==v2 ||
	       (tail!=v1 && tail!=v3 &&
		head!=v1 && head!=v2 && head!=v3)) continue;
	    double v32_old = sm[v3][v2];
	    double v12_old = sm[v1][v2];
	    double v32_new = GETNEWWTOLD_M(v3,v2,v32_old);
	    double v12_new = GETNEWWTOLD_M(v1,v2,v12_old);
	    if(v32_old>v31_old && v13_old>v12_old)
	      CHANGE_STAT[0]--;
	    if(v32_new>v31_new && v13_new>v12_new)
	      CHANGE_STAT[0]++;
	  }
	}
      }*/
}


WtS_CHANGESTAT_FN(s_deference){
  GET_AUX_STORAGE(0, double *, sm);
  for(Vertex v1=1; v1 <= N_NODES; v1++){
    for(Vertex v3=1; v3 <= N_NODES; v3++){
      if(v3==v1) continue;
      double v31 = sm[v3][v1], v13 = sm[v1][v3];
      for(Vertex v2=1; v2 <= N_NODES; v2++){
	if(v2==v1 || v3==v2) continue;
	if(sm[v3][v2]>v31 && v13>sm[v1][v2])
	  CHANGE_STAT[0]++;
      }
    }
  }
}

WtC_CHANGESTAT_FN(c_nodeicov_rank){
  GET_AUX_STORAGE(void, sm_raw);
  double **sm = (double **)sm_raw;
  unsigned int oshift = N_INPUT_PARAMS / N_CHANGE_STATS;
  Vertex v1=tail;
  Vertex v2=head;
  double v12_old = sm[tail][head];
  double v12_new = weight;
  for (Vertex v3=1; v3 <= N_NODES; v3++){
    if(v3==v2 || v3==v1) continue;
    double v13_old=sm[v1][v3];
    for(unsigned int j=0, o=0; j<N_CHANGE_STATS; j++, o+=oshift){
      double v23_covdiff = INPUT_PARAM[v2+o-1] - INPUT_PARAM[v3+o-1];
      if(v12_old>v13_old)
        CHANGE_STAT[j] -= v23_covdiff;
      if(v12_old<v13_old)
        CHANGE_STAT[j] += v23_covdiff;
      if(v12_new>v13_old)
        CHANGE_STAT[j] += v23_covdiff;
      if(v12_new<v13_old)
        CHANGE_STAT[j] -= v23_covdiff;
    }
  }
}

WtS_CHANGESTAT_FN(s_nodeicov_rank){
  GET_AUX_STORAGE(void, sm_raw);
  double **sm = (double **)sm_raw;
  unsigned int oshift = N_INPUT_PARAMS / N_CHANGE_STATS;
  for (Vertex v1=1; v1 <= N_NODES; v1++){
    for (Vertex v2=1; v2 <= N_NODES; v2++){
      if(v2==v1) continue;
      double v12=sm[v1][v2];
      for (Vertex v3=1; v3 <= N_NODES; v3++){
	if(v3==v2 || v3==v1) continue;
	if(v12>sm[v1][v3]){
	  for(unsigned int j=0, o=0; j<N_CHANGE_STATS; j++, o+=oshift){
	    CHANGE_STAT[j] += INPUT_PARAM[v2+o-1] - INPUT_PARAM[v3+o-1];
	  }
	}
      }
    }
  }
}

#define nonconform_perm(i, j, k, l, dir)\
  if (sm[l][j] > sm[l][k] && !(sm[i][j] > sm[i][k])) { CHANGE_STAT[0]--; }\
  if (GETNEWWTSM(l, j) > GETNEWWTSM(l, k) && !(GETNEWWTSM(i,j) > GETNEWWTSM(i, k))) { CHANGE_STAT[0]++;}

WtC_CHANGESTAT_FN(c_nonconformity) {
  GET_AUX_STORAGE(0, double *, sm);
  GET_AUX_STORAGE(1, Pair *, udsm);
  Vertex vth_old = sm[tail][head];
  Vertex vth_new = weight;
  if (vth_new > vth_old) { // New is above, so iterate upwards
      for (Vertex k : UpDownRange(tail, head, sm, udsm, vth_old, vth_new)) {
        for(Vertex l=1; l <= N_NODES; l++) {
        // i or l can be tail
        if (l != k && l != tail && l != head) {
          nonconform_perm(tail, head, k, l, up);
          nonconform_perm(tail, k, head, l, up);
          nonconform_perm(l, k, head, tail, up);
          nonconform_perm(l, head, k, tail, up);
        }
      }
    }
  } else { // New is below, so iterate downwards
      for (Vertex k : UpDownRange(tail, head, sm, udsm, vth_old, vth_new)) {
        for(Vertex l=1; l <= N_NODES; l++) {
        if (l != k && l != tail && l != head) {
          nonconform_perm(tail, head, k, l, down);
          nonconform_perm(tail, k, head, l, down);
          nonconform_perm(l, k, head, tail, down);
          nonconform_perm(l, head, k, tail, down);
        }
      }
    }
  }
  /*
  // GET_AUX_STORAGE(0, double *, sm);
      Vertex v1=tail;

      for(Vertex v2=1; v2 <= N_NODES; v2++){
	if(v2==v1) continue;

	for(Vertex v3=1; v3 <= N_NODES; v3++){
	  if(v3==v2 || v3==v1) continue;

	  double v13_old=sm[v1][v3];
	  double v23=sm[v2][v3];
	  double v13_new = GETNEWWTOLD_M(v1,v3,v13_old);

	  for(Vertex v4=1; v4 <= N_NODES; v4++){
	    if(v4==v3 || v4==v2 || v4==v1 ||
	       (head!=v3 && head!=v4)) continue;

	    double v14_old=sm[v1][v4];
	    double v24=sm[v2][v4];
	    double v14_new = GETNEWWTOLD_M(v1,v4,v14_old);

	    if((v13_old>v14_old)!=(v23>v24)) CHANGE_STAT[0]--;
	    if((v13_new>v14_new)!=(v23>v24)) CHANGE_STAT[0]++;
	  }
	}
      }*/
}

WtS_CHANGESTAT_FN(s_nonconformity){
  GET_AUX_STORAGE(0, double *, sm);
  for(Vertex v1=1; v1 <= N_NODES; v1++){
    for(Vertex v2=1; v2 < v1; v2++){
      for(Vertex v3=1; v3 <= N_NODES; v3++){
	if(v3==v2 || v3==v1) continue;
	double v13=sm[v1][v3], v23=sm[v2][v3];
	for(Vertex v4=1; v4 <= N_NODES; v4++){
	  if(v4==v3 || v4==v2 || v4==v1) continue;
	  unsigned int
	    v134 = v13>sm[v1][v4],
	    v234 = v23>sm[v2][v4];
	  if(v134!=v234) CHANGE_STAT[0]++;
	}
      }
    }
  }
}

// From Krivitsky and Butts paper, here, v1=i, v2=j, v3=l, v4=k.
WtC_CHANGESTAT_FN(c_local1_nonconformity){
  GET_AUX_STORAGE(0, double *, sm);
      Vertex v1=tail;

      for(Vertex v2=1; v2 <= N_NODES; v2++){
	if(v2==v1) continue;
	double v12_old=sm[v1][v2];
	double v12_new=GETNEWWTOLD_M(v1,v2,v12_old);

	for(Vertex v3=1; v3 <= N_NODES; v3++){
	  if(v3==v2 || v3==v1) continue;
	  double v13_old=sm[v1][v3];
	  double v13_new=GETNEWWTOLD_M(v1,v3,v13_old);

	  if(v13_old<=v12_old && v13_new<=v12_new) continue;

	  double v32_old=sm[v3][v2];
	  double v32_new=GETNEWWTOLD_M(v3,v2,v32_old);

	  for(Vertex v4=1; v4 <= N_NODES; v4++){
	    if(v4==v3 || v4==v2 || v4==v1 ||
	       (head!=v4 && head!=v3 && head!=v2)) continue;

	    double v14_old=sm[v1][v4];
	    double v14_new=GETNEWWTOLD_M(v1,v4,v14_old);

	    double v34_old=sm[v3][v4];
	    double v34_new=GETNEWWTOLD_M(v3,v4,v34_old);

	    if(v13_old>v12_old && v12_old<=v14_old && v32_old>v34_old) CHANGE_STAT[0]--;
	    if(v13_new>v12_new && v12_new<=v14_new && v32_new>v34_new) CHANGE_STAT[0]++;
	  }
	}
      }

      Vertex v3=tail;

      for(Vertex v1=1; v1 <= N_NODES; v1++){
	if(v1==v3) continue;
	double v13_old=sm[v1][v3];
	double v13_new=GETNEWWTOLD_M(v1,v3,v13_old);

	for(Vertex v2=1; v2 <= N_NODES; v2++){
	  if(v2==v3 || v2==v1) continue;
	    double v12_old=sm[v1][v2];
	    double v12_new=GETNEWWTOLD_M(v1,v2,v12_old);

	    if(v13_old<=v12_old && v13_new<=v12_new) continue;

	    double v32_old=sm[v3][v2];
	    double v32_new=GETNEWWTOLD_M(v3,v2,v32_old);

	    for(Vertex v4=1; v4 <= N_NODES; v4++){
	      if(v4==v3 || v4==v2 || v4==v1 ||
		 (head!=v4 && head!=v2 && head!=v1)) continue;

	      double v14_old=sm[v1][v4];
	      double v14_new=GETNEWWTOLD_M(v1,v4,v14_old);

	      double v34_old=sm[v3][v4];
	      double v34_new=GETNEWWTOLD_M(v3,v4,v34_old);

	      if(v13_old>v12_old && v12_old<=v14_old && v32_old>v34_old) CHANGE_STAT[0]--;
	      if(v13_new>v12_new && v12_new<=v14_new && v32_new>v34_new) CHANGE_STAT[0]++;
	    }
	}
      }
}


// From Krivitsky and Butts paper, here, v1=i, v2=j, v3=l, v4=k.
WtS_CHANGESTAT_FN(s_local1_nonconformity){
  GET_AUX_STORAGE(0, double *, sm);
  for(Vertex v1=1; v1 <= N_NODES; v1++){
    for(Vertex v2=1; v2 <= N_NODES; v2++){
      if(v2==v1) continue;
      double v12=sm[v1][v2];
      for(Vertex v3=1; v3 <= N_NODES; v3++){
	if(v3==v2 || v3==v1) continue;
	double v13=sm[v1][v3];
	if(v13<=v12) continue;
	double v32=sm[v3][v2];
	for(Vertex v4=1; v4 <= N_NODES; v4++){
	  if(v4==v3 || v4==v2 || v4==v1) continue;
	  double v14=sm[v1][v4];
	  double v34=sm[v3][v4];
	  if(v32>v34 && v12<=v14) CHANGE_STAT[0]++;
	}
      }
    }
  }
}


WtC_CHANGESTAT_FN(c_local2_nonconformity){
  GET_AUX_STORAGE(0, double *, sm);
      Vertex v1=tail;
      for(Vertex v2=1; v2 <= N_NODES; v2++){
	if(v2==v1) continue;
	double v12_old=sm[v1][v2];
	double v12_new=GETNEWWTOLD_M(v1,v2,v12_old);
	for(Vertex v3=1; v3 <= N_NODES; v3++){
	  if(v3==v2 || v3==v1) continue;
	  double v13_old=sm[v1][v3];
	  double v13_new=GETNEWWTOLD_M(v1,v3,v13_old);

	  if(v13_old<=v12_old && v13_new<=v12_new) continue;

	  double v32_old=sm[v3][v2];
	  double v32_new=GETNEWWTOLD_M(v3,v2,v32_old);

	  for(Vertex v4=1; v4 <= N_NODES; v4++){
	    if(v4==v3 || v4==v2 || v4==v1 ||
	       (head!=v2 && head!=v3 && head!=v4)) continue;
	    double v14_old=sm[v1][v4];
	    double v14_new=GETNEWWTOLD_M(v1,v4,v14_old);
	    double v34_old=sm[v3][v4];
	    double v34_new=GETNEWWTOLD_M(v3,v4,v34_old);

	    if(v13_old>v12_old && v14_old<=v12_old && v34_old>v32_old) CHANGE_STAT[0]--;
	    if(v13_new>v12_new && v14_new<=v12_new && v34_new>v32_new) CHANGE_STAT[0]++;
	  }
	}
      }

      for(Vertex v1=1; v1 <= N_NODES; v1++){
	for(Vertex v2=1; v2 <= N_NODES; v2++){
	  if(v2==v1) continue;
	  double v12_old=sm[v1][v2];
	  double v12_new=GETNEWWTOLD_M(v1,v2,v12_old);
	  Vertex v3=tail;
	  if(v3==v2 || v3==v1) continue;
	  double v13_old=sm[v1][v3];
	  double v13_new=GETNEWWTOLD_M(v1,v3,v13_old);

	  if(v13_old<=v12_old && v13_new<=v12_new) continue;

	  double v32_old=sm[v3][v2];
	  double v32_new=GETNEWWTOLD_M(v3,v2,v32_old);

	  for(Vertex v4=1; v4 <= N_NODES; v4++){
	    if(v4==v3 || v4==v2 || v4==v1 ||
	       (head!=v2 && head!=v3 && head!=v4)) continue;
	    double v14_old=sm[v1][v4];
	    double v14_new=GETNEWWTOLD_M(v1,v4,v14_old);
	    double v34_old=sm[v3][v4];
	    double v34_new=GETNEWWTOLD_M(v3,v4,v34_old);

	    if(v13_old>v12_old && v14_old<=v12_old && v34_old>v32_old) CHANGE_STAT[0]--;
	    if(v13_new>v12_new && v14_new<=v12_new && v34_new>v32_new) CHANGE_STAT[0]++;
	  }
	}
      }
}

WtS_CHANGESTAT_FN(s_local2_nonconformity){
  GET_AUX_STORAGE(0, double *, sm);
  for(Vertex v1=1; v1 <= N_NODES; v1++){
    for(Vertex v2=1; v2 <= N_NODES; v2++){
      if(v2==v1) continue;
      double v12=sm[v1][v2];
      for(Vertex v3=1; v3 <= N_NODES; v3++){
	if(v3==v2 || v3==v1) continue;
	double v13=sm[v1][v3];
	if(v13<=v12) continue;
	double v32=sm[v3][v2];
	for(Vertex v4=1; v4 <= N_NODES; v4++){
	  if(v4==v3 || v4==v2 || v4==v1) continue;
	  double v14=sm[v1][v4];
	  double v34=sm[v3][v4];
	  if(v14<=v12 && v34>v32) CHANGE_STAT[0]++;
	}
      }
    }
  }
}

#define nonconform_local_global_perm(i, j, k, l)                      \
  if ((sm[i][l] > sm[i][j]) && (sm[i][l] > sm[i][k]) && (sm[l][j] > sm[l][k]) && !(sm[i][j] > sm[i][k])) { CHANGE_STAT[0]--; Rprintf("Triggered: " # i "=%u " # j "=%u " # k "=%u " # l "=%u old\n", i, j, k, l);}\
  if ((GETNEWWTSM(i, l) > GETNEWWTSM(i,j)) && (GETNEWWTSM(i,l) > GETNEWWTSM(i, k)) && (GETNEWWTSM(l, j) > GETNEWWTSM(l, k)) && !(GETNEWWTSM(i, j) > GETNEWWTSM(i, k))) { CHANGE_STAT[0]++; Rprintf("Triggered: i=" # i "=%u j=" # j "=%u k=" # k "=%u l=" # l "=%u new\n", i, j, k, l);}

// From Krivitsky and Butts paper, here, v1=i, v2=j, v3=k, v4=l.
WtC_CHANGESTAT_FN(c_localAND_nonconformity){
  GET_AUX_STORAGE(0, double *, sm);
  GET_AUX_STORAGE(1, Pair *, udsm);
  Vertex vth_old = sm[tail][head];
  Vertex vth_new = weight;
  for (Vertex v3 : UpDownRange(tail, head, sm, udsm, vth_old, vth_new)) {
    for(Vertex v4 = 1; v4 <= N_NODES; v4++) {
      // i or l can be tail 
      if (v4 != v3 && v4 != tail && v4 != head) {
        // tail, head, v3, v4
        nonconform_local_global_perm(tail, head, v3, v4);
        // tail, head, v4, v3
        nonconform_local_global_perm(tail, head, v4, v3);
        // tail, v3, head, v4
        nonconform_local_global_perm(tail, v3, head, v4);
        // tail, v3, v4, head
        nonconform_local_global_perm(tail, v3, v4, head);
        // tail, v4, head, v3
        nonconform_local_global_perm(tail, v4, head, v3);
        // tail, v4, v3, head
        nonconform_local_global_perm(tail, v4, v3, head);

        // v4, tail, head, v3
        nonconform_local_global_perm(v4, tail, head, v3);
        // v4, tail, v3, head
        nonconform_local_global_perm(v4, tail, v3, head);
        // v4, head, tail, v3
        nonconform_local_global_perm(v4, head, tail, v3);
        // v4, head, v3, tail
        nonconform_local_global_perm(v4, head, v3, tail);
        // v4, v3, tail, head
        nonconform_local_global_perm(v4, v3, tail, head);
        // v4, v3, head, tail
        nonconform_local_global_perm(v4, v3, head, tail);


        
        // nonconform_local_global_perm(tail, head, v3, v4);
        // nonconform_local_global_perm(tail, v3, head, v4);
        // nonconform_local_global_perm(tail, v3, v4, head);
        // nonconform_local_global_perm(v3, head, tail, v4);

        
        // nonconform_local_global_perm(tail, head, v3, v4);
        // nonconform_local_global_perm(v3, head, tail, v4);
        // nonconform_local_global_perm(tail, v3, head, v4);
        // nonconform_local_global_perm(tail, v4, v3, head);

        // // --- v3 = tail orientation ---
        // nonconform_local_global_perm(v3, head, tail, v4);
        // nonconform_local_global_perm(tail, head, v3, v4);
        // nonconform_local_global_perm(v3, tail, head, v4);
        // nonconform_local_global_perm(v3, v4, tail, head);
      }
    }
  }
  /*
  GET_AUX_STORAGE(0, double *, sm);
      Vertex v1=tail;

      for(Vertex v2=1; v2 <= N_NODES; v2++){
	if(v2==v1) continue;
	double v12_old=sm[v1][v2];
	double v12_new=GETNEWWTOLD_M(v1,v2,v12_old);

	for(Vertex v3=1; v3 <= N_NODES; v3++){
	  if(v3==v2 || v3==v1) continue;
	  double v13_old=sm[v1][v3];
	  double v13_new=GETNEWWTOLD_M(v1,v3,v13_old);

	  if(v13_old<=v12_old && v13_new<=v12_new) continue;

	  double v32_old=sm[v3][v2];
	  double v32_new=GETNEWWTOLD_M(v3,v2,v32_old);

	  for(Vertex v4=1; v4 <= N_NODES; v4++){
	    if(v4==v3 || v4==v2 || v4==v1 ||
	       (head!=v4 && head!=v3 && head!=v2)) continue;

	    double v14_old=sm[v1][v4];
	    double v14_new=GETNEWWTOLD_M(v1,v4,v14_old);

	    if(v13_old<=v14_old && v13_new<=v14_new) continue;

	    double v34_old=sm[v3][v4];
	    double v34_new=GETNEWWTOLD_M(v3,v4,v34_old);

	    if(v13_old>v12_old && v13_old>v14_old && v12_old<=v14_old && v32_old>v34_old) CHANGE_STAT[0]--;
	    if(v13_new>v12_new && v13_new>v14_new && v12_new<=v14_new && v32_new>v34_new) CHANGE_STAT[0]++;
	  }
	}
      }

      Vertex v3=tail;

      for(Vertex v1=1; v1 <= N_NODES; v1++){
	if(v1==v3) continue;
	double v13_old=sm[v1][v3];
	double v13_new=GETNEWWTOLD_M(v1,v3,v13_old);

	for(Vertex v2=1; v2 <= N_NODES; v2++){
	  if(v2==v3 || v2==v1) continue;
	    double v12_old=sm[v1][v2];
	    double v12_new=GETNEWWTOLD_M(v1,v2,v12_old);

	    if(v13_old<=v12_old && v13_new<=v12_new) continue;

	    double v32_old=sm[v3][v2];
	    double v32_new=GETNEWWTOLD_M(v3,v2,v32_old);

	    for(Vertex v4=1; v4 <= N_NODES; v4++){
	      if(v4==v3 || v4==v2 || v4==v1 ||
		 (head!=v4 && head!=v2 && head!=v1)) continue;

	      double v14_old=sm[v1][v4];
	      double v14_new=GETNEWWTOLD_M(v1,v4,v14_old);

	      if(v13_old<=v14_old && v13_new<=v14_new) continue;

	      double v34_old=sm[v3][v4];
	      double v34_new=GETNEWWTOLD_M(v3,v4,v34_old);

	      if(v13_old>v12_old && v12_old<=v14_old && v32_old>v34_old) CHANGE_STAT[0]--;
	      if(v13_new>v12_new && v12_new<=v14_new && v32_new>v34_new) CHANGE_STAT[0]++;
	    }
	}
      }*/
}


// From Krivitsky and Butts paper, here, v1=i, v2=j, v3=l, v4=k.
WtS_CHANGESTAT_FN(s_localAND_nonconformity){
  GET_AUX_STORAGE(0, double *, sm);
  for(Vertex v1=1; v1 <= N_NODES; v1++){
    for(Vertex v2=1; v2 <= N_NODES; v2++){
      if(v2==v1) continue;
      double v12=sm[v1][v2];
      for(Vertex v3=1; v3 <= N_NODES; v3++){
	if(v3==v2 || v3==v1) continue;
	double v13=sm[v1][v3];
	if(v13<=v12) continue;
	double v32=sm[v3][v2];
	for(Vertex v4=1; v4 <= N_NODES; v4++){
	  if(v4==v3 || v4==v2 || v4==v1) continue;
	  double v14=sm[v1][v4];
	  if(v13<=v14) continue;
	  double v34=sm[v3][v4];
	  if(v32>v34 && v12<=v14) CHANGE_STAT[0]++;
	}
      }
    }
  }
}



WtD_FROM_S_FN(d_nonconformity_decay)

WtS_CHANGESTAT_FN(s_nonconformity_decay){
  GET_AUX_STORAGE(0, double *, sm);
  for(Vertex v1=1; v1 <= N_NODES; v1++){
    for(Vertex v2=1; v2 <= N_NODES; v2++){
      if(v2==v1) continue;
      double e = pow(INPUT_PARAM[0],INPUT_PARAM[1]-sm[v1][v2]);
      for(Vertex v3=1; v3 <= N_NODES; v3++){
	if(v3==v2 || v3==v1) continue;
	double v13=sm[v1][v3], v23=sm[v2][v3];
	for(Vertex v4=1; v4 <= N_NODES; v4++){
	  if(v4==v3 || v4==v2 || v4==v1) continue;
	  unsigned int
	    v134 = v13>sm[v1][v4],
	    v234 = v23>sm[v2][v4];
	  if(v134!=v234) CHANGE_STAT[0]+=e;
	}
      }
    }
  }
}

WtD_FROM_S_FN(d_nonconformity_thresholds)

WtS_CHANGESTAT_FN(s_nonconformity_thresholds){
  GET_AUX_STORAGE(0, double *, sm);
    for(Vertex v1=1; v1 <= N_NODES; v1++){
    for(Vertex v2=1; v2 <= N_NODES; v2++){
      if(v2==v1) continue;
      double v12 = sm[v1][v2];
      unsigned int i;
      for(i=0; i<N_CHANGE_STATS; i++) if(v12>=INPUT_PARAM[i]) break;
      if(i==N_CHANGE_STATS) continue;
      for(Vertex v3=1; v3 <= N_NODES; v3++){
	if(v3==v2 || v3==v1) continue;
	double v13=sm[v1][v3], v23=sm[v2][v3];
	for(Vertex v4=1; v4 <= N_NODES; v4++){
	  if(v4==v3 || v4==v2 || v4==v1) continue;
	  unsigned int
	    v134 = v13>sm[v1][v4],
	    v234 = v23>sm[v2][v4];
	  if(v134!=v234) CHANGE_STAT[i]+=1;
	}
      }
    }
  }
}

}
