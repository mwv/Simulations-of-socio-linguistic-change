/*
 * extension of nettle sims
 * maarten
 */

#ifndef _SOCIMPACTFUNCS_H
#define _SOCIMPACTFUNCS_H

typedef struct {
  int item;
  int status;
  int age;
} Agent;

typedef struct {
  Agent * grid;
  int size;
  int seed;
  int maxage;
  int agedistr;
  int nitems;
  int itemdistr;
  int statdistr;
  int learningmode;
  float bias;
  float murate;
  float normimpact;
  FILE * shortreportFP;
  FILE * longreportFP;
  FILE * finalreportFP;
  int nsteps;
  int currentstep;
  int mostfrequent;
  int nchanges;
  float tothomog;
} Simulation;

Simulation * init_sim(int size,
		      int nsteps,
		      int seed,
		      int maxage,
		      int agedistr, /* 0: RANDOM 1: COHORTS */
		      int nitems,
		      int itemdistr, /* 0: LOWEST 1: RANDOM */
		      int statdistr, /* 0: SAME 1: POISSON 2: HYPERS */
		      int learningmode, /* 0: M-M 1: M-S 2: S */
		      float bias,
		      float murate,
		      float normimpact,
		      char * path
		      );
void run(Simulation *);

#endif /* _SOCIMPACTFUNCS_H */
