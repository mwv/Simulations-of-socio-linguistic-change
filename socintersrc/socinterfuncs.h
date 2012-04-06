/*
 * socinterfuncs.h
 * mar 8, 2010
 * author maarten
 * re-implementation of sociolinguistic simulations
 */

#ifndef SOCINTERFUNCS_H_
#define SOCINTERFUNCS_H_

typedef struct {
  float status;
  float item;
  int age;
  float utility;
} Agent;

typedef struct {
  Agent * grid;
  int size;
  int nsteps;
  int maxage;
  float c;
  float deviationfactor;
  float driftfactor;
  float markpercentile;
  int distpower;
  int seed;
  int statusdistr; /* 0: uniform 1: normal */
  int agedistr; /* 0: random 1: cohorts */
  int itemdistr; /* 0: identical 1: uniform rand 2: bimodal */
  FILE * shortreportFP;
  FILE * longreportFP;
  FILE * finalreportFP;
  int mostfrequent;
  int numberofchanges;
  float tothomog;
  int currentstep;
} Simulation;


Simulation * init_sim(int size,
		      int nsteps,
		      int seed,
		      int distpower,
		      int maxage,
		      float c,
		      float deviationfactor,
		      float driftfactor,
		      int agedistr,
		      int statusdistr,
		      int itemdistr,
		      float markpercentile,
		      char * reportpath
		      ); 
void run(Simulation *);

#endif /* SOCINTERFUNCS_H_ */
