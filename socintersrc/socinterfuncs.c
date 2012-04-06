#include <stdio.h>
#include <stdlib.h>

#include <math.h>

#include "socinterfuncs.h"

#ifndef max
#define max(a , b) ( ((a) > (b)) ? (a) : (b) )
#endif /* max */

#ifndef min
#define min(a , b) ( ((a) < (b)) ? (a) : (b) )
#endif /* min */

#define DEBUG 0
#define LONGREPORT 0

#define NDEBUG
#include <assert.h>

#ifndef DPRINT
#define DPRINT(e) printf((DEBUG) ? ("DEBUG " e) : "");
#endif

#ifndef VPRINT
#define VPRINT(e) printf((DEBUG) ? ("DEBUG " #e ":\t%g\n", e) : "")
#endif

#define NAME_BUF_SIZE 300
#define HYPER_THRESH 0.025

/* prototypes */
static void step(Simulation *);
static void report(Simulation *);
static void reportfinal(Simulation *);
static void end_sim(Simulation *);
static int maxidx(int *, int);
static void sim_free(Simulation *);
static float rand01(void);
static char * makefilename(Simulation *, char *, char *);
static float distance(int, int, int);
static float convert(float, float, float, float, float);
static int cmp_stat(const void *, const void *);


/* functions */
Simulation * init_sim(
		      int size,
		      int nsteps,
		      int seed,
		      int distpower,
		      int maxage,
		      float c, /* scales differentiation */
		      float deviationfactor,
		      float driftfactor,
		      int agedistr, /* 0: random 1: cohorts */
		      int statusdistr, /* 0: uniform 1: normal */
		      int itemdistr, /* 0: uniform 1: random 2: bimodal */
		      float markpercentile, /* 0 < mpc < 1 */
		      char * reportpath
		      ) {
  // DPRINT("Entering init_sim...\n");
  Simulation * sim = (Simulation *) malloc (sizeof(Simulation));
  if (!sim) {
    fprintf(stderr,"Memory allocation failure: Simulation.\n");
    exit(EXIT_FAILURE);
  }

  // initialize bookkeeping fields in sim struct:
  sim -> nsteps = nsteps;
  sim -> seed = seed;
  sim -> currentstep = 0;

  // initialize simulation parameters
  sim -> distpower = distpower;
  sim -> size = size;
  sim -> c = c;
  sim -> deviationfactor = deviationfactor;
  sim -> driftfactor = driftfactor;
  sim -> distpower = distpower;
  sim -> maxage = maxage;
  sim -> agedistr = agedistr;
  sim -> statusdistr = statusdistr;
  sim -> itemdistr = itemdistr;
  sim -> markpercentile = markpercentile;
  
  // initialize file pointers
  char * shortreport = makefilename(sim, reportpath, "short");
  sim -> shortreportFP = fopen(shortreport,"w");
  assert(sim -> shortreportFP);
  free(shortreport);
  
  if (LONGREPORT) {
    char * longreport = makefilename(sim, reportpath, "long");
    sim -> longreportFP = fopen(longreport, "w");
    assert(sim -> longreportFP);
    free(longreport);
  } else {
    sim -> longreportFP = NULL;
  }
  

  char * finalreport = makefilename(sim, reportpath, "final");
  sim -> finalreportFP = fopen(finalreport,"w");
  assert(sim -> finalreportFP);
  free(finalreport);

  sim -> grid = malloc(size*size*sizeof(Agent));
  assert(sim -> grid);

  // initialize rand sequence
  srand(seed);

  // initialize population grid
  int i;
  for (i = 0; i < size * size; ++i) {
    // calculate x,y pos on grid
    int xpos = i / size;
    // determine status
    float status;
    if (statusdistr) 
      status = pow(rand01() * size, 2.0) / (size*size);
    else
      status = rand01();
    // determine age
    int age;
    if (agedistr) {
      if (xpos % maxage == xpos % (maxage * 2))
	age = xpos % maxage + 1;
      else
	age = maxage - (xpos % maxage);
    }
    else
      age = (int) (rand01() * maxage);
    // determine starting item
    float item = 0.5;
    switch(itemdistr) {
    case 0:
      item = 0.5;
      break;
    case 1:
      item = rand01();
      break;
    case 2:
      item = (status > 0.5) ? 0.75 : 0.25;
      break;
    }
    float utility = 0.0;
    Agent a = {status, item, age, utility};
    sim -> grid[i] = a;
  }

  return sim;
}

void run(Simulation * sim)
{
  printf("Starting run...\n");
  
  report(sim);
  while (sim -> currentstep++ < sim -> nsteps) {
    step(sim);
    report(sim);
  }
  end_sim(sim);
  printf("Run done.\n");
}

static void step(Simulation * sim)
{
  int size = sim -> size;

  // two copies, one the new, other sorted
  Agent * newgrid = (Agent *) malloc(size*size*sizeof(Agent));
  Agent * sortedgrid = (Agent *) malloc (size * size * sizeof(Agent));
  int i;
  for (i = 0; i < size * size; ++i) {
    Agent a = sim -> grid[i];
    sortedgrid[i] = a;
    newgrid[i] = a;
  }
  // sort the grid into sortedgrid
  qsort(sortedgrid, size*size, sizeof(Agent), cmp_stat);

  // determine high and lowmarks
  int lowpercentilemark = (int) ((float) (size * size) * (sim -> markpercentile));
  int highpercentilemark = size * size - lowpercentilemark;
  float itemsum = 0.0;
  for (i = 0; i < lowpercentilemark; ++i)
    itemsum += sortedgrid[i].item;

  float lowmark = itemsum / lowpercentilemark;

  itemsum = 0.0;
  for (i = highpercentilemark; i < size * size; ++i)
    itemsum += sortedgrid[i].item;
  float highmark = itemsum / lowpercentilemark;
  free(sortedgrid);

  // calculate all item statuses
  float itstats[size*size];
  for (i = 0; i < size * size; ++i) {
    if (lowmark != highmark) 
      itstats[i] = (sim -> grid[i].item - lowmark)/(highmark-lowmark);
    else
      itstats[i] = 0.0;
  }

  // calculate all utilitygains
  float utgains[size*size];
  // init to zero
  for (i = 0; i < size * size; ++i) 
    utgains[i] = 0.0;
  int j;
  for (i = 0; i < size * size; ++i) {
    for (j = i + 1; j < size * size; ++j) {
      Agent a1 = sim -> grid[i];
      Agent a2 = sim -> grid[j];
      // calculate conformity
      float socdist = fabs(a1.status - a2.status);
      float itdist = fabs(a1.item - a2.item);
      float conf;
      float confdev = fabs(socdist-itdist); // deviation from line (0,0)-(1,1)
      if (confdev < sim -> deviationfactor) 
	conf = convert(confdev,0.0,sim -> deviationfactor,1.0,0.0); 
      else
	conf = convert(confdev,sim -> deviationfactor,1.0,0.0,-1.0);
      float c = sim -> c;
      float eucldist = distance(i,j,size);
      float ut1 = (c*(itstats[i] - itstats[j]) + (1.0-c)*conf)/(pow(eucldist,(float) sim -> distpower));
      float ut2 = (c*(itstats[j] - itstats[i]) + (1.0-c)*conf)/(pow(eucldist,(float) sim -> distpower));
      utgains[i] += ut1;
      utgains[j] += ut2;
    }
  }

  // update the agents
  for (i = 0; i < size * size; ++i) {
    newgrid[i].utility += utgains[i];
    newgrid[i].age++;
    if (newgrid[i].age > sim -> maxage) {
      // determine drift
      float drift;
      if (newgrid[i].utility > 0) {
	drift = rand01() * 0.02 - 0.01;
      } else { // drift
	drift = (rand01()*2.0 -1.0) * sim -> driftfactor * fabs(newgrid[i].utility);
      }
      newgrid[i].item = min(1.0,max(0.0,drift + newgrid[i].item));
      newgrid[i].age = 0;
      newgrid[i].utility = 0.0;
    }  
  }

  free(sim -> grid);
  sim -> grid = newgrid;
}

static void report(Simulation * sim)
{
  // calculate bins and itemsum
  int size = sim -> size;
  Agent * sorted = (Agent *) malloc(size * size * sizeof(Agent));
  int i;
  int bin[11] = {0};
  float itemsum = 0.0;
  for (i = 0; i < size * size; ++i) {
    sorted[i] = sim -> grid[i];
    bin[(int) round(sorted[i].item * 10.0)]++;
    itemsum += sorted[i].item;
  }
  // calc mostfrequent and determine whether it has changed
  int mostfrequent = maxidx(bin, 11);
  if (sim -> mostfrequent != mostfrequent) {
    sim -> numberofchanges++;
    sim -> mostfrequent = mostfrequent;
  }
  // calc homog & avgitem
  float homogeneity = (float) bin[mostfrequent] / (size*size);
  sim -> tothomog += homogeneity;
  float avgitem = itemsum / (size*size);
  // calc lowmark & highmark
  
  // sort the grid into sortedgrid
  qsort(sorted, size*size, sizeof(Agent), cmp_stat);

  // determine high and lowmarks
  int lowpercentilemark = (int) ((float) (size * size) * (sim -> markpercentile));
  int highpercentilemark = size * size - lowpercentilemark;
  
  float itsum = 0.0;
  for (i = 0; i < lowpercentilemark; ++i)
    itsum += sorted[i].item;
  float lowmark = itsum / lowpercentilemark;
  itsum = 0.0;
  for (i = highpercentilemark; i < size * size; ++i)
    itsum += sorted[i].item;
  float highmark = itsum / lowpercentilemark;
  free(sorted);
  
  // write to short report
  fprintf(sim -> shortreportFP, 
	  "%d\t%d\t%.3f\t%.3f\t%.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", 
	  sim -> currentstep,
	  sim -> mostfrequent,
	  lowmark,
	  highmark,
	  avgitem,
	  bin[0],
	  bin[1],
	  bin[2],
	  bin[3],
	  bin[4],
	  bin[5],
	  bin[6],
	  bin[7],
	  bin[8],
	  bin[9],
	  bin[10]);

  if (LONGREPORT) {
    // write to long report
    for (i = 0; i < size * size; ++i)
      fprintf(sim -> longreportFP,
	      "%.3f %.3f %d %.3f ",
	      sim -> grid[i].status,
	      sim -> grid[i].item,
	      sim -> grid[i].age,
	      sim -> grid[i].utility);
    fprintf(sim -> longreportFP,"\n");
  }

}

static void reportfinal(Simulation * sim)
{
  fprintf(sim -> finalreportFP,"0. Simulation summary:\n");
  fprintf(sim -> finalreportFP,"1. Seed:\t%d\n",sim -> seed);
  fprintf(sim -> finalreportFP,"2. Size:\t%d\n",sim -> size);
  fprintf(sim -> finalreportFP,"3. Number of steps:\t%d\n",sim -> nsteps);
  fprintf(sim -> finalreportFP,"4. Max age:\t%d\n",sim -> maxage);
  fprintf(sim -> finalreportFP,"5. Scaler (c):\t%.5f\n",sim -> c);
  fprintf(sim -> finalreportFP,"6. Deviation factor:\t%.5f\n",sim -> deviationfactor);
  fprintf(sim -> finalreportFP,"7. Drift factor:\t%.5f\n",sim -> driftfactor);
  fprintf(sim -> finalreportFP,"8. Distance power:\t%d\n",sim -> distpower);
  fprintf(sim -> finalreportFP,"9. Status distribution:\t%s\n",(sim -> statusdistr) ? "NORMAL" : "UNIFORM");
  fprintf(sim -> finalreportFP,"10. Age distribution:\t%s\n",(sim -> agedistr) ? "RANDOM" : "COHORTS");
  char * stats;
  switch(sim -> itemdistr) {
  case 0:
    stats = "IDENTICAL";
    break;
  case 1:
    stats = "UNIFORM RAND";
    break;
  case 2:
    stats = "BIMODAL";
    break;
  default:
    exit(EXIT_FAILURE);
  }
  fprintf(sim -> finalreportFP,"11. Item distribution:\t%s\n",stats);
  fprintf(sim -> finalreportFP,"12. Number of changes:\t%d\n",sim -> numberofchanges);
  fprintf(sim -> finalreportFP,"13. Average homogeneity:\t%.5f\n",(sim -> tothomog) / (sim -> nsteps));
  fprintf(sim -> finalreportFP,"14. Mark percentile:\t%.5f\n",(sim -> markpercentile));
}

static void end_sim(Simulation * sim)
{
  reportfinal(sim);
  fclose(sim -> shortreportFP);
  if (LONGREPORT)
    fclose(sim -> longreportFP);
  fclose(sim -> finalreportFP);
  sim_free(sim);
}

static int maxidx(int * arr, int size)
{
  int i;
  int maxidx = 0;
  for (i = 1; i < size; ++i) {
    if (arr[i] > arr[maxidx])
      maxidx = i;
  }
  return maxidx;
}

static void sim_free(Simulation * sim)
{
  free(sim -> grid);
  free(sim);
}

/* helper functions */
static float rand01() 
{
  return (float) rand() / ((float) RAND_MAX + 1);
}

/* makefilename: prepare reportfilenames */
static char * makefilename(Simulation * sim, char * stem, char * type)
{
  char * buffer = (char *) malloc (NAME_BUF_SIZE * sizeof(char));
  sprintf(buffer,
	  "%sreport_%s_seed_%d_size_%d_c_%.2f_dev_%.1f_drift_%.2f_dist_%d_mark_%.2f_statd_%d_aged_%d_itd_%d.data",
	  stem,
	  type,
	  sim -> seed,
	  sim -> size,
	  sim -> c,
	  sim -> deviationfactor,
	  sim -> driftfactor,
	  sim -> distpower,
	  sim -> markpercentile,
	  sim -> statusdistr,
	  sim -> agedistr,
	  sim -> itemdistr
	  );
  return buffer;
}

/* distance: calculate euclidean distance from array indices */
static float distance(int pos1, int pos2, int size) 
{
  int x1 = pos1 / size;
  int y1 = pos1 % size;
  int x2 = pos2 / size;
  int y2 = pos2 % size;

  int xdiff = fabs(x1 - x2);
  int ydiff = fabs(y1 - y2);
  
  xdiff = min(xdiff, fabs(size - xdiff));
  ydiff = min(ydiff, fabs(size - ydiff));

  return sqrt(xdiff*xdiff + ydiff*ydiff);
}


/* range scaler */
static float convert(float x, float inmin, float inmax, float outmin, float outmax)
{
  return ((x-inmin)/(inmax-inmin))*(outmax-outmin) + outmin;
}

/* cmp agents by status */
static int cmp_stat(const void * vp1, const void * vp2)
{
  Agent * ap1 = (Agent *) vp1;
  Agent * ap2 = (Agent *) vp2;
  if (ap1 -> status < ap2 -> status)
    return -1;
  if (ap1 -> status > ap2 -> status)
    return 1;
  return 0;
}
  
  
