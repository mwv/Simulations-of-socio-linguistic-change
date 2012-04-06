#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "socimpactfuncs.h"

#ifndef max
#define max(a , b) ( ((a) > (b)) ? (a) : (b) )
#endif /* max */

#ifndef min
#define min(a , b) ( ((a) < (b)) ? (a) : (b) )
#endif /* min */

#define DEBUG 0

#ifndef DPRINT
#define DPRINT(e) printf((DEBUG) ? ("DEBUG " e) : "");
#endif

#ifndef VPRINT
#define VPRINT(e) printf((DEBUG) ? ("DEBUG " #e ":\t%g\n", e) : "")
#endif

#define LONGREPORT 0

#define NAME_BUF_SIZE 300
#define HYPER_THRESH 0.025
#define EPSILON 0.000001

/* prototypes */

static void end_sim(Simulation *);
static float rand01(void);
static char * makefilename(Simulation *, char *, char *);
static float distance(int, int, int);
static void step(Simulation *);
static void report(Simulation *);
static void reportfinal(Simulation *);
static void collectimpacts(Simulation *, int, float *);
static int sample(Simulation *, float *);
static int maxidx_float(float *, int);
static int maxidx_int(int *, int);

Simulation * init_sim(int size,
		      int nsteps,
		      int seed,
		      int maxage,
		      int agedistr,
		      int nitems,
		      int itemdistr,
		      int statdistr,
		      int learningmode,
		      float bias,
		      float murate,
		      float normimpact,
		      char * path
		      )
{
  Simulation * sim = (Simulation *) malloc (sizeof(Simulation));
  assert(sim);
  
  // init bookkeeping vars
  sim -> nsteps = nsteps;
  sim -> currentstep = 0;
  sim -> nchanges = 0;
  sim -> tothomog = 0.0;

  // carry sim params
  sim -> size = size;
  sim -> maxage = maxage;
  sim -> agedistr = agedistr;
  sim -> nitems = nitems;
  sim -> itemdistr = itemdistr;
  sim -> statdistr = statdistr;
  sim -> learningmode = learningmode;
  sim -> bias = bias;
  sim -> murate = murate;
  sim -> normimpact = normimpact;

  // init rand()
  sim -> seed = seed;
  srand(seed);
  
  // init file pointers
  char * shortreport = makefilename(sim, path, "short");
  sim -> shortreportFP = fopen(shortreport,"w");
  assert(sim -> shortreportFP);
  free(shortreport);

  if (LONGREPORT) {
    char * longreport = makefilename(sim, path, "long");
    sim -> longreportFP = fopen(longreport, "w");
    assert(sim -> longreportFP);
    free(longreport);
  } else 
    sim -> longreportFP = NULL;

  char * finalreport = makefilename(sim, path, "final");
  sim -> finalreportFP = fopen(finalreport, "w");
  assert(sim -> finalreportFP);
  free(finalreport);

  // allocate grid
  sim -> grid = (Agent *) malloc(size * size * sizeof(Agent));

  // populate grid while keeping track of item numbers
  int i;
  int itemsums[nitems];
  for (i = 0; i < nitems; ++i)
    itemsums[i] = 0;

  for (i = 0; i < size * size; ++i) {
    int xpos = i / size; // xposition on grid
    // determine status
    int status;
    switch (statdistr) {
    case 0: // all the same
      status = 1;
      break;
    case 2: // hypers
      if (rand01() < HYPER_THRESH) {
	status = size * size * 25;
	break;
      }
    case 1: // poisson approx
      status = (int) pow(rand01() * (size - 1) + 1, 2.0);
      break;
    default:
      printf("Illegal value for statdistr: %d\n",statdistr);
      exit(EXIT_FAILURE);
    }
    // determine age
    int age;
    if (agedistr) { // cohort distr
      if (xpos % maxage == xpos % (maxage * 2))
	age = xpos % maxage + 1;
      else
	age = maxage - (xpos % maxage);
    } else 
      age = (int) round(rand01() * (maxage-1)) + 1;
    // determine item
    int item;
    if (itemdistr)  // random items
      item = (int) floor(rand01() * nitems);
    else
      item = 0; // default to lowest item
    itemsums[item]++;
    Agent a = {item,status,age};
    sim -> grid[i] = a;
  }
  // determine most frequent item
  sim -> mostfrequent = maxidx_int(itemsums, nitems);

  return sim;
}

void run(Simulation * sim)
{
  report(sim);
  while (sim -> currentstep++ < sim -> nsteps) {
    step(sim);
    report(sim);
  }
  end_sim(sim);
}

static void step(Simulation * sim)
{
  int size = sim -> size;
  // copy the grid;
  Agent * newgrid = (Agent *) malloc(size * size * sizeof(Agent));
  int i;
  for (i = 0; i < size * size; ++i) {
    Agent a = {sim -> grid[i].item,sim -> grid[i].status,sim -> grid[i].age};
    sim -> grid[i] = a;
  }

  // update the agents
  for (i = 0; i < size * size; ++i) {
    // determine if item should be reset
    int item = sim -> grid[i].item;
    if (sim -> grid[i].age <= 2) {
      float impacts[sim -> nitems];
      collectimpacts(sim, i, impacts);
      switch(sim -> learningmode) {
      case 0: // maximize - maximize
	if (rand01() > sim -> murate) { // normal procedure: maximize
	  item = maxidx_float(impacts, sim -> nitems);
	} else { // exception procedure: maximize
	  impacts[maxidx_float(impacts, sim -> nitems)] = -1.0;
	  item = maxidx_float(impacts, sim -> nitems); // take second best
	}
	break;
      case 1: // maximize - sample
	if (rand01() > sim -> murate) { // normal procedure: maximize
	  item = maxidx_float(impacts, sim -> nitems);
	} else { // exception : sample
	  impacts[maxidx_float(impacts, sim -> nitems)] = -1.0;
	  item = sample(sim, impacts);
	}
	break;
      case 2: // sample - sample
	item = sample(sim, impacts);
	break;
      }
    }
    // determine status
    int status = sim -> grid[i].status;
    
    // determine age
    int age = 1 + sim -> grid[i].age;
    if (age > sim -> maxage) {
      age = 1;
      // determine new status
      switch (sim -> statdistr) {
      case 0: // all the same
	status = 1;
	break;
      case 2: // hypers
	if (rand01() < HYPER_THRESH) {
	  status = sim -> size * sim -> size * 25;
	  break;
	}
      case 1: // poisson approx
	status = (int) pow(rand01() * (sim -> size - 1) + 1, 2.0);
	break;
      }
    }
    Agent a = {item,status,age};
    newgrid[i] = a;
  }

  free(sim -> grid);
  sim -> grid = newgrid;
}

/* sample index from impacts, -1 indicates nonvalid index */
static int sample(Simulation * sim, float * impacts)
{
  float sum = 0.0;
  int nzeros = 0;
  int i;
  for (i = 0; i < sim -> nitems; ++i) {
    if (fabs(impacts[i]) <= EPSILON) // effectively 0
      nzeros++;
    if (impacts[i] > 0) // exclude -1
      sum += impacts[i];
  }
  
  float reservedmass = (sum / (1.0 - sim -> murate)) - sum;
  float massperzero = reservedmass / ((float) nzeros);
  sum += reservedmass;

  float normalized[sim -> nitems];
  if (fabs(sum) > EPSILON) {
    for (i = 0; i < sim -> nitems; ++i) {
      if (fabs(impacts[i] < EPSILON)) {
	normalized[i] = massperzero / sum;
      } else if (impacts[i] < 0) { //exclusion case
	normalized[i] = 0.0;
      } else {
	normalized[i] = impacts[i] / sum;
      }
    }
  } else {
    for (i = 0; i < sim -> nitems; ++i) 
      normalized[i] = 1.0/((float) sim -> nitems);
  }

  // now sample one index:
  while(1) {
    float r = rand01();
    for (i = 0; i < sim -> nitems; ++i) {
      if (normalized[i] > r)
	return i;
      else 
	r -= normalized[i];
    }
  }
}

static int maxidx_float(float * arr, int n)
{
  int maxidx = 0;
  int i;
  for (i = 1; i < n; ++i) {
    if (arr[i] > arr[maxidx])
      maxidx = i;
  }
  return maxidx;
}

static int maxidx_int(int * arr, int n)
{
  int maxidx = 0;
  int i;
  for (i = 1; i < n; ++i) {
    if (arr[i] > arr[maxidx])
      maxidx = i;
  }
  return maxidx;
}

static void collectimpacts(Simulation *  sim, int idx, float * arr)
{
  int sums[sim -> nitems];
  float status_over_dist_sums[sim -> nitems];
  int i;
  for (i = 0; i < sim -> nitems; ++i) {
    sums[i] = 0;
    status_over_dist_sums[i] = 0.0;
  }

  int j;
  for (j = 0; j < sim -> size * sim -> size; ++j) {
    if (idx != j) {
      if (sim -> grid[j].age > 1) {
	sums[sim -> grid[j].item]++;
	float dist = distance(idx,j,sim -> size);
	status_over_dist_sums[sim -> grid[j].item] += (float) sim -> grid[j].status / (dist * dist);
      }
    }
  }

  int specialitem = sim -> nitems - 1; // only item with bias
  if (sums[specialitem] != 0)
    arr[specialitem] = sim -> bias * pow(sums[specialitem],sim -> normimpact) * (status_over_dist_sums[specialitem] / ((float) sums[specialitem]));
  else
    arr[specialitem] = 0.0;
  for (i = 0; i < sim -> nitems - 1; ++i) {
    if (sums[i] != 0)
      arr[i] = pow(sums[i], sim -> normimpact) * (status_over_dist_sums[i] / ((float) sums[i]));
    else
      arr[i] = 0;
  }
}

static void report(Simulation * sim) 
{
  int i;
  int items[sim -> nitems];
  for (i = 0; i < sim -> nitems; ++i)
    items[i] = 0;
  for (i = 0; i < sim -> size * sim -> size; ++i) 
    items[sim -> grid[i].item]++;

  int mostfrequent = maxidx_int(items, sim -> nitems);
  if (sim -> mostfrequent != mostfrequent) {
    sim -> nchanges++;
    sim -> mostfrequent = mostfrequent;
  }

  float homogeneity = (float) items[mostfrequent] / ((float) sim -> size * sim -> size);
  sim -> tothomog += homogeneity;

  // write short report
  fprintf(sim -> shortreportFP,
	  "%d\t%d\t%.3f",
	  sim -> currentstep,
	  sim -> mostfrequent,
	  homogeneity);
  for (i = 0; i < sim -> nitems; ++i) 
    fprintf(sim -> shortreportFP,"\t%d",items[i]);
  fprintf(sim -> shortreportFP,"\n");

  // write long report
  if (LONGREPORT) {
    for (i = 0; i < sim -> size * sim -> size; ++i)
      fprintf(sim -> longreportFP,
	      "%d %d %d ",
	      sim -> grid[i].status,
	      sim -> grid[i].item,
	      sim -> grid[i].age
	      );
    fprintf(sim -> longreportFP,"\n");
  }
}

static void reportfinal(Simulation * sim)
{
  fprintf(sim -> finalreportFP, "0. Simulation summary:\n");
  fprintf(sim -> finalreportFP, "1. Seed:\t%d\n",sim -> seed);
  fprintf(sim -> finalreportFP, "2. Size:\t%d\n",sim -> size);
  fprintf(sim -> finalreportFP, "3. Number of steps:\t%d\n",sim -> nsteps);
  fprintf(sim -> finalreportFP, "4. Max age:\t%d\n",sim -> maxage);
  fprintf(sim -> finalreportFP, "5. Age distribution\t%s\n",(sim -> agedistr) ? "COHORTS" : "RANDOM");
  fprintf(sim -> finalreportFP, "6. Number of items:\t%d\n",sim -> nitems);
  fprintf(sim -> finalreportFP, "7. Item distribution:\t%s\n",(sim -> itemdistr) ? "RANDOM" : "SAME");
  char * stat;
  switch(sim -> statdistr) {
  case 0:
    stat = "SAME";
    break;
  case 1:
    stat = "POISSON";
    break;
  default:
    stat = "HYPER";
    break;
  }
  fprintf(sim -> finalreportFP, "8. Status distribution:\t%s\n",stat);
  char * lm;
  switch(sim -> learningmode) {
  case 0: 
    lm = "M-M";
    break;
  case 1:
    lm = "M-S";
    break;
  default:
    lm = "S";
    break;
  }
  fprintf(sim -> finalreportFP, "9. Learning mode:\t%s\n",lm);
  fprintf(sim -> finalreportFP, "10. Bias:\t%.1f\n",sim -> bias);
  fprintf(sim -> finalreportFP, "11. Mu:\t%.2f\n",sim -> murate);
  fprintf(sim -> finalreportFP, "12. Norm impact:\t%.2f\n", sim -> normimpact);
  fprintf(sim -> finalreportFP, "13. Number of changes:\t%d\n", sim -> nchanges);
  fprintf(sim -> finalreportFP, "14. Average homogeneity:\t%.3f\n",(sim -> tothomog / ((float) sim -> nsteps)));
}
  
  

/* helpers */
static void end_sim(Simulation * sim)
{
  reportfinal(sim);
  fclose(sim -> shortreportFP);
  if (LONGREPORT)
    fclose(sim -> longreportFP);
  fclose(sim -> finalreportFP);
  free(sim -> grid);
  free(sim);
}

static float rand01() 
{
  return (float) rand() / ((float) RAND_MAX + 1);
}


static char * makefilename(Simulation * sim, char * path, char * type)
{
  char * buffer = (char *) malloc(NAME_BUF_SIZE);
  sprintf(buffer,
	  "%sreport_%s_seed %d_size_%d_maxage_%d_aged_%d_nitems_%d_itemd_%d_statd_%d_learnm_%d_bias_%.1f_mu_%.2f_norm_%.2f.data",
	  path,
	  type,
	  sim -> seed,
	  sim -> size,
	  sim -> maxage,
	  sim -> agedistr,
	  sim -> nitems,
	  sim -> itemdistr,
	  sim -> statdistr,
	  sim -> learningmode,
	  sim -> bias,
	  sim -> murate,
	  sim -> normimpact
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


    
