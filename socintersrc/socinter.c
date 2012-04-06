/*
 * socinter.c
 *
 * implementation of a social interaction simulation
 * top level file
 * maarten
 */

#include <stdio.h>
#include <stdlib.h>
#include "socinterfuncs.h"

int main(int argc, char * argv[])
{
  if (argc != 14) {
    printf("\nsocinter: Simulation of linguistic change through social interaction.\nUsage:\n");
    printf("\t1. output path\n");
    printf("\t2. size of grid (20)\n");
    printf("\t3. number of steps (5000)\n");
    printf("\t4. seed for random (137)\n");
    printf("\t5. distance power (2)\n");
    printf("\t6. maximum age (5)\n");
    printf("\t7. scaler (c) [0.0..1.0]\n");
    printf("\t8. deviation factor [0.0..1.0]\n");
    printf("\t9. drift factor [0.0..]\n");
    printf("\t10. age distribution (0: cohorts 1: random)\n");
    printf("\t11. status distribution (0: uniform 1: normal)\n");
    printf("\t12. itemdistr (0: same for all 1: uniform 2: bimodal)\n");
    printf("\t13. markpercentile [0.0..1.0]\n\n");
  } else {
    int size = atoi(argv[2]);
    int nsteps = atoi(argv[3]);
    int seed = atoi(argv[4]);
    int dist = atoi(argv[5]);
    int maxage = atoi(argv[6]);
    float scaler = (float) atof( argv[7] );
    float dev = atof(argv[8]);
    float drift = atof(argv[9]);
    int agedistr = atoi(argv[10]);
    int statdistr = atoi(argv[11]);
    int itdistr = atoi(argv[12]);
    float markp = atof(argv[13]);
    char * path = argv[1];

    printf("\nSimulation of linguistic change through social interaction.\n");
    printf("\tSize:\t\t\t%d by %d\n",size,size);
    printf("\tNumber of steps:\t%d\n",nsteps);
    printf("\tSeed:\t\t\t%d\n",seed);
    printf("\nParameters:\n");
    printf("\tDistance power:\t\t%d\n",dist);
    printf("\tMaximum age:\t\t%d\n",maxage);
    printf("\tScaler (c):\t\t%.3f\n",scaler);
    printf("\tDeviation factor:\t%.3f\n",dev);
    printf("\tDrift factor:\t\t%.3f\n",drift);
    printf("\tMark percentile:\t%.3f\n",markp);
    printf("\nInitial settings:\n");
    printf("\tAge distribution:\t%s\n",(agedistr) ? "COHORT" : "RANDOM");
    char * itstring;
    switch(itdistr) {
    case 0:
      itstring = "SAME FOR ALL";
      break;
    case 1:
      itstring = "UNIFORM RANDOM";
      break;
    case 2:
      itstring = "BIMODAL";
      break;
    default:
      fprintf(stderr,"Illegal value for itemdistribution: %d\n",itdistr);
      exit(EXIT_FAILURE);
    }
    printf("\tItem distribution:\t%s\n",itstring);
    printf("\tStatus distribution:\t%s\n",(statdistr) ? "UNIFORM" : "NORMAL");
    printf("\nReporting to: %s\n\n",path);
    Simulation * sim = init_sim(size,
				nsteps,
				seed,
				dist,
				maxage,
				scaler,
				dev,
				drift,
				agedistr,
				statdistr,
				itdistr,
				markp,
				path
				);
				
    run(sim);
  }
  return 0;
}
  
