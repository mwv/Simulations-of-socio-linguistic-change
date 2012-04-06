/*
 * socimpact.c
 * implementation and extension of nettle's simulations
 * maarten
 */

#include <stdio.h>
#include <stdlib.h>
#include "socimpactfuncs.h"

int main(int argc, char * argv[])
{
  if (argc != 14) {
    printf("\nsocimpact: Social Impact simulation.\n");
    printf("\t1. Report path\n");
    printf("\t2. Size (grid is size by size)\n");
    printf("\t3. Number of steps\n");
    printf("\t4. Seed\n");
    printf("\t5. Maximum age\n");
    printf("\t6. Age distribution (0: random 1: cohorts)\n");
    printf("\t7. Number of items (> 1)\n");
    printf("\t8. Item distribution (0: same for all 1: random)\n");
    printf("\t9. Status distribution (0: same for all 1: poisson 2: hypers)\n");
    printf("\t10. Learning mode:\n\t\t0: maximize-maximize\n\t\t1. maximize-sample\n\t\t2. sample\n");
    printf("\t11. Bias [0.0-2.0]\n");
    printf("\t12. Mutation rate [0.0-1.0]\n");
    printf("\t13. Norm impact [1.0-2.0]\n\n");
  } else {
    char * path = argv[1];
    int size = atoi(argv[2]);
    int nsteps = atoi(argv[3]);
    int seed = atoi(argv[4]);
    int maxage = atoi(argv[5]);
    int agedistr = atoi(argv[6]);
    int nitems = atoi(argv[7]);
    int itemdistr = atoi(argv[8]);
    int statdistr = atoi(argv[9]);
    int learningmode = atoi(argv[10]);
    float bias = atof(argv[11]);
    float murate = atof(argv[12]);
    float normimpact = atof(argv[13]);

    printf("\nsocimpact: Social Impact simulation.\n");
   
    printf("\tNumber of steps:\t%d\n",nsteps);
    printf("\tSeed:\t\t\t%d\n",seed);

    printf("\nParameters:\n");
    printf("\tSize:\t\t\t%d by %d\n",size, size);
    printf("\tMaximum age:\t\t%d\n",maxage);
    printf("\tNumber of items:\t%d\n",nitems);
    printf("\tBias:\t\t\t%.1f\n",bias);
    printf("\tMu rate:\t\t%.2f\n",murate);
    printf("\tNorm impact:\t\t%.1f\n",normimpact);
    char * lm;
    switch(learningmode) {
    case 0:
      lm = "Maximize - maximize";
      break;
    case 1:
      lm = "Maximize - sample";
      break;
    default:
      lm = "Sample";
      break;
    }
    printf("\tLearning mode:\t\t%s\n",lm);

    printf("\nInitial settings:\n");
    printf("\tAge distribution:\t%s\n",(agedistr) ? "COHORTS" : "RANDOM");
    printf("\tItem distribution:\t%s\n",(itemdistr) ? "RANDOM" : "SAME FOR ALL");
    char * stat;
    switch(statdistr) {
    case 0:
      stat = "SAME FOR ALL";
      break;
    case 1:
      stat = "POISSON";
      break;
    default:
      stat = "HYPERS";
      break;
    }
    printf("\tStatus distribution:\t%s\n",stat);
    printf("\nReporting to: %s\n\n",path);
    Simulation * sim = init_sim(size,
				nsteps,
				seed,
				maxage,
				agedistr,
				nitems,
				itemdistr,
				statdistr,
				learningmode,
				bias,
				murate,
				normimpact,
				path);
    run(sim);
    printf("Simulation done.\n");
  }
  return 0;
}
    
    

    
