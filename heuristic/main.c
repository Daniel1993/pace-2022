#define _GNU_SOURCE
#include <sched.h>
#include <stdlib.h>

#include <sys/types.h>
#include <unistd.h>

#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <sys/resource.h>

#include "sa_lib.h"
#include "aux.h"

volatile int sigterm_called;

int
main(int argc,
     char** argv
     )
{
  /* Usage: ./project fileName n hotP hotD coldP coldD repOps */
  /*  */
  /* fileName: is the name of the file that contains the graph. */
  /* hotP: is the percentage for hot. */
  /* hotD: is the change for hot < 0. */
  /* coldP: is the percentage for cold. */
  /* coldD: is the change for cold < 0. */

  char fileNameBuffer[1<<10];
  char *fileName = NULL; // if set to NULL read from stdin
  // int exp;
  SA_parameters_s params = {
    .hot = 1.0e-5,
    .hotD = 1,
    .cold = 1.0e-15,
    .coldD = 1,
    .decay = 0.95
  };

  /* Default values are given in variable definition */
  switch(argc){
    // TODO: the compiler has a warning for this -Wimplicit-fallthrough=
  // case 9:
  //   sscanf(argv[8], "%lf", &params.decay);
  // case 8:
  //   sscanf(argv[7], "%d", &exp);
  //   params.repOps = 1<<exp;
  // case 7:
  //   sscanf(argv[6], "%d", &exp);
  //   params.n = 1<<exp;
  // case 6:
  //   sscanf(argv[5], "%d", &params.coldD);
  // case 5:
  //   sscanf(argv[4], "%lf", &params.cold);
  // case 4:
  //   sscanf(argv[3], "%d", &params.hotD);
  // case 3:
  //   sscanf(argv[2], "%lf", &params.hot);
  case 2:
    fileName = &fileNameBuffer[0];
    strcpy(fileName, argv[1]);
    break;
  case 1: default:
    break;
  }

  handle_sigterm();

  const rlim_t kStackSize = 1.0 * 1024 * 1024 * 1024; // 1GB (my system does not allow more than this)
  struct rlimit rl;
  // int r;
  /* r =  */getrlimit(RLIMIT_STACK, &rl); // TODO: check error 
  // fprintf(stderr, "Before first TarjanVisit (old stack size = %ld (ret=%i))\n", rl.rlim_cur, r);
  rl.rlim_cur = kStackSize;
  /* r =  */setrlimit(RLIMIT_STACK, &rl);
  // fprintf(stderr, "Before first TarjanVisit (new stack size = %ld (ret=%i))\n", rl.rlim_cur, r);

  // fprintf(stderr, "Before reading graph\n");
  SA_init_F(params, fileName);

  // fprintf(stderr, "Before SA_run\n");
  SA_run();

  SA_set_prof_file(stdout);
  SA_printFVS();

  // SA_destroy();

  return 0;
}
