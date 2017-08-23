/* A gettimeofday routine to give access to the wall
   clock timer on most UNIX-like systems. */

#include <sys/time.h>
#include <stdlib.h>

double time_of_day()
{
  struct timeval tp;
  int i;

  i = gettimeofday(&tp, NULL);
  return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

