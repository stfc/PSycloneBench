#include <stdio.h>
#include <stdint.h>

int rdtsc_available(void)
{
#if defined __INTEL_COMPILER
  return 1;
#else
  return 0;
#endif
}

/** Timer for use on Intel chips. Results are only for inter-comparison
   and not for conversion into some human measure of time.
   See http://en.wikipedia.org/wiki/Time_Stamp_Counter */
uint64_t getticks(void)
{
    uint32_t lo, hi;

    /* We can use the Intel intrinsic __rdtsc() when building with
       the Intel compiler */
#if defined __INTEL_COMPILER
    return __rdtsc();
#else
    fprintf(stderr, "TIMING: ERROR: attempting to use RDTSC timer when "
	    "not compiled with the Intel compiler\n");
    return 0;
#endif
}

/** Obtain the clock count using rdtscp instead of rdtsc. This is
  better than rdtsc since it ensures synchronisation of out-of-order
  instructions. This code provide by John D. McCalpin. */
unsigned long tacc_rdtscp(int *chip, int *core)
{
  unsigned a, d, c;

  __asm__ volatile("rdtscp" : "=a" (a), "=d" (d), "=c" (c));
  *chip = (c & 0xFFF000)>>12;
  *core = c & 0xFFF;

  return ((unsigned long)a) | (((unsigned long)d) << 32);;
}
