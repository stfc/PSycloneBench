#ifndef _TIMING_INCLUDE
#define _TIMING_INCLUDE

/** Returns the current time in seconds - for use in calculating 
    execution time */
#ifdef __cplusplus
extern "C" 
#endif
double CurrentTimeSeconds();

void TimerInit();

void TimerStart(const char *name);

void TimerStop();

void TimerReport();

#endif
