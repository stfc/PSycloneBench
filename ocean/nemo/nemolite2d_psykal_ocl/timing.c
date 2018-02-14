#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define USE_TIMING

/* If we were supporting MPI then idPE would hold the rank of
   the current process */
#define idPE 0
#define numPEs 1

/** Maximum length of string for timer region name */
#define NAME_LEN 64
/** Maximum number of timed regions in code */
#define NUM_TIMERS 40
/** Set to 1 to switch on timing output for each individual PE */
#define PER_PROCESS_TIMING 0

/*=====================================================*/

/** Object to hold timer information */
typedef struct TimerStruct
{
  /** Name of the region being timed */
  char name[NAME_LEN];
  /** No. of times we enter this region */
  long count;
  /** Sum of time spent in the region */
  double sum_time;
}
timer;

typedef struct TimerStackStruct
{
  /** Start time for this entry on the stack */
  double start_time;
  /** Pointer to the timer associated with this entry
      on the stack */
  timer *timer_entry;
}
timer_stack_entry;

/** List of timed regions in the code */
static timer timer_list[NUM_TIMERS];

/** Position in list of next free timer */
static int next_free_timer;

/** Stack of times. Allows timer calls to be nested
    without need to pass a tag around */
static timer_stack_entry timer_stack[NUM_TIMERS];

/** Pointer to top of stack of timed regions */
static int top_of_stack;

/** Number of timed regions that we've met so far */
static int num_active_timers;

/** Local routine to output multi-process timings */
#ifdef MPI_BUILD
static void TimerReportMPI();
#endif

/** Local routine to write individual PE timings to a file */
static void WriteSingleProcessTimings(char* fname);

/*=====================================================*/

double CurrentTimeSeconds()
{
#ifdef MPI_BUILD

  return MPI_Wtime();

#else /* not MPI_BUILD */

  struct timeval tv;

  if(gettimeofday(&tv, NULL)){
    return -1.0;
  }

  return (double)(tv.tv_sec) + 1.0e-6*(double)(tv.tv_usec);
#endif /* MPI_BUILD */
}

/*=====================================================*/

void TimerInit()
{
  int i;

#ifndef USE_TIMING
  return;
#endif

  num_active_timers = 0;

  for(i=0; i<NUM_TIMERS; i++){
    timer_list[i].name[0] = '\0';
    timer_list[i].count = 0;
    timer_list[i].sum_time = 0.0;

    timer_stack[i].start_time = 0.0;
    timer_stack[i].timer_entry = NULL;
  }

  next_free_timer = 0;
}

/*=====================================================*/

void TimerStart(const char* name)
{
  int i;

#ifndef USE_TIMING
  return;
#endif

  for(i=0; i<num_active_timers; i++){
    if(strcmp(timer_list[i].name, name) == 0)break;
  }

  /* This is a timed region we've not seen before */
  if(i == num_active_timers){
    i = next_free_timer;
    next_free_timer++;
    num_active_timers++;

    /* Check that we've not run out of timers */
    /** \todo Make number of timers dynamic. */
    if(i == NUM_TIMERS){
      fprintf(stderr, "TimerStart: ERROR: max no. of timers (%d) exceeded!\n"
     	              "          : label of timer that triggered this = %s\n",
	      NUM_TIMERS, name);
      exit(1);
    }

    strncpy(timer_list[i].name, name, NAME_LEN);
  }

  /* Increment count of no. of times we visit this region */
  timer_list[i].count++;

  if(top_of_stack > NUM_TIMERS){
    fprintf(stderr, 
	    "TimerStart: ERROR: max. no. of nested timers (%d) exceeded!\n"
     	    "          : label of timer that triggered this = %s\n",
            NUM_TIMERS, name);
    exit(1);
  }

  timer_stack[top_of_stack].timer_entry = &(timer_list[i]);
  timer_stack[top_of_stack].start_time = CurrentTimeSeconds();

  top_of_stack++;
}

/*=====================================================*/

void TimerStop()
{
  double tnow;

#ifndef USE_TIMING
  return;
#endif

  tnow = CurrentTimeSeconds();
  /* Take the start time from the top of the stack */
  top_of_stack--;
  if(top_of_stack < 0){
    fprintf(stderr, 
            "TimerStop: ERROR: no matching TimerStart() call?\n");
    exit(1);
  }
  tnow -= timer_stack[top_of_stack].start_time;
#ifndef USE_TIMING
  return;
#endif

  /* Add the time taken to the cumulative total 
     for this timed region */
  (timer_stack[top_of_stack].timer_entry)->sum_time += tnow;
}

/*=====================================================*/

void TimerReport()
{
  const int name_len = 30;
  char  fname[name_len];

#ifndef USE_TIMING
  return;
#endif

#ifdef MPI_BUILD
  TimerReportMPI();

  if(PER_PROCESS_TIMING){
    snprintf(fname, name_len, "timing_report_%d.dat", idPE);
    WriteSingleProcessTimings(fname);
  }
#else
  snprintf(fname, name_len, "timing_report_serial.dat");
  WriteSingleProcessTimings(fname);
#endif

}

/*=====================================================*/

void WriteSingleProcessTimings(char* fname)
{
  FILE *fp;
  int   i;
 
  fp = fopen(fname, "w");
  if(!fp){
    fprintf(stderr, "TimerReport: failed to open file for report.\n");
    fprintf(stderr, "TimerReport: no report will be generated.\n");
    return;
  }

  fprintf(fp, 
  /*       1         1         1         1         1          1          1 */
 	  "            Region              Count    Total time      Average\n");
  fprintf(fp, 
	  "                                             (s)           (s) \n");
  fprintf(fp, 
	  "=================================================================\n");
  for(i=0; i<num_active_timers; i++){

    /* Only output timing info for regions that have been visited 
       at least once */
    if( timer_list[i].count > 0 ){
      fprintf(fp,"%-30s  %5ld  %12.4le  %12.4le\n",
	      timer_list[i].name, 
	      timer_list[i].count,
	      timer_list[i].sum_time,
	      (timer_list[i].sum_time/(double)(timer_list[i].count)) );
    }
  }
  fprintf(fp, 
	  "=================================================================\n");

  fclose(fp);
}

/*=====================================================*/

void TimerReportMPI()
{
  /* Type needed for MPI reductions that obtain main/max value and
     its location over set of PEs */
  typedef struct {
    /* Value used in reduction operation */
    float val;
    /* Rank of PE with which this value associated */
    int   rank;
  } val_loc_type;
  
  val_loc_type *timings;
  val_loc_type *min_vals;
  val_loc_type *max_vals;
  float *raw_times;
  float *sum_times;
  int i;
  FILE *fp;

  /* For each of the timed regions we want the max and min time
     spent by any of the MPI processes */

  min_vals = (val_loc_type*)malloc(num_active_timers*sizeof(val_loc_type));
  max_vals = (val_loc_type*)malloc(num_active_timers*sizeof(val_loc_type));
  timings  = (val_loc_type*)malloc(num_active_timers*sizeof(val_loc_type));
  raw_times = (float*)calloc((size_t)num_active_timers,sizeof(float));
  sum_times = (float*)calloc((size_t)num_active_timers,sizeof(float));
  if(!timings || !min_vals || !max_vals || !raw_times || !sum_times){
    fprintf(stderr, "TimerReportMPI: failed to malloc buffers for stats\n");
    return;
  }

  /* Pack the times into an array of structs */
  for(i=0; i<num_active_timers; i++){
    timings[i].val  = timer_list[i].sum_time;
    raw_times[i]    = timer_list[i].sum_time;
    timings[i].rank = idPE;
  }

  /** \todo Ideally this MPI code should be in Parallel_MPI.c but
      this currently gives a circular dependency because some
      of the routines in that file call the timing API. This 
      implies that I should split out the truly MPI-specific stuff
      from the generic parallel code (and maybe put the latter in
      e.g. Parallel.c) */
#ifdef MPI_BUILD
  MPI_Reduce(timings, max_vals, num_active_timers, 
	     MPI_FLOAT_INT, MPI_MAXLOC, 0, MPI_COMM_WORLD);
  MPI_Reduce(timings, min_vals, num_active_timers, 
	     MPI_FLOAT_INT, MPI_MINLOC, 0, MPI_COMM_WORLD);
  MPI_Reduce(raw_times, sum_times, num_active_timers, 
	     MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
#endif

  /* Only PE 0 writes the timing report */
  if(idPE == 0){

    fp = fopen("timing_report_parallel.dat","w");
    if(!fp){
      fprintf(stderr, 
	      "TimerReportMPI: failed to open file to write timing data\n");
      return;
    }

    fprintf(fp, 
	    "                                                    Mean time in region (s)\n");
    fprintf(fp,
  	    "            Region              Count       Min[rank]          Mean         Max[rank]\n");
    fprintf(fp, 
	    "=========================================================================================\n");
    for(i=0; i<num_active_timers; i++){
      fprintf(fp, "%-30s  %5ld  %12.4e[%4d] %12.4e %12.4e[%4d]\n", 
              timer_list[i].name,
	      timer_list[i].count,
              min_vals[i].val/(float)(timer_list[i].count), 
	      min_vals[i].rank,
              sum_times[i]/(float)(numPEs*timer_list[i].count),
              max_vals[i].val/(float)(timer_list[i].count),
	      max_vals[i].rank);
    }
    fprintf(fp, 
	    "=========================================================================================\n");

    fclose(fp);
  }

  /* Free memory */
  free(min_vals);  min_vals=NULL;
  free(max_vals);  max_vals=NULL;
  free(timings);   timings=NULL;
  free(raw_times); raw_times=NULL;
  free(sum_times); sum_times=NULL;
}
