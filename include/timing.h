#ifndef TIMING_H
#define TIMING_H

#include"macro.h"

#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#ifdef OPENMP_MODE
#include<omp.h>
#endif

typedef struct Timer {
	char name[STD_STRING_LENGTH];
	double start_time;
	double stop_time;
	double elapsed_time;
	long count;
	double avg_elapsed_time;
	double max_elapsed_time;
	}	Timer;

typedef struct Time_Utils {
	Timer init_timer;
	Timer update_timer;
	Timer meas_timer;
	Timer step_timer;
	Timer prog_timer;
	double wall_time;
	}	Time_Utils;

double get_wtime();
void init_timer(Timer * const timer, char const * const name);
void start_timer(Timer * const timer);
void stop_timer(Timer * const timer);
void print_timer(FILE *fp, Timer const * const timer);

void init_time_utils(Time_Utils * const timers, double walltime);
void print_time_utils(FILE *fp, Time_Utils const * const timers);
int wall_time_check(Time_Utils const * const timers);


#endif
