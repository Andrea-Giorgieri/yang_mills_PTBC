#ifndef TIMING_C
#define TIMING_C

#include"../include/macro.h"
#include"../include/timing.h"

#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#ifdef OPENMP_MODE
#include<omp.h>
#endif

double get_wtime(void)
	{
	#ifdef OPENMP_MODE
	return omp_get_wtime();
	#else
	return (double) clock() / CLOCKS_PER_SEC;
	#endif
	}

void init_timer(Timer *const timer, char const *const name)
	{
	strcpy(timer->name, name);
	timer->count = 0;
	timer->avg_elapsed_time = 0;
	timer->max_elapsed_time = 0;
	}

void start_timer(Timer *const timer)
	{
	timer->start_time = get_wtime();
	}

void stop_timer(Timer *const timer)
	{
	timer->stop_time = get_wtime();
	timer->elapsed_time = (timer->stop_time - timer->start_time);
	if(timer->elapsed_time > timer->max_elapsed_time) timer->max_elapsed_time = timer->elapsed_time;
	timer->avg_elapsed_time = (timer->avg_elapsed_time * (double) timer->count + timer->elapsed_time) / ((double) (timer->count + 1));
	timer->count += 1;
	}

void print_timer(FILE *fp, Timer const *const timer)
	{
	if(timer->count > 1)
		fprintf(fp, "%15s : %-8.3g | %-8.3g\n", timer->name, timer->avg_elapsed_time, timer->max_elapsed_time);
	if(timer->count == 1)
		fprintf(fp, "%15s : %-8.3g \n", timer->name, timer->avg_elapsed_time);
	}

void init_time_utils(Time_Utils *const timers, double walltime)
	{
	char name[STD_STRING_LENGTH];
	strcpy(name, "Initialization");
	init_timer(&(timers->init_timer), name);
	strcpy(name, "Update");
	init_timer(&(timers->update_timer), name);
	strcpy(name, "Measure");
	init_timer(&(timers->meas_timer), name);
	strcpy(name, "Step");
	init_timer(&(timers->step_timer), name);
	strcpy(name, "Total");
	init_timer(&(timers->prog_timer), name);

	timers->wall_time = walltime;
	}

void print_time_utils(FILE *fp, Time_Utils const *const timers)
	{
	fprintf(fp, "Execution times in seconds (avg | max for repeated sections):\n");
	print_timer(fp, &(timers->init_timer));
	print_timer(fp, &(timers->update_timer));
	print_timer(fp, &(timers->meas_timer));
	print_timer(fp, &(timers->step_timer));
	print_timer(fp, &(timers->prog_timer));
	fprintf(fp, "Completed steps : %ld (%.3g steps/s)\n", timers->step_timer.count, (double) timers->step_timer.count / timers->prog_timer.avg_elapsed_time);
	}

int wall_time_check(Time_Utils const *const timers)
	{
	double const t = get_wtime();
	if((t - timers->prog_timer.start_time) + 2 * timers->step_timer.max_elapsed_time > 0.99 * timers->wall_time) return 1;
	return 0;
	}

#endif
