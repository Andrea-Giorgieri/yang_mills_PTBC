#ifndef RANDOM_H
#define RANDOM_H

double casuale(void);                    // random number in (0,1)

void initrand(unsigned int s);           // initialize random generator

double gauss1(void);                     // normal gaussian random number

void gauss2(double *res1, double *res2); // normal gaussian random numbers

#endif
