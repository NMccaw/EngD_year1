#include<stdio.h>
/* TIMING CODE BEGIN (We need the following lines to take the timings.) */
#include<stdlib.h>
#include<math.h>
#include <time.h>
clock_t startm, stopm;
#define RUNS 1

#define START if ( (startm = clock()) == -1) {printf("Error calling clock");exit(1);}
#define STOP if ( (stopm = clock()) == -1) {printf("Error calling clock");exit(1);}
#define PRINTTIME printf( "%8.5f seconds used .", (((double) stopm-startm)/CLOCKS_PER_SEC/RUNS));
/* TIMING CODE END */
double f(double x){
  return sqrt(1-pow(x,2));
}

double pi(long n){

  int i;
  double a =-1.,b=1.,h,pi_n,x,s;
  h = (b-a)/n;
  s = 0.5*f(a)+0.5*f(b);
  for(i=1; i<=n-1;i++){
    x = a + i*h;
    s = s+f(x);
  }
  pi_n = s*h*2;
  return pi_n;
}
int main(void) {
    /* Declarations */

    #define n 10000000

    /* Code */
    START;               /* Timing measurement starts here */
    /* Code to be written by student, calling functions from here is fine
       if desired
    */
    printf("pi = %f",pi(n));


    STOP;                /* Timing measurement stops here */
    PRINTTIME;           /* Print timing results */
    return 0;
}
