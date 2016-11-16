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
#define N 10000000.
/* TIMING CODE END */
double f(double x) {

return sqrt(1-pow(x,2));
}
double pi(long n){
double a,b,h,s,x,pi;
int i;
  a = -1.;
  b=1.;
  h=(b-a)/n;
  s=0.5* f(a) +0.5 *f(b);
  for(i=1; i<=n-1;i++){
    x = a +i*h;
    s = s+f(x);
  }
  pi = s*h*2;
  return pi;
}
int main(void) {

double pi_v;


    /* Code */
    START;               /* Timing measurement starts here */
    pi_v = pi(N);
    printf("pi = %f",pi_v);




    STOP;                /* Timing measurement stops here */
    PRINTTIME;           /* Print timing results */
    return 0;
}
