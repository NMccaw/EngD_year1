#include <stdio.h>
#include<limits.h>
#include<float.h>
#include<math.h>
long maxlong(void){
  /*LONG_MAX = %12ld\n, LONG_MAX);*/
  return LONG_MAX;
}
double upper_bound(long n){
double n_fact;

if(n<5){
  return 700;
}
else{
n_fact= pow(n/2.,n);
}


return n_fact;
}


int main(void) {
    long i;

    /* The next line should compile once "maxlong" is defined. */
    printf("maxlong()=%ld\n", maxlong());

    /* The next code block should compile once "upper_bound" is defined. */


    for (i=0; i<10; i++) {
        printf("upper_bound(%ld)=%f\n", i, upper_bound(i));
    }

    return 0;
}
