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

if(n<6){
  return 700;
}
else{
n_fact= pow(n/2.,n);
}


return n_fact;
}
long m=19;
long factorial(long m){
  double j, fact=1;
  if(m<0){
    return -2;
  }
else{
for(j=1.;j<=m;j++){
fact=fact*j;

  }
  if(fact>upper_bound(18)){
    return -1;
  }
  else{

return fact;
}
}


}

int main(void) {
    long i;

    /* The next line should compile once "maxlong" is defined. */
    printf("maxlong()=%ld\n", maxlong());

    /* The next code block should compile once "upper_bound" is defined. */


    for (i=0; i<25; i++) {
        printf("upper_bound(%ld)=%f\n", i, upper_bound(i));
    }
    printf("factorial()=%ld\n", factorial(m));
    return 0;
}
