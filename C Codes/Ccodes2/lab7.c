#include <stdio.h>
#include <stdlib.h>

long* make_fib_array(long n){

  long * pointer;
  int i;

  pointer = (long *)malloc(sizeof(long)*n);
  if (pointer == NULL) {
    printf("ERROR: Out of memory\n");
    return (long *)0;
  }


  *(pointer+0) = 0;
  *(pointer+1) = 1;
  for(i=2;i<=n;i++){
    *(pointer+i) = *(pointer+i-1) + *(pointer+i-2);
  }
  return pointer;
}

void use_fib_array(long N) {
  /* N is the maximum number for fibarray length */
  long n;      /* counter for fibarray length */
  long i;      /* counter for printing all elements of fibarray */
  long *fibarray;  /* pointer to long -- pointer to the fibarray itself*/

  /* Print one line for each fibarray length n*/
  for (n=2; n<=N; n++) {
    /* Obtain an array of longs with data */
    fibarray = make_fib_array(n);

    /* Print all elements in array */
    printf("fib(%2ld) : [",n);
    for (i=0; i<n; i++) {
      printf(" %ld", fibarray[i]);
    }
    printf(" ]\n");

    /* free array memory */
    free(fibarray);
  }

}

int main(void) {
  use_fib_array(10);
  return 0;
}
