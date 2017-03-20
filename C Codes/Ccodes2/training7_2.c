#include<stdio.h>
#include<stdlib.h>

long*  make_long_array(long n){
  long *pointer;
  pointer = (long *)malloc(sizeof(long)*n);
  if(pointer == NULL){
    printf("Error: Out of Memory");
    return NULL;
  }
  return pointer;
}

void use_make_long_array(long n) {
  int i;
  long *p = make_long_array(n);

  printf("In use_make_long_array(%ld)\n", n);

  /* if p is not NULL, we could allocate memory, and we proceed
     with testing: */
  if (p != NULL) {

    /* write some data to the array -- if the allocated memory
       is too short, this might trigger a segfault */
    for (i=0; i<n; i++) {
        p[i]=i+42;           /* just write some data */
    }

    /* free array -- if the allocated array is too short, we may
       have corrupted malloc/free metadata when writing the i+42 data
       above, and this may show when we call the free command: */

    free(p);

    /* if the program does not crash, it is a good sign [but
       no proof for correctness]. The other way round: if the
       program segfaults or crashes, this is not a good sign. */
  }
  else { /* we get here if memory allocation didn't work for
            some reason. */
    printf("Error - it seems that the memory allocation failed.\n");
  }
}

int main(void) {
  int n;
  for (n=0; n<20; n++) {
    use_make_long_array(n);
  }
  return 0;
}
