#include <stdio.h>

int main(void) {
    long i;

    /* The next line should compile once "maxlong" is defined. */
    printf("maxlong()=%ld\n", maxlong());

    /* The next code block should compile once "upper_bound" is defined. */

    /*
    for (i=0; i<10; i++) {
        printf("upper_bound(%ld)=%g\n", i, upper_bound(i));
    }
    */
    return 0;
}
