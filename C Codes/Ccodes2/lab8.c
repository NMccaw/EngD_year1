#include<stdio.h>  /* provides standard input/output tools */
#include<stdlib.h> /* provides RAND_MAX                    */

#define N 20       /* number of random numbers             */

/* Pseudo code reads:
"""procedure bubbleSort( A : list of sortable items )
  n := length(A)-1
  for(a=0; a<= n; a++)
     for(b=n; b>a; b--)
        if A[b-1] > A[b] then
           swap (A[b-1], A[b])
        end if
     end for
  end for
end procedure
"""
*/

/* Given an Array A of int, use bubble sort to sort elements in A (in
   place).*/
void bubble(int A[], int length) {
  int n,i,j;
  int temp;

  n = length - 1;
  for(i=0;i<=n;i++){
    for(j=n;j>i;j--){
      if (*(A+j-1) > *(A+j)) {
        temp = *(A+j-1);
        *(A+j-1)=*(A+j);
        *(A+j) = temp;
      }

    }

  }
}

/* Given an array of int 'a' of length 'length', print the first and
   last 'k' values */
void print_int_array(int a[], int length, int k) {
  int i;
  if (2*k < length) { /* longish array; only printing first and last k
			 values */
    for (i=0; i<k; i++)
      printf("a[%d]=%3d, ",i,a[i]);
    printf("   . . . .   ");
    for (i=length-k; i<length; i++)
      printf("a[%d]=%3d, ",i,a[i]);
  }
  else { /* for very short arrays, print all the data */
    for (i=0; i<length; i++ )
      printf("a[%d]=%3d, ",i,a[i]);
  }
  printf("\n");
}

int main(void) {
  int i;
  int data[N];
  /* initialises array with random integers between 0 and 999 */
  for (i=0; i<N; i++) {
    data[i] = (int) ((rand()+0.5)/(double) RAND_MAX * 999);
  }
  /* print data (at least beginning and end) */
  print_int_array(data,N,5);
  /* actual sorting: */
  bubble(data,N);
  printf("Data is now sorted:\n");
  /* print data (at least beginning and end) */
  print_int_array(data,N,5);

  return 0;
}
