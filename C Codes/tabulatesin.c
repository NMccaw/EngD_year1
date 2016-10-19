#include<stdio.h>
#include<math.h>

# define XMIN 1.0
# define XMAX 10.
# define N 10.

int main(void){
  double x,y;
  int i;
for(i=0;i<N;i++){
  x = XMIN + (XMAX - XMIN) / (N - 1.) * i;
  y=sin(x);

  printf( "%f %f\n" ,x,y);
}



  return 0;
}
