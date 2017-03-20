#include <stdio.h>

int find_greatest(int a, int b, int c){

  if(a>b && b>c){
    return a;
  }
    else if(b>a && b>c){
      return b;
    }
      else {
        return c;
      }


}
int main(void){
int a,b,c;
a=10;
b=8;
c=1;
printf("The largest of %d, %d and %d is: %d\n", a,
b, c, find_greatest(a, b, c));
 return 0;
}
