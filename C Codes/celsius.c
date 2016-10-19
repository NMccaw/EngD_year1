/*
celsius.c - a program to print a table of temperatures in celsius T_cel
and give their conversion in farnheight T_far

*/

#include<stdio.h>

int main(void){

  int T_cel ;
  double T_far ;
  for (T_cel=-30;T_cel<30;T_cel=T_cel+2){
  T_far=(T_cel*9.0/5.0) +32;
  printf("%3d = %5.1f\n",T_cel,T_far);
  }
  return 0;
}
