/*
interest.c- code which calculates the interest total interest and fraction of
original borrowing over a 24 month period
*/

#include<stdio.h>
int main(void){

  double s=1000,debt,rate,interest,total_interest,frac;
  int month;
  debt=s;
  rate=0.03;
  for(month =1 ;month<25;month++){
  interest=debt*rate;
  debt=interest+debt;
  total_interest=interest+total_interest;
  frac=(total_interest/s)*100;

  printf("month %2d: debt=%7.2f, interest=%.2f, total_interest=%7.2f, frac=%6.2f%%\n"
    ,month,debt,interest,total_interest,frac);
  }
}
return 0
