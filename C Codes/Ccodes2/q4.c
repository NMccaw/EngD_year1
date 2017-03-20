#include <stdio.h>
#include <math.h>
double standard_deviation(double a,double b,double c, double d){
  double N = 4.0;
  double sigma,mu;
  mu = (1.0/N)*(a+b+c+d);
  sigma = sqrt((1.0/N)*(pow((a-mu),2)+pow((b-mu),2)+pow((c-mu),2)+pow((d-mu),2)));
  return sigma;
}

int main(void) {
double a, b, c, d, stddev;
a = 16.3;
b = 24.2;
c = 733;
d = 12.27;
stddev = standard_deviation(a, b, c, d);
printf("Standard deviation of %.2f, %.2f, %.2f, %.2f is %.2f \n", a, b, c, d, stddev);
 return 0;
}
