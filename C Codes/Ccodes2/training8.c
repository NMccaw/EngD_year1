#include <stdio.h>
#include <stdlib.h>
#include <string.h>
char* mix(char *s1, char *s2){

   int N,i;
   char* mixed;

   N = strlen(s1);
   mixed = (char *)malloc(sizeof(char)*2*N);
   if (mixed == NULL) {
     printf("ERROR: Out of memory\n");
     return (char *)0;
   }
   for(i=0;i<=N;i++){
   *(mixed+0+i*2) = *(s1+i);
   *(mixed+1+i*2) = *(s2+i);
 }
   return mixed;
 }

 void use_mix(void) {
     char s1[] = "Hello World";
     char s2[] = "1234567890!";

     printf("s1 = %s\n", s1);
     printf("s2 = %s\n", s2);
     printf("r  = %s\n", mix(s1, s2));
 }

int main(void) {
   use_mix();
   return 0;
 }
