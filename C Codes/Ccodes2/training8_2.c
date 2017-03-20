#include<stdio.h>
#include<stdlib.h>
#include <string.h>

char* mix(char *s1, char *s2){
  int i;
  char *r;

  r = (char *)malloc(sizeof(char)*2*strlen(s1));

  for(i = 0;*(s1+i)!='\0';i++){
    *(r+(2*i))=*(s1+i);
    *(r+(2*i)+1) = *(s2+i);
  }
  return r;
}

void use_mix(void) {
    char s1[] = "Hello World";
    char s2[] = "1234567890!";

    printf("s1 = %s\n", s1);
    printf("s2 = %s\n", s2);
    printf("r  = %s\n", mix(s1, s2));
}

int main(void){
  use_mix();
}
