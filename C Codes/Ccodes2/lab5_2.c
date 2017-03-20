#include <stdio.h>
#define MAXLINE 1000 /* maximum length of string */

/* function prototype */
void reverse(char source[], char target[]);

long string_length(char s[]){
  long nc = 0; /* Number of Charachters */
  int i = 0;
  while(s[i] != '\0'){
    if(s[i]!= '\n'||s[i]!= '\t'||s[i]!='\0'){
      nc++;
    }
    i++;
  }
  return nc;
  }
int main(void) {
  char original[] = "This is a test: can you print me in reverse character order?";
  char reversed[MAXLINE];

  printf("%s\n", original);
  reverse(original, reversed);
  printf("%s\n", reversed);
  return 0;
}

/* reverse the order of characters in 'source', write to 'target'.
   Assume 'target' is big enough. */
void reverse(char source[], char target[]) {
  int i=0;
  long strlength;
  strlength = string_length(source);
  while(source[i]!= '\0'){
    target[i] = source[strlength - i-1];
    i++;
  }
}
