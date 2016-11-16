#include <stdio.h>
#define MAXLINE 1000 /* maximum length of string */

/* function prototype */
void reverse(char source[], char target[]);
long length;
int i;
int reversal;
int main(void) {
  char original[] = "This is a test: can you print me in reverse character order?";
  char reversed[MAXLINE];

  printf("%s\n", original);
  reverse(original, reversed);
  printf("%s\n", reversed);
  return 0;
}

long string_length(char s[]){
length=0;
for(i=0;s[i]!='\0';i++){
  length++;
}
return length;
}
/* reverse the order of characters in 'source', write to 'target'.
   Assume 'target' is big enough. */
void reverse(char source[], char target[]) {

length=string_length(source);
for(i=0; source[i]!='\0'; i++){
reversal=(length-1)-i;
target[i]=source[reversal];
}
return;
}
