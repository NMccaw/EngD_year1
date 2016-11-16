#include <stdio.h>
#include <string.h>
/* Function void rstrip(char s[])
modifies the string s: if at the end of the string s there are one or more spaces,
then remove these from the string.

The name rstrip stands for Right STRIP, trying to indicate that spaces at the 'right'
end of the string should be removed.
*/

void rstrip(char s[]) {
    int i;
    for ( i = 0; i < strlen(s); i++) {
      if (s[i] == ' ' && s[i+1]== ' '){
        s[i]='\0';
      }
    }
}


int main(void) {
  char test1[] = "la la la la       ";

  printf("Original string reads  : |%s|\n", test1);
  rstrip(test1);
  printf("r-stripped string reads: |%s|\n", test1);

  return 0;
}
