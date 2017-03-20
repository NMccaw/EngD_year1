#include <stdio.h>

/* Function void rstrip(char s[])
modifies the string s: if at the end of the string s there are one or more spaces,
then remove these from the string.

The name rstrip stands for Right STRIP, trying to indicate that spaces at the 'right'
end of the string should be removed.
*/

void rstrip(char s[]) {
    int i=0;
    while(s[i] != '\0'){
      if((s[i] == ' ' && s[i+1] ==' ') || (s[i] == ' ' && s[i+1] == '\0')){
        s[i] = '\0';
      }
      i++;
    }
}
void lstrip(char s[]) {
int i=0,j=0;
while (s[i] !='\0') {
  if((s[i] != ' ' && s[i+1] != ' ')){
    s[j] = s[i];
    j++;
    }
  i++;
}
s[j]='\0';
}


int main(void) {
  char test1[] = "     lalalala";

  printf("Original string reads  : |%s|\n", test1);
  lstrip(test1);
  printf("l-stripped string reads: |%s|\n", test1);

  return 0;
}
