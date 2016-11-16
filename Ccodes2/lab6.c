#include <stdio.h>
#include <string.h>
/* Function void rstrip(char s[])
modifies the string s: if at the end of the string s there are one or more spaces,
then remove these from the string.

The name rstrip stands for Right STRIP, trying to indicate that spaces at the 'right'
end of the string should be removed.
*/

void lstrip(char s[]) {
  int i,nospace, string_len,wspace;
  if(s[0]==' ' && s[1] !=' '){
    wspace=1;
  }
  else{
    wspace = 0;
  }
  string_len= strlen(s);
  for(i=0;i < string_len;i++){
    if((s[i] == ' ' && s[i+1] == ' ') || (s[i] == ' ' && s[i+1] != ' ' && s[i-1]==' ') ){
      ++wspace;
}
}


for(i=0;i <= string_len;i++){
   if(s[i] != ' ' || s[i] != '\n' || s[i] != '\t' || s[i] != '\0'){

     nospace = i;
     s[nospace-(wspace)] = s[nospace];
     /*s[string_len-wspace] = '\0';*/
   }
  }
}



int main(void) {
  char test1[] = " Hello World this is a test";

  printf("Original string reads  : |%s|\n", test1);
  lstrip(test1);
  printf("l-stripped string reads: |%s|\n", test1);

  return 0;
}
