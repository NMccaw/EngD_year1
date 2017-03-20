#include <stdio.h>

/* function prototype */
long length;

int i;
long string_length(char s[]){
length=0;

for(i=0;s[i]!='\0';i++){




  length++;

}

return length;
}

int main(void) {
  char s1[]="Hello";
  char s2[]="x";
  char s3[]="line 1\tline 2\n";

  printf("%20s | %s\n", "string_length(s)", "s");
  printf("%20ld | %s\n", string_length(s1), s1);
  printf("%20ld | %s\n", string_length(s2), s2);
  printf("%20ld | %s\n", string_length(s3), s3);
  return 0;
}
