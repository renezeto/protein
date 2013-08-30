#include <iostream>
using namespace std;
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct things {
  int mynum;
  int * pnum;
  char * mystring;
  string otherstring;
};

int main() {
  printf("Hello\n");
  things people;
  things * objects;
  people.mynum = 5;
  objects->mynum = 6;
  people.mystring = "this is my string";
  printf("%s",people.mystring);


return 0;
}
