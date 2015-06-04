#include <stdbool.h>

#include "symop.h"

int main(int argc, char **argv)
{
  SymOp s;

  symop_new(&s, "-x", "z+2/3", " -y -1/4");
  symop_dump(&s);

  toto("-F 2 2 -1d", true,
       3.2, 3.2, 3.2, 90., 90., 90., 0.0001);

  return 0;
}
