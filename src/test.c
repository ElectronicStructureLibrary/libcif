#include <stdbool.h>

#include "symop.h"
#include "uc.h"

int main(int argc, char **argv)
{
  SymOp s;
  double xred[] = {0., 0., 0.,
                   0.5, 0.15, 0.,
                   0.5, 0.15, 0.};
  Concentration conc[] = {{" B", 1.}, {" V", 0.05}, {" N", 0.95}};

  AtomicStructure *at;

  symop_new(&s, "-x", "z+2/3", " -y -1/4");
  symop_dump(&s);

  at = toto("-R 3 2\"", false,
            3.2, 3.2, 3.2, 90., 90., 60.,
            xred, conc, 3,
            0.0001);
  atomicStructureDump(at);
  atomicStructureFree(at);

  return 0;
}
