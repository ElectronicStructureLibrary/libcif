#pragma once

#include <stdbool.h>

#include "symop.h"

typedef struct Concentration_
{
  char symbol[3];
  double value;
} Concentration;

typedef struct AtomicStructure_ AtomicStructure;

AtomicStructure* toto(const char *hallsymbol, bool primitive,
                      double a, double b, double c,
                      double alpha, double beta, double gamma,
                      double *xred, Concentration *conc, unsigned int nat,
                      double eps);

void atomicStructureDump(AtomicStructure *at);

void atomicStructureFree(AtomicStructure *at);
