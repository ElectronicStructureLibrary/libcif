#pragma once

#include <stdbool.h>

#include "vector3.h"
#include "matrix3.h"

typedef struct SymOp_
{
  char ops[3][32];
  Matrix3 rot;
  Vector3 trans;
} SymOp;

bool symop_new(SymOp *symop, const char *opx, const char *opy, const char *opz);

void symop_dump(SymOp *symop);
