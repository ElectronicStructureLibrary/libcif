#include "symop.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include <math.h>

/* 'C 2y' : [['x',' y',' z'], */
/*           ['-x',' y',' -z'], */
/*           ['x+1/2',' y+1/2',' z'], */
/*           ['-x+1/2',' y+1/2',' -z']], */

static double dir2v(const char *op, char dir)
{
  char *pt;

  pt = strchr(op, dir);
  if (!pt)
    return 0.;

  /* Look for the sign. */
  for (; pt != op && *pt != '+' && *pt != '-'; pt--);
  if (*pt == '+') return +1.;
  else if (*pt == '-') return -1.;
  else return +1.;
}

static double tr2v(const char *op)
{
  char *pt;
  int num, denom;

  pt = strchr(op, '/');
  if (!pt)
    return 0.;

  denom = atoi(pt + 1);
  for (; pt != op && *pt != '+' && *pt != '-'; pt--);
  num = atoi(pt);

  return (double)num / (double)denom;
}

bool symop_new(SymOp *symop, const char *opx, const char *opy, const char *opz)
{
  strncpy(symop->ops[0], opx, sizeof(symop->ops[0]) / sizeof(char));
  strncpy(symop->ops[1], opy, sizeof(symop->ops[1]) / sizeof(char));
  strncpy(symop->ops[2], opz, sizeof(symop->ops[2]) / sizeof(char));
  /*
   * Return a rotation matrix from "x,y,z" representation of a symmetry operation
   * !!!With respect to cartesian axes!!!
   */
  m3_new(&symop->rot,
         dir2v(opx, 'x'), dir2v(opy, 'x'), dir2v(opz, 'x'),
         dir2v(opx, 'y'), dir2v(opy, 'y'), dir2v(opz, 'y'),
         dir2v(opx, 'z'), dir2v(opy, 'z'), dir2v(opz, 'z'));
  /* Return a translation vector from "x,y,z" representation of a symmetry operation */
  v3_new(&symop->trans,
         tr2v(opx), tr2v(opy), tr2v(opz));

  return fabs(m3_det(&symop->rot)) == 1.;
}

void symop_dump(SymOp *symop)
{
  fprintf(stdout, "symmetry operator: [%s, %s, %s]\n",
          symop->ops[0], symop->ops[1], symop->ops[2]);
  fprintf(stdout,
          "rotation: [%g, %g, %g\n"
          "           %g, %g, %g\n"
          "           %g, %g, %g]\n",
          symop->rot.x[0][0], symop->rot.x[0][1], symop->rot.x[0][2],
          symop->rot.x[1][0], symop->rot.x[1][1], symop->rot.x[1][2],
          symop->rot.x[2][0], symop->rot.x[2][1], symop->rot.x[2][2]);          
  fprintf(stdout,
          "translation: [%g, %g, %g]\n",
          symop->trans.x[0], symop->trans.x[1], symop->trans.x[2]);
}
