#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include <math.h>
#ifndef M_PI
#define M_PI acos(-1.0)
#endif

#include "uc.h"
#include "vector3.h"
#include "matrix3.h"

#include "hall_data.h"

/* Round value for close number to specific numbers. */
static const double magicVals[] =
  {1./3., 2./3., 1./2., 1./4., 1., 0., 1.4142135623730951,
   1./6., 5./6., 1.7320508075688772, 1.7320508075688772/2.};
static double roundToMagic(double x, double eps)
{
  int i;

  for (i = 0; i < sizeof(magicVals) / sizeof(double); i++)
    if (fabs(fabs(x) - magicVals[i]) <= eps)
      return (x < 0.) ? -magicVals[i] : magicVals[i];
      
  return x;
}

typedef enum {
  TRICLINIC,
  MONOCLINIC,
  ORTHORHOMBIC,
  TETRAGONAL,
  TRIGONAL,
  HEXAGONAL,
  CUBIC,
  N_CRYSTALSYSTEM
} CrystalSystem;

CrystalSystem crystal_system(unsigned int spacegroup)
{
  /* Determine crystal system */
  if (spacegroup == 0)
    return N_CRYSTALSYSTEM;
  else if (spacegroup <= 2)
    return TRICLINIC;
  else if (spacegroup <= 15)
    return MONOCLINIC;
  else if (spacegroup <= 74)
    return ORTHORHOMBIC;
  else if (spacegroup <= 142)
    return TETRAGONAL;
  else if (spacegroup <= 167)
    return TRIGONAL;
  else if (spacegroup <= 194)
    return HEXAGONAL;
  else if (spacegroup <= 230)
    return CUBIC;
  return N_CRYSTALSYSTEM;
}

Matrix3 conventionalCell(CrystalSystem crys,
                         double a, double b, double c,
                         double alpha, double beta, double gamma,
                         double eps)
{
  double coa, boa, alphar, betar, gammar;
  double angfac1, angfac2;
  Matrix3 cell;
  
  /* Set up Bravais lattice vectors of the conventional cell */
  coa = c / a;
  boa = b / a;
  alphar = alpha * M_PI / 180.;
  betar  = beta * M_PI / 180.;
  gammar = gamma * M_PI / 180.;

  switch (crys)
    {
    case (CUBIC):
      m3_new(&cell,
             1.0, 0.0, 0.0,
             0.0, 1.0, 0.0,
             0.0, 0.0, 1.0);
      break;
    case (HEXAGONAL):
      m3_new(&cell,
             sin(gammar), cos(gammar), 0.0,
             0.0, 1.0, 0.0,
             0.0, 0.0, coa);
      break;
    case (TETRAGONAL):
    case (ORTHORHOMBIC):
      m3_new(&cell,
             1.0, 0.0, 0.0,
             0.0, boa, 0.0,
             0.0, 0.0, coa);
      break;
    case (TRIGONAL):
      /* Hexagonal cell taken as conventional */
      if (fabs(gamma - 120.) < eps)
        gammar = 120. * M_PI / 180.;
      m3_new(&cell,
             sin(gammar), cos(gammar), 0.0,
             0.0, 1.0, 0.0,
             0.0, 0.0, coa);
      break;
    case (TRICLINIC):
    case (MONOCLINIC):
    case (N_CRYSTALSYSTEM):
      angfac1 = (cos(alphar) - cos(betar)*cos(gammar))/sin(gammar);
      angfac2 = sqrt(sin(gammar)*sin(gammar) - cos(betar)*cos(betar) -
                     cos(alphar)*cos(alphar) +
                     2*cos(alphar)*cos(betar)*cos(gammar)) / sin(gammar);
      m3_new(&cell,
             1.0, 0.0, 0.0,
             boa * cos(gammar), boa * sin(gammar), 0.0,
             coa * cos(betar), coa * angfac1, coa * angfac2);
      break;
    }
  return cell;
}

/**
 * rhomb2hex: a key / value matching between hexagonal setting of a space
 * group in Hall notation and corresponding rhombohedral setting.
 */
static const char *rhomb2hex[] = {"P 3*",     "R 3",
                                  "-P 3*",    "-R 3",
                                  "P 3* 2",   "R 3 2\"",
                                  "P 3* -2",  "R 3 -2\"",
                                  "P 3* -2n", "R 3 -2\"c",
                                  "-P 3* 2",  "-R 3 2\"",
                                  "-P 3* 2n", "-R 3 2\"c",
                                  NULL};
static bool hallIsRhombOrHex(const char *hall)
{
  int i;

  for (i = 0; rhomb2hex[i]; i++)
    if (!strcmp(hall, rhomb2hex[i]))
      return true;
  return false;
}

/**
 * Provide the spacegroupsetting from Hall symbol.
 */
static char sgSettingFromHall(const char *hallsymbol)
{
  return (hallsymbol[0] == '-') ? hallsymbol[1] : hallsymbol[0];
}

/**
 * Provide spacegroup from Hall symbol.
 */
static unsigned int sgFromHall(const char *hallsymbol)
{
  unsigned int i, j;

  for (i = 0; sg2Hall[i][0]; i++)
    for (j = 0; j < MAX_HALL_PER_SG && sg2Hall[i][j]; j++)
      if (!strcmp(hallsymbol, sg2Hall[i][j]))
        return i;
  return 0;
}

/**
 *  Choices of lattice vectors made to largely coincide with the choices
 *  made at http://cst-www.nrl.navy.mil/lattice/
 * 
 *  The induced translation vectors are from Zachariasen, "Theory of x-ray
 *  diffraction in crystals". 
 * 
 *  Relations between rhombohedral and hexagonal settings of trigonal
 *  space groups from Tilley, "Crystals and crystal structures"
 *  Bravais lattice vectors:
 * 
 *  a_r =  2/3 a_h + 1/3 b_h + 1/3 c_h
 *  b_r = -1/3 a_h + 1/3 b_h + 1/3 c_h
 *  c_r = -1/3 a_h - 2/3 b_h + 1/3 c_h
 * 
 *  a_h = a_r - b_r
 *  b_h = b_r - c_r
 *  c_h = a_r + b_r + c_r
 * 
 *  a, c and rhombohedral angle alpha:
 * 
 *  a_h = 2 * a_r * sin(alpha/2)
 *  c_h = a_r * sqrt(3 + 6*cos(alpha))
 *  
 *  a_r = sqrt(3*a_h^2 + c_h^2) / 3
 *  sin(alpha/2) = 3*a_h / (2 * sqrt(3*a_h^2 + c_h^2))
 * 
 * @param spacegroupsetting: 
 */
static void getLatticeTranslations(const char *hallsymbol,
                                   double gamma,
                                   double eps,
                                   Matrix3 *lattrans,
                                   Vector3 transvecs[4],
                                   size_t *tlen)
{
  switch (sgSettingFromHall(hallsymbol))
    {
    case ('I'):
      /* Body centered */
      *tlen = 2;
      v3_new(transvecs,     0.0, 0.0, 0.0);
      v3_new(transvecs + 1, 0.5, 0.5, 0.5);
      if (crystal_system(sgFromHall(hallsymbol)) == CUBIC)
        m3_new(lattrans,
                    -0.5,  0.5,  0.5,
                    0.5, -0.5,  0.5,
                    0.5,  0.5, -0.5);
      else
        m3_new(lattrans,
                    1.0,  0.0,  0.0,
                    0.0,  1.0,  0.0,
                    0.5,  0.5,  0.5);
      break;
    case ('F'):
      /* Face centered */
      *tlen = 4;
      v3_new(transvecs,     0.0, 0.0, 0.0);
      v3_new(transvecs + 1, 0.5, 0.5, 0.0);
      v3_new(transvecs + 2, 0.5, 0.0, 0.5);
      v3_new(transvecs + 3, 0.0, 0.5, 0.5);
      m3_new(lattrans,
                  0.5,  0.5,  0.0,
                  0.5,  0.0,  0.5,
                  0.0,  0.5,  0.5);
      break;
    case ('A'):
      /* A-centered */
      *tlen = 2;
      v3_new(transvecs,     0.0, 0.0, 0.0);
      v3_new(transvecs + 1, 0.0, 0.5, 0.5);
      m3_new(lattrans,
                  1.0,  0.0,  0.0,
                  0.0,  0.5, -0.5,
                  0.0,  0.5,  0.5);
      break;
    case ('B'):
      /* B-centered */
      *tlen = 2;
      v3_new(transvecs,     0.0, 0.0, 0.0);
      v3_new(transvecs + 1, 0.5, 0.0, 0.5);
      m3_new(lattrans,
                  0.5,  0.0, -0.5,
                  0.0,  1.0,  0.0,
                  0.5,  0.0,  0.5);
      break;
    case ('C'):
      /* C-centered */
      *tlen = 2;
      v3_new(transvecs,     0.0, 0.0, 0.0);
      v3_new(transvecs + 1, 0.5, 0.5, 0.0);
      m3_new(lattrans,
                  0.5, -0.5,  0.0,
                  0.5,  0.5,  0.0,
                  0.0,  0.0,  1.0);
      break;
    default:
      if (hallIsRhombOrHex(hallsymbol) && fabs(gamma - 120.) < eps)
        {
          /* Rhombohedral from hexagonal setting */
          *tlen = 3;
          v3_new(transvecs,     0.0, 0.0, 0.0);
          v3_new(transvecs + 1, 1./3., 2./3., 2./3.);
          v3_new(transvecs + 2, 2./3., 1./3., 1./3.);
          m3_new(lattrans,
                      2./3.,   1./3., 1./3.,
                      -1./3.,  1./3., 1./3.,
                      -1./3., -2./3., 1./3.);
        }
      else
        {
          *tlen = 1;
          v3_new(transvecs, 0.0, 0.0, 0.0);
          m3_new(lattrans,
                      1, 0, 0,
                      0, 1, 0,
                      0, 0, 1);
        }
    }
}

static void transformToPrimitive(Matrix3 *latticevectors,
                                 Matrix3 *lattrans, double eps)
{
  Matrix3 tmpmat;
  int i, j;
  
  /* Transform to primitive cell */
  tmpmat = m3_mm(latticevectors, lattrans);
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      latticevectors->x[i][j] = roundToMagic(tmpmat.x[i][j], eps);
}

static void getSymOpFromHall(const char *hallsymbol,
                             SymOp symops[MAX_HALL_SYMOPS],
                             size_t *nsymops)
{
  int i, j;

  for (i = 0; i < sizeof(hallSymOps) / sizeof(struct HallSymOp_); i++)
    if (!strcmp(hallSymOps[i].hallsymbol, hallsymbol))
      {
        for (j = 0; j < MAX_HALL_SYMOPS && hallSymOps[i].hallSymOps[j][0]; j++)
          symop_new(symops + j,
                    hallSymOps[i].hallSymOps[j][0],
                    hallSymOps[i].hallSymOps[j][1],
                    hallSymOps[i].hallSymOps[j][2]);
        *nsymops = j;
        return;
      }
  *nsymops = 0;
}

/**
 * Generate a mask of non-redundant symmetry operation with respect
 * to a set of lattice constant translations.
 */
static void maskSymOp(SymOp *symops, bool *mask, size_t nsymops,
                      Vector3 *transvects, size_t nvects, double eps)
{
  int i, j, t;
  Vector3 v;

  for (i = 0; i < nsymops; i++) mask[i] = true;

  /* Assume first translation vector being zero always. */
  if (nvects == 1) return;
  
  for (i = 0; i < nsymops; i++)
    for (j = 0; j < nsymops; j++)
      for (t = 1; t < nvects && mask[j]; t++)
        {
          v = v3_sum(&symops[i].trans, transvects + t);
          if (v3_eq(&v, &symops[j].trans, eps) &&
              m3_eq(&symops[i].rot, &symops[j].rot, eps) &&
              v3_len(&symops[i].trans) < v3_len(&symops[j].trans))
            mask[j] = false;
        }
}

static void fromSymOpToCartesians(SymOp *symops, bool *mask, size_t nsymops,
                                  Matrix3 *convCell, Matrix3 *lattrans)
{
  Matrix3 inv, tmp, inv2;
  int i;

  inv = m3_inv(convCell);
  inv2 = m3_inv(lattrans);
  for (i = 0; i < nsymops; i++)
    if (mask[i])
      {
        tmp = m3_mm(&symops[i].rot, convCell);
        symops[i].rot = m3_mm(&inv, &tmp);
        symops[i].trans = m3_mv(&inv2, &symops[i].trans);
      }
}

/**
 * Test that the lattice vectors are invariant under all space group operations
 * If not, the data is given in some non-standard representation that presently
 * can't be handled.
 */
static bool checkSymOp(SymOp *symops, bool *mask, size_t nsymops,
                       CrystalSystem crys, Matrix3 *convCell, double eps)
{
  int i, j;
  Vector3 v1, vm1, vrot, trans, lv[3], lvm[3];
  
  switch (crys)
    {
    case (HEXAGONAL):
    case (TRIGONAL):
      /* Hexagonal and trigonal as a special case...
         check that the hexagonal planes are in the ab plane */
      v3_new(&v1, 0., 0., 1.);
      v3_new(&vm1, 0., 0., -1.);
      for (i = 0; i < nsymops; i++)
        if (mask[i])
          {
            v3_from_array(&vrot, symops[i].rot.x[2]);
            if (!v3_eq(&vrot, &v1, eps) && !v3_eq(&vrot, &vm1, eps))
              return false;
          }
      return true;
    default:
      v3_from_array(lv + 0, convCell->x[0]);
      v3_from_array(lv + 1, convCell->x[1]);
      v3_from_array(lv + 2, convCell->x[2]);
      v3_new(lvm + 0, -convCell->x[0][0], -convCell->x[0][1], -convCell->x[0][2]);
      v3_new(lvm + 1, -convCell->x[1][0], -convCell->x[1][1], -convCell->x[1][2]);
      v3_new(lvm + 2, -convCell->x[2][0], -convCell->x[2][1], -convCell->x[2][2]);
      for (i = 0; i < nsymops; i++)
        if (mask[i] && v3_len(&symops[i].trans) < eps)
            for (j = 0; j < 3; j++)
              {
                v3_from_array(&vrot, lv[j].x);
                trans = m3_mv(&symops[i].rot, &vrot);
                if (!v3_eq(&trans, lv + 0, eps) &&
                    !v3_eq(&trans, lv + 1, eps) &&
                    !v3_eq(&trans, lv + 2, eps) &&
                    !v3_eq(&trans, lvm + 0, eps) &&
                    !v3_eq(&trans, lvm + 1, eps) &&
                    !v3_eq(&trans, lvm + 2, eps))
                  return false;
              }
    }
  return true;
}

typedef struct Node_
{
  Vector3 position;
  Concentration *species;
  unsigned int nspecies;
  char *label;
  bool mask;
} Node;

struct AtomicStructure_
{
  unsigned int nineq, nat;
  Node **ineq;
  unsigned int *ndup;
};

static AtomicStructure* generateCell(SymOp *symops, bool *mask, size_t nsymops,
                                     Vector3 *tvecs, unsigned int ntvec,
                                     double *xred, Concentration *conc, unsigned int nat,
                                     double eps)
{
  Node *ineq;
  int i, j, k, si, sj, n, nsym;
  AtomicStructure *at;
  Vector3 pos;
  bool add;

  /* Count symmetries for later. */
  nsym = 0;
  for (j = 0; j < nsymops; j++)
    if (mask[j])
      nsym += 1;

  /* Generate an array of Node based on given atoms. */
  ineq = malloc(sizeof(Node) * nat);
  for (i = 0; i < nat; i++)
    {
      v3_from_array(&ineq[i].position, xred + 3 * i);
      ineq[i].mask = true;
      ineq[i].nspecies = 1;
      ineq[i].species = malloc(sizeof(Concentration));
      ineq[i].species[0] = conc[i];
    }
  /* Clean up the given array in case of duplicates. */
  for (i = 0; i < nat; i++)
    for (j = nat -1; j > i; j--)
      if (v3_eq(&ineq[i].position, &ineq[j].position, eps))
        {
          /* Mixup species of j in i and mask j. */
          for (si = 0; si < ineq[i].nspecies; i++)
            if (!strcmp(ineq[i].species[si].symbol,
                        ineq[j].species[0].symbol))
              {
                ineq[i].species[si].value += ineq[j].species[0].value;
                break;
              }
          if (si == ineq[i].nspecies)
            {
              ineq[i].species = realloc(ineq[i].species, sizeof(Concentration) * (ineq[i].nspecies + 1));
              ineq[i].species[si] = ineq[j].species[0];
              ineq[i].nspecies += 1;
            }
          ineq[j].mask = false;
        }
  /* Count atoms and generate at. */
  n = 0;
  for (i = 0; i < nat; i++)
    if (ineq[i].mask) n += 1;
  at = malloc(sizeof(AtomicStructure));
  at->nat = n;
  at->nineq = n;
  at->ineq = malloc(sizeof(Node*) * n);
  at->ndup = malloc(sizeof(unsigned int) * n);
  /* Steal true inequivalent atoms, and clean the others. */
  for (i = 0, j = 0; i < nat; i++)
    if (ineq[i].mask)
      {
        at->ndup[j] = 1;
        /* Preallocate everyone. */
        at->ineq[j] = malloc(sizeof(Node) * (1 + nsym));
        at->ineq[j][0] = ineq[i];
        j += 1;
      }
    else
      free(ineq[i].species);
  free(ineq);
  /* Duplicate all atoms. */
  for (i = 0; i < at->nineq; i++)
    {
      for (j = 0; j < nsymops; j++)
        if (mask[j])
          {
            at->ineq[i][at->ndup[i]].position =
              symop_at(symops + j, &at->ineq[i][0].position);
            /* Only use this position if not equivalent by translation. */
            add = true;
            for (si = 0; si < at->ndup[i] && add; si++)
              for (k = 0; k < ntvec && add; k++)
                {
                  pos = v3_sum(&at->ineq[i][at->ndup[i]].position, tvecs + k);
                  add = !v3_eq(&at->ineq[i][si].position, &pos, eps);
                }
            if (add)
              {
                /* Shallow copy here, only first node store concentrations. */
                at->ineq[i][at->ndup[i]] = at->ineq[i][0];
                at->ndup[i] += 1;
              }
          }
      at->ineq[i] = realloc(at->ineq[i], sizeof(Node) * at->ndup[i]);
    }
  /* Mark duplicated nodes in full box. */
  for (i = 0; i < at->nineq; i++)
    for (si = 0; si < at->ndup[i]; si++)
      for (j = 0; j < at->nineq; j++)
        for (sj = 0; sj < at->ndup[j]; sj++)
          if (v3_eq(&at->ineq[i][si].position, &at->ineq[j][sj].position, eps))
            at->ineq[j][sj].mask = (i == j && si == sj);
  /* Count atoms. */
  at->nat = 0;
  for (i = 0; i < at->nineq; i++)
    for (si = 0; si < at->ndup[i]; si++)
      if (at->ineq[i][si].mask)
        at->nat += 1;

  return at;
}

void atomicStructureFree(AtomicStructure *at)
{
  int i;

  for (i = 0; i < at->nineq; i++)
    {
      free(at->ineq[i][0].species);
      /* for (j = 0; j < at->ndup[i]; j++) */
      /*   { */
      /*     free(at->ineq[i][j].label); */
      /*   } */
      free(at->ineq[i]);
    }
  free(at->ndup);
  free(at->ineq);
  free(at);
}

void atomicStructureDump(AtomicStructure *at)
{
  int i, si;
  
  fprintf(stdout, "atomic structure:\n");
  fprintf(stdout, "  number of ineq atoms: %d\n", at->nineq);
  fprintf(stdout, "  number of total atoms: %d\n", at->nat);
  fprintf(stdout, "  number of atoms per ineq:\n");
  for (i = 0; i < at->nineq; i++)
    fprintf(stdout, "  - %d\n", at->ndup[i]);
  fprintf(stdout, "  positions:\n");
  for (i = 0; i < at->nineq; i++)
    for (si = 0; si < at->ndup[i]; si++)
      if (at->ineq[i][si].mask)
        {
          fprintf(stdout,
                  "  - [%g, %g, %g]\n",
                  at->ineq[i][si].position.x[0],
                  at->ineq[i][si].position.x[1],
                  at->ineq[i][si].position.x[2]);
        }
}

AtomicStructure* toto(const char *hallsymbol, bool primitive,
                      double a, double b, double c,
                      double alpha, double beta, double gamma,
                      double *xred, Concentration *conc, unsigned int nat,
                      double eps)
{
  unsigned int sgnr;
  CrystalSystem crys;
  Matrix3 convCell, primCell;
  Matrix3 lattrans;
  Vector3 transvecs[4];
  size_t ntrans, nsymops;
  SymOp symops[MAX_HALL_SYMOPS];
  bool mask[MAX_HALL_SYMOPS];
  AtomicStructure *at;

  sgnr = sgFromHall(hallsymbol);
  crys = crystal_system(sgnr);

  /* Generate conventional and primitive cells. */
  primCell = convCell = conventionalCell(crys, a, b, c, alpha, beta, gamma, eps);
  fprintf(stdout,
          "conventional cell: [%g, %g, %g\n"
          "                 %g, %g, %g\n"
          "                 %g, %g, %g]\n",
          primCell.x[0][0], primCell.x[0][1], primCell.x[0][2],
          primCell.x[1][0], primCell.x[1][1], primCell.x[1][2],
          primCell.x[2][0], primCell.x[2][1], primCell.x[2][2]);          
  
  if (primitive)
    {
      getLatticeTranslations(hallsymbol,
                             gamma,
                             eps,
                             &lattrans,
                             transvecs,
                             &ntrans);
      transformToPrimitive(&primCell, &lattrans, eps);
      fprintf(stdout,
              "primitive cell: [%g, %g, %g\n"
              "                 %g, %g, %g\n"
              "                 %g, %g, %g]\n",
              primCell.x[0][0], primCell.x[0][1], primCell.x[0][2],
              primCell.x[1][0], primCell.x[1][1], primCell.x[1][2],
              primCell.x[2][0], primCell.x[2][1], primCell.x[2][2]);          
    }
  else
    {
      ntrans = 1;
      v3_new(transvecs, 0.0, 0.0, 0.0);
      m3_new(&lattrans,
             1, 0, 0,
             0, 1, 0,
             0, 0, 1);
    }

  /* Get unique symmetry operations in cartesians. */
  getSymOpFromHall(hallsymbol, symops, &nsymops);
  /* for (i = 0; i < nsymops; i++) */
  /*   symop_dump(symops + i); */

  maskSymOp(symops, mask, nsymops, transvecs, ntrans, eps);
  fromSymOpToCartesians(symops, mask, nsymops, &convCell, &lattrans);

  /* for (i = 0; i < nsymops; i++) */
  /*   if (mask[i]) */
  /*     symop_dump(symops + i); */

  if (!checkSymOp(symops, mask, nsymops, crys, &convCell, eps))
    return NULL;

  at = generateCell(symops, mask, nsymops, transvecs, ntrans, xred, conc, nat, eps);

  return at;
}
