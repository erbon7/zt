/*      
  *     zt - Simple and partial Mantel Test - version 1.1
  *     copyright (c) Eric Bonnet 2001 - 2007
  *
  *     This program is free software; you can redistribute it and/or modify
  *     it under the terms of the GNU General Public License as published by
  *     the Free Software Foundation; either version 2 of the License, or
  *     (at your option) any later version.
  *
  *     This program is distributed in the hope that it will be useful,
  *     but WITHOUT ANY WARRANTY; without even the implied warranty of
  *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *     GNU General Public License for more details.
  *
  *     You should have received a copy of the GNU General Public License
  *     along with this program; if not, write to the Free Software
  *     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
  *       
  *
  *     Here are all computation routines.
  *     The permutation procedure was written by Glenn C Rhoads,
  *     and is used here with the kind permission from the author.
  *     http://remus.rutgers.edu/~rhoads/Code/perm_lex.c
  *     Thanks, Glenn !
*/

#include <math.h>
#include <stdlib.h>
#include "rr.h"

/* 
  *       fact: compute factorial 
  *       input integer n
  *       return n! as long
*/
long fact (int n) {
  int i;
  long ret = 1;

  for (i = 1; i <= n; i++)
    ret *= i;
  return ret;
}

/* 
  *       somx:   compute simple sum of elements in a half-matrix
  *       input:  matrix pointer, size of the half-matrix
  *       return: sum as double
*/
double somx (double **a, long stop) {
  long i, j;
  double ret = 0;

  for (i = 0; i < stop; i++)
    for (j = 0; j <= i; j++)
      ret += a[i][j];
  return ret;
}

/*      
  *       somx2:  compute square sum of elements in a half-matrix
  *       input:  matrix pointer, size of the half-matrix
  *       return: square sum as double
*/
double somx2 (double **a, long stop) {
  long i, j;
  double ret = 0;

  for (i = 0; i < stop; i++)
    for (j = 0; j <= i; j++)
      ret += a[i][j] * a[i][j];
  return ret;
}

/*      
  *       moy:    compute mean for a half-matrix
  *       input:  matrix pointer, size of the half-matrix
  *       return: mean as double
*/
double moy (double **a, long stop) {
  long i, j, N;
  double ret = 0;
  N = ((stop * (stop - 1)) / 2) + stop;
  for (i = 0; i < stop; i++)
    for (j = 0; j <= i; j++)
      ret += a[i][j];
  ret = ret / N;
  return ret;
}

/* 
  *       ect:    compute standard deviation for a half-matrix
  *       input:  matrix pointer, size of the half-matrix
  *       return: standard deviation as double
*/
double ect (double **a, long stop) {
  long N;
  double ret = 0;
  double lsomx = somx (a, stop);
  double lsomx2 = somx2 (a, stop);
  N = ((stop * (stop - 1)) / 2) + stop;
  ret = sqrt ((lsomx2 - (lsomx * lsomx) / N) / (N - 1));
  return ret;
}

/* 
  *       shake:  shaking elements of a vector at random 
  *       input:  array pointer, size of the array
*/
void shake (long a[], long f) {
  long i;
  long aleat;
  long tmp;
  for (i = 0; i < (f - 1); i++) {
    aleat = i + (1 + rand () % (f - i - 1));
    tmp = a[i];
    a[i] = a[aleat];
    a[aleat] = tmp;
  }
}

/* 
  *       sompx: compute square sum of mean deviations for a half matrix
  *       input:  matrix pointer, size of the half-matrix
  *       return: sum as double
*/
double sompx (double **a, long stop) {
  double ret = 0;
  long N;

  double lsomx = somx (a, stop);
  double lsomx2 = somx2 (a, stop);
  N = ((stop * (stop - 1)) / 2) + stop;

  ret = (lsomx2 - ((lsomx * lsomx) / N));

  return ret;
}


/* 
  *       sompxy: sum x-mean(x) * y-mean(y) for two half matrices
  *       input:  matrix A pointer, matrix B pointer, size of the matrices, mean matrix A, mean matrix B
  *       return: sum as double
*/
double sompxy (double **a, double **b, long stop, double lmoyA, double lmoyB) {
  long i, j;
  double ret = 0;
  for (i = 0; i < stop; i++)
    for (j = 0; j <= i; j++)
      ret += (a[i][j] * b[i][j]) - (lmoyA * lmoyB);

  return ret;
}

/* 
  *       resid:  compute residuals of half-matrix A against half-matrix B
  *       input:  matrix A pointer, matrix B pointer, size of the matrices, mean matrix A, mean matrix B
*/
void resid (double **a, double **b, long stop, double lmoyA, double lmoyB) {
  double coef_b;
  double coef_a;
  long i, j;

  coef_b = sompxy (a, b, stop, lmoyA, lmoyB) / sompx (b, stop);
  coef_a = lmoyA - (coef_b * lmoyB);

  for (i = 0; i < stop; i++)
    for (j = 0; j <= i; j++)
      a[i][j] = a[i][j] - (coef_a + (coef_b * b[i][j]));
}

/* 
  *       norm: normalization of a half-matrix
  *       input:  matrix pointer, size of the matrix      
*/
void norm (double **a, long stop) {
  double lmoya;
  double lecta;
  long i, j;

  lmoya = moy (a, stop);
  lecta = ect (a, stop);

  for (i = 0; i < stop; i++)
    for (j = 0; j <= i; j++)
      a[i][j] = (a[i][j] - lmoya) / lecta;
}

/*
  *       pmt:    partial Mantel test
  *       input:  matrix A pointer, matrix B pointer, matrix C pointer, struct of parameters (see rr.h)
  *       return: 1 if ok
*/
int pmt (double **A, double **B, double **C, struct param *p) {
  double moyA;
  double moyC;
  long i, j, N;
  double r_abc = 0;
  double r_ab = 0;
  double r_ac = 0;
  double r_bc = 0;
  int ret = 0;

  N = p->matsize - 1;

  /* computing mean and standard deviation */
  moyA = moy (A, N);
  moyC = moy (C, N);

  /* computing residuals */
  if (p->raw == 0)
    resid (A, C, N, moyA, moyC);

  /* normalization */
  norm (A, N);
  norm (B, N);
  norm (C, N);

  /* computing initial r_ab */
  for (i = 0; i < N; i++)
    for (j = 0; j <= i; j++)
      r_ab += A[i][j] * B[i][j];
  r_ab = r_ab / (p->numelt - 1);

  /* computing initial r_ac */
  for (i = 0; i < N; i++)
    for (j = 0; j <= i; j++)
      r_ac += A[i][j] * C[i][j];
  r_ac = r_ac / (p->numelt - 1);

  /* computing initial r_bc */
  for (i = 0; i < N; i++)
    for (j = 0; j <= i; j++)
      r_bc += B[i][j] * C[i][j];
  r_bc = r_bc / (p->numelt - 1);

  /* computing initial r_abc */
  r_abc =
    (r_ab -
     (r_ac * r_bc)) / (sqrt (1 - (r_ac * r_ac)) * sqrt (1 - (r_bc * r_bc)));

  p->coef = r_abc;

  /* force exact permutation ? */
  if (p->exact == 0)
    ret = pmt_perm (A, B, C, &r_bc, p);
  else
    ret = pmt_perm_exact (A, B, C, &r_bc, p);

  return (ret);
}

/*
  *       pmt_perm:       randomization procedure for partial Mantel test
  *       input:          matrix A pointer, matrix B pointer, matrix C pointer, r_bc pointer, struct for parameters
  *       return:         1 if ok
*/
int pmt_perm (double **A, double **B, double **C, double *r_bc, struct param *p) {
  long i, j, r, cptega, cptinf, cptsup;
  double r_ab, r_ac, rrand, epsilon;
  long *ord;

  cptega = 1;
  cptsup = 0;
  cptinf = 0;
  epsilon = 0.0000001;

  ord = malloc (p->matsize * sizeof (long));
  if (ord == NULL)
    return (-1);

  for (i = 0; i < p->matsize; i++)
    ord[i] = i;

  for (r = 0; r < p->numrand; r++) {
    rrand = 0;
    shake (ord, p->matsize);

    r_ab = 0;
    for (i = 1; i < p->matsize; i++) {
      for (j = 0; j < i; j++) {
        if (ord[j] < ord[i])
          r_ab += B[i - 1][j] * A[ord[i] - 1][ord[j]];
        else
          r_ab += B[i - 1][j] * A[ord[j] - 1][ord[i]];
      }
    }
    r_ab = r_ab / (p->numelt - 1);

    r_ac = 0;
    for (i = 1; i < p->matsize; i++) {
      for (j = 0; j < i; j++) {
        if (ord[j] < ord[i])
          r_ac += C[i - 1][j] * A[ord[i] - 1][ord[j]];
        else
          r_ac += C[i - 1][j] * A[ord[j] - 1][ord[i]];
      }
    }
    r_ac = r_ac / (p->numelt - 1);

    rrand = (r_ab - (r_ac * (*r_bc))) / (sqrt (1 - (r_ac * r_ac)) * sqrt (1 - ((*r_bc) * (*r_bc))));

    if (fabs (rrand - p->coef) <= epsilon * fabs (rrand))
      cptega += 1;
    else {
      if (rrand > p->coef)
        cptsup += 1;
      if (rrand < p->coef)
        cptinf += 1;
    }
  }

  free (ord);

  if (p->coef < 0)
    p->proba = (double) (cptinf + cptega) / (p->numrand + 1);
  else
    p->proba = (double) (cptsup + cptega) / (p->numrand + 1);

  return (1);
}

/*
  *       smt:    simple Mantel test
  *       input:  matrix A pointer, matrix B pointer, struct for results  
  *       return  1 if ok
*/
int smt (double **A, double **B, struct param *p) {
  long i, j;
  double zini;
  long N = p->matsize - 1;
  int ret = 0;

  norm (A, N);
  norm (B, N);

  /*computing initial z */
  zini = 0;
  for (i = 0; i < N; i++)
    for (j = 0; j <= i; j++)
      zini += A[i][j] * B[i][j];

  p->coef = zini / (p->numelt - 1);

  if (p->exact == 0)
    ret = smt_perm (A, B, p);
  else
    ret = smt_perm_exact (A, B, p);

  return (ret);
}

/*
  *       smt_perm:       randomization procedure for simple Mantel test
  *       input:          matrix A pointer, matrix B pointer, struct for results
  *       return          1 if ok
*/
int smt_perm (double **A, double **B, struct param *p) {
  long i, j, r, cptega, cptinf, cptsup;
  double zrand, epsilon;
  long *ord;

  epsilon = 0.0000001;

  ord = malloc (p->matsize * sizeof (long));
  if (ord == NULL)
    return (-1);

  for (i = 0; i < p->matsize; i++)
    ord[i] = i;

  cptega = 1;
  cptsup = 0;
  cptinf = 0;
  zrand = 0;

  for (r = 0; r < p->numrand; r++) {
    zrand = 0;
    shake (ord, p->matsize);

    for (i = 1; i < p->matsize; i++) {
      for (j = 0; j < i; j++) {
        if (ord[j] < ord[i])
          zrand += A[i - 1][j] * B[ord[i] - 1][ord[j]];
        else
          zrand += A[i - 1][j] * B[ord[j] - 1][ord[i]];
      }
    }
    zrand = zrand / (p->numelt - 1);

    if (fabs (zrand - p->coef) <= epsilon * fabs (zrand))
      cptega += 1;
    else {
      if (zrand > (p->coef))
        cptsup += 1;
      if (zrand < (p->coef))
        cptinf += 1;
    }
  }

  if (p->coef < 0)
    p->proba = (double) (cptinf + cptega) / (p->numrand + 1);
  else
    p->proba = (double) (cptsup + cptega) / (p->numrand + 1);

  free (ord);

  return (1);
}

/*
  *       smt_perm_exact:         exact permutation procedure for simple Mantel test
  *       input:                  matrix pointer A, matrix pointer B, struct for results
  *       return                  1 if ok
*/
int smt_perm_exact (double **A, double **B, struct param *p) {
  long i, j, r, s, li, lj, temp, n, cpt;
  long cptega, cptinf, cptsup;
  long *ord;
  double zrand, epsilon;

  epsilon = 0.0000001;

  ord = malloc (p->matsize * sizeof (long));
  if (ord == NULL)
    return (-1);

  for (i = 0; i < p->matsize; i++)
    ord[i] = i;

  i = 1;
  n = p->matsize - 1;
  cptega = 0;
  cptsup = 0;
  cptinf = 0;
  cpt = 0;

  while (i >= 0) {
    cpt++;
    zrand = 0;

    for (li = 1; li < p->matsize; li++) {
      for (lj = 0; lj < li; lj++) {
        if (ord[lj] < ord[li])
          zrand += A[li - 1][lj] * B[ord[li] - 1][ord[lj]];
        else
          zrand += A[li - 1][lj] * B[ord[lj] - 1][ord[li]];
      }
    }
    zrand = zrand / (p->numelt - 1);

    if (fabs (zrand - p->coef) <= epsilon * fabs (zrand))
      cptega += 1;
    else {
      if (zrand > p->coef)
        cptsup += 1;
      if (zrand < p->coef)
        cptinf += 1;
    }

    i = n - 1;
    while (ord[i] > ord[i + 1])
      i--;

    j = n;
    while (ord[i] > ord[j])
      j--;

    if (i >= 0) {
      temp = ord[i];
      ord[i] = ord[j];
      ord[j] = temp;

      r = n;
      s = i + 1;

      while (r > s) {
        temp = ord[r];
        ord[r] = ord[s];
        ord[s] = temp;
        r--;
        s++;
      }
    }
  }

  if (p->coef < 0)
    p->proba = (double) (cptinf + cptega) / cpt;
  else
    p->proba = (double) (cptsup + cptega) / cpt;

  free (ord);

  return (1);
}

/*
  *       pmt_perm_exact: exact permutation procedure for partial Mantel test
  *       input:          matrix A pointer, matrix B pointer, matrix C pointer, r_bc pointer, struct for results
  *       return          1 if ok 
*/
int pmt_perm_exact (double **A, double **B, double **C, double *r_bc, struct param *p) {
  long i, j, r, s, li, lj, temp, n, cpt;
  long cptega, cptinf, cptsup;
  long *ord;
  double r_ab, r_ac, rrand, epsilon;

  epsilon = 0.0000001;

  ord = malloc (p->matsize * sizeof (long));
  if (ord == NULL)
    return (-1);

  for (i = 0; i < p->matsize; i++)
    ord[i] = i;

  i = 1;
  n = p->matsize - 1;
  cptega = 0;
  cptsup = 0;
  cptinf = 0;
  cpt = 0;

  while (i >= 0) {
    cpt++;
    rrand = 0;
    r_ab = 0;
    r_ac = 0;

    for (li = 1; li < p->matsize; li++) {
      for (lj = 0; lj < li; lj++) {
        if (ord[lj] < ord[li])
          r_ab += B[li - 1][lj] * A[ord[li] - 1][ord[lj]];
        else
          r_ab += B[li - 1][lj] * A[ord[lj] - 1][ord[li]];
      }
    }
    r_ab = r_ab / (p->numelt - 1);

    for (li = 1; li < p->matsize; li++) {
      for (lj = 0; lj < li; lj++) {
        if (ord[lj] < ord[li])
          r_ac += C[li - 1][lj] * A[ord[li] - 1][ord[lj]];
        else
          r_ac += C[li - 1][lj] * A[ord[lj] - 1][ord[li]];
      }
    }
    r_ac = r_ac / (p->numelt - 1);

    rrand = (r_ab - (r_ac * (*r_bc))) / (sqrt (1 - (r_ac * r_ac)) * sqrt (1 - ((*r_bc) * (*r_bc))));

    if (fabs (rrand - p->coef) <= epsilon * fabs (rrand))
      cptega += 1;
    else {
      if (rrand > p->coef)
        cptsup += 1;
      if (rrand < p->coef)
        cptinf += 1;
    }

    i = n - 1;
    while (ord[i] > ord[i + 1])
      i--;

    j = n;
    while (ord[i] > ord[j])
      j--;

    if (i >= 0) {
      temp = ord[i];
      ord[i] = ord[j];
      ord[j] = temp;

      r = n;
      s = i + 1;

      while (r > s) {
        temp = ord[r];
        ord[r] = ord[s];
        ord[s] = temp;
        r--;
        s++;
      }
    }
  }

  if (p->coef < 0)
    p->proba = (double) (cptinf + cptega) / cpt;
  else
    p->proba = (double) (cptsup + cptega) / cpt;

  free (ord);

  return (1);
}
