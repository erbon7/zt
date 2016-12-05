/*      
   *    zt - Simple and partial Mantel Test - version 1.1
   *    copyright (c) Eric Bonnet 2001 - 2007
   *
   *    This program is free software; you can redistribute it and/or modify
   *    it under the terms of the GNU General Public License as published by
   *    the Free Software Foundation; either version 2 of the License, or
   *    (at your option) any later version.
   *
   *    This program is distributed in the hope that it will be useful,
   *    but WITHOUT ANY WARRANTY; without even the implied warranty of
   *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   *    GNU General Public License for more details.
   *
   *    You should have received a copy of the GNU General Public License
   *    along with this program; if not, write to the Free Software
   *    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
   *       
   *
   *    Function prototypes.
 */

#define MIN_MAT_SIZE 5          /* mimimum size for matrices */
#define MAX_EXACT_SIZE 12       /* maximum size for exact permutation procedure */
#define EXACT_PROC_SIZE 8       /* size critical for automatic permutation procedure */

struct param {
    double coef;                /* reference statistic */
    double proba;               /* p-value */
    long numrand;               /* number of randomizations */
    long matsize;               /* size of matrices */
    long numelt;                /* number of elements in the half-matrix without diagonal values */
    int partial;                /* option partial 0|1 */
    int raw;                    /* option raw 0|1 */
    int help;                   /* option help 0|1 */
    int exact;                  /* option exact permutation 0|1 */
    int licence;                /* option licence terms 0|1 */
};

/* prototypes of functions */
long fact (int n);
double somx (double **a, long stop);
double somx2 (double **a, long stop);
double moy (double **a, long stop);
double ect (double **a, long stop);
void shake (long a[], long f);
double sompx (double **a, long stop);
double sompxy (double **a, double **b, long stop, double lmoyA, double lmoyB);
void resid (double **a, double **b, long stop, double lmoyA, double lmoyB);
void norm (double **a, long stop);
int pmt (double **A, double **B, double **C, struct param *p);
int pmt_perm (double **A, double **B, double **C, double *r_bc, struct param *p);
int pmt_perm_exact (double **A, double **B, double **C, double *r_bc, struct param *p);
int smt (double **A, double **B, struct param *p);
int smt_perm (double **A, double **B, struct param *p);
int smt_perm_exact (double **A, double **B, struct param *p);
