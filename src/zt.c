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
   *    Main body of the program:
   *    Verifications, data storage  and the launch of the test
 */

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include "rr.h"

int main (int argc, char *argv[]) {

  double **matA;                /* pointer to matrix A */
  double **matB;                /* pointer to matrix B */
  double **matC;                /* pointer to matrix C */
  FILE *fileA;                  /* pointer to file A */
  FILE *fileB;                  /* pointer to file B */
  FILE *fileC;                  /* pointer to file C */
  char fnameA[255];             /* array for filename A */
  char fnameB[255];             /* array for filename B */
  char fnameC[255];             /* array for filename C */
  long NA;                      /* matrix A size */
  long NB;                      /* matrix B size */
  long NC;                      /* matrix C size */
  long i, j;
  int c;
  int tmem, res;
  struct param myp;             /* struct for storing parameters and results see rr.h */

  tmem = 0;
  res = 0;

  myp.coef = 0;
  myp.proba = 0;
  myp.coef = 0;
  myp.numrand = 0;
  myp.matsize = 0;
  myp.numelt = 0;
  myp.partial = -1;
  myp.raw = 0;
  myp.help = 0;
  myp.exact = 0;
  myp.licence = 0;
  c = 0;

  /* seed for random function */
  srand ((unsigned) time (0));

  /* header */
  printf ("\nzt - version 1.1\n\n");
  printf ("copyright (c) Eric Bonnet 2001 - 2007\n\n");

  /* parse arguments a la Kernighan & Ritchie */
  if (argc == 1 || *argv[1] != '-') {
    myp.help = 1;
    printf ("Error: unknown option\n");
  }
  else {
    while (--argc > 0 && (*++argv)[0] == '-')
      while ((c = *++argv[0]))
        switch (c) {
        case 'p':
          myp.partial = 1;
          break;
        case 'r':
          myp.raw = 1;
          break;
        case 'h':
          myp.help = 1;
          break;
        case 'e':
          myp.exact = 1;
          break;
        case 's':
          myp.partial = 0;
          break;
        case 'l':
          myp.licence = 1;
          break;
        default:
          printf ("Error: unknown option %c\n", c);
          myp.help = 1;
          argc = 0;
          break;
        }
  }

  /* print help usage */
  if (myp.help == 1) {
    printf ("\nUsage:\n\n");
    printf
      ("Simple Mantel test:\nzt -s file1 file2 number_of_randomizations\n\n");
    printf
      ("Partial Mantel test:\nzt -p file1 file2 file3 number_of_randomizations\n");
    printf ("\nOptions:\n");
    printf ("-r partial Mantel test with raw option\n");
    printf ("-e force exact permutations procedure\n");
    printf ("-l print licence terms\n");
    printf ("-h display this help\n\n");
    return(1);
  }

/* print licence terms */
  if (myp.licence == 1) {
    printf
      ("\nThis program is free software; you can redistribute it and/or modify\n");
    printf
      ("it under the terms of the GNU General Public License as published by\n");
    printf
      ("the Free Software Foundation; either version 2 of the License, or\n");
    printf ("(at your option) any later version.\n\n");
    printf
      ("This program is distributed in the hope that it will be useful,\n");
    printf
      ("but WITHOUT ANY WARRANTY; without even the implied warranty of\n");
    printf
      ("MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n");
    printf ("GNU General Public License for more details.\n\n");
    printf
      ("You should have received a copy of the GNU General Public License\n");
    printf ("along with this program; if not, write to the Free Software\n");
    printf
      ("Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307\n\n");
    return(1);
  }

  /* test mandatory parameters */
  if (myp.partial == -1) {
    printf ("Error: -s or -p mandatory parameters.\n\n");
    return(1);
  }

  /* test number of arguments */
  if (myp.partial == 0) {
    if ((myp.exact == 0 && argc < 3) || (myp.exact == 1 && argc < 2)) {
      printf ("Error: bad number of arguments.\n\n");
      return(1);
    }
  }
  else {
    if ((myp.exact == 0 && argc < 4) || (myp.exact == 1 && argc < 3)) {
      printf ("Error: bad number of arguments.\n\n");
      return(1);
    }
  }

  /* get filenames */
  strcpy (fnameA, *argv);
  strcpy (fnameB, *++argv);
  if (myp.partial == 1)
    strcpy (fnameC, *++argv);

  /* open files and look for matrix size */
  fileA = fopen (fnameA, "r");
  if (fileA == NULL) {
    printf ("Error: can't open file %s.\n\n", fnameA);
    return(1);
  }
  fscanf (fileA, "%li", &NA);

  fileB = fopen (fnameB, "r");
  if (fileB == NULL) {
    printf ("Error: can't open file %s.\n\n", fnameB);
    return(1);
  }
  fscanf (fileB, "%li", &NB);

  if (myp.partial == 1) {
    fileC = fopen (fnameC, "r");
    if (fileC == NULL) {
      printf ("Error: can't open file %s.\n\n", fnameC);
      return(1);
    }
    fscanf (fileC, "%li", &NC);
  }

  /* test matrices size */
  if (myp.partial == 0) {
    if (NA != NB) {
      printf ("Error: matrices of different size.\n\n");
      return(1);
    }
  }
  else {
    if (NA != NB || NA != NC || NB != NC) {
      printf ("Error: matrices of different size.\n\n");
      return(1);
    }
  }

  /* define and test matrices size */
  /*myp.matsize = NA + 1; */
  myp.matsize = NA;
  myp.numelt = (myp.matsize * (myp.matsize - 1)) / 2;

  if (myp.matsize < MIN_MAT_SIZE) {
    printf ("Error: matrix size must be >= %i.\n\n", MIN_MAT_SIZE);
    return(1);
  }

  if (myp.exact == 1 && myp.matsize > MAX_EXACT_SIZE) {
    printf
      ("Error: matrix size is too big for exact permutations (should be <= %i).\n\n",
       MAX_EXACT_SIZE);
    return(1);
  }

  /* force exact permutation procedure if size too small */
  if (myp.matsize < EXACT_PROC_SIZE)
    myp.exact = 1;

  /* define and test number of randomizations */
  if (myp.exact == 0) {
    myp.numrand = atoi (*++argv);

    if (myp.numrand < 99 || myp.numrand > 999999999) {
      printf
        ("Error: Number of iterations must be between 99 and 999999999.\n\n");
      return(1);
    }
  }

  /* memory allocation */
  matA = malloc (NA * sizeof (double));
  matB = malloc (NA * sizeof (double));

  if (myp.partial == 1)
    matC = malloc (NA * sizeof (double));

  if (myp.partial == 0) {
    for (i = 0; i < NA; i++) {
      matA[i] = malloc ((i + 1) * sizeof (double));
      matB[i] = malloc ((i + 1) * sizeof (double));
      if (matA[i] == NULL || matB[i] == NULL) {
        tmem = 1;
        break;
      }
    }
  }
  else {
    for (i = 0; i < NA; i++) {
      matA[i] = malloc ((i + 1) * sizeof (double));
      matB[i] = malloc ((i + 1) * sizeof (double));
      matC[i] = malloc ((i + 1) * sizeof (double));
      if (matA[i] == NULL || matB[i] == NULL || matC[i] == NULL) {
        tmem = 1;
        break;
      }
    }
  }

  if (tmem == 1) {
    printf ("Error: not enough memory.\n\n");
    return(1);
  }

  /* get data into matrices */
  if (myp.partial == 0) {
    for (i = 0; i < NA-1; i++)
      for (j = 0; j <= i; j++) {
        fscanf (fileA, "%lf", &matA[i][j]);
        fscanf (fileB, "%lf", &matB[i][j]);
      }
  }
  else {
    for (i = 0; i < NA-1; i++)
      for (j = 0; j <= i; j++) {
        fscanf (fileA, "%lf", &matA[i][j]);
        fscanf (fileB, "%lf", &matB[i][j]);
        fscanf (fileC, "%lf", &matC[i][j]);
      }
  }

  if (myp.exact == 1)
    myp.numrand = fact (myp.matsize);

  /* print parameters of the test */
  printf ("File A:\t\t\t%s\n", fnameA);
  printf ("File B:\t\t\t%s\n", fnameB);
  if (myp.partial == 1)
    printf ("File C:\t\t\t%s\n", fnameC);
  printf ("Size of matrices:\t%li x %li\n", myp.matsize, myp.matsize);
  printf ("Number of iterations:\t%li\n", myp.numrand);
  printf ("Options:\t\t");
  if (myp.partial == 0)
    printf ("simple ");
  else {
    printf ("partial ");
    if (myp.raw == 1)
      printf ("raw ");
    else
      printf ("residuals ");
  }
  if (myp.exact == 1)
    printf ("exact ");
  printf ("\n\n");
  printf ("Randomizing...\n\n");

  /* launch the test */
  if (myp.partial == 1) {
    res = pmt (matA, matB, matC, &myp);
    if (res) {
      printf ("r =\t\t\t%f\n", myp.coef);
      printf ("p =\t\t\t%f (one-tailed)\n\n", myp.proba);
    }
    else
      printf
        ("An error has occurred during permutation procedure.\nPlease retry.\n\n");
  }
  else {
    res = smt (matA, matB, &myp);
    if (res) {
      printf ("r =\t\t\t%f\n", myp.coef);
      printf ("p =\t\t\t%f (one-tailed)\n\n", myp.proba);
    }
    else
      printf
        ("An error has occurred during permutation procedure.\nPlease retry.\n\n");
  }
  return(0);
}
