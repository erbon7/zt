# zt

zt is a software for testing correlation in distance matrices, using a permutation test. Simple correlation can be analysed for two distance matrices, and with three matrices partial correlation can be tested. This procedure is known as the Mantel test. For instance, genetic and geographical distances can be tested in group of species. 

For more informations, see [Bonnet and Van de Peer](https://www.jstatsoft.org/article/view/v007i10).

zt is a command line tool. The source files can be compiled easily, but pre-compiled binaries are available. See the readme.txt file for more instructions on how to compile, test and use the software.

Please cite the above reference if you use the software.

# Binaries
Pre-compiled binaries are available for Linux, Solaris, Windows and Mac OS X [here](https://www.jstatsoft.org/article/view/v007i10)

# Compilation

Compilation can be done with a c compiler such as gcc with the command:
```
gcc -o zt zt.c rr.c -lm
```

# Installation and usage 

First create a directory where you want to store the soft. Then copy the tar zt.tar in this directory and un-tar it:
```
tar xvf zt_linux.tar
```
This will expand the main program and the test files.

A simple Mantel test can be performed with the command:
```
  zt -s file1 file2 number_of_randomizations
```

while for a partial Mantel test, the command will be:

```
zt -p file1 file2 file3 number_of_randomizations

```
Options for the command line:
  -l display licence.
  -r partial Mantel test with raw option.
  -e force exact permutations procedure.
  -h display this help.

The number of randomizations are allowed between 99 and 999999999.

Examples with included test files:

Simple Mantel test:
```
./zt -s assoc.txt gond.txt 10000

zt - version 1.1

copyright (c) Eric Bonnet 2001 - 2007

File A:			assoc.txt
File B:			gond.txt
Size of matrices:	8 x 8
Number of iterations:	10000
Options:		simple 

Randomizing...

r =			-0.605379
p =			0.001900 (one-tailed)
```

Partial mantel test:
```
./zt -p gene.txt env.txt geo.txt 10000

zt - version 1.1

copyright (c) Eric Bonnet 2001 - 2007

File A:			gene.txt
File B:			env.txt
File C:			geo.txt
Size of matrices:	21 x 21
Number of iterations:	10000
Options:		partial residuals 

Randomizing...

r =			0.313143
p =			0.000700 (one-tailed)

```


For further details on the results and their interpretation, see the document zt.pdf.


# Data files format

Input files must be flat ASCII files. Matrix are half-matrix (lower part of the square matrix), without (!) diagonal terms. The first number is the size of the complete matrix. For Example if you have the following half-matrix data:
```
0.2 
0.3 0.4 
```
The size of this matrix is 3 x 3. So, the correct data file for zt will be:
```
3
0.2 
0.3 0.4 
```
Be careful, floating point must be dot(.) not comma(,). Try to avoid trailing spaces as well.
