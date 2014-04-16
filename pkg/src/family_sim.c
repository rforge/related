/*-------------------------------family_sim.c-------------------------------*/
/* This program takes allele frequencies from a given data set and          */
/* generates a user-defined number of pairs each of parent-offspring,       */
/* full-sib, half-sib, and unrelated individuals. These genotypes can then  */
/* be used for power analyses of various other applications.                */
/*                                                                          */
/* Copyright (C) 2009 Tim Frasier                                           */
/*                                                                          */
/* This program is free software; you can redistribute it and/or modify it  */
/* under the terms of the GNU General Public License as published by the    */
/* Free Software Foundation; either version 2 of the License, or (at your   */
/* option) any later version.                                               */
/*                                                                          */
/* This program is distributed in the hope that it will be useful, but      */
/* WITHOUT ANY WARRANTY; without even the implied warranty of               */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Public */
/* License for more details. (http://www.gnu.org/licenses/gpl.html)         */
/*--------------------------------------------------------------------------*/

/* Modifications (2013-08-12):
 *  replace the RNG (swapped GSL for R)
 *  replace Mendel arrays with a coin_flip_heads() function
 *  remove use of temporary file "freqs"
 *  replace variable length arrays (VLAs) with explicitly allocated arrays
 *    (the functions alloc_2d and free_2d help accomplish this)
 *    VLAs are an optional part of the C11 standard and should be avoided.
 *  some variables were renamed (count[1234] -> [ijk])
 *  freqs3 was removed and freqs2 used for cumulative frequencies
 *  code refactored into a family_sim function which takes integer arguments
 *    these arguements replace prompts for user imput
 *  the prompt for an input file was temporarily replaced with a hardcoded
 *    filename of "frequencies"
 *  the seed is now fixed as 42 rather than using the system clock
 *    this makes testing incremental changes possible
 *  indexing of individuals2 was changed to avoid writing beyond the
 *    allocated storage (was the change correct?)
 *  moved variable declarations to top of code to conform to ISO C90
 *    (this helps portability to systems with antiquated compilers)
 */
/* Further modifications (2013-08-15):
 *  created a header file to declare family_sim
 *  modified family_sim to take input via arrays rather than a file
 *  frequencies is now a 2D ragged array packed into a 1D array to
 *    facilitate passing this parameter in from R. there is one logical
 *    row per locus, which is indexed by f_rows. row i is length nalleles[i].
 *  family_sim now accepts integers to identify alleles, rather than
 *    using their index within the arrays. alleles is a packed 2D ragged
 *    array with the same shape as frequencies.
 *  freqs2 was replaced with a cfreqs, which has the same structure
 *    as the new frequencies parameter.
 *  nalleles was previously the number of alleles at the most polymorphic
 *    locus. it is now an array of #alleles at each locus.
 *    this array is used to index the ragged packed 2D arrays
 *    (which are alleles, frequencies and cfreqs).
 */
/* More modifications (2013-08-20):
 *  switched to using methods and macros from Rinternals.h for argument
 *    passing and memory allocation.
 *    note: R_alloc takes care of freeing memory, even when family_sim
 *    is interrupted by the user.
 *  individuals2 is no longer written to a file; it is instead copied to
 *    an R integer matrix and returned by family_sim through R's
 *    .External interface. The copying step involves a matrix transpose.
 *  the RNG is no longer seeded within this code; it shares state with R.
 *    (ie, calls to unif_rand() occur between GetRNGstate() and PutRNGstate())
 */
/* Compiling:
 *   gcc -Wall -Wextra -std=c90 -pedantic \
 *    -I/usr/include/R -fPIC -g family_sim.c -o fsim.so
 */

/*----------------*/
/* INCLUDED FILES */
/*----------------*/
#include <stdio.h>
#include <stdlib.h>

/* related.h includes R.h and Rinternals.h inside an ifdef _cpluscplus */
#include "related.h"

/* The inline keyword is not available in all compilers. */
#if __STDC_VERSION__ >= 199901L
#else
#    define inline  /* empty */
#endif

/*--------------------------------------------*/
/* This function returns a boolean value      */
/* which replaces the usage of Mendel arrays. */
/*--------------------------------------------*/
static inline int coin_flip_heads(void) {
    return (unif_rand() < 0.5);
}

/*---------------------------------------------*/
/* This function uses R's .External interface. */
/*---------------------------------------------*/
SEXP family_sim(SEXP args) {

    /* function parameters */
    SEXP arg_ninds;
    SEXP arg_nloci;
    SEXP arg_nalleles;
    SEXP arg_alleles;
    SEXP arg_frequencies;

    /* typed versions parameters */
    int ninds;             /* The number of simulated pairs to generate. */
    int nloci;             /* The number of loci analyzed. */
    int *nalleles;         /* Array of num alleles analysed at each locus. */
    int *alleles;          /* Packed 2D array of alleles,
                              length = Sum nalleles[i]. */
    double *frequencies;   /* Packed 2D array of frequencies of alleles. */

    /* returned value */
    SEXP ret_sims;          /* an integer matrix copy of individuals2 */

    /* indexing for returned value */
    int *sims;              /* typed pointer to ret_sims */
    int **sim_cols;         /* columns of sims */

    /* local variables */
    int **a_rows;           /* pointers for unpacking alleles */
    double **f_rows;        /* pointers for unpacking frequencies */
    double *cfreqs;         /* cumulative frequencies */
    double **cf_rows;       /* pointers for unpacking cumulative frequencies */
    int *individuals;       /* packed 2D array of alleles per locus for
                               intermediate new individuals */
    int *individuals2;      /* packed 2D array of alleles per locus for
                               simulated new individuals */
    int **ind_rows;         /* pointers for unpacking individuals */
    int **ind2_rows;        /* pointers for unpacking individuals2 */
    int new_ninds;          /* size of first dimension of individuals */
    int new_ninds2;         /* size of first dimension of individuals2 */
    int *pair;              /* packed 2D array for mating pairs */
    int *pair_rows[2];      /* pointers for indexing mating pairs */
    int *offspring;         /* array intermediate offspring */
    int *off_rows[1];       /* pointers for indexing offspring */
    int i, j, k;            /* loop indices */
    int len, z;             /* temporary counters */

    /* extract args (with Lispy macros) and do type checking */
    args = CDR(args);
    if (!isInteger(arg_ninds = CAR(args)))
        error("ninds should be an integer");
    args = CDR(args);
    if (!isInteger(arg_nloci = CAR(args)))
        error("nloci should be an integer");
    args = CDR(args);
    if (!isInteger(arg_nalleles = CAR(args)))
        error("nalleles should be an integer vector");
    args = CDR(args);
    if (!isInteger(arg_alleles = CAR(args)))
        error("alleles should be an integer vector");
    args = CDR(args);
    if (!isNumeric(arg_frequencies = CAR(args)))
        error("frequencies should be a numeric vector");

    /* do type casting on args for convenience */
    ninds       = *INTEGER(arg_ninds);
    nloci       = *INTEGER(arg_nloci);
    nalleles    =  INTEGER(arg_nalleles);
    alleles     =  INTEGER(arg_alleles);
    frequencies =  REAL(arg_frequencies);

    new_ninds  = ((2*ninds*4)+(3*ninds));
    new_ninds2 =  (2*ninds*4);

    /* some basic assertions on the input */
    if (ninds <= 0) error("require ninds > 0");
    if (nloci <= 0) error("require nloci > 0");

    /* find len = Sum nalleles[i] */
    for (i = len = 0; i < nloci; i++) {
        if (nalleles[i] <= 0) error("require all(nalleles > 0)");
        len += nalleles[i];
    }

    /* set pointers into packed 2D arrays for convenient access */
    /* alleles and frequencies have nloci "rows" packed into them */
    a_rows = (int**) R_alloc(nloci, sizeof(int*));
    f_rows = (double**) R_alloc(nloci, sizeof(int*));

    /* also allocate space for cumulative frequencies, which are */
    /* used in the random allele selection process.              */
    cfreqs = (double*) R_alloc(len, sizeof(double));
    cf_rows = (double**) R_alloc(nloci, sizeof(int*));

    a_rows[0]  = &alleles[0];
    f_rows[0]  = &frequencies[0];
    cf_rows[0] = &cfreqs[0];
    for (i = 1; i < nloci; i++) {
        a_rows[i]  = a_rows[i-1]  + nalleles[i-1];
        f_rows[i]  = f_rows[i-1]  + nalleles[i-1];
        cf_rows[i] = cf_rows[i-1] + nalleles[i-1];
    }

    /*************************************************/
    /* CONVERT TO CUMULATIVE, NORMALIZED FREQUENCIES */
    /*************************************************/
    for (i = 0; i < nloci; i++) {
        /* put cumulative frequencies in cfreqs, indexed by cf_rows */
        cf_rows[i][0] = f_rows[i][0];
        for (j = 1; j < nalleles[i]; j++) {
            cf_rows[i][j] = cf_rows[i][j-1] + f_rows[i][j];
        }

        /* normalize cfreqs using cf_rows[i][nalleles[i]-1] */
        for (j = 0; j < nalleles[i]; j++) {
            cf_rows[i][j] /= cf_rows[i][nalleles[i]-1];
        }
    }

#if 0
    /* some debug printing */
    printf("ninds = %d, nloci = %d\n", ninds, nloci);
    for (i = 0; i < nloci; i++) {
        printf("locus %d:\n", i+1);
        for (j = 0; j < nalleles[i]; j++) {
            printf("\t\t%d : %.2f,", a_rows[i][j], cf_rows[i][j]);
        }
        printf("\n");
    }

    return R_NilValue;
#endif

    /************************************/
    /* MAKE ARRAYS WITH NEW INDIVIDUALS */
    /************************************/
    individuals  = (int*) R_alloc(new_ninds  * (2*nloci), sizeof(int));
    individuals2 = (int*) R_alloc(new_ninds2 * (2*nloci), sizeof(int));
    ind_rows  = (int**) R_alloc(new_ninds, sizeof(int*));
    ind2_rows = (int**) R_alloc(new_ninds, sizeof(int*));

    ind_rows[0]  = &individuals[0];
    ind2_rows[0] = &individuals2[0];
    for (i = 1; i < new_ninds;  i++) ind_rows[i]  = ind_rows[i-1]  + (2*nloci);
    for (i = 1; i < new_ninds2; i++) ind2_rows[i] = ind2_rows[i-1] + (2*nloci);

    /* share RNG state with R */
    GetRNGstate();

    for (i = 0; i < new_ninds; i++) {
        for (j = 0; j < nloci; j++) {
            for (z = 0; z < 2; z++) {

                /* use a random number and cumulative */
                /* frequencies to pick an allele here */
                double rand = unif_rand();

                for (k = 0; k < nalleles[j]; k++) {
                    if (rand < cf_rows[j][k]) {
                        ind_rows[i][2*j+z] = a_rows[j][k];
                        break;
                    }
                }
            }
        }
    }

    /*****************************/
    /* GENERATE PARENT-OFFSPRING */
    /*****************************/
    pair = (int*) R_alloc(2 * (2*nloci), sizeof(int));
    pair_rows[0] = &pair[0];
    pair_rows[1] = &pair[2*nloci];

    offspring = (int*) R_alloc(1 * (2*nloci), sizeof(int));
    off_rows[0] = &offspring[0];

    for (k = 0; k < ninds; k++) {
        /*-----------------------------*/
        /* GET APPROPRIATE MATING PAIR */
        /*-----------------------------*/
        for (j = 0; j < 2*nloci; j++) {
            pair_rows[0][j] = ind_rows[(k*2)][j];
            pair_rows[1][j] = ind_rows[(k*2)+1][j];
        }

        /*------------------------*/
        /* CREATE A NEW OFFSPRING */
        /*------------------------*/
        for (i = 0; i < 2*nloci; i += 2) {
            if (coin_flip_heads())
                off_rows[0][i] = pair_rows[0][i];
            else
                off_rows[0][i] = pair_rows[0][i+1];
        }

        for (i = 1; i < 2*nloci; i += 2) {
            if (coin_flip_heads())
                off_rows[0][i] = pair_rows[1][i-1];
            else
                off_rows[0][i] = pair_rows[1][i];
        }

        /*--------------------------*/
        /* PLACE PAIR BACK IN ARRAY */
        /*--------------------------*/
        for (i = 0; i < 2*nloci; i++) {
            ind2_rows[k*2][i]     = pair_rows[0][i];
            ind2_rows[(k*2)+1][i] = off_rows[0][i];
        }
    }

    /**********************/
    /* GENERATE FULL-SIBS */
    /**********************/
    for (k = ninds; k < (2*ninds); k++) {
        /*-----------------------------*/
        /* GET APPROPRIATE MATING PAIR */
        /*-----------------------------*/
        for (j = 0; j < 2*nloci; j++) {
            pair_rows[0][j] = ind_rows[k*2][j];
            pair_rows[1][j] = ind_rows[(k*2)+1][j];
        }

        /*------------------------*/
        /* GENERATE TWO OFFSPRING */
        /*------------------------*/
        for (j = 0; j < 2; j++) {
            /*------------------------*/
            /* CREATE A NEW OFFSPRING */
            /*------------------------*/
            for (i = 0; i < 2*nloci; i += 2) {
                if (coin_flip_heads())
                    off_rows[0][i] = pair_rows[0][i];
                else
                    off_rows[0][i] = pair_rows[0][i+1];
            }

            for (i = 1; i < 2*nloci; i += 2) {
                if (coin_flip_heads())
                    off_rows[0][i] = pair_rows[1][i-1];
                else
                    off_rows[0][i] = pair_rows[1][i];
            }

            /*-------------------------------*/
            /* PLACE OFFSPRING BACK IN ARRAY */
            /*-------------------------------*/
            for (i = 0; i < 2*nloci; i++) {
                ind2_rows[(k*2)+j][i] = off_rows[0][i];
            }
        }

    }

    /**********************/
    /* GENERATE HALF-SIBS */
    /**********************/
    z = 0;

    for (k = (2*ninds); k < (3*ninds); k++) {
        /*------------------*/
        /* GET FIRST PARENT */
        /*------------------*/
        for (j = 0; j < 2*nloci; j++) {
            pair_rows[0][j] = ind_rows[(k*2)+z][j];
        }

        /*-----------------------------------------------------*/
        /* GENERATE TWO OFFSPRING EACH WITH A DIFFERENT SPOUSE */
        /*-----------------------------------------------------*/
        for (j = 0; j < 2; j++) {
            /*------------*/
            /* GET SPOUSE */
            /*------------*/
            for (i = 0; i < 2*nloci; i++) {
                pair_rows[1][i] = ind_rows[(k*2)+z+(j+1)][i];
            }

            /*-------------------------*/
            /* CREATE A NEW OFFSPRING  */
            /*-------------------------*/
            for (i = 0; i < 2*nloci; i += 2) {
                if (coin_flip_heads())
                    off_rows[0][i] = pair_rows[0][i];
                else
                    off_rows[0][i] = pair_rows[0][i+1];
            }

            for (i = 1; i < 2*nloci; i += 2) {
                if (coin_flip_heads())
                    off_rows[0][i] = pair_rows[1][i-1];
                else
                    off_rows[0][i] = pair_rows[1][i];
            }

            /*-------------------------------*/
            /* PLACE OFFSPRING BACK IN ARRAY */
            /*-------------------------------*/
            for (i = 0; i < 2*nloci; i++) {
                ind2_rows[(k*2)+j][i] = off_rows[0][i];
            }
        }
        z += 1;
    }

    /* share RNG state with R */
    PutRNGstate();

    /*****************************/
    /* GET UNRELATED INDIVIDUALS */
    /*****************************/
    z = (2*ninds*3);

    for (i = ((2*ninds*3)+(3*ninds)); i < ((2*ninds*4)+(3*ninds)); i++) {
        for (j = 0; j < 2*nloci; j++) {
            ind2_rows[z][j] = ind_rows[i][j];
        }
        z += 1;
    }

#if 0
    /* some debug printing */ {
      FILE *fp, *fp2;
      char ind_fname[] = "fam_sim.individuals";
      char ind2_fname[] = "fam_sim.individuals2";
      if (!(fp = fopen(ind_fname, "w")) || !(fp2 = fopen(ind2_fname, "w"))) {
        error("Can't open files for debug printing.");
        return R_NilValue;
      }
      warn("Created debug files %s and %s.", ind_fname, ind2_fname);

      for (i = 0; i < new_ninds; i++) {
        fprintf(fp, "%d\t", i+1);
        for (j = 0; j < 2*nloci; j++) {
          fprintf(fp, "%d\t", ind_rows[i][j]);
        }
        fprintf(fp, "\n");
      }
      for (i = 0; i < new_ninds2; i++) {
        fprintf(fp2, "%d\t", i+1);
        for (j = 0; j < 2*nloci; j++) {
          fprintf(fp2, "%d\t", ind2_rows[i][j]);
        }
        fprintf(fp2, "\n");
      }
      fclose(fp); fclose(fp2);
    }
#endif

    /**********************/
    /* LOAD OUTPUT MATRIX */
    /**********************/

    /* first, allocate and index an integer matrix */
    PROTECT(ret_sims = allocMatrix(INTSXP, new_ninds2, 2*nloci));
    sims = INTEGER(ret_sims);
    sim_cols = (int**) R_alloc(new_ninds2, sizeof(int*));
    sim_cols[0] = &sims[0];
    for (i = 1; i < 2*nloci; i++) sim_cols[i] = sim_cols[i-1] + new_ninds2;

    /* next, load a transpose of individuals2 into that matrix */
    for (i = 0; i < new_ninds2; i++) {
        for (j = 0; j < 2*nloci; j++) {
            sim_cols[j][i] = ind2_rows[i][j];
        }
    }

    UNPROTECT(1);
    return ret_sims;
}
