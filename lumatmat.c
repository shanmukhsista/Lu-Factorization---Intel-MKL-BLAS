//#include <solvers.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "mkl.h"

void lumatmat_(double *A, int *plda, int *pn, int *ipiv, int *perrflag, int *pnb) {
    //================================================================================
    //
    // Perform a BLAS-3 LU factorization of a square n x n matrix stored in the array
    // A. As a 2D array, A is stored in column-major order with leading declared
    // dimension lda. On exit, array A contains the L and U factors of the matrix A.
    //
    // Compute the LU factorization of a general n-by-n matrix A using partial
    // pivoting via row swaps. This implementation should use matrix-matrix products, a
    // BLAS level 3 operation and hence able to get and maintain a high fraction of a
    // machine's theoretical peak. The base version is the rank-1 update form, so the
    // matrix-matrix products are of the form of rank-nb updates.
    // 
    // The parameter small below determines when partial pivoting has failed from the
    // diagonal element A(k,k) being too small to safely divide by. When that happens
    // errflag is set to k and lumatmat returns.
    //
    //   A is an array of doubles containing the matrix A in column-major order.
    //   m is the number of rows in A 
    //   n is the number of columns in A 
    // ipiv is an array of length at least n which contains the pivot indices
    // errflag is what would normally be returned by a function in C. Values:
    // errflag = 0: success
    //         < 0: if errflag = -k, the k-th argument had an illegal value
    //         > 0: if errflag =  k, A(k,k) is too small to rely upon and the 
    //              factorization probably failed.
    // nb is the blocksize to use. Using nb = n or nb > n will lead to a BLAS-2 solve
    // using lumatvec.
    //
    //===================================================================================

    // value to use to bail out because pivot is too small
    const double small = ((double) 1000.0)*2.216e-16;
    int lda = *plda;
    int n = *pn;
    int nb = *pnb;
    int * tempCols;
    *perrflag = 0;
    double * b;
    int rows = 0;
    int j = 0;
    int cols = 0;
    int *tempRows;
    double alpha = -1.0, beta = 1.0; 
    char uplo, diag ; 
    int rowsb , colsb, rowsc ; 
    int colStartIndex = 0, iterationCount = 0;
    char side = 'L'; 
    // Sample invocation of lu4:
    //lumatvec_(A, &lda, &n, &n, ipiv, perrflag);
    char notrans = 'N'; 
    char trans = 'T';
    int i = 0, k = 0;
    //Create a loop for block sizes. 
    /* IN this case, i am creating two loops for the blocks and the remaining 
     elements separately..*/

    int remainingElements = n % nb;
    int nblocks = 0;
    int nElements = 0, incX = 0;
    /*
     * If the block size is greater than the total number of rows , then set 
     * the block size equal to the number of rows
        
     */
    if (n % nb == 0) {
        nblocks = (n - (n % nb)) / nb - 1;
    } else {
        nblocks = (n - (n % nb)) / nb;

    }



    //Print matrix before starting any operation
    if (nblocks < 1) {
        nblocks = 0;
        nb = n;
    }
    if (nb == 1 || nb >= n) {
        //Call Lumatvec_ directly without proceeding further. 
        rows = n;
        cols = n;
        lumatvec_(A, plda, &rows, &cols, ipiv, perrflag);
        return;
    }

    for (i = 0; i <= nblocks; i++) {
        /*This loop will perform all the block operations in the matrix.*/
        /*Perform a Matrix Vector Update for the elements*/
        //update the pivots before performing an lumatvec. 

        for (k = i * nb; k < nb * (i + 1); k++) {
            if (k < n) {
                ipiv[k] = k + 1;
            }
        }

        if (i > 0) {

            //Apply transformations for the previous values. 
            /*
             * C := alpha*op(A)*op(B) + beta*C,
             * alpha = -1 ; 
             * beta = 1 
             * A - previous block matrix 
             * C - Result 
             * B - Matrix above C 
             */
            rows = (n - (nb * i));

            if (nb * (i + 1) < n) {
                cols = nb;

            } else {
                cols = n - (nb * (i));


            }

            rows = (n - (nb * i));
            rowsb = (nb*i) ; 

            //rows = ((nb * i));

            dgemm(
                     
                        &notrans,
                        &notrans,
                        &rows, //rows in A 
                        &cols, //cols in B 
                        &rowsb, //cols in A  and rows in B 
                        &alpha,
                        A + (i * nb),
                        plda,
                        A + ((n * nb * (i))),
                        plda,
                        &beta,
                        A + ((n * nb * (i)) + (i * nb)),
                        plda
                        );

            rows = (n - (nb * i));


        }
        rows = (n - (nb * i));

        if (nb * (i + 1) > n) {
            cols = n - (nb * (i));

        } else {
            cols = nb;

        }

        //Do the lumatvec only if the block          
        lumatvec_(A + ((n * nb * i) + (i * nb)), plda, &rows, &cols, ipiv + (nb * i), perrflag);
        //Check if any error was encountered in lumatvec. 
        if ( *perrflag != 0 ){
            *perrflag = *perrflag + (nb*i) ; 
            return ; 
        }
        /*
         * Perform pivoting for the remaining columns of the column matrix.
         *  
         */
        //Loop through all elements in pivot vector
        /*
         since we have nblocks of size nb , we will loop through all the columns
         * of the current block and test the pivot elements. 
         */
        colStartIndex = colStartIndex + (i * nb);
        if (nb * (i + 1) > n) {
            iterationCount = n - (nb * (i));

        } else {
            iterationCount = nb;

        }
        for (j = 0; j < iterationCount; j++) {
            //update the pivots for each column : 
            //if a call to lumatvec is not made, then the pivot value will be 
            //the same as the column number starting from 1 
            //column value would be equal to (nb*i) + j 
            if (*(ipiv + (nb * i) + j) - 1 != j) {
                nElements = nb*i;
                incX = n;
                
                //Swap the previous columns and the next columns for the ipiv[j] value
                //ipiv[j] holds the row index for pivoting

                //Swap previous columns i = 0 to i = nb*i
                //Number of previous elements to swap = nb* i 
                if (i > 0) {
                    dswap(
                          &nElements, //num of elements in x and y 
                          A + (nb * i) + j,
                          &incX,
                          A + (nb * i) + *(ipiv + (nb * i) + j) - 1,
                          &incX
                          );
                }
                nElements = n - (nb * (i + 1));
                dswap(
                      &nElements, //num of elements in x and y 
                      A + (nb * (i + 1) * n) + (nb * i) + j,
                      &incX,
                      A + (nb * (i + 1) * n) + (nb * i) + *(ipiv + (nb * i) + j) - 1,
                      &incX

                      );

               
                //Now swap the remaining elements after the current block
                //nElements = n - (nb * (i+1))
                //Starting Element
                *(ipiv + (nb * i) + j) = *(ipiv + (i * nb) + j) + (i * nb);

            } else {
                *(ipiv + (nb * i) + j) = *(ipiv + (i * nb) + j) + (nb * i);
            }
        }
        if ((i + 1) * nb < n) {
            if (i > 0) {
                //C := alpha*op(A)*op(B) + beta*C,
                /*
             
                 */
                rows = nb;
                cols = (n - nb * (i + 1)); 
                rowsb = (nb*i) ; 
                dgemm(
                            
                            &notrans,
                            &notrans,
                            &rows, //Number of rows in A
                            &cols, //Number of columns in B 
                            &rowsb, //Number of columns in A and rows in B
                            &alpha,
                            A + (nb * (i)), //Matrix A
                            plda,
                            A + (nb * (i + 1) * n), //Matrix B
                            plda,
                            &beta,
                            A + ((n * nb * (i + 1)) + ((i) * nb)), //Matrix C
                            plda
                            );
            }
            //Update the trailing part of the matrix. 
            //Solve a triangular system of equations 
            //Ax = b - index = b = (nb * n *i) + nb 
            // 
            //Solve for all the remaining matrix 
            
           
            /*Use dtrsm to solve the matrix*/
            uplo = 'L'; 
                diag = 'U';
                rows = nb ;
                incX = 1;
                cols = (n - nb * (i + 1));
            dtrsm( 
                   &side, 
                   &uplo,
                   &notrans,
                   &diag, 
                   &rows,
                   &cols,
                   &beta,
                   A + ((n * nb * (i))) + (nb * i),
                  plda,
                  A + ((n * nb * (i + 1)) + (i * nb)) ,
                  plda
                    );
        }
    }
    return;

}
