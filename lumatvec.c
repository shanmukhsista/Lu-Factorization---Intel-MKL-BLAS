#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <mkl.h>
void PrintMatrix(double * mat, int lda, int nrows, int ncols) {
    int i = 0, j = 0;
    printf("\n");
    for (i = 0; i < nrows; i++) {
        for (j = 0; j < ncols; j++) {
            printf("%.04f\t", *(mat + (i + (lda * j))));
        }
        printf("\n");
    }
}

void PrintMatrixForMatlab(double * mat, int lda, int nrows, int ncols) {
    int i = 0, j = 0;
    printf("\n");
    printf("[");
    for (i = 0; i < nrows; i++) {
        for (j = 0; j < ncols; j++) {
            printf("%.04f", *(mat + (i + (lda * j))));
            if ( j != ncols -1 ){
                printf(",");
            }
        }
        printf(";");
    }
}

void lumatvec_(double *B, int *pldb, int *pm, int *pn, int *ipiv, int *perrflag) {
    //Initialize the pivot vectors 
    //===================================================================================
    //
    // LU factorization with partial pivoting on an m x n matrix B. This version should
    // be based on matrix-vector products, and calls BLAS-2 (BLAS level 2) routines.
    // The parameter small below determines when partial pivoting has failed from the
    // diagonal element B(k,k) being too small to safely divide by. When that happens
    // errflag is set to k and the function bails out. The array B can be rectangular
    // (viz., m /= n). This assumes that the indexing of the matrix B starts from
    // 1, not 0. Otherwise a failure on the first step is indistinguishable from a
    // no-error condition.
    //
    // B  is an array of doubles containing the matrix B in column-major order. It
    //    is (part of) an array of declared leading dimension ldb in the calling
    //    function. 
    //   m is the number of rows in B 
    //   n is the number of columns in B 
    // ipiv is an array of length at least m which contains the pivot indices
    // errflag is what would normally be returned by a function in C. Values:
    // errflag = 0: success
    //         < 0: if errflag = -k, the k-th argument had an illegal value
    //         > 0: if errflag =  k, B(k,k) is too small to rely upon and the 
    //              factorization probably failed.
    //
    // The BLAS requires each vector argument also has an "increment" giving how far
    // apart consecutive entries of the vector are in memory. E.g., if a matrix H is
    // stored in the upper 128 x 128 part of an array declared as G[128][256] in
    // row-major order, consecutive entries in a row of G have an increment of 1.
    // Accessing consecutive entries in a column of G will have an increment of 256.
    //
    // Similarly with column major ordering, if a matrix H stored in the first 128x128
    // entries of an array declared as G(256,128), then accessing elements along a
    // row of G have increment = 256, and elements along a column have increment = 1.
    //
    //===================================================================================

    // Value to use to bail out because pivot is too small
    double small = ((double) 1000.0)*2.216e-16;

    // To avoid messing with stars and ampersands, make local variables
    int ldb = *pldb;
    int n = *pn;
    int m = *pm;
    int debug = 1;
    double ONE = 1.0;
    double ZERO = 0.0;
    double MONE = -1.0;
    double constant = 1.0; 
    // Because quantities like m-k cannot be passed by address to the BLAS, use
    // some temporary variables, e.g., nrows = ldb-j and then pass in &nrows
    int nrows, ncols, whatever;
    *perrflag = 0;
    int i = 0, j = 0;
    int nrow = 0;
    int l = 0;
    int p = 0;
    int nElements = 0 ; 
    char trans = 'N';
    double alpha = -1.0;
    int incX = 0, incY = 0;
    double beta = 1.0;


    //print the matrix 
    /* For getting the A(i,j) - c indices we would use the following formula
     * array Index = (Column Index + (colCount*rowIndex))
     * 2nd row , 3rd column of 5x5 matrix with lda = 5 
     * (2 + ( 5*1)) = 7 B[8]b
     */

    /*Once pivoting is implemented, we would require to perform a few operations
     only after the 1st column*/
    //Matrix Multiply using blas
    //    cblas_dgemv (102,       //102 for column major order
    //                 'N',  //'N' for normal equation
    //                 m,               //Number of rows in Matrix
    //                 n,               //Number of Columns in Matrix    
    //                  -1,            //constant -> -1 in this case
    //                 B,
    //                 ldb,             //leading dimension of A
    //                 B,               // array to store results
    //                 m,            //
    //                 1 
    //                1);
    //        
    //Write all the array elements to the screen 
    for (i = 0; i < n; i++) {
        /*Outer loop is for looping through each columns 1 to n*/
        /*Matrix A has an index of i+1 to m*/
        /*Vector v has an index of i to n*/

        int index = 0;
        if (i > 0) {

            /* 
                 A(k:m, k)   = A(k:m, k) - A(k:m, 1:k-1)*A(1:k-1, k); 
             * 
             * for matrix A number of rows = (m - i )
             *      z          number of columns = (i)
             * for vector x number of rows = (i+1)
             *              number of columns =  1 ( ith Column)
             * for vector y - result
             *              number of rows = (m-i)) 
             *                        //Index : 
             *              number of columns = 1 (ith Column)
             */
            //Get the elements from i'th row using column major order
            /*y := alpha*A*x + beta*y Operation using dgemv*/
            nrows = m - i;
            ncols = i;
            incX = 1;
            incY = 1;
            dgemv(//102 for column major order
                  &trans, //'N' for normal equation
                  &nrows, //Number of rows in Matrix
                  &ncols, //Number of Columns in Matrix    
                  &alpha, //constant -> -1 in this case
                  B + (i), //A(k:n,1:k-1) x   
                  pldb,
                  B + (ldb * i), //vector x 
                  &incX, //increment x //leading dimension of A
                  &beta,
                  B + (i + (ldb * i)), // y - array to store results
                  &incY //
                  );

            //Find the largest element in the column for pivoting.
            /*CBLAS_INDEX cblas_idamax (const MKL_INT N, const double *X, const MKL_INT incX);*/
 
        }


        //Get the maximum value for the column vector using blas function
        //cblas_idmax

       
        nrows = m-i ; 
        incX = 1 ; 
        index = idamax(&nrows,
                             B + (i + (ldb * i)),
                             &incX) - 1;
         /*Get the value of the pivot element. If this value is less than small,
         then assign an error code and return. */
        if (fabs(B[i + (ldb * i)] + index) < small) {
            /* The indexing starts from one. So if there is an error in the 
             1st row 1st column pivot, *perrflag will hold a value of 1*/
            *perrflag = i + 1+index; //Error in the (k,K pivot) 
            return ;
        }
        
        if (index != 0) {
            //Swap both the rows in the matrix.
            /*Store the pivot elements in the ipiv array. 
             Index with respect to the entire matrix is stored in ipiv.
             */
            nElements = n ;  
            int nrow = 0;
    int l = 0;
    int p = 0;
                dswap(&nElements,
                        B + (i), 
                        pldb , 
                        B + (i) + index,
                        pldb);
            ipiv[i] = i + index + 1;
           
        } else {
            //Store the pivot index in the ipiv array. This has to be used in the 
            //matmat implementation for pivoting the remaining rows. 
           
            ipiv[i] = i + 1 ;
        }
        //printf("Pivot for column %d is %d", i+1 , ipiv[i]);
        if (i > 0) {
            //Step 2
            /* 
                 A(k:m, k)   = A(k:m, k) - A(k:m, 1:k-1)*A(1:k-1, k); 
             *   A(k, k+1:n) = A(k, k+1:n) - A(k, 1:k-1)*A(1:k-1, k+1:n);
             * 
             * for matrix A number of rows = i
             *              number of columns = ( n - i + 1 )
             * for vector x number of rows = 1 ( ith Column index)
             *              number of columns =  k - 1 ; 
             * for vector y - result
             *              number of rows = 1  
             *                        //Index : 
             *              number of columns = (n-i+1)
             */
            //Get the elements from i'th row using column major order
            /*y := alpha*A*x + beta*y Operation using dgemv*/
                nrows = i;
                ncols = (n - i - 1 );
            incX = ldb;
            incY = ldb;   
            char trans1 = 'T';
                dgemv( //102 for column major order
                        &trans1, //'N' for normal equation
                        &nrows, //Number of rows in Matrix
                        &ncols, //Number of Columns in Matrix    
                        &alpha, //constant -> -1 in this case
                        B + (ldb * (i + 1)), //A(k,k+1:n) x   
                        pldb,
                        B + i, //vector x 
                        &incX, //increment x //leading dimension of A
                        &beta,
                        B + (i + (ldb * (i + 1))), // y - array to store results
                        &incY //
                        );
        }
        /*
         *    A(k,k+1:n) = A(k,k+1:n) - A(k,1:k-1)*A(1:k-1,k+1:n)
         */
        //Perform the scaling operation
        /*
         *   A(k+1:n,k) = A(k+1:n,k)/A(k,k)
         */
        //Scaling using blas
        constant =   (1.0 / B[i + (ldb * i)]);
        incX = 1 ; 
        nrows  = m - i - 1 ; 
        dscal(&nrows, //Number of elements to scale.
                   &constant,
                    B + (i + 1 + (ldb * i)),
                    &incX);

        //        for (j = i + 1; j <= m - 1; j++) {
        //            B[j + (m * i)] = B[j + (m * i)] / B[i + (m * i)];
        //        }

    
    }
}
//Matrix Vector Multiplication 


