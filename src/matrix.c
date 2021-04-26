#include "matrix.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

// Include SSE intrinsics
#if defined(_MSC_VER)
#include <intrin.h>
#elif defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
#include <immintrin.h>
#include <x86intrin.h>
#endif

/* Below are some intel intrinsics that might be useful
 * void _mm256_storeu_pd (double * mem_addr, __m256d a)
 * __m256d _mm256_set1_pd (double a)
 * __m256d _mm256_set_pd (double e3, double e2, double e1, double e0)
 * __m256d _mm256_loadu_pd (double const * mem_addr)
 * __m256d _mm256_add_pd (__m256d a, __m256d b)
 * __m256d _mm256_sub_pd (__m256d a, __m256d b)
 * __m256d _mm256_fmadd_pd (__m256d a, __m256d b, __m256d c)
 * __m256d _mm256_mul_pd (__m256d a, __m256d b)
 * __m256d _mm256_cmp_pd (__m256d a, __m256d b, const int imm8)
 * __m256d _mm256_and_pd (__m256d a, __m256d b)
 * __m256d _mm256_max_pd (__m256d a, __m256d b)
*/

/* Generates a random double between low and high */
double rand_double(double low, double high) {
    double range = (high - low);
    double div = RAND_MAX / range;
    return low + (rand() / div);
}

/* Generates a random matrix */
void rand_matrix(matrix *result, unsigned int seed, double low, double high) {
    srand(seed);
    for (int i = 0; i < result->rows; i++) {
        for (int j = 0; j < result->cols; j++) {
            set(result, i, j, rand_double(low, high));
        }
    }
}

/*
 * Allocates space for a matrix struct pointed to by the double pointer mat with
 * `rows` rows and `cols` columns. You should also allocate memory for the data array
 * and initialize all entries to be zeros. `parent` should be set to NULL to indicate that
 * this matrix is not a slice. You should also set `ref_cnt` to 1.
 * You should return -1 if either `rows` or `cols` or both have invalid values, or if any
 * call to allocate memory in this function fails. If you don't set python error messages here upon
 * failure, then remember to set it in numc.c.
 * Return 0 upon success.
 */
int allocate_matrix(matrix **mat, int rows, int cols) {
    /* TODO: YOUR CODE HERE */
    if (rows < 1 || cols < 1) {
        PyErr_SetString(PyExc_RuntimeError, "Invalid matrix dimensions");
        return -1;
    }
    *mat = malloc(sizeof(matrix));
        if (!(*mat)) {
            PyErr_SetString(PyExc_RuntimeError, "Matrix structure allocation failed");
            return -1;
        }

    const int NUM_THREADS = 4;
    omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel
    {
        int id = omp_get_thread_num();
        switch(id) {
            case 0:
                (*mat)->rows = rows;
                break;
            case 1:
                (*mat)->cols = cols;
                break;
            case 2:
                (*mat)->data = calloc(rows * cols, sizeof(double));
                break;
            default:
                (*mat)->ref_cnt = 1;
                (*mat)->parent = NULL;
        }
    }
    if (!((*mat)->data)) {
        PyErr_SetString(PyExc_RuntimeError, "Data allocation failed");
        free(mat);
        return -1;
    }
    return 0;
}

/*
 * Allocates space for a matrix struct pointed to by `mat` with `rows` rows and `cols` columns.
 * Its data should point to the `offset`th entry of `from`'s data (you do not need to allocate memory)
 * for the data field. `parent` should be set to `from` to indicate this matrix is a slice of `from`.
 * You should return -1 if either `rows` or `cols` or both are non-positive or if any
 * call to allocate memory in this function fails.
 * If you don't set python error messages here upon failure, then remember to set it in numc.c.
 * Return 0 upon success.
 */
int allocate_matrix_ref(matrix **mat, matrix *from, int offset, int rows, int cols) {
    /* TODO: YOUR CODE HERE */
    if (rows < 1 || cols < 1) {
        PyErr_SetString(PyExc_RuntimeError, "Invalid matrix dimensions");
        return -1;
    }
    *mat = malloc(sizeof(matrix));
    if (!(*mat)) {
        PyErr_SetString(PyExc_RuntimeError, "Matrix structure allocation failed");
        return -1;
    }

    (*mat)->rows = rows;
    (*mat)->cols = cols;
    (*mat)->data = from -> data + offset;
    if (!((*mat)->data)) {
        PyErr_SetString(PyExc_RuntimeError, "Data pointer failed");
        return -1;
    }
    from->ref_cnt += 1;
    (*mat)->ref_cnt = 1;
    (*mat)->parent = from;
    return 0;
}

/*
 * You need to make sure that you only free `mat->data` if `mat` is not a slice and has no existing slices,
 * or if `mat` is the last existing slice of its parent matrix and its parent matrix has no other references
 * (including itself). You cannot assume that mat is not NULL.
 */
void deallocate_matrix(matrix *mat) {
    /* TODO: YOUR CODE HERE */
    if (mat) {
        mat->ref_cnt -= 1;
        if (mat->ref_cnt <= 0) {
            //has no children
            if(!mat->parent) {
                free(mat->data);
            } else {
                mat->parent->ref_cnt -= 1;
                if (mat->parent->ref_cnt <= 0) {
                    //if parent is dead and has no children, deallocate them
                    deallocate_matrix(mat->parent);
                }
            }
            free(mat);
        } //if there are children, nothing else happens
    }
}

/*
 * Returns the double value of the matrix at the given row and column.
 * You may assume `row` and `col` are valid.
 */
double get(matrix *mat, int row, int col) {
    /* TODO: YOUR CODE HERE */
    return mat->data[row * mat->cols + col];
}

/*
 * Sets the value at the given row and column to val. You may assume `row` and
 * `col` are valid
 */
void set(matrix *mat, int row, int col, double val) {
    /* TODO: YOUR CODE HERE */
    mat->data[row * mat->cols + col] = val;
}

/*
 * Sets all entries in mat to val
 */
void fill_matrix(matrix *mat, double val) {
    /* TODO: YOUR CODE HERE */
    #pragma omp parallel for
    for (int i = 0; i < mat->rows * mat->cols / 4 * 4; i+=4) {
        mat->data[i] = val;
        mat->data[i+1] = val;
        mat->data[i+2] = val;
        mat->data[i+3] = val;
    }
    #pragma omp parallel for
    for (int i = mat->rows * mat->cols / 4 * 4; i < mat->rows * mat->cols; i+=1) {
        mat->data[i] = val;
    }
}

/*
 * Store the result of adding mat1 and mat2 to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int add_matrix(matrix *result, matrix *mat1, matrix *mat2) {
    /* TODO: YOUR CODE HERE */
    if (mat1->rows == mat2->rows && result->rows == mat1->rows && mat1->cols == mat2->cols && mat2->cols == result->cols) {
        #pragma omp parallel for
        for(unsigned int i = 0; i < (mat1->rows * mat1->cols)/16 * 16; i += 16) {
			__m256d vector1 = _mm256_loadu_pd(&(mat1->data[i]));
			__m256d vector2 = _mm256_loadu_pd(&(mat2->data[i]));

			__m256d tempVector = _mm256_add_pd(vector1, vector2);
			_mm256_storeu_pd(&(result->data[i]), tempVector);
            
			vector1 = _mm256_loadu_pd(&(mat1->data[i+4]));
			vector2 = _mm256_loadu_pd(&(mat2->data[i+4]));
			tempVector = _mm256_add_pd(vector1, vector2);
			_mm256_storeu_pd(&(result->data[i+4]), tempVector);

            vector1 = _mm256_loadu_pd(&(mat1->data[i+8]));
			vector2 = _mm256_loadu_pd(&(mat2->data[i+8]));
			tempVector = _mm256_add_pd(vector1, vector2);
			_mm256_storeu_pd(&(result->data[i+8]), tempVector);

            vector1 = _mm256_loadu_pd(&(mat1->data[i+12]));
			vector2 = _mm256_loadu_pd(&(mat2->data[i+12]));
			tempVector = _mm256_add_pd(vector1, vector2);
			_mm256_storeu_pd(&(result->data[i+12]), tempVector);
		}

        for (int i = (mat1->rows * mat1->cols)/16 * 16; i < mat1->rows * mat1->cols; i++) {
            result->data[i] = mat1->data[i] + mat2->data[i];
        }
        return 0;
    }
    return -1;
}

/*
 * Store the result of subtracting mat2 from mat1 to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int sub_matrix(matrix *result, matrix *mat1, matrix *mat2) {
    /* TODO: YOUR CODE HERE */
    if (mat1->rows == mat2->rows && result->rows == mat1->rows && mat1->cols == mat2->cols && mat2->cols == result->cols) {
        #pragma omp parallel for
        for(unsigned int i = 0; i < (mat1->rows * mat1->cols)/16 * 16; i += 16) {
			__m256d vector1 = _mm256_loadu_pd(&(mat1->data[i]));
			__m256d vector2 = _mm256_loadu_pd(&(mat2->data[i]));

			__m256d tempVector = _mm256_sub_pd(vector1, vector2);
			_mm256_storeu_pd(&(result->data[i]), tempVector);
            
			vector1 = _mm256_loadu_pd(&(mat1->data[i+4]));
			vector2 = _mm256_loadu_pd(&(mat2->data[i+4]));
			tempVector = _mm256_sub_pd(vector1, vector2);
			_mm256_storeu_pd(&(result->data[i+4]), tempVector);

            vector1 = _mm256_loadu_pd(&(mat1->data[i+8]));
			vector2 = _mm256_loadu_pd(&(mat2->data[i+8]));
			tempVector = _mm256_sub_pd(vector1, vector2);
			_mm256_storeu_pd(&(result->data[i+8]), tempVector);

            vector1 = _mm256_loadu_pd(&(mat1->data[i+12]));
			vector2 = _mm256_loadu_pd(&(mat2->data[i+12]));
			tempVector = _mm256_sub_pd(vector1, vector2);
			_mm256_storeu_pd(&(result->data[i+12]), tempVector);
		}
        #pragma omp parallel for
        for (int i = (mat1->rows * mat1->cols)/16 * 16; i < mat1->rows * mat1->cols; i++) {
            result->data[i] = mat1->data[i] - mat2->data[i];
        }
        return 0;
    }
    return -1;
}

/*
 * Store the result of multiplying mat1 and mat2 to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 * Remember that matrix multiplication is not the same as multiplying individual elements.
 */
int mul_matrix(matrix *result, matrix *mat1, matrix *mat2) {
    /* TODO: YOUR CODE HERE */
    
    matrix *mat2t;
    allocate_matrix(&mat2t, mat2->cols, mat2->rows);

    const int NUM_THREADS = 8;
    omp_set_num_threads(NUM_THREADS);
    //transposes mat2 for faster multiplication
    #pragma omp parallel
    {
        int id = omp_get_thread_num();

        int lowerx = (id % 4) * mat2->cols / 4;
        int upperx = ((id % 4) + 1) * mat2->cols / 4;
        int lowery = (id / 4) * mat2->rows / 2;
        int uppery = ((id / 4) + 1) * mat2->rows / 2;
        for (int startx = lowerx; startx < upperx; startx += 16) {
            for (int starty = lowery; starty < uppery; starty += 16) {
                for (int x = startx; x < startx + 16; x++) {
                    for (int y = starty; y < starty + 16; y++) {
                        if (y < uppery && x < upperx) mat2t->data[x * mat2t->cols + y] = mat2->data[y*mat2->cols + x];
                    }
                }
            }
        }
        
    }
    fill_matrix(result, 0);

    //performs multiplication, column by column
    #pragma omp parallel
    {
        int id = omp_get_thread_num();

        int lowerx = (id) * mat1->rows/8;
        int upperx = (id + 1) * mat1->rows/8;
        //int lowerx = (id % 8) * mat2t->rows/8;
        //int upperx = ((id % 8) + 1) * mat2t->rows/8;
        //int lowery = (id / 4) * mat1->rows/2;
        //int uppery = ((id / 4) + 1) * mat1->rows/2;
        for (int part = 0; part < mat1->cols; part +=32) {
            for (int mat1row = lowerx; mat1row < upperx; mat1row++) {
                for (int mat2col = 0; mat2col < mat2t->rows; mat2col++) {
                    __m256d summedVector = _mm256_setzero_pd();
                    for (unsigned int i = part; i < part + 32; i+=16) {
                        if (i + 16 < mat1->cols) {
                            __m256d vectorVals = _mm256_loadu_pd(&(mat1->data[mat1row * mat1->cols + i]));
                            __m256d vector2 = _mm256_loadu_pd(&(mat2t->data[mat2col * mat2t->cols + i]));
                            summedVector = _mm256_fmadd_pd(vector2, vectorVals, summedVector);

                            vectorVals = _mm256_loadu_pd(&(mat1->data[mat1row * mat1->cols + i + 4]));
                            vector2 = _mm256_loadu_pd(&(mat2t->data[mat2col * mat2t->cols + i + 4]));
                            summedVector = _mm256_fmadd_pd(vector2, vectorVals, summedVector);

                            vectorVals = _mm256_loadu_pd(&(mat1->data[mat1row * mat1->cols + i + 8]));
                            vector2 = _mm256_loadu_pd(&(mat2t->data[mat2col * mat2t->cols + i + 8]));
                            summedVector = _mm256_fmadd_pd(vector2, vectorVals, summedVector);

                            vectorVals = _mm256_loadu_pd(&(mat1->data[mat1row * mat1->cols + i + 12]));
                            vector2 = _mm256_loadu_pd(&(mat2t->data[mat2col * mat2t->cols + i + 12]));
                            summedVector = _mm256_fmadd_pd(vector2, vectorVals, summedVector);
                        }
                        else if(i < mat1->cols) {
                            for (int j = 0; j < mat1->cols - i; j++) {
                                result->data[mat1row * result->cols + mat2col] += mat1->data[mat1row * mat1->cols + i + j]*mat2t->data[mat2col * mat2t->cols + i + j];
                            }
                            
                        }
                        
                    }
                    double *summedArray = malloc(32);
                    _mm256_storeu_pd(summedArray, summedVector);
                    for (unsigned int i = 0; i < 4; i++) {
                        result->data[mat1row * result->cols + mat2col] += summedArray[i];
                    }
                    free(summedArray);
                }
            }
        }
    }
    deallocate_matrix(mat2t);
    return 0;
}

/*
 * Store the result of raising mat to the (pow)th power to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 * Remember that pow is defined with matrix multiplication, not element-wise multiplication.
 */
int pow_matrix(matrix *result, matrix *mat, int pow) {
    /* TODO: YOUR CODE HERE */
    if (pow == 0) {
        #pragma omp parallel for
        for (int i = 0; i < mat->rows * mat->cols; i++) {
            if (i % (mat->cols + 1) == 0) {
                result->data[i] = 1;
            } else {
                result->data[i] = 0;
            }
        }
        return 0;
    }
    else if (pow == 1) {
        #pragma omp parallel for
        for (int i = 0; i < mat->rows * mat->cols; i++) {
            result->data[i] = mat->data[i];
        }
        return 0;
    }
    else if (pow % 2 == 0) {
        matrix *temp = NULL;
        allocate_matrix(&temp, result->rows, mat->cols);
        int a = pow_matrix(temp, mat, pow/2);
        mul_matrix(result, temp, temp);
        deallocate_matrix(temp);
        return a;
    }
    else {
        matrix *temp = NULL;
        allocate_matrix(&temp, result->rows, mat->cols);
        int a = pow_matrix(temp, mat, pow - 1);
        mul_matrix(result, temp, mat);
        deallocate_matrix(temp);
        return a;
    }
}

/*
 * Store the result of element-wise negating mat's entries to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int neg_matrix(matrix *result, matrix *mat) {
    /* TODO: YOUR CODE HERE */
    #pragma omp parallel for
    for (int i = 0; i < mat->rows * mat->cols / 4 * 4; i+=4) {
        result->data[i] = mat->data[i]*-1;
        result->data[i+1] = mat->data[i+1]*-1;
        result->data[i+2] = mat->data[i+2]*-1;
        result->data[i+3] = mat->data[i+3]*-1;
    }
    #pragma omp parallel for
    for (int i = mat->rows * mat->cols / 4 * 4; i < mat->rows * mat->cols; i++) {
        result->data[i] = mat->data[i]*-1;
    }
    return 0;
}

/*
 * Store the result of taking the absolute value element-wise to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int abs_matrix(matrix *result, matrix *mat) {
    /* TODO: YOUR CODE HERE */
    #pragma omp parallel for
    for (int i = 0; i < mat->rows * mat->cols / 4 * 4; i+=4) {
        result->data[i] = fabs(mat->data[i]);
        result->data[i+1] = fabs(mat->data[i+1]*-1);
        result->data[i+2] = fabs(mat->data[i+2]*-1);
        result->data[i+3] = fabs(mat->data[i+3]*-1);
    }
    #pragma omp parallel for
    for (int i = mat->rows * mat->cols / 4 * 4; i < mat->rows * mat->cols; i++) {
        result->data[i] = fabs(mat->data[i]);
    }
    return 0;
}

