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
    *mat = calloc((rows * cols) + 4, sizeof(double));
    if (!(*mat)) {
        PyErr_SetString(PyExc_RuntimeError, "Matrix structure allocation failed");
        return -1;
    }
    (*mat)->rows = rows;
    (*mat)->cols = cols;
    (*mat)->data = &(*mat)->placeholder;
    (*mat)->ref_cnt = 1;
    (*mat)->parent = NULL;
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
            if(mat->parent) {
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
    __m256d set = _mm256_set1_pd(val);
    double *matrix = mat->data;
    int size = mat->cols * mat->rows;

    #pragma omp parallel for
    for (int i = 0; i < size / 16 * 16; i+=16) {
        _mm256_storeu_pd(matrix + i, set);
		_mm256_storeu_pd(matrix + i + 4, set);
		_mm256_storeu_pd(matrix + i + 8, set);
		_mm256_storeu_pd(matrix + i + 12, set);
    }
    #pragma omp parallel for
    for (int i = size / 16 * 16; i < size / 4 * 4; i+=4) {
        _mm256_storeu_pd(matrix + i, set);
    }
    #pragma omp parallel for
    for (int i = size / 4 * 4; i < size; i++) {
        matrix[i] = val;
    }
}

/*
 * Store the result of adding mat1 and mat2 to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int add_matrix(matrix *result, matrix *mat1, matrix *mat2) {
    /* TODO: YOUR CODE HERE */
    if (mat1->rows == mat2->rows && result->rows == mat1->rows && mat1->cols == mat2->cols && mat2->cols == result->cols) {
        
        double *matrix1 = mat1->data;
        double *matrix2 = mat2->data;
        double *resultmatrix = result->data;
        
        int size = mat2->cols * mat2->rows;
        #pragma omp parallel for
        for(unsigned int i = 0; i < size/16 * 16; i += 16) {
			_mm256_storeu_pd(resultmatrix + i, _mm256_add_pd(_mm256_loadu_pd(matrix1 + i), _mm256_loadu_pd(matrix2 + i)));
			_mm256_storeu_pd(resultmatrix + i + 4, _mm256_add_pd(_mm256_loadu_pd(matrix1 + i + 4), _mm256_loadu_pd(matrix2 + i + 4)));
			_mm256_storeu_pd(resultmatrix + i + 8, _mm256_add_pd(_mm256_loadu_pd(matrix1 + i + 8), _mm256_loadu_pd(matrix2 + i + 8)));
			_mm256_storeu_pd(resultmatrix + i + 12, _mm256_add_pd(_mm256_loadu_pd(matrix1 + i + 12), _mm256_loadu_pd(matrix2 + i + 12)));
		}
        #pragma omp parallel for
        for (int i = size / 16 * 16; i < size / 4 * 4; i+=4) {
            _mm256_storeu_pd(resultmatrix + i, _mm256_add_pd(_mm256_loadu_pd(matrix1 + i), _mm256_loadu_pd(matrix2 + i)));
        }
        #pragma omp parallel for
        for (int i = size / 16 * 16; i < size; i++) {
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
        
        double *matrix1 = mat1->data;
        double *matrix2 = mat2->data;
        double *resultmatrix = result->data;
        
        int size = mat2->rows * mat2->cols;
        #pragma omp parallel for
        for(unsigned int i = 0; i < size/16 * 16; i += 16) {
			_mm256_storeu_pd(resultmatrix + i, _mm256_sub_pd(_mm256_loadu_pd(matrix1 + i), _mm256_loadu_pd(matrix2 + i)));
			_mm256_storeu_pd(resultmatrix + i + 4, _mm256_sub_pd(_mm256_loadu_pd(matrix1 + i + 4), _mm256_loadu_pd(matrix2 + i + 4)));
			_mm256_storeu_pd(resultmatrix + i + 8, _mm256_sub_pd(_mm256_loadu_pd(matrix1 + i + 8), _mm256_loadu_pd(matrix2 + i + 8)));
			_mm256_storeu_pd(resultmatrix + i + 12, _mm256_sub_pd(_mm256_loadu_pd(matrix1 + i + 12), _mm256_loadu_pd(matrix2 + i + 12)));
		}
        #pragma omp parallel for
        for (int i = size / 16 * 16; i < size / 4 * 4; i+=4) {
            _mm256_storeu_pd(resultmatrix + i, _mm256_sub_pd(_mm256_loadu_pd(matrix1 + i), _mm256_loadu_pd(matrix2 + i)));
        }
        #pragma omp parallel for
        for (int i = size / 16 * 16; i < size; i++) {
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
    
    double *mat2t = malloc(mat2->cols * mat2->rows * sizeof(double));
    if (!mat2t) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to allocate memory");
        return -1;
    }
    double *matrix1 = mat1->data;
    double *matrix2 = mat2->data;
    double *resultmatrix = result->data;
    int mat2cols = mat2->cols;
    int termsperentry = mat2->rows;
    int mat1rows = mat1->rows;

    //transpose mat2 for faster multiplication
    #pragma omp parallel for
    for (int startx = 0; startx < mat2->cols; startx += 64) {
        for (int starty = 0; starty < mat2->rows; starty += 64) {
            for (int x = startx; x < startx + 64 && x < mat2->cols; x++) {
                for (int y = starty; y < starty + 64 && y < mat2->rows; y++) {
                    mat2t[x * termsperentry + y] = matrix2[y * mat2cols + x];
                }
            }
        }
    }

    
    //performs multiplication, column by column
    #pragma omp parallel for
    for (int index = 0; index < mat1rows * mat2cols; index++) {
        int mat1row = index / mat2cols;
        int mat2col = index % mat2cols;
        int mat1rowindex = mat1row * termsperentry;
        int mat2colindex = mat2col * termsperentry;
        __m256d summedVector = _mm256_setzero_pd();
        for (int i = 0; i < termsperentry/16 * 16; i+=16) {
            summedVector = _mm256_fmadd_pd(_mm256_loadu_pd(matrix1 + mat1rowindex + i), _mm256_loadu_pd(mat2t + mat2colindex + i), summedVector);
            summedVector = _mm256_fmadd_pd(_mm256_loadu_pd(matrix1 + mat1rowindex + i + 4), _mm256_loadu_pd(mat2t + mat2colindex + i + 4), summedVector);
            summedVector = _mm256_fmadd_pd(_mm256_loadu_pd(matrix1 + mat1rowindex + i + 8), _mm256_loadu_pd(mat2t + mat2colindex + i + 8), summedVector);
            summedVector = _mm256_fmadd_pd(_mm256_loadu_pd(matrix1 + mat1rowindex + i + 12), _mm256_loadu_pd(mat2t + mat2colindex + i + 12), summedVector);
        }
        for (int i = termsperentry/16 * 16; i < termsperentry/4 * 4; i+=4) {
            summedVector = _mm256_fmadd_pd(_mm256_loadu_pd(matrix1 + mat1rowindex + i), _mm256_loadu_pd(mat2t + mat2colindex + i), summedVector);
        }

        resultmatrix[index] = summedVector[0] + summedVector[1] + summedVector[2] +summedVector[3];

        for (int i = termsperentry/4 * 4; i < termsperentry; i++) {
            resultmatrix[index] += matrix1[mat1rowindex + i]*mat2t[mat2colindex + i];
        }
    }
    free(mat2t);
    return 0;
}

/*
 * Store the result of raising mat to the (pow)th power to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 * Remember that pow is defined with matrix multiplication, not element-wise multiplication.
 */
int pow_matrix(matrix *result, matrix *mat, int pow) {
    /* TODO: YOUR CODE HERE */
    double *resultmatrix = result->data;
    int size = mat->rows;

    if (pow == 0) {
        #pragma omp parallel for
        for (int i = 0; i < size * size; i++) {
            if (i % (size + 1) == 0) {
                resultmatrix[i] = 1;
            } else {
                resultmatrix[i] = 0;
            }
        }
        return 0;
    }
    else if (pow % 2 == 0) {
        matrix *temp;
        allocate_matrix(&temp, size, size);
        int a = pow_matrix(temp, mat, pow/2);
        mul_matrix(result, temp, temp);
        deallocate_matrix(temp);
        return a;
    }
    else {
        matrix *temp;
        allocate_matrix(&temp, size, size);
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
    __m256d negativeone = _mm256_set1_pd(-1.0);
    double *matrix = mat->data;
    double *resultmatrix = result->data;
    int size = mat->cols * mat->rows;

    #pragma omp parallel for
    for (int i = 0; i < size / 16 * 16; i+=16) {
        _mm256_storeu_pd(resultmatrix + i, _mm256_mul_pd(_mm256_loadu_pd(matrix + i), negativeone));
		_mm256_storeu_pd(resultmatrix + i + 4, _mm256_mul_pd(_mm256_loadu_pd(matrix + i + 4), negativeone));
		_mm256_storeu_pd(resultmatrix + i + 8, _mm256_mul_pd(_mm256_loadu_pd(matrix + i + 8), negativeone));
		_mm256_storeu_pd(resultmatrix + i + 12, _mm256_mul_pd(_mm256_loadu_pd(matrix + i + 12), negativeone));
    }
    #pragma omp parallel for
    for (int i = size / 16 * 16; i < size / 4 * 4; i+=4) {
        _mm256_storeu_pd(resultmatrix + i, _mm256_mul_pd(_mm256_loadu_pd(matrix + i), negativeone));
    }
    #pragma omp parallel for
    for (int i = size / 4 * 4; i < size; i++) {
        resultmatrix[i] = matrix[i]*-1;
    }
    return 0;
}

/*
 * Store the result of taking the absolute value element-wise to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int abs_matrix(matrix *result, matrix *mat) {
    /* TODO: YOUR CODE HERE */
    __m256d negativeone = _mm256_set1_pd(-1.0);
    double *matrix = mat->data;
    double *resultmatrix = result->data;
    int size = mat->cols * mat->rows;

    #pragma omp parallel for
    for (int i = 0; i < size / 16 * 16; i+=16) {
        __m256d load1 = _mm256_loadu_pd(matrix + i);
        __m256d load2 = _mm256_loadu_pd(matrix + i + 4);
        __m256d load3 = _mm256_loadu_pd(matrix + i + 8);
        __m256d load4 = _mm256_loadu_pd(matrix + i + 12);
        _mm256_storeu_pd(resultmatrix + i, _mm256_max_pd(_mm256_mul_pd(load1, negativeone), load1));
		_mm256_storeu_pd(resultmatrix + i + 4, _mm256_max_pd(_mm256_mul_pd(load2, negativeone), load2));
		_mm256_storeu_pd(resultmatrix + i + 8, _mm256_max_pd(_mm256_mul_pd(load3, negativeone), load3));
		_mm256_storeu_pd(resultmatrix + i + 12, _mm256_max_pd(_mm256_mul_pd(load4, negativeone), load4));
    }
    #pragma omp parallel for
    for (int i = size / 16 * 16; i < size / 4 * 4; i+=4) {
        __m256d load1 = _mm256_loadu_pd(matrix + i);
        _mm256_storeu_pd(resultmatrix + i, _mm256_max_pd(_mm256_mul_pd(load1, negativeone), load1));
    }
    #pragma omp parallel for
    for (int i = size / 4 * 4; i < size; i++) {
        resultmatrix[i] = fabs(mat->data[i]);
    }
    return 0;
}

