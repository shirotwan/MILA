#pragma once
#include "varsetup.h"
#include "mathutils.h"

// Helpers
#ifndef EMB_SYS
template <class T>
inline void m_print(const T* RESTRICT mat, const __INT& rows, const __INT& cols){
    for(__INT i = 0; i != rows; i++){
    for(__INT j = 0; j != cols; j++){
        std::cout << mat[cols * i + j] << " ";
    }
        std::cout << "\n";
    }
}
#endif
// Functions - Operators of vectors

template <class T, __INT m, __INT n>
inline T __col_norm_from(const T* RESTRICT A, __INT col, __INT start) {
    T s = 0;
    for (int i = start; i != m; ++i) {
        T v = A[i * n + col];
        s = s + v * v;
    }
    return sqrt(s);
}

template <class T>
inline T __vec_norm(const T* RESTRICT v, __INT len) {
    T s = 0.0;
    for (__INT i = 0; i != len; ++i) s = s + v[i] * v[i];
    return sqrt(s);
}

// Functions - Operators of matrices

template <class T, __INT rows, __INT cols>
inline void __m_sum(const T* RESTRICT lhs, const T* RESTRICT rhs, T* RESTRICT res){
    for(__INT i = 0; i != rows; i++){
    for(__INT j = 0; j != cols; j++){
        res[cols * i + j] = lhs[cols * i + j] + rhs[cols * i + j];
    }
    }
}

template <class T, __INT rows_lhs, __INT common_dim, __INT cols_rhs>
inline void __m_prod(const T* RESTRICT lhs, const T* RESTRICT rhs, T* RESTRICT res){
    for(__INT i = 0; i != rows_lhs; i++){
    for(__INT j = 0; j != cols_rhs; j++){
        T sum = 0;
        for(__INT k = 0; k != common_dim; k++){
            sum = sum + lhs[common_dim * i + k] * rhs[cols_rhs * k + j];
        }
        res[cols_rhs * i + j] = sum;
    }
    }
}

template <class T, __INT n>
inline void __diag(T* RESTRICT ARR, T* RESTRICT RES){
    for(__INT i = 0; i != n; i++) RES[i*n+i] = ARR[i];
}

template <class T, __INT n>
inline T __trace(const T* RESTRICT A){
    T res = 0;
    for(__INT i = 0; i != n; i++) res = res + RES[i*n+i];
    return res;
}