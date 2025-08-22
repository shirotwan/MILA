#pragma once
#include "varsetup.h"

template <class T, __INT n>
void __inv_classic(const __INT* RESTRICT P, const T* RESTRICT L, const T* RESTRICT U, T* RESTRICT INV){
    T X[n] = {0}; T Y[n] = {0}; T B[n] = {0}; __INT swaps = 0;

    for(__INT j = 0; j != n; j++){
        for(__INT i = 0; i != n; i++) B[i] = T((P[i] == j)? 1.0 : 0.0);
        // Hacia adelante
        for(__INT i = 0; i != n; i++){
            T s = B[i];
            for(__INT k = 0; k != i; k++) s = s - L[i*n+k] * Y[k];
            Y[i] = s / L[i*n+i];
        }
        // Hacia atrás
        for(__SINT i = n-1; i >= 0; i--){
            T s = Y[i];
            for(__INT k = i+1; k < n; k++) s = s - U[i*n+k] * X[k];
            X[i] = s;
        }
        // Asignar
        for(__INT i = 0; i != n; i++) INV[i*n+j] = X[i];
    }
}

template <class T, __INT n>
T __det_classic(const T* RESTRICT L, const __INT& swaps){
    T res = static_cast<T>((swaps % 2 != 0)? 1.0 : -1.0);
    for(__INT i = 0; i != n; i++){
        res = res * L[i*n+i];
    }
    return res;
}