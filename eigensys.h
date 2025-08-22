#pragma once

#include "core.h"
#include "mcomplex.h"
#include "tol.h"

// More descs for complex numbers

template <class T, __INT n>
void __cplx_PLU(const CPLX<T>* RESTRICT M, __INT* RESTRICT P, CPLX<T>* RESTRICT L, CPLX<T>* RESTRICT U, __INT swaps){

    swaps = 0;
    CPLX<T> A[n*n] = {0};

    for(__INT i = 0; i != n*n; i++){
        __INT buf = (i % (n + 1) != 0) ? 0.0 : 1.0;
        U[i] = CPLX<T>(buf);
        A[i] = M[i];
    }

    for(__INT i = 0; i != n; i++){
        P[i] = i;
    }

    for(__INT k = 0; k != n; k++){
        __INT p = k;
        T maxv = abs( A[(n+1)*k] );

        for(__INT i = k+1; i != n; i++){
            T v = abs( A[i*n+k] );
            if(v > maxv){
                maxv = v;
                p = i;
            }
        }

        if(p != k){
            for(__INT j = 0; j != n; j++){
                __INT buf1 = k*n+j;
                __INT buf2 = p*n+j;
                CPLX<T> tmp = A[buf1];
                A[buf1] = A[buf2];
                A[buf2] = tmp;
            }
            __INT tmpi = P[k]; P[k] = P[p]; P[p] = tmpi;

            for(__INT j = 0; j != k; j++){
                __INT buf1 = k*n+j;
                __INT buf2 = p*n+j;
                CPLX<T> tmp = L[buf1];
                L[buf1] = L[buf2];
                L[buf2] = tmp;
            }

            swaps++;
        }

        for(__INT i = k; i != n; i++){
            CPLX<T> sum = 0;
            for(__INT s = 0; s != k; s++) sum = sum + L[i*n+s] * U[s*n+k];
            __INT buf = i*n+k;
            L[buf] = A[buf] - sum; 
        }

        for(__INT j = k+1; j != n; j++){
            CPLX<T> sum = 0;
            for(__INT s = 0; s != k; s++) sum = sum + L[k*n+s] * U[s*n+j];
            __INT buf1 = (n+1)*k; __INT buf2 = k*n+j;
            U[buf2] = (A[buf2]-sum)/L[buf1];
        }
    }
}

template <class T, __INT m, __INT n>
void __cplx_householderQR(const CPLX<T>* RESTRICT A, CPLX<T>* RESTRICT Q, CPLX<T>* RESTRICT R){
    // Copiar A en R
    for(__INT i = 0; i != m * n; i++) R[i] = A[i];
    // Inicializar Q como identidad (m x m)
    for(__INT i = 0; i != m * m; i++) Q[i] = (i % (m + 1) != 0) ? 0.0 : 1.0;

    for(__INT k = 0; k < n && k < m - 1; k++) {
        // Calcular norma de la columna k desde fila k hasta m-1
        T norm_x = 0;
        for(__INT i = k; i != m; i++) {
            CPLX<T> val = R[i * n + k];
            norm_x = norm_x + abs2(val);
        }
        norm_x = sqrt(norm_x);

        // Construir vector v (en 1D)
        CPLX<T> v[m] = {0}; // tamaño fijo para simplicidad
        //for(__INT i = 0; i != m; i++) v[i] = 0.0;
        for(__INT i = k; i != m; i++) v[i] = R[i * n + k];
        v[k] = v[k] + ((R[k * n + k] != 0) ? (norm_x * R[k * n + k] / abs(R[k * n + k])) : norm_x ); 
        // v[k] = v[k] + ((R[k * n + k] == 0) ? norm_x : (norm_x * R[k * n + k] / abs(R[k * n + k])) );

        // Normalizar v
        T norm_v = 0;
        for(__INT i = k; i != m; i++) norm_v = norm_v + abs2(v[i]);
        norm_v = sqrt(norm_v);
        for(__INT i = k; i != m; i++) v[i] = v[i]/norm_v;

        // Aplicar H a R → R = R - 2*v*(v^T*R)
        for(__INT j = k; j != n; j++) {
            CPLX<T> dot = 0;
            for(__INT i = k; i != m; i++) dot = dot + conj(v[i]) * R[i * n + j];
            for(__INT i = k; i != m; i++){
                __INT buf = i * n + j;
                R[buf] = R[buf] - 2.0 * v[i] * dot;
            } 
        }

        // Aplicar H a Q → Q = Q - 2*Q*v*v^T
        for(__INT j = 0; j != m; j++) {
            CPLX<T> dot = 0;
            for(__INT i = k; i != m; i++) dot = dot + Q[j * m + i] * v[i];
            for(__INT i = k; i != m; i++){
                __INT buf = j * m + i;
                Q[buf] = Q[buf] - 2.0 * dot * conj(v[i]);
            }
        }
    }
}

// Eigen vl and vv

template <class T, __INT n>
void __eigval(const CPLX<T>* RESTRICT A, CPLX<T>* RESTRICT EIGVL){
    CPLX<T> QMAT[n * n] = {0}; CPLX<T> RMAT[n * n] = {0}; CPLX<T> AT[n * n] = {0};
    for(__INT i = 0; i != n*n; i++) AT[i] = A[i];

    for(long int i = 0; i != __GLOBAL_MAX_ITERATIONS__; i++){
        __cplx_householderQR<T,n,n>(AT,QMAT,RMAT);
        __m_prod<CPLX<T>,n,n,n>(RMAT,QMAT,AT);
    }

    for(__INT i = 0; i != n; i++){
        if( i < n-1 && abs(AT[(i+1)*n+i])  > __GLOBAL_TOLERANCE__ ){
            __INT buf1 = (n+1)*i;
            __INT buf2 = buf1 + n + 1;
            CPLX<T> a = AT[buf1];
            CPLX<T> b = AT[buf1 + 1];
            CPLX<T> c = AT[buf2 - 1];
            CPLX<T> d = AT[buf2];

            CPLX<T> tr = a + d;
            CPLX<T> det = a * d - b * c;
            CPLX<T> disc = sqrt(tr*tr - 4.0 * det);

            EIGVL[i] = 0.5 * (tr + disc);
            EIGVL[i+1] = 0.5 * (tr - disc);
            i++;
        } else {
            EIGVL[i] = AT[(n+1)*i];
        }
    }
}

template <class T, __INT n>
void __eig(const CPLX<T>* RESTRICT A, CPLX<T>* RESTRICT EIGVL, CPLX<T>* RESTRICT EIGVV){
    CPLX<T> QMAT[n * n] = {0}; CPLX<T> RMAT[n * n] = {0}; CPLX<T> AT[n * n] = {0};
    for(__INT i = 0; i != n*n; i++) AT[i] = A[i];

    for(__LINT i = 0; i != __GLOBAL_MAX_ITERATIONS__; i++){
        __cplx_householderQR<T,n,n>(AT,QMAT,RMAT);
        __m_prod<CPLX<T>,n,n,n>(RMAT,QMAT,AT);
    }

    for(__INT i = 0; i != n; i++){
        if( i < n-1 && abs(AT[(i+1)*n+i])  > __GLOBAL_TOLERANCE__ ){
            __INT buf1 = (n+1)*i;
            __INT buf2 = buf1 + n + 1;
            CPLX<T> a = AT[buf1];
            CPLX<T> b = AT[buf1 + 1];
            CPLX<T> c = AT[buf2 - 1];
            CPLX<T> d = AT[buf2];

            CPLX<T> tr = a + d;
            CPLX<T> det = a * d - b * c;
            CPLX<T> disc = sqrt(tr*tr - 4.0 * det);

            EIGVL[i] = 0.5 * (tr + disc);
            EIGVL[i+1] = 0.5 * (tr - disc);
            i++;
        } else {
            EIGVL[i] = AT[(n+1)*i];
        }
    }

    for(__INT i = 0; i != n; i++){
        for(__INT k = 0; k != n*n; k++) QMAT[k] = A[k] - ( (k % (n+1) != 0)? 0.0 : EIGVL[i] );

        CPLX<T> LV[n*n] = {0}; CPLX<T> UV[n*n] = {0}; __INT PV[n] = {0}; __INT sw = 0;

        __cplx_PLU<T,n>(QMAT,PV,LV,UV,sw);

        CPLX<T> VV[n] = {0}; VV[n-1] = 1;

        for (__SINT j = n-2; j >= 0; j--) {
            CPLX<T> sum = 0;
            for (__INT k = j+1; k != n; k++) {
                sum = sum + UV[j*n + k] * VV[k];
            }

            if (abs(UV[j*n + j]) < __GLOBAL_TOLERANCE__) {
                // variable libre
                VV[j] = 1.0;
            } else {
                VV[j] = -1.0*sum / UV[j*n + j];
            }
        }

        T norm_c = 0;
        for(__INT j = 0; j != n; j++) norm_c = norm_c + abs2(VV[j]);
        norm_c = sqrt(norm_c);

        for(__INT j = 0; j != n; j++) VV[j] = VV[j]/norm_c;

        for(__INT j = 0; j != n; j++) EIGVV[j*n+i] = VV[j];
    }

}
