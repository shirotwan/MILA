#pragma once
#include "core.h"

template <class T, __INT n>
void __PLU(const T* RESTRICT M, __INT* RESTRICT P, T* RESTRICT L, T* RESTRICT U, __INT swaps){

    swaps = 0;
    T A[n*n] = {0};

    for(__INT i = 0; i != n*n; i++){
        __INT buf = (i % (n + 1) != 0) ? 0.0 : 1.0;
        U[i] = static_cast<T>(buf);
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
                T tmp = A[buf1];
                A[buf1] = A[buf2];
                A[buf2] = tmp;
            }
            __INT tmpi = P[k]; P[k] = P[p]; P[p] = tmpi;

            for(__INT j = 0; j != k; j++){
                __INT buf1 = k*n+j;
                __INT buf2 = p*n+j;
                T tmp = L[buf1];
                L[buf1] = L[buf2];
                L[buf2] = tmp;
            }

            swaps++;
        }

        for(__INT i = k; i != n; i++){
            T sum = 0;
            for(__INT s = 0; s != k; s++) sum = sum + L[i*n+s] * U[s*n+k];
            __INT buf = i*n+k;
            L[buf] = A[buf] - sum; 
        }

        for(__INT j = k+1; j != n; j++){
            T sum = 0;
            for(__INT s = 0; s != k; s++) sum = sum + L[k*n+s] * U[s*n+j];
            __INT buf1 = (n+1)*k; __INT buf2 = k*n+j;
            U[buf2] = (A[buf2]-sum)/L[buf1];
        }
    }
}

template <class T, __INT m, __INT n>
void __householderQR(const T* RESTRICT A, T* RESTRICT R, T* RESTRICT Q) {
    // Copiar A en R
    for(__INT i = 0; i != m * n; i++) R[i] = A[i];
    // Inicializar Q como identidad (m x m)
    for(__INT i = 0; i != m * m; i++) Q[i] = (i % (m + 1) == 0) ? 1.0 : 0.0;

    for(__INT k = 0; k < n && k < m - 1; k++) {
        // Calcular norma de la columna k desde fila k hasta m-1
        T norm_x = 0;
        for(__INT i = k; i != m; i++) {
            T val = R[i * n + k];
            norm_x = norm_x + val * val;
        }
        norm_x = sqrt(norm_x);

        // Construir vector v (en 1D)
        T v[m] = {0}; // tamaño fijo para simplicidad
        //for(__INT i = 0; i != m; i++) v[i] = 0.0;
        for(__INT i = k; i != m; i++) v[i] = R[i * n + k];
        v[k] += (R[k * n + k] >= 0) ? norm_x : -norm_x;

        // Normalizar v
        T norm_v = 0;
        for(__INT i = k; i != m; i++) norm_v = norm_v + v[i] * v[i];
        norm_v = sqrt(norm_v);
        for(__INT i = k; i != m; i++) v[i] = v[i]/norm_v;

        // Aplicar H a R → R = R - 2*v*(v^T*R)
        for(__INT j = k; j != n; j++) {
            T dot = 0;
            for(__INT i = k; i != m; i++) dot = dot + v[i] * R[i * n + j];
            for(__INT i = k; i != m; i++) R[i * n + j] = R[i * n + j] - 2.0 * v[i] * dot;
        }

        // Aplicar H a Q → Q = Q - 2*Q*v*v^T
        for(__INT j = 0; j != m; j++) {
            T dot = 0;
            for(__INT i = k; i != m; i++) dot = dot + Q[j * m + i] * v[i];
            for(__INT i = k; i != m; i++) Q[j * m + i] = Q[j * m + i] - 2.0 * dot * v[i];
        }
    }
}

template <class T, __INT m, __INT n>
void __pqr_householder(const T* RESTRICT A_in, T* RESTRICT Q, T* RESTRICT R, T* RESTRICT Pmat) {
    T v[m];            // vector reflector (uso local)
    T colNorms[n];
    __INT perm[n];

    // Copia entrada -> R (trabajamos sobre R)
    for (__INT i = 0; i != m * n; ++i) R[i] = A_in[i];

    // Inicializar Q = I_m
    for (__INT i = 0; i != m * m; ++i) Q[i] = 0.0;
    for (__INT i = 0; i != m; ++i) Q[i * m + i] = 1.0;

    // perm inicial
    for (__INT j = 0; j != n; ++j) perm[j] = j;

    // normas iniciales de columnas
    for (__INT j = 0; j != n; ++j) colNorms[j] = __col_norm_from<T,m,n>(R, j, 0);

    const __INT r = (m < n) ? m : n;

    for (__INT k = 0; k != r; ++k) {
        // --- Pivot: máxima norma (column pivoting) ---
        __INT pivot = k;
        T maxn = colNorms[k];
        for (__INT j = k + 1; j != n; ++j) {
            if (colNorms[j] > maxn) { maxn = colNorms[j]; pivot = j; }
        }
        if (pivot != k) {
            // __INTercambiar columnas k y pivot en R
            for (__INT i = 0; i != m; ++i) {
                T tmp = R[i * n + k];
                R[i * n + k] = R[i * n + pivot];
                R[i * n + pivot] = tmp;
            }
            T tmpd = colNorms[k]; colNorms[k] = colNorms[pivot]; colNorms[pivot] = tmpd;
            __INT tmpi = perm[k]; perm[k] = perm[pivot]; perm[pivot] = tmpi;
        }

        // --- Householder sobre la columna k (filas k..m-1) ---
        T normx = 0;
        for (__INT i = k; i != m; ++i) {
            T val = R[i * n + k];
            normx = normx + val * val;
        }
        normx = sqrt(normx);

        T alpha = (R[k * n + k] >= 0.0) ? -normx : normx;

        // construir v: v[k] = R[k,k] - alpha; v[i]=R[i,k] para i>k
        T vnorm2 = 0;
        for (__INT i = 0; i != k; ++i) v[i] = 0.0;
        v[k] = R[k * n + k] - alpha;
        vnorm2 = vnorm2 + v[k] * v[k];
        for (__INT i = k + 1; i != m; ++i) { v[i] = R[i * n + k]; vnorm2 = vnorm2 + v[i] * v[i]; }

        T beta = vnorm2; // v^T v

        // R := H * R  where H = I - 2 v v^T / beta
        for (__INT j = k; j != n; ++j) {
            T dot = 0;
            for (__INT i = k; i != m; ++i) dot = dot + v[i] * R[i * n + j];
            T coef = 2.0 * dot / beta;
            for (__INT i = k; i != m; ++i) R[i * n + j] = R[i * n + j] - coef * v[i];
        }
        // forzar forma superior en columna k
        R[k * n + k] = alpha;
        for (__INT i = k + 1; i != m; ++i) R[i * n + k] = 0.0;

        // Q := Q * H  (IMPORTANTE: post-multiplicar)
        // Primero u = Q * v  (u size m)
        T u[m];
        for (__INT i = 0; i != m; ++i) {
            T s = 0;
            // v[t]==0 for t<k, así que sumamos desde k
            for (__INT t = k; t < m; ++t) s = s + Q[i * m + t] * v[t];
            u[i] = s;
        }
        // Q -= 2 * u * v^T / beta
        for (__INT j = k; j != m; ++j) {
            T vj = v[j];
            T factor = 2.0 * vj / beta;
            for (__INT i = 0; i != m; ++i) Q[i * m + j] = Q[i * m + j] - u[i] * factor;
        }

        // actualizar normas columnas j>k (recalculo desde fila k+1)
        for (__INT j = k + 1; j != n; ++j) {
            T s = 0;
            for (__INT i = k + 1; i != m; ++i) {
                T vv = R[i * n + j];
                s = s + vv * vv;
            }
            colNorms[j] = sqrt(s);
        }
        colNorms[k] = 0.0;
    }

    // Construir matriz de permutación Pmat (n x n) tal que A * Pmat = columnas permutadas
    for (__INT i = 0; i != n * n; ++i) Pmat[i] = 0.0;
    for (__INT j = 0; j != n; ++j) {
        __INT orig = perm[j];
        Pmat[orig * n + j] = 1.0;
    }
}

template <class T, __INT n>
void __hess_householder(const T* RESTRICT A, T* RESTRICT Q, T* RESTRICT H){
    T v[n], u[n];
    for (__INT i = 0; i != n * n; ++i) {Q[i] = 0.0; H[i] = A[i];}
    for (__INT i = 0; i != n; ++i) Q[i * n + i] = 1.0;

    for (__INT k = 0; k < n - 2; ++k) {
        // longitud del subvector
        __INT m = n - (k + 1);            // elementos desde k+1 hasta n-1

        // construir x = A[k+1 : n-1, k]
        for (__INT i = 0; i < m; ++i) v[i] = H[(k + 1 + i)*n + k]; //H[(k + 1 + i)*n + k]

        // norma de x
        T sigma = __vec_norm<T>(v,m);
        if (sigma == 0.0) continue;

        // v = x + sign(x0)*||x|| * e1
        T sign = (v[0] >= 0.0) ? 1.0 : -1.0;
        v[0] = v[0] + sign * sigma;

        // beta = v^T v
        T beta = 0.0;
        for (__INT i = 0; i != m; ++i) beta = beta + v[i] * v[i];
        if (beta == 0.0) continue;

        // ----- Apply H on the left: A = H * A  for rows k+1..n-1 -----
        // For each column j=k .. n-1 compute dot = v^T * A[k+1.., j]
        for (__INT j = k; j != n; ++j) {
            T dot = 0.0;
            for (__INT i = 0; i != m; ++i) dot = dot + v[i] * H[(k + 1 + i)*n + j];
            T coef = 2.0 * dot / beta;
            for (__INT i = 0; i != m; ++i)
                H[(k + 1 + i)*n + j] = H[(k + 1 + i)*n + j] - coef * v[i];
        }

        // ----- Apply H on the right: A = A * H  for cols k+1..n-1 -----
        // For each row i=0..n-1 compute dot = A[i, k+1..] * v
        for (__INT i = 0; i != n; ++i) {
            T dot = 0.0;
            for (__INT j = 0; j != m; ++j) dot = dot + H[n*i + k + 1 + j] * v[j]; //H[n*i + k + 1 + j]
            T coef = 2.0 * dot / beta;
            for (__INT j = 0; j != m; ++j)
                H[n*i + k + 1 + j] = H[n*i + k + 1 + j] - coef * v[j];
        }

        // For exact Hessenberg enforce zeros under subdiagonal in column k
        // (numerical rounding may leave tiny values)
        for (__INT i = k + 2; i != n; ++i) H[n*i + k] = 0.0;

        // ----- Accumulate Q := Q * H  (post-multiply) -----
        // Compute u = Q * v_sub  where v_sub corresponds to positions k+1..n-1
        for (__INT i = 0; i != n; ++i) {
            T s = 0.0;
            for (__INT j = 0; j < m; ++j) s = s + Q[n*i + k + 1 + j] * v[j];
            u[i] = s;
        }
        // Q[:, k+1..] -= 2 * u * v^T / beta
        for (__INT j = 0; j != m; ++j) {
            T vj = v[j];
            T factor = 2.0 * vj / beta;
            for (__INT i = 0; i != n; ++i) {
                Q[n*i + k + 1 + j] = Q[n*i + k + 1 + j] - u[i] * factor;
            }
        }
    }
}