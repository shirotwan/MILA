#pragma once
#include "eigensys.h"
#include "msolver.h"

template <class T, __INT n>
void __spectralfun(CPLX<T> f(CPLX<T>), const CPLX<T>* RESTRICT M, CPLX<T>* RESTRICT S){
    CPLX<T> ZV[n*n] = {0}; CPLX<T> DV[n] = {0}; 
    __INT PM[n] = {0}; CPLX<T> LM[n*n] = {0}; CPLX<T> UM[n*n] = {0}; __INT sw = 0;
    CPLX<T> IZV[n*n] = {0};

    __eig<T,n>(M,DV,ZV);
    __cplx_PLU<T,n>(ZV,PM,LM,UM,sw);
    __inv_classic<CPLX<T>,n>(PM,LM,UM,IZV);

    // Reciclar L para optimizar
    for(__INT i = 0; i != n*n; i++) LM[i] = ZV[i] * f(DV[i%n]);

    __m_prod<CPLX<T>,n,n,n>(LM,IZV,S);
}