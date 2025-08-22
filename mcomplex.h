#pragma once
#include "mathutils.h"

#ifndef __USING_CPLX_VAL__
#define __USING_CPLX_VAL__
#endif

template <class T> class CPLX{
    T dcp[2] = {0,0};
    public:
    CPLX(){}
    CPLX(const T& real, const T& imag){dcp[0] = real; dcp[1] = imag;}
    CPLX(const T& real){dcp[0] = real;}
    inline T* data() noexcept {return dcp;} 
    inline T re() noexcept {return dcp[0];}
    inline T im() noexcept {return dcp[1];}

    friend T abs(CPLX<T> num) {return sqrt(num.dcp[0]*num.dcp[0] + num.dcp[1]*num.dcp[1]);}
    friend T abs2(CPLX<T> num) noexcept {return num.dcp[0]*num.dcp[0] + num.dcp[1]*num.dcp[1];}
    friend CPLX<T> conj(CPLX<T> num) noexcept {return CPLX<T>(num.dcp[0],-num.dcp[1]);}
    friend CPLX<T> sqrt(CPLX<T> num) {
        T mdl = abs(num), ct = num.dcp[0]/mdl, st = num.dcp[1]/mdl;
        return CPLX<T>( sqrt(0.5*mdl*(1 + ct)) , sgn(st) * sqrt(0.5*mdl*(1 - ct)) );
    }

    #ifndef EMB_SYS
    friend std::ostream& operator << (std::ostream& os, CPLX<T> num){ return os << "(" << num.dcp[0] << ";" << num.dcp[1] << ")"; }
    #endif
    
    friend CPLX<T> operator + (CPLX<T> num1, CPLX<T> num2) noexcept {return CPLX<T>(num1.dcp[0] + num2.dcp[0], num1.dcp[1] + num2.dcp[1]);}
    friend CPLX<T> operator - (CPLX<T> num1, CPLX<T> num2) noexcept {return CPLX<T>(num1.dcp[0] - num2.dcp[0], num1.dcp[1] - num2.dcp[1]);}
    friend CPLX<T> operator * (CPLX<T> num1, CPLX<T> num2) noexcept {return CPLX<T>(num1.dcp[0] * num2.dcp[0] - num1.dcp[1] * num2.dcp[1], num1.dcp[0] * num2.dcp[1] + num1.dcp[1] * num2.dcp[0]);}
    friend CPLX<T> operator / (CPLX<T> num1, CPLX<T> num2) 
    {
        return CPLX<T>( ( num1.dcp[0] * num2.dcp[0] + num1.dcp[1] * num2.dcp[1] )/( num2.dcp[0] *  num2.dcp[0] + num2.dcp[1] *  num2.dcp[1] ) , ( num1.dcp[1] * num2.dcp[0] - num1.dcp[0] * num2.dcp[1] )/( num2.dcp[0] *  num2.dcp[0] + num2.dcp[1] *  num2.dcp[1]) );
    }

    friend __INT operator == (CPLX<T> num1, CPLX<T> num2) noexcept { return (num1.dcp[0] == num2.dcp[0] && num1.dcp[1] == num2.dcp[1])?1:0; }
    friend __INT operator != (CPLX<T> num1, CPLX<T> num2) noexcept { return (num1.dcp[0] != num2.dcp[0] && num1.dcp[1] != num2.dcp[1])?1:0; }
};