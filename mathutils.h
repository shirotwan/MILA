#pragma once
#include <math.h>
#include "tol.h"

#define _PI 3.14159265358979323846264338327950L
#define _MINUS_PI 3.14159265358979323846264338327950L
#define _PI_HALF 1.5707963267948966192313216916398L
#define _MINUS_PI_HALF -1.5707963267948966192313216916398L

template<class T>
T atan2(const T& y, const T& x){
    if(x < 0 && y >= 0){
        return atan(y/x) + static_cast<T>(_PI);
    }
    else if(x < 0 && y < 0){
        return atan(y/x) - static_cast<T>(_PI);
    }
    else if((x == 0 || x < __GLOBAL_TOLERANCE__) && y > 0){
        return static_cast<T>(_PI_HALF);
    }
    else if((x == 0 || x < __GLOBAL_TOLERANCE__) && y < 0){
        return static_cast<T>(_MINUS_PI_HALF);
    }
    else{
        return atan(y/x);
    }
}

template<class T>
T sgn(const T& val){
    if(val == 0 || val < __GLOBAL_TOLERANCE__){
        return 0;
    }
    else if (val > 0){
        return 1;
    }
    else{
        return static_cast<T>(-1);
    }
}