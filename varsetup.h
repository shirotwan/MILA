#pragma once

#include <stdio.h>

#if defined(_WIN32)
    #define PLATFORM_NAME "windows" // Windows
#elif defined(_WIN64)
    #define PLATFORM_NAME "windows" // Windows
#elif defined(__CYGWIN__) && !defined(_WIN32)
    #define PLATFORM_NAME "windows" // Windows (Cygwin POSIX under Microsoft Window)
#elif defined(__ANDROID__)
    #define PLATFORM_NAME "android" // Android (implies Linux, so it must come first)
#elif defined(__linux__)
    #define PLATFORM_NAME "linux" // Debian, Ubuntu, Gentoo, Fedora, openSUSE, RedHat, Centos and other
#elif defined(__unix__) || !defined(__APPLE__) && defined(__MACH__)
    #include <sys/param.h>
    #if defined(BSD)
        #define PLATFORM_NAME "bsd" // FreeBSD, NetBSD, OpenBSD, DragonFly BSD
    #endif
#elif defined(__hpux)
    #define PLATFORM_NAME "hp-ux" // HP-UX
#elif defined(_AIX)
    #define PLATFORM_NAME "aix" // IBM AIX
#elif defined(__APPLE__) && defined(__MACH__) // Apple OSX and iOS (Darwin)
    #include <TargetConditionals.h>
    #if TARGET_IPHONE_SIMULATOR == 1
        #define PLATFORM_NAME "ios" // Apple iOS
    #elif TARGET_OS_IPHONE == 1
        #define PLATFORM_NAME "ios" // Apple iOS
    #elif TARGET_OS_MAC == 1
        #define PLATFORM_NAME "osx" // Apple OSX
    #endif
#elif defined(__sun) && defined(__SVR4)
    #define PLATFORM_NAME "solaris" // Oracle Solaris, Open Indiana
#else
    #define PLATFORM_NAME "embebbed"
    #define EMB_SYS
#endif

#ifndef EMB_SYS
    #include <iostream>
    typedef size_t __INT;
    typedef ptrdiff_t __SINT;
    typedef unsigned long int __LINT;
#else
    typedef unsigned int __INT;
    typedef int __SINT;
    typedef unsigned long int __LINT;
#endif

/**
* @brief Find the restrict keyword from "any" compiler possible
* @note If there is a better way to modify it to reach more compilers, please do it!
*/

#if defined(__clang__) || defined(__GNUC__)
  #define RESTRICT __restrict__
#elif defined(_MSC_VER)
  #define RESTRICT __restrict
#else
  #define RESTRICT
#endif

/**
* @brief Find the best int keyword for loops or space
* @note You can still modify the real type with ""
*/

template <class T, __INT rows, __INT cols>
class MTX{
    private:
        T __mtxbuff[rows * cols] = {0};
    public:
    MTX(){}
    MTX(const T* RESTRICT data_pointer){
        const __INT max_sized = rows*cols;
        for(__INT i = 0; i != max_sized; i++){
            __mtxbuff[i] = data_pointer[i];
        }
    }
    inline T* dat(){
        return __mtxbuff;
    }
};