#if defined TOOLCHAIN_MSVC

    #define LIB_PUBLIC __declspec(dllexport)
    #define LIB_PUBLIC_IMPL __declspec(dllexport)
    #define LIB_LOCAL

#elif defined TOOLCHAIN_GCC

    #define LIB_PUBLIC __attribute__ ((visibility ("default")))
    #define LIB_PUBLIC_IMPL
    #define LIB_LOCAL  __attribute__ ((visibility ("hidden")))

#elif defined TOOLCHAIN_CLANG

    #define LIB_PUBLIC __attribute__ ((visibility ("default")))
    #define LIB_PUBLIC_IMPL
    #define LIB_LOCAL  __attribute__ ((visibility ("hidden")))

#else

    #warning "No toolchain defined"

#endif
