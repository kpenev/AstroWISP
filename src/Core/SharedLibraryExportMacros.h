#if defined WINDOWS_VISIBILITY

    #define LIB_PUBLIC __declspec(dllexport)
    #define LIB_PUBLIC_IMPL __declspec(dllexport)
    #define LIB_LOCAL

#elif defined NIX_VISIBILITY

    #define LIB_PUBLIC __attribute__ ((visibility ("default")))
    #define LIB_PUBLIC_IMPL
    #define LIB_LOCAL  __attribute__ ((visibility ("hidden")))

#else

//    #warning "No toolchain defined"

#endif
