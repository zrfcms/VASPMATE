
#ifndef MSYM_EXPORT_H
#define MSYM_EXPORT_H

#ifdef MSYM_EXPORTS_BUILT_AS_STATIC
#  define MSYM_EXPORT
#  define MSYM_NO_EXPORT
#else
#  ifndef MSYM_EXPORT
#    ifdef msym_EXPORTS
        /* We are building this library */
#      define MSYM_EXPORT __attribute__((visibility("default")))
#    else
        /* We are using this library */
#      define MSYM_EXPORT __attribute__((visibility("default")))
#    endif
#  endif

#  ifndef MSYM_NO_EXPORT
#    define MSYM_NO_EXPORT __attribute__((visibility("hidden")))
#  endif
#endif

#ifndef MSYM_DEPRECATED
#  define MSYM_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef MSYM_DEPRECATED_EXPORT
#  define MSYM_DEPRECATED_EXPORT MSYM_EXPORT MSYM_DEPRECATED
#endif

#ifndef MSYM_DEPRECATED_NO_EXPORT
#  define MSYM_DEPRECATED_NO_EXPORT MSYM_NO_EXPORT MSYM_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef MSYM_NO_DEPRECATED
#    define MSYM_NO_DEPRECATED
#  endif
#endif

#endif
