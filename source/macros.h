/*
 * This macro was added to avoid compiler warnings as the compiler is unable
 * to recognise that ExitPython will cause the program to abort. This fixes
 * warnings where the compiler thinks an uninitialized value will be used or
 * returned from a function.
 */

#define Exit(error_code) \
{ \
  Exit_python(error_code);\
  exit(error_code);\
}
