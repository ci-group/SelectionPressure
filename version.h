/*
**  /Users/evert/projects/workspace/SelectionPressure/version.h -- Version Information for SelectionPressure (syntax: C/C++)
**  [automatically generated and maintained by GNU shtool]
*/

#ifdef _VERSION_H_AS_HEADER_

#ifndef __VERSION_H_
#define __VERSION_H_

#define VERSION 0x001207

typedef struct {
    const int   v_hex;
    const char *v_short;
    const char *v_long;
    const char *v_tex;
    const char *v_gnu;
    const char *v_web;
    const char *v_sccs;
    const char *v_rcs;
} version_t;

extern version_t version;

#endif /* __VERSION_H_ */

#else /* _VERSION_H_AS_HEADER_ */

#define _VERSION_H_AS_HEADER_
#include "version.h"
#undef  _VERSION_H_AS_HEADER_

version_t version = {
    0x001207,
    "0.1.7",
    "0.1.7 (01-Oct-2015)",
    "This is SelectionPressure, Version 0.1.7 (01-Oct-2015)",
    "SelectionPressure 0.1.7 (01-Oct-2015)",
    "SelectionPressure/0.1.7",
    "@(#)SelectionPressure 0.1.7 (01-Oct-2015)",
    "$Id: version.h 164 2015-10-12 08:12:34Z ehaasdi $"
};

#endif /* _VERSION_H_AS_HEADER_ */

