/*-----------------------------------------------------
 * error.cc
 * Implementation of error handling functions
 *
 * Author: P. Foggia
 * $Id: error.cc,v 1.3 2013/08/27 15:51:38 zfeng Exp $
 ----------------------------------------------------*/

/*----------------------------------------------------
 * REVISION HISTORY
 *   $Log: error.cc,v $
 *   Revision 1.3  2013/08/27 15:51:38  zfeng
 *   2013-08-27
 *
 *   Revision 1.2  2011/09/27 18:50:02  zfeng
 *   changed to exception handling
 *
 *   Revision 1.1.1.1  2010/11/25 02:01:42  zfeng
 *   import vflib
 *
 *   Revision 1.2  1998/12/12 12:18:07  foggia
 *   Now supports full printf syntax
 *
 *   Revision 1.1  1998/09/16 17:35:14  foggia
 *   Initial revision
 *
 ---------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string>
#include <stdexcept>

#include "error.h"


/*------------------------------------------
 * void error(msg, ...)
 * Prints an error message and exits 
 * the program.
 * The syntax is the same of printf, 
 * except that a trailing \n is automatically
 * appended.
 -----------------------------------------*/
void error(const char *msg, ...)
  { va_list ap;
// change to exception handling
    va_start(ap, msg);
    char buffer[1000];
    vsprintf(buffer, msg, ap);
    std::string cs = "ERROR: ";
    cs += buffer;
    cs += "\n";
    throw std::out_of_range(cs);
/*
    fprintf(stderr, "ERROR: ");
    vfprintf(stderr, msg, ap);
    fprintf(stderr, "\n");
    va_end(ap);
    exit(1);
*/
  }

