/* -------------------------------------------------------------

This file is a component of SDPA
Copyright (C) 2004 SDPA Project

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA

------------------------------------------------------------- */
/*--------------------------------------------------
  rsdpa_tool.h
  $Id: rsdpa_tool.h,v 1.2 2004/09/01 06:34:12 makoto Exp $
--------------------------------------------------*/

#ifndef __sdpa_tool_h__
#define __sdpa_tool_h__

#include <sdpa_right.h>

#include <iostream>
#include <sys/time.h>
#include <string>

#include <gmpxx.h>

namespace sdpa {

#if 1
#define rMessage(message) \
cout << message << " :: line " << __LINE__ \
  << " in " << __FILE__ << endl
#else
#define rMessage(message)
#endif

#define rError(message) \
cout << message << " :: line " << __LINE__ \
  << " in " << __FILE__ << endl; \
exit(false)

#if 0
#define rNewCheck() rMessage("new invoked");
#else
#define rNewCheck() ;
#endif

#define REVERSE_PRIMAL_DUAL 1


// These are constant. Do NOT change
extern int IZERO   ; // =  0;
extern int IONE    ; // =  1;
extern int IMONE   ; // = -1;
extern mpf_class MZERO; // =  0.0;
extern mpf_class MONE ; // =  1.0;
extern mpf_class MMONE; // = -1.0;

class Time
{
public:
  static double rGetUseTime();
  static void rSetTimeVal(struct timeval & targetVal);
  static double rGetRealTime(const struct timeval & start,
			     const struct timeval & end);
};

#if 1 // count time with process time
#define TimeStart(START__) \
   static double START__; START__ = Time::rGetUseTime()
#define TimeEnd(END__) \
   static double END__;   END__ = Time::rGetUseTime()
#define TimeCal(START__,END__) (END__ - START__)
#else // count time with real time
#define TimeStart(START__) \
   static struct timeval START__; Time::rSetTimeVal(START__)
#define TimeEnd(END__) \
   static struct timeval END__; Time::rSetTimeVal(END__)
#define TimeCal(START__,END__) Time::rGetRealTime(START__,END__)
#endif

}

#endif // __sdpa_tool_h__
