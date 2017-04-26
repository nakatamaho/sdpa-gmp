/* -------------------------------------------------------------

This file is a component of SDPA-C
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

/* $Id: rsdpa_chordal.h,v 1.1.1.1 2004/09/07 07:43:59 makoto Exp $ */
/*-----------------------------------------
  rsdpa_chordal.h
-----------------------------------------*/

#ifndef __sdpa_chordal_h__
#define __sdpa_chordal_h__

#include <sdpa_dataset.h>

#ifdef __cplusplus
extern "C" {
#endif

#include "InpMtx.h"
#include "IVL.h"
#include "Graph.h"
#include "ETree.h"
#include "misc.h"
#include "IV.h"
#include "SymbFac.h"

#define TIMES_PER_SECOND CLK_TCK

#ifdef __cplusplus
}
#endif

#define unsgnshrt unsigned int
#define MAXCHAR 60
#define UNITCLIQUE 30
#define UNITSON 20
#define FORDEBUG 1

//#define min(a,b) (a < b ? a : b)


namespace sdpa {

class Chordal {
public:

  // condition of sparse computation
  // m_threshold < mDim, 
  // b_threshold < nBlock, 
  // aggregate_threshold >= aggrigated sparsity ratio
  // extend_threshold    >= extended sparsity ratio
  int m_threshold;
  int b_threshold;
  double aggregate_threshold;
  double extend_threshold;

  int Method[5];
  int   best;                     /* indicates the best ordering method */
/* indicates the used ordering method */
  /* -1: dense computation */
  /* 0: METIS 4.0.1 - nested dissection <--- not support*/
  /* 1: Spooles 2.2 - mininum degree */
  /* 2: Spooles 2.2 - generalized nested dissection */
  /* 3: Spooles 2.2 - multisection */
  /* 4: Spooles 2.2 - better of 2 and 3 */

  ETree *etree;
  IVL *adjIVL;
  IVL *symbfacIVL_MMD, *symbfacIVL_NDMS, *symbfacIVL_ND, *symbfacIVL_MS;
  IV *newToOldIV_MMD, *newToOldIV_NDMS, *newToOldIV_ND,*newToOldIV_MS;
  Graph *graph;

  Chordal(void);
  ~Chordal();
  void initialize();
  void terminate();

  // marge array1 to array2
  void margeArray(int na1, int* array1, int na2, int* array2);
  void makeGraph(InputData& inputData, int m);
  int countNonZero(int m, IVL *symbfacIVL);
  int Spooles_MMD(int m);
  int Spooles_NDMS(int m);
  int Spooles_ND(int m);
  int Spooles_MS(int m);
  void ordering_bMat(int m, int nBlock,
                     InputData& inputData, 
                     FILE* fpOut);
  int   Best_Ordering(int *Method);
};

}

#endif // __sdpa_chordal_h__
