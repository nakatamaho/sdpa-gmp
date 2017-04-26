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

#define UseMETIS 0
#define PrintSparsity 0
#define OrderOnlyByMDO 1

#include <sdpa_chordal.h>

namespace sdpa {

Chordal::Chordal()
{
  initialize();
}

Chordal::~Chordal()
{
  //
}

void Chordal::initialize()
{
  // condition of sparse computation
  // m_threshold < mDim, 
  // b_threshold < nBlock, 
  // aggregate_threshold >= aggrigated sparsity ratio
  // extend_threshold    >= extended sparsity ratio
  m_threshold = 100;
  b_threshold = 5;
  aggregate_threshold = 0.25;
  extend_threshold = 0.4;
  
#if 0 // DENSE computation for debugging
  m_threshold = 10000000;
  b_threshold = 1000000; 
  aggregate_threshold = 0.0; 
  extend_threshold = 0.0; 
#endif
#if 0 // SPARSE computation for debugging
  m_threshold = 0;
  b_threshold = 0; 
  aggregate_threshold = 2.0; 
  extend_threshold = 2.0; 
#endif

/* indicates the used ordering method */
  /* 0: METIS 4.0.1 - nested dissection <--- not support */
  /* 1: Spooles 2.2 - mininum degree */
  /* 2: Spooles 2.2 - generalized nested dissection */
  /* 3: Spooles 2.2 - multisection */
  /* 4: Spooles 2.2 - better of 2 and 3 */
#if OrderOnlyByMDO
  Method[0] = 0;
  Method[1] = 1;
  Method[2] = 0;
  Method[3] = 0;
  Method[4] = 0;
#else
  Method[0] = 0;
  Method[1] = 1;
  Method[2] = 1;
  Method[3] = 1;
  Method[4] = 1;
#endif
  best = -1;
}

void Chordal::terminate()
{
  if (Method[0]) {
    rError("no support for METIS");
  }
  if (Method[1]>1)  {
    IV_free(newToOldIV_MMD);
    IVL_free(symbfacIVL_MMD);
  } 
  if (Method[2]>1)  {
    IV_free(newToOldIV_ND);
    IVL_free(symbfacIVL_ND);
  } 
  if (Method[3]>1)  {
    IV_free(newToOldIV_MS);
    IVL_free(symbfacIVL_MS);
  } 
  if (Method[4]>1)  {
    IV_free(newToOldIV_NDMS);
    IVL_free(symbfacIVL_NDMS);
  } 
}

// marge array1 to array2
void Chordal::margeArray(int na1, int* array1, int na2, int* array2)
{

  int ptr = na1 + na2 - 1;
  int ptr1 = na1-1;
  int ptr2 = na2-1;
  int idx1, idx2;
  
  while ((ptr1 >= 0) || (ptr2 >= 0)){

    if (ptr1 >= 0){
      idx1 = array1[ptr1];
    } else {
      idx1 = -1;
    }
    if (ptr2 >= 0 ){
      idx2 = array2[ptr2];
    } else {
      idx2 = -1;
    }
    if (idx1 > idx2){
      array2[ptr] = idx1;
      ptr1--;
    } else {
      array2[ptr] = idx2;
      ptr2--;
    }
    ptr--;

  }

  // error check
  if (ptr != -1){
    rMessage("Chordal::margeArray:: program bug");
  }

}

  // make aggrigate sparsity pattern
void Chordal::makeGraph(InputData& inputData, int m)
{

  int i,j,k,l;
  int SDP_nBlock = inputData.SDP_nBlock;
  int SOCP_nBlock = inputData.SOCP_nBlock;
  int LP_nBlock = inputData.LP_nBlock;

  int* counter;
  counter = new int[m];
  for (int i=0; i<m; i++){
    counter[i] = 0;
  }

  // count maximum mumber of index
  for (l = 0; l<SDP_nBlock; l++){
    int SDP_nConstraint = inputData.SDP_nConstraint[l];
    for (k=0; k<SDP_nConstraint; k++){
      i = inputData.SDP_constraint[l][k];
      counter[i] += SDP_nConstraint;
    }
  }
  for (l = 0; l<SOCP_nBlock; l++){
    int SOCP_nConstraint = inputData.SOCP_nConstraint[l];
    for (k=0; k<SOCP_nConstraint; k++){
      i = inputData.SOCP_constraint[l][k];
      counter[i] += SOCP_nConstraint;
    }
  }
  for (l = 0; l<LP_nBlock; l++){
    int LP_nConstraint = inputData.LP_nConstraint[l];
    for (k=0; k<LP_nConstraint; k++){
      i = inputData.LP_constraint[l][k];
      counter[i] += LP_nConstraint;
    }
  }

  // allocate temporaly workspace
  int** tmp;
  tmp = new int*[m];
  for (i=0; i<m; i++){
    tmp[i] = new int[counter[i]];
  }

  // merge index
  for (int i=0; i<m; i++){
    counter[i] = 0;
  }
  // marge index of for SDP
  for (l = 0; l<SDP_nBlock; l++){
    for (k=0; k<inputData.SDP_nConstraint[l]; k++){
      i = inputData.SDP_constraint[l][k];
      margeArray(inputData.SDP_nConstraint[l], inputData.SDP_constraint[l],
                 counter[i], tmp[i]);
      counter[i] += inputData.SDP_nConstraint[l];
    }
  }
  // marge index of for SOCP
  for (l = 0; l<SOCP_nBlock; l++){
    for (k=0; k<inputData.SOCP_nConstraint[l]; k++){
      i = inputData.SOCP_constraint[l][k];
      margeArray(inputData.SOCP_nConstraint[l], inputData.SOCP_constraint[l],
                 counter[i], tmp[i]);
      counter[i] += inputData.SOCP_nConstraint[l];
    }
  }
  // marge index of for LP
  for (l = 0; l<LP_nBlock; l++){
    for (k=0; k<inputData.LP_nConstraint[l]; k++){
      i = inputData.LP_constraint[l][k];
      margeArray(inputData.LP_nConstraint[l], inputData.LP_constraint[l],
                 counter[i], tmp[i]);
      counter[i] += inputData.LP_nConstraint[l];
    }
  }

  // construct adjacency list of SPOOLES  
  IVL_init1(adjIVL,IVL_CHUNKED,m);
  int isize, previous;
  int* ivec;
  ivec = new int[m];
  for (i=0; i<m; i++){
    isize = 0;
    previous = -1;
    for (j=0; j<counter[i]; j++){
      if (tmp[i][j] != previous){
        ivec[isize] = tmp[i][j];
        previous = ivec[isize];
        isize++;
      }
    }
    IVL_setList(adjIVL, i, isize, ivec);
  }

  // constract graph of SPOOLES
  Graph_init2(graph, 0, m, 0, IVL_tsize(adjIVL), m, IVL_tsize(adjIVL), 
              adjIVL, NULL, NULL);

  delete[] counter;
  for (int i=0; i<m; i++){
    delete[] tmp[i];
  }
  delete[] tmp;
  delete[] ivec;

}

int Chordal::countNonZero(int m, IVL *symbfacIVL)
{
  int nonzeros = 0;
  bool* bnode;

  // count non-zero element
  rNewCheck();
  bnode = new bool[m];
  if (bnode == NULL) {
    rError("Newton::initialize_sparse_bMat memory exhausted ");
  }
  for (int i=0; i<m; i++){
    bnode[i] = false;
  }

  int nClique = IVL_nlist(symbfacIVL);
  int psize;
  int* pivec;
  for (int l=nClique-1; l >= 0; l--){
    IVL_listAndSize(symbfacIVL,l,&psize,&pivec);
    for (int i=0; i<psize; i++){
      int ii = pivec[i];
      if (bnode[ii] == false){
        nonzeros += psize - i;
        bnode[ii] = true;
      }
    }
  }

  delete[] bnode;
  return nonzeros;
}

int Chordal::Spooles_MMD(int m)
{
  int seed = 0, msglvl = 0;
  FILE   *fp = NULL;

  //  rMessage("orderViaMMD:start");
  etree = orderViaMMD(graph, seed, msglvl, fp);
  //  rMessage("orderViaMMD:end");
  newToOldIV_MMD = ETree_newToOldVtxPerm(etree);
  symbfacIVL_MMD = SymbFac_initFromGraph(etree,graph);
  //  IVL_writeForHumanEye(symbfacIVL_MMD,stdout);

  int nonzeros = countNonZero(m, symbfacIVL_MMD);

  return nonzeros * 2 - m;
}

int Chordal::Spooles_NDMS(int m)
{
  int seed = 0, msglvl = 0;
  FILE   *fp = NULL;

  int maxdomainsize = m / 16 + 1; 
  int maxzeros = m / 10 + 1;
  int maxsize = 64;

  //  rMessage("orderViaMMD:start");
  etree = orderViaBestOfNDandMS(graph, maxdomainsize, maxzeros, maxsize, seed, msglvl, fp);
  //  rMessage("orderViaMMD:end");
  newToOldIV_NDMS = ETree_newToOldVtxPerm(etree);
  symbfacIVL_NDMS = SymbFac_initFromGraph(etree,graph);
  //  IVL_writeForHumanEye(symbfacIVL_NDMS,stdout);

  int nonzeros = countNonZero(m, symbfacIVL_NDMS);

  return nonzeros * 2 - m;
}

int Chordal::Spooles_ND(int m)
{
  int seed = 0, msglvl = 0;
  FILE   *fp = NULL;
  bool* bnode;

  int maxdomainsize = m / 16 + 1;

  //  rMessage("orderViaMMD:start");
  etree = orderViaND(graph, maxdomainsize, seed, msglvl, fp);
  //  rMessage("orderViaMMD:end");
  newToOldIV_ND = ETree_newToOldVtxPerm(etree);
  symbfacIVL_ND = SymbFac_initFromGraph(etree,graph);
  //  IVL_writeForHumanEye(symbfacIVL_ND,stdout);

  int nonzeros = countNonZero(m, symbfacIVL_ND);

  return nonzeros * 2 - m;
}

int Chordal::Spooles_MS(int m)
{
  int seed = 0, msglvl = 0;
  FILE   *fp = NULL;
  bool* bnode;

  int maxdomainsize = m / 16 + 1;
  
  //  rMessage("orderViaMMD:start");
  etree = orderViaMS(graph, maxdomainsize,seed, msglvl, fp);
  //  rMessage("orderViaMMD:end");
  newToOldIV_MS = ETree_newToOldVtxPerm(etree);
  symbfacIVL_MS = SymbFac_initFromGraph(etree,graph);
  //  IVL_writeForHumanEye(symbfacIVL_MS,stdout);

  int nonzeros = countNonZero(m, symbfacIVL_MS);

  return nonzeros * 2 - m;
}

void Chordal::ordering_bMat(int m, int nBlock,
                            InputData& inputData, 
                            FILE* fpOut)
{
  if ((m <= m_threshold)||(nBlock <= b_threshold)){
    best = -1;
    return;
  }
  for (int b=0; b<inputData.SDP_nBlock; b++){
    if (inputData.SDP_nConstraint[b] > m * sqrt(aggregate_threshold)){
      best = -1;
      return;
    }      
  }
  for (int b=0; b<inputData.SOCP_nBlock; b++){
    if (inputData.SOCP_nConstraint[b] > m * sqrt(aggregate_threshold)){
      best = -1;
      return;
    }      
  }
  for (int b=0; b<inputData.LP_nBlock; b++){
    if (inputData.LP_nConstraint[b] > m * sqrt(aggregate_threshold)){
      best = -1;
      return;
    }     
  }

  adjIVL = IVL_new();
  graph = Graph_new();

  makeGraph(inputData,m);

  if (IVL_tsize(adjIVL) > aggregate_threshold * m * m){
    best = -1;
    Graph_free(graph);
    return;
  }
#if PrintSparsity
  /* print sparsity information */
  printf("dense matrix               :\t\t\t%14d elements\n", m*m);
  fprintf(fpOut,"dense matrix               :\t\t\t%14d elements\n", m*m);
  printf("aggregate sparsity pattern :\t\t\t%14d elements\n",
	 IVL_tsize(adjIVL));
  fprintf(fpOut,"aggregate sparsity pattern :\t\t\t%14d elements\n",
          IVL_tsize(adjIVL));
#endif

  /* Uses METIS */
  if (Method[0]) {
    rError("no support for METIS");
  }
  
  /* Uses Spooles */
  if (Method[1])  {          /* Spooles 2.2 - minimum degree */
    Method[1] = Spooles_MMD(m);
    ETree_free(etree);
#if PrintSparsity
    printf("\tSpooles2.2 (minimum degree)\t\t%14d elements\n", Method[1]);
    fprintf(fpOut,"\tSpooles2.2 (minimum degree)\t\t%14d elements\n", Method[1]);
#endif
  }
  if (Method[2])  {         /* Spooles 2.2 - generalized nested dissection */
    Method[2] = Spooles_ND(m);
    ETree_free(etree);
#if PrintSparsity
    printf("\tSpooles2.2 (generalized nested dissection)%12d elements\n", Method[2]);
    fprintf(fpOut,"\tSpooles2.2 (generalized nested dissection)%12d elements\n", Method[2]);
#endif
  }
  if (Method[3])  {       /* Spooles 2.2 - multisection */
    Method[3] = Spooles_MS(m);
    ETree_free(etree);
#if PrintSparsity
    printf("\tSpooles2.2 (multisection)\t\t%14d elements\n", Method[3]);
    fprintf(fpOut,"\tSpooles2.2 (multisection)\t\t%14d elements\n", Method[3]);
#endif
  }
  if (Method[4])  {       /* Spooles 2.2 - best between nested
			     dissection and multisection */
    Method[4] = Spooles_NDMS(m);
    ETree_free(etree);
#if PrintSparsity
    printf("\tSpooles2.2 (best of ND and MS)\t\t%14d elements\n", Method[4]);
    fprintf(fpOut,"\tSpooles2.2 (best of ND and MS)\t\t%14d elements\n", Method[4]);
#endif
  }
  /* Select the best ordering */
  
  Graph_free(graph);
  
  best = Best_Ordering(Method);
  
  if (Method[best] > extend_threshold * m * m){
    best = -1;
  }

}

int Chordal::Best_Ordering(int *Method)
/************************************************************************
        Determine the best ordering.
************************************************************************/
{
        int   i, best;

        for (i = 0; Method[i] == 0; i++);
	best = i++;
	while (i < 5)
	{
	        for (; i < 5; i++)
		{
                        if (Method[i] != 0)  
                                break; 
                }
		if (i < 5)
		{
		        if (Method[i] < Method[best])
			        best = i;
			i++;
		}
        }
	return best;
}

}




