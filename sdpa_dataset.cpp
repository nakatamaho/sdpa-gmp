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

#include <sdpa_dataset.h>
#include <sdpa_parts.h>

namespace sdpa {

Solutions::Solutions()
{
  // Nothings needs.
}

Solutions::~Solutions()
{
  terminate();
}


Solutions::Solutions(int m,
		     int SDP_nBlock, int* SDP_blockStruct,
		     int SOCP_nBlock, int* SOCP_blockStruct,
		     int LP_nBlock, mpf_class lambda, ComputeTime& com)
{
  initialize(m,SDP_nBlock,SDP_blockStruct,
	     SOCP_nBlock,SOCP_blockStruct,
	     LP_nBlock, lambda, com);
}

void Solutions::initialize(int m,
			   int SDP_nBlock, int* SDP_blockStruct,
			   int SOCP_nBlock, int* SOCP_blockStruct,
			   int LP_nBlock, mpf_class lambda, ComputeTime& com)
{
  mDim = m;
  nDim = 0;
  for (int l=0; l<SDP_nBlock; ++l) {
    nDim += SDP_blockStruct[l];
  }
  for (int l=0; l<SOCP_nBlock; ++l) {
    nDim += SOCP_blockStruct[l];
  }
  nDim += LP_nBlock;

  xMat.initialize(SDP_nBlock,SDP_blockStruct,
		  SOCP_nBlock,SOCP_blockStruct,
		  LP_nBlock);
  xMat.setIdentity(lambda);
  zMat.initialize(SDP_nBlock,SDP_blockStruct,
		  SOCP_nBlock,SOCP_blockStruct,
		  LP_nBlock);
  zMat.setIdentity(lambda);
  yVec.initialize(m);
  yVec.setZero();

  invCholeskyX.initialize(SDP_nBlock,SDP_blockStruct,
			  SOCP_nBlock,SOCP_blockStruct,
			  LP_nBlock);
  invCholeskyX.setIdentity(1.0/sqrt(lambda));
  invCholeskyZ.initialize(SDP_nBlock,SDP_blockStruct,
			  SOCP_nBlock,SOCP_blockStruct,
			  LP_nBlock);
  invCholeskyZ.setIdentity(1.0/sqrt(lambda));
  invzMat.initialize(SDP_nBlock,SDP_blockStruct,
		     SOCP_nBlock,SOCP_blockStruct,
		     LP_nBlock);
  invzMat.setIdentity(1.0/lambda);
  //  invzMat.setIdentity(1.0/lambda);
}


void Solutions::terminate()
{
  xMat.terminate();
  zMat.terminate();
  yVec.terminate();
  invCholeskyX.terminate();
  invCholeskyZ.terminate();
  invzMat.terminate();
}

void Solutions::initializeZero(int m, 
			       int SDP_nBlock, int* SDP_blockStruct,
			       int SOCP_nBlock, int* SOCP_blockStruct,
			       int LP_nBlock, ComputeTime& com)
{
  // if we set initial point,
  // we malloc only the space of xMat, yVec, zMat
  xMat.initialize(SDP_nBlock,SDP_blockStruct,
		  SOCP_nBlock,SOCP_blockStruct,
		  LP_nBlock);
  xMat.setZero();
  zMat.initialize(SDP_nBlock,SDP_blockStruct,
		  SOCP_nBlock,SOCP_blockStruct,
		  LP_nBlock);
  zMat.setZero();
  yVec.initialize(m);
  yVec.setZero();
}


void Solutions::copyFrom(Solutions& other)
{
  if (this == &other) {
    return;
  }
  mDim = other.mDim;
  nDim = other.nDim;
  xMat.copyFrom(other.xMat);
  yVec.copyFrom(other.yVec);
  zMat.copyFrom(other.zMat);
  invCholeskyX.copyFrom(other.invCholeskyX);
  invCholeskyZ.copyFrom(other.invCholeskyZ);
  invzMat.copyFrom(other.invzMat);
}


bool Solutions::computeInverse(WorkVariables& work,
			       ComputeTime& com)
{

  bool total_judge = _SUCCESS;

  TimeStart(START1_3);
  if (Jal::getInvChol(invCholeskyX,xMat,work.DLS1) == false) {
    total_judge = FAILURE;
  }
  TimeEnd(END1_3);
  com.xMatTime += TimeCal(START1_3,END1_3);
  // rMessage(" xMat cholesky :: " << TimeCal(START1_3,END1_3));

  TimeStart(START1_4); 
  if (Jal::getInvCholAndInv(invCholeskyZ,invzMat,zMat,work.DLS2) == false) {
    total_judge = FAILURE;
  }
  TimeEnd(END1_4);
  // rMessage(" zMat cholesky :: " << TimeCal(START1_4,END1_4));
  com.zMatTime += TimeCal(START1_4,END1_4);
  
  xzMinEigenValue = 1.0;
  return total_judge;
}


bool Solutions::update(StepLength& alpha, Newton& newton,
		       WorkVariables& work,
		       ComputeTime& com)
{

  bool total_judge = _SUCCESS;

  TimeStart(START1_1);
  Lal::let(xMat,'=',xMat,'+',newton.DxMat,&alpha.primal);
  TimeEnd(END1_1);
  com.xMatTime += TimeCal(START1_1,END1_1);
  Lal::let(yVec,'=',yVec,'+',newton.DyVec,&alpha.dual);
  TimeStart(START1_2);
  Lal::let(zMat,'=',zMat,'+',newton.DzMat,&alpha.dual);
  TimeEnd(END1_2);
  com.zMatTime += TimeCal(START1_2,END1_2);

  const mpf_class cannot_move = 1.0e-4;
  if (alpha.primal < cannot_move && alpha.dual < cannot_move) {
    rMessage("Step length is too small. ");
    return FAILURE;
  }

  total_judge = computeInverse(work,com);

  return total_judge;
}





void Solutions::display(FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }
  fprintf(fpout,"dimension = \n",nDim);
  fprintf(fpout,"xMat = \n");
  xMat.display(fpout);
  fprintf(fpout,"yVec = \n");
  yVec.display(fpout);
  fprintf(fpout,"zMat = \n");
  zMat.display(fpout);
}


InputData::InputData()
{
  A = NULL;
  SDP_nBlock = 0;
  SDP_nConstraint = NULL;
  SDP_constraint = NULL;
  SDP_blockIndex = NULL;
  SOCP_nBlock = 0;
  SOCP_nConstraint = NULL;
  SOCP_constraint = NULL;
  SOCP_blockIndex = NULL;
  SDP_nBlock = 0;
  LP_nConstraint = NULL;
  LP_constraint = NULL;
  LP_blockIndex = NULL;
}

InputData::~InputData()
{
  terminate();
}

void InputData::initialize_bVec(int m)
{
  b.initialize(m);
}

  // not use         2008/02/29 kazuhide nakata
void InputData::initialize_CMat(int SDP_nBlock, int* SDP_blockStruct, 
				int* SDP_NonZeroNumber,
				int SOCP_nBlock, int* SOCP_blockStruct,
				int* SOCP_NonZeroNumber,
				int LP_nBlock, bool* LP_NonZeroNumber)
{
  C.initialize(SDP_nBlock, SDP_blockStruct, SDP_NonZeroNumber,
	       SOCP_nBlock, SOCP_blockStruct, SOCP_NonZeroNumber,
	       LP_nBlock, LP_NonZeroNumber);
}

  // not use         2008/02/29 kazuhide nakata
void InputData::initialize_AMat(int m)
{
  rNewCheck();
  A = new SparseLinearSpace[m];
  if (A==NULL) {
    rError("InputData:: Memory exhausted about blockIndex");
  }
}

  // not use         2008/02/29 kazuhide nakata
void InputData::initialize_AMat(int m, 
				int SDP_nBlock, int* SDP_blockStruct, 
				int* SDP_NonZeroNumber,
				int SOCP_nBlock, int* SOCP_blockStruct,
				int* SOCP_NonZeroNumber,
				int LP_nBlock, bool* LP_NonZeroNumber)
{
  rNewCheck();
  A = new SparseLinearSpace[m];
  if (A==NULL) {
    rError("InputData:: Memory exhausted about blockIndex");
  }
  for (int k=0 ; k<m; k++){
    A[k].initialize(SDP_nBlock, SDP_blockStruct, 
		    &SDP_NonZeroNumber[k*SDP_nBlock],
		    SOCP_nBlock, SOCP_blockStruct, 
		    &SOCP_NonZeroNumber[k*SOCP_nBlock],
		    LP_nBlock, &LP_NonZeroNumber[k*LP_nBlock]);
  }
}

void InputData::terminate()
{
  C.terminate();
  if (A){
    for (int k=0; k<b.nDim; ++k) {
      A[k].terminate();
    }
    delete[] A;
    A = NULL;
  }
  b.terminate();

  if (SDP_nConstraint && SDP_constraint && SDP_blockIndex){
    for (int k=0; k<SDP_nBlock; ++k) {
      delete[] SDP_constraint[k];SDP_constraint[k] = NULL;
      delete[] SDP_blockIndex[k];SDP_blockIndex[k] = NULL;
    }
    delete[] SDP_nConstraint;SDP_nConstraint = NULL;
    delete[] SDP_constraint;SDP_constraint = NULL;
    delete[] SDP_blockIndex;SDP_blockIndex = NULL;
  }
#if 0
  if (SOCP_nConstraint && SOCP_constraint && SOCP_blockIndex){
    for (int k=0; k<SOCP_nBlock; ++k) {
      delete[] SOCP_constraint[k];SOCP_constraint[k] = NULL;
      delete[] SOCP_blockIndex[k];SOCP_blockIndex[k] = NULL;
    }
    delete[] SOCP_nConstraint;SOCP_nConstraint = NULL;
    delete[] SOCP_constraint;SOCP_constraint = NULL;
    delete[] SOCP_blockIndex;SOCP_blockIndex = NULL;
  }
#endif
  if (LP_nConstraint && LP_constraint && LP_blockIndex){
    for (int k=0; k<LP_nBlock; ++k) {
      delete[] LP_constraint[k];LP_constraint[k] = NULL;
      delete[] LP_blockIndex[k];LP_blockIndex[k] = NULL;
    }
    delete[] LP_nConstraint;LP_nConstraint = NULL;
    delete[] LP_constraint;LP_constraint = NULL;
    delete[] LP_blockIndex;LP_blockIndex = NULL;
  }
}


void InputData::initialize_index_SDP(int SDP_nBlock, ComputeTime& com)
{
  int i,k;
  int mDim = b.nDim;
  int index;
  int* SDP_count;

  this->SDP_nBlock = SDP_nBlock;
  rNewCheck();
  SDP_nConstraint = new int[SDP_nBlock];
  if (SDP_nConstraint==NULL) {
    rError("InputData::initialize_index memory exhauseted ");
  }

  // count non-zero block matrix of A
  for (i=0; i<SDP_nBlock; i++){
      SDP_nConstraint[i] = 0;
  }
  for (i=0; i<mDim; i++){
    for (k=0; k<A[i].SDP_sp_nBlock; k++){
      index = A[i].SDP_sp_index[k];
      SDP_nConstraint[index]++;
    }
  }

  // malloc SDP_constraint, SDP_blockIndex
  rNewCheck();
  SDP_constraint = new int*[SDP_nBlock];
  if (SDP_constraint==NULL) {
    rError("InputData::initialize_index memory exhauseted ");
  }
  for (i=0; i<SDP_nBlock; i++){
    rNewCheck();
    SDP_constraint[i] = new int[SDP_nConstraint[i]];
    if (SDP_constraint[i]==NULL) {
      rError("InputData::initialize_index memory exhauseted ");
    }
  }
  rNewCheck();
  SDP_blockIndex = new int*[SDP_nBlock];
  if (SDP_blockIndex==NULL) {
    rError("InputData::initialize_index memory exhauseted ");
  }
  for (i=0; i<SDP_nBlock; i++){
    rNewCheck();
    SDP_blockIndex[i] = new int[SDP_nConstraint[i]];
    if (SDP_blockIndex[i]==NULL) {
      rError("InputData::initialize_index memory exhauseted ");
    }
  }

  // input index of non-zero block matrix of A
  rNewCheck()
  SDP_count = new int[SDP_nBlock];
  if (SDP_count == NULL){
      rError("InputData::initialize_index memory exhauseted ");
  }
  for (i=0; i<SDP_nBlock; i++){
      SDP_count[i] = 0;
  }
  for (i=0; i<mDim; i++){
    for (k=0; k<A[i].SDP_sp_nBlock; k++){
      index = A[i].SDP_sp_index[k];
      SDP_constraint[index][SDP_count[index]] = i;
      SDP_blockIndex[index][SDP_count[index]] = k;
      SDP_count[index]++;
    }
  }

  if (SDP_count){
    delete[] SDP_count;
    SDP_count = NULL;
  }
}


void InputData::initialize_index_SOCP(int SOCP_nBlock, ComputeTime& com)
{
  int i,k;
  int mDim = b.nDim;
  int index;
  int* SOCP_count;

  this->SOCP_nBlock = SOCP_nBlock;
  rNewCheck()
  SOCP_nConstraint = new int[SOCP_nBlock];
  if (SOCP_nConstraint==NULL) {
    rError("InputData::initialize_index memory exhauseted ");
  }

  // count non-zero block matrix of A
  for (i=0; i<SOCP_nBlock; i++){
      SOCP_nConstraint[i] = 0;
  }
  for (i=0; i<mDim; i++){
    for (k=0; k<A[i].SOCP_sp_nBlock; k++){
      index = A[i].SOCP_sp_index[k];
      SOCP_nConstraint[index]++;
    }
  }

  // malloc SOCP_constraint, SOCP_blockIndex
  rNewCheck()
  SOCP_constraint = new int*[SOCP_nBlock];
  if (SOCP_constraint==NULL) {
    rError("InputData::initialize_index memory exhauseted ");
  }
  for (i=0; i<SOCP_nBlock; i++){
    rNewCheck()
    SOCP_constraint[i] = new int[SOCP_nConstraint[i]];
    if (SOCP_constraint[i]==NULL) {
      rError("InputData::initialize_index memory exhauseted ");
    }
  }
  rNewCheck()
  SOCP_blockIndex = new int*[SOCP_nBlock];
  if (SOCP_blockIndex==NULL) {
    rError("InputData::initialize_index memory exhauseted ");
  }
  for (i=0; i<SOCP_nBlock; i++){
    rNewCheck()
    SOCP_blockIndex[i] = new int[SOCP_nConstraint[i]];
    if (SOCP_blockIndex[i]==NULL) {
      rError("InputData::initialize_index memory exhauseted ");
    }
  }

  // input index of non-zero block matrix of A
  rNewCheck()
  SOCP_count = new int[SOCP_nBlock];
    if (SOCP_count==NULL) {
      rError("InputData::initialize_index memory exhauseted ");
    }
  for (i=0; i<SOCP_nBlock; i++){
      SOCP_count[i] = 0;
  }
  for (i=0; i<mDim; i++){
    for (k=0; k<A[i].SOCP_sp_nBlock; k++){
      index = A[i].SOCP_sp_index[k];
      SOCP_constraint[index][SOCP_count[index]] = i;
      SOCP_blockIndex[index][SOCP_count[index]] = k;
      SOCP_count[index]++;
    }
  }

  if (SOCP_count){
    delete[] SOCP_count;
    SOCP_count = NULL;
  }
}


void InputData::initialize_index_LP(int LP_nBlock, ComputeTime& com)
{
  int i,k;
  int mDim = b.nDim;
  int index;
  int* LP_count;

  // for LP
  this->LP_nBlock = LP_nBlock;
  rNewCheck();
  LP_nConstraint = new int[LP_nBlock];
  if (LP_nConstraint==NULL) {
    rError("rInputData::initialize_index memory exhauseted ");
  }

  // count non-zero block matrix of A
  for (i=0; i<LP_nBlock; i++){
      LP_nConstraint[i] = 0;
  }
  for (i=0; i<mDim; i++){
    for (k=0; k<A[i].LP_sp_nBlock; k++){
      index = A[i].LP_sp_index[k];
      LP_nConstraint[index]++;
    }
  }

  // malloc LP_constraint, LP_blockIndex
  rNewCheck();
  LP_constraint = new int*[LP_nBlock];
  if (LP_constraint==NULL) {
    rError("InputData::initialize_index memory exhauseted ");
  }
  for (i=0; i<LP_nBlock; i++){
    rNewCheck();
    LP_constraint[i] = new int[LP_nConstraint[i]];
    if (LP_constraint[i]==NULL) {
      rError("InputData::initialize_index memory exhauseted ");
    }
  }
  rNewCheck();
  LP_blockIndex = new int*[LP_nBlock];
  if (LP_blockIndex==NULL) {
    rError("InputData::initialize_index memory exhauseted ");
  }
  for (i=0; i<LP_nBlock; i++){
    rNewCheck();
    LP_blockIndex[i] = new int[LP_nConstraint[i]];
    if (LP_blockIndex[i]==NULL) {
      rError("InputData::initialize_index memory exhauseted ");
    }
  }

  // input index of non-zero block matrix of A
  rNewCheck();
  LP_count = new int[LP_nBlock];
  if (LP_count==NULL) {
    rError("InputData::initialize_index memory exhauseted ");
  }
  for (i=0; i<LP_nBlock; i++){
      LP_count[i] = 0;
  }
  for (i=0; i<mDim; i++){
    for (k=0; k<A[i].LP_sp_nBlock; k++){
      index = A[i].LP_sp_index[k];
      LP_constraint[index][LP_count[index]] = i;
      LP_blockIndex[index][LP_count[index]] = k;
      LP_count[index]++;
    }
  }

  if (LP_count){
    delete[] LP_count;
    LP_count = NULL;
  }
}


void InputData::initialize_index(int SDP_nBlock,
				 int SOCP_nBlock,
				 int LP_nBlock,
				 ComputeTime& com)
{
  initialize_index_SDP(SDP_nBlock,com);
  //  initialize_index_SOCP(SOCP_nBlock,com);
  initialize_index_LP(LP_nBlock,com);
}

//   retVec_i := A_i bullet xMat (for i)
void InputData::multi_InnerProductToA(DenseLinearSpace& xMat, 
				      Vector& retVec)
{
  mpf_class ip;

  retVec.setZero();
  for (int i=0; i<retVec.nDim; i++){
    Lal::let(ip,'=',A[i],'.',xMat);
    retVec.ele[i] = ip;
  }    
}

  //   retMat := \sum_{i} A_i xVec_i
void InputData::multi_plusToA(Vector& xVec, DenseLinearSpace& retMat)
{
  retMat.setZero();
  for (int i=0; i<xVec.nDim; i++){
    Lal::let(retMat,'=',retMat,'+',A[i],&xVec.ele[i]);
  }    
}

void InputData::display(FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }

  fprintf(fpout,"b = \n");
  b.display(fpout);
  fprintf(fpout,"C = \n");
  C.display(fpout);
  for (int i=0; i<b.nDim; i++){
    fprintf(fpout,"A[%d] = \n",i);
    A[i].display(fpout);
  }
}

void InputData::display_index(FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }
  printf("display_index: %d %d %d\n",SDP_nBlock,SOCP_nBlock,LP_nBlock);

  for (int b=0; b<SDP_nBlock; b++){
    printf("SDP:%dth block\n",b);
    for (int i=0; i<SDP_nConstraint[b]; i++){
      printf("constraint:%d block:%d \n",
	     SDP_constraint[b][i],SDP_blockIndex[b][i]);
    }
  }

  for (int b=0; b<SOCP_nBlock; b++){
    printf("SOCP:%dth block\n",b);
    for (int i=0; i<SOCP_nConstraint[b]; i++){
      printf("constraint:%d block:%d \n",
	     SOCP_constraint[b][i],SOCP_blockIndex[b][i]);
    }
  }

  for (int b=0; b<LP_nBlock; b++){
    printf("LP:%dth block\n",b);
    for (int i=0; i<LP_nConstraint[b]; i++){
      printf("constraint:%d block:%d \n",
	     LP_constraint[b][i],LP_blockIndex[b][i]);
    }
  }


}




Residuals::Residuals()
{
  normPrimalVec = 0.0;
  normDualMat   = 0.0;
  centerNorm    = 0.0;
}

Residuals::Residuals(int m,
		     int SDP_nBlock, int* SDP_blockStruct,
		     int SOCP_nBlock, int* SOCP_blockStruct,
		     int LP_nBlock,
		     InputData& inputData, Solutions& currentPt)
{
  initialize(m,SDP_nBlock, SDP_blockStruct,
	     SOCP_nBlock ,SOCP_blockStruct,
	     LP_nBlock, inputData, currentPt);
}

Residuals::~Residuals()
{
  terminate();
}

void Residuals::initialize(int m, 
			   int SDP_nBlock, int* SDP_blockStruct,
			   int SOCP_nBlock, int* SOCP_blockStruct,
			   int LP_nBlock,
			   InputData& inputData, Solutions& currentPt)
{
  primalVec.initialize(m);
  dualMat.initialize(SDP_nBlock, SDP_blockStruct,
		     SOCP_nBlock ,SOCP_blockStruct,
		     LP_nBlock);
  compute(m, inputData, currentPt);
  
}

void Residuals::terminate()
{
  primalVec.terminate();
  dualMat.terminate();
}

void Residuals::copyFrom(Residuals& other)
{
  if (this==&other) {
    return;
  }
  primalVec.copyFrom(other.primalVec);
  dualMat.copyFrom(other.dualMat);
  normPrimalVec = other.normPrimalVec;
  normDualMat   = other.normDualMat;
  centerNorm    = other.centerNorm;
}

mpf_class Residuals::computeMaxNorm(Vector& primalVec)
{
  mpf_class ret = 0.0;
  #if 1
  for (int k=0; k<primalVec.nDim; ++k) {
    mpf_class tmp = abs(primalVec.ele[k]);
    if (tmp > ret) {
      ret = tmp;
    }
  }
  #else
  int index = gmp_idamax(primalVec.nDim,primalVec.ele,IONE);
  ret = fabs(primalVec.ele[index]);
  #endif
  return ret;
}

mpf_class Residuals::computeMaxNorm(DenseLinearSpace& dualMat)
{
  int SDP_nBlock = dualMat.SDP_nBlock;
  int SOCP_nBlock = dualMat.SOCP_nBlock;
  int LP_nBlock = dualMat.LP_nBlock;
  mpf_class ret = 0.0;
  mpf_class tmp;

  for (int l=0; l<SDP_nBlock; ++l) {
    mpf_class* target = dualMat.SDP_block[l].de_ele;
    int size = dualMat.SDP_block[l].nRow;
    for (int j=0; j<size*size; ++j) {
      tmp = abs(target[j]);
      if (tmp > ret) {
	ret = tmp;
      }
    }
  }

  for (int l=0; l<SOCP_nBlock; ++l) {
    rError("dataset:: current version do not support SOCP");
  }

  for (int l=0; l<LP_nBlock; ++l) {
    tmp = abs(dualMat.LP_block[l]);
    if (tmp > ret) {
      ret = tmp;
    }
  }

  return ret;
}

void Residuals::update(int m,
		       InputData& inputData,
		       Solutions& currentPt,
		       ComputeTime& com)
{
  TimeStart(UPDATE_START);
  compute(m,inputData,currentPt);
  TimeEnd(UPDATE_END);
  com.updateRes += TimeCal(UPDATE_START,UPDATE_END);
}

void Residuals::compute(int m,
			InputData& inputData,
			Solutions& currentPt)
{
  // p[k] = b[k] - A[k].X;
  inputData.multi_InnerProductToA(currentPt.xMat,primalVec);
  Lal::let(primalVec,'=',primalVec,'*',&MMONE);
  Lal::let(primalVec,'=',primalVec,'+',inputData.b);
  
  // D = C - Z - \sum A[k]y[k]
  inputData.multi_plusToA(currentPt.yVec, dualMat);
  Lal::let(dualMat,'=',dualMat,'*',&MMONE);
  Lal::let(dualMat,'=',dualMat,'+',inputData.C);
  Lal::let(dualMat,'=',dualMat,'-',currentPt.zMat);

  // rMessage("primal residual =");
  // primalVec.display();
  // rMessage("dual residual =");
  // dualMat.display();
  
  normPrimalVec = computeMaxNorm(primalVec);
  normDualMat   = computeMaxNorm(dualMat);
  centerNorm = 0.0;
}

void Residuals::display(FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }
  fprintf(fpout," currentRes.primalVec = \n");
  primalVec.display(fpout);
  fprintf(fpout," currentRes.dualMat = \n");
  dualMat.display(fpout);

  gmp_fprintf(fpout," currentRes.normPrimalVec = %8.3Fe\n",
	  normPrimalVec.get_mpf_t());
  gmp_fprintf(fpout," currentRes.normDualMat = %8.3Fe\n",
	  normDualMat.get_mpf_t());
}


WorkVariables::WorkVariables()
{
  // Nothings needs.
}

WorkVariables::WorkVariables(int m,
			     int SDP_nBlock, int* SDP_blockStruct,
			     int SOCP_nBlock, int* SOCP_blockStruct,
			     int LP_nBlock)
{
  initialize(m,SDP_nBlock,SDP_blockStruct,
	     SOCP_nBlock,SOCP_blockStruct,
	     LP_nBlock);
}

WorkVariables::~WorkVariables()
{
  terminate();
}


void WorkVariables::initialize(int m,
			       int SDP_nBlock, int* SDP_blockStruct,
			       int SOCP_nBlock, int* SOCP_blockStruct,
			       int LP_nBlock)
{
  DLS1.initialize(SDP_nBlock,SDP_blockStruct,
		  SOCP_nBlock,SOCP_blockStruct,
		  LP_nBlock);
  DLS2.initialize(SDP_nBlock,SDP_blockStruct,
		  SOCP_nBlock,SOCP_blockStruct,
		  LP_nBlock);
  DV1.initialize(m);
  DV1.initialize(m);

  if (SDP_nBlock > 0){
	SDP_BV1.initialize(SDP_nBlock,SDP_blockStruct);
	SDP_BV2.initialize(SDP_nBlock,SDP_blockStruct);
	SDP_BV3.initialize(SDP_nBlock,SDP_blockStruct);
	SDP_BV4.initialize(SDP_nBlock,SDP_blockStruct);
	SDP_BV5.initialize(SDP_nBlock,SDP_blockStruct);
	SDP_BV6.initialize(SDP_nBlock,SDP_blockStruct);
	SDP_BV7.initialize(SDP_nBlock,SDP_blockStruct);
	SDP_BV8.initialize(SDP_nBlock,SDP_blockStruct);
	SDP_BV9.initialize(SDP_nBlock,SDP_blockStruct);

	int* workStruct = NULL;
	rNewCheck();
	workStruct = new int[SDP_nBlock];
	if (workStruct == NULL) {
	  rMessage("WorkVariables :: memory exhausted");
	}

	for (int l=0; l<SDP_nBlock; ++l) {
	  workStruct[l] = max(1,3*SDP_blockStruct[l]-1);
	  //    workStruct[l] = max(1,2*SDP_blockStruct[l]-2);
	}
	SDP2_BV1.initialize(SDP_nBlock,workStruct);
	if (workStruct){
	  delete[] workStruct;
	  workStruct = NULL;
	}
  }
}

void WorkVariables::terminate()
{
  DLS1.terminate();
  DLS2.terminate();
  DV1.terminate();
  DV1.terminate();

  SDP_BV1.terminate();
  SDP_BV2.terminate();
  SDP_BV3.terminate();
  SDP_BV4.terminate();
  SDP_BV5.terminate();
  SDP_BV6.terminate();
  SDP_BV7.terminate();
  SDP_BV8.terminate();
  SDP_BV9.terminate();
  SDP2_BV1.terminate();
};



}
