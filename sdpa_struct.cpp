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

#include <sdpa_struct.h>

#define sdpa_dset(dset_length,dset_value,dset_pointer,dset_step) \
for (int dset_i=0,dset_index = 0; dset_i<dset_length; ++dset_i) { \
  dset_pointer[dset_index] = dset_value; \
  dset_index += dset_step; \
}

namespace sdpa{

Vector::Vector()
{
  nDim = 0;
  ele  = NULL;
}

Vector::Vector(int nDim, mpf_class value)
{
  ele  = NULL;
  initialize(nDim,value);
}

Vector::~Vector()
{
  terminate();
}

void Vector::initialize(int nDim,mpf_class value)
{
  // rMessage("Vector initialize");
  if (ele && this->nDim!=nDim) {
    if (ele) {
      delete[] ele;
      ele = NULL;
    }
    if (nDim<=0) {
      rError("Vector:: nDim is nonpositive");
    }
  }
  this->nDim = nDim;
  if (ele==NULL) {
    ele = NULL;
    rNewCheck();
    ele = new mpf_class[nDim];
    if (ele==NULL) {
      rError("Vector:: memory exhausted");
    }
  }
  sdpa_dset(nDim,value,ele,IONE);
}

void Vector::initialize(mpf_class value)
{
  if (nDim<=0) {
    rError("Vector:: nDim is nonpositive");
  }
  if (ele==NULL) {
    rNewCheck();
    ele = new mpf_class[nDim];
    if (ele==NULL) {
      rError("Vector:: memory exhausted");
    }
  }
  sdpa_dset(nDim,value,ele,IONE);
}

void Vector::terminate()
{
  if (ele) {
    delete[] ele;
  }
  ele = NULL;
}

void Vector::setZero()
{
  mpf_class zero=0.0;
  initialize(zero);
}

void Vector::display(FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }
  fprintf(fpout,"{");
  for (int j=0; j<nDim-1; ++j) {
    gmp_fprintf(fpout,P_FORMAT",",ele[j].get_mpf_t());
  }
  if (nDim>0) {
    gmp_fprintf(fpout,P_FORMAT"}\n",ele[nDim-1].get_mpf_t());
  } else {
    fprintf(fpout,"  }\n");
  }
}

void Vector::display(FILE* fpout,mpf_class scalar)
{
  if (fpout == NULL) {
    return;
  }
  fprintf(fpout,"{");
  for (int j=0; j<nDim-1; ++j) {
    mpf_class mtmp=ele[j]*scalar;
    gmp_fprintf(fpout,P_FORMAT",",mtmp.get_mpf_t());
  }
  if (nDim>0) {
    mpf_class mtmp=ele[nDim-1]*scalar;
    gmp_fprintf(fpout,P_FORMAT"}\n",mtmp.get_mpf_t());
  } else {
    fprintf(fpout,"  }\n");
  }
}

bool Vector::copyFrom(Vector& other)
{
  if (this == &other) {
    return _SUCCESS;
  }
  if (nDim != other.nDim && ele!=NULL) {
    delete[] ele;
    ele = NULL;
  }
  nDim = other.nDim;
  if (nDim<=0) {
    rError("Vector:: nDim is nonpositive");
  }
  if (ele==NULL) {
    rNewCheck();
    ele = new mpf_class[nDim];
    if (ele==NULL) {
      rError("Vector:: memory exhausted");
    }
  }
  Rcopy(nDim,other.ele,1,ele,1);
  return _SUCCESS;
}

BlockVector::BlockVector()
{
  nBlock = 0;
  blockStruct = NULL;
  ele = NULL;
}

BlockVector::BlockVector(int nBlock, int* blockStruct,
			   mpf_class value)
{
  initialize(nBlock,blockStruct,value);
}

BlockVector::~BlockVector()
{
  terminate();
}

void BlockVector::initialize(int nBlock, int* blockStruct,
			      mpf_class value)
{
  // rMessage("BlockVector initialize");
  this->nBlock = nBlock;
  if (nBlock<=0) {
    rError("BlockVector:: nBlock is nonpositive");
  }
  this->blockStruct = NULL;
  rNewCheck();
  this->blockStruct = new int[nBlock];
  if (this->blockStruct==NULL) {
    rError("BlockVector:: memory exhausted");
  }
  for (int l=0; l<nBlock; ++l) {
    this->blockStruct[l] = blockStruct[l];
  }

  ele = NULL;
  rNewCheck();
  ele = new Vector[nBlock];
  if (ele==NULL) {
    rError("BlockVector:: memory exhausted");
  }
  for (int l=0; l<nBlock; ++l) {
    int size = blockStruct[l];
    if (size<0) {
      size = -size;
    }
    ele[l].initialize(size,value);
  }
}

void BlockVector::initialize(mpf_class value)
{
  if (nBlock>0 && blockStruct && ele) {
    for (int l=0; l<nBlock; ++l) {
      ele[l].initialize(value);
    }
  }
}

void BlockVector::terminate()
{
  if (ele && blockStruct && nBlock>=0) {
    for (int l=0; l<nBlock; ++l) {
      ele[l].terminate();
    }
    delete[] ele;
    ele = NULL;

    delete[] blockStruct;
    blockStruct = NULL;
  }
}

void BlockVector::setZero()
{
  if (nBlock>0 && blockStruct && ele) {
    for (int l=0; l<nBlock; ++l) {
      ele[l].setZero();
    }
  }
}

void BlockVector::display(FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }
  fprintf(fpout,"{ ");
  if (nBlock>0 && blockStruct && ele) {
    for (int l=0; l<nBlock; ++l) {
      ele[l].display(fpout);
    }
  }
  fprintf(fpout,"} \n");
}

bool BlockVector::copyFrom(BlockVector& other)
{
  if (this == &other) {
    return _SUCCESS;
  }
  
  if (other.nBlock<=0) {
    rError("BlockVector:: nBlock is nonpositive");
  }
  if (nBlock!=other.nBlock && blockStruct) {
    delete[] blockStruct;
    blockStruct = NULL;
    delete[] ele;
    ele = NULL;
  }
  if (blockStruct==NULL) {
    nBlock = other.nBlock;
    rNewCheck();
    blockStruct = new int[nBlock];
    if (blockStruct==NULL) {
      rError("BlockVector:: memory exhausted");
    }
    for (int l=0; l<nBlock; ++l) {
      blockStruct[l] = other.blockStruct[l];
    }
  }
  if (ele==NULL) {
    rNewCheck();
    ele = new Vector[nBlock];
    if (ele==NULL) {
      rError("BlockVector:: memory exhausted");
    }
  }
  for (int l=0; l<nBlock; ++l) {
    ele[l].copyFrom(other.ele[l]);
  }
  return _SUCCESS;
}

SparseMatrix::SparseMatrix()
{
  nRow = 0;
  nCol = 0;
  type = SPARSE;

  NonZeroNumber = 0;
  
  de_ele = NULL;

  row_index     = NULL;
  column_index  = NULL;
  sp_ele        = NULL;
  NonZeroCount  = 0;
  NonZeroEffect = 0;
}

SparseMatrix::SparseMatrix(int nRow, int nCol,
			     SparseMatrix::Type type,
			     int NonZeroNumber)
{
  initialize(nRow, nCol, type, NonZeroNumber);
}

SparseMatrix::~SparseMatrix()
{
  terminate();
}

void SparseMatrix::
initialize(int nRow, int nCol,
	   SparseMatrix::Type type,
	   int NonZeroNumber)
{
  // rMessage("SparseMatrix initialize");

  SparseMatrix();
  if (nRow<=0 || nCol<=0) {
    rError("SparseMatrix:: Dimensions are nonpositive");
  }
  this->nRow          = nRow;
  this->nCol          = nCol;
  this->type          = type;

  int length;
  switch(type) {
  case SPARSE:
    this->NonZeroNumber  = NonZeroNumber;
    this->NonZeroCount   = 0;
    this->NonZeroEffect  = 0;
    if (NonZeroNumber > 0) {
      rNewCheck();
      row_index    = new int[NonZeroNumber];
      rNewCheck();
      column_index = new int[NonZeroNumber];
      rNewCheck();
      sp_ele       = new mpf_class[NonZeroNumber];
      if (row_index==NULL || column_index==NULL
	  || sp_ele==NULL) {
	rError("SparseMatrix:: memory exhausted");
      }
    }
    break;
  case DENSE:
    this->NonZeroNumber = nRow*nCol;
    this->NonZeroCount  = nRow*nCol;
    this->NonZeroEffect = nRow*nCol;
    rNewCheck();
    de_ele = new mpf_class[NonZeroNumber];
    if (de_ele==NULL) {
      rError("SparseMatrix:: memory exhausted");
    }
    length = nRow*nCol;
    sdpa_dset(length,MZERO,de_ele,IONE);
    // all elements are 0.
    break;
  }
}

void SparseMatrix::terminate()
{
  if (de_ele) {
    delete[] de_ele;
    de_ele = NULL;
  }
  if (row_index) {
    delete[] row_index;
    row_index = NULL;
  }
  if (column_index) {
    delete[] column_index;
    column_index = NULL;
  }
  if (sp_ele) {
    delete[] sp_ele;
    sp_ele = NULL;
  }
}

void SparseMatrix::display(FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }
  switch(type) {
  case SPARSE:
    fprintf(fpout,"{");
    for (int index=0; index<NonZeroCount; ++index) {
      int i        = row_index[index];
      int j        = column_index[index];
      mpf_class value = sp_ele[index];
      gmp_fprintf(fpout,"val[%d,%d] = "P_FORMAT"\n", i,j,value.get_mpf_t());
    }
    fprintf(fpout,"}\n");
    break;
  case DENSE:
    fprintf(fpout,"{\n");
    for (int i=0; i<nRow-1; ++i) {
      if (i==0) {
	fprintf(fpout," ");
      } else {
	fprintf(fpout,"  ");
      }
      fprintf(fpout,"{");
      for (int j=0; j<nCol-1; ++j) {
        gmp_fprintf(fpout, P_FORMAT",",de_ele[i+nCol*j].get_mpf_t());
      }
       gmp_fprintf(fpout,P_FORMAT" },\n",de_ele[i+nCol*(nCol-1)].get_mpf_t());
    }
    if (nRow>1) {
      fprintf(fpout,"  {");
    }
    for (int j=0; j<nCol-1; ++j) {
      gmp_fprintf(fpout,P_FORMAT",",de_ele[(nRow-1)+nCol*j].get_mpf_t());
    }
    gmp_fprintf(fpout,P_FORMAT" }",de_ele[(nRow-1)+nCol*(nCol-1)].get_mpf_t());
    if (nRow>1) {
      fprintf(fpout,"   }\n");
    } else {
      fprintf(fpout,"\n");
    }
    break;
  }
}

bool SparseMatrix::copyFrom(SparseMatrix& other)
{
  if (type != other.type || nRow != other.nRow
      || nCol != other.nCol) {
    this->~SparseMatrix();
    initialize(other.nRow,other.nCol,other.type,
	       NonZeroNumber);
    NonZeroCount  = other.NonZeroCount;
    NonZeroEffect = other.NonZeroEffect;
    int length;
    switch(type) {
    case SPARSE:
      for (int index = 0; index<NonZeroCount;++index) {
	row_index[index]    = other.row_index[index];
	column_index[index] = other.column_index[index];
	sp_ele[index]       = other.sp_ele[index];
      }
      break;
    case DENSE:
      length = nRow*nCol;
      Rcopy(length,other.de_ele,1,de_ele,1);
      break;
    }
  } else { // Sp_De_Di == other.Sp_De_Di
           // && nRow == other.nRow && nCol == other.nCol
    NonZeroCount  = other.NonZeroCount;
    NonZeroEffect = other.NonZeroEffect;
    int length;
    switch(type) {
    case SPARSE:
      if (NonZeroNumber!=other.NonZeroNumber) {
	delete[] row_index;
	delete[] column_index;
	delete[] sp_ele;
	row_index = column_index = NULL;
	sp_ele = NULL;
	rNewCheck();
	row_index    = new int[NonZeroNumber];
	rNewCheck();
	column_index = new int[NonZeroNumber];
	rNewCheck();
	sp_ele       = new mpf_class[NonZeroNumber];
	if (row_index==NULL || column_index==NULL
	    || sp_ele==NULL) {
	  rError("SparseMatrix:: memory exhausted");
	}
      }
      for (int index = 0; index<NonZeroCount;++index) {
	row_index[index]    = other.row_index[index];
	column_index[index] = other.column_index[index];
	sp_ele[index]       = other.sp_ele[index];
      }
      break;
    case DENSE:
      length = nRow*nCol;
      Rcopy(length,other.de_ele,1,de_ele,1);
      break;
    } // end of switch
  } // end of else
  return _SUCCESS;
}

void SparseMatrix::changeToDense(bool forceChange)
{
  if (type!=SPARSE) {
    return;
  }
  // if (false)
  // rMessage(" NonZeroCount " << NonZeroCount);
  // rMessage(" nRow*nCol*0.2 " << nRow*nCol*0.2);
  if (forceChange == false && NonZeroCount < (nRow*nCol) * 0.20) {
    // if the number of elements are less than 20 percent,
    // we don't change to Dense.
    return;
  }
  // rMessage("change");
  type = DENSE;
  de_ele = NULL;
  int length = nRow*nCol;
  rNewCheck();
  de_ele = new mpf_class[length];
  if (de_ele==NULL) {
    rError("SparseMatrix:: memory exhausted");
  }
  sdpa_dset(length,MZERO,de_ele,1);
  // all elements are set 0.
  for (int index=0; index<NonZeroCount; ++index) {
    int        i = row_index[index];
    int        j = column_index[index];
    mpf_class value = sp_ele[index];
    if (i==j) {
      de_ele[i+nCol*j] = value;
    } else {
      de_ele[i+nCol*j] = de_ele[j+nCol*i] = value;
    }
  }
  NonZeroCount = NonZeroNumber = NonZeroEffect = length;
  delete[] row_index;
  delete[] column_index;
  delete[] sp_ele;
  row_index = NULL;
  column_index = NULL;
  sp_ele = NULL;
}

void SparseMatrix::setZero()
{
  int length;
  switch(type) {
  case SPARSE:
    NonZeroCount  = 0;
    NonZeroEffect = 0;
    // No element is stored.
    break;
  case DENSE:
    length = nRow*nCol;
    sdpa_dset(length,MZERO,de_ele,IONE);
    break;
  }
}

void SparseMatrix::setIdentity(mpf_class scalar)
{
  if (nRow != nCol) {
    rError("SparseMatrix:: Identity matrix must be square matrix");
  }
  int length,step;
  switch(type) {
  case SPARSE:
    if (nCol > NonZeroNumber) {
      rError("SparseMatrix:: cannot store over NonZeroNumber");
      // the number of Diagonal elements equals nCol.
    }
    NonZeroCount  = nCol;
    NonZeroEffect = nCol;
    for (int index=0; index< NonZeroCount; ++index) {
      row_index[index]    = index;
      column_index[index] = index;
      sp_ele[index]       = scalar;
    }
    break;
  case DENSE:
    length = nRow*nCol;
    sdpa_dset(length,MZERO,de_ele,IONE);
    step = nCol+1;
    sdpa_dset(nCol,scalar,de_ele,step);
    // only diagonal elements are the value of scalar.
    break;
  }
}
    
bool SparseMatrix::sortSparseIndex(int& i, int& j)
{
  // if this matrix is not symmetric,
  // return the index(i,j) whose values are not symmetric.
  i = -1;
  j = -1;
  const mpf_class tolerance = 1.0e-8;
  switch(type) {
  case SPARSE:
    // Make matrix as Upper Triangluar
    for (int i1=0; i1<NonZeroCount; ++i1) {
      int tmpi = row_index[i1];
      int tmpj = column_index[i1];
      if (tmpi>tmpj) {
	row_index   [i1] = tmpj;
	column_index[i1] = tmpi;
      }
    }
    // simple sort
    for (int i1=0; i1<NonZeroCount; ++i1) {
      for (int i2=0; i2<i1; ++i2) {
	int index1 = row_index[i1]+nCol*column_index[i1];
	int index2 = row_index[i2]+nCol*column_index[i2];
	if (index1<index2) {
	  int         tmpi = row_index   [i2];
	  int         tmpj = column_index[i2];
	  mpf_class      tmpv = sp_ele      [i2];
	  row_index   [i2] = row_index   [i1];
	  column_index[i2] = column_index[i1];
	  sp_ele      [i2] = sp_ele      [i1];
	  row_index   [i1] = tmpi;
	  column_index[i1] = tmpj;
	  sp_ele      [i1] = tmpv;
	}
      }
    }
    // the process for the same index
    for (int i1=0; i1<NonZeroCount-1; ++i1) {
      int index1 = row_index[i1  ]+nCol*column_index[i1  ];
      int index2 = row_index[i1+1]+nCol*column_index[i1+1];
      if (index1 == index2) {
	if (abs(sp_ele[index1] - sp_ele[index2]) > tolerance) {
	  // Here must not be symmetric
	  if (i<0 || j<0) {
	    i = row_index   [i1];
	    j = column_index[i1];
	  }
	}

	// remove redudunt
	for (int i2 = i1+1; i2<NonZeroCount-2;++i2) {
	  row_index   [i2] = row_index   [i2+1];
	  column_index[i2] = column_index[i2+1];
	  sp_ele      [i2] = sp_ele      [i2+1];
	}
	NonZeroCount--;
	if (i==j) {
	  NonZeroEffect--;
	} else {
	  NonZeroEffect -= 2;
	}
      } // end of 'if (index1==index2)'
    }
    break;
  case DENSE:
    if (nRow!=nCol) {
      return FAILURE;
    }
    for (j=1; j<nCol; ++j) {
      for (i=0; i<j; ++i) {
	if (abs(de_ele[i+nCol*j]-de_ele[j+nCol*i]) > tolerance) {
	  return FAILURE;
	}
      }
    }
    break;
  }
  return _SUCCESS;
}

DenseMatrix::DenseMatrix()
{
  nRow = 0;
  nCol = 0;
  type = DENSE;

  de_ele = NULL;
}

DenseMatrix::DenseMatrix(int nRow, int nCol,
			   DenseMatrix::Type type)
{
  initialize(nRow, nCol, type);
}

DenseMatrix::~DenseMatrix()
{
  terminate();
}

void DenseMatrix::
initialize(int nRow, int nCol,
	   DenseMatrix::Type type)
{
  // rMessage("DenseMatrix::initialize");

  DenseMatrix();
  if (nRow<=0 || nCol<=0) {
    rError("DenseMatrix:: Dimensions are nonpositive");
  }
  int old_length = this->nRow*this->nCol;
  this->nRow  = nRow;
  this->nCol  = nCol;

  int length;
  switch(type) {
  case DENSE:
    length = nRow*nCol;
    if (de_ele && old_length!=length) {
      delete[] de_ele;
      de_ele = NULL;
    }
    if (de_ele==NULL) {
      rNewCheck();
      de_ele = new mpf_class[length];
      if (de_ele==NULL) {
	rError("DenseMatrix:: memory exhausted");
      }
    }
    sdpa_dset(length,MZERO,de_ele,IONE);
    break;
  case COMPLETION:
    rError("DenseMatrix:: no support for COMPLETION");
    break;
  }
}

void DenseMatrix::terminate()
{
  if (de_ele) {
    delete[] de_ele;
    de_ele = NULL;
  }
}

void DenseMatrix::display(FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }
  switch(type) {
  case DENSE:
    fprintf(fpout,"{");
    for (int i=0; i<nRow-1; ++i) {
      if (i==0) {
	fprintf(fpout," ");
      } else {
	fprintf(fpout,"  ");
      }
      fprintf(fpout,"{");
      for (int j=0; j<nCol-1; ++j) {
        gmp_fprintf(fpout, P_FORMAT",",de_ele[i+nCol*j].get_mpf_t());
      }
      gmp_fprintf(fpout,P_FORMAT" },\n",de_ele[i+nCol*(nCol-1)].get_mpf_t());
    }
    if (nRow>1) {
      fprintf(fpout,"  {");
    }
    for (int j=0; j<nCol-1; ++j) {
      gmp_fprintf(fpout,P_FORMAT",",de_ele[(nRow-1)+nCol*j].get_mpf_t());
    }
    gmp_fprintf(fpout,P_FORMAT" }",de_ele[(nRow-1)+nCol*(nCol-1)].get_mpf_t());
    if (nRow>1) {
      fprintf(fpout,"   }\n");
    } else {
      fprintf(fpout,"\n");
    }
    break;
  case COMPLETION:
    rError("DenseMatrix:: no support for COMPLETION");
    break;
  }
}

bool DenseMatrix::copyFrom(SparseMatrix& other)
{
  int length;
  switch(other.type) {
  case SparseMatrix::SPARSE:
    type = DENSE;
    if (de_ele) {
      delete[] de_ele;
    }
    de_ele = NULL;
    nRow = other.nRow;
    nCol = other.nCol;
    rNewCheck();
    de_ele = new mpf_class[nRow*nCol];
    if (de_ele==NULL) {
      rError("DenseMatrix:: memory exhausted");
    }
    length = nRow*nCol;
    sdpa_dset(length,MZERO,de_ele,IONE);
    for (int index = 0; index<other.NonZeroCount; ++index) {
      int i = other.row_index[index];
      int j = other.column_index[index];
      mpf_class value = other.sp_ele[index];
      de_ele[i+nCol*j] = de_ele[j+nCol*i] = value;
    }
    break;
  case SparseMatrix::DENSE:
    type = DENSE;
    if (de_ele && (other.nRow!=nRow || other.nCol!=nCol)) {
      delete[] de_ele;
      de_ele = NULL;
    }
    nRow = other.nRow;
    nCol = other.nCol;
    rNewCheck();
    de_ele = new mpf_class[nRow*nCol];
    if (de_ele==NULL) {
      rError("DenseMatrix:: memory exhausted");
    }
    length = nRow*nCol;
    Rcopy(length,other.de_ele,1,de_ele,1);
    break;
  }
  return _SUCCESS;
}

bool DenseMatrix::copyFrom(DenseMatrix& other)
{
  if (this == &other) {
    return _SUCCESS;
  }
  int length;
  switch(other.type) {
  case DENSE:
    type = DENSE;
    if (de_ele && (other.nRow!=nRow || other.nCol!=nCol)) {
      delete[] de_ele;
      de_ele = NULL;
    }
    nRow = other.nRow;
    nCol = other.nCol;
    if (de_ele==NULL) {
      rNewCheck();
      de_ele = new mpf_class[nRow*nCol];
      if (de_ele==NULL) {
	rError("DenseMatrix:: memory exhausted");
      }
    }
    length = nRow*nCol;
    Rcopy(length,other.de_ele,1,de_ele,1);
    break;
  case COMPLETION:
    rError("DenseMatrix:: no support for COMPLETION");
    break;
  }
  return _SUCCESS;
}


void DenseMatrix::setZero()
{
  int length;
  switch(type) {
  case DENSE:
    length = nRow*nCol;
    sdpa_dset(length,MZERO,de_ele,IONE);
    break;
  case COMPLETION:
    rError("DenseMatrix:: no support for COMPLETION");
    break;
  }
}

void DenseMatrix::setIdentity(mpf_class scalar)
{
  if (nRow != nCol) {
    rError("SparseMatrix:: Identity matrix must be square matrix");
  }
  int length,step;
  switch(type) {
  case DENSE:
    length = nRow*nCol;
    sdpa_dset(length,MZERO,de_ele,IONE);
    step = nCol+1;
    sdpa_dset(nCol,scalar,de_ele,step);
    break;
  case COMPLETION:
    rError("DenseMatrix:: no support for COMPLETION");
    break;
  }
}
    
SparseLinearSpace::SparseLinearSpace()
{
  SDP_sp_nBlock = 0;
  SDP_sp_index = NULL;
  SDP_sp_block = NULL;
  SOCP_sp_nBlock = 0;
  SOCP_sp_index = NULL;
  SOCP_sp_block = NULL;
  LP_sp_nBlock = 0;
  LP_sp_index = NULL;
  LP_sp_block = NULL;
}

SparseLinearSpace::SparseLinearSpace(int SDP_nBlock, 
				   int* SDP_blockStruct,
				   int* SDP_NonZeroNumber,
				   int SOCP_nBlock, 
				   int* SOCP_blockStruct,
				   int* SOCP_NonZeroNumber,
				   int LP_nBlock, 
				   bool* LP_NonZeroNumber)
{
  initialize(SDP_nBlock, SDP_blockStruct, SDP_NonZeroNumber,
	     SOCP_nBlock, SOCP_blockStruct, SOCP_NonZeroNumber,
	     LP_nBlock, LP_NonZeroNumber);
}

SparseLinearSpace::SparseLinearSpace(int SDP_sp_nBlock, 
                                     int* SDP_sp_index,
                                     int* SDP_sp_blockStruct, 
                                     int* SDP_sp_NonZeroNumber,
                                     int SOCP_sp_nBlock, 
                                     int* SOCP_sp_index,
                                     int* SOCP_sp_blockStruct,
                                     int* SOCP_sp_NonZeroNumber,
                                     int LP_sp_nBlock, 
                                     int* LP_sp_index)
{
  initialize(SDP_sp_nBlock, SDP_sp_index,
             SDP_sp_blockStruct, SDP_sp_NonZeroNumber,
	     SOCP_sp_nBlock, SOCP_sp_index,
             SOCP_sp_blockStruct, SOCP_sp_NonZeroNumber,
	     LP_sp_nBlock, LP_sp_index);
}

SparseLinearSpace::~SparseLinearSpace()
{
  terminate();
}

// dense form of block index
void SparseLinearSpace::initialize(int SDP_nBlock, 
				   int* SDP_blockStruct,
				   int* SDP_NonZeroNumber,
				   int SOCP_nBlock, 
				   int* SOCP_blockStruct,
				   int* SOCP_NonZeroNumber,
				   int LP_nBlock, 
				   bool* LP_NonZeroNumber)
{
  // rMessage("SparseLinearSpace::initialize");
  SDP_sp_nBlock = 0;
  SOCP_sp_nBlock = 0;
  LP_sp_nBlock = 0;
  int counter;

  // for SDP
  for (int l=0; l<SDP_nBlock; l++){
    if (SDP_NonZeroNumber[l] > 0){
      SDP_sp_nBlock++;
    }    
  }
  if (SDP_sp_nBlock > 0){
	SDP_sp_index = NULL;
	rNewCheck();
	SDP_sp_index = new int[SDP_sp_nBlock];
	if (SDP_sp_index==NULL) {
	  rError("SDP_index:: memory exhausted");
	}
	SDP_sp_block = NULL;
	rNewCheck();
	SDP_sp_block = new SparseMatrix[SDP_sp_nBlock];
	if (SDP_sp_block==NULL) {
	  rError("SparseLinearSpace:: memory exhausted");
	}
  }
  counter = 0;
  for (int l=0; l<SDP_nBlock; ++l) {
    if (SDP_NonZeroNumber[l] > 0){
      SDP_sp_index[counter] = l;
      int size = SDP_blockStruct[l];
      SDP_sp_block[counter].initialize(size,size,SparseMatrix::SPARSE,
				    SDP_NonZeroNumber[l]);
      counter++;
    }    
  }


  // for SOCP
#if 0
  for (int l=0; l<SOCP_nBlock; l++){
    if (SOCP_NonZeroNumber[l] > 0){
      SOCP_sp_nBlock++;
    }    
  }
  if (SOCP_sp_nBlock > 0){
	SOCP_sp_index = NULL;
	rNewCheck();
	SOCP_sp_index = new int[SOCP_sp_nBlock];
	if (SOCP_sp_index==NULL) {
	  rError("SOCP_index:: memory exhausted");
	}
	SOCP_sp_block = NULL;
	rNewCheck();
	SOCP_sp_block = new SparseMatrix[SOCP_sp_nBlock];
	if (SOCP_sp_block==NULL) {
	  rError("SparseLinearSpace:: memory exhausted");
	}
  }
  counter = 0;
  for (int l=0; l<SOCP_nBlock; ++l) {
    if (SOCP_NonZeroNumber[l] > 0){
      SOCP_sp_index[counter] = l;
      int size = SOCP_blockStruct[l];
      SOCP_sp_block[counter].initialize(size,size,SparseMatrix::SPARSE,
				    SOCP_NonZeroNumber[l]);
      counter++;
    }    
  }
#endif

  // for LP
  for (int l=0; l<LP_nBlock; l++){
    if (LP_NonZeroNumber[l] == true){
      LP_sp_nBlock++;
    }    
  }
  if (LP_sp_nBlock > 0){
	LP_sp_index = NULL;
	rNewCheck();
	LP_sp_index = new int[LP_sp_nBlock];
	if (LP_sp_index==NULL) {
	  rError("LP_index:: memory exhausted");
	}
	LP_sp_block = NULL;
	rNewCheck();
	LP_sp_block = new mpf_class[LP_sp_nBlock];
	if (LP_sp_block==NULL) {
	  rError("SparseLinearSpace:: memory exhausted");
	}
  }
  counter = 0;
  for (int l=0; l<LP_nBlock; ++l) {
    if (LP_NonZeroNumber[l] == true){
      LP_sp_index[counter] = l;
      counter++;
    }    
  }

}

// sparse form of block index      2008/02/27 kazuhide nakata
void SparseLinearSpace::initialize(int SDP_sp_nBlock, 
                                   int* SDP_sp_index,
                                   int* SDP_sp_blockStruct, 
                                   int* SDP_sp_NonZeroNumber,
                                   int SOCP_sp_nBlock, 
                                   int* SOCP_sp_index,
                                   int* SOCP_sp_blockStruct,
                                   int* SOCP_sp_NonZeroNumber,
                                   int LP_sp_nBlock, 
                                   int* LP_sp_index)
{
  // rMessage("SparseLinearSpace::initialize");

  // for SDP
  this->SDP_sp_nBlock = SDP_sp_nBlock;
  if (SDP_sp_nBlock > 0){
	this->SDP_sp_index = NULL;
	rNewCheck();
	this->SDP_sp_index = new int[SDP_sp_nBlock];
	if (this->SDP_sp_index==NULL) {
	  rError("SDP_index:: memory exhausted");
	}
	this->SDP_sp_block = NULL;
	rNewCheck();
	this->SDP_sp_block = new SparseMatrix[SDP_sp_nBlock];
	if (this->SDP_sp_block==NULL) {
	  rError("SparseLinearSpace:: memory exhausted");
	}
  }
  for (int l=0; l<SDP_sp_nBlock; ++l) {
    this->SDP_sp_index[l] = SDP_sp_index[l];
    int size = SDP_sp_blockStruct[l];
    SDP_sp_block[l].initialize(size,size,SparseMatrix::SPARSE,
                               SDP_sp_NonZeroNumber[l]);
  }

  // for SOCP
#if 0
  this->SOCP_sp_nBlock = SOCP_sp_nBlock;
  if (SOCP_sp_nBlock > 0){
	this->SOCP_sp_index = NULL;
	rNewCheck();
	this->SOCP_sp_index = new int[SOCP_sp_nBlock];
	if (this->SOCP_sp_index==NULL) {
	  rError("SOCP_index:: memory exhausted");
	}
	this->SOCP_sp_block = NULL;
	rNewCheck();
	this->SOCP_sp_block = new SparseMatrix[SOCP_sp_nBlock];
	if (this->SOCP_sp_block==NULL) {
	  rError("SparseLinearSpace:: memory exhausted");
	}
  }
  for (int l=0; l<SOCP_sp_nBlock; ++l) {
    this->SOCP_sp_index[l] = SOCP_sp_index[l];
    int size = SOCP_sp_blockStruct[l];
    SOCP_sp_block[l].initialize(size,size,SparseMatrix::SPARSE,
                               SOCP_sp_NonZeroNumber[l]);
  }
#endif

  // for LP
  this->LP_sp_nBlock = LP_sp_nBlock;
  if (LP_sp_nBlock > 0){
	this->LP_sp_index = NULL;
	rNewCheck();
	this->LP_sp_index = new int[LP_sp_nBlock];
	if (this->LP_sp_index==NULL) {
	  rError("LP_index:: memory exhausted");
	}
	this->LP_sp_block = NULL;
	rNewCheck();
	this->LP_sp_block = new mpf_class[LP_sp_nBlock];
	if (this->LP_sp_block==NULL) {
	  rError("SparseLinearSpace:: memory exhausted");
	}
  }
  for (int l=0; l<LP_sp_nBlock; ++l) {
    this->LP_sp_index[l] = LP_sp_index[l];
  }

}

void SparseLinearSpace::terminate()
{
  // for SDP
  if (SDP_sp_block && SDP_sp_index && SDP_sp_nBlock>=0) {
    for (int l=0; l<SDP_sp_nBlock; ++l) {
      SDP_sp_block[l].terminate();
    }
    delete[] SDP_sp_block;
    SDP_sp_block = NULL;
    
    delete[] SDP_sp_index;
    SDP_sp_index = NULL;
  }
  // for SOCP
#if 0
  if (SOCP_sp_block && SOCP_sp_index && SOCP_sp_nBlock>=0) {
    for (int l=0; l<SOCP_sp_nBlock; ++l) {
      SOCP_sp_block[l].terminate();
    }
    delete[] SOCP_sp_block;
    SOCP_sp_block = NULL;
    
    delete[] SOCP_sp_index;
    SOCP_sp_index = NULL;
  }
#endif 
  // for LP
  if (LP_sp_block && LP_sp_index && LP_sp_nBlock>=0) {
    delete[] LP_sp_block;
    LP_sp_block = NULL;
    delete[] LP_sp_index;
    LP_sp_index = NULL;
  }
}

void SparseLinearSpace::changeToDense(bool forceChange)
{
  if (SDP_sp_nBlock>0 && SDP_sp_index && SDP_sp_block) {
    for (int l=0; l<SDP_sp_nBlock; ++l) {
      SDP_sp_block[l].changeToDense(forceChange);
    }
  }
#if 0
  if (SOCP_nBlock>0 && SOCP_sp_index && SOCP_sp_block) {
    for (int l=0; l<SOCP_nBlock; ++l) {
      SOCP_sp_block[l].changeToDense(forceChange);
    }
  }
#endif

}

void SparseLinearSpace::display(FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }
  // SDP
  if (SDP_sp_nBlock>0 && SDP_sp_index && SDP_sp_block) {
    fprintf(fpout,"SDP part{\n");
    for (int l=0; l<SDP_sp_nBlock; ++l) {
      fprintf(fpout,"block %d\n",SDP_sp_index[l]);
      SDP_sp_block[l].display(fpout);
    }
    fprintf(fpout,"} \n");
  }
  // for SOCP
#if 0
  if (SOCP_sp_nBlock>0 && SOCP_sp_index && SOCP_sp_block) {
    fprintf(fpout,"SOCP part{\n");
    for (int l=0; l<SOCP_sp_nBlock; ++l) {
      fprintf(fpout,"block %d\n",SOCP_sp_index[l]);
      SOCP_sp_block[l].display(fpout);
    }
    fprintf(fpout,"} \n");
  }
#endif
  // LP
  if (LP_sp_nBlock>0 && LP_sp_index && LP_sp_block) {
    fprintf(fpout,"LP part{\n");
    for (int l=0; l<LP_sp_nBlock; ++l) {
      gmp_fprintf(fpout,"index: %d, element %Fe\n",
	      LP_sp_index[l] ,LP_sp_block[l].get_mpf_t());
    }
    fprintf(fpout,"} \n");
  }
}

bool SparseLinearSpace::copyFrom(SparseLinearSpace& other)
{
  bool total_judge;

  if (this == &other) {
    return _SUCCESS;
  }
  if (other.SDP_sp_nBlock+other.SOCP_sp_nBlock+LP_sp_nBlock<0) {
    rError("SparseLinearSpace:: nBlock is negative");
  }

  // for SDP
  if (other.SDP_sp_nBlock<0) {
    rError("SparseLinearSpace:: SDP_nBlock is negative");
  }
  if (SDP_sp_nBlock!=other.SDP_sp_nBlock && SDP_sp_index) {
    delete[] SDP_sp_index;
    SDP_sp_index = NULL;
    delete[] SDP_sp_block;
    SDP_sp_block = NULL;
  }
  SDP_sp_nBlock = other.SDP_sp_nBlock;
  if ((SDP_sp_nBlock > 0)&&(SDP_sp_index==NULL)) {
    rNewCheck();
    SDP_sp_index = new int[SDP_sp_nBlock];
    if (SDP_sp_index==NULL) {
      rError("SparseLinearSpace:: memory exhausted");
    }
    for (int l=0; l<SDP_sp_nBlock; ++l) {
      SDP_sp_index[l] = other.SDP_sp_index[l];
    }
  }
  if ((SDP_sp_nBlock > 0)&&(SDP_sp_block==NULL)) {
    rNewCheck();
    SDP_sp_block = new SparseMatrix[SDP_sp_nBlock];
    if (SDP_sp_block==NULL) {
      rError("SparseLinearSpace:: memory exhausted");
    }
  }
  total_judge = _SUCCESS;
  for (int l=0; l<SDP_sp_nBlock; ++l) {
    total_judge = SDP_sp_block[l].copyFrom(other.SDP_sp_block[l]);
  }
  if (total_judge==FAILURE) {
    rError("SparseLinearSpace:: copy miss");
  }


  // for SOCP
#if 0
  if (other.SOCP_sp_nBlock<0) {
    rError("SparseLinearSpace:: SOCP_nBlock is negative");
  }
  if (SOCP_sp_nBlock!=other.SOCP_sp_nBlock && SOCP_sp_index) {
    delete[] SOCP_sp_index;
    SOCP_sp_index = NULL;
    delete[] SOCP_sp_block;
    SOCP_sp_block = NULL;
  }
  SOCP_sp_nBlock = other.SOCP_sp_nBlock;
  if ((SOCP_sp_nBlock > 0)&&(SOCP_sp_index==NULL)) {
    rNewCheck();
    SOCP_sp_index = new int[SOCP_sp_nBlock];
    if (SOCP_sp_index==NULL) {
      rError("SparseLinearSpace:: memory exhausted");
    }
    for (int l=0; l<SOCP_sp_nBlock; ++l) {
      SOCP_sp_index[l] = other.SOCP_sp_index[l];
    }
  }
  if ((SOCP_sp_nBlock > 0)&&(SOCP_sp_block==NULL)) {
    rNewCheck();
    SOCP_sp_block = new SparseMatrix[SOCP_sp_nBlock];
    if (SOCP_sp_block==NULL) {
      rError("SparseLinearSpace:: memory exhausted");
    }
  }
  total_judge = _SUCCESS;
  for (int l=0; l<SOCP_sp_nBlock; ++l) {
    total_judge = SOCP_sp_block[l].copyFrom(other.SOCP_sp_block[l]);
  }
  if (total_judge==FAILURE) {
    rError("SparseLinearSpace:: copy miss");
  }
#endif

  // for LP
  if (other.LP_sp_nBlock<0) {
    rError("SparseLinearSpace:: LP_nBlock is negative");
  }
  if (LP_sp_nBlock!=other.LP_sp_nBlock && LP_sp_index) {
    delete[] LP_sp_index;
    LP_sp_index = NULL;
    delete[] LP_sp_block;
    LP_sp_block = NULL;
  }
  LP_sp_nBlock = other.LP_sp_nBlock;
  if ((LP_sp_nBlock > 0)&&(LP_sp_index==NULL)) {
    rNewCheck();
    LP_sp_index = new int[LP_sp_nBlock];
    if (LP_sp_index==NULL) {
      rError("SparseLinearSpace:: memory exhausted");
    }
    for (int l=0; l<LP_sp_nBlock; ++l) {
      LP_sp_index[l] = other.LP_sp_index[l];
    }
  }
  if ((LP_sp_nBlock > 0)&&(LP_sp_block==NULL)) {
    rNewCheck();
    LP_sp_block = new mpf_class[LP_sp_nBlock];
    if (LP_sp_block==NULL) {
      rError("SparseLinearSpace:: memory exhausted");
    }
  }
  total_judge = _SUCCESS;
  for (int l=0; l<LP_sp_nBlock; ++l) {
    LP_sp_block[l] = other.LP_sp_block[l];
  }
  if (total_judge==FAILURE) {
    rError("SparseLinearSpace:: copy miss");
  }


  return total_judge;
}

void SparseLinearSpace::setElement_SDP(int block, int i, int j, mpf_class ele)
{
  int k;

  // seek block
  for (k=0; k<SDP_sp_nBlock; k++){
    if (SDP_sp_index[k] == block){
      break;
    }
  }
  if (k == SDP_sp_nBlock){
    rError("SparseLinearSpace::setElement no block");
  }

  // check range
  if (SDP_sp_block[k].NonZeroCount >= SDP_sp_block[k].NonZeroNumber){
    rError("SparseLinearSpace::setElement NonZeroCount >= NonZeroNumber");
  }
  if ((i >= SDP_sp_block[k].nRow) || (j >= SDP_sp_block[k].nCol)){
    rError("out of range in input data");
  }

  // set element
  int count = SDP_sp_block[k].NonZeroCount;
  SDP_sp_block[k].row_index[count] = i;
  SDP_sp_block[k].column_index[count] = j;
  SDP_sp_block[k].sp_ele[count] = ele;
  SDP_sp_block[k].NonZeroCount++;
  if (i==j){
    SDP_sp_block[k].NonZeroEffect++;
  } else {
    SDP_sp_block[k].NonZeroEffect += 2;
  }
  
}

void SparseLinearSpace::setElement_SOCP(int block, int i, int j, mpf_class ele)
{
  rError("DenseLinearSpace:: current version does not support SOCP");
}

void SparseLinearSpace::setElement_LP(int block, mpf_class ele)
{
  int k;

  for (k=0; k<LP_sp_nBlock; k++){
    if (LP_sp_index[k] == block){
      break;
    }
  }
  if (k == LP_sp_nBlock){
    rError("SparseLinearSpace::setElement no block");
  }
  LP_sp_block[k] = ele;
}

void SparseLinearSpace::setZero()
{
  // for SDP
  if (SDP_sp_nBlock>0 && SDP_sp_index && SDP_sp_block) {
    for (int l=0; l<SDP_sp_nBlock; ++l) {
      SDP_sp_block[l].setZero();
    }
  }
  // for SOCP
#if 0
  if (SOCP_sp_nBlock>0 && SOCP_sp_index && SOCP_sp_block) {
    for (int l=0; l<SOCP_sp_nBlock; ++l) {
      SOCP_sp_block[l].setZero();
    }
  }
#endif
  // for LP
  if (LP_sp_nBlock>0 && LP_sp_index && LP_sp_block) {
    for (int l=0; l<LP_sp_nBlock; ++l) {
      LP_sp_block[l] = 0;
    }
  }
}

void SparseLinearSpace::setIdentity(mpf_class scalar)
{
  rError("SparseLinearSpace::setIdentity   no support");
  if (SDP_sp_nBlock>0 && SDP_sp_index && SDP_sp_block) {
    for (int l=0; l<SDP_sp_nBlock; ++l) {
      SDP_sp_block[l].setIdentity(scalar);
    }
  }
#if 0
  if (SOCP_sp_nBlock>0 && SOCP_sp_index && SOCP_sp_block) {
    for (int l=0; l<SOCP_sp_nBlock; ++l) {
      SOCP_sp_block[l].setIdentity(scalar);
    }
  }
  if (LP_sp_nBlock>0 && LP_sp_index && LP_sp_block) {
    for (int l=0; l<LP_sp_nBlock; ++l) {
      LP_sp_block[l].setIdentity(scalar);
    }
  }
#endif
}

bool SparseLinearSpace::sortSparseIndex(int& l, int& i, int& j)
{
  bool total_judge = _SUCCESS;
  l = -1;
  int i_in,j_in; 
  // for SDP
  if (SDP_sp_nBlock>0 && SDP_sp_index && SDP_sp_block) {
    for (int l_in=0; l_in<SDP_sp_nBlock; ++l_in) {
      total_judge = SDP_sp_block[l_in].sortSparseIndex(i_in,j_in);
      if (total_judge==FAILURE && l<0) {
	l = l_in;
	i = i_in;
	j = j_in;
      }
    }
  }
  // for SOCP
  l = -1;
  if (SOCP_sp_nBlock>0 && SOCP_sp_index && SOCP_sp_block) {
    for (int l_in=0; l_in<SOCP_sp_nBlock; ++l_in) {
      total_judge = SOCP_sp_block[l_in].sortSparseIndex(i_in,j_in);
      if (total_judge==FAILURE && l<0) {
	l = l_in;
	i = i_in;
	j = j_in;
      }
    }
  }

  return total_judge;
}


DenseLinearSpace::DenseLinearSpace()
{
  SDP_nBlock = 0;
  SDP_block = NULL;
  SOCP_nBlock = 0;
  SOCP_block = NULL;
  LP_nBlock = 0;
  LP_block = NULL;
}

DenseLinearSpace::DenseLinearSpace(int SDP_nBlock, int* SDP_blockStruct,
				   int SOCP_nBlock,  int* SOCP_blockStruct,
				   int LP_nBlock)
{
  initialize(SDP_nBlock, SDP_blockStruct,
	     SOCP_nBlock, SOCP_blockStruct,
	     LP_nBlock);
}

DenseLinearSpace::~DenseLinearSpace()
{
  terminate();
}

void DenseLinearSpace::initialize(int SDP_nBlock, int* SDP_blockStruct,
				  int SOCP_nBlock,  int* SOCP_blockStruct,
				  int LP_nBlock)
{
  // rMessage("DenseLinearSpace::initialize");
  if (SDP_nBlock + SOCP_nBlock + LP_nBlock <= 0) {
    rError("DenseLinearSpace:: SDP + SOCP + LP Block is nonpositive");
  }

  // for SDP
  this->SDP_nBlock = SDP_nBlock;
  if (SDP_nBlock<0) {
    rError("DenseLinearSpace:: SDP_nBlock is negative");
  }
  if ((SDP_nBlock > 0) && (SDP_block==NULL)) {
    rNewCheck();
    SDP_block = new DenseMatrix[SDP_nBlock];
    if (SDP_block==NULL) {
      rError("DenseLinearSpace:: memory exhausted");
    }
  }
  for (int l=0; l<SDP_nBlock; ++l) {
    int size = SDP_blockStruct[l];
    if (size>0) {
      SDP_block[l].initialize(size,size,DenseMatrix::DENSE);
    } else {
      rError("DenseLinearSpace:: SDP size is nonpositive");
    }
  }

  // for SOCP
  this->SOCP_nBlock = 0;
  this->SOCP_block = NULL;
#if 0
  this->SOCP_nBlock = SOCP_nBlock;
  if (SOCP_nBlock<0) {
    rError("DenseLinearSpace:: SOCP_nBlock is negative");
  }
  if ((SOCP_nBlock > 0) && (SOCP_block==NULL)) {
    rNewCheck();
    SOCP_block = new DenseMatrix[SOCP_nBlock];
    if (SOCP_block==NULL) {
      rError("DenseLinearSpace:: memory exhausted");
    }
  }
  for (int l=0; l<SOCP_nBlock; ++l) {
    int size = SOCP_blockStruct[l];
    if (size>0) {
      SOCP_block[l].initialize(size,size,DenseMatrix::DENSE);
    } else {
      rError("DenseLinearSpace:: SOCP size is nonpositive");
    }
  }
#endif

  // for LP
  this->LP_nBlock = LP_nBlock;
  if (LP_nBlock<0) {
    rError("DenseLinearSpace:: LP_nBlock is negative");
  }
  if ((LP_nBlock > 0) && (LP_block==NULL)) {
    rNewCheck();
    LP_block = new mpf_class[LP_nBlock];
    if (LP_block==NULL) {
      rError("DenseLinearSpace:: memory exhausted");
    }
  }
  for (int l=0; l<LP_nBlock; ++l) {
    LP_block[l] = 0.0;
  }

}

void DenseLinearSpace::terminate()
{
  // for SDP
  if (SDP_block && SDP_nBlock>0) {
    for (int l=0; l<SDP_nBlock; ++l) {
      SDP_block[l].terminate();
    }
    delete[] SDP_block;
    SDP_block = NULL;
  }

  // SOCP
#if 0
  if (SOCP_block && SOCP_nBlock>0) {
    for (int l=0; l<SOCP_nBlock; ++l) {
      SOCP_block[l].terminate();
    }
    delete[] SOCP_block;
    SOCP_block = NULL;
  }
#endif

  // LP
  if (LP_block && LP_nBlock>0) {
    delete[] LP_block;
    LP_block = NULL;
  }

}

void DenseLinearSpace::display(FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }

  if (SDP_nBlock>0 && SDP_block) {
    fprintf(fpout,"SDP part{\n");
    for (int l=0; l<SDP_nBlock; ++l) {
      SDP_block[l].display(fpout);
    }
    fprintf(fpout,"} \n");
  }

#if 0
  if (SOCP_nBlock>0 && SOCP_block) {
    fprintf(fpout,"SOCP part{\n");
    for (int l=0; l<SOCP_nBlock; ++l) {
      SOCP_block[l].display(fpout);
    }
    fprintf(fpout,"} \n");
  }
#endif

  if (LP_nBlock>0 && LP_block) {
    fprintf(fpout,"LP part{\n");
    for (int l=0; l<LP_nBlock; ++l) {
      gmp_fprintf(fpout,"%Fe, ",LP_block[l].get_mpf_t());
    }
    fprintf(fpout,"} \n");
  }
}

bool DenseLinearSpace::copyFrom(DenseLinearSpace& other)
{
  if (this == &other) {
    return _SUCCESS;
  }

  if (other.SDP_nBlock+other.SOCP_nBlock+other.LP_nBlock<=0) {
    rError("DenseLinearSpace:: SDP + SOCP + LP Block is nonpositive");
  }
  bool total_judge = _SUCCESS;

  // for SDP
  if (other.SDP_nBlock<0) {
    rError("DenseLinearSpace:: SDP_nBlock is negative");
  }
  if (SDP_nBlock!=other.SDP_nBlock) {
    delete[] SDP_block;
    SDP_block = NULL;
  }
  SDP_nBlock = other.SDP_nBlock;
  if ((SDP_nBlock > 0) && (SDP_block == NULL)) {
    rNewCheck();
    SDP_block = new DenseMatrix[SDP_nBlock];
    if (SDP_block==NULL) {
      rError("DenseLinearSpace:: memory exhausted");
    }
  }
  for (int l=0; l<SDP_nBlock; ++l) {
    total_judge = SDP_block[l].copyFrom(other.SDP_block[l]);
  }
  if (total_judge==FAILURE) {
    rError("DenseLinearSpace:: copy miss");
  }

  
  // for SOCP
#if 0
  if (other.SOCP_nBlock<0) {
    rError("DenseLinearSpace:: SOCP_nBlock is negative");
  }
  if (SOCP_nBlock!=other.SOCP_nBlock) {
    delete[] SOCP_block;
    SOCP_block = NULL;
  }
  SOCP_nBlock = other.SOCP_nBlock;
  if ((SOCP_nBlock > 0) && (SOCP_block == NULL)) {
    rNewCheck();
    SOCP_block = new DenseMatrix[SOCP_nBlock];
    if (SOCP_block==NULL) {
      rError("DenseLinearSpace:: memory exhausted");
    }
  }
  for (int l=0; l<SOCP_nBlock; ++l) {
    total_judge = SOCP_block[l].copyFrom(other.SOCP_block[l]);
  }
  if (total_judge==FAILURE) {
    rError("DenseLinearSpace:: copy miss");
  }
#endif

  // for LP
  if (other.LP_nBlock<0) {
    rError("DenseLinearSpace:: LP_nBlock is negative");
  }
  if (LP_nBlock!=other.LP_nBlock) {
    delete[] LP_block;
    LP_block = NULL;
  }
  LP_nBlock = other.LP_nBlock;
  if ((LP_nBlock > 0) && (LP_block == NULL)) {
    LP_block = new mpf_class[LP_nBlock];
    if (LP_block==NULL) {
      rError("DenseLinearSpace:: memory exhausted");
    }
  }
  for (int l=0; l<LP_nBlock; ++l) {
    LP_block[l] = other.LP_block[l];
  }


  return total_judge;
}

void DenseLinearSpace::setElement_SDP(int block, int i, int j, mpf_class ele)
{

  // check range
  if (block >= SDP_nBlock){
    rError("out of range in input data");
  }
  if ((i >= SDP_block[block].nRow) || (j >= SDP_block[block].nCol)){
    rError("out of range in input data");
  }

  int nCol = SDP_block[block].nCol;
  SDP_block[block].de_ele[i + j * nCol] = ele;
  SDP_block[block].de_ele[j + i * nCol] = ele;
}

void DenseLinearSpace::setElement_SOCP(int block, int i, int j, mpf_class ele)
{
  rError("DenseLinearSpace:: current version does not support SOCP");
}

void DenseLinearSpace::setElement_LP(int block, mpf_class ele)
{
  // check range
  if (block >= LP_nBlock){
    rError("out of range in input data");
  }
  LP_block[block] = ele;
}

void DenseLinearSpace::setZero()
{
  // for SDP
  if (SDP_nBlock>0 && SDP_block) {
    for (int l=0; l<SDP_nBlock; ++l) {
      SDP_block[l].setZero();
    }
  }

  // for SOCP
#if 0  
  if (SOCP_nBlock>0 && SOCP_block) {
    for (int l=0; l<SOCP_nBlock; ++l) {
      SOCP_block[l].setZero();
    }
  }
#endif

  // for LP
  if (LP_nBlock>0 && LP_block) {
    for (int l=0; l<LP_nBlock; ++l) {
      LP_block[l] = 0.0;
    }
  }

}

void DenseLinearSpace::setIdentity(mpf_class scalar)
{
  // for SDP
  if (SDP_nBlock>0 && SDP_block) {
    for (int l=0; l<SDP_nBlock; ++l) {
      SDP_block[l].setIdentity(scalar);
    }
  }
  // for SOCP
#if 0
  if (SOCP_nBlock>0 && SOCP_block) {
    for (int l=0; l<SOCP_nBlock; ++l) {
      SOCP_block[l].setIdentity(scalar);
    }
  }
#endif
  // for LP
  if (LP_nBlock>0 && LP_block) {
    for (int l=0; l<LP_nBlock; ++l) {
      LP_block[l] = scalar;
    }
  }


}


//print vector as matlab compat format.
void gmp_printvec (int N, mpf_class *A)
{
  printf(" [ ");
  for (int i=0; i<N-1; i++){
    gmp_printf( P_FORMAT "; ", A[i].get_mpf_t());
  }
  gmp_printf( P_FORMAT " ] ", A[N-1].get_mpf_t());
}

//print vector as matlab compat format.
void gmp_printveci (int N, mpf_class *A, int inc)
{
  int ix=0;
  printf(" [ ");
  for (int i=0; i<N-1; i++){
    gmp_printf( P_FORMAT "; ", A[ix].get_mpf_t());
    ix = ix + inc;
  }
  gmp_printf( P_FORMAT " ] ", A[ix].get_mpf_t());
}

//print matrix as matlab compat format.
void gmp_printmat (int N, int M, mpf_class *A, int LDA)
{
  mpf_class mtmp;
  printf("[ ");
  for (int i=0; i<N; i++){
    printf("[ ");
    for (int j=0; j<M; j++){
      mtmp=A[i+j*LDA];
      gmp_printf( P_FORMAT, mtmp.get_mpf_t());
      if (j<M-1) printf(", ");
    }
    if (i<N-1) printf("]; "); else printf ("] ");
  }
  printf("]");
}

//print symmetry matrix as matlab compat format.
void gmp_printsymmmat (int N, mpf_class *A, int LDA)
{
  mpf_class mtmp;
  printf("[ ");
  for (int i=0; i<N; i++){
    printf("[ ");
    for (int j=0; j<N; j++){
      mtmp=A[i+j*LDA];
      gmp_printf(P_FORMAT, mtmp.get_mpf_t());
      if (j<N-1) printf(", ");
    }
    if (i<N-1) printf("]; "); else printf ("] ");
  }
  printf("]");
}

void gmp_printmat3 (int N, mpf_class *p, mpf_class *q)
{
  mpf_class mtmp;
  mpf_class zero;
  zero = 0.0;
  printf(" [ ");
  for (int i=0; i<N; i++){
    printf(" [ ");
    for (int j=0; j<N; j++){
      if (i==j){
        gmp_printf(P_FORMAT, p[i].get_mpf_t());
      }
      if(abs(i-j)==1){
	gmp_printf(P_FORMAT, q[(i<j?i:j)].get_mpf_t());
      }
      if(abs(i-j)>1){
        gmp_printf(P_FORMAT, zero.get_mpf_t());
      }
      if (j<N-1)  printf (", ");
    }
    if (i<N-1) printf("]; "); else printf ("] ");
  }
  printf(" ] ");
}

}
