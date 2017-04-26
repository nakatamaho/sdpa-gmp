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

#define DIMACS_PRINT 0

#include <sdpa_io.h>
#include <vector>
#include <algorithm>

namespace sdpa {

// 2008/02/27  kazuhide nakata 
#if 0  // not use
void IO::read(FILE* fpData, int m,
	      int SDP_nBlock,int* SDP_blockStruct,
	      int SOCP_nBlock,int* SOCP_blockStruct,
	      int LP_nBlock, 
		  int nBlock, int* blockStruct, int* blockType, int* blockNumber,
	      InputData& inputData, bool isDataSparse)
{
  inputData.initialize_bVec(m);
  read(fpData,inputData.b);
  long position = ftell(fpData);
  // C,A must be accessed "twice".

  // count numbers of elements of C and A
  int* SDP_CNonZeroCount = NULL;
  SDP_CNonZeroCount = new int[SDP_nBlock];
  if (SDP_CNonZeroCount==NULL) {
    rError("Memory exhausted about blockStruct");
  }
  int* SDP_ANonZeroCount = NULL;
  SDP_ANonZeroCount = new int[m*SDP_nBlock];
  if (SDP_ANonZeroCount==NULL) {
    rError("Memory exhausted about blockStruct");
  }

  // count numbers of elements of C and A
  int* SOCP_CNonZeroCount = NULL;
  SOCP_CNonZeroCount = new int[SOCP_nBlock];
  if (SOCP_CNonZeroCount==NULL) {
    rError("Memory exhausted about blockStruct");
  }
  int* SOCP_ANonZeroCount = NULL;
  SOCP_ANonZeroCount = new int[m*SOCP_nBlock];
  if (SOCP_ANonZeroCount==NULL) {
    rError("Memory exhausted about blockStruct");
  }

  // count numbers of elements of C and A
  bool* LP_CNonZeroCount = NULL;
  LP_CNonZeroCount = new bool[LP_nBlock];
  if (LP_CNonZeroCount==NULL) {
    rError("Memory exhausted about blockStruct");
  }
  bool* LP_ANonZeroCount = NULL;
  LP_ANonZeroCount = new bool[m*LP_nBlock];
  if (LP_ANonZeroCount==NULL) {
    rError("Memory exhausted about blockStruct");
  }

  //   initialize C and A
  read(fpData,m,
       SDP_nBlock, SDP_blockStruct, SDP_CNonZeroCount, SDP_ANonZeroCount,
       SOCP_nBlock, SOCP_blockStruct, SOCP_CNonZeroCount, SOCP_ANonZeroCount,
       LP_nBlock, LP_CNonZeroCount, LP_ANonZeroCount,
	   nBlock, blockStruct, blockType, blockNumber,
       isDataSparse);
  //   rMessage(" C and A count over");
  inputData.initialize_CMat(SDP_nBlock, SDP_blockStruct,
			    SDP_CNonZeroCount,
			    SOCP_nBlock,  SOCP_blockStruct,
			    SOCP_CNonZeroCount,
			    LP_nBlock, LP_CNonZeroCount);
  inputData.initialize_AMat(m,SDP_nBlock, SDP_blockStruct,
			    SDP_ANonZeroCount,
			    SOCP_nBlock,  SOCP_blockStruct,
			    SOCP_ANonZeroCount,
			    LP_nBlock, LP_ANonZeroCount);
  delete[] SDP_CNonZeroCount;
  SDP_CNonZeroCount = NULL;
  delete[] SDP_ANonZeroCount;
  SDP_ANonZeroCount = NULL;
  delete[] SOCP_CNonZeroCount;
  SOCP_CNonZeroCount = NULL;
  delete[] SOCP_ANonZeroCount;
  SOCP_ANonZeroCount = NULL;
  delete[] LP_CNonZeroCount;
  LP_CNonZeroCount = NULL;
  delete[] LP_ANonZeroCount;
  LP_ANonZeroCount = NULL;
    
  //   rMessage(" C and A initialize over");
  read(fpData, inputData, m, 
       SDP_nBlock, SDP_blockStruct, 
       SOCP_nBlock, SOCP_blockStruct, 
       LP_nBlock, 
	   nBlock, blockStruct, blockType, blockNumber,
       position, isDataSparse);
  //   rMessage(" C and A have been read");
}
#endif

void IO::read(FILE* fpData, FILE* fpout, int& m, char* str)
{
  while (true) {
  volatile int dummy=0; //for gcc-3.3 bug
    fgets(str,lengthOfString,fpData);
    if (str[0]=='*' || str[0]=='"') {
      fprintf(fpout,"%s",str);
    } else {
      sscanf(str,"%d",&m);
      break;
    }
  }
}

void IO::read(FILE* fpData, int& nBlock)
{
  fscanf(fpData,"%d",&nBlock);
}

void IO::read(FILE* fpData, 
			  int nBlock, int* blockStruct)
{
  for (int l=0; l<nBlock; ++l) {
    fscanf(fpData,"%*[^0-9+-]%d",&blockStruct[l]);
  }
}

void IO::read(FILE* fpData, Vector& b)
{
  for (int k=0; k<b.nDim; ++k) {
    gmp_fscanf(fpData,"%*[^0-9+-]%Fe",&b.ele[k]);
  }
}

void IO::read(FILE* fpData, DenseLinearSpace& xMat,
	      Vector& yVec, DenseLinearSpace& zMat,
	      bool inputSparse)
{
  // read initial point 
  int SDP_nBlock = xMat.SDP_nBlock;
  int SOCP_nBlock = xMat.SOCP_nBlock;
  int LP_nBlock = xMat.LP_nBlock;

  // yVec is opposite sign
  for (int k=0; k<yVec.nDim; ++k) {
    mpf_class tmp;
    gmp_fscanf(fpData,"%*[^0-9+-]%Fe",&tmp);
    yVec.ele[k] = -tmp;
    //     rMessage("yVec.ele[" << k << "] = " << tmp);
  }
  
  if (inputSparse) {
    // sparse case , zMat , xMat in this order
    int i,j,l,target;
    mpf_class value;
    while (true) {
      if (fscanf(fpData,"%*[^0-9+-]%d",&target)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&l)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&i)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&j)<=0) {
	break;
      }
      if (gmp_fscanf(fpData,"%*[^0-9+-]%Fe",&value)<=0) {
	break;
      }
      #if 0
      rMessage("target = " << target
	       << ": l " << l
	       << ": i " << i
	       << ": j " << j
	       << ": value " <<value);
      #endif

      if (l <= SDP_nBlock){
	// SDP part
	if (target==1) {
	  zMat.setElement_SDP(l-1,i-1,j-1,value);
	} else {
	  xMat.setElement_SDP(l-1,i-1,j-1,value);
	}
      } else if (l <= SDP_nBlock + SOCP_nBlock){
	// SOCP part
	rError("io:: current version does not support SOCP");
	int ll = l - SDP_nBlock;
	if (target==1) {
	  zMat.setElement_SOCP(ll-1,i-1,j-1,value);
	} else {
	  xMat.setElement_SOCP(ll-1,i-1,j-1,value);
	}
      } else {
	// LP part
	if (i != j){
	  rError("io:: LP part  3rd elemtn != 4th elemnt");
	}
	if (target==1) {
	  zMat.setElement_LP(i-1,value);
	} else {
	  xMat.setElement_LP(i-1,value);
	}
      }
    } // end of 'while (true)'
  } else {
    // dense case , zMat , xMat in this order
    // for SDP
    for (int l=0; l<SDP_nBlock; ++l) {
      int size = zMat.SDP_block[l].nRow;
      for (int i=0; i<size; ++i) {
	for (int j=0; j<size; ++j) {
	  mpf_class tmp;
	  gmp_fscanf(fpData,"%*[^0-9+-]%Fe",&tmp);
	  if (i<=j && tmp!=0.0) {
	    zMat.setElement_SDP(l,i,j,tmp);
	  }
	}
      }
    }
    // for SOCP
    for (int l=0; l<SOCP_nBlock; ++l) {
	rError("io:: current version does not support SOCP");
    }
    // for LP
    for (int j=0; j<LP_nBlock; ++j) {
      mpf_class tmp;
      gmp_fscanf(fpData,"%*[^0-9+-]%Fe",&tmp);
      if (tmp!=0.0) {
	zMat.setElement_LP(j,tmp);
      }
    }

    // for SDP
    for (int l=0; l<SDP_nBlock; ++l) {
      int size = xMat.SDP_block[l].nRow;
      for (int i=0; i<size; ++i) {
	for (int j=0; j<size; ++j) {
	  mpf_class tmp;
	  gmp_fscanf(fpData,"%*[^0-9+-]%Fe",&tmp);
	  if (i<=j && tmp!=0.0) {
	    xMat.setElement_SDP(l,i,j,tmp);
	  }
	}
      }
    }
    // for SOCP
    for (int l=0; l<SOCP_nBlock; ++l) {
	rError("io:: current version does not support SOCP");
    }
    // for LP
    for (int j=0; j<LP_nBlock; ++j) {
      mpf_class tmp;
      gmp_fscanf(fpData,"%*[^0-9+-]%Fe",&tmp);
      if (tmp!=0.0) {
	xMat.setElement_LP(j,tmp);
      }
    }
  } // end of 'if (inputSparse)'
}

// 2008/02/27 kazuhide nakata   
// not use
void IO::read(FILE* fpData, int m, 
	      int SDP_nBlock, int* SDP_blockStruct,
	      int* SDP_CNonZeroCount, int* SDP_ANonZeroCount,
	      int SOCP_nBlock, int* SOCP_blockStruct,
	      int* SOCP_CNonZeroCount, int* SOCP_ANonZeroCount,
	      int LP_nBlock,
	      bool* LP_CNonZeroCount, bool* LP_ANonZeroCount,
		  int nBlock, int* blockStruct, int* blockType, int* blockNumber,
	      bool isDataSparse)
{
  // only count the numbers of C,A[k]
  for (int l=0; l<SDP_nBlock; ++l) {
    SDP_CNonZeroCount[l] = 0;
  }
  for (int k=0; k<m; ++k) {
    for (int l=0; l<SDP_nBlock; ++l) {
      SDP_ANonZeroCount[k*SDP_nBlock + l] = 0;
    }
  }
  for (int l=0; l<SOCP_nBlock; ++l) {
    SOCP_CNonZeroCount[l] = 0;
  }
  for (int k=0; k<m; ++k) {
    for (int l=0; l<SOCP_nBlock; ++l) {
      SOCP_ANonZeroCount[k*SOCP_nBlock + l] = 0;
    }
  }
  for (int l=0; l<LP_nBlock; ++l) {
    LP_CNonZeroCount[l] = false;
  }
  for (int k=0; k<m; ++k) {
    for (int l=0; l<LP_nBlock; ++l) {
      LP_ANonZeroCount[k*LP_nBlock + l] = false;
    }
  }
  if (isDataSparse) {
    int i,j,k,l;
    mpf_class value;
    while (true) {
      if (fscanf(fpData,"%*[^0-9+-]%d",&k)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&l)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&i)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&j)<=0) {
	break;
      }
      if (gmp_fscanf(fpData,"%*[^0-9+-]%Fe",&value)<=0) {
	break;
      }

      if (blockType[l-1] == 1){	// SDP part
		int l2 = blockNumber[l-1];
		if (k==0) {
		  SDP_CNonZeroCount[l2]++;
		} else {
		  SDP_ANonZeroCount[(k-1)*SDP_nBlock+l2]++;
		}
      } else if (blockType[l-1] == 2){	// SOCP part
		rError("io:: current version does not support SOCP");
		int l2 = blockNumber[l-1];;
		if (k==0) {
		  SOCP_CNonZeroCount[l2]++;
		} else {
		  SOCP_ANonZeroCount[(k-1)*SOCP_nBlock+l2]++;
		}
      } else if (blockType[l-1] == 3){ // LP part
		int l2 =blockNumber[l-1];
		if (k==0) {
		  LP_CNonZeroCount[l2+i-1] = true;
		} else {
		  LP_ANonZeroCount[(k-1)*LP_nBlock+l2+i-1] = true;
		}
      } else {
		rError("io::read not valid blockType");
	  }
    }// end of 'while (true)'

  } else { // isDataSparse == false

    // k==0
	for (int l2=0; l2<nBlock; ++l2){
	  if (blockType[l2] == 1) { // SDP part
		int l = blockNumber[l2];
		int size = SDP_blockStruct[l];
		for (int i=0; i<size; ++i) {
		  for (int j=0; j<size; ++j) {
			mpf_class tmp;
			gmp_fscanf(fpData,"%*[^0-9+-]%Fe",&tmp);
			if (i<=j && tmp!=0.0) {
			  SDP_CNonZeroCount[l]++;
			}
		  }
		}
	  } else if (blockType[l2] == 2) { // SOCP part
		  rError("io:: current version does not support SOCP");
	  } else if (blockType[l2] == 3) { // LP part
		for (int j=0; j<blockStruct[l2]; ++j) {
		  mpf_class tmp;
		  gmp_fscanf(fpData,"%*[^0-9+-]%Fe",&tmp);
		  if (tmp!=0.0) {
			LP_CNonZeroCount[blockNumber[l2]+j]++;
		  }
		}
	  } else {
		rError("io::read not valid blockType");
	  }
	}

    for (int k=0; k<m; ++k) {
    // k>0
	for (int l2=0; l2<nBlock; ++l2){
	  if (blockType[l2] == 1) { // SDP part
		int l = blockNumber[l2];
		int size = SDP_blockStruct[l];
		for (int i=0; i<size; ++i) {
		  for (int j=0; j<size; ++j) {
			mpf_class tmp;
			gmp_fscanf(fpData,"%*[^0-9+-]%Fe",&tmp);
			if (i<=j && tmp!=0.0) {
			  SDP_ANonZeroCount[k*SDP_nBlock+l]++;
			}
		  }
		}
	  } else if (blockType[l2] == 2) { // SOCP part
		  rError("io:: current version does not support SOCP");
	  } else if (blockType[l2] == 3) { // LP part
		for (int j=0; j<blockStruct[l2]; ++j) {
		  mpf_class tmp;
		  gmp_fscanf(fpData,"%*[^0-9+-]%Fe",&tmp);
		  if (tmp!=0.0) {
			LP_ANonZeroCount[k*LP_nBlock+blockNumber[l2]+j] = true;
		  }
		}
	  } else {
		rError("io::read not valid blockType");
	  }
	}
	}

  } // end of 'if (isDataSparse)'

}


// 2008/02/27 kazuhide nakata   
// not use
void IO::read(FILE* fpData,  InputData& inputData, int m, 
	      int SDP_nBlock, int* SDP_blockStruct, 
	      int SOCP_nBlock, int* SOCP_blockStruct, 
	      int LP_nBlock, 
		  int nBlock, int* blockStruct, int* blockType, int* blockNumber,
	      long position, bool isDataSparse)

{
  // in Sparse, read C,A[k]

  // seed the positon of C in the fpData
  fseek(fpData, position, 0);

  if (isDataSparse) {
    int i,j,k,l;
    mpf_class value;
    while (true) {
      if (fscanf(fpData,"%*[^0-9+-]%d",&k)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&l)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&i)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&j)<=0) {
	break;
      }
      if (gmp_fscanf(fpData,"%*[^0-9+-]%Fe",&value)<=0) {
	break;
      }
#if 0
      rMessage("input k:" << k <<
	       " l:" << l <<
	       " i:" << i <<
	       " j:" << j);
#endif     

      if (blockType[l-1] == 1){	// SDP part
		int l2 = blockNumber[l-1];
		if (k==0) {
		  inputData.C.setElement_SDP(l2,i-1,j-1,-value);
		} else {
		  inputData.A[k-1].setElement_SDP(l2,i-1,j-1,value);
		}
      } else if (blockType[l-1] == 2){	// SOCP part
		rError("io:: current version does not support SOCP");
		int l2 = blockNumber[l-1];
		if (k==0) {
		  inputData.C.setElement_SOCP(l2,i-1,j-1,-value);
		} else {
		  inputData.A[k-1].setElement_SOCP(l2,i-1,j-1,value);
		}
      } else if (blockType[l-1] == 3){ // LP part
		if (i != j){
		  rError("io:: LP part  3rd elemtn != 4th elemnt");
		}
		if (k==0) {
		  inputData.C.setElement_LP(blockNumber[l-1]+i-1,-value);
		} else {
		  inputData.A[k-1].setElement_LP(blockNumber[l-1]+i-1,value);
		}
      } else {
		rError("io::read not valid blockType");
	  }
	} 
  } else {  // dense

    // k==0
	for (int l2=0; l2<nBlock; ++l2){
	  if (blockType[l2] == 1) { // SDP part
		int l = blockNumber[l2];
		int size = SDP_blockStruct[l];
		for (int i=0; i<size; ++i) {
		  for (int j=0; j<size; ++j) {
			mpf_class tmp;
			gmp_fscanf(fpData,"%*[^0-9+-]%Fe",&tmp);
			if (i<=j && tmp!=0.0) {
			  inputData.C.setElement_SDP(l,i,j,-tmp);
			}
		  }
		}
	  } else if (blockType[l2] == 2) { // SOCP part
		rError("io:: current version does not support SOCP");
	  } else if (blockType[l2] == 3) { // LP part
		for (int j=0; j<blockStruct[l2]; ++j) {
		  mpf_class tmp;
		  gmp_fscanf(fpData,"%*[^0-9+-]%Fe",&tmp);
		  if (tmp!=0.0) {
			inputData.C.setElement_LP(blockNumber[l2]+j,-tmp);
		  }
		}
	  } else {
		rError("io::read not valid blockType");
	  }
	}

	// k > 0
    for (int k=0; k<m; ++k) {
	  
	for (int l2=0; l2<nBlock; ++l2){
	  if (blockType[l2] == 1) { // SDP part
		int l = blockNumber[l2];
		int size = SDP_blockStruct[l];
		for (int i=0; i<size; ++i) {
		  for (int j=0; j<size; ++j) {
			mpf_class tmp;
			gmp_fscanf(fpData,"%*[^0-9+-]%Fe",&tmp);
			if (i<=j && tmp!=0.0) {
			  inputData.A[k].setElement_SDP(l,i,j,tmp);
			}
		  }
		}
	  } else if (blockType[l2] == 2) { // SOCP part
		rError("io:: current version does not support SOCP");
	  } else if (blockType[l2] == 3) { // LP part
		for (int j=0; j<blockStruct[l2]; ++j) {
		  mpf_class tmp;
		  gmp_fscanf(fpData,"%*[^0-9+-]%Fe",&tmp);
		  if (tmp!=0.0) {
			inputData.A[k].setElement_LP(blockNumber[l2]+j,tmp);
		  }
		}
	  } else {
		rError("io::read not valid blockType");
	  }
	}

	} // for k

  } // end of 'if (isDataSparse)'

}

// 2008/02/27 kazuhide nakata
// without LP_ANonZeroCount
#if 1
void IO::read(FILE* fpData, int m,
	      int SDP_nBlock,
              int* SDP_blockStruct,
	      int SOCP_nBlock,
              int* SOCP_blockStruct,
	      int LP_nBlock, 
              int nBlock, 
              int* blockStruct,
              int* blockType, 
              int* blockNumber,
	      InputData& inputData, bool isDataSparse)
{
  inputData.initialize_bVec(m);
  read(fpData,inputData.b);
  long position = ftell(fpData);

  // C,A must be accessed "mpf_class".

  //   initialize block struct of C and A
  setBlockStruct(fpData, inputData, m,
                 SDP_nBlock, 
                 SDP_blockStruct,
                 SOCP_nBlock, 
                 SOCP_blockStruct,
                 LP_nBlock,
                 nBlock, blockStruct, blockType, blockNumber,
                 position, isDataSparse);
  //   rMessage(" C and A initialize over");
    
  setElement(fpData, inputData, m, 
             SDP_nBlock, 
             SDP_blockStruct, 
             SOCP_nBlock, 
             SOCP_blockStruct, 
             LP_nBlock, 
             nBlock, blockStruct, blockType, blockNumber,
             position, isDataSparse);
  //   rMessage(" C and A have been read");
}
#endif

  // 2008/02/27 kazuhide nakata   
  // without LP_ANonZeroCount
void IO::setBlockStruct(FILE* fpData, InputData& inputData, int m,
                        int SDP_nBlock, 
                        int* SDP_blockStruct,
                        int SOCP_nBlock, 
                        int* SOCP_blockStruct,
                        int LP_nBlock,
                        int nBlock, int* blockStruct, 
                        int* blockType, int* blockNumber,
                        long position, bool isDataSparse)
{
  // seed the positon of C in the fpData
  fseek(fpData, position, 0);

  vector<int>* SDP_index = NULL;
  SDP_index = new vector<int>[m+1];
  vector<int>* SOCP_index = NULL;
  SOCP_index = new vector<int>[m+1];
  vector<int>* LP_index = NULL;
  LP_index = new vector<int>[m+1];

  // for SDP
  int SDP_sp_nBlock;
  int* SDP_sp_index = NULL;
  int* SDP_sp_blockStruct = NULL;
  int* SDP_sp_NonZeroNumber = NULL;
  SDP_sp_index = new int[SDP_nBlock];
  SDP_sp_blockStruct = new int[SDP_nBlock];
  SDP_sp_NonZeroNumber = new int[SDP_nBlock];
  // for SOCP
  int SOCP_sp_nBlock;
  int* SOCP_sp_blockStruct;
  int* SOCP_sp_index;
  int* SOCP_sp_NonZeroNumber;
  // for LP
  int LP_sp_nBlock;
  int* LP_sp_index;
  LP_sp_index = new int[LP_nBlock];


  if (isDataSparse) {
    int i,j,k,l;
    mpf_class value;
    while (true) {
      if (fscanf(fpData,"%*[^0-9+-]%d",&k)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&l)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&i)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&j)<=0) {
	break;
      }
      if (gmp_fscanf(fpData,"%*[^0-9+-]%Fe",&value)<=0) {
	break;
      }
      
      if (blockType[l-1] == 1){	// SDP part
        int l2 = blockNumber[l-1];
        SDP_index[k].push_back(l2);
      } else if (blockType[l-1] == 2){	// SOCP part
        rError("io:: current version does not support SOCP");
        int l2 = blockNumber[l-1];;
        SOCP_index[k].push_back(l2);
      } else if (blockType[l-1] == 3){ // LP part
        if (i!=j){
          printf("invalid data file k:%d, l:%d, i:%d, j:%d, "
                 ,k,l,i,j);
          printf("value:%9.1Felf\n",
                 value.get_mpf_t());
          rError("IO::initializeLinearSpace");
        }
        int l2 =blockNumber[l-1];
        LP_index[k].push_back(l2+i-1);
      } else {
        rError("io::read not valid blockType");
      }
    }// end of 'while (true)'
    
  } else { // isDataSparse == false
    
    // k==0
    for (int l2=0; l2<nBlock; ++l2){
      if (blockType[l2] == 1) { // SDP part
        int l = blockNumber[l2];
        int size = SDP_blockStruct[l];
        for (int i=0; i<size; ++i) {
          for (int j=0; j<size; ++j) {
            mpf_class tmp;
            gmp_fscanf(fpData,"%*[^0-9+-]%Fe",&tmp);
            if (i<=j && tmp!=0.0) {
              SDP_index[0].push_back(l);
            }
          }
        }
      } else if (blockType[l2] == 2) { // SOCP part
        rError("io:: current version does not support SOCP");
      } else if (blockType[l2] == 3) { // LP part
        for (int j=0; j<blockStruct[l2]; ++j) {
          mpf_class tmp;
          gmp_fscanf(fpData,"%*[^0-9+-]%Fe",&tmp);
          if (tmp!=0.0) {
              LP_index[0].push_back(blockNumber[l2]+j);
          }
        }
      } else {
        rError("io::read not valid blockType");
      }
    }
    
    for (int k=0; k<m; ++k) {
      // k>0
      for (int l2=0; l2<nBlock; ++l2){
        if (blockType[l2] == 1) { // SDP part
          int l = blockNumber[l2];
          int size = SDP_blockStruct[l];
          for (int i=0; i<size; ++i) {
            for (int j=0; j<size; ++j) {
              mpf_class tmp;
              gmp_fscanf(fpData,"%*[^0-9+-]%Fe",&tmp);
              if (i<=j && tmp!=0.0) {
                SDP_index[k+1].push_back(l);
              }
            }
          }
        } else if (blockType[l2] == 2) { // SOCP part
          rError("io:: current version does not support SOCP");
        } else if (blockType[l2] == 3) { // LP part
          for (int j=0; j<blockStruct[l2]; ++j) {
            mpf_class tmp;
            gmp_fscanf(fpData,"%*[^0-9+-]%Fe",&tmp);
            if (tmp!=0.0) {
              LP_index[k+1].push_back(blockNumber[l2]+j);
            }
          }
        } else {
          rError("io::read not valid blockType");
        }
      }
    }
    
  } // end of 'if (isDataSparse)'
  

  inputData.A = new SparseLinearSpace[m];
  for (int k=0 ; k<m+1; k++){
    sort(SDP_index[k].begin(),SDP_index[k].end());
    SDP_sp_nBlock = 0;
    int previous_index = -1;
    int index;
    for (int i=0; i<SDP_index[k].size(); i++){
      index = SDP_index[k][i];
      if (previous_index != index){
        SDP_sp_index[SDP_sp_nBlock] = index;
        SDP_sp_blockStruct[SDP_sp_nBlock] = SDP_blockStruct[index];
        SDP_sp_NonZeroNumber[SDP_sp_nBlock] = 1;
        previous_index = index;
        SDP_sp_nBlock++;
      } else {
        SDP_sp_NonZeroNumber[SDP_sp_nBlock-1]++;
      }
    }
    sort(LP_index[k].begin(),LP_index[k].end());
    LP_sp_nBlock=0;
    previous_index = -1;
    for (int i=0; i<LP_index[k].size(); i++){
      index = LP_index[k][i];
      if (previous_index != index){
        LP_sp_index[LP_sp_nBlock] = index;
        previous_index = index;
        LP_sp_nBlock++;
      }
    }

    if (k==0){
      inputData.C.initialize(SDP_sp_nBlock, 
                             SDP_sp_index,
                             SDP_sp_blockStruct, 
                             SDP_sp_NonZeroNumber,
                             SOCP_sp_nBlock, 
                             SOCP_sp_blockStruct, 
                             SOCP_sp_index,
                             SOCP_sp_NonZeroNumber,
                             LP_sp_nBlock, 
                             LP_sp_index);
    } else {
      inputData.A[k-1].initialize(SDP_sp_nBlock, 
                                  SDP_sp_index,
                                  SDP_sp_blockStruct, 
                                  SDP_sp_NonZeroNumber,
                                  SOCP_sp_nBlock, 
                                  SOCP_sp_blockStruct, 
                                  SOCP_sp_index,
                                  SOCP_sp_NonZeroNumber,
                                  LP_sp_nBlock, 
                                  LP_sp_index);
    }
  }

  delete[] SDP_index;
  SDP_index = NULL;
  delete[] SOCP_index;
  SOCP_index = NULL;
  delete[] LP_index;
  LP_index = NULL;

  delete[] SDP_sp_index;
  SDP_sp_index = NULL;
  delete[] SDP_sp_blockStruct;
  SDP_sp_blockStruct = NULL;
  delete[] SDP_sp_NonZeroNumber;
  SDP_sp_NonZeroNumber = NULL;
#if 0
  delete[] SOCP_sp_index;
  SOCP_sp_index = NULL;
  delete[] SOCP_sp_blockStruct;
  SOCP_sp_blockStruct = NULL;
  delete[] SOCP_sp_NonZeroNumber;
  SOCP_sp_NonZeroNumber = NULL;
#endif
  delete[] LP_sp_index;
  LP_sp_index = NULL;
}


  // 2008/02/27 kazuhide nakata   
  // without LP_ANonZeroCount
void IO::setElement(FILE* fpData, InputData& inputData, int m, 
                    int SDP_nBlock, int* SDP_blockStruct, 
                    int SOCP_nBlock, int* SOCP_blockStruct, 
                    int LP_nBlock, 
                    int nBlock, int* blockStruct, 
                    int* blockType, int* blockNumber,
                    long position, bool isDataSparse)
{
  // in Sparse, read C,A[k]

  // seed the positon of C in the fpData
  fseek(fpData, position, 0);

  if (isDataSparse) {
    int i,j,k,l;
    mpf_class value;
    while (true) {
      if (fscanf(fpData,"%*[^0-9+-]%d",&k)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&l)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&i)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&j)<=0) {
	break;
      }
      if (gmp_fscanf(fpData,"%*[^0-9+-]%Fe",&value)<=0) {
	break;
      }
#if 0
      rMessage("input k:" << k <<
	       " l:" << l <<
	       " i:" << i <<
	       " j:" << j);
#endif     

      if (blockType[l-1] == 1){	// SDP part
		int l2 = blockNumber[l-1];
		if (k==0) {
		  inputData.C.setElement_SDP(l2,i-1,j-1,-value);
		} else {
		  inputData.A[k-1].setElement_SDP(l2,i-1,j-1,value);
		}
      } else if (blockType[l-1] == 2){	// SOCP part
		rError("io:: current version does not support SOCP");
		int l2 = blockNumber[l-1];
		if (k==0) {
		  inputData.C.setElement_SOCP(l2,i-1,j-1,-value);
		} else {
		  inputData.A[k-1].setElement_SOCP(l2,i-1,j-1,value);
		}
      } else if (blockType[l-1] == 3){ // LP part
		if (i != j){
		  rError("io:: LP part  3rd elemtn != 4th elemnt");
		}
		if (k==0) {
		  inputData.C.setElement_LP(blockNumber[l-1]+i-1,-value);
		} else {
		  inputData.A[k-1].setElement_LP(blockNumber[l-1]+i-1,value);
		}
      } else {
		rError("io::read not valid blockType");
	  }
	} 
  } else {  // dense

    // k==0
	for (int l2=0; l2<nBlock; ++l2){
	  if (blockType[l2] == 1) { // SDP part
		int l = blockNumber[l2];
		int size = SDP_blockStruct[l];
		for (int i=0; i<size; ++i) {
		  for (int j=0; j<size; ++j) {
			mpf_class tmp;
			gmp_fscanf(fpData,"%*[^0-9+-]%Fe",&tmp);
			if (i<=j && tmp!=0.0) {
			  inputData.C.setElement_SDP(l,i,j,-tmp);
			}
		  }
		}
	  } else if (blockType[l2] == 2) { // SOCP part
		rError("io:: current version does not support SOCP");
	  } else if (blockType[l2] == 3) { // LP part
		for (int j=0; j<blockStruct[l2]; ++j) {
		  mpf_class tmp;
		  gmp_fscanf(fpData,"%*[^0-9+-]%Fe",&tmp);
		  if (tmp!=0.0) {
			inputData.C.setElement_LP(blockNumber[l2]+j,-tmp);
		  }
		}
	  } else {
		rError("io::read not valid blockType");
	  }
	}

	// k > 0
    for (int k=0; k<m; ++k) {
	  
	for (int l2=0; l2<nBlock; ++l2){
	  if (blockType[l2] == 1) { // SDP part
		int l = blockNumber[l2];
		int size = SDP_blockStruct[l];
		for (int i=0; i<size; ++i) {
		  for (int j=0; j<size; ++j) {
			mpf_class tmp;
			gmp_fscanf(fpData,"%*[^0-9+-]%Fe",&tmp);
			if (i<=j && tmp!=0.0) {
			  inputData.A[k].setElement_SDP(l,i,j,tmp);
			}
		  }
		}
	  } else if (blockType[l2] == 2) { // SOCP part
		rError("io:: current version does not support SOCP");
	  } else if (blockType[l2] == 3) { // LP part
		for (int j=0; j<blockStruct[l2]; ++j) {
		  mpf_class tmp;
		  gmp_fscanf(fpData,"%*[^0-9+-]%Fe",&tmp);
		  if (tmp!=0.0) {
			inputData.A[k].setElement_LP(blockNumber[l2]+j,tmp);
		  }
		}
	  } else {
		rError("io::read not valid blockType");
	  }
	}

	} // for k

  } // end of 'if (isDataSparse)'

}

void IO::printHeader(FILE* fpout, FILE* Display)
{
  if (fpout) {
    fprintf(fpout,"   mu      thetaP  thetaD  objP      objD "
	    "     alphaP  alphaD  beta \n");
  }
  if (Display) {
    fprintf(Display,"   mu      thetaP  thetaD  objP      objD "
	    "     alphaP  alphaD  beta \n");
  }
}

void IO::printOneIteration(int pIteration,
			    AverageComplementarity& mu,
			    RatioInitResCurrentRes& theta,
			    SolveInfo& solveInfo,
			    StepLength& alpha,
			    DirectionParameter& beta,
			    FILE* fpout,
			    FILE* Display)
{
  #if REVERSE_PRIMAL_DUAL
  if (Display) {
    mpf_class mtmp1=-solveInfo.objValDual;
    mpf_class mtmp2=-solveInfo.objValPrimal;
    gmp_fprintf(Display,"%2d %4.1Fe %4.1Fe %4.1Fe %+7.2Fe %+7.2Fe"
	    " %4.1Fe %4.1Fe %4.F2e\n", pIteration, mu.current.get_mpf_t(),
	    theta.dual.get_mpf_t(), theta.primal.get_mpf_t(),
	    mtmp1.get_mpf_t(), mtmp2.get_mpf_t(),
	    alpha.dual.get_mpf_t(), alpha.primal.get_mpf_t(), beta.value.get_mpf_t());
  }
  if (fpout) {
    mpf_class mtmp1=-solveInfo.objValDual;
    mpf_class mtmp2=-solveInfo.objValPrimal;
    gmp_fprintf(fpout,"%2d %4.1Fe %4.1Fe %4.1Fe %+7.2Fe %+7.2Fe"
	    " %4.1Fe %4.1Fe %4.2Fe\n", pIteration, mu.current.get_mpf_t(),
	    theta.dual.get_mpf_t(), theta.primal.get_mpf_t(),
	    mtmp1.get_mpf_t(),mtmp2.get_mpf_t(),
	    alpha.dual.get_mpf_t(), alpha.primal.get_mpf_t(), beta.value.get_mpf_t());
  }
  #else
  if (Display) {
    gmp_fprintf(Display,"%2d %4.1Fe %4.1Fe %4.1Fe %+7.2Fe %+7.2Fe"
	    " %4.1Fe %4.1Fe %4.2Fe\n", pIteration.get_mpf_t(), mu.current.get_mpf_t(),
	    theta.primal.get_mpf_t(), theta.dual.get_mpf_t(),
	    solveInfo.objValPrimal.get_mpf_t(), solveInfo.objValDual.get_mpf_t(),
	    alpha.primal.get_mpf_t(), alpha.dual.get_mpf_t(), beta.value.get_mpf_t());
  }
  if (fpout) {
    gmp_fprintf(fpout,"%2d %4.1Fe %4.1Fe %4.1Fe %+7.2Fe %+7.2Fe"
	    " %4.1Fe %4.1Fe %4.2Fe\n", pIteration.get_mpf_t(), mu.current.get_mpf_t(),
	    theta.primal.get_mpf_t(), theta.dual.get_mpf_t(),
	    solveInfo.objValPrimal.get_mpf_t(), solveInfo.objValDual.get_mpf_t(),
	    alpha.primal.get_mpf_t(), alpha.dual.get_mpf_t(), beta.value.get_mpf_t());
  }
  #endif
}

void IO::printLastInfo(int pIteration,
		       AverageComplementarity& mu,
		       RatioInitResCurrentRes& theta,
		       SolveInfo& solveInfo,
		       StepLength& alpha,
		       DirectionParameter& beta,
		       Residuals& currentRes,
		       Phase & phase,
		       Solutions& currentPt,
		       double cputime,
		       InputData& inputData,
                       WorkVariables& work,
		       ComputeTime& com,
		       Parameter& param,
		       FILE* fpout,
		       FILE* Display,
		       bool printTime)
{
  int nDim = currentPt.nDim;

  printOneIteration(pIteration,mu,theta,solveInfo,alpha,
		    beta, fpout, Display);

  mpf_class mean = (abs(solveInfo.objValPrimal)
		 + abs(solveInfo.objValDual)) / 2.0;
  mpf_class PDgap = abs(solveInfo.objValPrimal
		      - solveInfo.objValDual);
  // mpf_class dominator;
  mpf_class relgap;
  if (mean < 1.0) {
    relgap = PDgap;
  } else {
    relgap = PDgap/mean;
  }

  mpf_class gap    = mu.current*nDim; 
  mpf_class digits = -mpf_approx_log10(abs(PDgap/mean));

  #if DIMACS_PRINT
  mpf_class tmp = 0.0;
  mpf_class b1 = 0.0;
  for (int k=0; k<inputData.b.nDim; ++k) {
    tmp= abs(inputData.b.ele[k]);
    b1 = max(b1, tmp);
  }
  mpf_class c1 = 0.0;
  for (int l=0; l<inputData.C.SDP_sp_nBlock; ++l) {
    SparseMatrix& Cl = inputData.C.SDP_sp_block[l];
    if (Cl.type == SparseMatrix::SPARSE) {
      for (int i=0; i<Cl.NonZeroCount; ++i) {
	tmp= abs(Cl.sp_ele[i]);
	c1 = max(c1, tmp);
      }
    } else if (Cl.type == SparseMatrix::DENSE) {
      for (int i=0; i<Cl.nRow*Cl.nCol; ++i) {
	tmp= abs(Cl.de_ele[i]);
	c1 = max(c1, tmp);
      }
    }
  }
  for (int l=0; l<inputData.C.SOCP_sp_nBlock; ++l) {
    rError("io:: current version does not support SOCP");
  }
  for (int l=0; l<inputData.C.LP_sp_nBlock; ++l) {
    tmp = abs(inputData.C.LP_sp_block[l]);
    c1 = max(c1, tmp);
  }
  mpf_class p_norm;
  Lal::let(tmp,'=',currentRes.primalVec,'.',currentRes.primalVec);
  p_norm = sqrt(tmp);
  mpf_class d_norm = 0.0;
  for (int l=0; l<currentRes.dualMat.SDP_nBlock; ++l) {
    Lal::let(tmp,'=',currentRes.dualMat.SDP_block[l],'.',currentRes.dualMat.SDP_block[l]);
    d_norm += sqrt(tmp);
  }
  for (int l=0; l<currentRes.dualMat.SOCP_nBlock; ++l) {
    rError("io:: current version does not support SOCP");
  }
  tmp = 0.0;
  for (int l=0; l<currentRes.dualMat.LP_nBlock; ++l) {
    tmp += currentRes.dualMat.LP_block[l] * currentRes.dualMat.LP_block[l];
  }
  d_norm += sqrt(tmp);
  mpf_class x_min =  Jal::getMinEigen(currentPt.xMat,work);
  mpf_class z_min =  Jal::getMinEigen(currentPt.zMat,work);
					
  // gmp_printf("b1:%Fe\n",b1);
  // gmp_printf("c1:%Fe\n",c1);
  // gmp_printf("p_norm:%Fe\n",p_norm);
  // gmp_printf("d_norm:%Fe\n",d_norm);
  // gmp_printf("x_min:%Fe\n",x_min);
  // gmp_printf("z_min:%Fe\n",z_min);
  
  mpf_class ctx = solveInfo.objValPrimal;
  mpf_class bty = solveInfo.objValDual;
  mpf_class xtz = 0.0;
  Lal::let(xtz,'=',currentPt.xMat,'.',currentPt.zMat);

  mpf_class mzero = 0.0;
  mpf_class err1 = p_norm / (1+b1);
  mpf_class err2 = max( mzero, - x_min / (1+b1));
  mpf_class err3 = d_norm / (1+c1);
  mpf_class err4 = max( mzero, - z_min / (1+c1));
  mpf_class err5 = (ctx - bty) / (1 + abs(ctx) + abs(bty));
  mpf_class err6 = xtz / (1 + abs(ctx) + abs(bty));
    
  #endif
  if (Display) {
    fprintf(Display, "\n");
    phase.display(Display);
        fprintf(Display, "   Iteration = %d\n",       pIteration);
    gmp_fprintf(Display, "          mu = %4.16Fe\n",  mu.current.get_mpf_t());
    gmp_fprintf(Display, "relative gap = %4.16Fe\n",  relgap.get_mpf_t());
    gmp_fprintf(Display, "         gap = %4.16Fe\n",  gap.get_mpf_t());
    gmp_fprintf(Display, "      digits = %4.16Fe\n",  digits.get_mpf_t());

    #if REVERSE_PRIMAL_DUAL
    mpf_class mtmp1 = -solveInfo.objValDual;
    mpf_class mtmp2 = -solveInfo.objValPrimal;
    gmp_fprintf(Display, "objValPrimal = %10.16Fe\n",
	    mtmp1.get_mpf_t());
    gmp_fprintf(Display, "objValDual   = %10.16Fe\n",
	    mtmp2.get_mpf_t());
    gmp_fprintf(Display, "p.feas.error = %10.16Fe\n",
	    currentRes.normDualMat.get_mpf_t());
    gmp_fprintf(Display, "d.feas.error = %10.16Fe\n",
	    currentRes.normPrimalVec.get_mpf_t());
    gmp_fprintf(Display, "relative eps = %10.16Fe\n",
            Rlamch_gmp("E").get_mpf_t());
    #else
    gmp_fprintf(Display, "objValPrimal = %10.16Fe\n",
	    solveInfo.objValPrimal.get_mpf_t());
    gmp_fprintf(Display, "objValDual   = %10.16Fe\n",
	    solveInfo.objValDual.get_mpf_t());
    gmp_fprintf(Display, "p.feas.error = %10.16Fe\n",
	    currentRes.normPrimalVec.get_mpf_t());
    gmp_fprintf(Display, "d.feas.error = %10.16Fe\n",
	    currentRes.normDualMat.get_mpf_t());
    gmp_fprintf(Display, "relative eps = %10.16Fe\n",
            gmp_dlamchE().get_mpf_t());
    #endif
    if (printTime == true) {
      fprintf(Display, "total time   = %.3f\n",cputime);
    }
    #if DIMACS_PRINT
    fprintf(Display, "\n");
    gmp_fprintf(Display, "* DIMACS_ERRORS * \n");
    gmp_fprintf(Display, "err1 = %4.16Fe  [%40s]\n",
	    err1.get_mpf_t(), "||Ax-b|| / (1+||b||_1) ");
    gmp_fprintf(Display, "err2 = %4.16Fe  [%40s]\n",
	    err2.get_mpf_t(), "max(0, -lambda(x) / (1+||b||_1))");
    gmp_fprintf(Display, "err3 = %4.16Fe  [%40s]\n",
	    err3.get_mpf_t(), "||A^Ty + z - c || / (1+||c||_1) ");
    gmp_fprintf(Display, "err4 = %4.16Fe  [%40s]\n",
	    err4.get_mpf_t(), "max(0, -lambda(z) / (1+||c||_1))");
    gmp_fprintf(Display, "err5 = %4.16Fe  [%40s]\n",
	    err5.get_mpf_t(),"(<c,x> - by) / (1 + |<c,x>| + |by|)");
    gmp_fprintf(Display, "err6 = %4.16Fe  [%40s]\n",
	    err6.get_mpf_t(),"<x,z> / (1 + |<c,x>| + |by|)");
    fprintf(Display, "\n");
    #endif
  }
  if (fpout) {
    fprintf(fpout, "\n");
    phase.display(fpout);
        fprintf(fpout, "   Iteration = %d\n",  pIteration);
    gmp_fprintf(fpout, "          mu = %4.16Fe\n",  mu.current.get_mpf_t());
    gmp_fprintf(fpout, "relative gap = %4.16Fe\n",  relgap.get_mpf_t());
    gmp_fprintf(fpout, "         gap = %4.16Fe\n",  gap.get_mpf_t());
    gmp_fprintf(fpout, "      digits = %4.16Fe\n",  digits.get_mpf_t());

    #if REVERSE_PRIMAL_DUAL
    mpf_class mtmp1=-solveInfo.objValDual;
    mpf_class mtmp2=-solveInfo.objValPrimal;
    gmp_fprintf(fpout, "objValPrimal = %10.16Fe\n",
	    mtmp1.get_mpf_t());
    gmp_fprintf(fpout, "objValDual   = %10.16Fe\n",
	    mtmp2.get_mpf_t());
    gmp_fprintf(fpout, "p.feas.error = %10.16Fe\n",
	    currentRes.normDualMat.get_mpf_t());
    gmp_fprintf(fpout, "d.feas.error = %10.16Fe\n",
	    currentRes.normPrimalVec.get_mpf_t());
    gmp_fprintf(fpout, "relative eps = %10.16Fe\n",
            Rlamch_gmp("E").get_mpf_t());
    #else
    gmp_fprintf(fpout, "objValPrimal = %10.16Fe\n",
	    solveInfo.objValPrimal.get_mpf_t());
    gmp_fprintf(fpout, "objValDual   = %10.16Fe\n",
	    solveInfo.objValDual.get_mpf_t());
    gmp_fprintf(fpout, "p.feas.error = %10.16Fe\n",
	    currentRes.normPrimalVec.get_mpf_t());
    gmp_fprintf(fpout, "d.feas.error = %10.16Fe\n",
	    currentRes.normDualMat.get_mpf_t());
    gmp_fprintf(fpout, "relative eps = %10.16Fe\n",
            Rlamch_gmp("E").get_mpf_t());
    #endif
    fprintf(fpout, "total time   = %.3f\n",cputime);
    #if DIMACS_PRINT
    fprintf(fpout, "\n");
    gmp_fprintf(fpout, "* DIMACS_ERRORS * \n");
    gmp_fprintf(fpout, "err1 = %4.16Fe  [%40s]\n",
	    err1.get_mpf_t(), "||Ax-b|| / (1+||b||_1) ");
    gmp_fprintf(fpout, "err2 = %4.16Fe  [%40s]\n",
	    err2.get_mpf_t(), "max(0, -lambda(x) / (1+||b||_1))");
    gmp_fprintf(fpout, "err3 = %4.16Fe  [%40s]\n",
	    err3.get_mpf_t(), "||A^Ty + z - c || / (1+||c||_1) ");
    gmp_fprintf(fpout, "err4 = %4.16Fe  [%40s]\n",
	    err4.get_mpf_t(), "max(0, -lambda(z) / (1+||c||_1))");
    gmp_fprintf(fpout, "err5 = %4.16Fe  [%40s]\n",
	    err5.get_mpf_t(),"(<c,x> - by) / (1 + |<c,x>| + |by|)");
    gmp_fprintf(fpout, "err6 = %4.16Fe  [%40s]\n",
	    err6.get_mpf_t(),"<x,z> / (1 + |<c,x>| + |by|)");
    fprintf(fpout, "\n");
    #endif

    gmp_fprintf(fpout, "\n\nParameters are\n");
    param.display(fpout);
    com.display(fpout);

    #if 1
    #if REVERSE_PRIMAL_DUAL
    fprintf(fpout,"xVec = \n");
    currentPt.yVec.display(fpout,-1.0);
    fprintf(fpout,"xMat = \n");
    currentPt.zMat.display(fpout);
    fprintf(fpout,"yMat = \n");
    currentPt.xMat.display(fpout);
    #else
    currentPt.display(fpout);
    #endif
    #endif
  }
}


void IO::printLastInfo(int pIteration,
		       AverageComplementarity& mu,
		       RatioInitResCurrentRes& theta,
		       SolveInfo& solveInfo,
		       StepLength& alpha,
		       DirectionParameter& beta,
		       Residuals& currentRes,
		       Phase & phase,
		       Solutions& currentPt,
		       double cputime,
		       int nBlock,
		       int* blockStruct,
		       int* blockType,
		       int* blockNumber,
		       InputData& inputData,
                       WorkVariables& work,
		       ComputeTime& com,
		       Parameter& param,
		       FILE* fpout,
		       FILE* Display,
		       bool printTime)
{
  int nDim = currentPt.nDim;

  printOneIteration(pIteration,mu,theta,solveInfo,alpha,
		    beta, fpout, Display);

  mpf_class mean = (abs(solveInfo.objValPrimal)
		 + abs(solveInfo.objValDual)) / 2.0;
  mpf_class PDgap = abs(solveInfo.objValPrimal
		      - solveInfo.objValDual);
  // mpf_class dominator;
  mpf_class relgap;
  if (mean < 1.0) {
    relgap = PDgap;
  } else {
    relgap = PDgap/mean;
  }
  mpf_class gap    = mu.current*nDim;
  mpf_class digits = 0.0;
  if (PDgap != 0.0 && mean != 0.0) {
    digits = -mpf_approx_log10(abs(PDgap/mean));
  }

  #if DIMACS_PRINT
  mpf_class tmp = 0.0;
  mpf_class b1 = 0.0;
  for (int k=0; k<inputData.b.nDim; ++k) {
    tmp = abs(inputData.b.ele[k]);
    b1 = max(b1, tmp);
  }
  mpf_class c1 = 0.0;
  for (int l=0; l<inputData.C.SDP_sp_nBlock; ++l) {
    SparseMatrix& Cl = inputData.C.SDP_sp_block[l];
    if (Cl.type == SparseMatrix::SPARSE) {
      for (int i=0; i<Cl.NonZeroCount; ++i) {
	tmp = abs(Cl.sp_ele[i]);
	c1 = max(c1, tmp);
      }
    } else if (Cl.type == SparseMatrix::DENSE) {
      for (int i=0; i<Cl.nRow*Cl.nCol; ++i) {
	tmp = abs(Cl.de_ele[i]);
	c1 = max(c1, tmp);
      }
    }
  }
  for (int l=0; l<inputData.C.SOCP_sp_nBlock; ++l) {
    rError("io:: current version does not support SOCP");
  }
  for (int l=0; l<inputData.C.LP_sp_nBlock; ++l) {
    tmp = abs(inputData.C.LP_sp_block[l]);
    c1 = max(c1, tmp);
  }
  mpf_class p_norm;
  Lal::let(tmp,'=',currentRes.primalVec,'.',currentRes.primalVec);
  p_norm = sqrt(tmp);
  mpf_class d_norm = 0.0;
  for (int l=0; l<currentRes.dualMat.SDP_nBlock; ++l) {
    Lal::let(tmp,'=',currentRes.dualMat.SDP_block[l],'.',currentRes.dualMat.SDP_block[l]);
    d_norm += sqrt(tmp);
  }
  for (int l=0; l<currentRes.dualMat.SOCP_nBlock; ++l) {
    rError("io:: current version does not support SOCP");
  }
  tmp = 0.0;
  for (int l=0; l<currentRes.dualMat.LP_nBlock; ++l) {
    tmp += currentRes.dualMat.LP_block[l] * currentRes.dualMat.LP_block[l];
  }
  d_norm += sqrt(tmp);
  mpf_class x_min =  Jal::getMinEigen(currentPt.xMat,work);
  mpf_class z_min =  Jal::getMinEigen(currentPt.zMat,work);
					
  //  gmp_printf("b1:%Fe\n",b1);
  //  gmp_printf("c1:%Fe\n",c1);
  //  gmp_printf("p_norm:%Fe\n",p_norm);
  //  gmp_printf("d_norm:%Fe\n",d_norm);
  //  gmp_printf("x_min:%Fe\n",x_min);
  //  gmp_printf("z_min:%Fe\n",z_min);
  
  mpf_class ctx = solveInfo.objValPrimal;
  mpf_class bty = solveInfo.objValDual;
  mpf_class xtz = 0.0;
  Lal::let(xtz,'=',currentPt.xMat,'.',currentPt.zMat);

  mpf_class mzero = 0.0;
  mpf_class err1 = p_norm / (1+b1);
  mpf_class err2 = max( mzero, - x_min / (1+b1));
  mpf_class err3 = d_norm / (1+c1);
  mpf_class err4 = max( mzero, - z_min / (1+c1));
  mpf_class err5 = (ctx - bty) / (1 + abs(ctx) + abs(bty));
  mpf_class err6 = xtz / (1 + abs(ctx) + abs(bty));
    
  #endif
  
  if (Display) {
    fprintf(Display, "\n");
    phase.display(Display);
    gmp_fprintf(Display, "   Iteration = %d\n",       pIteration);
    gmp_fprintf(Display, "          mu = %4.16Fe\n",  mu.current.get_mpf_t());
    gmp_fprintf(Display, "relative gap = %4.16Fe\n",  relgap.get_mpf_t());
    gmp_fprintf(Display, "         gap = %4.16Fe\n",  gap.get_mpf_t());
    gmp_fprintf(Display, "      digits = %4.16Fe\n",  digits.get_mpf_t());

    #if REVERSE_PRIMAL_DUAL
    mpf_class mtmp1 = -solveInfo.objValDual;
    mpf_class mtmp2 = -solveInfo.objValPrimal;
    gmp_fprintf(Display, "objValPrimal = %10.16Fe\n",
	    mtmp1.get_mpf_t());
    gmp_fprintf(Display, "objValDual   = %10.16Fe\n",
	    mtmp2.get_mpf_t());
    gmp_fprintf(Display, "p.feas.error = %10.16Fe\n",
	    currentRes.normDualMat.get_mpf_t());
    gmp_fprintf(Display, "d.feas.error = %10.16Fe\n",
	    currentRes.normPrimalVec.get_mpf_t());
    gmp_fprintf(Display, "relative eps = %10.16Fe\n",
            Rlamch_gmp("E").get_mpf_t());
    #else
    gmp_fprintf(Display, "objValPrimal = %10.16Fe\n",
	    solveInfo.objValPrimal.get_mpf_t());
    gmp_fprintf(Display, "objValDual   = %10.16Fe\n",
	    solveInfo.objValDual.get_mpf_t());
    gmp_fprintf(Display, "p.feas.error = %10.16Fe\n",
	    currentRes.normPrimalVec.get_mpf_t());
    gmp_fprintf(Display, "d.feas.error = %10.16Fe\n",
	    currentRes.normDualMat.get_mpf_t());
    gmp_fprintf(Display, "relative eps = %10.16Fe\n",
            Rlamch_gmp("E").get_mpf_t());
    #endif
    if (printTime == true) {
      fprintf(Display, "total time   = %.3f\n",cputime);
    }
    #if DIMACS_PRINT
    fprintf(Display, "\n");
    gmp_fprintf(Display, "* DIMACS_ERRORS * \n");
    gmp_fprintf(Display, "err1 = %4.16Fe  [%40s]\n",
	    err1.get_mpf_t(), "||Ax-b|| / (1+||b||_1) ");
    gmp_fprintf(Display, "err2 = %4.16Fe  [%40s]\n",
	    err2.get_mpf_t(), "max(0, -lambda(x) / (1+||b||_1))");
    gmp_fprintf(Display, "err3 = %4.16Fe  [%40s]\n",
	    err3.get_mpf_t(), "||A^Ty + z - c || / (1+||c||_1) ");
    gmp_fprintf(Display, "err4 = %4.16Fe  [%40s]\n",
	    err4.get_mpf_t(), "max(0, -lambda(z) / (1+||c||_1))");
    gmp_fprintf(Display, "err5 = %4.16Fe  [%40s]\n",
	    err5.get_mpf_t(),"(<c,x> - by) / (1 + |<c,x>| + |by|)");
    gmp_fprintf(Display, "err6 = %4.16Fe  [%40s]\n",
	    err6.get_mpf_t(),"<x,z> / (1 + |<c,x>| + |by|)");
    fprintf(Display, "\n");
    #endif
  }
  if (fpout) {
    fprintf(fpout, "\n");
    phase.display(fpout);
    fprintf(fpout, "   Iteration = %d\n",  pIteration);
    gmp_fprintf(fpout, "          mu = %4.16Fe\n",  mu.current.get_mpf_t());
    gmp_fprintf(fpout, "relative gap = %4.16Fe\n",  relgap.get_mpf_t());
    gmp_fprintf(fpout, "         gap = %4.16Fe\n",  gap.get_mpf_t());
    gmp_fprintf(fpout, "      digits = %4.16Fe\n",  digits.get_mpf_t());

    #if REVERSE_PRIMAL_DUAL
    mpf_class mtmp1=-solveInfo.objValDual;
    mpf_class mtmp2=-solveInfo.objValPrimal;
    gmp_fprintf(fpout, "objValPrimal = %10.16Fe\n",
	    mtmp1.get_mpf_t());
    gmp_fprintf(fpout, "objValDual   = %10.16Fe\n",
	    mtmp2.get_mpf_t());
    gmp_fprintf(fpout, "p.feas.error = %10.16Fe\n",
	    currentRes.normDualMat.get_mpf_t());
    gmp_fprintf(fpout, "d.feas.error = %10.16Fe\n",
	    currentRes.normPrimalVec.get_mpf_t());
    gmp_fprintf(fpout, "relative eps = %10.16Fe\n",
            Rlamch_gmp("E").get_mpf_t());
    #else
    gmp_fprintf(fpout, "objValPrimal = %10.16Fe\n",
	    solveInfo.objValPrimal.get_mpf_t());
    gmp_fprintf(fpout, "objValDual   = %10.16Fe\n",
	    solveInfo.objValDual.get_mpf_t());
    gmp_fprintf(fpout, "p.feas.error = %10.16Fe\n",
	    currentRes.normPrimalVec.get_mpf_t());
    gmp_fprintf(fpout, "d.feas.error = %10.16Fe\n",
	    currentRes.normDualMat.get_mpf_t());
    gmp_fprintf(fpout, "relative eps = %10.16Fe\n",
            Rlamch_gmp("E").get_mpf_t());
    #endif
    fprintf(fpout, "total time   = %.3f\n",cputime);
    #if DIMACS_PRINT
    fprintf(fpout, "\n");
    gmp_fprintf(fpout, "* DIMACS_ERRORS * \n");
    gmp_fprintf(fpout, "err1 = %4.16Fe  [%40s]\n",
	    err1.get_mpf_t(), "||Ax-b|| / (1+||b||_1) ");
    gmp_fprintf(fpout, "err2 = %4.16Fe  [%40s]\n",
	    err2.get_mpf_t(), "max(0, -lambda(x) / (1+||b||_1))");
    gmp_fprintf(fpout, "err3 = %4.16Fe  [%40s]\n",
	    err3.get_mpf_t(), "||A^Ty + z - c || / (1+||c||_1) ");
    gmp_fprintf(fpout, "err4 = %4.16Fe  [%40s]\n",
	    err4.get_mpf_t(), "max(0, -lambda(z) / (1+||c||_1))");
    gmp_fprintf(fpout, "err5 = %4.16Fe  [%40s]\n",
	    err5.get_mpf_t(),"(<c,x> - by) / (1 + |<c,x>| + |by|)");
    gmp_fprintf(fpout, "err6 = %4.16Fe  [%40s]\n",
	    err6.get_mpf_t(),"<x,z> / (1 + |<c,x>| + |by|)");
    fprintf(fpout, "\n");
    #endif

    fprintf(fpout, "\n\nParameters are\n");
    param.display(fpout);
    com.display(fpout);

    #if 1
    #if REVERSE_PRIMAL_DUAL
    fprintf(fpout,"xVec = \n");
    currentPt.yVec.display(fpout,-1.0);
    fprintf(fpout,"xMat = \n");
    displayDenseLinarSpaceLast(currentPt.zMat,
							   nBlock, blockStruct, blockType, blockNumber,
							   fpout);
    fprintf(fpout,"yMat = \n");
    displayDenseLinarSpaceLast(currentPt.xMat,
							   nBlock, blockStruct, blockType, blockNumber,
							   fpout);
    #else
	fprintf(fpout,"xMat = \n");
    displayDenseLinarSpaceLast(currentPt.xMat,
							   nBlock, blockStruct, blockType, blockNumber,
							   fpout);
	fprintf(fpout,"yVec = \n");
	currentPt.yVec.display(fpout);
	fprintf(fpout,"zMat = \n");
    displayDenseLinarSpaceLast(currentPt.zMat,
							   nBlock, blockStruct, blockType, blockNumber,
							   fpout);
    #endif
    #endif
  }
}

void IO::displayDenseLinarSpaceLast(DenseLinearSpace& aMat,
									int nBlock,
									int* blockStruct,
									int* blockType,
									int* blockNumber,
									FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }

  fprintf(fpout,"{\n");
  for (int i=0; i<nBlock; i++){
	if (blockType[i] == 1){
	  int l = blockNumber[i];
	  aMat.SDP_block[l].display(fpout);
	} else if (blockType[i] == 2){
	  rError("io:: current version does not support SOCP");
	  int l = blockNumber[i];
	  aMat.SOCP_block[l].display(fpout);
	} else if (blockType[i] == 3){
	  fprintf(fpout,"{");
	  for (int l=0; l<blockStruct[i]-1; ++l) {
		gmp_fprintf(fpout,P_FORMAT",",aMat.LP_block[blockNumber[i]+l].get_mpf_t());
	  }
	  if (blockStruct[i] > 0) {
		gmp_fprintf(fpout,P_FORMAT"}\n",aMat.LP_block[blockNumber[i]+blockStruct[i]-1].get_mpf_t());
	  } else {
		fprintf(fpout,"  }\n");
	  }
	} else {
		rError("io::displayDenseLinarSpaceLast not valid blockType");
	}
  }
  fprintf(fpout,"}\n");
}


}
