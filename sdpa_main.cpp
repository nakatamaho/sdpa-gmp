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

#ifndef _MAIN_
#define _MAIN_
#endif

#define LengthOfBuffer 1024
static double KAPPA = 1.2;

#include <sdpa_io.h>

#if UseMETIS
#ifdef __cplusplus
extern "C" {
#endif
#include "metis.h"
#ifdef __cplusplus
}
#endif
#endif

namespace sdpa {

bool pinpal(char* dataFile, char* initFile, char* outFile,
	    char* paraFile, bool isInitFile, bool isInitSparse,
	    bool isDataSparse, bool isParameter,
	    Parameter::parameterType parameterType,
	    FILE* Display)
{

  TimeStart(TOTAL_TIME_START1);
  TimeStart(FILE_READ_START1);
  ComputeTime com;

  FILE* fpData      = NULL;
  FILE* fpOut       = NULL;

  if ((fpOut=fopen(outFile,"w"))==NULL) {
    rError("Cannot open out file " << outFile);
  }
  Parameter param;
  param.setDefaultParameter(parameterType);
  if (isParameter) {
    FILE* fpParameter = NULL;
    if ((fpParameter=fopen(paraFile,"r"))==NULL) {
      fprintf(Display,"Cannot open parameter file %s \n",
	      paraFile);
      exit(0);
    } else {
      param.readFile(fpParameter);
      fclose(fpParameter);
    }
  }
  // param.display(Display);

  if ((fpData=fopen(dataFile,"r"))==NULL) {
    rError("Cannot open data file " << dataFile);
  }
  char titleAndComment[LengthOfBuffer];
  int m;
  time_t ltime;
  time( &ltime );
  fprintf(fpOut,"SDPA-GMP start at %s",ctime(&ltime));
  IO::read(fpData,fpOut,m,titleAndComment);
  fprintf(fpOut,"data      is %s\n",dataFile);
  if (paraFile) {
    fprintf(fpOut,"parameter is %s\n",paraFile);
  }
  if (initFile) {
    fprintf(fpOut,"initial   is %s\n",initFile);
  }
  fprintf(fpOut,"out       is %s\n",outFile);


#if 1 // 2007/11/28 nakata    for multi LP block

  int SDP_nBlock, SOCP_nBlock,LP_nBlock, nBlock;
  IO::read(fpData,nBlock);

  int* blockStruct = NULL;
  int* blockType = NULL;
  int* blockNumber = NULL;
  int* SDP_blockStruct = NULL;
  int* SOCP_blockStruct = NULL;
  blockStruct = new int[nBlock];
  if (blockStruct==NULL) {
    rError("Memory exhausted about blockStruct");
  }
  blockType = new int[nBlock];
  if (blockType==NULL) {
    rError("Memory exhausted about blockType");
  }
  blockNumber = new int[nBlock];
  if (blockNumber==NULL) {
    rError("Memory exhausted about blockNumber");
  }
  IO::read(fpData,nBlock,blockStruct);

  SDP_nBlock = 0;
  SOCP_nBlock = 0;
  LP_nBlock = 0;
  for (int i=0; i<nBlock; i++){
	if (blockStruct[i] >= 2) {
	  blockType[i] = 1;
	  blockNumber[i] = SDP_nBlock;
	  SDP_nBlock++;
	} else if (blockStruct[i] < 0) {
	  blockType[i] = 3;
	  blockStruct[i] = - blockStruct[i];
	  blockNumber[i] = LP_nBlock;
	  LP_nBlock += blockStruct[i];
	} else if (blockStruct[i] == 1) {
	  blockType[i] = 3;
	  blockNumber[i] = LP_nBlock;
	  LP_nBlock += blockStruct[i];
	} else {
	  rError("block struct");
	}
  }

  SDP_blockStruct = new int[SDP_nBlock];
  if (SDP_blockStruct==NULL) {
    rError("Memory exhausted about SDP blockStruct");
  }
  SOCP_blockStruct = new int[SOCP_nBlock];
  if (SOCP_blockStruct==NULL) {
    rError("Memory exhausted about SOCP blockStruct");
  }

  SDP_nBlock = 0;
  for (int i=0; i<nBlock; i++){
	if (blockType[i] == 1) {
	  SDP_blockStruct[SDP_nBlock] = blockStruct[i];
	  SDP_nBlock++;
	} 
  }

  InputData inputData;
  //  rMessage("read input data: start");
  IO::read(fpData, m, SDP_nBlock, SDP_blockStruct,
	   SOCP_nBlock, SOCP_blockStruct, LP_nBlock,
	   nBlock, blockStruct, blockType, blockNumber,
	   inputData,isDataSparse);
  //  rMessage("read input data: end");
  inputData.initialize_index(SDP_nBlock,SOCP_nBlock,LP_nBlock,com);
#else 

  int SDP_nBlock, SOCP_nBlock,LP_nBlock;
  IO::read(fpData,SDP_nBlock,SOCP_nBlock,LP_nBlock);
  int* SDP_blockStruct = NULL;
  int* SOCP_blockStruct = NULL;
  SDP_blockStruct = new int[SDP_nBlock];
  if (SDP_blockStruct==NULL) {
    rError("Memory exhausted about SDP blockStruct");
  }
  SOCP_blockStruct = new int[SOCP_nBlock];
  if (SOCP_blockStruct==NULL) {
    rError("Memory exhausted about SOCP blockStruct");
  }
  IO::read(fpData,SDP_nBlock,SDP_blockStruct,
	   SOCP_nBlock,SOCP_blockStruct, LP_nBlock);


  for (int i=0; i<SDP_nBlock-1; i++){
    if (SDP_blockStruct[i] < 0){
      rError("LP block must be in last block");
    }
  }
  // muriyari nyuuryoku saseru
  if (SDP_blockStruct[SDP_nBlock-1] < 0){
    LP_nBlock = - SDP_blockStruct[SDP_nBlock-1];
    SDP_nBlock--;
  }

  InputData inputData;
  IO::read(fpData, m, SDP_nBlock, SDP_blockStruct,
	   SOCP_nBlock, SOCP_blockStruct, LP_nBlock,
	   inputData,isDataSparse);
  inputData.initialize_index(SDP_nBlock,SOCP_nBlock,LP_nBlock,com);

#endif

  fclose(fpData);
  
#if 0
  inputData.display();
#endif

#if 1
  TimeStart(FILE_CHANGE_START1);
  // if possible , change C and A to Dense
  inputData.C.changeToDense();
  for (int k=0; k<m; ++k) {
    inputData.A[k].changeToDense();
  }
  TimeEnd(FILE_CHANGE_END1);
  com.FileChange += TimeCal(FILE_CHANGE_START1,
			     FILE_CHANGE_END1);
#endif

  // rMessage("C = ");
  // inputData.C.display(Display);
  // for (int k=0; k<m; ++k) {
  //   rMessage("A["<<k<<"] = ");
  //   inputData.A[k].display(Display);
  //   }
  
  // the end of initialization of C and A

  Newton newton(m, SDP_nBlock, SDP_blockStruct,
		SOCP_nBlock, SOCP_blockStruct,	LP_nBlock);
  int nBlock2 = SDP_nBlock+SOCP_nBlock+LP_nBlock;
  // 2008/03/12 kazuhide nakata
  Chordal chordal;
  // rMessage("ordering bMat: start");
  chordal.ordering_bMat(m, nBlock2, inputData, fpOut);
  // rMessage("ordering bMat: end");
  newton.initialize_bMat(m, chordal,inputData, fpOut);
  chordal.terminate();

  //  rMessage("newton.computeFormula_SDP: start");
  newton.computeFormula_SDP(inputData,0.0,KAPPA);
  //  rMessage("newton.computeFormula_SDP: end");

  // set initial solutions.
  Solutions currentPt;
  WorkVariables work;
  DenseLinearSpace initPt_xMat;
  DenseLinearSpace initPt_zMat;

  currentPt.initialize(m, SDP_nBlock, SDP_blockStruct,
		       SOCP_nBlock, SOCP_blockStruct, LP_nBlock,
		       param.lambdaStar,com);
  work.initialize(m, SDP_nBlock, SDP_blockStruct,
		  SOCP_nBlock, SOCP_blockStruct, LP_nBlock);

  if (isInitFile) {
    FILE* fpInit = NULL;
    if ((fpInit=fopen(initFile,"r"))==NULL) {
      rError("Cannot open init file " << initFile);
    }
    IO::read(fpInit,currentPt.xMat,currentPt.yVec,currentPt.zMat,
	     isInitSparse);
    fclose(fpInit);
    currentPt.computeInverse(work,com);

    initPt_xMat.initialize(SDP_nBlock, SDP_blockStruct,
			   SOCP_nBlock, SOCP_blockStruct, 
			   LP_nBlock);
    
    initPt_zMat.initialize(SDP_nBlock, SDP_blockStruct,
			   SOCP_nBlock, SOCP_blockStruct, 
			   LP_nBlock);
    initPt_xMat.copyFrom(currentPt.xMat);
    initPt_zMat.copyFrom(currentPt.zMat);

  }
  //  rMessage("initial xMat = "); initPt_xMat.display(Display);
  //  rMessage("initial yVec = "); currentPt.yVec.display(Display);
  //  rMessage("initial zMat = "); initPt_zMat.display(Display);
  //  rMessage("current pt = "); currentPt.display(Display);
  
  TimeEnd(FILE_READ_END1);
  com.FileRead += TimeCal(FILE_READ_START1,
			   FILE_READ_END1);
  // -------------------------------------------------------------
  // the end of file read
  // -------------------------------------------------------------
  
  Residuals initRes(m, SDP_nBlock, SDP_blockStruct,
		    SOCP_nBlock, SOCP_blockStruct,
		    LP_nBlock, inputData, currentPt);
  Residuals currentRes;
  currentRes.copyFrom(initRes);
  // rMessage("initial currentRes = ");
  // currentRes.display(Display);


  StepLength alpha;
  DirectionParameter beta(param.betaStar);
  Switch reduction(Switch::ON);
  AverageComplementarity mu(param.lambdaStar);

  // rMessage("init mu"); mu.display();

  if (isInitFile) {
    mu.initialize(currentPt);
  }

  RatioInitResCurrentRes theta(param, initRes);
  SolveInfo solveInfo(inputData, currentPt, 
			mu.initial, param.omegaStar);
  Phase phase(initRes, solveInfo, param, currentPt.nDim);

  int pIteration = 0;
  IO::printHeader(fpOut, Display);
  // -----------------------------------------------------
  // Here is MAINLOOP
  // -----------------------------------------------------

  TimeStart(MAIN_LOOP_START1);

  // explicit maxIteration
  // param.maxIteration = 2;
  while (phase.updateCheck(currentRes, solveInfo, param)
	 && pIteration < param.maxIteration) {
    // rMessage(" turn hajimari " << pIteration );
    // Mehrotra's Predictor
    TimeStart(MEHROTRA_PREDICTOR_START1);
    // set variable of Mehrotra
    reduction.MehrotraPredictor(phase);
    beta.MehrotraPredictor(phase, reduction, param);

    // rMessage("reduction = "); reduction.display();
    // rMessage("phase = "); phase.display();
    // rMessage("beta.predictor.value = " << beta.value);
    // rMessage(" mu = " << mu.current);
    // rMessage("currentPt = "); currentPt.display();

    bool isSuccessCholesky;
    isSuccessCholesky = newton.Mehrotra(Newton::PREDICTOR,
					inputData, currentPt,
					currentRes,
					mu, beta, reduction,
					phase,work,com);
    if (isSuccessCholesky == false) {
      break;
    }
    // rMessage("newton predictor = "); newton.display();

    TimeEnd(MEHROTRA_PREDICTOR_END1);
    com.Predictor += TimeCal(MEHROTRA_PREDICTOR_START1,
			      MEHROTRA_PREDICTOR_END1);
    
    TimeStart(STEP_PRE_START1);
    alpha.MehrotraPredictor(inputData, currentPt, phase, newton,
			      work, com);
    // rMessage("alpha predictor = "); alpha.display();

    TimeStart(STEP_PRE_END1);
    com.StepPredictor += TimeCal(STEP_PRE_START1,STEP_PRE_END1);

    // rMessage("alphaStar = " << param.alphaStar);
    // Mehrotra's Corrector
    // rMessage(" Corrector ");

    TimeStart(CORRECTOR_START1);
    beta.MehrotraCorrector(phase,alpha,currentPt, newton,mu,param);

    // rMessage("beta corrector = " << beta.value);

#if 1 // 2007/08/29 kazuhide nakata
	// add stopping criteria: objValPrimal < ObjValDual
	//	if ((pIteration > 10) &&
	if ((phase.value == SolveInfo::pdFEAS) &&
		((beta.value > 5)||(solveInfo.objValPrimal < solveInfo.objValDual))){
	  break;
	}
#endif

    newton.Mehrotra(Newton::CORRECTOR,
		    inputData, currentPt, currentRes,
		    mu, beta, reduction, phase,work,com);

    // rMessage("currentPt = "); currentPt.display();
    // rMessage("newton corrector = "); newton.display();

    TimeEnd(CORRECTOR_END1);
    com.Corrector += TimeCal(CORRECTOR_START1,
			      CORRECTOR_END1);
    TimeStart(CORRECTOR_STEP_START1);
    alpha.MehrotraCorrector(inputData, currentPt, phase,
			    reduction, newton, mu, theta,
			    work, param, com);
    // rMessage("alpha corrector = "); alpha.display();
    TimeEnd(CORRECTOR_STEP_END1);
    com.StepCorrector += TimeCal(CORRECTOR_STEP_START1,
				  CORRECTOR_STEP_END1);
    // the end of Corrector
    
    IO::printOneIteration(pIteration, mu, theta, solveInfo,
			   alpha, beta, fpOut, Display);

    if (currentPt.update(alpha,newton,work,com)==false) {
      // if step length is too short,
      // we finish algorithm
      rMessage("cannot move");
      //   memo by kazuhide nakata
      //   StepLength::MehrotraCorrector
      //   thetaMax*mu.initial -> thetamax*thetaMax*mu.initial
      break;
    }

    // rMessage("currentPt = "); currentPt.display();
    // rMessage("updated");

    theta.update(reduction,alpha);
    mu.update(currentPt);
    currentRes.update(m,inputData, currentPt, com);
    theta.update_exact(initRes,currentRes);

    if (isInitFile) {
      solveInfo.update(inputData, initPt_xMat, initPt_zMat, currentPt,
		       currentRes, mu, theta, param);
    } else {
      solveInfo.update(param.lambdaStar,inputData, currentPt,
		       currentRes, mu, theta, param);
    }
	// 2007/09/18 kazuhide nakata
	// print information of ObjVal, residual, gap, complementarity
	//	solveInfo.check(inputData, currentPt, currentRes, mu, theta, param);
    pIteration++;
  } // end of MAIN_LOOP

  TimeEnd(MAIN_LOOP_END1);

  com.MainLoop = TimeCal(MAIN_LOOP_START1,
			  MAIN_LOOP_END1);
  currentRes.compute(m,inputData,currentPt);
  TimeEnd(TOTAL_TIME_END1);
  
  com.TotalTime = TimeCal(TOTAL_TIME_START1,
			   TOTAL_TIME_END1);
  #if REVERSE_PRIMAL_DUAL
  phase.reverse();
  #endif
#if 1
  IO::printLastInfo(pIteration, mu, theta, solveInfo, alpha, beta,
		    currentRes, phase, currentPt, com.TotalTime,
			nBlock, blockStruct, blockType, blockNumber,
		    inputData, work, com, param, fpOut, Display);
#else
  IO::printLastInfo(pIteration, mu, theta, solveInfo, alpha, beta,
		    currentRes, phase, currentPt, com.TotalTime,
		    inputData, work, com, param, fpOut, Display);
#endif
  // com.display(fpOut);

  if (SDP_blockStruct) {
    delete[] SDP_blockStruct;
    SDP_blockStruct = NULL;
  }
  if (SOCP_blockStruct) {
    delete[] SOCP_blockStruct;
    SOCP_blockStruct = NULL;
  }
  if (blockStruct) {
    delete[] blockStruct;
    blockStruct = NULL;
  }
  if (blockType) {
    delete[] blockType;
    blockType = NULL;
  }
  if (blockNumber) {
    delete[] blockNumber;
    blockNumber = NULL;
  }
  
  fprintf(Display,   "  main loop time = %.6f\n",com.MainLoop);
  fprintf(fpOut,   "    main loop time = %.6f\n",com.MainLoop);
  fprintf(Display,   "      total time = %.6f\n",com.TotalTime);
  fprintf(fpOut,   "        total time = %.6f\n",com.TotalTime);
  #if 0
  fprintf(Display,   "file  check time = %.6f\n",com.FileCheck);
  fprintf(fpOut,   "  file  check time = %.6f\n",com.FileCheck);
  fprintf(Display,   "file change time = %.6f\n",com.FileChange);
  fprintf(fpOut,   "  file change time = %.6f\n",com.FileChange);
  #endif
  fprintf(Display,   "file   read time = %.6f\n",com.FileRead);
  fprintf(fpOut,   "  file   read time = %.6f\n",com.FileRead);
  fclose(fpOut);


#if 0
  rMessage("memory release");
  currentRes.terminate();
  initRes.terminate();
  currentPt.terminate();
  initPt_xMat.terminate();
  initPt_zMat.terminate();
  newton.terminate();
  work.terminate();
  inputData.terminate();
  com.~ComputeTime();
  param.~Parameter();
  alpha.~StepLength();
  beta.~DirectionParameter();
  reduction.~Switch();
  mu.~AverageComplementarity();
  theta.~RatioInitResCurrentRes();
  solveInfo.~SolveInfo();
  phase.~Phase();
#endif

  return true;
}

static void message(char* argv0)
{
  cout << endl;
  cout << "*** Please assign data file and output file.***" << endl;
  cout << endl;
  cout << "---- option type 1 ------------" << endl;
  cout << argv0 <<" DataFile OutputFile [InitialPtFile]"
    " [-pt parameters]"<< endl;
  cout << "parameters = 0 default, 1 fast (unstable),"
    " 2 slow (stable)" << endl;
  cout << "example1-1: " << argv0
       << " example1.dat example1.result" << endl;
  cout << "example1-2: " << argv0
       << " example1.dat-s example1.result" << endl;
  cout << "example1-3: " << argv0
       << " example1.dat example1.result example1.ini" << endl;
  cout << "example1-4: " << argv0
       << " example1.dat example1.result -pt 2" << endl;

  cout << endl;
  cout << "---- option type 2 ------------" << endl;
  cout << argv0 << " [option filename]+ " << endl;
  cout << "  -dd : data dense :: -ds : data sparse     " << endl;
  cout << "  -id : init dense :: -is : init sparse     " << endl;
  cout << "  -o  : output     :: -p  : parameter       " << endl;
  cout << "  -pt : parameters , 0 default, 1 fast (unstable)" << endl;
  cout << "                     2 slow (stable)         " << endl;
  // cout << "  -k  : Kappa(RealValue)" << endl;
  cout << "example2-1: " << argv0
       << " -o example1.result -dd example1.dat" << endl;
  cout << "example2-2: " << argv0
       << " -ds example1.dat-s -o example2.result "
       << "-p param.sdpa" << endl;
  cout << "example2-3: " << argv0
       << " -ds example1.dat-s -o example3.result "
       << "-pt 2" << endl;
  exit(1);
}
  
}


using namespace sdpa;

int main(int argc, char** argv)
{
  FILE* Display = stdout;
  setbuf(Display,NULL);

  time_t ltime;
  time( &ltime );
  cout << "SDPA-GMP start at    " << ctime(&ltime);
  // << "... (built at "<< __DATE__ << " " <<__TIME__ ")" << endl;
  // cout << "let me see your ..." << endl;

  bool isInitFile   = false;
  bool isInitSparse = false;
  bool isOutFile    = false;
  bool isDataSparse = false;
  bool isParameter  = false;
  
  char* dataFile = NULL;
  char* initFile = NULL;
  char* outFile  = NULL;
  char* paraFile = NULL;

  Parameter::parameterType parameterType =
    Parameter::PARAMETER_DEFAULT;
  
  if (argc == 1) {
    message(argv[0]);
  }
  if (argv[1][0] == '-') {
    // rsdpa argument
    
    for (int index = 0; index < argc; ++index) {
      char* target = argv[index];
      if (strcmp(target,"-dd")==0 && index+1 < argc) {
	dataFile = argv[index+1];
	index++;
	continue;
      }
      if (strcmp(target,"-ds")==0 && index+1 < argc) {
	dataFile = argv[index+1];
	index++;
	isDataSparse = true;
	continue;
      }
      if (strcmp(target,"-id")==0 && index+1 < argc) {
	initFile = argv[index+1];
	index++;
	isInitFile = true;
	continue;
      }
      if (strcmp(target,"-is")==0 && index+1 < argc) {
	initFile = argv[index+1];
	index++;
	isInitFile   = true;
	isInitSparse = true;
	continue;
      }
      if (strcmp(target,"-o")==0 && index+1 < argc) {
	outFile = argv[index+1];
	index++;
	isOutFile = true;
	continue;
      }
      if (strcmp(target,"-p")==0 && index+1 < argc) {
	paraFile = argv[index+1];
	index++;
	isParameter = true;
	continue;
      }
      if (strcmp(target,"-k")==0 && index+1 < argc) {
	KAPPA = atof(argv[index+1]);
	rMessage("Kappa = " << KAPPA);
	index++;
	continue;
      }
      if (strcmp(target,"-pt")==0 && index+1 < argc) {
	int tmp = atoi(argv[index+1]);
	switch (tmp) {
	case 0:
	  parameterType = Parameter::PARAMETER_DEFAULT;
	  break;
	case 1:
	  parameterType = Parameter::PARAMETER_UNSTABLE_BUT_FAST;
	  break;
	case 2:
	  parameterType = Parameter::PARAMETER_STABLE_BUT_SLOW;
	  break;
	default:
	  parameterType = Parameter::PARAMETER_DEFAULT;
	}
	index++;
	paraFile = NULL;
	isParameter = false;
	continue;
      }
    }
  }
  else { // SDPA argument
    dataFile = argv[1];
    int len = strlen(dataFile);
    if (dataFile[len-1] == 's'
	&& dataFile[len-2] == '-') {
      isDataSparse = true;
    }
	
    outFile = argv[2];

    paraFile = (char *)"./param.sdpa";
    isParameter = true;

    for (int index=3; index<argc; ++index) {
      if (strcmp(argv[index],"-pt")==0 && index+1 < argc) {
	int tmp = atoi(argv[index+1]);
	switch (tmp) {
	case 0:
	  parameterType = Parameter::PARAMETER_DEFAULT;
	  break;
	case 1:
	  parameterType = Parameter::PARAMETER_UNSTABLE_BUT_FAST;
	  break;
	case 2:
	  parameterType = Parameter::PARAMETER_STABLE_BUT_SLOW;
	  break;
	default:
	  parameterType = Parameter::PARAMETER_DEFAULT;
	}
	index++;
	paraFile = NULL;
	isParameter = false;
      } // end of "-pt"
      else {
	initFile = argv[index];
	isInitFile = true;
	int len = strlen(initFile);
	if (initFile[len-1] == 's'
	    && initFile[len-2] == '-') {
	  isInitSparse = true;
	}
      }
    } // end of 'for'
    
  }
  
  if (dataFile == NULL || outFile == NULL) {
    message(argv[0]);
  }

  cout << "data      is " << dataFile;
  if (isDataSparse) {
    cout << " : sparse" << endl;
  } else {
    cout << " : dense" << endl;
  }
  if (paraFile) {
    cout << "parameter is " << paraFile << endl;
  }
  if (outFile) {
    cout << "out       is " << outFile <<endl;
  }
  if (initFile) {
    cout << "initial   is " << initFile;
  }
  if (isInitFile) {
    if (isInitSparse) {
      cout << " : sparse" << endl;
    } else {
      cout << " : dense" << endl;
    }
  } else {
    cout << endl;
  }
  if (paraFile == NULL) {
    if (parameterType == Parameter::PARAMETER_DEFAULT) {
      cout << "set       is DEFAULT" << endl;
    }
    else if (parameterType == Parameter::PARAMETER_UNSTABLE_BUT_FAST) {
      cout << "set       is UNSTABLE_BUT_FAST" << endl;
    }
    else if (parameterType == Parameter::PARAMETER_STABLE_BUT_SLOW) {
      cout << "set       is STABLE_BUT_SLOW" << endl;
    }
  }
  pinpal(dataFile, initFile, outFile, paraFile, isInitFile,
  	 isInitSparse, isDataSparse, isParameter,
	 parameterType, Display);
  return 0;
}

