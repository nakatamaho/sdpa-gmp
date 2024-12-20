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

#ifndef __sdpa_parts_h__
#define __sdpa_parts_h__

#include <sdpa_newton.h>

namespace sdpa {

class Newton;

class Solutions;
class InputData;
class Residuals;
class WorkVariables;

class ComputeTime;
class Parameter;
class StepLength;
class DirectionParameter;
class Switch;
class RatioInitResCurrentRes;
class SolveInfo;
class Phase;
class AverageComplementarity;

class ComputeTime {
  public:
    double Predictor;
    double Corrector;
    double StepPredictor;
    double StepCorrector;
    double xMatTime;
    double zMatTime;
    double invzMatTime;
    double xMatzMatTime;
    double EigxMatTime;
    double EigzMatTime;
    double EigxMatzMatTime;
    double makerMat;
    double makebMat;
    double B_DIAG;
    double B_F1;
    double B_F2;
    double B_F3;
    double B_PRE;
    double makegVecMul;
    double makegVec;
    double choleskybMat;
    double solve;
    double sumDz;
    double makedX;
    double symmetriseDx;
    double makedXdZ;
    double updateRes;
    double MainLoop;
    double FileRead;
    double FileCheck;
    double FileChange;
    double TotalTime;
    ComputeTime();
    ~ComputeTime();
    void display(FILE *fpout = stdout);
};

class Parameter {
  public:
    enum parameterType { PARAMETER_DEFAULT, PARAMETER_UNSTABLE_BUT_FAST, PARAMETER_STABLE_BUT_SLOW };
    int maxIteration;
    int precision;
    double epsilonStar;
    double lambdaStar;
    double omegaStar;
    double lowerBound;
    double upperBound;
    double betaStar;
    double betaBar;
    double gammaStar;
    double epsilonDash;
#define PRINT_DEFAULT_LENGTH 30
    static char xPRINT_DEFAULT[PRINT_DEFAULT_LENGTH];
    static char XPRINT_DEFAULT[PRINT_DEFAULT_LENGTH];
    static char YPRINT_DEFAULT[PRINT_DEFAULT_LENGTH];
    static char infPRINT_DEFAULT[PRINT_DEFAULT_LENGTH];
    char xPrint[PRINT_DEFAULT_LENGTH];
    char XPrint[PRINT_DEFAULT_LENGTH];
    char YPrint[PRINT_DEFAULT_LENGTH];
    char infPrint[PRINT_DEFAULT_LENGTH];
    Parameter();
    Parameter(FILE *parameterFile);
    ~Parameter();
    void setDefaultParameter(parameterType type = PARAMETER_DEFAULT);
    void readFile(FILE *parameterFile);
    void display(FILE *fpout = stdout, const char *printFormat = infPRINT_DEFAULT);
};

class StepLength {
  public:
    mpf_class primal;
    mpf_class dual;
    StepLength();
    StepLength(mpf_class alphaP, mpf_class alphaD, int nBlock, int *blockStruct);
    ~StepLength();
    void initialize(mpf_class alphaP, mpf_class alphaD);
    void terminate();

    static mpf_class minBlockVector(BlockVector &aVec);

    void computeStepLength(Solutions &currentPt, Newton &newton, WorkVariables &work, ComputeTime &com);
    void MehrotraPredictor(InputData &inputData, Solutions &currentPt, Phase &phase, Newton &newton, WorkVariables &work, ComputeTime &com);
    void MehrotraCorrector(InputData &inputData, Solutions &currentPt, Phase &phase, Switch &reduction, Newton &newton, AverageComplementarity &mu, RatioInitResCurrentRes &theta, WorkVariables &work, Parameter &param, ComputeTime &com);
    void display(FILE *fpout = stdout);
};

class DirectionParameter {
  public:
    mpf_class value;
    DirectionParameter(mpf_class betaStar = 0.0);
    ~DirectionParameter();
    void initialize(mpf_class betaStar = 0.0);

    void MehrotraPredictor(Phase &phase, Switch &reduction, Parameter &param);
    void MehrotraCorrector(Phase &phase, StepLength &alpha, Solutions &currentPt, Newton &newton, AverageComplementarity &mu, Parameter &param);
    void display(FILE *fpout = stdout);
};

class Switch {
  public:
    enum SwitchType { ON, OFF };
    SwitchType switchType;

    Switch(SwitchType switchType = ON);
    ~Switch();
    void initialize(SwitchType switchType = ON);

    void MehrotraPredictor(Phase &phase);
    void display(FILE *fpout = stdout);
};

class AverageComplementarity {
  public:
    mpf_class initial;
    mpf_class current;
    AverageComplementarity(mpf_class lambdaStar = 0.0);
    ~AverageComplementarity();
    void initialize(mpf_class lambdaStar = 0.0);
    void initialize(Solutions &initPt);
    void update(Solutions &currentPt);
    void display(FILE *fpout = stdout);
};

class RatioInitResCurrentRes {
  public:
    mpf_class primal;
    mpf_class dual;

    RatioInitResCurrentRes();
    RatioInitResCurrentRes(Parameter &param, Residuals &initRes);
    ~RatioInitResCurrentRes();

    void initialize(Parameter &param, Residuals &initRes);

    void update(Switch &reduction, StepLength &alpha);
    void update_exact(Residuals &initRes, Residuals &currentRes);
    void display(FILE *fpout = stdout);
};

class SolveInfo {
  public:
    enum phaseType { noINFO, pFEAS, dFEAS, pdFEAS, pdINF, pFEAS_dINF, pINF_dFEAS, pdOPT, pUNBD, dUNBD };

    mpf_class rho;
    mpf_class etaPrimal;
    mpf_class etaDual;
    mpf_class objValPrimal;
    mpf_class objValDual;

    SolveInfo();
    SolveInfo(InputData &inputData, Solutions &currentPt, mpf_class mu0, mpf_class omegaStar);
    ~SolveInfo();

    void initialize(InputData &inputData, Solutions &currentPt, mpf_class mu0, mpf_class omegaStar);

    void update(InputData &inputData, DenseLinearSpace &initPt_xMat, DenseLinearSpace &initPt_zMat, Solutions &currentPt, Residuals &currentRes, AverageComplementarity &mu, RatioInitResCurrentRes &theta, Parameter &param);
    // assume   initPt = lambda * I
    void update(double &lambda, InputData &inputData, Solutions &currentPt, Residuals &currentRes, AverageComplementarity &mu, RatioInitResCurrentRes &theta, Parameter &param);
    // check mu, gap, feasibility   2007/08/27
    void check(InputData &inputData, Solutions &currentPt, Residuals &currentRes, AverageComplementarity &mu, RatioInitResCurrentRes &theta, Parameter &param);
    void display(FILE *fpout = stdout);
};

class Phase {
  public:
    int nDim;
    SolveInfo::phaseType value;

    Phase();
    Phase(Residuals &initRes, SolveInfo &solveInfo, Parameter &param, int nDim);
    ~Phase();

    bool initialize(Residuals &initRes, SolveInfo &solveInfo, Parameter &param, int nDim);
    bool updateCheck(Residuals &currentRes, SolveInfo &solveInfo, Parameter &param);
    void reverse();
    void display(FILE *fpout = stdout);
};

} // namespace sdpa

#endif // __sdpa_parts_h__
