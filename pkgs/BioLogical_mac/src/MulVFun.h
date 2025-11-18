/*
  Copyright (c) 2023-2024, YAO Yuxiang
    date: 2024-01-16
    version: 1.0

   Generate multi-valued function.
*/
#ifndef MULTIVALUED_FUNCTION_H_INCLUDED
#define MULTIVALUED_FUNCTION_H_INCLUDED

#include <iostream>
#include <cstdio>
#include <vector>
#include <cstring>
#include <string>
#include <algorithm>
#include <cstdarg>
#include <z3++.h>
//#include "RandomSetting.h"
#include <Rcpp.h>
class mulvfun {
private:
    int k;                      // Number of input variables
    int L;                      // L-level discrete value
    int length;      // Length of maptable
public:
    char Type;                  // Label of BF type
    double *bias;               // Bias of maptable (Only use external data, no need input)
    short *ttt;                 // Pointer of table (Only use external data, no need input)
    std::vector<int> argus;
    mulvfun(int a,int b,int c,char type):k(a),L(b),length(c),Type(type), 
        bias(nullptr),ttt(nullptr) {}// Set the class.
    ~mulvfun(){}// Destructor

// Functions & Operations
    mulvfun& Reset();                       // Reset.
    mulvfun& Configuration(double *BiasSet, short *maps,
        std::vector<int> &pars);// Configuration the parameters.
    
// Special & Ordered Boolean functions (Generation)
    mulvfun& Gen_Rand();
    mulvfun& Gen_Cana_Free(Rcpp::IntegerVector &CanaVarID, Rcpp::IntegerVector &EachCanaNum);
    mulvfun& Gen_Cana_Config(Rcpp::IntegerVector &CanaVarID, 
        std::vector<int> &EachCanaNum, 
        std::vector<std::vector<short>> &CanaIns, 
        std::vector<std::vector<short>> &CanaOut);
    std::vector<std::vector<std::vector<short>>> NestCana(
        std::vector<int> &CanalizingIDs, bool Quick);

    mulvfun& Gen_MulVF_Thre();// Multi-valued linear threshold function
    Rcpp::List is_MulVF_Thre();// Judge it is a multi-valued linear threshold function
    mulvfun& Gen_MulVF_Domi();// Multi-valued dominant function
    Rcpp::List is_MulVF_Domi();// Judge it is a multi-valued dominant function
    
    Rcpp::List is_MulVF_Sign();// Judge it is a multi-valued broadly defined-signed function
    double Sensitivity();// Function's sensitivity

private:
    int OneVariableActivity(std::vector<unsigned short> &Vec, int &VarID, std::vector<unsigned int> &LL_system);
    std::vector<std::vector<short>> OnceCanalizedChecker(int PointedID, std::vector<bool> &UnCana);
    short FindDomiantComponent(std::vector<unsigned short> &vec1,
        std::vector<int> &wights, std::vector<double> &theta);
    
// Some auxiliary functions
    void aux_RMF();// Return a Bernoulli distribution results with bias (p).
};

// Some auxiliary functions
void MulVFun_Adder(std::vector<unsigned short> &Slots, int k, int L);
void Replaced_Shuffle(std::vector<short> &aVec);
int VecDotTimes(std::vector<unsigned short> &Vec, std::vector<unsigned int> &LL_system, int lens);
int VecDotTimes_pn(std::vector<unsigned short> &vec1, std::vector<int> &vec2, int lens);
bool AnyStateTRUE(std::vector<bool> &vec);

Rcpp::List MulVFun_PolynomialFunction(Rcpp::IntegerMatrix &VariableMat, Rcpp::IntegerVector &maptab, int k, int L_sys);


double MeansOfSet(std::vector<int> &xx);
std::vector<int> MulV2Bool(std::vector<int> &mulvfun, int k, int L, std::vector<int> &threshold);
std::vector<int> Bool2MulV(std::vector<int> &boolfun, int k, int L, std::vector<int> &threshold);

// Generate a list of Boolean funcitons with pointed OBFs.
// void BatchGenerationOBF(std::vector<std::vector<bool>> *BF_List, 
//     int *RandIndex, int *InDegs, char OBF_Type, int SysSize, int part,
//     double P_bias, double OBFv1, double OBFv2);
// Analyze possible polynomial expression of one Boolean function.
//Rcpp::List PolynomialFunction(Rcpp::IntegerMatrix &VariableMat, Rcpp::IntegerVector &maptab,int IsSpinLike);
#endif // MULTIVALUED_FUNCTION_H_INCLUDED
