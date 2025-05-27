#ifndef BOOLEAN_FUNCTION_H_INCLUDED
#define BOOLEAN_FUNCTION_H_INCLUDED

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
class boolfun {
public:
    int k;                      // Number of input variables
    int length;                 // Length of maptable
    double bias;                // Bias of maptable
    bool *ttt;                  // Bool Pointer of table 
    char Type;                  // Label of BF type
    int argus[4];               // 0,1,9 for ...; 0~2^15 for Cana.
    boolfun(){// default-value
        k=0; 
        length=1;
        bias=0.5;
        ttt=nullptr;
        Type='X';
        argus[0]=argus[1]=argus[2]=argus[3]=-1;}
    ~boolfun() {}
// Functions & Operations
    boolfun& Reset();                       // Reset.
    double NumSize1();                      // Count number of 1 in truth table.
    void PrintMappingTable();               // Print out the truth table.
    //boolfun& Configuration(char BF_Type, int indeg, double bias, bool *maps, ...);    // Configuration the parameters.
    boolfun& Configuration(char BF_Type, int indeg, double Bias, bool *maps, int par1,int par2,int par3,int par4);// Configuration the parameters.
    double Sensitivity();
    double Energy();
    double RelativeEnergy();
    boolfun& Gen_BF();                      // Generate a specified type Boolean function.
    bool is_PointedType(char types,bool Showit);
    std::vector<std::vector<int>> NestCana();

private:
// Special & Ordered Boolean functions (Generation)
    boolfun& Gen_Rand();
    boolfun& Gen_Cana();
    boolfun& Gen_Post();
    boolfun& Gen_Mone();
    boolfun& Gen_Spin();
    boolfun& Gen_Domi();
// Special & Ordered Boolean functions (Judgement)
    bool is_Cana();// Canalized function
    bool is_Post();// Post-2/Clique function
    bool is_Mone();// Monotonic function (strict constraint; all variables)
    bool is_Sign();// Signed function (Non-strict constraint; independent variables)
    bool is_Spin();// W_ij=±1, theta=0.
    bool is_Thre(bool PlusMinus,bool ShowIt);// W_ij, theta ∈R
    bool is_Effc();
    bool is_Domi();
// Some auxiliary functions
    bool aux_MoneSign(bool independent);// Monotone pattern.
    void aux_RBF(double p);// Return a Bernoulli distribution results with bias (p).
    void aux_ShuffleTTT();
};

// Generate a list of Boolean funcitons with pointed OBFs.
void BatchGenerationOBF(std::vector<std::vector<bool>> *BF_List, 
    int *RandIndex, int *InDegs, char OBF_Type, int SysSize, int part,
    double P_bias, double OBFv1, double OBFv2);
// Analyze possible polynomial expression of one Boolean function.
Rcpp::List PolynomialFunction(Rcpp::IntegerMatrix &VariableMat, Rcpp::IntegerVector &maptab,int IsSpinLike);
#endif // BOOLEAN_FUNCTION_H_INCLUDED
