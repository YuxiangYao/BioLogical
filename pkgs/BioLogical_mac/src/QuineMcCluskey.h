#ifndef MUL_QC_QUICK_H_INCLUDED
#define MUL_QC_QUICK_H_INCLUDED
#include <iostream>
#include <cstdio>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <cmath>
#include <set>
#include <vector>
#include <Rcpp.h>
// Quine-McCluskey method (Only consider Boolean system)
struct Prime{
    int id;
    int mapout;
    std::vector<int> slot;
    std::vector<int> blank;
    int cell_id;
    bool Covered;
    std::string Name;
    // Default
    Prime(): id(-1), mapout(-1),
        slot(0),blank(0),cell_id(-1),Covered(false),Name(""){}
    ~Prime() {}
    // Operator overloading.
    bool operator< (const Prime &a) const{
        return Name<a.Name;}
    bool operator== (const Prime &a) {
        if(Name==a.Name)return true;
        else return false;}
};

// Return number of invalid elements.
int qmc_aux_InvalidCount(std::vector<int> &Slot);
// Can pointed bit/dit of Prime A be covered by Prime B's one?
std::string qmc_aux_ExactMatchOnce(struct Prime &a,const struct Prime &b,int label);
// Prime A can be covered by Prime B but inversely.
std::string qmc_aux_ExactMatchOnceInverse(struct Prime &a,const struct Prime &b,int label);
// Return ID of a discrete vector belonging to which [cell].
int qmc_auc_Non0CellID(std::vector<int> &Counts,int *Order,int num);
// Count all non-zero numbers. [Boolean means "1"][Other means N"1" + N"2"+ ...]
int qmc_auc_Non0Summation(std::vector<int> &Counts, int Level);
// Can Prime A convered by Prime B?
bool qmc_aux_ConveredByPrime(struct Prime &a,struct Prime &b,int k);
// Self-addation operation in vector.
void qmc_aux_Adder1(std::vector <int> &num, int BitLeng, int Level);
// Count different elements in one input vector.
void qmc_aux_Counter(std::vector<int> &Seqs, std::vector<int> &Counter,int lengths);
// Change "the logic false -> true" of item named NNaammee.
void qmc_aux_ChangeLogic(std::set<struct Prime> &ss,std::string NNaammee);
// Remove invalid labels (count the number).
int qmc_aux_RemoveInvalidLabels(std::vector<int> &a,int target);
// Find the location or index. 
int qmc_aux_FirstFind(std::vector<int> &a,int target);

// Class of Quine-McCluskey Method
class QuineMcCluskeyMethod {
private:
    int k;// k-input function
    int L;// L-level discrete system (Boolean is 2)
    bool Incomparability;   // Each element is incomparable.
    int lengths;
    int *maptab;
    int *LL_system;
    struct Prime *implicants;                       // Original items
    std::vector<struct Prime> FinalPrime;           // Final primes should be presented.
    std::set<struct Prime> TmpImplicant;            // Used in each loop of checking prime items, current pending items.
    std::vector<std::set<struct Prime>> Records;    // Save the irreducible items, set[i]: remaining terms in i-th check.
    std::vector<std::set<struct Prime>> dicts;      // Save the sequencial primes in each repeat.
public:
    QuineMcCluskeyMethod(int a,int b,bool inc):k(a),L(b),Incomparability(inc),lengths(0),
        maptab(nullptr),LL_system(nullptr),implicants(nullptr) {}// Set the class.
    ~QuineMcCluskeyMethod(){}// Destructor
    // Functions & Operations
    QuineMcCluskeyMethod& Inits(int *p);
    QuineMcCluskeyMethod& Reset();          
    // void ShowImplicants();
    void ShowFinalPrime_B();// Boolean 
    Rcpp::IntegerMatrix ShowImplicants_M();// multi-valued
    std::vector<std::vector<int>> ReturnFinalPrimeInfo();
    QuineMcCluskeyMethod& Do_QuineMcCluskey();
    double Effective();
    int NumberFinalPrime();
    void SingleEdgeConnect_bool(std::vector<double> &WeightConnect);
private:
    void BuildCandidate(int Cell_ID,std::set<struct Prime> &a);
    void HyperCube();
    void PrepareAnalysis();
    std::vector<int> AnalysisMode(struct Prime item,std::set<struct Prime> &a);
    std::vector<std::vector<std::string>> Inverse_AnalysisMode(struct Prime item,
        std::set<struct Prime> &a,std::vector<int> &mapto,std::vector<int> &WhichVar);
    int RecursionImplicant();
    void RecursionProcess();
    std::vector<struct Prime> SetCoveredMatrix();
};

// Get the effective connections of nodes.
double OnceEffective(int *maptab,int k, int L);
double toR_MultipleEffective(int *maptab,int k,int L,int lens);
// Show the disjunctive normal form of Booelan function.
void toR_ShowBoolFunDNF(int *maptab,int k,Rcpp::CharacterVector VarNames);
// Analyze each connect effective. (Boolean or multiple)
void OnceEdgeConnect(int *maptab, int k, int L, std::vector<double> &Slot);
std::vector<double> toR_MultipleEdgeConnect(int *maptab,int k,int L, int lens);
// The complexity of ternary/multiple expressions (Independent & Comparable)
double toR_BoolMulComplexity(int *maptab,int k,int L,int lens, bool incomparable);

Rcpp::IntegerMatrix toR_ShowMulVFunDNF(int *maptab,int k, int L);

#endif // MUL_QC_QUICK_H_INCLUDED
