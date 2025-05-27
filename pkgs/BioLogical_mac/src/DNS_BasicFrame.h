#ifndef DNS_BASICFRAME_H_INCLUDED
#define DNS_BASICFRAME_H_INCLUDED
#include "BoolFun.h"
#include "MulVFun.h"
#include "NetGraphFrame.h"
#include <set>
#include <Rcpp.h>
#include <cmath>

void DNS_Aux_GenerateTopo(char NetType, double NetPara, int Size, 
    std::vector<std::vector<int>> *inedge,// [!!warning] Memory allocated over.
    std::vector<std::vector<int>> *otedge, int *InDegs);

void DNS_Aux_GenerateFunc(int sys_size, int ll_sys, int *InDeg, int *IDsss,
    double *bias, char of_type, double of_part,
    int of_i_para1, int of_i_para2, std::vector<std::vector<int>> *bnsss,
    std::vector<int> &ControlNodes, std::vector<int> &ControlValues);

void DNS_Aux_Bias(short *slot, int length, double *p, int Ls);


// Class of Boolean Network System basic framework.
class DNS_Basic {
protected:
    int N=0; // System size
    int L=2; // Discrete level
    std::vector<std::vector<int>> *inedge=nullptr;      // Info of in-degree
    std::vector<std::vector<int>> *otedge=nullptr;      // Info of out-degree
    std::vector<std::vector<int>> *bns=nullptr;         // Info of function
    short *sss=nullptr;     // State slot, double size for analyzing.
    int *LL_sys=nullptr;    // order: 1, L, L^2, ... L^12
    int i_argus[16];        
    double f_argus[16];     
public:
    int *Labels=nullptr;    // Label: all zeros to control
    DNS_Basic();
    DNS_Basic(int size, int L);
    ~DNS_Basic();
    DNS_Basic& SystemEvolution_Syn(int steps, short *ss0, std::vector<bool> *flag);  // Synchronous update
    DNS_Basic& SystemEvolution_Asy1(int steps, short *ss0, std::vector<bool> *flag); // Allow invalid update
    DNS_Basic& SystemEvolution_Asy2(int steps, short *ss0, std::vector<bool> *flag); // Must select one update
    DNS_Basic& LoadedModel(std::vector<std::vector<int>> *load_InDeg, 
      std::vector<std::vector<int>> *load_OtDeg, std::vector<std::vector<int>> *Bnsss);
protected:
    
};

// Derived classes: Derrida damage spread
class DNS_Derrida: public DNS_Basic {
public:
    DNS_Derrida();
    DNS_Derrida(int size, int L);
    ~DNS_Derrida();
    DNS_Derrida& DynamicParaConfig(double *s_bias,double s_delta);// Set dynamic parameter, can be used repeatedly.
    DNS_Derrida& DerridaDamageSpread(int Steps, int UpdateType);
    double FinalDistance();
};

// Derived classes: Percolation & Phase transition
class DNS_Percolation: public DNS_Basic {
protected:
    int Lattice=4;
    std::vector<bool> *flag=nullptr;
    std::vector<int> MaxCluster_ID;
public:
    DNS_Percolation();
    DNS_Percolation(int size, int L, int LaType);
    ~DNS_Percolation();
    DNS_Percolation& DynamicParaConfig(double *s_bias);
    DNS_Percolation& PercolationModel(int Steps, int Windows, int UpdateType);
    DNS_Percolation& PercolationWhetherNot();
    DNS_Percolation& OutputFinalLattice(int *IsFixed, int *MaxCluster);
    DNS_Percolation& PercolationStableFraction(double *Returner);
};

// Derived classes: Scaling pattern.
class DNS_Engaged: public DNS_Basic {
protected:
    // int NumSys=2;               // The system of number 
    // int Order[16];              // Number system
    int *PossibleLocal=nullptr; // Remain possible locals of unstable values
    int *PossibleValue=nullptr; // Remain possible the number of discrete values
    int *Pruned=nullptr;        // Label: Whether the node is invalid due to no useful output instead of stable.public: // Scaling Pattern >>>
    int **Mapping=nullptr;      // Transient mapping tables [dynamical]
    int *InDeg_temp=nullptr;    // Transient in-degree [dynamical]
    int *OtDeg_temp=nullptr;    // Transient out-degree [dynamical]
    int **Parents_temp=nullptr; // Transient parents [dynamical]
    std::vector<std::set<int>> *Children_temp=nullptr;
    std::set<int> GlobalCandidate;
    std::set<int> StableNode;
    std::set<int> UselessNode;
public: 
    DNS_Engaged();
    DNS_Engaged(int size, int L,
        std::vector<std::vector<int>> *ipt1, 
        std::vector<std::vector<int>> *ipt2, 
        std::vector<std::vector<int>> *ipt3);
    ~DNS_Engaged();
    void OnlyScalingPattern(int *OutputVec, int *SUW);
    DNS_Engaged& Export2ResidualNetwork(std::vector<std::vector<int>> &Exp_InEdge,
        std::vector<std::vector<int>> &Exp_OtEdge, std::vector<std::vector<int>> &BoolFn);
    std::vector<int> Export2AllNodeState();// Whether nodes are stable (1 or 0) or not?
    std::vector<int> Export2AllNodeType();// Nodes belong to stable, useless, or engaged?
protected:
    int NodeShrunk(int Code);
    void RemoveInvalidInputs(int Code);
    DNS_Engaged& ClampedVertex();// Initial & subsequent clamping
    DNS_Engaged& PrunedVertex();// Initial & subsequent pruning
    DNS_Engaged& CuttingLeavesVertex();// Cutting all stable leaf-nodes (Part A of RemoveStable.....)
};

class DNS_CoreDyn: public DNS_Engaged {
public: // Now only suit for Boolean system
    DNS_CoreDyn();
    DNS_CoreDyn(int size, int L,
        std::vector<std::vector<int>> *ipt1, 
        std::vector<std::vector<int>> *ipt2, 
        std::vector<std::vector<int>> *ipt3);
    ~DNS_CoreDyn();
    void OnceCoreDynamic(int *OutputVec,int Times);
protected:
    DNS_CoreDyn& RemoveStableNonZeroOutput2();// Simplify topology 
    DNS_CoreDyn& SingleInputCouple(int IDxx);
    DNS_CoreDyn& DoubleInputCoupleType01(int IDxx, int Successor, int Par_H, int Par_L);
    DNS_CoreDyn& DoubleInputCoupleType02(int IDxx, int Successor, int Par_H, int Par_L);
    DNS_CoreDyn& DoubleInputCouple(int IDxx);
};

#endif // BNS_BASICFRAME_H_INCLUDED