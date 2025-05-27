#ifndef NETWORK_GRAPHTHEORY_H_INCLUDED
#define NETWORK_GRAPHTHEORY_H_INCLUDED
#include <stdio.h>
#include <stdlib.h>// for itoa
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cstdarg>
#include <vector>
#include <Rcpp.h>

// Structure of regulations.
struct Regulation{
    int code;// Node code.
    struct Regulation *prev;// Preceder.
    struct Regulation *next;// Successor.
};
void net_chooseKfromN(int N,int K,int *Selected,int *Label,int flag);
// Insert a parent-type node.
void InsertParentsNode(struct Regulation *p, int parent);
// Insert a child-type node.
void InsertChildsNode(struct Regulation *p,int child);
// Build a regulation edge/link.
void BuildRegulationship(struct Regulation *net,int child,int parent);
// Build non-directional edge (Only use child, not contain the parent-type)
void BuildNeighborEdge(struct Regulation *net,int node1,int node2);
// Remove one specified child from parent node.
int DeleteChilds(struct Regulation *p,int child);
// Remove one specified parent from child node.
int DeleteParents(struct Regulation *p,int parent);
// Node has No.xx child/parent.
int ExistChild(struct Regulation *p,int child);
int ExistParent(struct Regulation *p,int parent);
// Build a mutual loop of all nodes.
void LinkMutualLoop(struct Regulation *Loops,int total);
void LinkRomoveCurrent(struct Regulation *rpt);
// First-in first-out to queue of distance.
void ShortestFind_FIFO(struct Regulation *ParentNode,int Code,int logi);
// Replaced Poisson distribution
int Replaced_Poisson(double lambda);

// Define network class.
class NetGraphFrame {
public:
    int *InDeg=nullptr;
    int *OtDeg=nullptr;
    struct Regulation *Network=nullptr; // Slot for edge's relations.
private:
    int total;              // Size of network
    char Type;
    int *Coding=nullptr;    // Code: 0,1,2,3,4,...
    int **Address=nullptr;  // Pointing address of [code] (Double of total+1).
    std::string NetPara;
    int tmp_i_arg[4];
    double tmp_f_arg[4];
public:
    NetGraphFrame(){// default-value
        total=-1;
        Type='X';
        InDeg=nullptr;
        OtDeg=nullptr;
        Network=nullptr;
        Coding=nullptr;
        Address=nullptr;
        NetPara.clear();}
    ~NetGraphFrame(){
        NetGraphFrame::Reset_Network();// First reset.
        free(InDeg);
        InDeg=OtDeg=Coding=nullptr;
        free(Network);
        Network=nullptr;
        free(Address);
        Address=nullptr;
        NetPara.clear();
    }
    NetGraphFrame& ConfigurationBuildNet(char Net_Type, int Size, int ipar1, int ipar2, double fpar1);
    //NetGraphFrame& Delete_Network();
    NetGraphFrame& Build_GraphNet();
    NetGraphFrame& Out2VecVecIntFrame(
        std::vector<std::vector<int>> *InEdge, std::vector<std::vector<int>> *OutEdge);
    NetGraphFrame& LoadFromVecVecIntFrame(
        std::vector<std::vector<int>> &InEdge, std::vector<std::vector<int>> &OutEdge);
    NetGraphFrame& IsolatedPointedNode_D(int ID);
    NetGraphFrame& BreakDownPointedEdge_D(int source, int target);
    NetGraphFrame& Tarjon(std::vector<std::vector<int>> &SCC_List);
    void Count_InOt_Deg();
    void Show_Network();
    void Show_Degree();
    void ShortestPath(int source);
private:
    void Initial_Network();
    void Reset_Network();
// Generate various types of network
    void Net_ER_D(double AvDeg);
    void Net_SF_U(int MaxDeg,int MinDeg,double Gamma);
    void Net_BA_U(int Core);
    void Net_RR_D(int K);
    void Net_LT_D(int K);
    void Net_NK_D(int K);
// Some auxiliary functions
    void Address_Reset();// Reset temporary address slot.
    void Auxiliary_ER_D_sub(int *OtLab,int child,int parent);
    void Auxiliary_ER_D();
    void Auxiliary_SF_U();
    void Auxiliary_RR_D();
    void Auxiliary_RR_D_sub(int *OtLab,int child,int parent);
    void Auxiliary_LT_D(int ids,int *Neighbor);
    void Auxiliary_NetworkTarjan(int ThisID,std::vector<bool> &hold_on, std::vector<int> &finding,
        std::vector<int> &lowrank, std::vector<int> &sequnces, std::vector<std::vector<int>> &scc,
        int &step, int &index);
};

#endif // NETWORK_GRAPHTHEORY_H_INCLUDED
