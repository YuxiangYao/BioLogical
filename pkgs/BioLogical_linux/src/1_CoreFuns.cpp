#include <Rcpp.h>
#include "BoolFun.h"
#include "QuineMcCluskey.h"
#include "NetGraphFrame.h"
#include "DNS_BasicFrame.h"
#include "MulVFun.h"
#include <map>

using namespace Rcpp;

// Inner function: Load networks from R space. (User invisible)
void c_inner_Robj2Cobj(Rcpp::List &aRealNet, IntegerVector &PointedGene, IntegerVector &PointValues,
    IntegerVector &InD, IntegerVector &OtD,
    std::vector<std::vector<int>> &ot_edge,// Outer1 
    std::vector<std::vector<int>> &in_edge,// Outer2
    std::vector<std::vector<int>> &bn_list,
    IntegerVector &Len5Int32, bool AlsoBN, int ll_sys){// Outer3
    // 3: External (No input node), 4: Terminal (No output node)
    Len5Int32[3]=0;Len5Int32[4]=0;
    // Convert List into appropriate data_frame.
    Rcpp::IntegerVector tmpdata;
    Rcpp::List tmplist1=aRealNet[1], tmplist2=aRealNet[2], tmplist3=aRealNet[3];
    int sys_size=tmplist1.size();
    in_edge.resize(sys_size);
    ot_edge.resize(sys_size);
    if(AlsoBN)bn_list.resize(sys_size);
    for(int ii=0; ii<sys_size; ++ii){
        in_edge[ii].resize(0);
        ot_edge[ii].resize(0);}
    if(AlsoBN)for(int ii=0; ii<sys_size; ++ii)bn_list[ii].resize(0);
    for(int ii=0; ii<sys_size; ++ii){
        if(OtD[ii]>0){// Outs not empty // if(Rcpp::internal::Rcpp_IsNA(tmpdata[0]))
            tmpdata=tmplist2[ii];
            for(int jj=0; jj<tmpdata.size(); ++jj){
                ot_edge[ii].push_back(tmpdata[jj]);}}
        else {
            Len5Int32[4]=Len5Int32[4]+1;}
        if(InD[ii]>0){// Ins not empty // if(Rcpp::internal::Rcpp_IsNA(tmpdata[0]))
            tmpdata=tmplist1[ii];
            for(int jj=0; jj<tmpdata.size(); ++jj){
                in_edge[ii].push_back(tmpdata[jj]);}
            if(AlsoBN){
                tmpdata=tmplist3[ii];
                for(int jj=0; jj<tmpdata.size(); ++jj){
                    bn_list[ii].push_back(tmpdata[jj]);}}
        }
        else {// No input, should be set one
            if(AlsoBN)bn_list[ii].push_back(PointValues[ii]);// Actually this value could be from Exponents or Controller.
            Len5Int32[3]=Len5Int32[3]+1;}}
    // Check whether exist non-exponenet controller?
    if(AlsoBN){
        if(PointedGene[0]>=0){// Exist!! Should reset all corresponding mapping tables.
            int po_id, po_va;
            for(int ii=0; ii<PointedGene.size(); ++ii){
                po_id=PointedGene[ii];
                po_va=PointValues[po_id];
                //std::vector<int> tmps(1<<InD[po_id],po_va);
                std::vector<int> tmps((int)(std::pow(ll_sys,InD[po_id])),po_va);
                bn_list[po_id]=tmps;}}}  
    Len5Int32[0]=sys_size;
} 

// Inner funciton: format conversion from C++ to R. (User invisible)
Rcpp::List Cobj2Robj_Net(int SysSize, std::vector<std::vector<int>> *IN_edge, 
    std::vector<std::vector<int>> *OT_edge, std::vector<std::vector<int>> *BN_list){
    Rcpp::List Returner;
    Rcpp::List tmp1,tmp2,tmp3;
    Returner.push_back(NA_INTEGER);
    for(int ii=0; ii<SysSize; ++ii){
        tmp1.push_back(
            Rcpp::IntegerVector((*IN_edge)[ii].begin(), (*IN_edge)[ii].end()) );
        tmp2.push_back(
            Rcpp::IntegerVector((*OT_edge)[ii].begin(), (*OT_edge)[ii].end()) );
        tmp3.push_back(
            Rcpp::IntegerVector((*BN_list)[ii].begin(), (*BN_list)[ii].end()) );}
    Returner.push_back(tmp1);
    Returner.push_back(tmp2);
    Returner.push_back(tmp3);
    return Returner;
}

// Inner function: Random generation

//' @title Determine whether a function belongs to a special type
//' @description C++ prototype function of \link{BoolFun_Type}.
//' @param boolfuns [\code{LogicalVector}] the Boolean function logical vector.
//' @param k [\code{int}] the number of input variables of the function.
//' @param Type [\code{char}] function type.
//' @param Showit [\code{bool}] whether show the function?
//[[Rcpp::export]]
bool c_BF_isPointed(LogicalVector boolfuns,int k,char Type,bool Showit){
    boolfun abf;
    bool *bitmap=(bool*)malloc((1<<k)*sizeof(bool));
    for(int ii=0;ii<(1<<k);ii++){
        bitmap[ii]=boolfuns[ii];}
    bool xx=false;
    xx=abf.Configuration('X',k,0.5,bitmap,-1,-1,-1,-1).is_PointedType(Type,Showit);
    free(bitmap);
    return xx;
}

//' @title Generate a specific type of Boolean function
//' @description C++ prototype function of \link{BoolFun_Generator}.
//' @param Type [\code{char}] function type.
//' @param k [\code{k}] the number of input variables of the function.
//' @param bias [\code{double}] function bias.
//' @param Vars [\code{IntegerVector}] some configuration parameters.
//[[Rcpp::export]]
Rcpp::IntegerVector c_BF_Generator(char Type, int k, double bias, IntegerVector Vars){
    boolfun abf;
    bool *bitmap=(bool*)malloc((1<<k)*sizeof(bool));
    abf.Configuration(Type,k,bias,bitmap,Vars[0],Vars[1],Vars[2],Vars[3]).Gen_BF();
    Rcpp::IntegerVector xx(1<<k);
    for(int ii=0; ii<(1<<k); ii++){
        xx[ii]=bitmap[ii]?1:0;}
    free(bitmap);
    return xx;
}

//' @title Calculate the sensitivity of a Boolean function
//' @description C++ prototype function of \link{BoolFun_Sensitivity}.
//' @param boolfuns [\code{LogicalVector &}] the Boolean function logical vector.
//' @param k [\code{int}] the number of input variables of the function.
//[[Rcpp::export]]
double c_BF_Sensitivity(LogicalVector &boolfuns, int k){
    boolfun abf;
    bool *bitmap=(bool*)malloc((1<<k)*sizeof(bool));
    for(int ii=0; ii<(1<<k); ++ii){
        bitmap[ii]=boolfuns[ii];}
    double xx=abf.Configuration('X',k,0.5,bitmap,-1,-1,-1,-1).Sensitivity();
    free(bitmap);
    return xx;
}

//' @title Calculate the effective edges of a Boolean function
//' @description C++ prototype function of \link{BoolFun_Effectiveness}.
//' @param boolfuns [\code{IntegerVector &}] the Boolean function.
//' @param k [\code{int}] the number of input variables of the function.
//[[Rcpp::export]]
double c_BF_Effective(IntegerVector &boolfuns,int k){
    bool constant=true;
    double xx=0;
    int *ttt=boolfuns.begin();
    int First=boolfuns[0];
    for(int ii=0; ii<(1<<k); ++ii){
        constant=constant&&(First==boolfuns[ii]);
        if(!constant)break;}
    if(!constant){
        xx=toR_MultipleEffective(ttt, k, 2, (1<<k));}
    return xx;
}

//' @title Calculate the complexity of a Boolean function
//' @description C++ prototype function of \link{BoolFun_Complexity}.
//' @param boolfuns [\code{IntegerVector &}] the Boolean function.
//' @param k [\code{int}] the number of input variables of the function.
//[[Rcpp::export]]
double c_BF_Complexity(IntegerVector &boolfuns, int k){
    int *ttt=boolfuns.begin();
    double xx=toR_BoolMulComplexity(ttt, k, 2, (1<<k), true);// Fixed L:2 , lens:1<<k
    return xx;
}

//' @title Calculate the effective edges of a Boolean function
//' @description C++ prototype function of \link{BoolFun_Effectiveness}.
//' @param boolfuns [\code{IntegerVector &}] the Boolean function.
//' @param k [\code{int}] the number of input variables of the function.
//[[Rcpp::export]]
Rcpp::NumericVector c_BF_EffectiveEdges(IntegerVector &boolfuns, int k){
    int *ttt=boolfuns.begin();
    int First=boolfuns[0];
    bool constant=true;
    std::vector<double> tmp_res;
    for(int ii=0; ii<(1<<k); ++ii){
        constant=constant&&(First==boolfuns[ii]);
        if(!constant)break;}
    if(!constant){
        tmp_res=toR_MultipleEdgeConnect(ttt, k, 2, (1<<k));}
    else {
        tmp_res.resize(k,0);}
    Rcpp::NumericVector Res=Rcpp::wrap(tmp_res);
    return Res;
}

//' @title Display the Quine-McCluskey form of a Boolean function
//' @description C++ prototype function of \link{BoolFun_QMForm}.
//' @param boolfuns [\code{IntegerVector &}] the Boolean function.
//' @param k [\code{int}] the number of input variables of the function.
//' @param VarsName [\code{CharacterVector}] variable names.
//[[Rcpp::export]]
void c_BF_QuineMcCluskey(IntegerVector &boolfuns, int k, CharacterVector VarsName){
    int *ttt=boolfuns.begin();
    toR_ShowBoolFunDNF(ttt,k,VarsName);
    
}

//' @title Enumerate all possible input vectors
//' @description C++ prototype function for enumerate all possible input vectors
//' under the given length.
//' @param VarNum [\code{int}] the number of input variables of the function
//[[Rcpp::export]]
Rcpp::IntegerMatrix c_FrameTruthTable(int VarNum){
    int length=1<<VarNum;
    int ii,jj,sumer;
    Rcpp::IntegerMatrix res(length, VarNum);
    for(ii=0; ii<length; ii++){
        sumer=ii;
        for(jj=0; jj<VarNum; jj++){
            res(ii,VarNum-jj-1)=sumer&1;
            sumer>>=1;}}
    return (res);
}

//' @title Decode the structure of a nested canalized Boolean function
//' @description C++ prototype function of \link{BoolFun_NestedCana}.
//' @param boolfuns [\code{IntegerVector &}] the Boolean function logical vector.
//' @param k [\code{int}] the number of input variables of the function.
//[[Rcpp::export]]
Rcpp::List c_B_NestedCanalized(IntegerVector &boolfuns,int k){
    boolfun abf;
    bool *bitmap=(bool*)malloc((1<<k)*sizeof(bool));
    for(int ii=0;ii<(1<<k);ii++){
        bitmap[ii]=boolfuns[ii];}
    std::vector<std::vector<int>> xx;
    xx=abf.Configuration('X',k,0.5,bitmap,-1,-1,-1,-1).NestCana();
    free(bitmap);
    Rcpp::List Res;
    Rcpp::IntegerVector id1(xx[0].begin(), xx[0].end());Res.push_back(id1);
    Rcpp::IntegerVector id2(xx[1].begin(), xx[1].end());Res.push_back(id2);
    Rcpp::IntegerVector id3(xx[2].begin(), xx[2].end());Res.push_back(id3);
    return (Res);
}

//' @title Analyze damage spread within system
//' @description C++ prototype function of \link{DNS_DamageSpread}
//' @param sys_size [\code{int}] system size
//' @param ll_system [\code{int}] level number of discrete system
//' @param sim_step [\code{int}] simulating steps
//' @param obf_type [\code{char}] ordered Boolean function type
//' @param bias_rf [\code{NumericVector}] bias in random function 
//' @param obf_ratio [\code{double}] proportion of ordered Boolean function 
//' @param init_dis [\code{double}] initial normalized Hamming distance
//' @param init_1_ratio [\code{NumericVector}] "1" proportion of initial state
//' @param net_type [\code{char}] system topological type
//' @param net_f_para [\code{double}] topological configured parameters
//' @param obf_i_para1 [\code{int}] configuration parameter for OBF (1)
//' @param obf_i_para2 [\code{int}] configuration parameter for OBF (2)
//' @param RuleType [\code{int}] update strategy
//[[Rcpp::export]]
double c_Derrida_Simualtion(int sys_size, int ll_system, int sim_step, char obf_type, 
    Rcpp::NumericVector bias_rf, double obf_ratio,
    double init_dis, Rcpp::NumericVector init_1_ratio, char net_type, // Only simulate pseudo-networks.
    double net_f_para, int obf_i_para1, int obf_i_para2, int RuleType){
    // Set a pseudo-model.
    std::vector<int> con_ids, con_val; //con_id.clear(); con_val.clear();
    std::vector<std::vector<int>> in_edge(sys_size);
    for(int ii=0; ii<sys_size; ++ii){
        in_edge[ii].resize(0);}
    std::vector<std::vector<int>> ot_edge=in_edge;
    std::vector<std::vector<int>> bn_list=in_edge;
    std::vector<int> iTmpSlot(sys_size);
    std::vector<double> bias_rbf=as<std::vector<double>>(bias_rf);
    std::vector<double> bias_int=as<std::vector<double>>(init_1_ratio);
    DNS_Derrida amodel(sys_size, ll_system);
    DNS_Aux_GenerateTopo(net_type, net_f_para, sys_size, &in_edge, &ot_edge, iTmpSlot.data());
    DNS_Aux_GenerateFunc(sys_size, ll_system, iTmpSlot.data(),// <- int* in_deg
        amodel.Labels, bias_rbf.data(), obf_type, obf_ratio, obf_i_para1,
        obf_i_para2, &bn_list, con_ids, con_val);
    amodel.LoadedModel(&in_edge, &ot_edge, &bn_list);
    amodel.DynamicParaConfig(bias_int.data(), init_dis);// Set dynamic parameters.
    amodel.DerridaDamageSpread(sim_step, RuleType);
    return amodel.FinalDistance();
}

// Inner function: Derrida Loaded mode.
//[[Rcpp::export]]
double c_Derrida_Load(Rcpp::List &aRealNet, IntegerVector &PointedGene, IntegerVector &PointValues,
    IntegerVector &InD, IntegerVector &OtD, double init_dis, Rcpp::NumericVector init_1_ratio,
    int sim_step, int RuleType){
    Rcpp::List Returner;
    Rcpp::IntegerVector res(5);//res[3]=0;res[4]=0;// 3: External (No input node), 4: Terminal (No output node)
    std::vector<std::vector<int>> ot_edge, in_edge, bn_list;
    c_inner_Robj2Cobj(aRealNet, PointedGene, PointValues,InD, OtD,// In !!
        ot_edge, in_edge, bn_list, res, true, 2L);//Out!! 
    std::vector<double> bias_int=as<std::vector<double>>(init_1_ratio);
    DNS_Derrida amodel(res[0], 2);
    //Rcpp::List tmplist=Cobj2Robj_Net(res[0], &in_edge, &ot_edge, &bn_list);
    amodel.LoadedModel(&in_edge, &ot_edge, &bn_list);
    amodel.DynamicParaConfig(bias_int.data(), init_dis);
    amodel.DerridaDamageSpread(sim_step, RuleType);
    return amodel.FinalDistance();
}

//' @title Analyze percolation within system
//' @description C++ prototype function of \link{DNS_Percolation}
//' @param sys_size [\code{int}] system size
//' @param ll_system [\code{int}] level number of discrete system
//' @param sim_step [\code{int}] simulating steps
//' @param lat_type [\code{int}] lattice type
//' @param obr_window [\code{int}] observe window
//' @param obf_type [\code{char}] ordered Boolean function type
//' @param bias_rf [\code{NumericVector}] biases in random function 
//' @param obf_ratio [\code{double}] proportion of ordered Boolean function 
//' @param init_bias [\code{NumericVector}] "1" proportion of initial state 
//' @param net_type [\code{char}] system topological type
//' @param net_f_para [\code{char}] topological configured parameters
//' @param obf_i_para1 [\code{int}] configuration parameter for OBF (1)
//' @param obf_i_para2 [\code{int}] configuration parameter for OBF (2)
//' @param OutPut [\code{bool}] whehther output all states?
//' @param RuleType [\code{int}] update strategy
//[[Rcpp::export]]
Rcpp::List c_Percolation_Simualtion(int sys_size, int ll_system, int sim_step, 
    int lat_type, int obr_window, char obf_type, 
    Rcpp::NumericVector bias_rf,
    double obf_ratio,
    Rcpp::NumericVector init_bias, 
    char net_type, // Only simulate pseudo-networks.
    double net_f_para, int obf_i_para1, int obf_i_para2, bool OutPut,
    int RuleType){
    // Return values.
    Rcpp::List perco_sim_case(3);
    Rcpp::NumericVector MaxCluster(2);
    // Set a pseudo- model & configure it.
    //mt.seed(RandSeed);
    std::vector<int> con_ids, con_val; //con_id.clear(); con_val.clear();
    std::vector<std::vector<int>> in_edge(sys_size);
    for(int ii=0; ii<sys_size; ++ii){
        in_edge[ii].resize(0);}
    std::vector<std::vector<int>> ot_edge=in_edge;
    std::vector<std::vector<int>> bn_list=in_edge;
    std::vector<int> iTmpSlot(sys_size);// Record in_dege & unstable type
    std::vector<double> bias_rbf=as<std::vector<double>>(bias_rf);
    std::vector<double> Init_Bias=as<std::vector<double>>(init_bias);
    DNS_Percolation amodel(sys_size, ll_system, lat_type);
    DNS_Aux_GenerateTopo(net_type, net_f_para, sys_size, &in_edge, &ot_edge, iTmpSlot.data());
    DNS_Aux_GenerateFunc(sys_size, ll_system, iTmpSlot.data(),// <- int* in_deg
        amodel.Labels, bias_rbf.data(), obf_type, obf_ratio, obf_i_para1,
        obf_i_para2, &bn_list, con_ids, con_val);
    amodel.LoadedModel(&in_edge, &ot_edge, &bn_list);
    amodel.DynamicParaConfig(Init_Bias.data());// Initial bais of each values.
    amodel.PercolationModel(sim_step, obr_window, RuleType);
    // Calculate it.
    amodel.PercolationWhetherNot().PercolationStableFraction(MaxCluster.begin());
    if(OutPut){// Output the current states.
        Rcpp::IntegerVector StableLattice_rINT(sys_size);
        Rcpp::IntegerVector Stable_in_MaxClust(sys_size);
        int *ipt1=StableLattice_rINT.begin();
        int *ipt2=Stable_in_MaxClust.begin();
        amodel.OutputFinalLattice(ipt1,ipt2);
        perco_sim_case[1]=StableLattice_rINT;
        perco_sim_case[2]=Stable_in_MaxClust;
    }
    perco_sim_case[0]=MaxCluster;
    return perco_sim_case;
}

//' @title Analyze potential engaged nodes within Boolean/multi-valued system
//' @description C++ prototype function of \link{DNS_Engaged}
//' @param sys_size [\code{int}] system size
//' @param ll_system [\code{int}] level number of discrete system
//' @param obf_type [\code{char}] ordered Boolean function type
//' @param bias_rf [\code{NumericVector}] biases in random function 
//' @param obf_ratio [\code{double}] proportion of ordered Boolean function 
//' @param net_type [\code{char}] system topological type
//' @param net_f_para [\code{char}] topological configured parameters
//' @param obf_i_para1 [\code{int}] configuration parameter for OBF (1)
//' @param obf_i_para2 [\code{int}] configuration parameter for OBF (2)
//' @param PointedNode [\code{IntegerVector}] which nodes should be fixed
//' @param PointValues [\code{IntegerVector}] pointed values of free/controlled nodes
//' @param NodeDetailInfor [\code{bool}] should return detail informations of node
//' @param ReturnResidualNetwork [\code{bool}] should return residual networks?
//[[Rcpp::export]]
Rcpp::List c_ScalingLaw_Simualtion(int sys_size, int ll_system, char obf_type, 
    Rcpp::NumericVector bias_rf, double obf_ratio,
    char net_type, double net_f_para, int obf_i_para1, int obf_i_para2, 
    Rcpp::IntegerVector PointedNode, Rcpp::IntegerVector PointValues,
    int NodeDetailInfor, int ReturnResidualNetwork){
    // Return values.
    Rcpp::List Returner;
    Rcpp::IntegerVector res(5);
    int ii, tmp_res[3];
    // Check whether exist non-exponenet controller?
    std::vector<int> con_ids, con_val;
    if(PointedNode[0]>=0){// Exist!! Convert into cpp-format; unlike BioNetAnalysis(), 
        //  ... constant value's setting in the inner function.
        con_ids=as<std::vector<int>>(PointedNode);
        con_val=as<std::vector<int>>(PointValues);}
    std::vector<double> Biasss=as<std::vector<double>>(bias_rf);
    // Set a pseudo-model.
    //mt.seed(RandSeed);
    std::vector<std::vector<int>> in_edge(sys_size);
    for(ii=0; ii<sys_size; ++ii){
        in_edge[ii].resize(0);}
    std::vector<std::vector<int>> ot_edge=in_edge;
    std::vector<std::vector<int>> bn_list=in_edge;
    std::vector<int> UnstableType(sys_size);// Record in_dege & unstable type
    std::vector<int> Num_ID(sys_size);
    for(ii=0; ii<sys_size; ++ii){
        Num_ID[ii]=ii;}
    DNS_Aux_GenerateTopo(net_type, net_f_para, sys_size, &in_edge, &ot_edge, UnstableType.data());
    DNS_Aux_GenerateFunc(sys_size, ll_system, UnstableType.data(),// <- int* in_deg
        Num_ID.data(), Biasss.data(), obf_type, obf_ratio, obf_i_para1,
        obf_i_para2, &bn_list, con_ids, con_val);
    DNS_Engaged amodel(sys_size, ll_system, &in_edge, &ot_edge, &bn_list);
    amodel.OnlyScalingPattern(tmp_res, UnstableType.data());
    // Set return 1: integrate vector (Stable, useless, unstable node).
    res[0]=tmp_res[0]; res[1]=tmp_res[1]; res[2]=tmp_res[2]; // Generally, no external nodes.
    Returner.push_back(res);
    // Set return 2: node's attribute
    if(NodeDetailInfor>0){// Return detail node's attributes.
        Returner.push_back(Rcpp::wrap(amodel.Export2AllNodeState()));// Node: stable (0|1), unstable
        Returner.push_back(Rcpp::wrap(amodel.Export2AllNodeType()));
        Returner.push_back(Rcpp::wrap(UnstableType));
    }// Node: stable, useless, engaged.
    else {
        Returner.push_back(NA_LOGICAL);// or NA_INTEGER
        Returner.push_back(NA_LOGICAL);
        Returner.push_back(NA_LOGICAL);}
    // Set return 3: residual network.
    if((ReturnResidualNetwork>0)&&(res[2]>0)){// Should return and exist engaged nodes.
        // std::vector<std::vector<int>> in_edge, ot_edge, bn_list;
        amodel.Export2ResidualNetwork(in_edge, ot_edge, bn_list);
        Rcpp::List ins, ots, bns, res_net;
        for(ii=0; ii<sys_size; ++ii){
            if(in_edge[ii].size()){// Has inputs.
                ins.push_back(Rcpp::wrap(in_edge[ii]));
                bns.push_back(Rcpp::wrap(bn_list[ii]));}
            else {
                ins.push_back(NA_LOGICAL);//R_NilValue
                bns.push_back(NA_LOGICAL);}
            if(ot_edge[ii].size()){// Has outputs.
                ots.push_back(Rcpp::wrap(ot_edge[ii]));}
            else {
                ots.push_back(NA_LOGICAL);}}
        //res_net.push_back(aRealNet[0],"AllMember");
        Rcpp::IntegerVector codenumber(sys_size);
        res_net.push_back(codenumber,"AllMember");
        res_net.push_back(ins,"InEdge");
        res_net.push_back(ots,"OutEdge");
        res_net.push_back(bns,"BoolFun");
        Returner.push_back(res_net);}
    else {
        Returner.push_back(NA_LOGICAL);}
    return Returner;
}

//' @title Analyze core dynamic components of Boolean networks (2)
//' @description C++ prototype function of \link{BoolBioNet_CoreDyn} (\code{scaling} part). 
//' @param aRealNet [\code{List}] a real/silico genetic network (See \link{BoolGRN_CellCollective}).
//' @param PointedGene [\code{IntegerVector}] which nodes should be fixed.
//' @param PointValues [\code{IntegerVector}] pointed values of free/controlled nodes.
//' @param InD [\code{IntegerVector}] the vector of in-degrees.
//' @param OtD [\code{IntegerVector}] the vector of out-degrees.
//' @param NodeDetailInfor [\code{int}] should return detail information of node (>0, yes).
//' @param ReturnResidualNetwork [\code{int}] recursive number of analysis (>0).
//[[Rcpp::export]]
Rcpp::List c_ScalingLaw_RealNet(Rcpp::List &aRealNet, IntegerVector &PointedGene, IntegerVector &PointValues,
    IntegerVector &InD, IntegerVector &OtD, int NodeDetailInfor, int ReturnResidualNetwork){
    // Return values.
    Rcpp::List Returner;
    Rcpp::IntegerVector res(5);// 3: External (No input node), 4: Terminal (No output node)
    int tmp_res[3];
    std::vector<std::vector<int>> in_edge, ot_edge, bn_list;
    c_inner_Robj2Cobj(aRealNet, PointedGene, PointValues, InD, OtD,// In !!
        ot_edge, in_edge, bn_list, res, true, 2L);//Out!!
    int sys_size=res[0];
    DNS_Engaged amodel(res[0], 2L, &in_edge, &ot_edge, &bn_list);
    std::vector<int> UnstaleType(res[0]);
    amodel.OnlyScalingPattern(tmp_res, UnstaleType.data());
    // Set return 1: integrate vector (Stable, useless, unstable node).
    res[0]=tmp_res[0]; res[1]=tmp_res[1]; res[2]=tmp_res[2];
    Returner.push_back(res);
    // Set return 2: node's attribute
    if(NodeDetailInfor>0){// Return detail node's attributes.
        Returner.push_back(Rcpp::wrap(amodel.Export2AllNodeState()));// Node: stable (0|1), unstable
        Returner.push_back(Rcpp::wrap(amodel.Export2AllNodeType()));// Node: stable, useless, engaged
        //Returner.push_back(Rcpp::wrap(UnstaleType));// Node: Unstable Types
    }
    else {
        Returner.push_back(NA_LOGICAL);// NA_LOGICAL
        Returner.push_back(NA_LOGICAL);}
    // Set return 3: residual network.
    if((ReturnResidualNetwork>0)&&(res[2]>0)){// Should return and exist engaged nodes.
        amodel.Export2ResidualNetwork(in_edge, ot_edge, bn_list);// Reuse slots.
        Rcpp::List ins, ots, bns, res_net;
        for(int ii=0; ii<sys_size; ++ii){
            if(in_edge[ii].size()){// Has inputs.
                ins.push_back(Rcpp::wrap(in_edge[ii]));
                bns.push_back(Rcpp::wrap(bn_list[ii]));}
            else {
                ins.push_back(NA_LOGICAL);
                bns.push_back(NA_LOGICAL);}
            if(ot_edge[ii].size()){// 
                ots.push_back(Rcpp::wrap(ot_edge[ii]));}
            else {
                ots.push_back(NA_LOGICAL);}}
        res_net.push_back(aRealNet[0],"AllMember");
        res_net.push_back(ins,"InEdge");
        res_net.push_back(ots,"OutEdge");
        res_net.push_back(bns,"BoolFun");
        Returner.push_back(res_net);}
    else {
        Returner.push_back(NA_INTEGER);}
    return Returner;
}

//' @title Analyze core dynamic components of Boolean networks (1)
//' @description C++ prototype function of \link{BoolBioNet_CoreDyn} (\code{coredyn} part).
//' @param aRealNet [\code{List}] a real/silico genetic network (See \link{BoolGRN_CellCollective}).
//' @param PointedGene [\code{IntegerVector}] which nodes should be fixed.
//' @param PointValues [\code{IntegerVector}] pointed values of free/controlled nodes.
//' @param InD [\code{IntegerVector}] the vector of in-degrees.
//' @param OtD [\code{IntegerVector}] the vector of out-degrees.
//' @param NodeDetailInfor [\code{int}] should return detail information of node (>0, yes).
//' @param ReturnResidualNetwork [\code{int}] should return Residual Network (>0, yes).
//' @param Times [\code{int}] recursive number of analysis (>0).
// [[Rcpp::export]]
Rcpp::List c_CoreDynamicNode(Rcpp::List aRealNet, IntegerVector PointedGene, IntegerVector PointValues,
    IntegerVector InD, IntegerVector OtD, int NodeDetailInfor, int ReturnResidualNetwork, int Times){
    // Return values.
    Rcpp::List Returner;
    Rcpp::IntegerVector res(5);//res[3]=0;res[4]=0;// 3: External (No input node), 4: Terminal (No output node)
    int tmp_res[3];
    std::vector<std::vector<int>> in_edge, ot_edge, bn_list;
    c_inner_Robj2Cobj(aRealNet, PointedGene, PointValues,InD, OtD,// In !!
        ot_edge, in_edge, bn_list, res, true, 2L);//Out!!
    int sys_size=res[0];
    DNS_CoreDyn amodel(res[0], 2L, &in_edge, &ot_edge, &bn_list);
    amodel.OnceCoreDynamic(tmp_res, Times);
    // Set return 1: integrate vector (Stable, useless, unstable node).
    res[0]=tmp_res[0]; res[1]=tmp_res[1]; res[2]=tmp_res[2];
    Returner.push_back(res);
    // Set return 2: node's attribute
    if(NodeDetailInfor>0){// Return detail node's attributes.
        Returner.push_back(Rcpp::wrap(amodel.Export2AllNodeState()));// Node: stable (0|1), unstable
        Returner.push_back(Rcpp::wrap(amodel.Export2AllNodeType()));}// Node: stable, useless, engaged.
    else {
        Returner.push_back(NA_LOGICAL);// NA_LOGICAL
        Returner.push_back(NA_LOGICAL);}
    // Set return 3: residual network.
    if((ReturnResidualNetwork>0)&&(res[2]>0)){// Should return and exist engaged nodes.
        amodel.Export2ResidualNetwork(in_edge, ot_edge, bn_list);
        Rcpp::List ins, ots, bns, res_net;
        for(int ii=0; ii<sys_size; ++ii){
            if(in_edge[ii].size()){// Has inputs.
                ins.push_back(Rcpp::wrap(in_edge[ii]));
                bns.push_back(Rcpp::wrap(bn_list[ii]));}
            else {
                ins.push_back(NA_LOGICAL);
                bns.push_back(NA_LOGICAL);}
            if(ot_edge[ii].size()){
                ots.push_back(Rcpp::wrap(ot_edge[ii]));}
            else {
                ots.push_back(NA_LOGICAL);}}
        res_net.push_back(aRealNet[0],"AllMember");
        res_net.push_back(ins,"InEdge");
        res_net.push_back(ots,"OutEdge");
        res_net.push_back(bns,"BoolFun");
        Returner.push_back(res_net);}
    else {
        Returner.push_back(NA_INTEGER);}
    return Returner;
}

//' @title Convert a Boolean function into threshold-based or polynomial forms
//' @description C++ prototype function of \link{BoolFun_Polynomial}.
//' @param VariableMat [\code{IntegerMatrix &}] a matrix of variable combination 
//' @param MapTab [\code{IntegerMatrix &}] the truth table of Boolean function.
//' @param LogiSpin [\code{int}] Should show the detail polynomial form?
// [[Rcpp::export]]
Rcpp::List c_BoolFun2Polynomial(Rcpp::IntegerMatrix &VariableMat, Rcpp::IntegerVector &MapTab, int LogiSpin){
    return (PolynomialFunction(VariableMat, MapTab, LogiSpin));
}

//' @title Caluclate the network's strongly connected components
//' @description C++ prototype function of \link{BoolBioNet_StrConComp}
//' @param aRealNet [\code{List}] a real/silico genetic network (See \link{BoolGRN_CellCollective}).
//' @param InD [\code{IntegerVector}] the vector of in-degrees.
//' @param OtD [\code{IntegerVector}] the vector of out-degrees.
// [[Rcpp::export]]
Rcpp::List c_StrongConnectComponent(Rcpp::List &aRealNet, Rcpp::IntegerVector &InD, Rcpp::IntegerVector &OtD){
    // Return values.
    Rcpp::List Returner;
    Rcpp::IntegerVector res(5);
    std::vector<std::vector<int>> in_edge, ot_edge, bn_list;
    c_inner_Robj2Cobj(aRealNet, res, res, InD, OtD,// In !! (2nd, 3rd no use here)
        ot_edge, in_edge, bn_list, res, false, 2L);//Out!!
    NetGraphFrame aNet;// No need {DNS_Basic}, {NetGraphFrame} enough!
    aNet.ConfigurationBuildNet('N',res[0],-1,-1,-1.0).Build_GraphNet();
    aNet.LoadFromVecVecIntFrame(in_edge,ot_edge);
    std::vector<std::vector<int>> scc;
    aNet.Tarjon(scc);
    Rcpp::List res_scc;
    Rcpp::IntegerVector non_trivial(scc.size(),0);
    for(int ii=0; ii<((int)(scc.size())); ++ii){
        res_scc.push_back(scc[ii]);
        if(scc[ii].size()>1){
            non_trivial[ii]=1;}
        else {// Only one node.
            int id=scc[ii][0];
            non_trivial[ii]=(std::find(in_edge[id].begin(), in_edge[id].end(),
                id)!=in_edge[id].end());}}
    Returner.push_back(res_scc);
    Returner.push_back(non_trivial);
    return Returner;
}

//' @title Calculate the complexity of a multi-valued function
//' @description C++ prototype function of \link{MulVFun_Complexity}
//' @param avec [\code{IntegerVector}] the function.
//' @param k [\code{int}] the number of input variables.
//' @param L [\code{int}] the level of function.
//' @return [\code{double}] the complexity of given function.
// [[Rcpp::export]]
double c_MulF_Complexity(IntegerVector &avec, int k, int L){
    int Lens=(int)pow(L,k), First=avec[0];
    int *ttt=avec.begin();
    bool constant=true;
    double xx=0;
    for(int ii=0; ii<Lens; ++ii){
        constant=constant&&(First==avec[ii]);// Deal with 1/0-constant function.
        if(!constant)break;}
    if(!constant){
        xx=toR_BoolMulComplexity(ttt, k, L, Lens, true);}
    return xx;
}

//' @title Calculate the effetive edges loading of a multi-valued function
//' @description C++ prototype function (1) of \link{MulVFun_Effectiveness}
//' @param avec [\code{IntegerVector}] the function.
//' @param k [\code{int}] the number of input variables.
//' @param L [\code{int}] the level of function.
//' @return [\code{double}] the global effetive or redundant edges of given function.
//[[Rcpp::export]]
double c_MulF_Effective(IntegerVector &avec, int k, int L){
    int Lens=(int)pow(L,k), First=avec[0];
    int *ttt=avec.begin();
    bool constant=true;
    double xx=0;
    for(int ii=0; ii<Lens; ++ii){
        constant=constant&&(First==avec[ii]);
        if(!constant)break;}
    if(!constant){
        xx=toR_MultipleEffective(ttt, k, L, Lens);}
    return xx;
}

//' @title Calculate the effetive edges loading of a multi-valued function
//' @description C++ prototype function (2) of \link{MulVFun_Effectiveness}
//' @param avec [\code{IntegerVector}] the function.
//' @param k [\code{int}] the number of input variables.
//' @param L [\code{int}] the level of function.
//' @return [\code{NumericVector}] the effetive or redundant value of each edge.
//[[Rcpp::export]]
Rcpp::NumericVector c_MulF_EffectiveEdges(IntegerVector &avec, int k, int L){
    int Lens=(int)pow(L,k);
    int *ttt=avec.begin();
    bool constant=true;
    std::vector<double> tmp_res;
    for(int ii=0; ii<Lens; ++ii){
        constant=constant&&(avec[0]==avec[ii]);
        if(!constant)break;}
    if(!constant){
        tmp_res=toR_MultipleEdgeConnect(ttt, k, L, Lens);}
    else {
        tmp_res.resize(k,0);}
    Rcpp::NumericVector Res=Rcpp::wrap(tmp_res);
    return Res;
}

//' @title Show the Quine-McCluskey form of a multi-valued function
//' @description C++ prototype function of \link{MulVFun_QMForm}
//' @param avec [\code{IntegerVector &}] that represents a mapping table of multi-valued function.
//' @param k [\code{int}] that denotes the number of input variables.
//' @param L [\code{int}] that represents the level of multi-valued system.
//[[Rcpp::export]]
IntegerMatrix c_MulF_QuineMcCluskey(IntegerVector &avec, int k, int L){
    int *ttt=avec.begin();
    IntegerMatrix res=toR_ShowMulVFunDNF(ttt,k,L);
    return res;
}

//' @title Generate a specific type of multi-valued function
//' @description C++ prototype function of \link{MulVFun_Generator}
//' @param FunType [\code{char}] function type.
//' @param k [\code{int}] the number of input variables.
//' @param L [\code{int}] the level of function.
//' @param CanaDeep [\code{int}] the layer of canalization.
//' @param CanaVar [\code{IntegerVector &}] canalized variable indexes.
//' @param CanaVarNum [\code{IntegerVector &}] the canalizing number of each variable.
//' @param CanaInfo1 [\code{List &}] canalizing information.
//' @param CanaInfo2 [\code{List &}] canalized information.
//' @param bias [\code{NumericVector &}] function bias.
//' @param Cana_Free [\code{bool}] is a parameter-free canalized function?
//[[Rcpp::export]]
Rcpp::IntegerVector c_MulVF_Generator(char FunType,int k,int L, int CanaDeep, 
    Rcpp::IntegerVector &CanaVar, Rcpp::IntegerVector &CanaVarNum,
    Rcpp::List &CanaInfo1, Rcpp::List &CanaInfo2,
    Rcpp::NumericVector &bias, bool Cana_Free){
    unsigned int Lens=(unsigned int)(pow(L,k));
    short *bitmap=(short*)malloc(Lens*sizeof(short));
    char fun_type=FunType;
    double *BiasSet=bias.begin();
    // Set config slot.
    std::vector<int> Argus(4,0); Argus[0]=CanaDeep;
    mulvfun amvf(k, L, Lens, FunType);
    amvf.Configuration(BiasSet, bitmap, Argus);
    if(fun_type=='C'){
        // Rprintf("Here okay 1\n");
        if(Cana_Free){// User not offer detail information
            amvf.Gen_Cana_Free(CanaVar,CanaVarNum);
        }
        else {// User have provide information
            // Set 2 std::vector<std::vector<int>> slots.
            // Rprintf("Here okay 2\n");
            std::vector<std::vector<short>> CanaIN(CanaDeep);
            std::vector<std::vector<short>> CanaOT(CanaDeep);
            std::vector<int> EachVarNum(CanaDeep);
            for(int ii=0; ii<CanaDeep; ++ii){
                Rcpp::IntegerVector c_in=CanaInfo1[ii];
                Rcpp::IntegerVector c_ot=CanaInfo2[ii];
                EachVarNum[ii]=c_in.size();
                for(int jj=0; jj<EachVarNum[ii]; ++jj){
                    // Rcpp::IntegerVector c_in_sub=c_in[ii];
                    // Rcpp::IntegerVector c_ot_sub=CanaInfo2[ii];
                    CanaIN[ii].push_back((short)c_in[jj]);
                    CanaOT[ii].push_back((short)c_ot[jj]);
                }
            }
            amvf.Gen_Cana_Config(CanaVar,EachVarNum,CanaIN,CanaOT);
        }
    }
    else {       
        switch(fun_type){
            case 'D':// Domainted-valued
                amvf.Gen_MulVF_Domi();break;
            case 'T':// Linear threshold
                amvf.Gen_MulVF_Thre();break;
            default:// Random
                amvf.Gen_Rand();break;
        }
    }
    amvf.Reset();
    Rcpp::IntegerVector xx(Lens);
    for(unsigned int ii=0; ii<Lens; ++ii){
        xx[ii]=bitmap[ii];
    }
    free(bitmap);
    return xx;
}

//' @title Calculate the sensitivity of a multi-valued function
//' @description C++ prototype function of \link{MulVFun_Sensitivity}
//' @param amulfun [\code{IntegerVector &}] the function.
//' @param k [\code{int}] the number of input variables.
//' @param L [\code{int}] the level of function.
//' @param Lens [\code{int}] the length of mapping table (\eqn{L^k}).
//[[Rcpp::export]]
double c_MulVF_Sensitivity(Rcpp::IntegerVector &amulfun, int k, int L, int Lens){
    mulvfun amvf(k, L, Lens, 'X');
    short *maptab=(short*)malloc(Lens*sizeof(short));
    for(int ii=0; ii<Lens; ++ii){
        maptab[ii]=amulfun[ii];}
    std::vector<int> NoUse(2);
    double xx=amvf.Configuration(nullptr, maptab,NoUse).Sensitivity();
    free(maptab);
    amvf.Reset();
    return xx;
}

//' @title Analyze the structure of a multi-valued nested canalized function
//' @description C++ prototype function of \link{MulVFun_is_NestedCana}
//' @param amulfun [\code{IntegerVector &}] the function.
//' @param k [\code{int}] the number of input variables.
//' @param L [\code{int}] the level of function.
//' @param Lens [\code{int}] the length of mapping table (\eqn{L^k}).
//[[Rcpp::export]]
Rcpp::List c_M_NestedCanalized(Rcpp::IntegerVector &amulfun, int k, int L, int Lens){
    mulvfun amvf(k, L, Lens, 'X');
    short *maptab=(short*)malloc(Lens*sizeof(short));
    for(int ii=0; ii<Lens; ++ii){
        maptab[ii]=amulfun[ii];}
    std::vector<int> Canalizing_V;
    std::vector<int> NoUse(2,-1);
    std::vector<std::vector<std::vector<short>>> tmplist =
        amvf.Configuration(nullptr, maptab, NoUse).NestCana(Canalizing_V, false);
    amvf.Reset();
    free(maptab);
    Rcpp::List FinRes;
    FinRes.push_back(-1);
    if(Canalizing_V.size()>0){
        FinRes[0]=1;
        FinRes.push_back(Canalizing_V);
        Rcpp::List C_in, C_ot;
        for(int ii=0; ii<((int)Canalizing_V.size()); ++ii){
            Rcpp::IntegerVector id1(tmplist[ii][0].begin(), tmplist[ii][0].end());
            C_in.push_back(id1);
            Rcpp::IntegerVector id2(tmplist[ii][1].begin(), tmplist[ii][1].end());
            C_ot.push_back(id2);}
        FinRes.push_back(C_in);
        FinRes.push_back(C_ot);}
    return (FinRes);
}

//' @title Decode the information of a multi-valued linear threshold function
//' @description C++ prototype function of \link{MulVFun_is_Threshold}
//' @param amulfun [\code{IntegerVector &}] the function.
//' @param k [\code{int}] the number of input variables.
//' @param L [\code{int}] the level of function.
//' @param Lens [\code{int}] the length of mapping table (\eqn{L^k}).
//[[Rcpp::export]]
Rcpp::List c_M_Threshold(Rcpp::IntegerVector &amulfun, int k, int L, int Lens){
    mulvfun amvf(k, L, Lens, 'T');
    short *maptab=(short*)malloc(Lens*sizeof(short));
    for(int ii=0; ii<Lens; ++ii){
        maptab[ii]=amulfun[ii];}
    std::vector<int> NoUse(2,-1);
    Rcpp::List FinRes=amvf.Configuration(nullptr, maptab, NoUse).is_MulVF_Thre();
    amvf.Reset();
    free(maptab);
    return (FinRes);
}

//' @title Decode the information of a multi-valued domaineted function
//' @description C++ prototype function of \link{MulVFun_is_Domainted}
//' @param amulfun [\code{IntegerVector &}] the function.
//' @param k [\code{int}] the number of input variables of original function.
//' @param L [\code{int}] the level of original function.
//' @param Lens [\code{int}] the length of mapping table (\eqn{L^k}).
//[[Rcpp::export]]
Rcpp::List c_M_Domainted(Rcpp::IntegerVector &amulfun, int k, int L, int Lens){
    mulvfun amvf(k, L, Lens, 'D');
    short *maptab=(short*)malloc(Lens*sizeof(short));
    for(int ii=0; ii<Lens; ++ii){
        maptab[ii]=amulfun[ii];}
    std::vector<int> NoUse(2,-1);// Useless, just to meet function parameters
    Rcpp::List FinRes=amvf.Configuration(nullptr, maptab, NoUse).is_MulVF_Domi();
    amvf.Reset();
    free(maptab);
    return (FinRes);
}

//' @title Decode the information of a multi-valued signed function
//' @description C++ prototype function of \link{MulVFun_is_Signed}
//' @param amulfun [\code{IntegerVector &}] the function.
//' @param k [\code{int}] the number of input variables of original function.
//' @param L [\code{int}] the level of original function.
//' @param Lens [\code{int}] the length of mapping table (\eqn{L^k}).
//[[Rcpp::export]]
Rcpp::List c_M_Signed(Rcpp::IntegerVector &amulfun, int k, int L, int Lens){
    mulvfun amvf(k, L, Lens, 'S');
    short *maptab=(short*)malloc(Lens*sizeof(short));
    for(int ii=0; ii<Lens; ++ii){
        maptab[ii]=amulfun[ii];}
    std::vector<int> NoUse(2,-1);
    Rcpp::List FinRes=amvf.Configuration(nullptr, maptab, NoUse).is_MulVF_Sign();
    amvf.Reset();
    free(maptab);
    return (FinRes);
}

//' @title Convert a multi-valued function as polynomial like forms
//' @description C++ prototype function of \link{MulVFun_Polynomial}
//' @param VariableMat [\code{IntegerMatrix &}] a matrix that record input vectors.
//' @param MapTab [\code{IntegerVector &}] the function.
//' @param k [\code{int}] the number of input variables of original function.
//' @param L [\code{int}] the level of original function.
// [[Rcpp::export]]
Rcpp::List c_MulVFun2Polynomial(Rcpp::IntegerMatrix &VariableMat, Rcpp::IntegerVector &MapTab, int k, int L){
    return (MulVFun_PolynomialFunction(VariableMat, MapTab, k, L));
}

//' @title Boolean to multi-valued and multi-valued to Boolean transformation
//' @description C++ prototype function of \link{MulV2Bool_Bool2MulV}
//' @param OriMapTab, [\code{IntegerVector}] original mapping table.
//' @param k, [\code{int}] the number of input variables.
//' @param L, [\code{int}] the level of multi-valued system.
//' @param Thresholds, [\code{IntegerVector}] thresholds for a multi-valued system. 
//' @param b2m, [\code{int}] Is from Boolean to Multi-valued system? (>0, yes)
// [[Rcpp::export]]
Rcpp::IntegerVector c_MulV2Bool_Bool2MulV(Rcpp::IntegerVector OriMapTab, int k,
    int L, Rcpp::IntegerVector Thresholds, int b2m){
    std::vector<int> in_fun=as<std::vector<int>>(OriMapTab);
    std::vector<int> thres=as<std::vector<int>>(Thresholds);
    if(b2m){// Bool 2 MulV
        return Rcpp::wrap(Bool2MulV(in_fun, k, L, thres));}
    else {
        return Rcpp::wrap(MulV2Bool(in_fun, k, L, thres));}
}

// [[Rcpp::export]]
Rcpp::List cin_OnlyTempState(int sys_size, int ll_system, int sim_step, char obf_type, 
    Rcpp::NumericVector bias_rf, double obf_ratio,
    double init_dis, Rcpp::NumericVector init_1_ratio, char net_type, // Only simulate pseudo-networks.
    double net_f_para, int obf_i_para1, int obf_i_para2, int RuleType){
    std::vector<int> con_ids, con_val; //con_id.clear(); con_val.clear();
    std::vector<std::vector<int>> in_edge(sys_size);
    for(int ii=0; ii<sys_size; ++ii){
        in_edge[ii].resize(0);}
    std::vector<std::vector<int>> ot_edge=in_edge;
    std::vector<std::vector<int>> bn_list=in_edge;
    std::vector<int> iTmpSlot(sys_size);
    std::vector<double> bias_rbf=as<std::vector<double>>(bias_rf);
    std::vector<double> bias_int=as<std::vector<double>>(init_1_ratio);
    DNS_JustSim amodel(sys_size, ll_system, sim_step);
    DNS_Aux_GenerateTopo(net_type, net_f_para, sys_size, &in_edge, &ot_edge, iTmpSlot.data());
    DNS_Aux_GenerateFunc(sys_size, ll_system, iTmpSlot.data(),// <- int* in_deg
        amodel.Labels, bias_rbf.data(), obf_type, obf_ratio, obf_i_para1,
        obf_i_para2, &bn_list, con_ids, con_val);
    // Return values.
    Rcpp::List Returner;
    Rcpp::List tmplist=Cobj2Robj_Net(sys_size, &in_edge, &ot_edge, &bn_list);
    amodel.LoadedModel(&in_edge, &ot_edge, &bn_list);
    amodel.DynamicParaConfig(bias_int.data());
    amodel.SelfEvolution(sim_step, RuleType);
    Returner.push_back( amodel.Out2Robj(sim_step) );
    Returner.push_back( tmplist );
    return Returner;
}

// [[Rcpp::export]]
Rcpp::List cin_OnlyTempState_Load(Rcpp::List &aRealNet, IntegerVector &PointedGene, IntegerVector &PointValues,
    IntegerVector &InD, IntegerVector &OtD, Rcpp::NumericVector init_1_ratio,
    int sim_step,int RuleType){
    // Return values.
    Rcpp::List Returner;
    Rcpp::IntegerVector res(5);//res[3]=0;res[4]=0;// 3: External (No input node), 4: Terminal (No output node)
    std::vector<std::vector<int>> ot_edge, in_edge, bn_list;
    c_inner_Robj2Cobj(aRealNet, PointedGene, PointValues,InD, OtD,// In !!
        ot_edge, in_edge, bn_list, res, true, 2L);//Out!! 
    std::vector<double> bias_int=as<std::vector<double>>(init_1_ratio);
    DNS_JustSim amodel(res[0], 2, sim_step);
    // Return values.
    Rcpp::List tmplist=Cobj2Robj_Net(res[0], &in_edge, &ot_edge, &bn_list);
    amodel.LoadedModel(&in_edge, &ot_edge, &bn_list);
    amodel.DynamicParaConfig(bias_int.data());
    amodel.SelfEvolution(sim_step, RuleType);
    Returner.push_back( amodel.Out2Robj(sim_step) );
    Returner.push_back( tmplist );
    return Returner;
}

// Inner function: analyze edge's activity.
double c_inner_Activity(int *MapTab, int WhichID, int k){
    int sumer=0;
    int Bits=k-1-WhichID;
    int lower=(1<<Bits), upper=1<<(k-1-Bits);
    int *ipt0=MapTab, *ipt1=MapTab+lower;
    for(int jj=0; jj<upper; ++jj){
        for(int kk=0; kk<lower; ++kk){
            sumer+=(ipt0[kk])!=(ipt1[kk]);}
        ipt0+=(lower<<1);
        ipt1+=(lower<<1);}
    return (double)(sumer<<1)/(double)(1<<k);
}
// Inner function: analyze edge's effectiveness (redundancy).
double c_inner_EffecRedu(int *MapTab, int WhichID, int k){
    double xx=0;
    bool constant=true;
    int First=MapTab[0];
    for(int ii=0; ii<(1<<k); ++ii){
        constant=constant&&(First==MapTab[ii]);
        if(!constant)break;}// Constant
    if(!constant){
        std::vector<double> Results(k,0);
        int *InverseBit=(int*)malloc((1<<k)*sizeof(int));
        for(int ii=0; ii<(1<<k); ++ii){
            InverseBit[ii]=!MapTab[ii];}// 1-pattern + 0-pattern
        OnceEdgeConnect(MapTab, k, 2, Results);// From QuineMcCluskey.h
        OnceEdgeConnect(InverseBit, k, 2, Results);
        free(InverseBit);
        xx=Results[k-1-WhichID];
        xx=xx/(1<<k);}
    return xx;// <-- It is effective not redundancy (1-xx)
}
// Inner function: analyze edge's sign: Activation(0), Inhibition(1), Hybird(2).
int c_inner_EdgeSigned(int *MapTab, int WhichID, int k){
    bool logi0=true, logi1=true, exit=false;
    int Type;
    int Bits=k-1-WhichID;
    int lower=1<<Bits;
    int upper=1<<(k-1-Bits);
    int *ipt0=MapTab, *ipt1=MapTab+lower;
    for(int jj=0; jj<upper; ++jj){
        for(int kk=0; kk<lower; ++kk){
            if(ipt0[kk]>ipt1[kk])logi1=false;// Contradiction to increase
            if(ipt0[kk]<ipt1[kk])logi0=false;// Contradiction to decrease
            if(!(logi1||logi0)){
                exit=true;// Exit repeat judgment, not increase/decrease monotone.
                break;}}
        if(exit)break;
        ipt0+=(lower<<1);
        ipt1+=(lower<<1);}
    if(exit){
        Type=2;}
    else {
        if(logi1)Type=0;// Positive!
        else if(logi0)Type=1;// Negative!
        else Type=-1;}//HALT!!
    return Type;
} 
// Inner function: analyze edge's loading of coupled loops.
Rcpp::List c_inner_EdgeLoadLoop(std::vector<std::vector<int>> &Loops, bool SelfKeep,
    Rcpp::List &BoolFun_List,
    std::vector<std::vector<int>> &InTopo,
    Rcpp::IntegerVector &InD){
    auto edge_key=[](int u, int v) -> int64_t {// Edge's id
        return (static_cast<int64_t>(u) << 32) | (static_cast<int64_t>(v) & 0xFFFFFFFF);};
    auto get_u=[](int64_t key) -> int {// Edge's source 
        return static_cast<int>((key>>32) & 0xFFFFFFFF);};
    auto get_v=[](int64_t key) -> int {// Edge's target 
        return static_cast<int>(key & 0xFFFFFFFF);};
    std::map<int64_t,int> edge_support;
    for(const auto& cycle:Loops){
        int nn=cycle.size()-1;
        if(nn>1||SelfKeep){// more than 2 nodes or allow self-node.
            for(int ii=0; ii<nn; ++ii){
                int64_t key=edge_key(cycle[ii], cycle[ii+1]);
                edge_support[key]++;}}}
    Rcpp::IntegerMatrix result1(edge_support.size(), 4);// Return statistic matrix 
    Rcpp::NumericMatrix result2(edge_support.size(), 2);// Return statistic matrix 
    int row=0;
    for(const auto& [key,sup]:edge_support){
        int ins_=get_u(key), out_=get_v(key);
        int id=std::distance(InTopo[out_].begin(),
            std::find(InTopo[out_].begin(), InTopo[out_].end(), ins_));
        result1(row,0)=ins_;
        result1(row,1)=out_;
        result1(row,2)=sup;
        Rcpp::IntegerVector bf=BoolFun_List[out_];
        result1(row,3)=c_inner_EdgeSigned(bf.begin(), id, InD[out_]);
        result2(row,0)=c_inner_Activity(bf.begin(), id, InD[out_]);
        result2(row,1)=c_inner_EffecRedu(bf.begin(), id, InD[out_]);
        ++row;}
    Rcpp::colnames(result1)=Rcpp::CharacterVector::create("from","to","support","SignType");
    Rcpp::colnames(result2)=Rcpp::CharacterVector::create("activity","effectiveness");
    Rcpp::List Returner;
    Returner.push_back(result1);
    Returner.push_back(result2);
    return Returner;
}

// Inner function: analyze loop's sign: Activation(0), Inhibition(1), Hybird(2).
Rcpp::IntegerMatrix c_inner_LoopSign(std::vector<std::vector<int>> &Loops, 
    Rcpp::List &BoolFun_List,
    std::vector<std::vector<int>> &InTopo,
    Rcpp::IntegerVector &InD){
    Rcpp::IntegerMatrix result(Loops.size(),3);// Return statistic matrix
    size_t row=0;
    for(const auto& aCycle:Loops){
        for(size_t ii=0; ii<(aCycle.size()-1); ++ii){
            int in_=aCycle[ii];
            int ot_=aCycle[ii+1];
            int id=std::distance(InTopo[ot_].begin(), 
                std::find(InTopo[ot_].begin(), InTopo[ot_].end(), in_));
            Rcpp::IntegerVector bf=BoolFun_List[ot_];
            int types=c_inner_EdgeSigned(bf.begin(), id, InD[ot_]);
            result(row,types)=result(row,types)+1;}
        row++;}
    Rcpp::colnames(result)=Rcpp::CharacterVector::create("Positive","Negative","Hybrid");
    return result;
}

// Inner fucntion: Multipurpose functions involving feedback loops.
// [[Rcpp::export]]
Rcpp::List c_ObtainFeedBackLoop(Rcpp::List &aRealNet, Rcpp::IntegerVector &InD, Rcpp::IntegerVector &OtD,
    int Isolated, int Type){// Return network's feedback loops or other features.
    Rcpp::List Returner;
    Rcpp::IntegerVector res(5);
    std::vector<std::vector<int>> in_edge, ot_edge, bn_list;
    c_inner_Robj2Cobj(aRealNet, res, res, InD, OtD,// In!! (2nd, 3rd no use here)
        ot_edge, in_edge, bn_list, res, false,2L);// Out!!
    NetGraphFrame aNet;// No need using class of {DNS_Basic}, {NetGraphFrame} enough!
    aNet.ConfigurationBuildNet('N',res[0],-1,-1,-1.0).Build_GraphNet();
    aNet.LoadFromVecVecIntFrame(in_edge,ot_edge);
    std::vector<std::vector<int>> scc, fbl;
    std::vector<std::vector<int>> scc_valid;
    aNet.Tarjon(scc);
    scc_valid.reserve(scc.size());
    for(int ii=0; ii<((int)(scc.size())); ++ii){
        if(scc[ii].size()>1){
            scc_valid.push_back(scc[ii]);}
        else {// Only one node.
            int id=scc[ii][0];
            if(std::find(in_edge[id].begin(), in_edge[id].end(),id)!=in_edge[id].end()){
                scc_valid.push_back(scc[ii]);}}}
    if(1==Type){// Only list of SCC's loops.
        for(int ii=0; ii<((int)(scc_valid.size())); ++ii){
            fbl.clear();
            aNet.FeedBackLoop(scc_valid[ii],fbl);
            if(scc_valid[ii].size()>1||Isolated){// >1 or allow isolated node
                Rcpp::List result(fbl.size());// vector<vector<int>> to List<IntegerVector>
                for(size_t jj=0; jj<fbl.size(); ++jj){
                    if(fbl[jj].size()>1||Isolated){// >1 or allow isolated node (A <-> A -> B -> A)
                        result[jj]=Rcpp::IntegerVector(fbl[jj].begin(), fbl[jj].end());
                    } else {
                        result[jj]=Rcpp::NumericVector::get_na();}}
                Returner.push_back(result);}}}
    else if(2==Type){// Loop's length distribution.
        std::map<int,int> LoopLengDistr;
        for(int ii=0; ii<((int)(scc_valid.size())); ++ii){
            fbl.clear();
            aNet.FeedBackLoop(scc_valid[ii],fbl);
            for(const auto& _fbl:fbl){
                LoopLengDistr[_fbl.size()-1]++;}}
        int nn=(int)LoopLengDistr.size();
        Rcpp::IntegerMatrix result(nn, 2);
        int rows=0;
        for(const auto& pair:LoopLengDistr){
            result(rows,0)=pair.first;
            result(rows,1)=pair.second;
            ++rows;}
        Returner.push_back(result);
    }
    else if(3==Type||4==Type){
        std::vector<std::vector<int>> tmp_loops;
        Rcpp::List bf_list=aRealNet[3];
        for(int ii=0; ii<((int)(scc_valid.size())); ++ii){
            fbl.clear();
            if(scc_valid[ii].size()>1||Isolated){// >1 or allow isolated node
                aNet.FeedBackLoop(scc_valid[ii],fbl);
                tmp_loops.insert(tmp_loops.end(),fbl.begin(),fbl.end());}}
        if(3==Type){// Loop's length distribution.
            Returner.push_back(c_inner_EdgeLoadLoop(tmp_loops,Isolated, bf_list, in_edge, InD));} 
        else {// Feedback types (signed symbols):
            Returner.push_back(c_inner_LoopSign(tmp_loops, bf_list, in_edge, InD));}
    }
    return Returner;// Is a 3-order list
}

// Inner function: Coarse graining analysis of attractors
// [[Rcpp::export]]
Rcpp::IntegerMatrix c_AttractorSim(Rcpp::List &aRealNet, int SimNum,
    IntegerVector &PointedGene, IntegerVector &PointValues, 
    IntegerVector &InD, IntegerVector &OtD, 
    int sim_step, int RuleType, int L_sys){
    // Return values.
    Rcpp::IntegerVector res(5);//res[3]=0;res[4]=0;// 3: External (No input node), 4: Terminal (No output node)
    std::vector<std::vector<int>> ot_edge, in_edge, bn_list;
    c_inner_Robj2Cobj(aRealNet, PointedGene, PointValues,InD, OtD,// In !!
        ot_edge, in_edge, bn_list, res, true, L_sys);//Out!! 
    int SysSize=res[0];
    DNS_OnceSim amodel(SysSize, L_sys);
    int *tmp_int=(int *)malloc(SimNum*SysSize*sizeof(int));
    Rcpp::IntegerMatrix Returner(SysSize,SimNum);
    std::vector<double> bias_int(L_sys, 1.0/L_sys);
    int *iptr=tmp_int;
    amodel.LoadedModel(&in_edge, &ot_edge, &bn_list);
    for(int ii=0; ii<SimNum; ++ii){
        amodel.DynamicParaConfig(bias_int.data());
        amodel.SelfEvolution(sim_step, RuleType);
        amodel.FinalState(iptr);
        iptr=iptr+SysSize;}
    std::copy(tmp_int, tmp_int + SimNum*SysSize, Returner.begin());
    free(tmp_int);
    iptr=tmp_int=nullptr;
    return Returner;
}

// Inner function: Multi-valued dynamic simulation
// [[Rcpp::export]]
Rcpp::IntegerMatrix c_AttractorSim_MulVal(Rcpp::List &aRealNet_MulVal, int SimNum,
                                          IntegerVector &PointedGene, IntegerVector &PointValues, 
                                          IntegerVector &InD, IntegerVector &OtD, IntegerVector &AllNegative,
                                          int sim_step, int L_sys){
    Rcpp::IntegerVector res(5);
    std::vector<std::vector<int>> ot_edge, in_edge, bn_list;
    c_inner_Robj2Cobj(aRealNet_MulVal, PointedGene, PointValues,InD, OtD,// In !!
                    ot_edge, in_edge, bn_list, res, true, L_sys);//Out!!
    int SysSize=res[0];
    int *tmp_int=(int *)malloc(SimNum*SysSize*sizeof(int));
    Rcpp::IntegerMatrix Returner(SysSize, SimNum);
    int *ipt1, *ipt2, *iptc;
    int *tmp_int2=(int *)malloc(SysSize*sizeof(int));
    std::mt19937 InnerRNG( (int)(INT32_MAX*unif_rand()) );// Generate a seed from R.
    std::uniform_int_distribution<int> State(0,L_sys);
    std::uniform_int_distribution<int> Nodes(0,SysSize-1);
    //Rprintf("Here-2\n");
    for(auto ii=0; ii<SimNum; ++ii){
        ipt1=tmp_int+ii*SysSize;
        ipt2=ipt1;
        for(int jj=0; jj<SysSize; ++jj){// Random initial.
            ipt1[jj]=State(InnerRNG);}
        for(int jj=0; jj<sim_step; ++jj){
            int kk=Nodes(InnerRNG);
            if(InD[kk]>0){
                int tmpx=ipt1[kk];
                double Targets=0;
                for(int pp=0; pp<InD[kk]; ++pp){
                    Targets+=ipt1[in_edge[kk][pp]]*bn_list[kk][pp];}
                if(AllNegative[kk]){// All is negative
                    Targets=Targets+L_sys;}
                else {
                    if(Targets<0){// Max(0, Act-Inh)
                        Targets=0;}}
                if(tmpx>Targets){
                    ipt2[kk]=tmpx-1;
                    if(ipt2[kk]<0)ipt2[kk]=0;}
                else if(tmpx<Targets){
                    ipt2[kk]=tmpx+1;
                    if(ipt2[kk]>L_sys)ipt2[kk]=L_sys;}
                else {
                    ipt2[kk]=tmpx;}}   
            else {
                ipt2[kk]=ipt1[kk];}}}
    std::copy(tmp_int, tmp_int + SimNum*SysSize, Returner.begin());
    free(tmp_int);
    free(tmp_int2);
    return Returner;
}
// Code is over!