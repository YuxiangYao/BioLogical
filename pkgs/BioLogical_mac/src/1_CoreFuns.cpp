#include <Rcpp.h>
#include "BoolFun.h"
#include "QuineMcCluskey.h"
#include "NetGraphFrame.h"
#include "DNS_BasicFrame.h"
#include "MulVFun.h"

using namespace Rcpp;

//' C++ prototype function for BoolFun_Type()
//' @param boolfunss, bool function logical vector.
//' @param k, input variable size.
//' @param Leixing, function's type.
//' @param Showit, whether show the function?
//' @details Please see document of \link{BoolFun_Type}
//[[Rcpp::export]]
bool c_BF_isPointed(LogicalVector boolfunss,int k,char Leixing,bool Showit){
    boolfun abf;
    bool *bitmap=(bool*)malloc((1<<k)*sizeof(bool));
    for(int ii=0;ii<(1<<k);ii++){
        bitmap[ii]=boolfunss[ii];}
    bool xx=false;
    xx=abf.Configuration('X',k,0.5,bitmap,-1,-1,-1,-1).is_PointedType(Leixing,Showit);
    free(bitmap);
    return xx;
}

//' C++ prototype function for BoolFun_Generator()
//' @param Leixing, function's type.
//' @param k, input variable size.
//' @param bias, function's bias.
//' @param Vars, some configuration parameters.
//' @details Please see document of \link{BoolFun_Generator}
//[[Rcpp::export]]
Rcpp::IntegerVector c_BF_Generator(char Leixing, int k, double bias, IntegerVector Vars){
    boolfun abf;
    bool *bitmap=(bool*)malloc((1<<k)*sizeof(bool));
    abf.Configuration(Leixing,k,bias,bitmap,Vars[0],Vars[1],Vars[2],Vars[3]).Gen_BF();
    Rcpp::IntegerVector xx(1<<k);
    for(int ii=0; ii<(1<<k); ii++){
        xx[ii]=bitmap[ii]?1:0;}
    free(bitmap);
    return xx;
}

//' C++ prototype function for BoolFun_Sensitivity()
//' @param boolfunss, bool function logical vector.
//' @param k, input variable size.
//' @details Please see document of \link{BoolFun_Sensitivity}
//[[Rcpp::export]]
double c_BF_Sensitivity(LogicalVector &boolfunss, int k){
    boolfun abf;
    bool *bitmap=(bool*)malloc((1<<k)*sizeof(bool));
    for(int ii=0; ii<(1<<k); ++ii){
        bitmap[ii]=boolfunss[ii];}
    double xx=abf.Configuration('X',k,0.5,bitmap,-1,-1,-1,-1).Sensitivity();
    free(bitmap);
    return xx;
}

//' C++ prototype function for BoolFun_EffectiveEdges()
//' @param boolfunss, bool function logical vector.
//' @param k, input variable size.
//' @details Please see document of \link{BoolFun_EffectiveEdges}
//[[Rcpp::export]]
double c_BF_Effective(IntegerVector &boolfunss,int k){
    bool constant=true;
    double xx=0;
    int *ttt=boolfunss.begin();
    int First=boolfunss[0];
    for(int ii=0; ii<(1<<k); ++ii){
        constant=constant&&(First==boolfunss[ii]);
        if(!constant)break;}
    if(!constant){
        xx=toR_MultipleEffective(ttt, k, 2, (1<<k));}
    return xx;
}

//' C++ prototype function for BoolFun_Complexity()
//' @param boolfunss, bool function logical vector.
//' @param k, input variable size.
//' @details Please see document of \link{BoolFun_Complexity}
//[[Rcpp::export]]
double c_BF_Complexity(IntegerVector &boolfunss, int k){
    int *ttt=boolfunss.begin();
    // int *ttt=(int*)malloc((1<<k)*sizeof(int));
    // for(int ii=0; ii<(1<<k); ++ii){
    //     ttt[ii]=boolfunss[ii];}
    double xx=toR_BoolMulComplexity(ttt, k, 2, (1<<k), true);// Fixed L:2 , lens:1<<k
    // free(ttt);ttt=nullptr;
    return xx;
}

//' C++ prototype function for BoolFun_EffectiveEdges()
//' @param boolfunss, truth table.
//' @param k, input variable size.
//' @details Please see document of \link{BoolFun_EffectiveEdges}
//[[Rcpp::export]]
Rcpp::NumericVector c_BF_EffectiveEdges(IntegerVector &boolfunss, int k){
    // int *ttt=(int*)malloc((1<<k)*sizeof(int));
    int *ttt=boolfunss.begin();
    int First=boolfunss[0];
    bool constant=true;
    std::vector<double> tmp_res;
    for(int ii=0; ii<(1<<k); ++ii){
        constant=constant&&(First==boolfunss[ii]);
        if(!constant)break;}
    if(!constant){
        tmp_res=toR_MultipleEdgeConnect(ttt, k, 2, (1<<k));}
    else {
        tmp_res.resize(k,0);}
    Rcpp::NumericVector Res=Rcpp::wrap(tmp_res);
    return Res;
}

//' C++ prototype function for BoolFun_QMForm()
//' @param boolfunss, bool function logical vector.
//' @param k, input variable size.
//' @param VarsName, variables' names.
//' @details Please see document of \link{BoolFun_QMForm}
//[[Rcpp::export]]
void c_BF_QuineMcCluskey(IntegerVector &boolfunss, int k, CharacterVector VarsName){
    int *ttt=boolfunss.begin();
    toR_ShowBoolFunDNF(ttt,k,VarsName);
}

//' C++ prototype function for FrameTruthTable()
//' @param VarNum, input variable size.
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

//' C++ prototype function for BoolFun_NestedCanalized()
//' @param boolfunss, truth table.
//' @param k, input variable size.
//' @details Please see document of \link{BoolFun_NestedCanalized}
//[[Rcpp::export]]
Rcpp::List c_B_NestedCanalized(IntegerVector boolfunss,int k){
    boolfun abf;
    bool *bitmap=(bool*)malloc((1<<k)*sizeof(bool));
    for(int ii=0;ii<(1<<k);ii++){
        bitmap[ii]=boolfunss[ii];}
    std::vector<std::vector<int>> xx;
    xx=abf.Configuration('X',k,0.5,bitmap,-1,-1,-1,-1).NestCana();
    free(bitmap);
    Rcpp::List Res;
    Rcpp::IntegerVector id1(xx[0].begin(), xx[0].end());Res.push_back(id1);
    Rcpp::IntegerVector id2(xx[1].begin(), xx[1].end());Res.push_back(id2);
    Rcpp::IntegerVector id3(xx[2].begin(), xx[2].end());Res.push_back(id3);
    return (Res);
}

//' C++ prototype function for DNS_DamageSpread()
//' @param sys_size, system size
//' @param ll_system, level number of discrete system
//' @param sim_step, simulating steps
//' @param obf_type, ordered Boolean function type
//' @param bias_rf, biases in random function 
//' @param obf_ratio, proportion of ordered Boolean function 
//' @param init_dis, initial normalized Hamming distance
//' @param init_1_ratio, "1" proportion of initial vector
//' @param net_type, system topological type
//' @param net_f_para, topological configured parameters
//' @param obf_i_para1, configuration parameter for OBF (1)
//' @param obf_i_para2, configuration parameter for OBF (2)
//' @param RuleType, update strategy 
//' @details Please see document of \link{DNS_DamageSpread}
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

//' C++ prototype function for DNS_Percolation()
//' @param sys_size, system size
//' @param ll_system, level number of discrete system
//' @param sim_step, simulating steps
//' @param lat_type, lattice type
//' @param obr_window, Observe window
//' @param obf_type, ordered Boolean function type
//' @param bias_rf, biases in random function 
//' @param obf_ratio, proportion of ordered Boolean function 
//' @param init_bias, each 
//' @param net_type, system topological type
//' @param net_f_para, topological configured parameters
//' @param obf_i_para1, configuration parameter for OBF (1)
//' @param obf_i_para2, configuration parameter for OBF (2)
//' @param OutPut, whehther output all states?
//' @param RuleType, update strategy
//' @details Please see document of \link{DNS_Percolation}
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

//' C++ prototype function for DNS_Engaged()
//' @param sys_size, system size
//' @param ll_system, level number of discrete system
//' @param obf_type, ordered Boolean function type
//' @param bias_rf, biases in random function 
//' @param obf_ratio, proportion of ordered Boolean function 
//' @param net_type, system topological type
//' @param net_f_para, topological configured parameters
//' @param obf_i_para1, configuration parameter for OBF (1)
//' @param obf_i_para2, configuration parameter for OBF (2)
//' @param PointedNode, set which nodes should be controlled.
//' @param PointValues, set corresponding values of controlled genes.
//' @param NodeDetailInfor, should return detail informations of node.
//' @param ReturnResidualNetwork, should return residual networks?
//' @details Please see document of \link{DNS_Engaged}
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

//' C++ prototype function for BoolBioNet_CoreDyn()
//' @param aRealNet, list-object of real/silico genetic network
//' @param PointedGene, integer or name vector, which nodes should be fixed.
//' @param PointValues, pointed values, for free or controlled nodes.
//' @param InD, in-degree vector.
//' @param OtD, out-degree vector.
//' @param NodeDetailInfor, should return detail informations of node.
//' @param ReturnResidualNetwork, should return Residual Network >0: yes.
//' @details Please see document of \link{BoolBioNet_CoreDyn}
//[[Rcpp::export]]
Rcpp::List c_ScalingLaw_RealNet(Rcpp::List aRealNet, IntegerVector PointedGene, IntegerVector PointValues,
    IntegerVector InD, IntegerVector OtD, int NodeDetailInfor, int ReturnResidualNetwork){
    // Return values.
    Rcpp::List Returner;
    Rcpp::IntegerVector res(5);res[3]=0;res[4]=0;// 3: External (No input node), 4: Terminal (No output node)
    int tmp_res[3];
    // Convert List into appropriate data_frame.
    Rcpp::IntegerVector tmpdata;
    Rcpp::List tmplist1=aRealNet[1], tmplist2=aRealNet[2], tmplist3=aRealNet[3];
    int ii, jj, sys_size=tmplist1.size();
    std::vector<std::vector<int>> in_edge(sys_size);
    for(ii=0; ii<sys_size; ++ii){
        in_edge[ii].resize(0);}
    std::vector<std::vector<int>> ot_edge=in_edge;
    std::vector<std::vector<int>> bn_list=in_edge;
    for(ii=0; ii<sys_size; ++ii){
        if(OtD[ii]>0){// Outs not empty // if(Rcpp::internal::Rcpp_IsNA(tmpdata[0]))
            tmpdata=tmplist2[ii];
            for(jj=0; jj<tmpdata.size(); ++jj){
                ot_edge[ii].push_back(tmpdata[jj]);}}
        else {
            res[4]=res[4]+1;}
        if(InD[ii]>0){// Ins not empty // if(Rcpp::internal::Rcpp_IsNA(tmpdata[0]))
            tmpdata=tmplist1[ii];
            for(jj=0; jj<tmpdata.size(); ++jj){
                in_edge[ii].push_back(tmpdata[jj]);}
            tmpdata=tmplist3[ii];
            for(jj=0; jj<tmpdata.size(); ++jj){
                bn_list[ii].push_back(tmpdata[jj]);}}
        else {// No input, should be set one
            bn_list[ii].push_back(PointValues[ii]);// Actually this value could be from Exponents or Controller.
            res[3]=res[3]+1;}}
    // Check whether exist non-exponenet controller?
    if(PointedGene[0]>=0){// Exist!! Should reset all corresponding mapping tables.
        int po_id, po_va;
        for(ii=0; ii<PointedGene.size(); ++ii){
            po_id=PointedGene[ii];
            po_va=PointValues[po_id];
            std::vector<int> tmps(1<<InD[po_id],po_va);
            bn_list[po_id]=tmps;}}
    // Set a model for analyzing.
    // Default Boolean system
    DNS_Engaged amodel(sys_size, 2L, &in_edge, &ot_edge, &bn_list);
    std::vector<int> UnstaleType(sys_size);
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
        for(ii=0; ii<sys_size; ++ii){
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

//' C++ prototype function for BoolBioNet_CoreDyn()
//' @param aRealNet, list-object of real/silico genetic network
//' @param PointedGene, integer or name vector, which nodes should be fixed.
//' @param PointValues, pointed values, for free or controlled nodes.
//' @param InD, in-degree vector.
//' @param OtD, out-degree vector.
//' @param NodeDetailInfor, should return detail informations of node.
//' @param ReturnResidualNetwork, should return Residual Network >0: yes.
//' @param Times, recursive number of analysis.
//' @details Please see document of \link{BoolBioNet_CoreDyn}
// [[Rcpp::export]]
Rcpp::List c_CoreDynamicNode(Rcpp::List aRealNet, IntegerVector PointedGene, IntegerVector PointValues,
    IntegerVector InD, IntegerVector OtD, int NodeDetailInfor, int ReturnResidualNetwork, int Times){
    // Return values.
    Rcpp::List Returner;
    Rcpp::IntegerVector res(5);res[3]=0;res[4]=0;// 3: External (No input node), 4: Terminal (No output node)
    int tmp_res[3];
    // Convert List into appropriate data_frame.
    Rcpp::IntegerVector tmpdata;
    Rcpp::List tmplist1=aRealNet[1], tmplist2=aRealNet[2], tmplist3=aRealNet[3];
    int ii, jj, sys_size=tmplist1.size();
    std::vector<std::vector<int>> in_edge(sys_size);
    for(ii=0; ii<sys_size; ++ii){
        in_edge[ii].resize(0);}
    std::vector<std::vector<int>> ot_edge=in_edge;
    std::vector<std::vector<int>> bn_list=in_edge;
    for(ii=0; ii<sys_size; ++ii){
        if(OtD[ii]>0){// Outs not empty // if(Rcpp::internal::Rcpp_IsNA(tmpdata[0]))
            tmpdata=tmplist2[ii];
            for(jj=0; jj<tmpdata.size(); ++jj){
                ot_edge[ii].push_back(tmpdata[jj]);}}
        else {
            res[4]=res[4]+1;}
        if(InD[ii]>0){// Ins not empty // if(Rcpp::internal::Rcpp_IsNA(tmpdata[0]))
            tmpdata=tmplist1[ii];
            for(jj=0; jj<tmpdata.size(); ++jj){
                in_edge[ii].push_back(tmpdata[jj]);}
            tmpdata=tmplist3[ii];
            for(jj=0; jj<tmpdata.size(); ++jj){
                bn_list[ii].push_back(tmpdata[jj]);}}
        else {// No input, should be set one
            bn_list[ii].push_back(PointValues[ii]);// Actually this value could be from Exponents or Controller.
            res[3]=res[3]+1;}}
    // Check whether exist non-exponenet controller?
    if(PointedGene[0]>=0){// Exist!! Should reset all corresponding mapping tables.
        int po_id, po_va;
        for(ii=0; ii<PointedGene.size(); ++ii){
            po_id=PointedGene[ii];
            po_va=PointValues[po_id];
            std::vector<int> tmps(1<<InD[po_id],po_va);
            bn_list[po_id]=tmps;}}
    // Set a model for analyzing.
    // Default Boolean system
    DNS_CoreDyn amodel(sys_size, 2L, &in_edge, &ot_edge, &bn_list);
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
        for(ii=0; ii<sys_size; ++ii){
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

//' C++ prototype function for BoolFun_Polynomial()
//' @param VariableMat, a matrix that record input variables.
//' @param MapTab, the truth table of Boolean function.
//' @param LogiSpin, an integer, set the detail form of Polynomial transformation.
//' @details Please see document of \link{BoolFun_Polynomial}
//' @return a List, [[1]] is.SAT ture or not? [[2]] Weights of coupled variable; [[3]] threshold.
// [[Rcpp::export]]
Rcpp::List c_BoolFun2Polynomial(Rcpp::IntegerMatrix &VariableMat, Rcpp::IntegerVector &MapTab, int LogiSpin){
    return (PolynomialFunction(VariableMat, MapTab, LogiSpin));
}

//' C++ prototype function for BoolBioNet_StrConComp()
//' @param aRealNet, list-object of real/silico genetic network
//' @param InD, IntegerVector of in-degree
//' @param OtD, IntegerVector of out-degree
//' @details Please see document of \link{BoolBioNet_StrConComp}
//' @return a List, each element is a strong connect component.
// [[Rcpp::export]]
Rcpp::List c_StrongConnectComponent(Rcpp::List aRealNet, Rcpp::IntegerVector &InD, Rcpp::IntegerVector &OtD){
    // Return values.
    Rcpp::List Returner;
    // Convert List into appropriate data_frame.
    Rcpp::IntegerVector tmpdata;
    Rcpp::List tmplist1=aRealNet[1], tmplist2=aRealNet[2], tmplist3=aRealNet[3];
    int ii, jj, sys_size=tmplist1.size();
    std::vector<std::vector<int>> in_edge(sys_size);
    for(ii=0; ii<sys_size; ++ii){
        in_edge[ii].resize(0);}
    std::vector<std::vector<int>> ot_edge=in_edge;
    for(ii=0; ii<sys_size; ++ii){
        if(OtD[ii]>0){// Outs not empty // if(Rcpp::internal::Rcpp_IsNA(tmpdata[0]))
        //if(Rcpp::is_na(tmplist2[ii][0])){
            tmpdata=tmplist2[ii];
            for(jj=0; jj<tmpdata.size(); ++jj){
                ot_edge[ii].push_back(tmpdata[jj]);}}
        if(InD[ii]>0){// Ins not empty // if(Rcpp::internal::Rcpp_IsNA(tmpdata[0]))
        //if(Rcpp::is_na(tmplist1[ii][0])){
            tmpdata=tmplist1[ii];
            for(jj=0; jj<tmpdata.size(); ++jj){
                in_edge[ii].push_back(tmpdata[jj]);}}
    }
    NetGraphFrame aNet;// No need using class of {BNS_Basic}, {NetGraphFrame} enough!
    aNet.ConfigurationBuildNet('N',sys_size,-1,-1,-1.0).Build_GraphNet();
    aNet.LoadFromVecVecIntFrame(in_edge,ot_edge);
    std::vector<std::vector<int>> scc;
    aNet.Tarjon(scc);
    for(ii=0; ii<((int)(scc.size())); ++ii){
        Returner.push_back(scc[ii]);
    }
    // aNet.Delete_Network();// Clear network 
    return Returner;
}

//' C++ prototype function for MulVFun_Complexity()
//' @param avec, IntegerVector, each element in {0,1,...L-1}
//' @param k, input number (in-degree)
//' @param L, system of L:{v0,v1,...v_{L-1}}
//' @details Please see document of \link{MulVFun_Complexity}
//' @return numeric, the complexity of given function.
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

//' C++ prototype function for MulVFun_EffectiveEdges()
//' @param avec, IntegerVector, each element in {0,1,...L-1}
//' @param k, input number (in-degree)
//' @param L, system of L:{v0,v1,...v_{L-1}}
//' @details Please see document of \link{MulVFun_EffectiveEdges}
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

//' C++ prototype function for MulVFun_EffectiveEdges()
//' @param avec, IntegerVector, each element in {0,1,...L-1}
//' @param k, input number (in-degree)
//' @param L, system of L:{v0,v1,...v_{L-1}}
//' @details Please see document of \link{MulVFun_EffectiveEdges}
//[[Rcpp::export]]
Rcpp::NumericVector c_MulF_EffectiveEdges(IntegerVector &avec, int k, int L){
    int Lens=(int)pow(L,k), First=avec[0];
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

//' C++ prototype function for MulVFun_QMForm()
//' @param avec, IntegerVector, each element in {0,1,...L-1}
//' @param k, input number (in-degree)
//' @param L, system of L:{v0,v1,...v_{L-1}}
//' @details Please see document of \link{MulVFun_QMForm}
//' @return a matrix, rows are multi-valued prime implicants, columns are
//' each variable's state [0~L-1], -1 means "*". Last column is mapped results.
//[[Rcpp::export]]
IntegerMatrix c_MulF_QuineMcCluskey(IntegerVector &avec, int k, int L){
    int *ttt=avec.begin();
    IntegerMatrix res=toR_ShowMulVFunDNF(ttt,k,L);
    return res;
}

//' C++ prototype function for MulVFun_Generator()
//' @param FunType, function's type.
//' @param k, input variable size.
//' @param L, level of discrete system.
//' @param CanaDeep, layer of canalization.
//' @param CanaVar, Canalized varibale's ID.
//' @param CanaVarNum, Each canalizing number.
//' @param CanaInfo1, Rcpp::List, canalizing information.
//' @param CanaInfo2, Rcpp::List, canalized information.
//' @param bias, function's bias.
//' @param Cana_Free, is a parameter-free canalized function?
//' @details Please see document of \link{MulVFun_Generator}
//[[Rcpp::export]]
Rcpp::IntegerVector c_MulVF_Generator(char FunType,int k,int L, int CanaDeep, 
    Rcpp::IntegerVector &CanaVar, Rcpp::IntegerVector &CanaVarNum,
    Rcpp::List &CanaInfo1, Rcpp::List &CanaInfo2,
    Rcpp::NumericVector &bias, bool Cana_Free){
    unsigned short Lens=(unsigned short)(pow(L,k));
    short *bitmap=(short*)malloc(Lens*sizeof(short));
    char fun_type=FunType;
    double *BiasSet=bias.begin();
    // Set config slot.
    std::vector<int> Argus(4,0); Argus[0]=CanaDeep;
    mulvfun amvf(k, L, Lens, FunType);
    amvf.Configuration(BiasSet, bitmap, Argus);
    if(fun_type=='C'){
        if(Cana_Free){// User not offer detail information
            amvf.Gen_Cana_Free(CanaVar,CanaVarNum);
        }
        else {// User have provide information
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
                    CanaIN[ii].push_back(c_in[jj]);
                    CanaOT[ii].push_back(c_ot[jj]);
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
    for(int ii=0; ii<Lens; ++ii){
        xx[ii]=bitmap[ii];
    }
    free(bitmap);
    return xx;
}

//' C++ prototype function for MulVFun_Sensitivity()
//' @param amulfunss, a multi-valued function vector.
//' @param k, input variable size.
//' @param L, L-level discrete value.
//' @param Lens, length of mapping table: L^k.
//' @details Please see document of \link{MulVFun_Sensitivity}
//[[Rcpp::export]]
double c_MulVF_Sensitivity(Rcpp::IntegerVector &amulfunss, int k, int L, int Lens){
    mulvfun amvf(k, L, Lens, 'X');
    short *maptab=(short*)malloc(Lens*sizeof(short));
    for(int ii=0; ii<Lens; ++ii){
        maptab[ii]=amulfunss[ii];}
    std::vector<int> NoUse(2);
    double xx=amvf.Configuration(nullptr, maptab,NoUse).Sensitivity();
    free(maptab);
    amvf.Reset();
    return xx;
}

//' C++ prototype function for MulVFun_is_NestedCanalized()
//' @param amulfunss, a multi-valued function vector.
//' @param k, input variable size.
//' @param L, L-level discrete value.
//' @param Lens, length of mapping table: L^k.
//' @details Please see document of \link{MulVFun_is_NestedCanalized}
//[[Rcpp::export]]
Rcpp::List c_M_NestedCanalized(Rcpp::IntegerVector &amulfunss, int k, int L, int Lens){
    mulvfun amvf(k, L, Lens, 'X');
    short *maptab=(short*)malloc(Lens*sizeof(short));
    for(int ii=0; ii<Lens; ++ii){
        maptab[ii]=amulfunss[ii];}
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

//' C++ prototype function for MulVFun_is_Threshold()
//' @param amulfunss, a multi-valued function vector.
//' @param k, input variable size.
//' @param L, L-level discrete value.
//' @param Lens, length of mapping table: L^k.
//' @details Please see document of \link{MulVFun_is_Threshold}
//[[Rcpp::export]]
Rcpp::List c_M_Threshold(Rcpp::IntegerVector &amulfunss, int k, int L, int Lens){
    mulvfun amvf(k, L, Lens, 'T');
    short *maptab=(short*)malloc(Lens*sizeof(short));
    for(int ii=0; ii<Lens; ++ii){
        maptab[ii]=amulfunss[ii];}
    std::vector<int> NoUse(2,-1);
    Rcpp::List FinRes=amvf.Configuration(nullptr, maptab, NoUse).is_MulVF_Thre();
    amvf.Reset();
    free(maptab);
    return (FinRes);
}
//' C++ prototype function for MulVFun_is_Domainted()
//' @param amulfunss, a multi-valued function vector.
//' @param k, input variable size.
//' @param L, L-level discrete value.
//' @param Lens, length of mapping table: L^k.
//' @details Please see document of \link{MulVFun_is_Domainted}
//[[Rcpp::export]]
Rcpp::List c_M_Domainted(Rcpp::IntegerVector &amulfunss, int k, int L, int Lens){
    mulvfun amvf(k, L, Lens, 'D');
    short *maptab=(short*)malloc(Lens*sizeof(short));
    for(int ii=0; ii<Lens; ++ii){
        maptab[ii]=amulfunss[ii];}
    std::vector<int> NoUse(2,-1);
    Rcpp::List FinRes=amvf.Configuration(nullptr, maptab, NoUse).is_MulVF_Domi();
    amvf.Reset();
    free(maptab);
    return (FinRes);
}
//' C++ prototype function for MulVFun_is_Signed()
//' @param amulfunss, a multi-valued function vector.
//' @param k, input variable size.
//' @param L, L-level discrete value.
//' @param Lens, length of mapping table: L^k.
//' @details Please see document of \link{MulVFun_is_Signed}
//[[Rcpp::export]]
Rcpp::List c_M_Signed(Rcpp::IntegerVector &amulfunss, int k, int L, int Lens){
    mulvfun amvf(k, L, Lens, 'S');
    short *maptab=(short*)malloc(Lens*sizeof(short));
    for(int ii=0; ii<Lens; ++ii){
        maptab[ii]=amulfunss[ii];}
    std::vector<int> NoUse(2,-1);
    Rcpp::List FinRes=amvf.Configuration(nullptr, maptab, NoUse).is_MulVF_Sign();
    amvf.Reset();
    free(maptab);
    return (FinRes);
}

//' C++ prototype function for MulVFun_Polynomial()
//' @param VariableMat, a matrix that record input variables.
//' @param MapTab, the truth table of multi-valued function.
//' @param k, input variable size.
//' @param L, L-level discrete value.
//' @details Please see document of \link{MulVFun_Polynomial}
//' @return a List, [[1]] is.SAT ture or not? [[2]] Weights of coupled variable; [[3]] thresholds of L-levels.
// [[Rcpp::export]]
Rcpp::List c_MulVFun2Polynomial(Rcpp::IntegerMatrix &VariableMat, Rcpp::IntegerVector &MapTab, int k, int L){
    return (MulVFun_PolynomialFunction(VariableMat, MapTab, k, L));
}

//' C++ prototype function for MulV2Bool_Bool2MulV()
//' @param OriMapTab, IntegerVector of original mapping table.
//' @param k, input variable size.
//' @param L, L-level discrete value.
//' @param Thresholds, IntegerVector of binary thresholds for multi-valued system.
//' @param b2m, Boolean value, Is from Boolean to Multi-valued system?
//' @details Please see document of \link{MulV2Bool_Bool2MulV}
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

// Code is over!
