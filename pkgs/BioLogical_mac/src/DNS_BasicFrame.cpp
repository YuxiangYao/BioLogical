#include "DNS_BasicFrame.h"

// DNS auxiliary function: generate self-defined systems
void DNS_Aux_GenerateTopo(char NetType, double NetPara, int Size, 
    std::vector<std::vector<int>> *inedge,// [!!warning] Memory allocated over.
    std::vector<std::vector<int>> *otedge, int *InDegs){
    NetGraphFrame aNet;
    if('E'==NetType){//Erdös-Rényi network
        aNet.ConfigurationBuildNet(NetType, Size, -666, -666, NetPara);
    }
    else {// Kauffman, lattice, and regular random networks
        aNet.ConfigurationBuildNet(NetType, Size, (int)(NetPara), -666, NetPara);
    }
    // Build network and output into two vector<vector<int>>
    aNet.Build_GraphNet().Out2VecVecIntFrame(inedge, otedge);
    memcpy(InDegs, aNet.InDeg, Size*sizeof(int));
}

void DNS_Aux_GenerateFunc(int sys_size, int ll_sys, int *InDeg, int *IDsss,
    double *bias, char of_type, double of_part,
    int of_i_para1, int of_i_para2, std::vector<std::vector<int>> *bnsss,
    std::vector<int> &ControlNodes, std::vector<int> &ControlValues){
    (*bnsss).resize(sys_size);
    //std::shuffle(IDsss, IDsss+sys_size, mt);
    int idss, vals;
    for(int ii=sys_size-1; ii>=0; ii=ii-1){
        idss=(int)(unif_rand()*ii);
        vals=IDsss[ii];
        IDsss[ii]=IDsss[idss];
        IDsss[idss]=vals;}
    int part=(int)(sys_size*of_part);
    if(2==ll_sys){// Boolean system 
        int ii, jj, bn_lens, tmp_ind;
        bool *bitmap=(bool*)malloc((1<<17)*sizeof(bool));// Enough large. (Actually k-input ≤ 16)
        boolfun a_bf;
        // Ordered part.
        for(ii=0; ii<part; ++ii){// OBF part.
            tmp_ind=InDeg[IDsss[ii]];
            if(0<tmp_ind){// Has input-edges
                a_bf.Configuration(of_type, tmp_ind, 0.50, bitmap, of_i_para1, of_i_para2, -9, -9);
                a_bf.Gen_BF().Reset();
                (*bnsss)[IDsss[ii]].resize(1<<tmp_ind);
                for(jj=0; jj<(1<<tmp_ind); ++jj){// Bool to int
                    (*bnsss)[IDsss[ii]][jj]=(bitmap[jj]>0);
                }//(*bnsss)[IDsss[ii]].assign(bitmap,bitmap+(1<<tmp_ind));
            }
            else {// No input, should be only one fixed value.
                (*bnsss)[IDsss[ii]].push_back(unif_rand()>0.50);}}
        free(bitmap); bitmap=nullptr;
        // Random part.
        for(ii=part; ii<sys_size; ++ii){// Random part.
            bn_lens=1<<InDeg[IDsss[ii]];
            (*bnsss)[IDsss[ii]].resize(bn_lens);
            for(jj=0; jj<bn_lens; ++jj){
                (*bnsss)[IDsss[ii]][jj]=unif_rand()<bias[0];
            }
        }
    }
    else {// Multi-valued system,only allow for c('R','C','D','T')
        int ii, jj, bn_lens, tmp_ind;
        std::vector<short> bitmap;
        // Ordered part.
        std::vector<int> Argus(4,0);
        Argus[0]=of_i_para1;
        Argus[1]=of_i_para2;
        Rcpp::IntegerVector cana_config1(1,-1);
        Rcpp::IntegerVector cana_config2(Argus[0],Argus[1]);
        for(ii=0; ii<part; ++ii){// OBF part.
            tmp_ind=InDeg[IDsss[ii]];
            if(0<tmp_ind){// Has input-edges
                bn_lens=(int)pow(ll_sys,tmp_ind);
                bitmap.resize(bn_lens);
                mulvfun amvf(tmp_ind, ll_sys, bn_lens, of_type);
                amvf.Configuration(bias, bitmap.data(), Argus);
                switch(of_type){
                    case 'C':amvf.Gen_Cana_Free(cana_config1, cana_config2);break;
                    case 'T':amvf.Gen_MulVF_Thre();break;
                    case 'D':amvf.Gen_MulVF_Domi();break;
                    case 'R':amvf.Gen_Rand();break;
                    default:break;
                }
                (*bnsss)[IDsss[ii]].resize(bn_lens);
                for(jj=0; jj<bn_lens; ++jj){// Bool to int
                    (*bnsss)[IDsss[ii]][jj]=(int)bitmap[jj];
                }
            }
            else {// No input, should be only one fixed value. (can be random)
                (*bnsss)[IDsss[ii]].push_back((int)(unif_rand()*ll_sys));
            }
        }
        // Random part.
        for(ii=part; ii<sys_size; ++ii){// Random part.
            tmp_ind=InDeg[IDsss[ii]];
            if(0<tmp_ind){// Has input-edges
                bn_lens=(int)pow(ll_sys,tmp_ind);
                bitmap.resize(bn_lens);
                mulvfun amvf(tmp_ind, ll_sys, bn_lens, of_type);
                amvf.Configuration(bias, bitmap.data(), Argus).Gen_Rand();
                (*bnsss)[IDsss[ii]].resize(bn_lens);
                for(jj=0; jj<bn_lens; ++jj){// Bool to int
                    (*bnsss)[IDsss[ii]][jj]=(int)bitmap[jj];
                }
            }
            else {// No input, should be only one fixed value.
                short tmp_f;
                DNS_Aux_Bias(&tmp_f, 1L, bias, ll_sys);
                (*bnsss)[IDsss[ii]].push_back((int)(tmp_f));
            }
        }
    }
    if(ControlNodes.size()>0){// Due to index>=0, a minus number means nonexisted controller.
        for(int ii=0; ii<((int)ControlNodes.size()); ++ii){
            int code=ControlNodes[ii];
            std::vector<int> tmps((int)pow(ll_sys, InDeg[code]), ControlValues[code]);
            (*bnsss)[code]=tmps;
        }
    }
}

// DNS auxiliary function: assign states according to configured biases//[v]
void DNS_Aux_Bias(short *slot, int length, double *p, int Ls){
    if(2==Ls){// Boolean system, quick setting.
        for(int ii=0; ii<length; ++ii){
            slot[ii]=(short)(unif_rand()<p[0]);
        }
    }
    else {// Multi-valued setting.
        short jj;
        double Rander;
        bool flag;
        for(int ii=0; ii<length; ++ii){
            Rander=unif_rand();
            flag=true;
            while(flag){
                for(jj=0; jj<(short)Ls; ++jj){
                    if(Rander<p[jj]){
                        slot[ii]=jj; 
                        flag=false;
                        break;
                    }
                    else if(jj==(Ls-1)){// <- As a safety setting, avoid sum<1
                        slot[ii]=jj; 
                        flag=false;
                        break;
                    }
                    else {
                        Rander-=p[jj];
                    }
                }
            }
        }
    }
}

// DNS auxiliary function: Set a different state (MOD method)//[v]
void DNS_Aux_BitDiff(short *ss, int L){
    if(2==L){// only 0 and 1, inverse them
        *ss=1-(*ss);
    }
    else {
        *ss=(*ss+(short)(unif_rand()*(L-1)))%L;
    }
}

// DNS_Basic Constructor 1: Default//[v]
DNS_Basic::DNS_Basic() {// N(0), L(2), inedge(nullptr), otedge(nullptr), 
}// bns(nullptr), sss(nullptr), Labels(nullptr), LL_sys(nullptr)

// DNS_Basic Constructor 2: with arguments//[v]
DNS_Basic::DNS_Basic(int size, int L): N(size), L(L){
    Labels=(int*)malloc(size*sizeof(int));// randomly configure node's IDs.
    for(int ii=0; ii<size; ++ii){
        Labels[ii]=ii;
    }
    LL_sys=(int*)malloc(20*sizeof(int));
    LL_sys[0]=1;
    for(int ii=0; ii<19; ++ii){// [!!warning] Large L may Overflow int32.
        LL_sys[ii+1]=LL_sys[ii]*L;// Here not check their rationality.
    }
}

// DNS_Basic Destructor//[v]
DNS_Basic::~DNS_Basic() {
    inedge=nullptr;
    otedge=nullptr;
    bns=nullptr;
    free(Labels);
    Labels=nullptr;
    free(LL_sys);
    LL_sys=nullptr;
}

// Synchronous update//[v]
DNS_Basic& DNS_Basic::SystemEvolution_Syn(int steps, short *ss0, std::vector<bool> *flag){
    int ii, jj, kk, index;
    short *olds=ss0, *news=ss0+N, *changer;
    for(ii=0; ii<steps; ++ii){
        for(jj=0; jj<N; ++jj){
            index=0;
            for(kk=0; kk<((int)(*inedge)[jj].size()); ++kk){
                index+=(olds[(*inedge)[jj][kk]]*LL_sys[kk]);
            }// index+=(olds[(*inedge)[jj][kk]]<<kk);// only for Boolean system
            news[jj]=(*bns)[jj][index];
            (*flag)[jj]=(*flag)[jj]||(olds[jj]!=news[jj]);
        }
        changer=olds; olds=news; news=changer;
    }
    return *this;
}

// Asynchronous update 1: random select//[v]
DNS_Basic& DNS_Basic::SystemEvolution_Asy1(int steps, short *ss0, std::vector<bool> *flag){
    int ii, jj, kk, index;
    short tmp;
    for(ii=0; ii<steps; ++ii){
        jj=(int)(unif_rand()*N);// Random select one.
        index=0;
        for(kk=0; kk<((int)(*inedge)[jj].size()); ++kk){
            index+=(ss0[(*inedge)[jj][kk]]*LL_sys[kk]);
        }
        tmp=ss0[jj];
        ss0[jj]=(*bns)[jj][index];// Just update this.
        (*flag)[jj]=(*flag)[jj]||(ss0[jj]!=tmp);
    }
    return *this;
}

// Asynchronous update 2: random select valid one//[v]
DNS_Basic& DNS_Basic::SystemEvolution_Asy2(int steps, short *ss0, std::vector<bool> *flag){
    int ii, jj, kk, index, changed;
    std::set<int> candidate;     // Record all unstable nodes' ID.
    for(ii=0; ii<N; ++ii){
        index=0;
        for(jj=0; jj<((int)(*inedge)[ii].size()); ++jj){
            index+=(ss0[(*inedge)[ii][jj]]*LL_sys[jj]);
        }
        if(ss0[ii]!=(*bns)[ii][index]){// Not as stable states.
            candidate.insert(ii);
        }
    }
    for(ii=0; (ii<steps)&&(candidate.size()>0); ++ii){
        jj=(int)(unif_rand()*candidate.size());
        auto it = candidate.begin();
        std::advance(it, jj);
        jj=(*it);
        index=0;
        for(kk=0; kk<((int)(*inedge)[jj].size()); ++kk){
            index+=(ss0[(*inedge)[jj][kk]]*LL_sys[kk]);
        }
        ss0[jj]=(*bns)[jj][index];// Update this jj node.
        (*flag)[jj]=(true)||(*flag)[jj];
        // jj's succssors maybe turn from unstble to stable or inverse.
        changed=jj;
        for(int pp=0; pp<((int)(*otedge)[changed].size()); ++pp){
            jj=(*otedge)[changed][pp];
            index=0;
            for(kk=0; kk<((int)(*inedge)[jj].size()); ++kk){
                index+=(ss0[(*inedge)[jj][kk]]*LL_sys[kk]);
            }
            if(ss0[jj]==(*bns)[jj][index]){// Keep stable
                candidate.erase(jj);
            }
            else {// Still unstable
                candidate.insert(jj);
            }
        }
    }
    return *this;
}

// Load the topological structure and mapping tables//[v]
DNS_Basic& DNS_Basic::LoadedModel(std::vector<std::vector<int>> *load_InDeg, 
  std::vector<std::vector<int>> *load_OtDeg, std::vector<std::vector<int>> *Bnsss){ 
    inedge=load_InDeg;      // Info of in-degree
    otedge=load_OtDeg;      // Info of out-degree
    bns=Bnsss;              // Info of truth/mappin tables
    return *this;
}

// DNS_Derrida's Constructor & Destructor//[V]
DNS_Derrida::DNS_Derrida(): DNS_Basic(){
}
DNS_Derrida::DNS_Derrida(int size, int L): DNS_Basic(size, L){
    sss=(short *)malloc(4*size*sizeof(short)); // Store temporary state
}
DNS_Derrida::~DNS_Derrida(){
    free(sss);
    sss=nullptr;
}

// Execute once simulation of damage-spread//[v]
DNS_Derrida& DNS_Derrida::DerridaDamageSpread(int Steps, int UpdateType){
    std::vector<bool> tmp(N,false);// Useless here
    std::vector<bool> *pb=&tmp;
    switch(UpdateType){
        case 1: // Synchronous
            DNS_Basic::SystemEvolution_Syn(Steps,sss,pb);
            DNS_Basic::SystemEvolution_Syn(Steps,sss+N+N,pb);
            break;
        case 2: // Asynchronous random one
            DNS_Basic::SystemEvolution_Asy1(Steps,sss,pb);
            DNS_Basic::SystemEvolution_Asy1(Steps,sss+N+N,pb);
            break;
        case 3: // Asynchronous must valid one
            DNS_Basic::SystemEvolution_Asy2(Steps,sss,pb);
            DNS_Basic::SystemEvolution_Asy2(Steps,sss+N+N,pb);
            break;
        default:;
    }
    return *this;
}

// Configure some dynamic parameters//[v]
DNS_Derrida& DNS_Derrida::DynamicParaConfig(double *s_bias,double s_delta){
    f_argus[9]=s_delta;// [9] s_delta (normal Hamming distance)
    short *ipt1=sss, *ipt2=sss+N+N;
    DNS_Aux_Bias(ipt1, N, s_bias, L);       // Set original system states
    memcpy(ipt2, ipt1, N*sizeof(short));    // Copy patterns
    //std::shuffle(Labels, Labels+N, mt);     // Shuffle codes to set "delta"
    int idss, vals;
    for(int ii=N-1; ii>=0; ii=ii-1){
        idss=(int)(unif_rand()*ii);
        vals=Labels[ii];
        Labels[ii]=Labels[idss];
        Labels[idss]=vals;}
    for(int ii=0; ii<(int)(f_argus[9]*N); ++ii){
        DNS_Aux_BitDiff(ipt2+Labels[ii], L);
    }
    return *this;
}

// Return the final Hamming distance//[v]
double DNS_Derrida::FinalDistance(){
    int ii, sum=0;
    short *ipt1=sss, *ipt2=sss+N+N;
    for(ii=0; ii<N; ii++){
        sum+=(ipt1[ii]!=ipt2[ii]);}
    return sum/(double)N;
}


// DNS_Percolation's Constructor & Destructor//[V]
DNS_Percolation::DNS_Percolation(): DNS_Basic(){
}
DNS_Percolation::DNS_Percolation(int size, int L, int LaType): 
    DNS_Basic(size, L), Lattice(LaType) {
    sss=(short *)malloc(2*size*sizeof(short)); // Store temporary state information
    flag=new std::vector<bool>(size);
}
DNS_Percolation::~DNS_Percolation(){
    free(sss);
    sss=nullptr;
    delete flag;
    flag=nullptr;
}

// Configure some dynamic parameters//[v]
DNS_Percolation& DNS_Percolation::DynamicParaConfig(double *s_bias){
    DNS_Aux_Bias(sss, N, s_bias, L);       // Set original system states
    return *this;
}

// Execute once simulation of damage-spread in lattice//[v]
DNS_Percolation& DNS_Percolation::PercolationModel(int Steps, int Windows, int UpdateType){
    std::vector<bool> tmp(N,false);// Useless here
    std::vector<bool> *pb=&tmp;// [REF] J. theor. Biol. (1988) 135, 255-261]
    switch(UpdateType){
        case 1: // Synchronous
            DNS_Basic::SystemEvolution_Syn(Steps,sss,pb);
            DNS_Basic::SystemEvolution_Syn(Windows,sss,flag);
            break;
        case 2: // Asynchronous random one
            DNS_Basic::SystemEvolution_Asy1(Steps,sss,pb);
            DNS_Basic::SystemEvolution_Asy1(Windows,sss,flag);
            break;
        case 3: // Asynchronous must valid one
            DNS_Basic::SystemEvolution_Asy2(Steps,sss,pb);
            DNS_Basic::SystemEvolution_Asy2(Windows,sss,flag);
            break;
        default:;
    }
    return *this;
}

// Judge the max cluster exist? Percolation path through boundaries//[v]
DNS_Percolation& DNS_Percolation::PercolationWhetherNot(){
    int LengHigh=(int)(sqrt(N));// Length of (only) tetragonal lattice
    int ii, jj, Boundary=N-LengHigh, tmpLOC;
    // Split the top+bottom boundary.
    NetGraphFrame aNet;
    aNet.ConfigurationBuildNet('L', N, Lattice,-666, -666.0).Build_GraphNet();
    if(4==Lattice){// Square
        for(ii=0; ii<LengHigh; ++ii){
            aNet.BreakDownPointedEdge_D(ii,Boundary+ii).BreakDownPointedEdge_D(Boundary+ii,ii);
        }
    }
    else if(3==Lattice){// Triangle
        for(ii=1; ii<LengHigh; ii+=2){
            aNet.BreakDownPointedEdge_D(ii,Boundary+ii).BreakDownPointedEdge_D(Boundary+ii,ii);
        }
    }
    else {// Hexagon
        for(ii=0; ii<LengHigh; ++ii){
            aNet.BreakDownPointedEdge_D(ii,Boundary+ii).BreakDownPointedEdge_D(Boundary+ii,ii);
            tmpLOC=ii-1;
            if(tmpLOC<0){
                tmpLOC=LengHigh-1;
            }
            aNet.BreakDownPointedEdge_D(ii,Boundary+tmpLOC).BreakDownPointedEdge_D(Boundary+tmpLOC,ii);
        }
    }
    // Delete unstable nodes
    for(ii=0; ii<N; ++ii){
        if((*flag)[ii]){// Isolate unstable nodes.
            aNet.IsolatedPointedNode_D(ii);
        }
    }
    // Analyze whether percolation and largest cluster.
    int *routers=aNet.InDeg+N+N+N;      // 
    int *topline=routers+N-LengHigh;    // [NOTE] Length≠High can lead to error
    bool passthrough;
    std::vector<bool> ReduceLabel(LengHigh);
    for(ii=0; ii<LengHigh; ++ii){// Label the bottom nodes.
        ReduceLabel[ii]=!((*flag)[ii]); // flag's TRUE is "Oscillation"
    }
    std::set<int> MaxCluster; MaxCluster.clear();
    std::set<int> tmpCluster;
    for(ii=0; ii<LengHigh; ++ii){// Check nodes within bottom-line
        if(ReduceLabel[ii]){// It is a stable node
            aNet.ShortestPath(ii);  // The routing saved in [aNet.InDeg+N*3]
            passthrough=false;      // Find the shortest way from node-ii
            for(jj=0; jj<LengHigh; ++jj){
                if(topline[jj]<N){  // Has a path? (Only check top line)
                    passthrough=true;
                    break;
                }
            }
            tmpCluster.clear();
            if(passthrough){// Find all nodes in this path/cluster
                for(jj=0; jj<N; ++jj){
                    if(routers[jj]<N){
                        tmpCluster.insert(jj);
                    }
                }
                for(jj=0; jj<LengHigh; ++jj){
                    if(routers[jj]<N){// Avoid repetitive computation
                        ReduceLabel[jj]=false;
                    }
                }
            }
            if(tmpCluster.size()>MaxCluster.size()){
                MaxCluster=tmpCluster;
            }
        }
    }
    MaxCluster_ID.clear();
    MaxCluster_ID.insert(MaxCluster_ID.end(), MaxCluster.begin(), MaxCluster.end());
    return *this;
}

// Output fractions of stable nodes & size of max stable-cluster//[v]
DNS_Percolation& DNS_Percolation::PercolationStableFraction(double *Returner){
    int sum=0;
    for(int ii=0; ii<N; ++ii){
        sum+=(*flag)[ii];
    }
    i_argus[8]=Returner[1]=N-sum;
    i_argus[9]=Returner[0]=(int)MaxCluster_ID.size();// MaxCluster size
    return *this;
}

// Output final lattice status (stable-0, stable-1, unstable)//[v]
DNS_Percolation& DNS_Percolation::OutputFinalLattice(int *IsFixed, int *MaxCluster){
    for(int ii=0; ii<N; ++ii){
        if((*flag)[ii]){// Unstable: TRUE ~~-> Oscillation
            IsFixed[ii]=-1;
        }
        else {// Stable
            IsFixed[ii]=(int)sss[ii];
        }
    }
    for(int ii=0; ii<(int)MaxCluster_ID.size(); ++ii){
        MaxCluster[MaxCluster_ID[ii]]=1;
    }
    return *this;
}


// Auxiliary function: count the number of "1" in a int32 variable.
int Aux_Num1inInt32(int aInt32Num){
    int count=0, tmp=aInt32Num;
    while(tmp){
        count++;
        tmp&=(tmp-1);}
    return count;
}

// Auxiliary function: scalar product.
int Aux_DotProductSummation(int *Vec1,int *Vec2,int Length){
    int ii,sum=0;
    for(ii=0; ii<Length; ++ii){
        sum+=Vec1[ii]*Vec2[ii];}
    return sum;
}

// Auxiliary function: Simplify mapping table due to some fixed nodes. (By bit-comparasion)
void Aux_SimplifyMappingTable(int *oldmap, int inDeg, std::vector<int> &Bit_0stable, std::vector<int> &Bit_1stable){
    int valid_ind=inDeg-Bit_0stable.size()-Bit_1stable.size();
    int ii, num=0, checker0=0, checker1=0;
    int *newmap=(int *)malloc((1<<valid_ind)*sizeof(int));// newmap[1<<valid_ind];
    for(ii=0; ii<((int)Bit_0stable.size()); ii++){// Set 0-stable check
        checker0+=(1<<Bit_0stable[ii]);}
    for(ii=0; ii<((int)Bit_1stable.size()); ii++){// Set 1-stable check
        checker1+=(1<<Bit_1stable[ii]);}
    for(ii=0; ii<(1<<inDeg); ii++ ){
        if((checker1==(checker1&ii))&&(checker0==(checker0&(~ii)))){
            newmap[num++]=oldmap[ii];}}
    memcpy(oldmap,newmap,(1<<valid_ind)*sizeof(int));
    free(newmap);
}

// Auxiliary function: Simplify one node's unstable parents due to some fixed nodes.
void Aux_SimplifyParents(int *temp_par, int inDeg, int del_code){
    for(int ii=del_code; ii<(inDeg-1); ii++){
        temp_par[ii]=temp_par[ii+1];}
}

// Auxiliary function: Two upstream points coupling, which lead to some stable node, invalid edges.
void Aux_TwoUpperNodeCoupling(int K_var, int ID1, int ID2, bool IsInCoherent, int *maptab){
    // Note that only simplify variables not involve ramaining ones.
    // {A,B,C} -> D & A->B: bf: 1X1X0X0X, ABC as AC via 1100 (This means C is invalidt that can be evalued at next loop)
    bool logi=false;
    int ii, tar1=1<<(K_var-1-ID1), tar2=1<<(K_var-1-ID2), lengths=(1<<K_var);
    int *bn_list0=(int *)malloc((lengths>>1)*sizeof(int));//int bn_list0[lengths>>1];
    int count=0;
    // Normal order or delete lower bits. The turth table can directly copy.
    if(ID1<ID2){// Delete the lower bit (ID2); Smaller ID is higher bit.
        for(ii=0; ii<lengths; ++ii){
            logi=(((ii&tar1)==0)&&((ii&tar2)==0))||(((ii&tar1)!=0)&&((ii&tar2)!=0));
            if(logi!=IsInCoherent){// Is Incoherent check. A --> B or A --| B ?
                bn_list0[count++]=maptab[ii];}}
        // if(count!=(lengths>>1)){
        //     Rcpp::stop("\nError!\nHere.\n(%d,%d)\n(InDeg: %d,ID1: %d, ID2: %d)\n",
        //         count,lengths>>1,K_var,ID1,ID2);}
        memcpy(maptab, bn_list0, (lengths>>1)*sizeof(int));}
    else {// Delete 
        int Fuller=(1<<K_var)-1;
        int LowerBitGet=tar2-1;
        int HigherBitGet=Fuller^(tar2+LowerBitGet);
        int *bn_list_loc=(int*)malloc((lengths>>1)*sizeof(int));//int bn_list_loc[lengths>>1];
        for(ii=0; ii<lengths; ++ii){
            logi=(((ii&tar1)==0)&&((ii&tar2)==0))||(((ii&tar1)!=0)&&((ii&tar2)!=0));
            if(logi!=IsInCoherent){// Is Incoherent check. A --> B or A --| B ?
                bn_list0[count]=maptab[ii];
                bn_list_loc[count]=((ii&HigherBitGet)>>1)+(ii&LowerBitGet);// New locations
                count++;}}
        for(ii=0; ii<(lengths>>1); ++ii){
            maptab[bn_list_loc[ii]]=bn_list0[ii];}
        free(bn_list_loc);}
    free(bn_list0);
}

// Auxiliary function: autoadder of a numerical system 
void Aux_NumSysAdder(int *vec, int k, int L){
    int loc=0, Ceiling=L-1;
    while(loc<k){
        if(vec[loc]<Ceiling){vec[loc]++; break;}
        else {vec[loc]=0; loc++;}}
}

// Auxiliary function: an int32 that records possible states are converted into a vector.
void Aux_Convert2(int MulitState, int *Vec, int NumSys){
    int ii, flag=MulitState;
    for(ii=0; ii<NumSys; ++ii){
        Vec[ii]=(flag&1);
        flag>>=1;}
}

// Auxiliary function: possible states under pointed mapping tables.
int Aux_Map2PossibleValue(int *maps, int length, int OriMulState, int NumSys){
    int ii,logi=0;
    std::vector<int> nx_bit(NumSys,1);
    for(ii=0; ii<NumSys; ++ii){
        nx_bit[ii]=(nx_bit[ii]<<ii);}
    for(ii=0; ii<length; ++ii){
        logi=logi|nx_bit[maps[ii]];
        if(logi==OriMulState)break;}// Have same possible states.
    return logi;
}


// DNS_Engaged's Constructor & Destructor//[V]
DNS_Engaged::DNS_Engaged(): DNS_Basic(){
}
DNS_Engaged::DNS_Engaged(int size, int L, std::vector<std::vector<int>> *ipt1,
    std::vector<std::vector<int>> *ipt2, std::vector<std::vector<int>> *ipt3): 
    DNS_Basic(size, L){
    inedge=ipt1;
    otedge=ipt2;
    bns=ipt3;
    PossibleLocal=(int*)malloc(3*sizeof(int)*N);// Remain possible locals   
    PossibleValue=PossibleLocal+N;  // Remain possible number of discrete values
    Pruned=PossibleValue+N;         // Label: Whether the node is invalid
    InDeg_temp=(int*)malloc(N*sizeof(int));
    OtDeg_temp=(int*)malloc(N*sizeof(int));
    Parents_temp=(int**)malloc(N*sizeof(int *));
    Children_temp=new std::vector<std::set<int>>(N);
    Mapping=(int**)malloc(N*sizeof(int *));
    GlobalCandidate.clear();
    StableNode.clear();
    UselessNode.clear();
    // Configure all information
    int ii, jj, tmp=(1<<L)-1, temp_bn_len;
    for(ii=0; ii<N; ++ii){
        InDeg_temp[ii]=(int)(*inedge)[ii].size();
        OtDeg_temp[ii]=(int)(*otedge)[ii].size();
        PossibleLocal[ii]=tmp;
        PossibleValue[ii]=L;
        // Should be set, not redundant steps, avoid zero-input or zero-output.
        Mapping[ii]=nullptr;
        Parents_temp[ii]=nullptr;
        GlobalCandidate.insert(ii);// Add all node for checking.
        (*Children_temp)[ii].clear();
        // Set temp_parents
        if(InDeg_temp[ii]>0){
            Parents_temp[ii]=(int*)malloc(InDeg_temp[ii]*sizeof(int));
            for(jj=0; jj<InDeg_temp[ii]; jj++){
                Parents_temp[ii][jj]=(*inedge)[ii][jj];
            }
        }
        // Set temp_children
        for(jj=0; jj<OtDeg_temp[ii]; jj++){
            (*Children_temp)[ii].insert((*otedge)[ii][jj]);
        }
        // Set temp_mapping_table
        temp_bn_len=LL_sys[InDeg_temp[ii]];// [WARNING] Not larger than 15.
        Mapping[ii]=(int*)malloc(temp_bn_len*sizeof(int));
        for(jj=0; jj<temp_bn_len; jj++){
            Mapping[ii][jj]=(*bns)[ii][jj];
        }
    }
}
DNS_Engaged::~DNS_Engaged(){
    free(PossibleLocal);
    PossibleLocal=nullptr;
    PossibleValue=nullptr;
    Pruned=nullptr;
    if(nullptr!=InDeg_temp){
        free(InDeg_temp);
        InDeg_temp=nullptr;
    }
    if(nullptr!=OtDeg_temp){
        free(OtDeg_temp);
        OtDeg_temp=nullptr;
    }
    if(Mapping==nullptr||Parents_temp==nullptr||Children_temp==nullptr);
    else {
        for(int ii=0; ii<N; ++ii){
            free(Mapping[ii]);
            free(Parents_temp[ii]);
            (*Children_temp)[ii].clear();
        }
        free(Mapping);
        Mapping=nullptr;
        free(Parents_temp);
        Parents_temp=nullptr;
        delete Children_temp;
        Children_temp=nullptr;
    }
}

// Remove invalid input due to stable point or fixed signals//[v]
void DNS_Engaged::RemoveInvalidInputs(int Code){
    int ii,jj,kk,gg;
    int *Standard, *Cursor, *Map=Mapping[Code];
    int Upper, Lower, logi;
    int step_inner, step_outer;
    int FreeK=0, RemainK=0;
    int indeg=InDeg_temp[Code];//int Independ[indeg],Depends[indeg],new_Pars[indeg];
    int *Independ=(int*)malloc(3*indeg*sizeof(int));
    int *Depends=Independ+indeg;
    int *new_Pars=Depends+indeg;
    for(ii=0; ii<indeg; ++ii){// The tmpI-th variable.
        Standard=Map;
        step_inner=Lower=LL_sys[ii];
        step_outer=LL_sys[ii+1];
        Upper=LL_sys[indeg-ii-1];// Upper=Order[indeg]/step_outer;
        logi=0;
        for(jj=0; jj<Upper; ++jj){// Current variable.
            for(kk=0; kk<Lower; ++kk){// Lower variable.
                Cursor=Standard+step_inner;
                for(gg=1; gg<L; ++gg){// 
                    if(Standard[kk]!=Cursor[kk]){
                        logi=1;break;
                    }
                    Cursor+=step_inner;
                }
            }
            if(logi){
                break;}
            Standard+=step_outer;}
        if(logi){// logi=1: depend this variable
            Depends[RemainK++]=ii;
        }
        else {// These variables 
            Independ[FreeK++]=ii;
        }
    }
    // Simplify the mapping table, in/out-degree, and temporary Parents/Children.
    if(FreeK){// Have "FreeK" independent variables, the indexes not the real codes.
        // Only change the *_temp slots.
        int *loc1=(int*)malloc(LL_sys[RemainK]*sizeof(int));// Slot for new mapping table.
        InDeg_temp[Code]-=FreeK;// Update code's in-degree.
        // Which parents can be free (Index) ?
        for(ii=0; ii<FreeK; ++ii){
            //jj=Parents_temp[Code][Independ[ii]];
            jj=Parents_temp[Code][indeg-1-Independ[ii]];// <-- Inserve the bit!
            OtDeg_temp[jj]--;
            (*Children_temp)[jj].erase(Code);// Delete child
        }
        // Some parents should be fixed.
        if(RemainK>0){
            for(ii=0; ii<RemainK; ++ii){
                //new_Pars[ii]=Parents_temp[Code][Depends[ii]];//[X] Due to [0]->[k] from high bit to low bit
                new_Pars[RemainK-1-ii]=Parents_temp[Code][indeg-1-Depends[ii]];// <-- Inserve the bit!
                Depends[ii]=LL_sys[Depends[ii]];
            }
        }// NumSys Order of ii-th remain variable.
        else {
            Depends[0]=0;}
        // Update mapping table.
        memset(Independ,0,indeg*sizeof(int));// Set 0 for updating mapping table.
        if(RemainK>0){
            for(ii=0; ii<LL_sys[RemainK]; ++ii){
                loc1[ii]=Map[Aux_DotProductSummation(Depends,Independ,RemainK)];
                Aux_NumSysAdder(Independ,indeg,L);}}// Automatic adder.
        else {
            loc1[0]=Map[0];}
        memcpy(Map, loc1, LL_sys[RemainK]*sizeof(int));// Insert into old mapping slot
        // Update new parents.
        if(RemainK>0){
            memcpy(Parents_temp[Code],new_Pars,RemainK*sizeof(int));
        }
        free(loc1);
    }
    free(Independ);
}

// Predecessor's states can lead to [Code] be shrunk//[v]
int DNS_Engaged::NodeShrunk(int Code){
    int ii, jj, flag, New_State, Old_State=PossibleLocal[Code];
    int indeg=InDeg_temp[Code];// int magnitude=Order[InDeg_temp[Code]];
    int *TempMappingVec=(int*)malloc(indeg*sizeof(int));//int TempMappingVec[indeg];
    memset(TempMappingVec,0,indeg*sizeof(int));
    int *ValidVariable=(int*)malloc(L*indeg*sizeof(int));
    memset(ValidVariable, 0, L*indeg*sizeof(int));
    int *ipt_map=Mapping[Code];
    int **local=(int**)malloc(indeg*sizeof(int*));// int *local[indeg];
    // Set the address of each variable.
    local[0]=ValidVariable;
    for(ii=1; ii<indeg; ++ii){
        local[ii]=local[ii-1]+L;}
    // Set the valid value markers of each variable.
    for(ii=0; ii<indeg; ++ii){
        //Aux_Convert2(PossibleLocal[Parents_temp[Code][ii]], local[ii], NumSys);
        Aux_Convert2(PossibleLocal[Parents_temp[Code][indeg-1-ii]], local[ii], L);}
    // Search mapping table.
    New_State=0;
    for(ii=0; ii<LL_sys[InDeg_temp[Code]]; ++ii){
        flag=1;
        for(jj=0;jj<indeg;jj++){
            flag&=local[jj][TempMappingVec[jj]];
            if(0==flag)break;
        }
        if(flag){// This mapping input vector is valid.
            New_State|=(1<<(*ipt_map));
        }
        if(New_State==Old_State){break;}// Still contain all old possible state.
        Aux_NumSysAdder(TempMappingVec, indeg, L);// Next mapping vector.
        ipt_map++;
    }
    free(TempMappingVec);
    free(ValidVariable);
    free(local);
    return New_State;
}

// Clamp nodes. (Fixed to next fixed)//[v]
DNS_Engaged& DNS_Engaged::ClampedVertex(){
    int ii,jj;
    std::set<int> NotExamine, NextLoop;
    std::set<int>::iterator it;
    // Check fixed point and invalid inputting.
    NotExamine=GlobalCandidate;
    for(auto it=NotExamine.begin(); it!=NotExamine.end(); ++it){
        ii=(*it);
        DNS_Engaged::RemoveInvalidInputs(ii);
        jj=InDeg_temp[ii];// The in-degree has been updated.
        if(jj){// Have valid inputs.
            PossibleLocal[ii]=Aux_Map2PossibleValue(Mapping[ii], LL_sys[jj], PossibleLocal[ii], L);
            PossibleValue[ii]=Aux_Num1inInt32(PossibleLocal[ii]);}
        else {// Is a fixed point.
            PossibleLocal[ii]=(1<<Mapping[ii][0]);
            PossibleValue[ii]=1;}
        if(1==PossibleValue[ii]){// Remove initial node.
            StableNode.insert(ii);
            GlobalCandidate.erase(ii);}}
    // Set all nodes for being checked.
    NotExamine=GlobalCandidate;
    // Repeating check. (Espeacailly for multi-value system)
    while(NotExamine.size()){
        NextLoop.clear();
        for(it=NotExamine.begin(); it!=NotExamine.end(); ++it){
            jj=DNS_Engaged::NodeShrunk(*it);
            if(jj<PossibleLocal[*it]){// Mulitple-state's shrinking.
                // Note: for Boolean system, they are directly shrunk as stable node.
                //StableNode.insert(*it);// Only for Boolean
                //GlobalCandidate.erase(*it);// Only for Boolean 
                PossibleLocal[*it]=jj;
                PossibleValue[*it]=Aux_Num1inInt32(jj);
                if(1==PossibleValue[*it]){// Can be fixed!
                    StableNode.insert(*it);
                    GlobalCandidate.erase(*it);
                }
                for(auto xx=(*Children_temp)[*it].begin();xx!=(*Children_temp)[*it].end();++xx){
                    if(PossibleValue[*xx]>1){// Only unstable nodes can be added.
                        NextLoop.insert(*xx);}
                    }
                }
            }
        NotExamine=NextLoop;
    }
    return *this;
}

// Prune edges//[v]
DNS_Engaged& DNS_Engaged::PrunedVertex(){
    int ii, jj, tmp;
    std::set<int> candidate,NextLoop; candidate.clear();
    std::set<int>::iterator it;
    // Stage A: Remove the stable vertexes.
    for(it=StableNode.begin(); it!=StableNode.end(); ++it){
        for(jj=0; jj<InDeg_temp[*it]; ++jj){
            tmp=Parents_temp[*it][jj];
            (*Children_temp)[tmp].erase(*it);
            OtDeg_temp[tmp]--;
        }
        InDeg_temp[*it]=0;
    }
    candidate=GlobalCandidate;
    // Stage B: Remove the non-regulating vertexes. (Also Including merely self-loop out-links)
    while(candidate.size()){
        NextLoop.clear();
        for(it=candidate.begin(); it!=candidate.end(); ++it){
            ii=(*it);
            if(0==OtDeg_temp[ii]){// No output!!
                GlobalCandidate.erase(ii);
                if(StableNode.find(ii)==StableNode.end()){// Priority Stable > Useless, unless repeat counting.
                    UselessNode.insert(ii);
                }
                for(jj=0; jj<InDeg_temp[ii]; ++jj){
                    tmp=Parents_temp[ii][jj];
                    (*Children_temp)[tmp].erase(ii);
                    OtDeg_temp[tmp]--;
                    NextLoop.insert(tmp);// Changed, need next check.
                }
                InDeg_temp[ii]=0;
            }
            else if(1==OtDeg_temp[ii]&&((*Children_temp)[ii].find(ii)!=(*Children_temp)[ii].end())){// Merely its self-pointed node.
                GlobalCandidate.erase(ii);
                if(StableNode.find(ii)==StableNode.end()){// Priority Stable > Useless, unless repeat counting.
                    UselessNode.insert(ii);}
                for(jj=0; jj<InDeg_temp[ii]; ++jj){
                    tmp=Parents_temp[ii][jj];
                    (*Children_temp)[tmp].erase(ii);
                    OtDeg_temp[tmp]--;
                    if(tmp!=ii){// No need check this self-node.
                        NextLoop.insert(tmp);
                    }
                }// Changed, need next check.
                InDeg_temp[ii]=0;
            }
            else ;// Other cases, no operation.
        }
        candidate=NextLoop;
    }
    return *this;
}

// Get once analysis of clamp + prune //[v]
void DNS_Engaged::OnlyScalingPattern(int *OutputVec,int *SUW){
    DNS_Engaged::ClampedVertex();
    DNS_Engaged::PrunedVertex();
    OutputVec[0]=(int)(StableNode.size());// Length of OutputVec >= 3
    OutputVec[1]=(int)(UselessNode.size());
    OutputVec[2]=(int)(GlobalCandidate.size());
    memcpy(SUW, PossibleValue, N*sizeof(int));
}

// Whether nodes are stable (1 or 0) or not?//[v]
std::vector<int> DNS_Engaged::Export2AllNodeState(){
    std::vector<int> res(N);
    memcpy(&res[0], PossibleLocal, N*sizeof(int));
    return res;
}

// Nodes belong to stable, useless, or engaged?//[v]
std::vector<int> DNS_Engaged::Export2AllNodeType(){
    std::vector<int> res(N,666);// 666 for invalid cases.
    for(auto iitt=StableNode.begin(); iitt!=StableNode.end(); ++iitt){
        res[*iitt]=1;}// for stable
    for(auto iitt=UselessNode.begin(); iitt!=UselessNode.end(); ++iitt){
        res[*iitt]=0;}// for useless
    for(auto iitt=GlobalCandidate.begin(); iitt!=GlobalCandidate.end(); ++iitt){
        res[*iitt]=-1;}// for engaged
    return res;
}

// Export residual Network into vec<vec<int>> objects.
DNS_Engaged& DNS_Engaged::Export2ResidualNetwork(
    std::vector<std::vector<int>> &Exp_InEdge, std::vector<std::vector<int>> &Exp_OtEdge, 
    std::vector<std::vector<int>> &BoolFn){
    // Before using this function, "CuttingLeavesVertex" should be executed in advanced.
    DNS_Engaged::CuttingLeavesVertex();// to remove "s-->u" useless edge.
    int n_residual=GlobalCandidate.size();
    Exp_InEdge.clear(); Exp_InEdge.resize(N);
    Exp_OtEdge.clear(); Exp_OtEdge.resize(N);
    BoolFn.clear();     BoolFn.resize(N);
    // Configure the network.
    if(n_residual>0){
        int ii,jj,id;
        for(auto iitt1=GlobalCandidate.begin(); iitt1!=GlobalCandidate.end(); ++iitt1){
            // [NOTE] Due to Clamped -> Pruned -> Cutting Leaves, all remaining nodes are
            // unstable, both having in-edge and out-edge.
            id=(*iitt1);
            // Set in-edges
            Exp_InEdge[id].resize(InDeg_temp[id]);
            for(ii=0; ii<InDeg_temp[id]; ++ii){
                Exp_InEdge[id][ii]=Parents_temp[id][ii];}
            // Set out-edges
            Exp_OtEdge[id].resize(OtDeg_temp[id]); jj=0;
            for(auto iitt2=(*Children_temp)[id].begin(); iitt2!=(*Children_temp)[id].end(); ++iitt2){
                Exp_OtEdge[id][jj++]=(*iitt2);}
            // Set Boolean function.
            BoolFn[id].resize(1<<InDeg_temp[id]);
            std::copy(Mapping[id], Mapping[id] + (1<<InDeg_temp[id]), BoolFn[id].begin());
        }
    }
    return *this;
}

// Cutting all stable leaf-nodes//[v]
DNS_Engaged& DNS_Engaged::CuttingLeavesVertex(){
    int ii, jj, tmp_id;
    std::vector<int> ss1,ss0;
    for(auto it=StableNode.begin(); it!=StableNode.end(); ++it){
        if(OtDeg_temp[*it]>0){
            ss0.clear(); ss1.clear(); ii=(*it);
            // [NOTE] Only for Boolean system. Possible Local, 2=>1, 1=>0
            if(PossibleLocal[ii]>>1){// Stable value is 1 
                for(auto iitt=(*Children_temp)[ii].begin(); iitt!=(*Children_temp)[ii].end(); ++iitt){
                    tmp_id=(*iitt);
                    for(jj=0; jj<InDeg_temp[tmp_id]; jj++){
                        if(ii==Parents_temp[tmp_id][jj]){
                            break;}}
                    //ss1.push_back(jj);// jj-th parent should be removed.
                    ss1.push_back(InDeg_temp[tmp_id]-1-jj);// << --inverse! From high to low!
                    Aux_SimplifyMappingTable(Mapping[tmp_id],InDeg_temp[tmp_id],ss0,ss1);
                    Aux_SimplifyParents(Parents_temp[tmp_id],InDeg_temp[tmp_id],jj);
                    InDeg_temp[tmp_id]--;
                    ss1.clear();}}
            else {// Stable value is 0
                for(auto iitt=(*Children_temp)[ii].begin(); iitt!=(*Children_temp)[ii].end(); ++iitt){
                    tmp_id=(*iitt);
                    for(jj=0; jj<InDeg_temp[tmp_id]; jj++){
                        if(ii==Parents_temp[tmp_id][jj]){
                            break;}}
                    //ss0.push_back(jj);
                    ss0.push_back(InDeg_temp[tmp_id]-1-jj);// << --inverse!
                    Aux_SimplifyMappingTable(Mapping[tmp_id],InDeg_temp[tmp_id],ss0,ss1);
                    Aux_SimplifyParents(Parents_temp[tmp_id],InDeg_temp[tmp_id],jj);
                    InDeg_temp[tmp_id]--;
                    ss0.clear();}}
            OtDeg_temp[ii]=0;
            (*Children_temp)[ii].clear();}}
    return *this;
}


// DNS_CoreDyn's Constructor & Destructor//[V]
DNS_CoreDyn::DNS_CoreDyn(): DNS_Engaged(){
}
DNS_CoreDyn::DNS_CoreDyn(int size, int L, std::vector<std::vector<int>> *ipt1,
    std::vector<std::vector<int>> *ipt2, std::vector<std::vector<int>> *ipt3): 
    DNS_Engaged(size, L, ipt1, ipt2, ipt3){
}
DNS_CoreDyn::~DNS_CoreDyn(){
}

// Decode feed-forward loops & Remove some redundant edges.
DNS_CoreDyn& DNS_CoreDyn::SingleInputCouple(int IDxx){
    int ii, tmp_id, tmp_par, tmp_loc;
    bool isFFL, Inverse;
    std::set<int> tmp_child;
    if(InDeg_temp[IDxx]==1&&Parents_temp[IDxx][0]!=(IDxx)){
        tmp_id=(IDxx);
        tmp_par=Parents_temp[tmp_id][0];
        Inverse=(1==Mapping[tmp_id][0]);// Should inverse the mapping table.
        tmp_child=(*Children_temp)[tmp_id];// Surrogate set, avoid intervenct
        for(auto iipp=tmp_child.begin(); iipp!=tmp_child.end(); iipp++){
            // Check whether tri-nodes form a reducible regulatory pattern.
            isFFL=((*Children_temp)[tmp_par].find(*iipp)!=(*Children_temp)[tmp_par].end());
            if(isFFL){// Yes, a triangle pattern (Forward feedback loop): A -> B -> C & A -> C
                // It also contains the case: A -> B -> A, added or combined as A -> A.
                // Check two upper can be integrated!
                int id1, id2;// 1: top node, 2: intermediate node
                bool found1=false, found2=false;
                for(ii=0; ii<InDeg_temp[*iipp]; ++ii){
                    if(tmp_par==Parents_temp[*iipp][ii]){
                        id1=ii; found1=true;}
                    if(tmp_id==Parents_temp[*iipp][ii]){
                        id2=ii; found2=true;}
                    if(found1&&found2)break;}
                if(found1==false||found2==false)Rcpp::stop("Somthing worng!!\n");
                if(id1==id2){// May have a self-loop node
                    Rcpp::stop("%d ==>(this node:%d)~~~%d\n",
                        tmp_par,tmp_id,InDeg_temp[IDxx]);}
                Aux_TwoUpperNodeCoupling(InDeg_temp[*iipp],id1,id2,Inverse,Mapping[*iipp]);
                Aux_SimplifyParents(Parents_temp[*iipp],InDeg_temp[*iipp],id2);
                // OtDeg_temp[tmp_par]++;// Note that top node originally exist.
                OtDeg_temp[tmp_id]--;
                InDeg_temp[*iipp]--;// B is removed!
                (*Children_temp)[tmp_id].erase(*iipp);
            }
            else {// No, only a straight-like patterns: A -> B -> C
                for(ii=0; ii<InDeg_temp[*iipp]; ++ii){// Replace node [B]
                    if(tmp_id==Parents_temp[*iipp][ii]){
                        tmp_loc=ii;
                        Parents_temp[*iipp][ii]=tmp_par;
                        (*Children_temp)[tmp_par].insert(*iipp);
                        OtDeg_temp[tmp_par]++;
                        break;}}
                // Change some elements in C's mapping table.
                if(Inverse){// Should inverse some bits.
                    int *tmp_bn=(int*)malloc(sizeof(int)*(1<<InDeg_temp[*iipp]));// int tmp_bn[(1<<InDeg_temp[*iipp])];
                    tmp_loc=InDeg_temp[*iipp]-1-tmp_loc;// First element is highest bit (inverse location!).
                    tmp_loc=1<<tmp_loc;
                    for(ii=0; ii<(1<<InDeg_temp[*iipp]); ii++){// Only check half.
                        tmp_bn[ii]=Mapping[*iipp][ii^tmp_loc];}
                    memcpy(Mapping[*iipp],tmp_bn,sizeof(int)*(1<<InDeg_temp[*iipp]));
                    free(tmp_bn);}
                // Handle the intermediate node [B]
                OtDeg_temp[tmp_id]--;
                (*Children_temp)[tmp_id].erase(*iipp);
            }
        }
        if(0==OtDeg_temp[tmp_id]){// Has no child.
            InDeg_temp[IDxx]=0;
            GlobalCandidate.erase(IDxx);// Remove intermediate nodes.
            UselessNode.insert(IDxx);// Record these nodes.
            (*Children_temp)[IDxx].clear();// Actually, Redundant opertion.
            (*Children_temp)[tmp_par].erase(IDxx);// Remove parent -->> current
            OtDeg_temp[tmp_par]--;// Remember this.
        }
    }
    return *this;
}
DNS_CoreDyn& DNS_CoreDyn::DoubleInputCoupleType01(int IDxx, int Successor, int Par_H, int Par_L){// Single couple
    int ii;
    int tmp_child_id=InDeg_temp[Successor];
    int tmp_child_ln=1<<(tmp_child_id);
    int *new_bn_val=(int*)malloc(sizeof(int)*tmp_child_ln);// int new_bn_val[tmp_child_ln];
    int Replace=(1<<tmp_child_id)-1;
    int tmp_match;
    int loc_md, loc_cm;
    int loc_re=666, loc_pl=666, loc_ph=666;
    for(ii=0; ii<tmp_child_id; ++ii){
        if(IDxx==Parents_temp[Successor][ii]){
            loc_re=tmp_child_id-1-ii;}
        else if(Par_L==Parents_temp[Successor][ii]){
            loc_pl=tmp_child_id-1-ii;}
        else if(Par_H==Parents_temp[Successor][ii]){
            loc_ph=tmp_child_id-1-ii;}
        else ;}
    if(loc_pl==666){// Coupled the higher bit
        loc_md=1<<loc_re;
        loc_cm=1<<loc_ph;
        Replace=Replace^loc_md;
        for(ii=0; ii<tmp_child_ln; ++ii){
            tmp_match=Mapping[IDxx][(((ii&loc_cm)>0)<<1)+((ii&loc_md)>0)];
            new_bn_val[ii]=Mapping[Successor][(ii&Replace)+(tmp_match<<loc_re)];
        }
        // Replace the parent.
        Parents_temp[Successor][tmp_child_id-1-loc_re]=Par_L;
        (*Children_temp)[Par_L].insert(Successor);
        OtDeg_temp[Par_L]++;
    }
    else {// Coupled the lower bit
        loc_md=1<<loc_re;
        loc_cm=1<<loc_pl;
        Replace=Replace^loc_md;
        for(ii=0; ii<tmp_child_ln; ++ii){
            tmp_match=Mapping[IDxx][(((ii&loc_md)>0)<<1)+((ii&loc_cm)>0)];
            new_bn_val[ii]=Mapping[Successor][(ii&Replace)+(tmp_match<<loc_re)];
        }
        // Replace the parent.
        Parents_temp[Successor][tmp_child_id-1-loc_re]=Par_H;
        (*Children_temp)[Par_H].insert(Successor);
        OtDeg_temp[Par_H]++;
    }
    memcpy(Mapping[Successor],new_bn_val,tmp_child_ln*sizeof(int));
    OtDeg_temp[IDxx]--;
    (*Children_temp)[IDxx].erase(Successor);
    free(new_bn_val);
    return *this;
}
DNS_CoreDyn& DNS_CoreDyn::DoubleInputCoupleType02(int IDxx, int Successor, int Par_H, int Par_L){// Double couple
    int ii;// int Upper_l=Par_L, Upper_h=Par_H;// Inverse!
    int tmp_child_id=InDeg_temp[Successor];
    int tmp_child_ln=1<<(tmp_child_id);
    int loc_de=666, loc_p1=666, loc_p2=666;
    // Begin analysis: {A+B+X --> C} & {A --> C} & {B --> C} ==> {A+B --> C}
    int *new_bn_val=(int*)malloc(sizeof(int)*(tmp_child_ln>>1));// int new_bn_val[tmp_child_ln>>1];
    int *new_bn_loc=(int*)malloc(sizeof(int)*(tmp_child_ln>>1));// int new_bn_loc[tmp_child_ln>>1];
    for(ii=0; ii<tmp_child_id; ++ii){
        if(IDxx==Parents_temp[Successor][ii]){
            loc_de=tmp_child_id-1-ii;}
        else if(Par_L==Parents_temp[Successor][ii]){
            loc_p1=tmp_child_id-1-ii;}
        else if(Par_H==Parents_temp[Successor][ii]){
            loc_p2=tmp_child_id-1-ii;}
        else ;
    }
    // Judge new location and shrunk truth tables (TTs).
    int Fuller=(1<<tmp_child_id)-1;
    int LowerBitGet=(1<<loc_de)-1;
    int HigherBitGet=Fuller^((1<<loc_de)+LowerBitGet);
    int tmp1=1<<loc_p1, tmp2=1<<loc_p2, tmp_de=1<<loc_de, tmp_match;
    int count=0;
    for(ii=0; ii<tmp_child_ln; ++ii){
        tmp_match=Mapping[IDxx][(((ii&tmp2)>0)<<1)+((ii&tmp1)>0)];
        if(((ii&tmp_de)>0)==tmp_match){// Meet condition in child's TT.
            new_bn_val[count]=Mapping[Successor][ii];
            new_bn_loc[count]=((ii&HigherBitGet)>>1)+(ii&LowerBitGet);// 
            count++;}}
    for(ii=0; ii<(tmp_child_ln>>1); ++ii){
        Mapping[Successor][new_bn_loc[ii]]=new_bn_val[ii];}
    // Adjust the network's topology.
    Aux_SimplifyParents(Parents_temp[Successor],tmp_child_id,tmp_child_id-1-loc_de);
    InDeg_temp[Successor]--;
    OtDeg_temp[IDxx]--;
    (*Children_temp)[IDxx].erase(Successor);
    free(new_bn_val);free(new_bn_loc);
    return *this;
}
DNS_CoreDyn& DNS_CoreDyn::DoubleInputCouple(int IDxx){
    int par_h=Parents_temp[IDxx][0], par_l=Parents_temp[IDxx][1];// NOTE: Inverse order in memory.
    bool coup_h,coup_l;
    std::set<int> tmp_child=(*Children_temp)[IDxx];// Surrogate set, avoid intervenct
    for(auto iipp=tmp_child.begin(); iipp!=tmp_child.end(); iipp++){
        coup_h=((*Children_temp)[par_h].find(*iipp)!=(*Children_temp)[par_h].end());
        coup_l=((*Children_temp)[par_l].find(*iipp)!=(*Children_temp)[par_l].end());
        if(coup_h||coup_l){
            if(coup_h&&coup_l){// {A+B+X --> C} & {A --> C} & {B --> C}
                DoubleInputCoupleType02(IDxx, *iipp, par_h, par_l);
            }
            else {// // {A+B+X --> C} & ( {A --> C} | {B --> C} )
                DoubleInputCoupleType01(IDxx, *iipp, par_h, par_l);
            }
        }
    }
    return *this;
}
DNS_CoreDyn& DNS_CoreDyn::RemoveStableNonZeroOutput2(){
    // Part1: Stable-leaves cutting. (StableNode, {UselessNode} not change)
    DNS_CoreDyn::CuttingLeavesVertex();
    // Part2: Single transmition remove.
    std::set<int> tmp_child, tmp_can=GlobalCandidate;
    for(auto iitt=tmp_can.begin(); iitt!=tmp_can.end(); iitt++){
        if(InDeg_temp[*iitt]==1&&Parents_temp[*iitt][0]!=(*iitt)){
            DNS_CoreDyn::SingleInputCouple(*iitt);}
        else if(2==InDeg_temp[*iitt]){
            if((*Children_temp)[*iitt].find(*iitt)==(*Children_temp)[*iitt].end()){// !NOT contain self-looping cricuit.
                DNS_CoreDyn::DoubleInputCouple(*iitt);}}
        else ;}
    return *this;
}
void DNS_CoreDyn::OnceCoreDynamic(int *OutputVec, int Times){
    int tmpI[3];
    bool GoOn=true;int coun=0;// &&coun<2
    while(GoOn&&(coun<Times)){
        // Rprintf("A loop:\n");
        tmpI[0]=(int)(StableNode.size());
        tmpI[1]=(int)(UselessNode.size());
        tmpI[2]=(int)(GlobalCandidate.size());
        DNS_CoreDyn::ClampedVertex();// Due to triangle-loop merging, some node can be stable.
        DNS_CoreDyn::PrunedVertex();
        DNS_CoreDyn::RemoveStableNonZeroOutput2();
        // for(int ii=0; ii<N; ii++){
        //     Rprintf("(%d,%d,%d)",ii,InDeg_temp[ii],OtDeg_temp[ii]);}Rprintf("\n");
        GoOn=(tmpI[0]!=(int)(StableNode.size()))||
            (tmpI[1]!=(int)(UselessNode.size()))||
            (tmpI[2]!=(int)(GlobalCandidate.size()));
        coun++;
    }
    OutputVec[0]=tmpI[0];
    OutputVec[1]=tmpI[1];
    OutputVec[2]=tmpI[2];
}
// Code is over!
