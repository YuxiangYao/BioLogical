#include "MulVFun.h"

void MulVFun_Adder(std::vector<unsigned short> &Slots, int k, int L){
    int loc=0, Ceiling=L-1;
    while(loc<k){
        if(Slots[loc]<Ceiling){Slots[loc]++; break;}
        else {Slots[loc]=0; loc++;}}
}

void Replaced_Shuffle(std::vector<short> &aVec){
    int idss, ii;
    short vals;
    for(ii=(int)(aVec.size())-1; ii>=0; ii=ii-1){
        idss=(int)(unif_rand()*ii);
        vals=aVec[ii];
        aVec[ii]=aVec[idss];
        aVec[idss]=vals;
    }
}

int VecDotTimes(std::vector<unsigned short> &Vec, std::vector<unsigned int> &LL_system, int lens){
    int ii, sum=0;
    for(ii=0; ii<lens; ++ii){
        sum+=Vec[ii]*LL_system[ii];}
    return sum;
}

int VecDotTimes_pn(std::vector<unsigned short> &vec1, std::vector<int> &vec2, int lens){
    int ii, sum=0;
    for(ii=0; ii<lens; ++ii){
        sum+=vec1[ii]*vec2[ii];}
    return sum;
}

// Basic operations mulvfun
mulvfun& mulvfun::Reset(){
    k=0;
    L=3;
    length=1;
    bias=nullptr;
    ttt=nullptr;
    Type='X';
    argus.resize(0);
    return *this;
}
void mulvfun::aux_RMF(){
    int jj;
    double Rander;
    bool flag;
    for(int ii=0; ii<length; ++ii){
        Rander=unif_rand();
        flag=true;
        while(flag){
            for(jj=0; jj<L; ++jj){
                if(Rander<bias[jj]){
                    ttt[ii]=jj; flag=false;
                    break;}
                else if(jj==(L-1)){
                    ttt[ii]=jj; flag=false;break;}
                else {
                    Rander-=bias[jj];
                }
            }
        }
    }
}

mulvfun& mulvfun::Configuration(double *BiasSet, short *maps, std::vector<int> &pars){
    bias=BiasSet;
    ttt=maps;
    switch(Type){
        case 'R':case 'X'://break;
        case 'C':case 'D':case 'T':case 'S':
            argus=pars;
            break;
        default: Rcpp::stop("Invalid type input.\n");break;}
    return *this;
}

mulvfun& mulvfun::Gen_Rand(){
    mulvfun::aux_RMF();
    return *this;
}

mulvfun& mulvfun::Gen_Cana_Free(Rcpp::IntegerVector &CanaVarID, Rcpp::IntegerVector &EachCanaNum){// std::vector<int> &Config
    int ii, jj, kk, par1=argus[0];// par1: Deep-level: par1<=k
    std::vector<short> k_ID(k), L_ID(L);// short k_ID[k], L_ID[L];
    std::vector<int> nEachVar(par1);
    std::vector<std::vector<short>> CanaIns(k, std::vector<short> (L,-1) );
    std::vector<std::vector<short>> CanaOut(k, std::vector<short> (L,-1) );
    for(ii=0; ii<L; ++ii){// Set slot of tags of discrete system.
        L_ID[ii]=(short)ii;}
    // Check arguments from user-setting.
    if(CanaVarID[0]<=-1){// Not provided by user.
        for(ii=0; ii<k; ++ii){// Set slot of tags of discrete system.
            k_ID[ii]=(short)ii;}
        //std::shuffle(k_ID.begin(), k_ID.end(), mt);
        Replaced_Shuffle(k_ID);
    }
    else {
        for(ii=0; ii<par1; ++ii){
            k_ID[ii]=(short)CanaVarID[ii];}}
    if(EachCanaNum[0]<=0){// Not provided by user.
        for(ii=0; ii<par1; ++ii){// Cana: [1,L]
            nEachVar[ii]=(int)(unif_rand()*L)+1;}}
    else {
        for(ii=0; ii<par1; ++ii){
            nEachVar[ii]=EachCanaNum[ii];}}
    for(ii=0; ii<par1; ++ii){
        // Config canalizing values (in): randomly choose Config[ii][0] from [0,1,2,...,L-1]
        //std::shuffle(L_ID.begin(), L_ID.end(), mt);
        Replaced_Shuffle(L_ID);
        std::copy(L_ID.begin(), L_ID.begin()+nEachVar[ii], CanaIns[ii].begin());
        // Config canalized values (out): randomly map to any one of [0,1,2,...,L-1]
        for(jj=0; jj<nEachVar[ii]; ++jj){// Config[ii][1] ≤ Config[ii][0]
            CanaOut[ii][jj]=(short)(unif_rand()*L);
        }
    }
    // Set input vector [0,0, ...] -> [L-1,L-1, ...]
    std::vector<unsigned short> InputCombined(k,0);
    int Upper=par1-1;
    short tmp_map_to, code_id;
    bool tar1;
    for(ii=0; ii<length; ++ii){// Each input vector.
        tmp_map_to=(short)(unif_rand()*L);// Firstly, set a random number [0,L-1].
        for(jj=Upper; jj>-1; --jj){// Each canalizing/canalized layer.
            code_id=k_ID[jj];
            tar1=true;
            for(kk=0; kk<nEachVar[jj]&&tar1; ++kk){// Each canalizing values.
                if(InputCombined[code_id]==CanaIns[jj][kk]){
                    tmp_map_to=CanaOut[jj][kk];
                    tar1=false;// One input vector, only one "code_id" can meet.
                }
            }
        }
        ttt[ii]=tmp_map_to;
        // Rprintf("{%hu,%hu,%hu}\n",InputCombined[2],InputCombined[1],InputCombined[0]);
        MulVFun_Adder(InputCombined,k,L);
    }
    return *this;
}

mulvfun& mulvfun::Gen_Cana_Config(Rcpp::IntegerVector &CanaVarID, std::vector<int> &EachCanaNum, 
    std::vector<std::vector<short>> &CanaIns, std::vector<std::vector<short>> &CanaOut){
    // par1: Deep-level: par1<=k
    int ii, jj, kk;int par1=argus[0];
    // Set input vector [0,0, ...] -> [L-1,L-1, ...]
    // short k_ID[k];// L_ID[L];
    std::vector<short> k_ID(k);
    std::vector<unsigned short> InputCombined(k,0);
    int Upper=argus[0]-1;
    short tmp_map_to, code_id;
    bool tar1;
    // Not provided by user.
    if(CanaVarID[0]==NA_INTEGER){
        for(ii=0; ii<k; ++ii){// Set slot of tags of discrete system.
            k_ID[ii]=(short)ii;}
        //std::shuffle(k_ID.begin(), k_ID.end(), mt);
        Replaced_Shuffle(k_ID);
    }
    else {
        for(ii=0; ii<par1; ++ii){
            k_ID[ii]=(short)CanaVarID[ii];}}
    // EachCanaNum can be inferred by &CanaIns + &CanaOut in upstream code.
    // Fill the mapping table.
    for(ii=0; ii<length; ++ii){// Each input vector.
        tmp_map_to=(short)(unif_rand()*L);// Firstly, set a random number [0,L-1].
        //tmp_map_to=666;
        for(jj=Upper; jj>-1; --jj){// Each canalizing/canalized layer.
            code_id=k_ID[jj];
            tar1=true;
            for(kk=0; kk<EachCanaNum[jj]&&tar1; ++kk){// Each canalizing values.
                if(InputCombined[code_id]==CanaIns[jj][kk]){
                    tmp_map_to=CanaOut[jj][kk];
                    tar1=false;
                }
            }
        }
        ttt[ii]=tmp_map_to;
        MulVFun_Adder(InputCombined,k,L);
    }
    return *this;
}

// Calculate sensitivity of multi-valued functions.
int mulvfun::OneVariableActivity(std::vector<unsigned short> &Vec, int &VarID, std::vector<unsigned int> &LL_system){
    int ii, different=0, Pointed=Vec[VarID], Multiper=LL_system[VarID];
    // Remove the value.
    Vec[VarID]=0;
    int RemainingSum=VecDotTimes(Vec,LL_system,k);
    int id_std=ttt[RemainingSum+Pointed*Multiper];
    // Its forward data.
    for(ii=0; ii<Pointed; ++ii){
        different+=(id_std!=ttt[RemainingSum]);
        RemainingSum+=Multiper;}
    RemainingSum+=Multiper;
    // Its subsequent data.
    for(ii=Pointed+1; ii<L; ++ii){
        different+=(id_std!=ttt[RemainingSum]);
        RemainingSum+=Multiper;}
    Vec[VarID]=Pointed;// Recover values.
    return different;
}
double mulvfun::Sensitivity(){
    int ii,jj;
    std::vector<unsigned short> input_vec(k,0);
    std::vector<unsigned int> LL_sys(k); LL_sys[0]=1;
    for(ii=0; ii<k-1; ++ii){
        LL_sys[ii+1]=LL_sys[ii]*L;}
    unsigned int sum=0;
    for(ii=0; ii<length; ++ii){
        for(jj=0; jj<k; ++jj){
            sum+=OneVariableActivity(input_vec, jj, LL_sys);}
        MulVFun_Adder(input_vec,k,L);}
    return (double)sum/(double)((L-1)*length);
}

// Return a multi_valued canalized function's hierarchical structure.
// Any one element is TRUE can keep search solves of table.
bool AnyStateTRUE(std::vector<bool> &vec){
    bool tag=true; 
    int n=vec.size();
    for(int ii=0; ii<n&&(tag); ++ii){
        tag=!vec[ii];}
    return (!tag);
}
// Check one canalized value.
std::vector<std::vector<short>> mulvfun::OnceCanalizedChecker(int PointedID, std::vector<bool> &UnCana){
    int ii, code;
    std::vector<unsigned short> input_vec(k,0);
    short *MapTab=ttt;
    std::vector<short> Tag(L,-1);// Record canalized values.
    std::vector<bool> ForAll(L,true);// If this value is canalizing?
    bool StillCheck=true;
    for(ii=0; ii<length&&StillCheck; ++ii){
        if(UnCana[ii]){// Not be canalized!
            code=input_vec[PointedID];
            if(Tag[code]<0){
                Tag[code]=(*MapTab);// Record the first canalized valued.
            }
            else {
                if(ForAll[code]){
                    ForAll[code]=ForAll[code]&&(Tag[code]==(*MapTab));// Keep same value?
                    StillCheck=AnyStateTRUE(ForAll);
                }
            }
        }
        MulVFun_Adder(input_vec,k,L);
        MapTab++;// Move pointer forward. 
    }
    if(StillCheck){// Record information & update UnCana-slot
        std::vector<std::vector<short>> Results(2);
        for(ii=0; ii<L; ++ii){
            if(ForAll[ii]){
                Results[0].push_back(ii);
                Results[1].push_back(Tag[ii]);
            }
        }
        // Update UnCana-slot.
        std::vector<unsigned short> input_vec2(k,0);
        //MapTab=ttt;
        for(ii=0; ii<length; ++ii){
            if(UnCana[ii]){
                code=input_vec2[PointedID];
                UnCana[ii]=!ForAll[code];// Canalized/true ==> false for next loop check.
            }
            MulVFun_Adder(input_vec2,k,L);
            //MapTab++;
        }
        //Transfer=Results;
        return Results;//std::move(Results);
    }
    else {
        std::vector<std::vector<short>> Results(2, std::vector<short>(1,-1));
        //Transfer=Results;
        return Results;//std::move(Results);
    }
}
// Main function of checking caanalizing/canalized correlation.
std::vector<std::vector<std::vector<short>>> mulvfun::NestCana(std::vector<int> &CanalizingIDs, bool Quick){
    std::vector<bool> checker(length,true);
    std::vector<bool> Canalized(k,true);
    int ii;
    int unstable=k;
    if(Quick){// Only need check outermost layer 
        unstable=1;}
    bool FullCana=true, inner;
    std::vector<std::vector<short>> Transfer;
    // Loop through each variable. 
    CanalizingIDs.clear();
    std::vector<std::vector<std::vector<short>>> Final_Cana;
    while(unstable>0&&FullCana){
        inner=true;
        FullCana=false;
        for(ii=0; (ii<k)&&inner; ++ii){
            if(Canalized[ii]){// Not been canalized.
                Transfer=OnceCanalizedChecker(ii, checker);//Rprintf("insert-1,(%d)\n",ii);
                if(0<=Transfer[0][0]){// Can be canalized!
                    CanalizingIDs.push_back(ii);
                    Final_Cana.push_back(Transfer);//Rprintf("insert-2,(%d)\n",ii);
                    unstable--;
                    Canalized[ii]=false;
                    inner=false;
                    FullCana=((int)Transfer[0].size())<L;
                }
            }
        }
    }
    return Final_Cana;// return std::move(Final_Cana);
}

// Multi-valued linear threshold function
mulvfun& mulvfun::Gen_MulVF_Thre(){
    int ii;
    std::vector<int> Adj(k);
    int MAXer=(k<<1)+1;
    int Upper_boundary=0, Down_boundary=0;
    for(ii=0; ii<k; ++ii){
        // Adj[ii]=(int)(MAXer*unif_rand())-k;// Boundary of [-MAX, +MAX]
        Adj[ii]=(int)(MAXer*unif_rand())+1;// Only limited to [-2k,-1] & [+1,+2k]
        if(unif_rand()<0.5){
            Adj[ii]=-Adj[ii];}
        if(Adj[ii]>0){// Positive
            Upper_boundary+=Adj[ii]*(L-1);}
        else if(Adj[ii]<0){// Negative
            Down_boundary+=Adj[ii]*(L-1);} }
    double delta=(double)(Upper_boundary-Down_boundary)/(double)(L);
    int tmp;
    std::vector<unsigned short> InputCombined(k,0);
    for(ii=0; ii<length; ++ii){
        tmp=VecDotTimes_pn(InputCombined, Adj, k)-Down_boundary;
        ttt[ii]=(short)(tmp/delta);// Localed in which bin?
        if(ttt[ii]==L){// Sometime can exactly in +MAX.
            ttt[ii]=L-1;}
        MulVFun_Adder(InputCombined, k, L);}
    return *this;
}
Rcpp::List mulvfun::is_MulVF_Thre(){
    int ii, jj, tmp_id;
    std::string TmpName;
    z3::context MapThreshold;
    z3::solver solve_it(MapThreshold);
    std::vector <z3::expr> Thetas, Weights;
    Rcpp::List Res;
    for(ii=0; ii<=L; ++ii){// Set all (L+1) thresholds 
        TmpName.clear();
        TmpName.append("yuzhi_").append(std::to_string(ii));// yuzhi_ii
        Thetas.push_back(MapThreshold.int_const(TmpName.c_str()));}
    for(ii=0; ii<k; ++ii){// Set all k weights
        TmpName.clear();
        TmpName.append("k_").append(std::to_string(ii));// yuzhi_ii
        Weights.push_back(MapThreshold.int_const(TmpName.c_str()));}
    std::vector<unsigned short> InputCombined(k,0);
    for(ii=0; ii<length; ++ii){        
        z3::expr express_tmp=((int)(InputCombined[0]))*Weights[0];//MapThreshold.int_const("0");
        tmp_id=(int)(ttt[ii]);
        for(jj=1; jj<k; ++jj){
            // This is key to judge: \Theta(\sum_j a_{ij}x_j) belongs to which interval
            express_tmp=express_tmp+((int)(InputCombined[jj]))*Weights[jj];}
        MulVFun_Adder(InputCombined, k, L);
        solve_it.add(Thetas[tmp_id]<express_tmp);// Threshold set: 
        solve_it.add(express_tmp<=Thetas[tmp_id+1]);// lower <= ii < upper
    }
    for(ii=0; ii<L; ++ii){
        solve_it.add(Thetas[ii]<Thetas[ii+1]);
    }
    if(solve_it.check()==z3::sat){// Sum zero set as 1.
        z3::model Models=solve_it.get_model();
        Res.push_back(1);
        Rcpp::IntegerVector wight(k);
        Rcpp::IntegerVector theta(L-1);
        for(ii=0; ii<k; ++ii){
            wight[ii]=Models.eval(Weights[ii]).get_numeral_int();}
        for(ii=0; ii<L-1; ++ii){
            theta[ii]=Models.eval(Thetas[ii+1]).get_numeral_int();}
        Res.push_back(wight);
        Res.push_back(theta);}
    else {
        Res.push_back(-1);}
    return Res;
}

// Multi-valued dominant function
short mulvfun::FindDomiantComponent(std::vector<unsigned short> &vec1,
    std::vector<int> &wights, std::vector<double> &theta){
    short Max_ID;
    std::vector<double> Candi(L);
    int tmp;
    for(int ii=0; ii<L; ++ii){
        Candi[ii]=theta[ii];}
    for(int ii=0; ii<k; ++ii){
        tmp=vec1[ii];// This bit is which value?
        Candi[tmp]+=wights[ii];}// Corresponding slot add the weights.
    // Check which ID of L-level is okay.
    Max_ID=0;
    double MAXer=Candi[0];
    for(int ii=1; ii<L; ++ii){
        if(MAXer<Candi[ii]){// Set this new ID!!
            MAXer=Candi[ii];
            Max_ID=ii;}}
    return Max_ID;
}   
mulvfun& mulvfun::Gen_MulVF_Domi(){
    // Boundary of [-MAX, +MAX]
    int ii;
    std::vector<int> Adj(k);
    // int MAXer=(argus[0]<<1)+1;
    for(ii=0; ii<k; ++ii){
        //Adj[ii]=(int)(MAXer*unif_rand())-argus[0];
        Adj[ii]=(int)(k*unif_rand())+1;// Only limited to [-k,-1] & [+1,+k]
        if(unif_rand()<0.5){
            Adj[ii]=-Adj[ii];}
    }
    std::vector<double> l_delta(L);
    for(ii=0; ii<L; ++ii){
        l_delta[ii]=0.1*unif_rand();}
    std::vector<unsigned short> InputCombined(k,0);
    for(ii=0; ii<length; ++ii){
        ttt[ii]=FindDomiantComponent(InputCombined, Adj, l_delta);
        MulVFun_Adder(InputCombined, k, L);}
    return *this;
}
Rcpp::List mulvfun::is_MulVF_Domi(){
    int ii, jj, mm;
    std::string TmpName;
    z3::context MapThreshold;
    z3::solver solve_it(MapThreshold);
    std::vector <z3::expr> Weights,Thetas;
    Rcpp::List Res;
    // Set all weights.
    for(ii=0; ii<k; ++ii){// Set all k weights 
        TmpName.clear();
        TmpName.append("k_").append(std::to_string(ii));// quanzhong_ii
        Weights.push_back(MapThreshold.int_const(TmpName.c_str()));
    }
    for(ii=0; ii<L; ++ii){// Set all k weights 
        TmpName.clear();
        TmpName.append("t_").append(std::to_string(ii));// yuzhi_ii
        Thetas.push_back(MapThreshold.int_const(TmpName.c_str()));
    }
    std::vector<unsigned short> InputCombined(k,0);
    for(ii=0; ii<length; ++ii){
        // New setting.
        std::vector <z3::expr> Accum;
        for(jj=0; jj<L; ++jj){
            Accum.push_back(Thetas[jj]);}
        // Check
        for(jj=0; jj<k; ++jj){
            mm=InputCombined[jj];
            Accum[mm]=Accum[mm]+Weights[jj];}
        // Set MAX!
        mm=(int)(ttt[ii]);
        for(jj=0; jj<mm; ++jj){
            solve_it.add(Accum[mm]>Accum[jj]);}
        for(jj=mm+1; jj<L; ++jj){
            solve_it.add(Accum[mm]>Accum[jj]);}
        MulVFun_Adder(InputCombined, k, L);
    }
    if(solve_it.check()==z3::sat){// Sum zero set as 1.
        z3::model Models=solve_it.get_model();
        Res.push_back(1);
        Rcpp::IntegerVector wight(k);
        Rcpp::IntegerVector theta(L);
        for(ii=0; ii<k; ++ii){
            wight[ii]=Models.eval(Weights[ii]).get_numeral_int();}
        for(ii=0; ii<L; ++ii){
            theta[ii]=Models.eval(Thetas[ii]).get_numeral_int();}
        Res.push_back(wight);
        Res.push_back(theta);}
    else {
        Res.push_back(-1);}
    return Res;
}

// C++ Prototype function of PolynomialFunction(), all parameters are integer
Rcpp::List MulVFun_PolynomialFunction(Rcpp::IntegerMatrix &VariableMat, Rcpp::IntegerVector &maptab, int k, int L_sys){
    Rcpp::List res;
    int ii,jj,mm;
    int nVar=(int) VariableMat.ncol();// { (K,1), or (K,1)+(K,2), or (K,1)+(K,2)+(K,3), or ..., or 2^K-1 }
    int nLen=(int) VariableMat.nrow();// L^K length
    z3::context TT_T;
    std::vector<z3::expr> inputs;
    std::vector<z3::expr> thetas;
    // z3::expr exprss=theta;
    z3::solver solve_it(TT_T);
    for(ii=0; ii<nVar; ++ii){// Set variable number of VAR slots.
        inputs.push_back(TT_T.int_const(("a_" + std::to_string(ii)).c_str()));}
    for(ii=0; ii<L_sys; ++ii){
        thetas.push_back(TT_T.int_const(("b_" + std::to_string(ii)).c_str()));}
    
    // std::vector<unsigned short> InputCombined(k,0);
    for(ii=0; ii<nLen; ++ii){
        // New setting.
        std::vector <z3::expr> Accum;
        for(jj=0; jj<L_sys; ++jj){
            Accum.push_back(thetas[jj]);}
        // Check
        for(jj=0; jj<nVar; ++jj){
            mm=VariableMat(ii,jj);
            Accum[mm]=Accum[mm]+inputs[jj];}
        // Set MAX! Pointed is the domianted one!
        mm=maptab[ii];
        for(jj=0; jj<mm; ++jj){
            solve_it.add(Accum[mm]>Accum[jj]);}
        for(jj=mm+1; jj<L_sys; ++jj){
            solve_it.add(Accum[mm]>Accum[jj]);}
    }
    if(solve_it.check()==z3::sat){// Find out appropriate weights.
        res.push_back(true);// Return the logical tag.
        z3::model Models=solve_it.get_model();
        Rcpp::IntegerVector Weights(nVar,-666);
        Rcpp::IntegerVector Thetass(L_sys,-888);
        // Insert the variables' weights.
        for(ii=0; ii<nVar; ++ii){
            Weights[ii]=Models.eval(inputs[ii]).get_numeral_int();}
        res.push_back(Weights);
        // Return the threshold value.
        for(ii=0; ii<L_sys; ++ii){
            Thetass[ii]=Models.eval(thetas[ii]).get_numeral_int();}
        res.push_back(Thetass);
    }
    else {// Failure 
        res.push_back(false);}
    return (res);
}

// Multi-valued broadly defined-signed function
Rcpp::List mulvfun::is_MulVF_Sign(){
    int ii, jj, kk, mm, res=-1;
    Rcpp::List Res;
    std::string VarName, TmpName;
    z3::context MapCotext;// z3's context
    z3::solver solve_it(MapCotext);// Solver of mapping context.
    char letters[10]={'a','b','c','d','e','f','g','h','i','j'};// MAX is ten!
    std::vector <std::vector<z3::expr> > wights;// The k*L*L weight tensor.
    std::vector<z3::expr> thetas;// L-level benchmarks (threshold)
    wights.resize(k);
    for(ii=0; ii<k; ii++){// Each subsequence has L*L A->A, A->B, ...
        for(jj=0; jj<L; ++jj){
            TmpName.clear();
            TmpName.append("X").append(std::to_string(ii)).append("_");
            TmpName.push_back(letters[jj]); 
            TmpName.append("_to_");
            for(kk=0; kk<L; ++kk){
                VarName.clear();
                VarName.append(TmpName).push_back(letters[kk]);
                wights[ii].push_back(MapCotext.int_const(VarName.c_str()));}}}
    // Build basic solver.
    for(ii=0; ii<L; ++ii){
        thetas.push_back(MapCotext.int_const(("v_" + std::to_string(ii)).c_str()));}
    // Set adder.
    std::vector<unsigned short> InputCombined(k,0);
    for(ii=0; ii<length; ++ii){
        // New setting.
        std::vector <z3::expr> Accum;
        for(jj=0; jj<L; ++jj){
            Accum.push_back(thetas[jj]);}
        // Check them.
        for(jj=0; jj<k; ++jj){
            mm=((int)InputCombined[jj])*L;
            for(kk=0; kk<L; ++kk, ++mm){
                Accum[kk]=Accum[kk]+wights[jj][mm];
            }
        }
        // Set MAX! Pointed is the domianted one!
        mm=(int)ttt[ii];
        for(jj=0; jj<mm; ++jj){
            solve_it.add(Accum[mm]>Accum[jj]);}
        for(jj=mm+1; jj<L; ++jj){
            solve_it.add(Accum[mm]>Accum[jj]);}
        // Auto-AddOne.
        MulVFun_Adder(InputCombined, k, L);
    }  
    if(solve_it.check()==z3::sat){// Sum zero set as 1.
        res=1;
        z3::model Models=solve_it.get_model();
        // Rcpp::IntegerMatrix Wight(k,L*L);
        Rcpp::IntegerMatrix Wight(k*L,L);
        mm=0;
        for(ii=0; ii<k; ++ii){
            int ww=0;
            for(jj=0; jj<L; ++jj){
                for(kk=0; kk<L; ++kk){
                    Wight(mm,kk)=Models.eval(wights[ii][ww]).get_numeral_int();
                    ww++;}
                mm++;
            }
        }
        // for(ii=0; ii<k; ++ii){
        //     for(jj=0; jj<(L*L); ++jj){
        //         Wight(ii,jj)=Models.eval(wights[ii][jj]).get_numeral_int();}}
        Rcpp::IntegerVector Theta(L);
        for(ii=0; ii<L; ++ii){
            Theta[ii]=Models.eval(thetas[ii]).get_numeral_int();}
        // Return valid values of weight-matrix and theta-vector.
        Res.push_back(res);
        Res.push_back(Wight);
        Res.push_back(Theta);}
    else if(solve_it.check()==z3::unknown){// To many to slove it.
        res=0;
        Res.push_back(res);}
    else {
        Res.push_back(res);}
    return Res;
}


// Obtain the means of int vector.
double MeansOfSet(std::vector<int> &xx){// xx.size() always >=1 
    double sum=0;
    if(1==xx.size()){
        sum=(*(xx.begin()));
    }
    else {
        for(auto iitt=xx.begin(); iitt!=xx.end(); ++iitt){
            sum+=(*iitt);
        }
        sum=sum/(double)(xx.size());
    }
    return sum;
}

std::vector<int> MulV2Bool(std::vector<int> &mulvfun, int k, int L, std::vector<int> &threshold){
    int ii, jj, tmp, lens=(1<<k), lengths=(int)pow(L,k);
    std::vector<unsigned short> InputVecs(k,0);
    std::vector<int> NewMapTab(lens);
    std::vector<std::vector<int>> candis(lens);
    for(ii=0; ii<lengths; ++ii){
        tmp=0;
        for(jj=0; jj<k; ++jj){
            tmp+=((InputVecs[jj]>=threshold[jj])<<jj);
        }
        candis[tmp].push_back(mulvfun[ii]);
        MulVFun_Adder(InputVecs, k, L);
    }
    for(ii=0; ii<lens; ++ii){
        NewMapTab[ii]=(int)(threshold[k]<=MeansOfSet(candis[ii]));
    }
    return NewMapTab;
}

std::vector<int> Bool2MulV(std::vector<int> &boolfun, int k, int L, std::vector<int> &threshold){
    int ii, jj, tmp, lens=(int)pow(L,k);// Randomly assign multi2one cases.
    int flag=threshold[k], remain=L-flag;
    std::vector<unsigned short> InputVecs(k,0);
    std::vector<int> NewMapTab(lens);
    for(ii=0; ii<lens; ++ii){
        tmp=0;
        for(jj=0; jj<k; ++jj){
            tmp+=((InputVecs[jj]>=threshold[jj])<<jj);
        }
        if(boolfun[tmp]>0){// Original Booelan vectors Mapto 1.
            NewMapTab[ii]=flag+(int)(remain*unif_rand());
        }
        else {// Original Booelan vectors Mapto 0.
            NewMapTab[ii]=(int)(flag*unif_rand());
        }
        MulVFun_Adder(InputVecs, k, L);
    }
    return NewMapTab;
}
// Code is over.