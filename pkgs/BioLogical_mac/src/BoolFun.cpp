#include "BoolFun.h"
// Basic operations
boolfun& boolfun::Reset(){
    k=0;
    length=1;
    bias=0.5;
    ttt=nullptr;
    Type='X';
    argus[0]=argus[1]=argus[2]=argus[3]=-1;
    return *this;
}
void boolfun::aux_RBF(double p){
    for(int ii=0; ii<length; ++ii){
        ttt[ii]=(unif_rand()<p);
    }
}
double boolfun::NumSize1(){
    int sum=0;
    for(int ii=0; ii<length; ++ii){
        sum+=ttt[ii];
    }
    return (double)(sum);
}
void boolfun::PrintMappingTable(){
    int ii,jj,index;
    for(ii=0; ii<length; ++ii){
        Rcpp::Rcout << '[';
        index=1<<(k-1);
        for(jj=0; jj<k; jj++, index>>=1){
            Rcpp::Rcout << ((ii&index)>0);
        }
        Rcpp::Rcout << "] ["<< ttt[ii] << ']' << std::endl;
    }
}
boolfun& boolfun::Configuration(char BF_Type, int indeg, double Bias, bool *maps, int par1, int par2, int par3, int par4){
    // Old version: Configuration(char BF_Type, int indeg, double bias, bool *maps, ...)
    Type=BF_Type;
    bias=Bias;
    k=indeg;
    length=1<<indeg;
    ttt=maps;
    switch(Type){
        case 'R':break;
        case 'C':case 'P':case 'M':case 'D':case 'T':
            argus[0]=par1; argus[1]=par2; argus[2]=par3; argus[3]=par4;
            break;
        case 'X':break;
        default: Rcpp::stop("Invalid type input.\n");break;
    }
    return *this;
}
// Calculate sensitivity of BFs.
double boolfun::Sensitivity(){
    // [NOTE] Ergodic checking only half part, should be doubled.
    bool *ipt1,*ipt2;
    int ii,jj,kk,boundary_l=1,boundary_u=length>>1,sumer=0;
    for(ii=0; ii<k; ++ii){
        ipt1=ttt;
        ipt2=ttt+boundary_l;
        for(jj=0; jj<boundary_u; ++jj){
            for(kk=0; kk<boundary_l; ++kk){
                sumer+=(ipt1[kk])!=(ipt2[kk]);
            }
            ipt1+=(boundary_l<<1);
            ipt2+=(boundary_l<<1);
        }
        boundary_l=boundary_l<<1;
        boundary_u=boundary_u>>1;
    }
    return (sumer<<1)/(double)(length);
}
void boolfun::aux_ShuffleTTT(){// Random exchange mapping table's values.
    bool exchanger;
    int ii,id1,id2,Less=length-1;
    for(ii=0; ii<(length<<1); ++ii){// Should exchange "0.5*length" times.
        id1=(int)(unif_rand()*length);
        id2=(int)(unif_rand()*Less)+1;
        id2=(id1+id2)%length;// Get the loop.
        exchanger=ttt[id1];
        ttt[id1]=ttt[id2];
        ttt[id2]=exchanger;
    }
}
double boolfun::Energy(){// Calculate function energy.
    int ii,jj,energy=0;// spin[length],bits[k];
    std::vector<int> spin(length), bits(k);
    for(ii=0;ii<k;ii++)bits[ii]=1<<ii;// Set bit-xor slot.
    for(ii=0;ii<length;ii++)spin[ii]=(ttt[ii]<<1)-1;
    for(ii=0;ii<length;ii++){
        for(jj=0;jj<k;jj++){
            energy+=spin[ii]*spin[ii^bits[jj]];
        }
    }
    return 0.5*((double)energy);
}
double boolfun::RelativeEnergy(){
    double eng1,eng2,ReNe=0;
    if('R'==Type);
    else {// Only analyze ordered Boolean functions
        std::vector<bool> tmp(length);// bool tmp[length];
        std::copy(ttt, ttt+length, tmp.begin());
        eng1=boolfun::Energy();
        if(fabs(eng1)>1e-6){
            boolfun::aux_ShuffleTTT();
            eng2=(boolfun::Energy());
            std::copy(tmp.begin(), tmp.begin()+length, ttt);
            ReNe=(fabs(eng1)-fabs(eng2))/fabs(eng1);
        }
        else ReNe=0;
    }
    return ReNe;
}
boolfun& boolfun::Gen_BF(){
    switch(Type){
        case 'R': boolfun::Gen_Rand();break;
        case 'C': boolfun::Gen_Cana();break;
        case 'M': boolfun::Gen_Mone();break;
        case 'P': boolfun::Gen_Post();break;
        case 'T': boolfun::Gen_Spin();break;
        case 'D': boolfun::Gen_Domi();break;
        case 'X':break;
        default:Rcpp::stop("Error BF type.\n");break;
    }
    return *this;
}
bool boolfun::is_PointedType(char types,bool Showit){
    bool Yes_or_No=false;
    switch(types){
        case 'C': Yes_or_No=boolfun::is_Cana();break;
        case 'P': Yes_or_No=boolfun::is_Post();break;
        case 'M': Yes_or_No=boolfun::is_Mone();break;
        case 'S': Yes_or_No=boolfun::is_Sign();break;
        case 'N': Yes_or_No=boolfun::is_Spin();break;
        case 'D': Yes_or_No=boolfun::is_Domi();break;
        case 'T': Yes_or_No=boolfun::is_Thre(false, Showit);break;
        case 'E': Yes_or_No=boolfun::is_Effc();break;
        default:Rcpp::stop("Illegal type!\n");break;
    }
    return Yes_or_No;
}

boolfun& boolfun::Gen_Rand(){
    boolfun::aux_RBF(bias);
    return *this;
}
boolfun& boolfun::Gen_Cana(){
    int ii,jj,code;// int CanaIN[k],CanaOT[k],VarsID[k];
    std::vector<int> CanaIN(k), CanaOT(k), VarsID(k);
    for(ii=0; ii<k; ii++){
        VarsID[ii]=ii;}
    // std::shuffle(VarsID.begin(),VarsID.end(),mt);// For all k, but only choose former layers.
    int idss, vals;
    for(ii=k-1; ii>=0; ii=ii-1){
        idss=(int)(unif_rand()*ii);
        vals=VarsID[ii];
        VarsID[ii]=VarsID[idss];
        VarsID[idss]=vals;}
    if(0==argus[1]){
        for(ii=0; ii<k; ii++){
            CanaIN[ii]=0;
            CanaOT[ii]=0;
        }
        argus[2]=argus[3]=0;
    } else if(1==argus[1]){
        argus[2]=0;
        for(ii=0; ii<k; ii++){
            CanaIN[ii]=1;
            CanaOT[ii]=1;
            argus[2]+=(1<<ii);
        }
        argus[3]=argus[2];
    } else {// Random ~~
        argus[2]=0;argus[3]=0;
        for(ii=0; ii<k; ii++){
            CanaIN[ii]=unif_rand()>0.5;
            CanaOT[ii]=unif_rand()>0.5;
            argus[2]+=(CanaIN[ii]<<ii);
            argus[3]+=(CanaOT[ii]<<ii);
        }
    }
    boolfun::aux_RBF(0.5);// Random set mapping tables.
    for(jj=0; jj<argus[0]; jj++){
        code=VarsID[jj];
        for(ii=0; ii<length; ii++){
            if(CanaIN[jj]==((ii>>code)&1)){
                ttt[ii]=CanaOT[jj];
            }
        }
    }
    return *this;
}
boolfun& boolfun::Gen_Mone(){
    //  Set some random seed vectors, and set their successor vectors as candidate.
    //  Generate an increasing patterns, if (0110010) as one candidate map to one, then all successors should be configurated as one, like (1110010)=>1, (0111010)=>1, ..., (1111111)=>1 (totally 2^4-1=15 successors). */
    int ii,jj,tmp,inverse;// n_seed[k];
    std::vector<int> n_seed(k);
    for(ii=0; ii<argus[1]; ++ii){// Candidate: 1 ~ length-1, can be repeated.
        n_seed[ii]=(int)((length-1)*unif_rand())+1;
    }
    for(ii=0; ii<length; ++ii){
        ttt[ii]=false;
    }
    for(ii=0;ii<argus[1]; ++ii){
        tmp=n_seed[ii];
        ttt[tmp]=true;
        for(jj=tmp+1; jj<length; ++jj){// Forward analysis.
            if((jj&tmp)==tmp){
                ttt[jj]=true;
            }
        }
    }
    if(0==argus[0]||1==argus[0]){// Have been set.
        inverse=argus[0];
    } else {// Randomly set one.
        inverse=unif_rand()>0.5;
    }
    if(0==inverse){// Inverse the relation.
        for(ii=0; ii<length; ++ii){
            ttt[ii]=!ttt[ii];
        }
    }
    return *this;
}
boolfun& boolfun::Gen_Post(){
    // [NOTE] Only consider first element as seed, successor must be intersect each other.
    int ii,jj,tmp=0,flag, intersect, Selected; 
    std::vector<int> Slots(length);// int Slots[length];
    if(argus[1]>0){// Selected nodes can be repeated?
        while(0==tmp){// Get first element, intersecting seed, not zero
            tmp=(int)(length*unif_rand());
        }
        Selected=1;
        Slots[0]=tmp;
        for(ii=0; ii<(length>>1); ii++){// Get 2^k/2 attempts, regardless of whether code is suitable.
            tmp=(int)(length*unif_rand());
            flag=1;intersect=tmp;// Get seed.
            for(jj=0;jj<Selected;jj++){
                intersect&=Slots[jj];
                if(0==intersect){
                    flag=0;break;
                }
            }
            if(flag){// Find a suitable slot.
                Slots[Selected++]=tmp;
            }
        }
    } else {// Selected nodes can be non-repeated!!
        std::vector<bool> labels(length,false);
        while(0==tmp){// Get first element, intersecting seed, not zero
            tmp=(int)(length*unif_rand());
        }
        Selected=1;
        labels[tmp]=true;
        Slots[0]=tmp;
        for(ii=0; ii<(length>>1); ii++){// Get 2^k/2 attempts, regardless of whether code is suitable.
            do{tmp=(int)(length*unif_rand());}while(labels[tmp]);// Has been selected.
            flag=1;intersect=tmp;// Get seed.
            for(jj=0;jj<Selected;jj++){
                intersect&=Slots[jj];
                if(0==intersect){
                    flag=0;break;
                }
            }
            if(flag){// Find a suitable slot.
                labels[tmp]=true;
                Slots[Selected++]=tmp;
            }
        }
    }
    // Need to judge {A^2} or {a^2}.
    if(argus[0]==1||argus[0]==0){// {A^2} or {a^2} class.
        flag=argus[0];
    } else {// Not set Post_Type.
        flag=unif_rand()>0.5;
    }
    if(1==flag){// {A^2} subclass
        for(ii=0; ii<length; ii++){ttt[ii]=false;}
        for(ii=0; ii<Selected; ii++){// Set "one" position.
            ttt[Slots[ii]]=true;}}
    else {// {a^2} subclass
        tmp=length-1;
        for(ii=0; ii<length; ii++){ttt[ii]=true;}
        for(ii=0; ii<Selected; ii++){
            ttt[(~Slots[ii])&tmp]=false;}}// Should get inverse bits.
    return *this;
}
boolfun& boolfun::Gen_Spin(){
    /* [NOTE] Bias p can regulate positive/nagetive inputs. */
    int ii,jj,sum=0; std::vector<int> pn_in(k);// int pn_in[k];
    bool rand_fill;
    for(ii=0;ii<k;ii++){// Set +/- inputs. Here is the core of threshold Boolean function.
        pn_in[ii]=((unif_rand()<bias)<<1)-1;
        sum+=pn_in[ii];}// Balance the positive/nagetive inputs.
    if(0==argus[0]||1==argus[0]){
        rand_fill=(argus[0]>0);}
    else {
        if(sum!=0){// Balance the positive and negative values.
            rand_fill=(sum<0);}
        else rand_fill=(unif_rand()>0.5);}
    for(ii=0;ii<length;ii++){
        sum=0;
        for(jj=0;jj<k;jj++){
            sum+=pn_in[jj]*((ii>>jj)&1);}
        if(sum>0){
            ttt[ii]=true;}
        else if(sum<0){
            ttt[ii]=false;}
        else ttt[ii]=rand_fill;}
    return *this;
}
boolfun& boolfun::Gen_Domi(){
    bool rander=false;
    int ii,jj,nn0,nn1,code,Domin_Type;
    if(0==argus[0]||1==argus[0]){
        rander=false;
        Domin_Type=argus[0];}
    else {
        rander=true;}
    for(ii=0; ii<length; ii++){
        nn0=nn1=0;
        code=ii;
        for(jj=0; jj<k; jj++){
            if(code&1)nn1++;
            else nn0++;
            code>>=1;}
        if(nn0>nn1){
            ttt[ii]=false;}
        else if(nn0<nn1){
            ttt[ii]=true;}
        else {
            if(rander)ttt[ii]=unif_rand()>0.5;
            else ttt[ii]=Domin_Type;}}
    return *this;
}

bool boolfun::is_Cana(){// [NOTE] Only consider outermost layer.
    bool cana_0,cana_1,results=false,*ipt1,*ipt2;
    int ii,jj,kk,logi,boundary_l=1,boundary_u=length>>1;
    int flag_0=1,flag_1=1;
    for(ii=0; ii<k; ii++){
        ipt1=ttt;
        ipt2=ttt+boundary_l;
        cana_0=*ipt1;// Set 0-input canalized values.
        cana_1=*ipt2;// Set 1-input canalized values.
        logi=flag_0=flag_1=0;
        for(jj=0; jj<boundary_u; jj++){
            for(kk=0; kk<boundary_l; kk++){
                flag_0+=(ipt1[kk]!=cana_0);// Keep same, in-0.
                flag_1+=(ipt2[kk]!=cana_1);// keep same, in-1.
                if(flag_0&&flag_1){// Both 0/1 input are not canalized.
                    logi=1;break;}}
            if(logi)break;
            ipt1+=(boundary_l<<1);
            ipt2+=(boundary_l<<1);}
        if(logi){// Check next variable.
            boundary_l=boundary_l<<1;
            boundary_u=boundary_u>>1;}
        else {
            results=true;break;}}
    return results;
}
std::vector<std::vector<int>> boolfun::NestCana(){
    std::vector<bool> checker(length,true);
    std::vector<bool> global_logic(length,true);
    std::vector<int> Canalized(k,1);
    std::vector<std::vector<int>> res(3);
    res[0].clear();res[1].clear();res[2].clear();
    int ii, jj, masker, tag0, tag1;// GlobalID[length];
    int unstable=k;// layer=length>>1;
    bool lab1, lab0, failure, inner;
    //for(ii=0; ii<length; ++ii){GlobalID[ii]=ii;}
    failure=true;
    while(unstable>0&&failure){
        inner=true;
        failure=false;
        for(ii=0; (ii<k)&&inner; ++ii){
            if(0<Canalized[ii]){// Not been canalized.
                masker=1<<ii;
                tag1=-1; tag0=-1;
                lab1=true; lab0=true;
                for(jj=0; (jj<length)&&(lab1||lab0); ++jj){
                    if(checker[jj]){// Not been canalized upper layer?
                        if(0<(jj&masker)){// Tag:1
                            if(tag1<0){
                                tag1=ttt[jj];}
                            else {
                                lab1&=(ttt[jj]==tag1);}}
                        else {// Tag:0
                            if(tag0<0){
                                tag0=ttt[jj];}
                            else {
                                lab0&=(ttt[jj]==tag0);}}}}
                if(lab1||lab0){
                    unstable--;
                    inner=false;
                    failure=(lab0+lab1)<2;
                    if(lab0){
                        res[0].push_back(ii);// varaiable
                        res[1].push_back(0);// Canalizing
                        res[2].push_back(tag0);}// Canalized
                    if(lab1){
                        res[0].push_back(ii);// varaiable
                        res[1].push_back(1);// Canalizing
                        res[2].push_back(tag1);}// Canalized
                    Canalized[ii]=-1;
                    masker=1<<ii;
                    for(jj=0; jj<length; ++jj){
                        if(lab1&&((jj&masker)>0) )checker[jj]=false;
                        if(lab0&&((jj&masker)==0))checker[jj]=false;}}}}}
    return res;
}
bool boolfun::is_Post(){// [NOTE] Only consider the form of clique function.
    bool results=true,*ipt=ttt;
    int ii,jj,slot0,slot1,fuller=(1<<k)-1,flag_0=1,flag_1=1;
    std::vector<int> sets0(0),sets1(0);
    for(ii=0; ii<length; ii++, ipt++){
        if(*ipt==true&&flag_1){// is 1
            slot1=ii;
            for(jj=0; jj<((int)sets1.size()); jj++){
                if(slot1&sets1[jj]);// Have common "1"
                else flag_1=0;}
            if(flag_1)sets1.push_back(ii);}
        else if(*ipt==false&&flag_0){// is 0, note the inverse bits!!
            slot0=(~ii)&fuller;
            for(jj=0; jj<((int)sets0.size()); jj++){
                if(slot0&((~sets0[jj])&fuller));// Have common "0"
                else flag_0=0;}
            if(flag_0)sets0.push_back(ii);}// Record the actual code, not inverse code.
        else ;
        if(flag_0+flag_1==0){
            results=false;break;}}
    return results;
}
bool boolfun::is_Mone(){
    return boolfun::aux_MoneSign(false);
}
bool boolfun::is_Sign(){
    return boolfun::aux_MoneSign(true);
}
bool boolfun::aux_MoneSign(bool independent){
    int ii,jj,kk,logi0=1,logi1=1;
    int exit=0,lower=1,upper=length>>1;
    bool *ipt0,*ipt1;
    for(ii=0; ii<k; ii++){// Should ensure all inputs are SAME monotone type.
        ipt0=ttt;
        ipt1=ttt+lower;
        if(independent)logi0=logi1=1;
        for(jj=0; jj<upper; jj++){
            for(kk=0; kk<lower; kk++){
                if(ipt0[kk]>ipt1[kk])logi1=0;// Contradiction to increase
                if(ipt0[kk]<ipt1[kk])logi0=0;// Contradiction to decrease
                if(logi1+logi0==0){
                    exit=1;// Exit repeat judgment, not increase/decrease monotone.
                    break;}}
            if(exit)break;
            ipt0+=(lower<<1);
            ipt1+=(lower<<1);}
        if(exit)break;
        lower<<=1;// Shrinking range of analysis.
        upper>>=1;}
    return (logi0+logi1>0);
}
bool boolfun::is_Effc(){
    int results=0,logi;
    bool *ipt1,*ipt2;
    int ii,jj,kk,lower=1,upper=length>>1;
    for(ii=0; ii<k; ii++){
        ipt1=ttt;
        ipt2=ttt+lower;
        logi=0;
        for(jj=0; jj<upper; jj++){
            for(kk=0; kk<lower; kk++){
                if((ipt1[kk])!=(ipt2[kk])){// Just one case can meet condition.
                    logi=1;break;}}
            if(logi)break;
            ipt1+=(lower<<1);
            ipt2+=(lower<<1);}
        if(logi){// Meet condition
            lower=lower<<1;
            upper=upper>>1;
            results++;}
        else break;}
    return (results==k);
}
bool boolfun::is_Spin(){// [WARNING] require Z3 solver library
    bool res;
    int ii,jj;
    z3::context TT_V,TT_W;// Set z3's context.
    std::vector <z3::expr> vars,wars;
    for(ii=0; ii<k; ii++){// Set variable number of VAR slots.
        vars.push_back(TT_V.int_const(("a_" + std::to_string(ii)).c_str()));
        wars.push_back(TT_W.int_const(("b_" + std::to_string(ii)).c_str()));}
    // Build solver.
    z3::expr theta=TT_V.int_const("theta");
    z3::expr exprss1=theta;
    z3::solver solve_it1(TT_V);
    z3::expr gamma=TT_W.int_const("gamma");
    z3::expr exprss2=gamma;
    z3::solver solve_it2(TT_W);
    for(ii=0; ii<length; ii++){
        exprss1=theta;exprss2=gamma;
        for(jj=0; jj<k; jj++){
            if((ii&(1<<jj))){// This bit is one (effective).
                exprss1=exprss1+vars[jj];
                exprss2=exprss2+wars[jj];}}
        if(ttt[ii]){// TURE/1
            solve_it1.add(exprss1>=0);
            solve_it2.add(exprss2>0);}
        else {
            solve_it1.add(exprss1<0);
            solve_it2.add(exprss2<=0);}
        for(jj=0; jj<k; jj++){// Regulatory weight value are limited to ±1 (a_{xy}=±1).
            solve_it1.add(1==vars[jj]||-1==vars[jj]);
            solve_it2.add(1==wars[jj]||-1==wars[jj]);}
        solve_it1.add(0==theta);
        solve_it2.add(0==gamma);}// Thershold is zero.
    if(solve_it1.check()==z3::sat){// Zero-summation is set as 1.
        res=true;}
    else if(solve_it2.check()==z3::sat){// Zero-summation is set as 0.
        res=true;}
    else res=false;
    return res;
}
bool boolfun::is_Thre(bool PlusMinus,bool ShowIt){// [WARNING] require Z3 solver library
    bool res=false;
    int ii,jj;
    z3::context TT_T;
    std::vector <z3::expr> vars;
    for(ii=0; ii<k; ii++){// Set variable number of VAR slots.
        vars.push_back(TT_T.int_const(("a_" + std::to_string(ii)).c_str()));}
    // Build solver.
    z3::expr theta=TT_T.int_const("theta");
    z3::expr exprss=theta;
    z3::solver solve_it(TT_T);
    for(jj=0; jj<k; jj++){// a_{xy}!=0
        solve_it.add(vars[jj]!=0);}
    for(ii=0; ii<length; ii++){
        exprss=theta;
        for(jj=0; jj<k; jj++){
            if((ii&(1<<jj))){// This bit is one (effective).
                exprss=exprss+vars[jj];}
            else if(PlusMinus)exprss=exprss-vars[jj];}
        if(ttt[ii]){
            solve_it.add(exprss>0);}// The equation is .
        else {
            solve_it.add(exprss<=0);}}
    if(solve_it.check()==z3::sat){
        if(ShowIt){
            z3::model Models=solve_it.get_model();
            for(ii=0;ii<k;ii++){
                Rprintf("a_%d=%ld, ",ii,Models.eval(vars[ii]).get_numeral_int64());}
            Rprintf("theta=%ld\n",Models.eval(theta).get_numeral_int64());}
        res=true;}// theta: threshold.
    return res;
}
bool boolfun::is_Domi(){
    bool logi=true;// *ipt1=ttt, 
    int ii,jj,code,nn0,nn1;
    for(ii=0;ii<length;ii++){
        nn0=nn1=0;
        code=ii;
        for(jj=0; jj<k; jj++){
            if(code&1)nn1++;
            else nn0++;
            code>>=1;}
        if(((nn0>nn1)&&(ttt[ii]==true))||((nn0<nn1)&&(ttt[ii]==false))){
            logi=false;break;}}// The state not meet the state-major condition.
    return logi;
}

// Generate a list of Boolean funcitons with pointed OBFs.
void BatchGenerationOBF(std::vector<std::vector<bool>> *BF_List, 
    int *RandIndex, int *InDegs, char OBF_Type, int SysSize, int part,
    double P_bias, double OBFv1, double OBFv2){
    (*BF_List).resize(SysSize);
    // Set boolean functions
    boolfun a_bf;
    int ii, jj, bn_lens, tmp_ind;
    bool *bitmap=(bool*)malloc((1<<17)*sizeof(bool));// Enough large. (Actually k-input ≤ 16)
    // Ordered part.
    for(ii=0; ii<part; ++ii){// OBF part.
        tmp_ind=InDegs[RandIndex[ii]];
        if(0<tmp_ind){// Has input-edges
            a_bf.Configuration(OBF_Type, tmp_ind, 0.50, bitmap, OBFv1, OBFv2,-9,-9).Gen_BF().Reset();
            (*BF_List)[RandIndex[ii]].assign(bitmap,bitmap+(1<<tmp_ind));}
        else {// No input, should be only one fixed value.
            (*BF_List)[RandIndex[ii]].push_back(unif_rand()>0.50);}}
    free(bitmap); bitmap=nullptr;
    // Random part.
    for(ii=part; ii<SysSize; ++ii){// Random part.
        bn_lens=1<<InDegs[RandIndex[ii]];
        (*BF_List)[RandIndex[ii]].resize(bn_lens);
        for(jj=0; jj<bn_lens; ++jj){
            (*BF_List)[RandIndex[ii]][jj]=unif_rand()<P_bias;}}
}

// Auxiliary function of PolynomialFunction(), all parameters are integer
Rcpp::List PolynomialFunction_Int32(Rcpp::IntegerMatrix &VariableMat, Rcpp::IntegerVector &TTT, int Type){
    Rcpp::List res;
    int ii,jj;
    int nVar=VariableMat.ncol();// { (K,1), or (K,1)+(K,2), or (K,1)+(K,2)+(K,3), or ..., or 2^K-1 }
    int nLen=VariableMat.nrow();// 2^K length
    z3::context TT_T;
    std::vector<z3::expr> vars;
    z3::expr theta=TT_T.int_const("theta");
    z3::expr exprss=theta;
    z3::solver solve_it(TT_T);
    for(ii=0; ii<nVar; ++ii){// Set variable number of VAR slots.
        vars.push_back(TT_T.int_const(("a_" + std::to_string(ii)).c_str()));}
    switch (Type){
        case 0:// Boolean-based threshold.
            for(ii=0; ii<nLen; ++ii){
            exprss=theta;
            for(jj=0; jj<nVar; ++jj){
                    if(VariableMat(ii,jj)>0){// The variable matrix is ±1.
                        exprss=exprss+vars[jj];}}
                if(TTT[ii]>0){
                    solve_it.add(exprss>0);}
                else {
                    solve_it.add(exprss<=0);}}
                break;
        case 1:// Boolean-based polynomial form.
            for(ii=0; ii<nLen; ++ii){
                exprss=theta;
                for(jj=0; jj<nVar; ++jj){
                    if(VariableMat(ii,jj)>0){// The variable matrix is ±1.
                        exprss=exprss+vars[jj];}}
                if(TTT[ii]>0){
                    solve_it.add(1==exprss);}
                else {
                    solve_it.add(0==exprss);}}
                break;
        case 2:// Spin-like threshold.
        default:
            for(ii=0; ii<nLen; ++ii){
                exprss=theta;
                for(jj=0; jj<nVar; ++jj){
                    if(VariableMat(ii,jj)>0){// The variable matrix is ±1.
                        exprss=exprss+vars[jj];}
                    else {
                        exprss=exprss-vars[jj];}}
                if(TTT[ii]>0){
                    solve_it.add(exprss>0);}
                else {
                    solve_it.add(exprss<=0);}}
                break;}
    if(solve_it.check()==z3::sat){// Find out appropriate weights.
        res.push_back(true);// Return the logical tag.
        z3::model Models=solve_it.get_model();
        Rcpp::IntegerVector Weights(nVar,-666);
        for(ii=0; ii<nVar; ++ii){
            Weights[ii]=Models.eval(vars[ii]).get_numeral_int();}
        // Insert the variables' weights.
        res.push_back(Weights);
        // Return the threshold value.
        res.push_back(Models.eval(theta).get_numeral_int());
    }
    else {// Failure 
        res.push_back(false);}
    return (res);
}
// Auxiliary function of PolynomialFunction(), all parameters are double/float
Rcpp::List PolynomialFunction_Float(Rcpp::IntegerMatrix &VariableMat, Rcpp::IntegerVector &TTT, int Type){
    Rcpp::List res;
    int ii,jj;
    int nVar=VariableMat.ncol();// { (K,1), or (K,1)+(K,2), or (K,1)+(K,2)+(K,3), or ..., or 2^K-1 }
    int nLen=VariableMat.nrow();// 2^K length
    z3::context TT_T;
    std::vector<z3::expr> vars;
    z3::expr theta=TT_T.real_const("theta");
    z3::expr exprss=theta;
    z3::solver solve_it_r(TT_T);
    for(ii=0; ii<nVar; ++ii){// Set variable number of VAR slots.
        vars.push_back(TT_T.real_const(("a_" + std::to_string(ii)).c_str()));}
    for(ii=0; ii<nLen; ++ii){
        exprss=theta;
        for(jj=0; jj<nVar; ++jj){
            if(VariableMat(ii,jj)>0){// The variable matrix is ±1.
                exprss=exprss+vars[jj];}
            else {
                exprss=exprss-vars[jj];}}
        if(TTT[ii]>0){
            solve_it_r.add(TT_T.real_val("1.0")==exprss);}
        else {
            solve_it_r.add(TT_T.real_val("-1.0")==exprss);}}
    if(solve_it_r.check()==z3::sat){// Find out appropriate weights.
        res.push_back(true);// Return the logical tag.
        z3::model Models=solve_it_r.get_model();
        Rcpp::NumericVector Weights(nVar,-666.0);
        int numerator, denominator;
        for(ii=0; ii<nVar; ++ii){
            numerator=Models.eval(vars[ii]).numerator().get_numeral_int();
            denominator=Models.eval(vars[ii]).denominator().get_numeral_int();
            Weights[ii]=(double)(numerator)/denominator;}
            res.push_back(Weights);
            res.push_back(
                (double)(Models.eval(theta).numerator().get_numeral_int())/
                    Models.eval(theta).denominator().get_numeral_int());
    }
    else {// Failure 
        res.push_back(false);}
    return res;
}
// Analyze possible polynomial expressions of given Boolean function.
Rcpp::List PolynomialFunction(Rcpp::IntegerMatrix &VariableMat, Rcpp::IntegerVector &maptab, int IsSpinLike){
    if(3!=IsSpinLike){
        return(PolynomialFunction_Int32(VariableMat,maptab,IsSpinLike));}
    else {
        return(PolynomialFunction_Float(VariableMat,maptab,IsSpinLike));}
}

//Code is over!