#include "QuineMcCluskey.h"
// QuineMcCluskey Method Auxiliary Functions:
// Return number of invalid elements.
int qmc_aux_InvalidCount(std::vector<int> &Slot){
    int ii,sum=0;
    for(ii=0;ii<((int)Slot.size());ii++){
        sum+=(Slot[ii]==-1);}
    return sum;
}

// Can pointed bit/dit of Prime A be covered by Prime B's one?
std::string qmc_aux_ExactMatchOnce(struct Prime &a,const struct Prime &b,int label){
    int ii,sum=1;
    std::string res="Failed!";
    if(a.mapout<=b.mapout){// (2&a&b)=(1&a&b)|(2&a&b): expanding for multiple levels.
        for(ii=0;ii<((int)a.slot.size());ii++){
            if(ii==label){// Objected bit: Pointed/Non-pointed one should be different/same value.
                sum*=(a.slot[ii]!=b.slot[ii]);}
            else {
                sum*=(a.slot[ii]==b.slot[ii]);}
            if(!sum)break;}}
    else sum=0;
    if(sum){res.erase(0).append(b.Name);}
    return res;
}

// Prime A can be covered by Prime B but inversely.
std::string qmc_aux_ExactMatchOnceInverse(struct Prime &a,const struct Prime &b,int label){
    int ii,sum=1;
    std::string res="Failed!";
    if(a.mapout>=b.mapout){// Inverse form: 2:A0B1 ==> [1:A0B1] + 1:A1B1 + 1:A2B1, latter three can integrat as 1:A*B1  
        for(ii=0;ii<((int)a.slot.size());ii++){
            if(ii==label){// Objected bit.
                sum*=(a.slot[ii]!=b.slot[ii]);}
            else {
                sum*=(a.slot[ii]==b.slot[ii]);}
            if(!sum)break;}}
    else sum=0;
    if(sum){res.erase(0).append(b.Name);}
    return res;
}

// Return ID of a discrete vector belonging to which [cell].
int qmc_auc_Non0CellID(std::vector<int> &Counts,int *Order,int num){
    int sum=0;
    for(int ii=0;ii<num;ii++){
        sum=sum+Counts[ii+1]*Order[ii];}
    return sum;
}

// Count all non-zero numbers. [Boolean means "1"][Other means N"1" + N"2"+ ...]
int qmc_auc_Non0Summation(std::vector<int> &Counts, int Level){
    if(2==Level){// Boolean system.
        return Counts[1];}
    else {
        int sum=0;
        for(int ii=1; ii<((int)Counts.size()); ++ii)sum+=Counts[ii];
        return sum;}
}

// Can Prime A convered by Prime B?
bool qmc_aux_ConveredByPrime(struct Prime &a,struct Prime &b,int k){
    bool Convered=true;
    for(int ii=0;ii<k;++ii){
        if(-1!=b.slot[ii]&&a.slot[ii]!=b.slot[ii]){Convered=false;break;}}
    return Convered;
}

// Self-addation operation in vector.
void qmc_aux_Adder1(std::vector <int> &num,int BitLeng, int Level){
    int loc=0, Ceiling=Level-1;
    while(loc<BitLeng){
        if(num[loc]<Ceiling){num[loc]++; break;}
        else {num[loc]=0; loc++;}}
}

// Count different elements in one input vector.
void qmc_aux_Counter(std::vector<int> &Seqs, std::vector<int> &Counter,int lengths){
    for(int ii=0;ii<((int)Counter.size());++ii){
        Counter[ii]=0;}
    for(int ii=0;ii<lengths;ii++){
        if(Seqs[ii]>-1){// Covered signals
            Counter[Seqs[ii]]++;}}
}

// Change "the logic false -> true" of item named NNaammee.
void qmc_aux_ChangeLogic(std::set<struct Prime> &ss,std::string NNaammee){
    struct Prime Deleter,Keeper;
    Deleter.Name=NNaammee;
    auto where=ss.find(Deleter);
    Keeper=(*where);
    Keeper.Covered=true;
    ss.erase(Deleter);
    ss.insert(Keeper);
}

// Remove invalid labels (count the number). 
int qmc_aux_RemoveInvalidLabels(std::vector<int> &a,int target){
    int sum=0;
    for(int ii=0;ii<((int)a.size());ii++){
        if(a[ii]<target){
            a[ii]=-1;}
        else sum++;}
    return sum;
}

// Find the location or index. 
int qmc_aux_FirstFind(std::vector<int> &a,int target){
    int ii;
    for(ii=0;ii<((int)a.size());ii++){
        if(a[ii]==target)break;}
    return ii;
}


// Member of class [QuineMcCluskeyMethod]
// Initial a class and assign memory.
QuineMcCluskeyMethod& QuineMcCluskeyMethod::Inits(int *p){
    //lengths=1<<k;// 2^k
    int ii,jj;
    LL_system=(int*)malloc(L*sizeof(int));
    LL_system[0]=1;
    if(2==L){
        lengths=1<<k;// 2^k
        LL_system[1]=k+1;}
    else {
        lengths=(int)(pow(L,k));
        for(ii=0; ii<L-1; ++ii){// Set the slot of conter's mullter: {1,k+1,(k+1)^2,...,(k+1)^(L-1)}
            LL_system[ii+1]=LL_system[ii]*(k+1);}}
    FinalPrime.clear();
    TmpImplicant.clear();
    Records.clear();
    dicts.clear();
    maptab=(int*)malloc(lengths*sizeof(int));
    implicants=(struct Prime*)malloc(lengths*sizeof(struct Prime));
    //int tmp_scale[k];
    //for(ii=0; ii<k; ++ii)tmp_scale[ii]=0;
    std::vector<int> tmp_scale(k,0);
    for(ii=0; ii<lengths; ++ii){
        new (&implicants[ii]) struct Prime ();
        maptab[ii]=p[ii];
        implicants[ii].mapout=p[ii];
        implicants[ii].id=ii;
        implicants[ii].Covered=false;
        implicants[ii].slot=tmp_scale;
        implicants[ii].blank.resize(L,0);
        qmc_aux_Counter(implicants[ii].slot,implicants[ii].blank,k);
        implicants[ii].cell_id=qmc_auc_Non0CellID((&implicants[ii])->blank,LL_system,L-1);
        if(implicants[ii].Name.size()!=0)implicants[ii].Name.erase(0);
        for(jj=0;jj<k;jj++){// Set item's name!!
            implicants[ii].Name.append(std::to_string(tmp_scale[jj]));}
        qmc_aux_Adder1(tmp_scale,k,L);}
    dicts.resize((int)(pow(k+1,L-1)));
    for(ii=0; ii<((int)dicts.size()); ++ii){
        dicts[ii].clear();}
    QuineMcCluskeyMethod::PrepareAnalysis();
    return *this;
}
// Free memory and reset class.
QuineMcCluskeyMethod& QuineMcCluskeyMethod::Reset(){
    free(LL_system);
    LL_system=nullptr;
    free(maptab);
    maptab=nullptr;
    for(int ii=0;ii<lengths;ii++){
        implicants[ii].slot.clear();
        implicants[ii].blank.clear();
        implicants[ii].Name.clear();
        implicants[ii].~Prime();}
    free(implicants);
    implicants=nullptr;
    std::vector<struct Prime>().swap(FinalPrime);// Clear FinalPrime.
    TmpImplicant.clear();
    for(int ii=0;ii<((int)Records.size());ii++){
        Records[ii].clear();}
    for(int ii=0;ii<((int)dicts.size());ii++){
        dicts[ii].clear();}
    std::vector<std::set<struct Prime>>().swap(Records);// Clear FinalPrime.
    std::vector<std::set<struct Prime>>().swap(dicts);// Clear FinalPrime.
    return *this;
}

// Show a possible final DNF for Boolean system.
void QuineMcCluskeyMethod::ShowFinalPrime_B(){
    char letters[16]={'a','b','c','d', 'e','f','g','h', 'i','j','l','m', 'n','o','p','q'};
    char LETTERS[16]={'A','B','C','D', 'E','F','G','H', 'I','J','L','M', 'N','O','P','Q'};
    // C++'s output stream >>>
    /*printf("{%d:",FinalPrime[0].mapout);
    for(int jj=0;jj<k;jj++){
        if(FinalPrime[0].slot[jj]>-1){
            printf("%c=%d,",Letters[jj],FinalPrime[0].slot[jj]);}
    }printf("}\n");
    for(int ii=1;ii<FinalPrime.size();ii++){
        printf(" & {%d:",FinalPrime[ii].mapout);
        for(int jj=0;jj<k;jj++){
            if(FinalPrime[ii].slot[jj]>-1){
                printf("%c=%d,",Letters[jj],FinalPrime[ii].slot[jj]);}
        }printf("}\n");
    }*/
    // R's output stream >>>
    Rprintf("{ 1= ");
    for(int jj=0; jj<k; ++jj){
        if(1==FinalPrime[0].slot[jj]){// Constant functions have been checked in outer layer.
            Rprintf("%c",LETTERS[jj]);}
        else if(0==FinalPrime[0].slot[jj]){
            Rprintf("%c",letters[jj]);}}
    for(int ii=1; ii<(int)FinalPrime.size(); ++ii){
        Rprintf(" + ");
        for(int jj=0; jj<k; ++jj){
            if(1==FinalPrime[ii].slot[jj]){
                Rprintf("%c",LETTERS[jj]);}
            else if(0==FinalPrime[ii].slot[jj]){
                Rprintf("%c",letters[jj]);}}}
    Rprintf(" }\n");
}
// Return information of multi-valued QMC results.
Rcpp::IntegerMatrix QuineMcCluskeyMethod::ShowImplicants_M(){
    // ShowImplicantSet_M(FinalPrime,k,L);
    // char LETTERS[16]={'A','B','C','D', 'E','F','G','H', 'I','J','L','M', 'N','O','P','Q'};
    // Output table title >>>
    // for(ii=0; ii<k; ++ii){
    //     Rprintf(" %c",LETTERS[ii]);}
    // Rprintf(" f(~)\n");
    // Output final prime implicants >>>
    // for(auto iitt=FinalPrime.begin(); iitt!=FinalPrime.end(); ++iitt){
    //     for(ii=0; ii<k; ++ii){
    //         ((iitt->slot)[ii]>-1)?Rprintf("%2d",(iitt->slot)[ii]):Rprintf(" *");}
    //     Rprintf("%3d  \n",iitt->mapout);}
    int ii;
    Rcpp::IntegerMatrix Res((int)(FinalPrime.size()),k+1);
    int count=0;
    for(auto iitt=FinalPrime.begin(); iitt!=FinalPrime.end(); ++iitt){
        for(ii=0; ii<k; ++ii){
            Res(count,ii)=(iitt->slot)[k-1-ii];}
        Res(count,k)=iitt->mapout;
        count++;}
    return Res;
}
// Return final implicants' #/* location information.
std::vector<std::vector<int>> QuineMcCluskeyMethod::ReturnFinalPrimeInfo(){
    // [WARNING!!] The first index begins from "1" instead of "0" (this simplifies sign's labels).
    std::vector<std::vector<int>> Res;
    Res.resize(FinalPrime.size());
    std::vector<int> tmp_vec;
    for(int ii=0; ii<(int)FinalPrime.size(); ++ii){
        tmp_vec.clear();
        for(int jj=0; jj<k; ++jj){
            if(1==FinalPrime[ii].slot[jj]){// Constant functions have been checked in outer layer.
                tmp_vec.push_back(jj+1);}
            else if(0==FinalPrime[ii].slot[jj]){
                tmp_vec.push_back(-jj-1);}
            else ;}
        Res[ii]=(tmp_vec);}
    return Res;
}
// Set (n1=a)'s next plane (n1=a+1)
void QuineMcCluskeyMethod::BuildCandidate(int Cell_ID,std::set<struct Prime> &a){
    int ii,NeighborCell;
    std::set<Prime>::iterator iitt,inner_;
    std::vector<int> Orders;
    iitt=dicts[Cell_ID].begin();// Check non-zero in advance
    Orders.assign(iitt->blank.begin(),iitt->blank.end());
    a.clear();
    for(ii=1; ii<L; ++ii){// Only check the non-zero counters.
        if(k>Orders[ii]){// Should not over the boundary.
            Orders[ii]++;
            NeighborCell=qmc_auc_Non0CellID(Orders,LL_system,L-1);
            if(qmc_auc_Non0Summation(Orders,L)<=k){// If reach the hyper-plane?
                for(inner_=dicts[NeighborCell].begin();inner_!=dicts[NeighborCell].end();++inner_){
                    a.insert(*inner_);}}
            Orders[ii]--;}}
}
// Build a hyper cube with different ones.
void QuineMcCluskeyMethod::HyperCube(){
    for(int ii=0;ii<(int)dicts.size();ii++){// Clear old record.
        dicts[ii].clear();}
    for(auto iitt=TmpImplicant.begin();iitt!=TmpImplicant.end();++iitt){
        dicts[iitt->cell_id].insert(*iitt);}// Has calculate the cell's ID.
}
// Only consider input vectors with non-zero value.
void QuineMcCluskeyMethod::PrepareAnalysis(){
    TmpImplicant.clear();
    for(int ii=0;ii<lengths;ii++){
        if(implicants[ii].mapout>0){
            TmpImplicant.insert(implicants[ii]);}}
}
// Check each prime in current plane can be merged.
std::vector<int> QuineMcCluskeyMethod::AnalysisMode(struct Prime item,std::set<struct Prime> &a){
    int ii,jj,indicator;
    std::vector<int> ress;
    std::string marker;
    std::vector<std::string> NameSet;
    for(ii=0;ii<k;ii++){
        if(-1<item.slot[ii]){// "-1" means i-th variable has be reduced.
            NameSet.clear();
            indicator=0;// Each variable should reset the vector.
            for(auto iitt=a.begin();iitt!=a.end();++iitt){
                marker=qmc_aux_ExactMatchOnce(item,*iitt,ii);// item covered by *iitt ?
                if(marker!="Failed!"){
                    indicator++;
                    if(iitt->mapout==item.mapout){// 1:A0B1~(2:A1B1+1:A2B1) => 1AXB1, but 2:A1B1 should be kept.
                        NameSet.push_back(marker);}}
                if(indicator==(L-1)){// Break advanced.
                    break;}}
            if(indicator==(L-1)){// If it meets the condition
                ress.push_back(ii);// ii-th variable can be ignored.
                for(jj=0;jj<(int)NameSet.size();jj++){// Should immediately updating all TmpImplicant.
                    qmc_aux_ChangeLogic(TmpImplicant,NameSet[jj]);}}}}
    return ress;
}
// Check each prime in current plane can be merged but inversely!
std::vector<std::vector<std::string>> QuineMcCluskeyMethod::Inverse_AnalysisMode(
    struct Prime item,std::set<struct Prime> &a,std::vector<int> &mapto,std::vector<int> &WhichVar){
    int ii,indicator,min_er;
    std::string marker;
    std::vector<std::string> NameSet;
    std::vector<std::vector<std::string>> ress;
    mapto.clear();WhichVar.clear();
    // Clear the return.
    for(ii=0; ii<k; ++ii){// ii: variable's 
        if(-1<item.slot[ii]){// "-1" means i-th variable has be reduced.
            NameSet.clear();    indicator=0;    min_er=item.mapout;// Each variable should reset the vector.
            for(auto iitt=a.begin();iitt!=a.end();++iitt){
                marker=qmc_aux_ExactMatchOnceInverse(item,*iitt,ii);
                if(marker!="Failed!"){
                    indicator++;
                    if(iitt->mapout<min_er){// Only save items with current minimum mapto.
                        min_er=iitt->mapout;
                        NameSet.clear();
                        NameSet.push_back(marker);}
                    else if(iitt->mapout==min_er){// "iitt->mapout" is identical to current minimum.
                        NameSet.push_back(marker);}
                    else ;// Current iitt->mapout larger than min_er, no need to record.
                }
                if(indicator==(L-1))break;}// Break advanced.
            if(indicator==(L-1)){// Meet the condition ~~
                if(min_er<item.mapout){// Only save some minimum mapto results, to avoid duplicate.
                    ress.push_back(NameSet);
                    mapto.push_back(min_er);
                    WhichVar.push_back(ii);
                    for(int jj=0;jj<(int)NameSet.size();jj++){
                        qmc_aux_ChangeLogic(TmpImplicant,NameSet[jj]);}}}}}
    return ress;
}
// Once analyze the merging implicants.
int QuineMcCluskeyMethod::RecursionImplicant(){
    int ii,valid=0;
    std::vector<int> remove,Maptos,IDs;
    struct Prime Newer;
    // Set extra one slot for suitable index.
    std::set<struct Prime> pools;
    pools.clear();
    // Slot for next loop of TmpImplicant.
    std::set<struct Prime> New_TmpImplicant;
    New_TmpImplicant.clear();
    std::vector<std::vector<std::string>> Inverses;
    // Build a hyper-cube, divide the each cell's componetns (Refer to dicts[ii]).
    QuineMcCluskeyMethod::HyperCube();
    // Check each cell.
    for(ii=0;ii<(int)dicts.size();ii++){
        if(dicts[ii].size()>0){// If Cell_(ii) is non-empty.
            // Build the Candidate-pool for Cell_(ii)
            QuineMcCluskeyMethod::BuildCandidate(ii,pools);
            // Analyze each element of dicts[ii].
            for(auto iitt=dicts[ii].begin();iitt!=dicts[ii].end();iitt++){
                // Forward merging (In Boolean context, no inverse merging).
                remove=AnalysisMode(*iitt,pools);
                if(remove.size()){// One lower item can be merged in (single or multiple) higher items
                    qmc_aux_ChangeLogic(TmpImplicant,iitt->Name);
                    for(int jj=0;jj<(int)remove.size();++jj){
                        Newer= *iitt;
                        Newer.id=valid;
                        Newer.slot[remove[jj]]=-1;
                        // Newer.blank.resize(2);
                        Newer.blank.resize(L);
                        qmc_aux_Counter(Newer.slot,Newer.blank,k);
                        // Newer.cell_id=(Newer.blank)[1];
                        Newer.cell_id=qmc_auc_Non0CellID(Newer.blank,LL_system,L-1);
                        Newer.Covered=false;
                        Newer.Name.erase(0);
                        for(int kk=0;kk<k;kk++){
                            Newer.Name.append(std::to_string(Newer.slot[kk]));}
                        New_TmpImplicant.insert(Newer);
                        valid++;}
                }
                // Inverse merging.
                if((L>2)&&(!Incomparability)){// Pattern {0/1} or {max, min, delta}?
                    Inverses=Inverse_AnalysisMode(*iitt,pools,Maptos,IDs);
                    if(!(Inverses.size()==Maptos.size()&&Maptos.size()==IDs.size())){
                        Rcpp::stop("Maybe some error in checking inverse cover.\n");}
                    for(int kk=0; kk<(int)Inverses.size(); ++kk){
                        // In fact, only condsider the first item enough.
                        Newer.Name=Inverses[kk][0];
                        auto iivv=TmpImplicant.find(Newer);
                        Newer= *iivv;
                        Newer.id=valid;
                        Newer.slot[IDs[kk]]=-1;
                        qmc_aux_Counter(Newer.slot,Newer.blank,k);
                        Newer.cell_id=qmc_auc_Non0CellID(Newer.blank,LL_system,L-1);
                        Newer.Covered=false;
                        Newer.Name.erase(0);
                        for(int pp=0;pp<k;pp++){
                            Newer.Name.append(std::to_string(Newer.slot[pp]));}
                        New_TmpImplicant.insert(Newer);
                        valid++;}
                }
            }
        }
    }
    // Update echo record~~.
    pools.clear();// "pools" is reused as temporarily save non-convered items.
    if(valid>0){// Has new emerging items, should open next loop.
        for(auto iitt=TmpImplicant.begin();iitt!=TmpImplicant.end();++iitt){
            if(!(iitt->Covered))pools.insert(*iitt);}// Not covered, can record them.
        Records.push_back(pools);
        TmpImplicant.clear();
        TmpImplicant=New_TmpImplicant;}
    else {
        for(auto iitt=TmpImplicant.begin();iitt!=TmpImplicant.end();++iitt){
            if(!(iitt->Covered))pools.insert(*iitt);}
        Records.push_back(pools);}
    return valid;
}
// Recursively check potential merging items.
void QuineMcCluskeyMethod::RecursionProcess(){
    int flag=1;
    while(flag>0){
        flag=QuineMcCluskeyMethod::RecursionImplicant();}
}
// Analyze the covering situation.
std::vector<struct Prime> QuineMcCluskeyMethod::SetCoveredMatrix(){
    int ii,jj=0;
    std::vector<struct Prime> prime_implicant;// for convenience, modify the bool.
    // Covert each echo Record into all prime implicant;
    for(ii=0;ii<(int)Records.size();++ii){
        for(auto iitt=Records[ii].begin();iitt!=Records[ii].end();++iitt){
            prime_implicant.push_back(*iitt);
            prime_implicant[jj].id=jj;jj++;}}// Reset id for idenfying items.
    std::vector<std::set<int>> CoverMatrix_r(lengths);
    std::vector<std::vector<int>> CorssMatrix_b(lengths, std::vector<int> (prime_implicant.size(),-1));
    std::vector<int> CurrentMax(lengths), Non0(lengths);// int CurrentMax[lengths],Non0[lengths];
    memcpy(CurrentMax.data(),maptab,lengths*sizeof(int));
    for(ii=0;ii<lengths;ii++){
        CoverMatrix_r[ii].clear();
        if(implicants[ii].mapout){// non-zero
            for(jj=0;jj<(int)prime_implicant.size();jj++){
                if(qmc_aux_ConveredByPrime(implicants[ii],prime_implicant[jj],k)){// Original item can be covered.
                    CoverMatrix_r[ii].insert(jj);
                    CorssMatrix_b[ii][jj]=prime_implicant[jj].mapout;
                    if(CurrentMax[ii]<prime_implicant[jj].mapout){
                        CurrentMax[ii]=prime_implicant[jj].mapout;}}}
            Non0[ii]=qmc_aux_RemoveInvalidLabels(CorssMatrix_b[ii],CurrentMax[ii]);}}
    // First loop of Must Keep:
    std::set<int> iidd;// Original IDs.
    std::set<int> MustKeep;
    for(ii=0;ii<lengths;ii++){
        if(implicants[ii].mapout>0){
            if(Non0[ii]==1){
                MustKeep.insert(qmc_aux_FirstFind(CorssMatrix_b[ii],CurrentMax[ii]));}
            else {
                iidd.insert(ii);}}}
    // Repeat the checking.
    int doit=1;
    std::set<int> Deler;// Need delete elements.
    while(doit&&iidd.size()>0){
        doit=0;Deler.clear();
        for(auto iitt=iidd.begin();iitt!=iidd.end();++iitt){
            for(auto it2=CoverMatrix_r[*iitt].begin();it2!=CoverMatrix_r[*iitt].end();++it2){
                if(MustKeep.end()!=MustKeep.find(*it2)){
                    Deler.insert(*iitt);
                    doit=1;break;}}}
        for(auto iitt=Deler.begin();iitt!=Deler.end();++iitt){
            iidd.erase(*iitt);}}
    // Some coupled items.
    while(iidd.size()>0){
        // Just pointed one item.
        auto it0=iidd.begin();
        auto it1=CoverMatrix_r[*it0].begin();// Only keep the first, and continue analyze.
        MustKeep.insert(*it1);
        iidd.erase(*it0);
        doit=1;
        while(doit&&iidd.size()>0){
            doit=0;Deler.clear();
            for(auto iitt=iidd.begin();iitt!=iidd.end();++iitt){
                for(auto it2=CoverMatrix_r[*iitt].begin();it2!=CoverMatrix_r[*iitt].end();++it2){
                    if(MustKeep.end()!=MustKeep.find(*it2)){
                        Deler.insert(*iitt);
                        doit=1;break;}}}
            for(auto iitt=Deler.begin();iitt!=Deler.end();++iitt){
                iidd.erase(*iitt);}}}
    std::vector<struct Prime> Res;
    for(auto it=MustKeep.begin();it!=MustKeep.end();++it){
        Res.push_back(prime_implicant[*it]);}
    return Res;
}
// Do once Quine McCluskey analysis.
QuineMcCluskeyMethod& QuineMcCluskeyMethod::Do_QuineMcCluskey(){
    QuineMcCluskeyMethod::RecursionProcess();
    FinalPrime=QuineMcCluskeyMethod::SetCoveredMatrix();
    return *this;
}
// Calculate the effective based mapping to one vectors (Obtain redundant information, indirect methods)
double QuineMcCluskeyMethod::Effective(){// Only pick pointed item's "#" symbol (-1 in slot)
    double avg,counter,sum=0;
    for(int jj=0;jj<lengths;jj++){
        avg=0;counter=0;
        for(int ii=0;ii<(int)FinalPrime.size();ii++){
            if(qmc_aux_ConveredByPrime(implicants[jj],FinalPrime[ii],k)){
                counter=counter+1;
                avg+=qmc_aux_InvalidCount(FinalPrime[ii].slot);}}
        if(counter>1e-7){
            sum+=(avg/counter);}}
    return sum;
}
// Return the number of final primes.
int QuineMcCluskeyMethod::NumberFinalPrime(){
    return (int)(FinalPrime.size());
}
// Return each input's average "*/#" (Obtain effective information, direct methods)
void QuineMcCluskeyMethod::SingleEdgeConnect_bool(std::vector<double> &WeightConnect){
    double total;
    for(int jj=0;jj<lengths;jj++){
        if(maptab[jj]>0){// 1-type and 0-type both set 1, execute twice.
            total=0;
            std::vector<int> counters(k,0);
            for(int ii=0; ii<(int)FinalPrime.size(); ++ii){
                if(qmc_aux_ConveredByPrime(implicants[jj],FinalPrime[ii],k)){
                    total+=1.0;
                    for(int pp=0; pp<k; ++pp){
                        if(-1!=FinalPrime[ii].slot[pp]){// Not "*" symbols.
                            counters[pp]++;
                        }
                    }
                }
            }
            for(int pp=0; pp<k; ++pp){
                WeightConnect[pp]+=counters[pp]/total;
            }
        }
    }
}

// Get the effective connections of nodes.
double OnceEffective(int *maptab,int k, int L){
    QuineMcCluskeyMethod mpt(k,L,true);
    double sum=mpt.Inits(maptab).Do_QuineMcCluskey().Effective();
    mpt.Reset();
    return sum;
}
double toR_MultipleEffective(int *maptab, int k, int L, int lens){
    int ii;// tempint1[lens],tempint0[lens];
    std::vector<std::vector<int>> Mapto(L,std::vector<int>(lens,0));
    std::vector<double> Results(k,0);
    for(ii=0; ii<lens; ++ii){
        Mapto[maptab[ii]][ii]=1;}
    double sum=0;
    for(ii=0; ii<L; ++ii){
        sum+=OnceEffective(Mapto[ii].data(),k,L);}
    //double sum=OnceEffective(tempint1,k,L)+OnceEffective(tempint0,k,L);
    return k-sum/(double)(lens);
}
// Show the disjunctive normal form of Booelan function.
void toR_ShowBoolFunDNF(int *maptab,int k, Rcpp::CharacterVector VarNames){
    QuineMcCluskeyMethod mpt(k,2,true);
    mpt.Inits(maptab).Do_QuineMcCluskey();
    if(VarNames[0]==NA_STRING){// Not provided variable's name.
        mpt.ShowFinalPrime_B();}
    else {
        std::vector<std::vector<int>> tmp_id=mpt.ReturnFinalPrimeInfo();
        std::vector<std::string> Var_Name(VarNames.begin(), VarNames.end());
        Rprintf("{ 1= ");
        if(tmp_id[0][0]>0){
            Rprintf("%s",Var_Name[abs(tmp_id[0][0])-1].c_str());}
        else {
            Rprintf("~%s",Var_Name[abs(tmp_id[0][0])-1].c_str());}
        for(int jj=1; jj<(int)tmp_id[0].size(); ++jj){
            Rprintf("*");
            if(tmp_id[0][jj]>0){
                Rprintf("%s",Var_Name[abs(tmp_id[0][jj])-1].c_str());}
            else {
                Rprintf("~%s",Var_Name[abs(tmp_id[0][jj])-1].c_str());}}
        for(int ii=1; ii<(int)tmp_id.size(); ++ii){
            Rprintf(" + ");
            if(tmp_id[ii][0]>0){
                Rprintf("%s",Var_Name[abs(tmp_id[ii][0])-1].c_str());}
            else {
                Rprintf("~%s",Var_Name[abs(tmp_id[ii][0])-1].c_str());}
            for(int jj=1; jj<(int)tmp_id[ii].size(); ++jj){
                Rprintf("*");
                if(tmp_id[ii][jj]>0){
                    Rprintf("%s",Var_Name[abs(tmp_id[ii][jj])-1].c_str());}
                else {
                    Rprintf("~%s",Var_Name[abs(tmp_id[ii][jj])-1].c_str());}}}
       Rprintf("}\n");}
    mpt.Reset();
}
// Analyze each connect effective. (Boolean or multiple)
void OnceEdgeConnect(int *maptab, int k, int L,std::vector<double> &Slot){
    QuineMcCluskeyMethod mpt(k,L,true);// All should be incomparable.
    mpt.Inits(maptab).Do_QuineMcCluskey().SingleEdgeConnect_bool(Slot);
    mpt.Reset();
}
std::vector<double> toR_MultipleEdgeConnect(int *maptab,int k,int L,int lens){
    int ii;// tempint1[lens],tempint0[lens];
    std::vector<std::vector<int>> Mapto(L,std::vector<int>(lens,0));
    std::vector<double> Results(k,0);
    for(ii=0; ii<lens; ++ii){
        Mapto[maptab[ii]][ii]=1;}
    for(ii=0; ii<L; ++ii){
        OnceEdgeConnect(Mapto[ii].data(),k,L,Results);}
    double fLen=(double)(lens);
    for(ii=0; ii<k; ++ii){
        Results[ii]=Results[ii]/fLen;
    }
    return Results;
}
// The complexity of ternary/multiple expressions (Independent & Comparable)
double toR_BoolMulComplexity(int *maptab,int k,int L,int lens, bool incomparable){
    double Res=0;
    if(incomparable){
        int ii,sum=0;
        std::vector<std::vector<int>> Mapto(L,std::vector<int>(lens,0));
        for(ii=0; ii<lens; ++ii){
            Mapto[maptab[ii]][ii]=1;}
        for(ii=0; ii<L; ++ii){
            QuineMcCluskeyMethod mpt(k,L,incomparable);// Treated as binary.
            sum+=mpt.Inits(Mapto[ii].data()).Do_QuineMcCluskey().NumberFinalPrime();
            mpt.Reset();}
        Res=(double)(sum);
    }
    else {
        QuineMcCluskeyMethod mpt(k,L,incomparable);// Elements are comparable.
        mpt.Inits(maptab).Do_QuineMcCluskey();
        Res=(double)(mpt.NumberFinalPrime());
        mpt.Reset();}
    return Res;
}
// Show the disjunctive normal form of Booelan function.
Rcpp::IntegerMatrix toR_ShowMulVFunDNF(int *maptab,int k, int L){
    QuineMcCluskeyMethod mpt(k,L,false);// Should use comparable elements.
    mpt.Inits(maptab).Do_QuineMcCluskey();
    Rcpp::IntegerMatrix Res=mpt.ShowImplicants_M();
    mpt.Reset();
    return Res;
}
// Code is over.
