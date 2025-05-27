#include "NetGraphFrame.h"
#include <Rcpp.h>

void net_chooseKfromN(int N,int K,int *Selected,int *Label,int flag){
    int ii;
    if(K<N/2){// Little selected.
        int tmp;
        if(flag>=0)Label[flag]=1;
        for(ii=0; ii<K; ++ii){
            do{
                tmp=(int)(N*unif_rand());
            }while(Label[tmp]);
            Selected[ii]=tmp;
            Label[tmp]=1;}
        if(flag>=0)Label[flag]=0;// Recover the slot Label
        for(ii=0; ii<K; ++ii){
            Label[Selected[ii]]=0;}}
    else {// Large selection.
        std::vector<int> elements(N);
        int vals, idss;
        for(ii=0; ii<N; ++ii)elements[ii]=ii;
        //std::shuffle(elements.begin(), elements.end(), mt);
        for(ii=N-1; ii>=0; ii=ii-1){
            idss=(int)(unif_rand()*ii);
            vals=elements[ii];
            elements[ii]=elements[idss];
            elements[idss]=vals;}
        memcpy(Selected,elements.data(),sizeof(int)*K);}// Only choose former K's.
}
// Insert a parent-type node.
void InsertParentsNode(struct Regulation *p, int parent){
    struct Regulation *tmp=(struct Regulation *)malloc(sizeof(struct Regulation));
    tmp->code=parent;
    if(p->prev==nullptr){// No parent
        p->prev=tmp;
        tmp->prev=nullptr;
        tmp->next=p;}
    else {// Exist parents
        p->prev->next=tmp;
        tmp->prev=p->prev;
        tmp->next=p;
        p->prev=tmp;}
}
// Insert a child-type node.
void InsertChildsNode(struct Regulation *p,int child){
    struct Regulation *tmp=(struct Regulation *)malloc(sizeof(struct Regulation));
    tmp->code=child;
    if(p->next==nullptr){// No child
        p->next=tmp;
        tmp->prev=p;
        tmp->next=nullptr;}
    else {// Exist children
        p->next->prev=tmp;
        tmp->next=p->next;
        p->next=tmp;
        tmp->prev=p;}
}
// Build a regulation edge/link.
void BuildRegulationship(struct Regulation *net,int child,int parent){
    InsertParentsNode(&net[child],parent);
    InsertChildsNode(&net[parent],child);
}
// Build non-directional edge (Only use child, not contain the parent-type)
void BuildNeighborEdge(struct Regulation *net,int node1,int node2){
    InsertChildsNode(&net[node1],node2);
    InsertChildsNode(&net[node2],node1);
}
// Remove one specified child from parent node.
int DeleteChilds(struct Regulation *p,int child){
    struct Regulation *npt=p,*Dels;
    int GotIt=0;
    while(npt->next!=nullptr){// Delete child target in this parents node
        if(npt->next->code==child){// Find where target node
            GotIt=1;
            if(npt->next->next==nullptr){// At the last location
                free(npt->next);
                npt->next=nullptr;}
            else {// Not at last
                npt->next->next->prev=npt;
                Dels=npt->next;
                npt->next=npt->next->next;
                free(Dels);}
            break;}
        else npt=npt->next;}
    return GotIt;
}
// Remove one specified parent from child node.
int DeleteParents(struct Regulation *p,int parent){
    struct Regulation *npt=p,*Dels;
    int GotIt=0;
    while(npt->prev!=nullptr){// Delete parent target in this son node
        if(npt->prev->code==parent){// Find where target vertex
            GotIt=1;
            if(npt->prev->prev==nullptr){// At the first location
                free(npt->prev);
                npt->prev=nullptr;}
            else {// Not at first
                npt->prev->prev->next=npt;
                Dels=npt->prev;
                npt->prev=npt->prev->prev;
                free(Dels);}
            break;}
        else npt=npt->prev;}
    return GotIt;
}
// Node has No.xx child/parent.
int ExistChild(struct Regulation *p,int child){
    struct Regulation *npt=p; int inhere=0;
    while(npt->next!=nullptr){
        if(child==npt->next->code){inhere=1;break;}
        else npt=npt->next;}
    return inhere;
}
int ExistParent(struct Regulation *p,int parent){
    struct Regulation *npt=p; int inhere=0;
    while(npt->prev!=nullptr){
        if(parent==npt->prev->code){inhere=1;break;}
        else npt=npt->prev;}
    return inhere;
}
// Build a mutual loop of all nodes.
void LinkMutualLoop(struct Regulation *Loops,int total){
    int ii;struct Regulation *rpt=Loops;
    for(ii=0;ii<total;ii++,rpt++)rpt->code=ii;
    Loops[0].next=Loops+1;
    Loops[0].prev=Loops+total-1;
    Loops[total-1].next=Loops;
    Loops[total-1].prev=Loops+total-2;
    for(ii=1;ii<total-1;ii++){
        Loops[ii].next=Loops+ii+1;
        Loops[ii].prev=Loops+ii-1;}
}
void LinkRomoveCurrent(struct Regulation *rpt){
    rpt->prev->next=rpt->next;
    rpt->next->prev=rpt->prev;
    rpt->next=rpt->prev=rpt;// Current node points itself!
}
// First-in first-out to queue of distance.
void ShortestFind_FIFO(struct Regulation *ParentNode,int Code,int logi){// Priority queue for shortest path algorithm.
    struct Regulation *tmp=ParentNode;
    if(logi>0){// Insert candidate nodes at end.
        struct Regulation *Newer=(struct Regulation *)malloc(sizeof(struct Regulation));
        Newer->next=Newer->prev=nullptr;Newer->code=Code;
        while(tmp->next!=nullptr){// Follow the last member.
            tmp=tmp->next;}
        tmp->next=Newer;Newer->prev=tmp;}
    else {// Delete the queue's first.
        if(ParentNode->next->next==nullptr){// Only one member.
            tmp=ParentNode->next;ParentNode->next=nullptr;}
        else {tmp=tmp->next;
            ParentNode->next->next->prev=ParentNode;
            ParentNode->next=ParentNode->next->next;}
        free(tmp);}
}

int Replaced_Poisson(double lambda){
    double SampleR=unif_rand();
    int k=0;
    double p=std::exp(-lambda);
    double F=p;
    while(F<SampleR){
        k+=1;
        p*=(lambda/k);
        F+=p;}
    return k;
}

NetGraphFrame& NetGraphFrame::ConfigurationBuildNet(char Net_Type, int Size, int ipar1, int ipar2, double fpar1){
    total=Size;
    Type=Net_Type;
    NetGraphFrame::Initial_Network();
    char tmpBuffer[32];
    tmpBuffer[0]=Net_Type;tmpBuffer[1]='\0';
    NetPara.append(tmpBuffer).append("-type, Size=");
    sprintf(tmpBuffer, "%d", Size);
    NetPara.append(tmpBuffer);
    switch(Net_Type){
        case 'N':// a null-model
            NetPara.append(" Null Model").push_back('\n');
            break;
        case 'K':case 'L':case 'R':case 'B':// Kauffman, Lattice, Regular-Random, Barabasi-Albert
            tmp_i_arg[0]=ipar1;// "K" parameter or "Core" of BA.
            sprintf(tmpBuffer, "%d", tmp_i_arg[0]);
            NetPara.append(", iPara=").append(tmpBuffer).push_back('\n');
            break;
        case 'E':// Erdos-Renyi
            tmp_f_arg[0]=fpar1;// "AvDeg" of ER.
            sprintf(tmpBuffer, "%.4f", tmp_f_arg[0]);
            NetPara.append(", fPara=").append(tmpBuffer).push_back('\n');
            break;
        case 'S':// Scale-free
            tmp_i_arg[0]=ipar1;// MaxDeg
            tmp_i_arg[1]=ipar2;// MinDeg
            tmp_f_arg[0]=fpar1;// Gamma, scaling degree.
            sprintf(tmpBuffer, "%d", tmp_i_arg[0]);
            NetPara.append(", iPara1=").append(tmpBuffer);
            sprintf(tmpBuffer, "%d", tmp_i_arg[1]);
            NetPara.append(", iPara2=").append(tmpBuffer);
            sprintf(tmpBuffer, "%.4f", tmp_f_arg[0]);
            NetPara.append(", fPara=").append(tmpBuffer).push_back('\n');
            break;
        default: Rcpp::stop("Error Type of Network.\n");}
    return *this;
}
NetGraphFrame& NetGraphFrame::Build_GraphNet(){
    switch(Type){
        case 'K':
            NetGraphFrame::Net_NK_D(tmp_i_arg[0]);break;
        case 'E':
            NetGraphFrame::Net_ER_D(tmp_f_arg[0]);break;
        case 'S':
            NetGraphFrame::Net_SF_U(tmp_i_arg[0],tmp_i_arg[1],tmp_f_arg[0]);break;
        case 'B':
            NetGraphFrame::Net_BA_U(tmp_i_arg[0]);break;
        case 'R':
            NetGraphFrame::Net_RR_D(tmp_i_arg[0]);break;
        case 'L':
            NetGraphFrame::Net_LT_D(tmp_i_arg[0]);break;
        case 'N':
            break;
        default: 
            Rcpp::stop("Error Type of Network.\n");}
    return *this;
}
NetGraphFrame& NetGraphFrame::Out2VecVecIntFrame(std::vector<std::vector<int>> *InEdge,
    std::vector<std::vector<int>> *OutEdge){// Output in-deg/out-deg data.
    struct Regulation *npt,*net=Network;
    for(int ii=0;ii<total;ii++,net++){
        npt=net;
        while(npt->prev!=nullptr){
            (*InEdge)[ii].push_back(npt->prev->code);
            npt=npt->prev;}
        npt=net;
        while(npt->next!=nullptr){
            (*OutEdge)[ii].push_back(npt->next->code);
            npt=npt->next;}}
    return *this;
}
NetGraphFrame& NetGraphFrame::LoadFromVecVecIntFrame(std::vector<std::vector<int>> &InEdge,
    std::vector<std::vector<int>> &OutEdge){// Output in-deg/out-deg data.
    // [NOTE] Null model has been built over; all memory space have been created. 
    for(int ii=0; ii<total; ++ii){
        InDeg[ii]=InEdge[ii].size();
        OtDeg[ii]=OutEdge[ii].size();
        for(int jj=0; jj<((int)InEdge[ii].size()); ++jj){// [NOTE] only add itself to avoid repeated 
           InsertParentsNode(Network+ii, InEdge[ii][jj]);}// or set &((Network)[ii])
        for(int jj=0; jj<((int)OutEdge[ii].size()); ++jj){// [NOTE] only add itself to avoid repeated 
           InsertChildsNode(Network+ii, OutEdge[ii][jj]);}}
    return *this;
}

NetGraphFrame& NetGraphFrame::IsolatedPointedNode_D(int ID){
    struct Regulation *tmp=Network+ID,*del;
    int loc,ii;// [NOTE] Some self-loop edges can be
    // Remove the parents & changed the degrees
    tmp=Network+ID;
    while(tmp->prev!=nullptr){
        loc=tmp->prev->code;
        ii=DeleteChilds(&Network[loc],ID);
        OtDeg[loc]--;
        tmp=tmp->prev;}
    // Remove the children & changed the degrees
    tmp=Network+ID;
    while(tmp->next!=nullptr){
        loc=tmp->next->code;
        ii=DeleteParents(&Network[loc],ID);
        InDeg[loc]--;
        tmp=tmp->next;}
    // Remove itself structure
    InDeg[ID]=OtDeg[ID]=0;// [NOTE] Is or not self-loop node both can be set.
    tmp=Network+ID;
    while(tmp->prev!=nullptr){// Delete parents
        del=tmp->prev;tmp->prev=tmp->prev->prev;free(del);}
    while(tmp->next!=nullptr){// Delete children
        del=tmp->next;tmp->next=tmp->next->next;free(del);}
    return *this;
}
NetGraphFrame& NetGraphFrame::BreakDownPointedEdge_D(int source, int target){
    int ii;
    ii=DeleteChilds(&Network[source],target);
    ii=DeleteParents(&Network[target],source);
    return *this;
}

// Initial one network (Allocate memory and initialization).
void NetGraphFrame::Initial_Network(){
    Network=(struct Regulation *)malloc(total*sizeof(struct Regulation));
    InDeg=(int *)malloc(4*total*sizeof(int));// InDeg + OtDeg + Coding + ReserveSlot
    OtDeg=InDeg+total;
    Coding=OtDeg+total;
    Address=(int **)malloc(2*(total+1)*sizeof(int *));// For quick build, exchange, &
    int **iptt=Address+total+1;
    struct Regulation *npt=Network;
    for(int ii=0;ii<total;ii++,npt++){
        InDeg[ii]=OtDeg[ii]=0;
        Coding[ii]=ii;
        iptt[ii]=Coding+ii;
        npt->code=ii;
        npt->prev=npt->next=nullptr;}
    iptt[total]=nullptr;
}
// Reset one network (Only delete and free the relations).
void NetGraphFrame::Reset_Network(){
    struct Regulation *del,*npt=Network;
    if(npt!=nullptr){
        for(int ii=0;ii<total;ii++,npt++){
            InDeg[ii]=OtDeg[ii]=0;
            while(npt->prev!=nullptr){// Delete parents
                del=npt->prev;npt->prev=npt->prev->prev;free(del);}
            while(npt->next!=nullptr){// Delete children
                del=npt->next;npt->next=npt->next->next;free(del);}}}
}
// Show one network.
void NetGraphFrame::Show_Network(){
    struct Regulation *npt,*net=Network;
    for(int ii=0;ii<total;ii++,net++){
        npt=net;
        while(npt->prev!=nullptr){
            Rprintf("%d->",npt->prev->code);npt=npt->prev;}
        npt=net;
        Rprintf("[%d]",npt->code);
        while(npt->next!=nullptr){
            Rprintf("->%d",npt->next->code);npt=npt->next;}
        Rprintf("\n");}
}
// Show each in/out degree.
void NetGraphFrame::Show_Degree(){
    int *ind=InDeg,*otd=OtDeg;
    for(int ii=0;ii<total;ii++,ind++,otd++){
        Rprintf("[%d]  %d,%d\n",ii,*ind,*otd);}
}
// Count in/out degree of each vertex.
void NetGraphFrame::Count_InOt_Deg(){
    struct Regulation *net=Network,*npt;
    for(int ii=0;ii<total;ii++,net++){
        npt=net;
        while(npt->next!=nullptr){// Count "npt->prev" is same for directional net.
            InDeg[npt->next->code]++;
            OtDeg[ii]++;
            npt=npt->next;}}
}
// Distance table of one node to all (Dijkstra, single source).
void NetGraphFrame::ShortestPath(int source){// Priority queue, not save routes
    int tmpI,tmp_source,tmp_node;
    int *Distance=OtDeg+(total<<1);
    struct Regulation *pt=nullptr,*head=(struct Regulation *)malloc(sizeof(struct Regulation));
    head->next=head->prev=nullptr;head->code=-1;
    for(tmpI=0;tmpI<total;tmpI++){Distance[tmpI]=total<<1;}// Initialization
    Distance[source]=0;
    ShortestFind_FIFO(head,source,666);// Insert source node
    if(Network[source].next!=nullptr){
        while(1){
            if(head->next==nullptr)break;
            tmp_source=head->next->code;
            ShortestFind_FIFO(head,-888,-666);// Out of queue (2nd Argu useless)
            pt=&Network[tmp_source];
            while(pt->next!=nullptr){
                tmp_node=pt->next->code;
                if(Distance[tmp_source]+1<Distance[tmp_node]){// A shorter path
                    Distance[tmp_node]=Distance[tmp_source]+1;
                    ShortestFind_FIFO(head,tmp_node,666);}// Insert candidate node
                pt=pt->next;}}}
    free(head);
}
// Generate a directed Erdos-Renyi graph.
void NetGraphFrame::Net_ER_D(double AvDeg){
    //std::poisson_distribution<int> PoissonDist((int)AvDeg);
    int ii,tmp,tmp1=0,tmp2=0;
    int *ipt_din=InDeg,*ipt_dot=OtDeg;
    // Set Poisson Distribution of IN/OUT.
    for(ii=0;ii<total;ii++){
        //tmp=PoissonDist(mt);
        tmp=Replaced_Poisson(AvDeg);
        if(0==tmp)tmp=1;
        if(tmp>=total)tmp=total-1;// Should not be larger than total-1.
        *ipt_dot++=tmp;
        //tmp=PoissonDist(mt);
        tmp=Replaced_Poisson(AvDeg);
        if(0==tmp)tmp=1;
        if(tmp>=total)tmp=total-1;
        *ipt_din++=tmp;}
    // Balance the in-out degree.
    ipt_din=InDeg;ipt_dot=OtDeg;
    for(ii=0;ii<total;ii++){
        tmp1+=*ipt_din++;tmp2+=*ipt_dot++;}
    if(tmp1>tmp2){// IN is larger than OUT.
        while(tmp1>tmp2){
            OtDeg[(int)(unif_rand()*total)]++;tmp2++;}}
    else if(tmp1<tmp2){// OUT is larger than IN.
        while(tmp1<tmp2){
            InDeg[(int)(unif_rand()*total)]++;tmp1++;}}
    else ;// Exactly balance
    NetGraphFrame::Address_Reset();
    NetGraphFrame::Auxiliary_ER_D();
}
// Generate a Barabasi-Albert network (Core+ Degree priority).
void NetGraphFrame::Net_BA_U(int Core){
    int ii,jj,tmp1,tmp2,TotalDeg;
    for(ii=0;ii<Core;ii++){
        InDeg[ii]=Core;
        for(jj=ii+1;jj<=Core;jj++){
            BuildNeighborEdge(Network,ii,jj);}}
    InDeg[Core]=Core;
    TotalDeg=(Core+1)*Core;
    for(ii=Core+1;ii<total;ii++){
        for(jj=0;jj<Core;jj++){
            while(1){
                tmp1=-1;
                tmp2=(int)(TotalDeg*unif_rand())+1;
                while(tmp2>0){// Large degree first.
                    tmp1++;
                    tmp2=tmp2-InDeg[tmp1];}
                if(ExistChild(&Network[tmp1],ii));// Has existed in neighbor table.
                else {// Can be as neighbors.
                    InDeg[tmp1]++;
                    BuildNeighborEdge(Network,tmp1,ii);
                    TotalDeg++;
                    break;}}}
        TotalDeg+=Core;// Now add Core, because of ignoring self-loop.
        InDeg[ii]=Core;}
}
// Generate a scale free network (assigning power law distribution).
void NetGraphFrame::Net_SF_U(int MaxDeg,int MinDeg,double Gamma){
    int tmpI,tmpJ,max_deg=0,*DegSlot;
    int *ipt_din=InDeg,*ipt_dot=OtDeg;
    double gamma=1-Gamma,exp_gamma=1.0/gamma;
    double min_bin=pow(MinDeg,gamma),max_min=pow(MaxDeg+1,gamma)-min_bin;
    for(tmpI=0;tmpI<total;tmpI++){
        tmpJ=(int)floor(pow(max_min*unif_rand()+min_bin,exp_gamma));
        if(max_deg<tmpJ)max_deg=tmpJ;// Find the max degree.
        ipt_dot[tmpI]=tmpJ;}
    DegSlot=(int *)malloc((max_deg+1)*sizeof(int));
    memset(DegSlot,0,(max_deg+1)*sizeof(int));// Set all zeros.
    for(tmpI=0;tmpI<total;tmpI++)DegSlot[ipt_dot[tmpI]]++;// Count degree distribution.
    tmpI=max_deg;
    while(tmpI>0){// Order degree from out-degree to in-degree.
        if(DegSlot[tmpI]>0){
            tmpJ=DegSlot[tmpI];
            while(tmpJ>0){
                *ipt_din++=tmpI;*ipt_dot++=tmpI;tmpJ--;}}
        tmpI--;}// Out-degree for assigning neighbors.
    // Build each edge by assigned degree.
    NetGraphFrame::Address_Reset();
    NetGraphFrame::Auxiliary_SF_U();
    free(DegSlot);
}
// Generate a directed regular random network.[V]
void NetGraphFrame::Net_RR_D(int K){
    int *ipt_din=InDeg,*ipt_dot=OtDeg;
    for(int ii=0;ii<total;ii++){
        *ipt_dot++=*ipt_din++=K;}
    NetGraphFrame::Address_Reset();
    NetGraphFrame::Auxiliary_RR_D();
}
// Generate a lattice network quickly.
void NetGraphFrame::Net_LT_D(int K){
    int code,ii;// Periodic boundary condition
    int *ipt_din=InDeg,*ipt_dot=OtDeg;
    std::vector<int> Neighbor(K);//int Neighbor[K];// *ipt_par=(int *)malloc(K*sizeof(int));
    tmp_i_arg[1]=(int)(sqrt(total));// [Note] Length=Height=sqrt(total);
    for(code=0; code<total; code++, ipt_din++, ipt_dot++){
        *ipt_din=K;*ipt_dot=K;
        NetGraphFrame::Auxiliary_LT_D(code,Neighbor.data());
        for(ii=0; ii<K; ii++){
            BuildRegulationship(Network,code,Neighbor[ii]);}}
}
// Generate a Kuaffman network quickly.
void NetGraphFrame::Net_NK_D(int K){
    int *ipt,*ipt_din=InDeg,*ipt_dot=OtDeg,*ipt_par=(int *)malloc(K*sizeof(int));
    int ii,code;
    int *labels=(int *)malloc(total*sizeof(int));
    ipt=labels;
    for(code=0;code<total;code++){// Degree & label initialization.
        *ipt_din++=K;*ipt_dot++=0;*ipt++=0;}
    for(code=0;code<total;code++){
        net_chooseKfromN(total,K,ipt_par,labels,code);// Assign each node's parents.
        ipt=ipt_par;
        for(ii=0;ii<K;ii++){OtDeg[*ipt++]++;}// Update the out-degree.
        for(ii=0;ii<K;ii++){// Regulations are directed.
            BuildRegulationship(Network,code,ipt_par[ii]);}}
    free(labels);
}
// Rest address.
void NetGraphFrame::Address_Reset(){
    memcpy(Address,Address+total+1,(total+1)*sizeof(int *));
}
// Rewiring of ER direct network.
void NetGraphFrame::Auxiliary_ER_D_sub(int *OtLab,int child,int parent){
    int tmp1,tmp2; struct Regulation *npt;
    while(1){
        tmp1=child;// Find a third node that not child & parent, and must be done-node.
        while(tmp1==child||tmp1==parent||OtLab[tmp1])tmp1=(int)(unif_rand()*total);
        tmp2=(int)(unif_rand()*OtDeg[tmp1]);// Should the find the replaced child of node [tmp1].
        npt=&Network[tmp1];
        while(tmp2>0){npt=npt->next;tmp2--;}
        tmp2=npt->next->code;
        if(tmp2!=parent&&!ExistChild(&Network[parent],tmp2)){// Find it! [Mod_20230709: add the !Exist...
            DeleteParents(&Network[tmp2],tmp1);
            DeleteChilds(&Network[tmp1],tmp2);
            BuildRegulationship(Network,child,tmp1);
            BuildRegulationship(Network,tmp2,parent);
            break;}}
}
// Assign neighbor with directed network.
void NetGraphFrame::Auxiliary_ER_D(){
    int Re_In,Re_Ot=total,N_INT=total*sizeof(int);
    int *tmp_ind=(int *)malloc(N_INT);
    int *tmp_otd=(int *)malloc(N_INT);
    int tmpI,tmp1,tmp2;
    memcpy(tmp_ind,InDeg,N_INT);// Copy temporary degrees.
    memcpy(tmp_otd,OtDeg,N_INT);// Copy temporary degrees.
    for(tmpI=0;tmpI<total;tmpI++){// Loop for assigning each node's parents.
        if(tmp_ind[tmpI]>=Re_Ot){break;}// Not enough to set parents, the equal also.
        while(tmp_ind[tmpI]){
            tmp1=(int)(unif_rand()*(Re_Ot));// Search.
            tmp2=**(Address+tmp1);// Candidate
            if(tmp2==tmpI||ExistChild(&Network[tmp2],tmpI));// Failure! (tmp2 -> tmpI) or self-loop
            else {
                BuildRegulationship(Network,tmpI,tmp2);
                tmp_ind[tmpI]--;// Child minus one.
                tmp_otd[tmp2]--;// Parent minus one.
                if(0==tmp_otd[tmp2]){
                    Re_Ot--;
                    memcpy(Address+tmp1,Address+tmp1+1,(Re_Ot-tmp1)*sizeof(int *));}}}}
    Re_In=total-tmpI;// Some unsatisfied nodes exist.
    while(Re_In){
        while(tmp_ind[tmpI]){
            tmp1=(int)(unif_rand()*Re_Ot);// Search.
            tmp2=**(Address+tmp1);// Candidate
            if(tmp2==tmpI||ExistChild(&Network[tmp2],tmpI)){
                NetGraphFrame::Auxiliary_ER_D_sub(tmp_otd,tmpI,tmp2);}
            else {
                BuildRegulationship(Network,tmpI,tmp2);}
            tmp_ind[tmpI]--;// Child minus one.
            tmp_otd[tmp2]--;// Parent minus one.
            if(0==tmp_otd[tmp2]){
                Re_Ot--;
                memcpy(Address+tmp1,Address+tmp1+1,(Re_Ot-tmp1)*sizeof(int *));}}
        Re_In--;tmpI++;}
    free(tmp_ind);free(tmp_otd);
}
// Assign suitable neighbors for SF network without degree priority.
void NetGraphFrame::Auxiliary_SF_U(){
    int ii,tmp1,tmp2,Remain=total;
    for(ii=0;ii<total;ii++){
        if(OtDeg[ii]){
            if(OtDeg[ii]>=Remain)break;// Remain is not enough to be assigned.
            memcpy(Address,Address+1,Remain*sizeof(int *));
            while(OtDeg[ii]){
                tmp1=(int)(unif_rand()*(Remain-1));// Forward search.
                tmp2=**(Address+tmp1);// Candidate
                if(ExistChild(&Network[ii],tmp2)){;}// Failure!
                else {
                    BuildNeighborEdge(Network,ii,tmp2);
                    OtDeg[ii]--;OtDeg[tmp2]--;
                    if(0==OtDeg[tmp2]){
                        Remain--;
                        memcpy(Address+tmp1,Address+tmp1+1,(Remain-tmp1)*sizeof(int *));}}}
            Remain--;}}
    while(Remain){// Some unsatisfied nodes exist.
        for(tmp1=ii;tmp1<total;tmp1++){
            if(OtDeg[tmp1]){// Unsatisfied from large to small degree.
                tmp2=0;
                while(OtDeg[tmp1]){
                    if(ExistChild(&Network[tmp1],tmp2));
                    else {
                        BuildNeighborEdge(Network,tmp1,tmp2);
                        OtDeg[tmp1]--;InDeg[tmp2]++;}
                    tmp2++;}
                Remain--;}}}
}
void NetGraphFrame::Auxiliary_RR_D_sub(int *OtLab,int child,int parent){
    int tmp1,tmp2;
    struct Regulation *rpt;
    while(1){
        tmp1=child;// Find a third node that not child & parent, and must be done-node.
        while(tmp1==child||tmp1==parent||OtLab[tmp1])tmp1=(int)(unif_rand()*total);
        tmp2=(int)(unif_rand()*OtDeg[tmp1]);// Should the find the replaced child of node [tmp1].
        rpt=&Network[tmp1];
        while(tmp2>0){rpt=rpt->next;tmp2--;}
        tmp2=rpt->next->code;
        if(tmp2!=parent){// Find it.
            DeleteParents(&Network[tmp2],tmp1);
            DeleteChilds(&Network[tmp1],tmp2);
            BuildRegulationship(Network,child,tmp1);
            BuildRegulationship(Network,tmp2,parent);
            break;}}
}
void NetGraphFrame::Auxiliary_RR_D(){
    int Re_In,Re_Ot=total,N_INT=total*sizeof(int),K=*InDeg;
    int *tmp_ind=(int *)malloc(N_INT);
    int *tmp_otd=(int *)malloc(N_INT);
    int tmpI,tmp1,tmp2;
    memcpy(tmp_ind,InDeg,N_INT);// Copy temporary degrees.
    memcpy(tmp_otd,OtDeg,N_INT);// Copy temporary degrees.
    for(tmpI=0;tmpI<total;tmpI++){// Loop for assigning each node's parents.
        if(K>=Re_Ot)break;// Not enough to set parents, the equal also.
        while(tmp_ind[tmpI]){
            tmp1=(int)(unif_rand()*(Re_Ot));// Search.
            tmp2=**(Address+tmp1);// Candidate
            if(tmp2==tmpI||ExistChild(&Network[tmp2],tmpI));// Failure! (tmp2 -> tmpI) or self-loop
            else {
                BuildRegulationship(Network,tmpI,tmp2);
                tmp_ind[tmpI]--;// Child minus one.
                tmp_otd[tmp2]--;// Parent minus one.
                if(0==tmp_otd[tmp2]){
                    Re_Ot--;
                    memcpy(Address+tmp1,Address+tmp1+1,(Re_Ot-tmp1)*sizeof(int *));}}}}
    Re_In=total-tmpI;// Some unsatisfied nodes exist.
    while(Re_In){
        while(tmp_ind[tmpI]){
            tmp1=(int)(unif_rand()*Re_Ot);// Search.
            tmp2=**(Address+tmp1);// Candidate
            if(tmp2==tmpI||ExistChild(&Network[tmp2],tmpI)){
                NetGraphFrame::Auxiliary_RR_D_sub(tmp_otd,tmpI,tmp2);}
            else {
                BuildRegulationship(Network,tmpI,tmp2);}
            tmp_ind[tmpI]--;// Child minus one.
            tmp_otd[tmp2]--;// Parent minus one.
            if(0==tmp_otd[tmp2]){
                Re_Ot--;
                memcpy(Address+tmp1,Address+tmp1+1,(Re_Ot-tmp1)*sizeof(int *));}}
        Re_In--;tmpI++;}
    free(tmp_ind);free(tmp_otd);
}

// Auxiliary function for lattice undirected.
void NetGraphFrame::Auxiliary_LT_D(int ids,int *Neighbor){// Note: ids only 3, 4, 6; total is a square number (even*even).
    int Length=tmp_i_arg[1],tmp1,tmp2,tmp3;
    int Height=total/Length;
    int ii=ids/Length,jj=ids%Length;
    if(tmp_i_arg[0]==4){// Square
        Neighbor[0]=jj-1;if(Neighbor[0]<0)Neighbor[0]=Length-1;Neighbor[0]=Neighbor[0]+ii*Length;
        Neighbor[1]=jj+1;if(Neighbor[1]==Length)Neighbor[1]=0;Neighbor[1]=Neighbor[1]+ii*Length;
        Neighbor[2]=ii-1;if(Neighbor[2]<0)Neighbor[2]=Height-1;Neighbor[2]=Neighbor[2]*Length+jj;
        Neighbor[3]=ii+1;if(Neighbor[3]==Height)Neighbor[3]=0;Neighbor[3]=Neighbor[3]*Length+jj;}
    else if(tmp_i_arg[0]==3){// Triangle
        Neighbor[0]=jj-1;if(Neighbor[0]<0)Neighbor[0]=Length-1;Neighbor[0]=Neighbor[0]+ii*Length;
        Neighbor[1]=jj+1;if(Neighbor[1]==Length)Neighbor[1]=0;Neighbor[1]=Neighbor[1]+ii*Length;
        if((ii+jj)%2==0){// even
            Neighbor[2]=ii+1;if(Neighbor[2]==Height)Neighbor[2]=0;Neighbor[2]=Neighbor[2]*Length+jj;}
        else {// odd
            Neighbor[2]=ii-1;if(Neighbor[2]<0)Neighbor[2]=Height-1;Neighbor[2]=Neighbor[2]*Length+jj;}}
    else if(tmp_i_arg[0]==6){// Hexagon
        Neighbor[0]=jj-1;if(Neighbor[0]<0)Neighbor[0]=Length-1;Neighbor[0]=Neighbor[0]+ii*Length;
        Neighbor[1]=jj+1;if(Neighbor[1]==Length)Neighbor[1]=0;Neighbor[1]=Neighbor[1]+ii*Length;
        tmp1=ii-1;tmp2=ii+1;
        if(tmp1<0){tmp1=Height-1;}if(tmp2==Height){tmp2=0;}
        if(1==(ii&1)){// Odd lines
            tmp3=jj+1;if(tmp3==Length)tmp3=0;}
        else {// Even lines
            tmp3=jj-1;if(tmp3<0)tmp3=Length-1;}
        Neighbor[2]=tmp1*Length+jj;
        Neighbor[3]=tmp1*Length+tmp3;
        Neighbor[4]=tmp2*Length+jj;
        Neighbor[5]=tmp2*Length+tmp3;
    }
}

void NetGraphFrame::Auxiliary_NetworkTarjan(int ThisID,std::vector<bool> &hold_on, std::vector<int> &finding,
    std::vector<int> &lowrank, std::vector<int> &sequnces, std::vector<std::vector<int>> &scc, int &step, int &index){
    finding[ThisID]=lowrank[ThisID]=step;
    step++;
    sequnces.push_back(ThisID);
    hold_on[ThisID]=true;
    int AdjNode;
    struct Regulation *tmp=Network+ThisID;// Get node's address;
    while(nullptr!=tmp->next){
        AdjNode=tmp->next->code;
        if(finding[AdjNode]<0){
            NetGraphFrame::Auxiliary_NetworkTarjan(AdjNode, hold_on,
                finding, lowrank, sequnces, scc, step, index);
            lowrank[ThisID]=(lowrank[ThisID]<lowrank[AdjNode])?lowrank[ThisID]:lowrank[AdjNode];
        }
        else if(hold_on[AdjNode]){
            lowrank[ThisID]=(lowrank[ThisID]<finding[AdjNode])?lowrank[ThisID]:finding[AdjNode];
        }
        tmp=tmp->next;
    }
    if(finding[ThisID]==lowrank[ThisID]){// In the SCC
        //scc[index]
        std::vector<int> tmp_scc;
        int tmp_id;
        while(true){
            tmp_id=sequnces.back();
            sequnces.pop_back();
            hold_on[tmp_id]=false;
            tmp_scc.push_back(tmp_id);
            if(ThisID==tmp_id)break;
        }
        scc.push_back(tmp_scc);
        index++;
    }
}

// Return the list of strong connected components of a network.
NetGraphFrame& NetGraphFrame::Tarjon(std::vector<std::vector<int>> &SCC_List){
    SCC_List.clear();
    std::vector<bool> HoldOn(total,false);
    std::vector<int> finding(total,-666);
    std::vector<int> lowrank(total,-666);
    std::vector<int> Sequen;Sequen.reserve(total);
    int Timer=0, scc_id=0; 
    //struct Regulation *tmp=nullptr;
    for(int ii=0; ii<total; ++ii){
        if(finding[ii]<0){
            NetGraphFrame::Auxiliary_NetworkTarjan(ii, HoldOn,
                finding, lowrank, Sequen, SCC_List, Timer, scc_id);}}
    return *this;
}
// Code is over.
