#pragma once
#include "Defines.h"
using namespace std;

extern int k,bestSz,*bestSol, maxsec, optimal;
extern long long int branches1,branches2;
extern int rule1, rule2, rule3, rule4;
#include "MVarSet.h"
#include "MStack.h"
#include "MGraph.h"

class HeuriSearcher{
public:
	MCSRGraph &g;
	int k;
	int *core;
	int *seq;
	int *pos;
	int* sol;
	int solSz;

	HeuriSearcher(MCSRGraph &_g, int _k,	int *_core, int *_seq, int *_pos,int *_sol)
	:g(_g), k(_k), core(_core), seq(_seq), pos(_pos),sol(_sol),solSz(0)
	{
	}

	int coreHeuristic(){    
		int *deg = new int[g.nbvtx];
		MBitSet64 *isSol = new MBitSet64(g.nbvtx);
		int maxcore = 0;
		bool hit = false;

		for (int i = 0; i < g.nbvtx; i++) {		
			deg[i] = g.degree(i);
			seq[i] = i;
			core[i] = -1;
		}
		ListLinearHeap *linear_heap = new ListLinearHeap(g.nbvtx, g.nbvtx);
		linear_heap->init(g.nbvtx, g.nbvtx, seq, deg);
		for (int i = 0; i < g.nbvtx; i++) {
			int u, key;
			linear_heap->pop_min(u, key);
			if (key > maxcore) 
				maxcore = key;	
			seq[i] = u;
			core[u] = maxcore;
			pos[u] = i;
			if ( i + key + k >= g.nbvtx && !hit){//这一步为什么??
				int sz = i+1;
				linear_heap->get_ids_of_larger_keys(seq, sz, key);
				solSz=g.nbvtx-i;
				int cnt=0;							
				for (int j = i ; j < g.nbvtx; j++){//vertices after i is not ordered，but the key is larger or equal than i
					sol[cnt++] = seq[j];				
					isSol->set(seq[j]);
				}
				hit = true;
			}
			for (int j = g.pstart[u]; j < g.pstart[u + 1]; j++){
				if (core[g.edges[j]] == -1)//重置u的邻居的core值
					linear_heap->decrement(g.edges[j]);//减少u的第j个邻居的度数，因此每个点的key是在不断变化的，整个过程就是求degenracy order的过程
			}
		}
		//extend the heuristic solution.
		//compute degrees
		memset(deg, 0, sizeof(int) * g.nbvtx);	
		int totSat = 0;
		const int lim=solSz-k;
		for (int i = 0; i < solSz; i++){
			const int u = sol[i];
			for (int j = g.pstart[u]; j < g.pstart[u+1]; j++){
				const int nei = g.edges[j];
				if (isSol->test(nei)){
					deg[u]++;
				}
			}
			if (deg[u] == lim){
				totSat++;//记录饱和点数量
			}
		}
		
		//extending the remaining vertices
		for (int i = g.nbvtx - 1 - solSz; i >= 0; i--){
			const int lim=solSz-k;
			const int u = seq[i];
			if (core[u] + k <= solSz) break;
			int cntNei = 0, cntSat = 0;
			for (int j = g.pstart[u]; j < g.pstart[u+1]; j++){
				int nei = g.edges[j];
				if (isSol->test(nei)){
					cntNei++;
					if (deg[nei] == lim){
						cntSat++;
					}
				}
			}
			//extend u
			if (cntNei > lim && cntSat == totSat){ 
				deg[u] = cntNei;
				sol[solSz++] = u;
				isSol->set(u);
				//update number of saturated vertices.			
				for (int j = g.pstart[u]; j < g.pstart[u+1]; j++){
					const int nei = g.edges[j];
					if (isSol->test(nei)){
						deg[nei]++;
					}
				}
				//recount satu vertices
				totSat = 0;
				for (int j = 0; j < solSz; j++){				
					if (deg[sol[j]] == lim+1){
						totSat++;
					}
				}
			}
		}
		delete linear_heap;
		delete isSol;
		delete[] deg;
		return maxcore;
	}
};

using PlexStack = Stack<int>;
using VtxSet = VarSet<int>;
struct ExactSearcher
{
	const MCSRGraph &subg;
	MCSRGraph &cog;
	int *neiInP;
	int *neiInG;
	bool found;
	bool stopable;
	PlexStack *plex;
	VtxSet *cand1;
	VtxSet *cand2;
	VtxSet *candT;
	VtxSet *candTCo;
	const uint8_t* const adjMtx;
	int* const commMtx;
	int depth;
	int depthLim;
	clock_t startClk;
	ExactSearcher(const MCSRGraph &_subg, MCSRGraph &_cog, PlexStack *_plex, VtxSet *_cand1, VtxSet *_cand2, VtxSet *_candT, VtxSet *_candTCo,int *_neiInP, int *_neiInG, uint8_t* _adjMtx, int* _commMtx, clock_t _startClk)
	:subg(_subg),cog(_cog),plex(_plex),cand1(_cand1),cand2(_cand2),candT(_candT),candTCo(_candTCo),depth(0),depthLim(k/2),found(false),stopable(true)
	,neiInP(_neiInP),neiInG(_neiInG),adjMtx(_adjMtx),commMtx(_commMtx),startClk(_startClk)
	{
	}

	int interrupt() {
		if ((double)(clock() - startClk) / CLOCKS_PER_SEC > maxsec) {
			optimal = 0;
			return 1;
		}
		return 0;
	}

	bool isAdjMtx(const int v1,const int v2)
	{
		return adjMtx[v1*subg.nbvtx+v2];
	}

	bool isAdjMtxCo(const int v1,const int v2)
	{
		return !isAdjMtx(v1,v2);
	}

	void addG(const int u){
		for (int i = subg.pstart[u]; i < subg.pstart[u+1]; i++)
		{
			const int nei = subg.edges[i];
			neiInG[nei]++;
			if(stopable&&depth<depthLim){
				for (int j = i+1; j < subg.pstart[u+1]; j++)
				{
					const int nei2 = subg.edges[j];
					commMtx[nei*subg.nbvtx+nei2]++;
					commMtx[nei2*subg.nbvtx+nei]++;
				}
			}
		}
	}

	void subG(const int u){
		for (int i = subg.pstart[u]; i < subg.pstart[u+1]; i++)
		{
			const int nei = subg.edges[i];
			neiInG[nei]--;
			if(stopable&&depth<depthLim){
				for (int j = i+1; j < subg.pstart[u+1]; j++)
				{
					const int nei2 = subg.edges[j];
					commMtx[nei*subg.nbvtx+nei2]--;
					commMtx[nei2*subg.nbvtx+nei]--;
				}
			}
		}
	}

	void addGCo(const int u){
		for (int i = cog.pstart[u]; i < cog.pstart[u+1]; i++)
		{
			const int nei = cog.edges[i];
			neiInG[nei]++;
		}
	}

	void subGCo(const int u){
		for (int i = cog.pstart[u]; i < cog.pstart[u+1]; i++)
		{
			const int nei = cog.edges[i];
			neiInG[nei]--;
		}
	}

	void addP(const int u){
		for (int i = subg.pstart[u]; i < subg.pstart[u+1]; i++)
		{
			const int nei = subg.edges[i];
			neiInP[nei]++;
		}
	} 

	void subP(const int u){
		for (int i = subg.pstart[u]; i < subg.pstart[u+1]; i++)
		{
			const int nei = subg.edges[i];
			neiInP[nei]--;
		}
	}

	void addPCo(const int u){
		for (int i = cog.pstart[u]; i < cog.pstart[u+1]; i++)
		{
			const int nei = cog.edges[i];
			neiInP[nei]++;
		}
	} 

	void subPCo(const int u){
		for (int i = cog.pstart[u]; i < cog.pstart[u+1]; i++)
		{
			const int nei = cog.edges[i];
			neiInP[nei]--;
		}
	}

	ExactSearcher* cand1ToPlex(const int v)
	{
		plex->push(v);
		cand1->remove(v);
		addP(v);
		return this;
	}

	ExactSearcher* cand1ToPlexCo(const int v)
	{
		plex->push(v);
		cand1->remove(v);
		addPCo(v);
		return this;
	}

	ExactSearcher* plexToCand1()
	{
		assert(plex->sz > 0);
		const int u = plex->top();
		cand1->add(u);
		plex->pop();
		subP(u);
		return this;
	}

	ExactSearcher* plexToCand1Co()
	{
		assert(plex->sz > 0);
		const int u = plex->top();
		cand1->add(u);
		plex->pop();
		subPCo(u);
		return this;
	}

	ExactSearcher* cand2ToPlex(const int v)
	{
		cand2->remove(v);
		plex->push(v);
		addP(v);
		return this;
	}

	ExactSearcher* plexToCand2()
	{
		assert(plex->sz > 0);
		const int v = plex->top();
		cand2->add(v);
		plex->pop();
		subP(v);
		return this;
	}

	ExactSearcher* cand1Move(const int v)
	{
		cand1->remove(v);
		subG(v);
		return this;
	}

	ExactSearcher* cand1Add(const int v)
	{
		cand1->add(v);
		addG(v);
		return this;
	}

	ExactSearcher* cand1MoveCo(const int v)
	{
		assert(cand1->contains(v));
		cand1->remove(v);
		subGCo(v);
		return this;
	}

	ExactSearcher* cand1AddCo(const int v)
	{
		assert(!cand1->contains(v));
		cand1->add(v);
		addGCo(v);
		return this;
	}

	ExactSearcher* cand2Move(const int v)
	{
		cand2->remove(v);
		subG(v);
		return this;
	}

	ExactSearcher* cand2Add(const int v)
	{
		cand2->add(v);
		addG(v);
		return this;
	}

	ExactSearcher* updateCand1(int& recCand1){
		int rec = cand1->sz;
		for(int i = 0; i < cand1->sz;)
		{
			int ele = cand1->members[i];
			if((!canFormPlex(ele))){
				cand1->remove(ele);
				subG(ele);
			}
			else ++i;
		}
		rec -= cand1->sz;
		recCand1 += rec;
		return this;
	}

	ExactSearcher* updateCand1Co(int& recCand1){
		int rec = cand1->sz;
		for(int i = 0; i < cand1->sz;)
		{
			int ele = cand1->members[i];
			if((!canFormPlexCo(ele))){
				cand1->remove(ele);
				subGCo(ele);
			}
			else ++i;
		}
		rec -= cand1->sz;
		recCand1 += rec;
		return this;
	}

	ExactSearcher* updateCand2(int& recCand2){
		int rec = cand2->sz;
		for(int i = 0; i < cand2->sz;)
		{
			int ele = cand2->members[i];
			if((!canFormPlex(ele))){
				cand2->remove(ele);
				subG(ele);
			}
			else ++i;
		}
		rec -= cand2->sz;
		recCand2 += rec;
		return this;
	}

	ExactSearcher* recoverCand1(const int recCand1){
        int* cursor=cand1->members+cand1->sz;
        for(int i=0;i<recCand1;++i){
            const int u=*cursor;
            addG(u);
            cursor++;
        }
		cand1->recover(recCand1);
		return this;
	}

	ExactSearcher* recoverCand1Co(const int recCand1){
        int* cursor=cand1->members+cand1->sz;
        for(int i=0;i<recCand1;++i){
            const int u=*cursor;
            addGCo(u);
            cursor++;
        }
		cand1->recover(recCand1);
		return this;
	}

	ExactSearcher* recoverCand2(const int recCand2){
        int* cursor=cand2->members+cand2->sz;
        for(int i=0;i<recCand2;++i){
            const int u=*cursor;
            addG(u);
            cursor++;
        }
		cand2->recover(recCand2);
		return this;
	}

	//Check if plex+u is a k-plex
	bool canFormPlex(const int u)
	{
		if (neiInP[u]+k<plex->sz+1||neiInG[u]+k<=bestSz) return false;
		for (int i = 0; i < plex->sz; i++)
		{
			const int v = plex->members[i];
			if (stopable && commMtx[u*subg.nbvtx+v] < bestSz+1-2*k+2-2*isAdjMtx(v, u)) return false;
			if (neiInP[v]+k==plex->sz&&!isAdjMtx(v, u)) return false;
		}
		return true;
	}

	//Check if plex+u is a k-1-plexco
	bool canFormPlexCo(const int u)
	{
		if (neiInP[u] >= k || plex->sz+cand1->sz+k-1-neiInG[u]<=bestSz) return false;
		for (int i = 0; i < plex->sz; i++)
		{
			const int v = plex->members[i];
			if (neiInP[v] == k-1 && isAdjMtxCo(v, u))
			{
				return false;
			}
		}
		return true;
	}

	bool checkNei(int u) {
		int notSat = 0;
		for(int i = 0; i < cand1->sz; ++i) {
			int ele = cand1->members[i];
			if(!isAdjMtx(u, ele) && u != ele) {
				if(subg.nbvtx-neiInG[ele] > k) {
					notSat++;
					if(notSat > 1)
						return false;
				}
			}
		}
		for(int i = 0; i < cand2->sz; ++i) {
			int ele = cand2->members[i];
			if(!isAdjMtx(u, ele) && u != ele) {
				if(subg.nbvtx-neiInG[ele] > k) {
					notSat++;
					if(notSat > 1)
						return false;
				}
			}
		}
		return true;
	}

	bool checkNeiCo(int u) {
		int notSat = 0;
		for(int i = 0; i < cand1->sz; ++i) {
			int ele = cand1->members[i];
			if(isAdjMtxCo(u, ele) && u != ele) {
				if(neiInG[ele] >= k) {
					notSat++;
					if(notSat > 1)
						return false;
				}
			}
		}
		return true;
	}

	bool checkTrangle(int u) {
		vector<int> NotNei;
		for(int i = 0; i < cand1->sz; ++i) {
			int ele = cand1->members[i];
			if(!isAdjMtx(u, ele) && u != ele) 
				NotNei.emplace_back(ele);
		}
		for(int i = 0; i < cand2->sz; ++i) {
			int ele = cand2->members[i];
			if(!isAdjMtx(u, ele) && u != ele) 
				NotNei.emplace_back(ele);
		}
		if(NotNei.size() == 2) {
			int ele1 = NotNei[0], ele2 = NotNei[1];
			if(!isAdjMtx(ele1, ele2)) 
				return true;
		}
		return false;
	}

	bool checkTrangleCo(int u) {
		vector<int> nei;
		for(int i = 0; i < cand1->sz; ++i) {
			int ele = cand1->members[i];
			if(isAdjMtxCo(u, ele) && u != ele) 
				nei.emplace_back(ele);
		}
		if(nei.size() == 2) {
			int ele1 = nei[0], ele2 = nei[1];
			if(isAdjMtxCo(ele1, ele2)) 
				return true;
		}
		return false;
	}

	bool checkPlex(int u) {
		for(int i = 0; i < plex->sz; i++) {
			int ele = plex->members[i];
			if(!isAdjMtx(u, ele)) {
				if(plex->sz-neiInP[ele] < k)
					return true;
				else
					return false;
			}
		}
		return false;
	}

	bool checkPlexCo(int u) {
		for(int i = 0; i < plex->sz; i++) {
			int ele = plex->members[i];
			if(!isAdjMtx(u, ele)) {
				if(neiInP[ele] < k-1)
					return true;
				else
					return false;
			}
		}
		return false;
	}

	bool bound() {
		candT->clear();
		for(int i = 0; i < cand1->sz; i++) candT->add(cand1->members[i]);
		for(int i = 0; i < cand2->sz; i++) candT->add(cand2->members[i]);
		auto cmp= [](pair<int, int> a, pair<int, int> b) {
			if(a.second == b.second)
				return a.first < b.first;
			return a.second < b.second;
		};
		int UB = plex->sz;
		vector<pair<int, int> > ord_P;
		for(int i = 0; i < plex->sz; i++) {
			int u = plex->members[i];
			ord_P.push_back(make_pair(u, neiInG[u]));
		}
		sort(ord_P.begin(), ord_P.end(), cmp);
		for(auto PP : ord_P) {
			int u = PP.first;
			int NonAdj = 0;
			for(int j = candT->sz-1; j >= 0; j--) {
				int v = candT->members[j];
				if(!isAdjMtx(u, v)) {
					candT->remove(v);
					NonAdj++;
				}
			}
			UB += min(NonAdj, neiInP[u]+k-plex->sz);
			if(UB > bestSz) return true;
		}
		if(UB + candT->sz > bestSz) return true;
		return false;
	}

	bool boundCo() {
		candTCo->clear();
		for(int i = 0; i < cand1->sz; i++) candTCo->add(cand1->members[i]);
		auto cmp= [](pair<int, int> a, pair<int, int> b) {
			if(a.second == b.second)
				return a.first > b.first;
			return a.second > b.second;
		};
		int UB = plex->sz;
		vector<pair<int, int> > ord_P;
		for(int i = 0; i < plex->sz; i++) {
			int u = plex->members[i];
			ord_P.push_back(make_pair(u, neiInG[u]));
		}
		sort(ord_P.begin(), ord_P.end(), cmp);
		for(auto PP : ord_P) {
			int u = PP.first;
			int adj = 0;
			for(int j = candTCo->sz-1; j >= 0; j--) {
				int v = candTCo->members[j];
				if(isAdjMtxCo(u, v)) {
					candTCo->remove(v);
					adj++;
				}
			}
			UB += min(adj, k-1-neiInP[u]);
			if(UB > bestSz) return true;
		}
		if(UB + candTCo->sz > bestSz) return true;
		return false;
	}

	void reduce(int &popcnt1, int &popcnt2, int &recReCand1, int &recReCand2) {
		bool flag;
		do{
			flag=false;
			for(int i = 0; i < cand1->sz; ) {
				const int ele = cand1->members[i];
				const int NotNei = cand1->sz+cand2->sz+plex->sz-neiInG[ele];
				if(neiInP[ele] == plex->sz && NotNei <= 2) {
					cand1ToPlex(ele);
					popcnt1++;
					rule1++;
					flag=true;
				}
				else if(neiInP[ele] == plex->sz && NotNei <= k+1 && checkNei(ele)) {
					cand1ToPlex(ele);
					popcnt1++;
					rule2++;
					flag=true;
				}
				else if(neiInP[ele] == plex->sz && NotNei == 3 && checkTrangle(ele)) {
					cand1ToPlex(ele);
					popcnt1++;
					rule3++;
					flag=true;
				}
				else if(neiInP[ele] == plex->sz-1 && NotNei == 2 && checkPlex(ele)) {
					cand1ToPlex(ele);
					popcnt1++;
					rule4++;
					flag=true;
				}
				else
					++i;
			}

			for(int i = 0; i < cand2->sz; ) {
				const int ele = cand2->members[i];
				const int NotNei = cand1->sz+cand2->sz+plex->sz-neiInG[ele];
				if(neiInP[ele] == plex->sz-1 && NotNei == 2 && checkPlex(ele)) {
					cand2ToPlex(ele);
					popcnt2++;
					rule4++;
					flag=true;
				}
				else
					++i;
			}
		}while(flag);
		updateCand1(recReCand1);
		updateCand2(recReCand2);
	}

	void reduceCo(int &popcnt1, int &recReCand1) {
		bool flag;
		do{
			flag=false;
			for(int i = 0; i < cand1->sz; ) {
				const int ele = cand1->members[i];
				const int nei = neiInG[ele];
				if(neiInP[ele] == 0 && nei <= 1) {
					cand1ToPlexCo(ele);
					popcnt1++;
					rule1++;
					flag=true;
				}
				else if(neiInP[ele] == 0 && nei <= k && checkNeiCo(ele)) {
					cand1ToPlexCo(ele);
					popcnt1++;
					rule2++;
					flag=true;
				}
				else if(neiInP[ele] == 0 && nei == 2 && checkTrangleCo(ele)) {
					cand1ToPlexCo(ele);
					popcnt1++;
					rule3++;
					flag=true;
				}
				else if(neiInP[ele] == 1 && nei == 1 && checkPlexCo(ele)) {
					cand1ToPlexCo(ele);
					popcnt1++;
					rule4++;
					flag=true;
				}
				else
					++i;
			}
		}
		while(flag);
		updateCand1Co(recReCand1);
	}

	void recover(int &popcnt1, int &popcnt2, int &recReCand1, int &recReCand2) {
		recoverCand2(recReCand2);
		recoverCand1(recReCand1);
		while(popcnt2--) plexToCand2();
		while(popcnt1--) plexToCand1();
	}

	void recoverCo(int &popcnt1, int &recReCand1) {
		recoverCand1Co(recReCand1);
		while(popcnt1--) plexToCand1Co();
	}
	
	void hop2Search(int res){
		if(interrupt()) return ;
		if(stopable&&found) return;
		branches2++;

		//reduce
		int popcnt1 = 0, popcnt2 = 0;
		int recReCand1 = 0, recReCand2 = 0;
		int recCand1=0; int recCand2=0;
		reduce(popcnt1, popcnt2, recReCand1, recReCand2);
		res -= popcnt2;

		//bound
		int minnei=INT_MAX-k; int pivot;
		if(plex->sz+cand1->sz+std::min(res,cand2->sz)<=bestSz) goto RECOVER2;
		for(int i=0;i<plex->sz;++i){
			const int ele=plex->members[i];
			const int nei=neiInG[ele];
			if(nei<minnei) minnei=nei;
		}
		if (minnei + k <= bestSz) goto RECOVER2;
		if(!cand2->sz){depth++,hop1Search(),depth--;goto RECOVER2;}
		for(int i=0;i<cand2->sz;++i){
			const int ele=cand2->members[i];
			const int nei=neiInG[ele];
			if(nei<minnei) minnei=nei;
		}
		for(int i=0;i<cand1->sz;++i){
			const int ele=cand1->members[i];
			const int nei=neiInG[ele];
			if(nei<minnei) minnei=nei;
		}
		if (minnei + k >= plex->sz + cand1->sz + cand2->sz)
		{
			bestSz = plex->sz + cand1->sz + cand2->sz;
			found=true;
			memcpy(bestSol,plex->members,plex->sz*sizeof(int));
			memcpy(&bestSol[plex->sz],cand1->members,cand1->sz*sizeof(int));
			memcpy(&bestSol[plex->sz+cand1->sz],cand2->members,cand2->sz*sizeof(int));
			cout<<bestSz<<endl;
			goto RECOVER2;
		}
		if(!bound()) goto RECOVER2;
		//branch
		pivot=cand2->members[cand2->sz-1];
		cand2ToPlex(pivot);
		updateCand1(recCand1),updateCand2(recCand2);
		depth++,hop2Search(res-1),depth--;
		recoverCand1(recCand1);
		recoverCand2(recCand2);
		plexToCand2();
		cand2Move(pivot);
		depth++,hop2Search(res),depth--;
		cand2Add(pivot);

	RECOVER2:
		recover(popcnt1, popcnt2, recReCand1, recReCand2);
	}

	void hop1Search(){
		candT->clear();
		for(int i=0;i<cand1->sz;++i) candT->add(cand1->members[i]);
		for(int i=0;i<plex->sz;++i) candT->add(plex->members[i]);
		int nbedges=subg.nbvtx*(subg.nbvtx-1)-subg.nbedges;
		cog.edges = new int[nbedges];
		cog.nbvtx=subg.nbvtx;
		cog.nbedges=0;
		cog.pstart = new int[cog.nbvtx+1];
		cog.pstart[0] = 0;
		for(int i=0;i<cog.nbvtx;++i){
			cog.pstart[i+1] = cog.pstart[i];
			if(!candT->contains(i))continue;
			for(int j=0; j<cog.nbvtx;++j ){
				if( j!=i && candT->contains(j) && isAdjMtxCo(i,j)){
					cog.edges[cog.pstart[i+1]++] = j;
					cog.nbedges++;
				}
			}
		}
		// cout<<(plex->sz+cand1->sz)<<" "<<(plex->sz+cand1->sz)-(bestSz+1)<<endl;
		for(int i=0;i<plex->sz;++i){
			neiInP[plex->members[i]]=plex->sz-1-neiInP[plex->members[i]];
			neiInG[plex->members[i]]=candT->sz-1-neiInG[plex->members[i]];
		}
		for(int i=0;i<cand1->sz;++i){
			neiInP[cand1->members[i]]=plex->sz-neiInP[cand1->members[i]];
			neiInG[cand1->members[i]]=candT->sz-1-neiInG[cand1->members[i]];
		}
		hop1SearchCo();
		for(int i=0;i<plex->sz;++i){
			neiInP[plex->members[i]]=plex->sz-1-neiInP[plex->members[i]];
			neiInG[plex->members[i]]=candT->sz-1-neiInG[plex->members[i]];
		}
		for(int i=0;i<cand1->sz;++i){
			neiInP[cand1->members[i]]=plex->sz-neiInP[cand1->members[i]];
			neiInG[cand1->members[i]]=candT->sz-1-neiInG[cand1->members[i]];
		}
		cog.del();
	}

	void hop1SearchCo(){
		if(interrupt())return;
		if(stopable && found) return;
		branches1++;

		//reduce
		int popcnt1 = 0, recReCand1 = 0;
		reduceCo(popcnt1, recReCand1);

		//bound
		int maxnei=INT_MIN;int maxneiPlex;int pivot;
		auto maxswap = [&](int u) {
			if (neiInG[u] > maxnei)
			{
				maxnei = neiInG[u];
				pivot = u;
			}
		};
		if (plex->sz + cand1->sz <= bestSz) goto RECOVER1;
		if (cand1->sz == 0)
		{
			bestSz=plex->sz + cand1->sz;
			found=true;
			memcpy(bestSol,plex->members,plex->sz*sizeof(int));
			cout<<bestSz<<endl;
			goto RECOVER1;
		}
		plex->for_each(maxswap);
		if (plex->sz + cand1->sz-max(0,(maxnei-(k-1))) <= bestSz) goto RECOVER1; 
		maxneiPlex=maxnei;
		cand1->for_each(maxswap);
		if (maxnei < k)
		{
			bestSz = plex->sz + cand1->sz;
			found=true;
			memcpy(bestSol,plex->members,plex->sz*sizeof(int));
			memcpy(bestSol+plex->sz,cand1->members,cand1->sz*sizeof(int));
			cout<<bestSz<<endl;
			goto RECOVER1;
		}
		if(!boundCo()) goto RECOVER1;
		
		//branch
		{
			int addList[K_LIMIT];
			int recCand1[K_LIMIT];
			memset(recCand1,0,sizeof recCand1);
			if(maxnei!=maxneiPlex){
				cand1MoveCo(pivot);
				hop1SearchCo();
				cand1AddCo(pivot);
				cand1ToPlexCo(pivot)->updateCand1Co(recCand1[0]);
				hop1SearchCo();
				recoverCand1Co(recCand1[0])->plexToCand1Co();
				goto RECOVER1;
			}

			int tot = k-1 - neiInP[pivot];
			int idx=0;
			for(int i=0;i<cand1->sz;++i){
				int ele=cand1->members[i];
				if(isAdjMtxCo(pivot,ele)) addList[idx++]=ele;
				if(idx==tot)break;
			}
			//There are at most tol+1 branches.
			//br0;
			int v2del = addList[0];
			cand1MoveCo(v2del)->hop1SearchCo();
			cand1AddCo(v2del);

			//Proceed br1 in [1,tol-1); there are at least 2 branches
			//add in addList[0,...,br), remove addList[br].
			int v2add;
			int movs = 1; // recording the number of vertices pushed

			//start from the 1st br.
			int br = 1;
			for (; br < tot; br++)
			{
				v2add = addList[br - 1];
				v2del = addList[br];
				if (cand1->contains(v2add) && plex->sz + cand1->sz > bestSz) 
				{
					cand1ToPlexCo(v2add)->updateCand1Co(recCand1[movs++]);
					if (cand1->contains(v2del))
					{
						cand1MoveCo(v2del);
						hop1SearchCo();
						cand1AddCo(v2del);
					}
					else hop1SearchCo();
				}
				else break;
			}
			//The last (tol-th) branch.
			if (br == tot && cand1->contains(addList[br - 1]) && plex->sz + cand1->sz > bestSz)
			{
				v2add = addList[tot - 1];
				cand1ToPlexCo(v2add)->updateCand1Co(recCand1[movs++]); //all vertices in addList has moved to plex        
				hop1SearchCo();
			}
			//recover
			movs--;
			while(movs){
				recoverCand1Co(recCand1[movs]);
				plexToCand1Co();
				movs--;
			}
		}

	RECOVER1:
		recoverCo(popcnt1, recReCand1);
	}
};