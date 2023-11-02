// #include "Defines.h"
// #include "MGraph.h"
// #include "LinearHeap.h"
// //grpah algorithm: core decomposition
// int coreDecomposition(const MCSRGraph &g, int *core, int* seq, int* pos){//这里core数组应该不是递增的??
//     assert(core!= nullptr && seq!= nullptr && pos!= nullptr);

//     int *degree = new int[g.nbvtx];	
// 	char *isout = new char[g.nbvtx];
// 	for (int i = 0; i < g.nbvtx; i++) {		
// 		degree[i] = g.degree(i);		
//         seq[i] = i;
// 		isout[i] = 0;
// 	}

// 	int max_core = 0;
// 	ListLinearHeap *linear_heap = new ListLinearHeap(g.nbvtx, g.nbvtx);
// 	linear_heap->init(g.nbvtx, g.nbvtx, seq, (int*)degree);	
// 	for (int i = 0; i < g.nbvtx; i++) {
// 		int u, key;
// 		linear_heap->pop_min(u, key);
// 		if (key > max_core)
// 			max_core = key;		
// 		seq[i] = u;
// 		core[u] = max_core;
//         pos[u] = i;
// 		isout[u] = 1;
// 		for (int j = g.pstart[u]; j < g.pstart[u + 1]; j++)
// 			if (!isout[g.edges[j]])
// 				linear_heap->decrement(g.edges[j]);
// 	}	
//     delete[] degree;
//     delete[] isout;
//     delete linear_heap;
//     return max_core;
// }

// using namespace std;
// //The number of vertices will not change
// int *cg_pstart;
// int *cg_pend;
// int *cg_edges; //unordered
// int *cg_valid;
// int cg_nbvtx;
// int cg_nbedges;
// int cg_maxvtx;
// int cg_maxedge;
// int cg_MAXNBEDGES;
// //TODO:keep vertice 
// map<int, vector<int> > cg_vmp;

// inline int cg_deg(int u){
// 	return cg_pend[u] - cg_pstart[u];
// }

// /** Counting triangles for each edge cg_edge[j], j=cg_pstart[u] + i
//  *  traingle[j] is avalaible only if u < cg_edges[j];
//  */
// void cg_countingTriangles(int *triangles){
// 	//  triangle counting
// 	for (int u = 0; u < cg_maxvtx; u++){
// 		if (cg_valid[u]){
// 			for (int j = cg_pstart[u]; j < cg_pend[u]; j++){
// 				int neiu = cg_edges[j];
// 				if (u < neiu){
// 					triangles[j] = Utility::countIntersect(cg_edges + cg_pstart[u], cg_deg(u), cg_edges + cg_pstart[neiu], cg_deg(neiu));//计算以边j为一条边的三角形个数
// 				}
// 			}
// 		}
// 	}
// }

// /** Remove vertex u from cg graph.
//  * ensure the neighbors of u are sorted in ascending order.
//  * */
// void cg_removeVtx(int u){
// 	assert(cg_valid[u] && u < cg_maxvtx);
// 	for (int i = cg_pstart[u]; i < cg_pend[u]; i++){
// 		int neiu = cg_edges[i];
// 		//Find u from neighbors of neiu.
// 		int posu = distance(cg_edges, find(cg_edges + cg_pstart[neiu], cg_edges + cg_pend[neiu], u));//这里用distance函数找到u在neiu邻居中的位置
// 		//remove u from N(neiu)
// 		for (int p = posu; p < cg_pend[neiu] - 1; p++) cg_edges[p] = cg_edges[p+1];
// 		cg_pend[neiu]--;		
// 	}
// 	cg_nbedges -= 2*(cg_pend[u] - cg_pstart[u]);
// 	cg_valid[u] = 0;	//mark u as removed
// 	cg_nbvtx--;
// 	cg_pend[u] = cg_pstart[u]; //TODO:remove
// }

// /** Remove edge (u, v) from cg
//  */
// void cg_removeEdge(int u, int v){
// 	assert(cg_valid[u] && cg_valid[v]);
// 	int posv = distance(cg_edges, find(cg_edges + cg_pstart[u], cg_edges + cg_pend[u], v));
// 	for (int p = posv; p < cg_pend[u] - 1; p++) cg_edges[p] = cg_edges[p+1];
// 	cg_pend[u]--;
// 	int posu = distance(cg_edges, find(cg_edges + cg_pstart[v], cg_edges + cg_pend[v], u));
// 	for (int p = posu; p < cg_pend[v] - 1; p++) cg_edges[p] = cg_edges[p+1];	
// 	cg_pend[v]--;
// 	cg_nbedges -= 2;
// }

// /**
//  * make a copy of g into current graph : cg
//  * In cg: the neigbors of u is located from cg_edges[cg_pstart[u]] to cg_edges[cg_pend[u]]. 
//  * These vertices in [cg_pstart[u]] to cg_edges[cg_pend[u]] sorted in ascending order, all of them are valid.
//  * 
//  * cg_nbvtx: is the number of vertex in cg. cg_maxvtx is the maximum vertex id in cg. A vertex id x is a valid vertex in cg 
//  * 	if only if cg_valid[x]==true.
//  */
// void cg_prepareTempGraph(const MCSRGraph &g){
// 	cg_pstart = new int[g.nbvtx+1];
// 	cg_pend = new int[g.nbvtx+1];

// 	cg_MAXNBEDGES = g.nbedges;
// 	cg_edges = new int[cg_MAXNBEDGES]; // extending the list
// 	cg_valid = new int[g.nbvtx];
// 	cg_nbvtx = g.nbvtx;
// 	cg_maxvtx = g.nbvtx;// change
// 	cg_nbedges = g.nbedges;
// 	cg_maxedge = g.nbedges;
// 	memcpy(cg_edges, g.edges, sizeof(int) * g.nbedges);
// 	for (int u= 0; u < g.nbvtx; u++){
// 		cg_vmp[u].push_back(u);
// 		cg_pstart[u] = g.pstart[u];
// 		cg_pend[u] = g.pstart[u+1];
// 		cg_valid[u] = 1;
// 		sort(cg_edges + cg_pstart[u], cg_edges + cg_pend[u]);
// 	}
// }

// void cg_back2CSRGraph(MCSRGraph &g, int lb, int* new2ori){
// 	if(cg_maxvtx <= lb || cg_nbvtx <= lb) {
// 		cg_maxedge = 0;
// 		cg_maxvtx = 0;
// 		cg_nbvtx = 0;
// 		cg_nbedges = 0;
// 	}
// 	int *newid = new int[cg_maxvtx];
// 	int *orid = new int[cg_nbvtx];
// 	//reorder the vertices
// 	int cnt = 0;
// 	for (int u = 0; u < cg_maxvtx; u++){
// 		if (cg_valid[u]) {
// 			newid[u] = cnt;
// 			orid[cnt] = new2ori[u];
// 			++cnt;
// 		}else{
// 			newid[u] = cg_nbvtx;
// 		}
// 	}
// 	g.nbvtx = cg_nbvtx;
// 	g.nbedges = cg_nbedges;
// 	g.pstart = new int[g.nbvtx + 1];
// 	g.edges = new int[g.nbedges];
// 	g.pstart[0] = 0;
// 	int ecnt = 0;
// 	for (int u = 0; u < cg_maxvtx; u++)	{
// 		if (cg_valid[u]){
// 			int neu = newid[u];
// 			g.pstart[neu+1] = g.pstart[neu];
// 			for (int j = cg_pstart[u]; j < cg_pend[u]; j++){
// 				assert(newid[cg_edges[j]] != cg_nbvtx);
// 				g.edges[g.pstart[neu+1]++] = newid[cg_edges[j]];
// 				ecnt++;
// 			}
// 			assert(cg_deg(u) == g.degree(neu));
// 		}
// 	}
// 	assert(g.pstart[g.nbvtx] == ecnt);
// 	assert(g.nbvtx == cg_nbvtx);
// 	memcpy(new2ori,orid,sizeof(int)*g.nbvtx);
// 	delete[] newid;
// 	delete[] orid;
// }

// void cg_edgeReduction(int k, int lb){
// 	const int lim=lb-k;
// 	int *triangles = new int[cg_MAXNBEDGES];
// 	set<int> Q1;
// 	bool first=true;
// 	while(true){
// 		if(!Q1.empty()){						
// 			int u = *(Q1.begin());			
// 			Q1.erase(Q1.begin());
// 			cg_removeVtx(u);
// 			int st = cg_pstart[u], end=cg_pend[u];
// 			for (int i =  st; i < end; i++){
// 				int v = cg_edges[i];								
// 				if (Q1.find(v) == Q1.end() && cg_deg(v) <= lim ){
// 					Q1.insert(v);						
// 				}
// 			}	
// 		}
// 		else{
// 			cg_countingTriangles(triangles);	
// 			set<pair<int, int> > egs;
// 			for (int u = 0; u < cg_maxvtx; u++){
// 				if (cg_valid[u]){
// 					for (int j = cg_pstart[u]; j < cg_pend[u]; j++){
// 						if (u < cg_edges[j] && triangles[j] <= lb-2*k){
// 							assert(triangles[j] != -1);
// 							//cg_removeEdge(u, j); //WARNING: inconsistency
// 							egs.insert(make_pair(u, cg_edges[j]));							
// 						}
// 					}
// 				}
// 			}
// 			set<int> vts;
// 			for (auto e: egs){ 
// 				cg_removeEdge(e.first, e.second);
// 				vts.insert(e.first); vts.insert(e.second);
// 			}
// 			for (auto v: vts){
// 				if (cg_deg(v) + k <= lb){
// 					Q1.insert(v);
// 				}
// 			}
// 			if (Q1.empty())
// 				break;
// 		}
// 	}
// 	delete[] triangles;
// }

// void strongReduction(const MCSRGraph &g, MCSRGraph &newG, int K, int lb, int* new2ori){
// 	//triangle counting
// 	cg_prepareTempGraph(g);	
// 	cg_edgeReduction(K, lb);
// 	cg_back2CSRGraph(newG, lb, new2ori);
// }

// bool peelReduction(const MCSRGraph &orG, MCSRGraph &coreG, const int k, const int lb, int *core, int *dseq, int *dpos, int *new2ori){
// 	const int lim=lb-k;
// 	const int n=orG.nbvtx;
// 	int idx=0;
// 	while(core[dseq[idx]]<=lim)idx++;
// 	if(n-idx<=lb)return true;
// 	int coreSz = 0;
// 	for (int i = idx; i < n; i++){
// 		dseq[coreSz++] = dseq[i];
// 	}

// 	int *ori2core = new int[n];
//     coreG.nbvtx = coreSz;
//     coreG.pstart = new int[coreG.nbvtx+1];
//     coreG.nbedges = 0;
//     coreG.pstart[0] = 0;
// 	coreG.edges = new int[orG.nbedges];
//     int * const edges = coreG.edges;
//     MBitSet64 mark(orG.nbvtx);
//     for (int i = 0; i < coreSz;i++){
// 		const int ele=dseq[i];
//         assert(!(mark.test(ele)));  
//         mark.set(ele);
//         ori2core[ele] = i;
//         new2ori[i] = ele;
//     } 
//     for (int i = 0; i < coreSz; i++){
//         const int u = dseq[i];
//         coreG.pstart[i+1] = coreG.pstart[i];	
//         for (int j = orG.pstart[u]; j < orG.pstart[u+1]; j++){
//             const int nei = orG.edges[j];
//             if (mark.test(nei)){
//                 edges[coreG.pstart[i+1]++] = ori2core[nei]; 
//                 coreG.nbedges++;
//             }
//         }
// 	}
//     delete[] ori2core;
// 	return false;
// }