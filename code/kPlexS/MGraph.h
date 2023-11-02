#pragma once 
#include <algorithm>
#include "Utility.h"
#include "MBitSet.h"
class MCSRGraph{
public:
	int nbedges;
	int nbvtx;
	int *pstart;
	int *edges;
    
    MCSRGraph():nbvtx(0),nbedges(0),pstart(nullptr),edges(nullptr){}
    MCSRGraph(const MCSRGraph& csrg): nbvtx(csrg.nbvtx), nbedges(csrg.nbedges){
        pstart = new int[nbvtx+1];
        edges = new int[nbedges];
        memcpy(pstart, csrg.pstart, sizeof(int) * (nbvtx+1));
        memcpy(edges, csrg.edges, sizeof(int) * nbedges);
    }
    MCSRGraph& operator=(const MCSRGraph &g){
        if (pstart != nullptr) delete[] pstart;
        if (edges != nullptr) delete[] edges;
        nbvtx = g.nbvtx;
        nbedges = g.nbedges;
        pstart = new int[nbvtx+1];
        edges = new int[nbedges];
        memcpy(pstart, g.pstart, sizeof(int) * (nbvtx+1));
        memcpy(edges, g.edges, sizeof(int) * nbedges);
        return *this;
    }

    void del(){
        delete[] pstart;
        delete[] edges;
    }

    void fromBinaryFile(const char* filepath) {
        FILE *f = Utility::open_file(filepath, "rb");
        int tt;
        fread(&tt, sizeof(int), 1, f);
        if (tt != sizeof(int)) {
            printf("sizeof unsigned int is different: file %u, machine %lu\n", tt, sizeof(int));
        }
        fread(&nbvtx, sizeof(int), 1, f);	// the number of vertices
        fread(&nbedges, sizeof(int), 1, f); // the number of edges (twice the acutal number).
 
        int *degree = new int[nbvtx];
        fread(degree, sizeof(int), nbvtx, f);
        if (pstart != nullptr) delete[] pstart;
        pstart = new int[nbvtx + 1];
        if (edges != nullptr) delete[] edges;
        edges = new int[nbedges];

        pstart[0] = 0;
        for (int i = 0; i < nbvtx; i++) {
            if (degree[i] > 0){
                fread(edges + pstart[i], sizeof(int), degree[i], f);
                //std::sort(edges+pstart[i], edges + pstart[i] + degree[i]);
            }
            pstart[i + 1] = pstart[i] + degree[i];
        }
        fclose(f);        
        delete[] degree; 
    }
    
    void toBinaryFile(const char* filepath){
        FILE *f = Utility::open_file(filepath, "wb");
        int tt = sizeof(int);
        fwrite(&tt, sizeof(int), 1, f); //length of ui
        fwrite(&nbvtx, sizeof(int), 1, f);
        fwrite(&nbedges, sizeof(int), 1, f);
        int *d = new int[nbvtx];
        for (int i = 0; i < nbvtx; i++) d[i] =degree(i);
        fwrite(d, sizeof(int), nbvtx, f);
        fwrite(edges, sizeof(int), nbedges, f);
        fclose(f);
    }
    
    void reverseGraph(){
        fprintf(stderr, "function not implemented\n" );
        exit(1);
    }

    int degree(int u) const{
        assert(u >= 0 && u < nbvtx);
        return pstart[u+1]-pstart[u];
    }
};

class MBitGraph{
public:
    int nbvtx;
    int nbedges;
    MBitSet64 **rows;//用邻接矩阵存储图，若存在边(i,j)则将矩阵i行j列赋值为1，求两点的公共邻居可以将row[i]和row[j]按位与在求点1的个数

    MBitGraph():nbvtx(0), nbedges(0), rows(nullptr){}
    
    MBitGraph(const MBitGraph &mbg){
        nbvtx = mbg.nbvtx;
        nbedges = mbg.nbedges;
        rows = new MBitSet64*[nbvtx];
        for (int i = 0; i < nbvtx; i++){
            rows[i] = new MBitSet64(*(mbg.rows[i]));
        }
    }
    
    MBitGraph& operator=(const MBitGraph &mbg) = delete;
    MBitGraph(const MCSRGraph& cg){
        nbvtx = cg.nbvtx;
        nbedges = cg.nbedges;
        rows = new MBitSet64*[nbvtx];
        assert(rows!=nullptr);
        for (int i = 0; i < nbvtx; i++){
            rows[i] = new MBitSet64(nbvtx);
            for (int j = cg.pstart[i]; j < cg.pstart[i+1]; j++){
                rows[i]->set(cg.edges[j]);
            }
        }
    }
    
    void reverseGraph(){
        for (int i = 0; i < nbvtx; i++){
            rows[i]->flip();
            rows[i]->set(i); 
        }
        nbedges = nbvtx * (nbvtx - 1) - nbedges;
    }
    
    void writeToFile(char *filename){
        fprintf(stderr, "function  not implemented\n");
        exit(1);
    }

    int degree(int u){
        return rows[u]->size();
    }

    ~MBitGraph(){
        for (int i = 0; i < nbvtx; i++)
            delete rows[i];
        delete[] rows;
        nbvtx = nbedges = 0;
    }

};

extern int coreDecomposition(const MCSRGraph &g, int *core, int* coreIncSeq, int *pos);
extern bool peelReduction(const MCSRGraph &orG, MCSRGraph &coreG, const int k, const int lb, int *core, int *dseq, int *dpos, int *new2ori);
extern void strongReduction(const MCSRGraph &g, MCSRGraph &newG, int K, int lb, int* new2ori);