#pragma once
#include <cassert>
#include <cstdint>
#include <cstring>
#include <algorithm>
#include <climits>

template <typename T=int> 
class VarSet{    
public:
    using TIndex = uint16_t;
    int sz;
    int cap;
    T *members;
    TIndex *pos; 
    
public:
    VarSet() = delete;
    VarSet(const VarSet &_set) = delete;
    VarSet& operator=(const VarSet &_set) = delete;
    
    VarSet(int _cap):cap(_cap),sz(0){       
        members = new T[_cap];
        pos = new TIndex[_cap];
        memset(pos,0x3f,sizeof(TIndex)*_cap);
    }

    ~VarSet(){
        delete[] members;
        delete[] pos;
    }

    bool contains(const T& v){
        return pos[v]<sz;
    }

    void add(const T& v){
        pos[v] = sz; 
        members[sz++]=v;
    }

    void remove(const T& v){
        TIndex idx = pos[v];            
        pos[members[sz-1]] = idx;
        pos[v] = sz-1;
        swap(members[idx],members[--sz]);
    }

    void recover(const int len){
        sz+=len;
    }

    void clear(){
        sz = 0;
    }

    template<typename Fnc>
    void for_each(const Fnc &f){        
        for (TIndex i = 0; i < sz; i++){
            f(members[i]);
        }
    }
};

