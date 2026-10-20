// Optional BGK_INCREMENTAL_GAIN: exact gain reuse only; frozen original in bgk_campaign_reference/.
// basegen_k.cpp - generic-K port of lrc14/v65/basegen_v65.cpp (LRC(K+1 runners)
// level-1 cover generator, I(K,p,1) family, ST framework).
//
// Algorithm is the audited v6.5 search tree with every hard-coded 13/14/12
// replaced by K / K+1 / K-1.  With -DK=13 the outputs must be byte-identical
// to basegen_v65 (validation gate: p=239 irsubrootrawc(0,0) -> 31,198 rows,
// p=43 both -> direct=decomp=163507).  With -DK=14 it generates the
// I(14,p,1) family for the 15-runner campaign.
//
// Build (K=14, p <= 975 so that n=(p-1)/2 <= 487 <= 512 and m <= 64):
//   g++ -O3 -march=native -DK=14 -DNW=8 -DMAXN=512 -o bgk14 basegen_k.cpp
// Validation build (K=13, matches basegen_v65 wide build):
//   g++ -O3 -march=native -DK=13 -DNW=7 -DMAXN=448 -o bgk13 basegen_k.cpp
//
// Modes: irroots | irsubroots r | irsubrootrawc r s out | irsubrootraw r s out
//        km1roots | km1root r out | km1all out | both | direct | decomp
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <unordered_set>
#include <algorithm>
#include <functional>
#include <utility>
#include <climits>
#include <chrono>
#include <stdexcept>
#if defined(__AVX512VPOPCNTDQ__)
#include <immintrin.h>
#endif
#ifdef BGK_INCREMENTAL_GAIN
#include <memory>
#include "bgk_incremental_gain.h"
#if defined(EXP_MINPRIV) || defined(EXP_SYM)
#error "Gain comparison excludes unrelated experimental modes."
#endif
#endif
using namespace std;

#ifndef K
#define K 14
#endif
#ifndef NW
#define NW 8
#endif
#ifndef MAXN
#define MAXN (64*NW)
#endif
static_assert(MAXN <= 64 * NW, "MAXN must fit in NW words");
static_assert(K >= 3 && K <= 20, "K range");
#define KM1 (K-1)
#define KP1 (K+1)
// every class covers exactly m folded times; m <= 64 for the 64-bit quota mask
#define MAXM 64

#if defined(__AVX512VPOPCNTDQ__) && NW <= 8
#define VEC 1
#define VMASK ((__mmask8)((1u<<NW)-1))
#else
#define VEC 0
#endif

struct Bits {
    uint64_t w[NW] = {};
    inline void set(int i){ w[i>>6] |= 1ULL<<(i&63); }
    inline void reset(int i){ w[i>>6] &= ~(1ULL<<(i&63)); }
    inline bool test(int i)const{ return (w[i>>6]>>(i&63))&1; }
    inline int count()const{ int s=0; for(int j=0;j<NW;j++) s+=__builtin_popcountll(w[j]); return s; }
    inline bool any()const{ for(int j=0;j<NW;j++) if(w[j]) return true; return false; }
};
inline Bits operator|(Bits a,Bits b){ Bits r; for(int j=0;j<NW;j++) r.w[j]=a.w[j]|b.w[j]; return r; }
inline Bits operator&(Bits a,Bits b){ Bits r; for(int j=0;j<NW;j++) r.w[j]=a.w[j]&b.w[j]; return r; }
inline Bits andnot(Bits a,Bits b){ Bits r; for(int j=0;j<NW;j++) r.w[j]=a.w[j]&~b.w[j]; return r; }
inline bool covers_all(const Bits&c,const Bits&all){ for(int j=0;j<NW;j++) if((c.w[j]&all.w[j])!=all.w[j]) return false; return true; }

struct Generator {
#ifdef BGK_INCREMENTAL_GAIN
    std::unique_ptr<BgkIncrementalGain<MAXN>> gain_kernel;
#ifdef BGK_AUDIT_GAINS
    mutable uint64_t gain_values_checked=0,gain_delta_calls=0;
    ~Generator(){cerr<<"GAIN_AUDIT values="<<gain_values_checked<<" deltas="<<gain_delta_calls<<"\n";}
#endif
#endif
    int p,n,m;
    vector<Bits> cov;   // 0-based classes
    vector<Bits> cand;  // per point, covering classes
#if VEC
    alignas(64) array<array<uint64_t,8>,MAXN> cov8{},cand8{};
#endif
    vector<vector<uint16_t>> trans; // g class x v class => class
    vector<uint16_t> invfold;       // folded inverse multiplier of each class
    vector<uint64_t> qm;            // which C1 points class v covers (bit i = i-th C1 point)
    Bits ALL, C1;
    uint64_t nodes=0, leaves=0, canon_calls=0;
#ifdef EXP_SYM
    // dilation-redundancy measurement: per depth, visits vs distinct orbit keys of the partial set
    vector<unordered_set<uint64_t>> sym_sets; vector<uint64_t> sym_visits;
    void sym_record(const vector<int>&sel){
        int d=(int)sel.size(); if((int)sym_sets.size()<=d){sym_sets.resize(d+1);sym_visits.resize(d+1);}
        sym_visits[d]++; if(d<EXP_SYM||d>K-3)return;
        uint64_t best=~0ull; int img[K+2];
        for(int s:sel){ int g=invfold[s-1];           // 1-based folded inverse
            for(int i=0;i<d;i++)img[i]=trans[g-1][sel[i]-1];
            sort(img,img+d);
            uint64_t h=1469598103934665603ull; for(int i=0;i<d;i++){h^=(uint64_t)img[i]+1;h*=1099511628211ull;h^=h>>29;}
            if(h<best)best=h; }
        sym_sets[d].insert(best);
    }
#endif
    unordered_set<string> out;
    bool irredundant=false;
    int quota_limit=-1;
    int coordK=K;
    int hexw, wb;

    explicit Generator(int P):p(P),n((P-1)/2),m((P-1)/KP1),cov(n),cand(n),trans(n,vector<uint16_t>(n)),invfold(n){
        if(n>MAXN){cerr<<"p too large for MAXN\n";exit(2);}
        hexw = (n>255) ? 4 : 2;
        wb = (hexw==4) ? 9 : 8;
        for(int t=0;t<n;t++) ALL.set(t);
        for(int v=1;v<=n;v++){
            for(int t=1;t<=n;t++){
                int r=(int)(1LL*v*t%p); int d=min(r,p-r);
                if((long long)KP1*d<p){cov[v-1].set(t-1);cand[t-1].set(v-1);}
            }
        }
        for(int v=0;v<n;v++) if(cov[v].count()!=m){cerr<<"class "<<v+1<<" covers "<<cov[v].count()<<" != m="<<m<<"\n";exit(2);}
        for(int t=0;t<n;t++) if(cand[t].count()!=m){cerr<<"point "<<t+1<<" has "<<cand[t].count()<<" != m covering classes\n";exit(2);}
        C1=cov[0];
#if VEC
        for(int v=0;v<n;v++)for(int j=0;j<NW;j++)cov8[v][j]=cov[v].w[j];
        for(int t=0;t<n;t++)for(int j=0;j<NW;j++)cand8[t][j]=cand[t].w[j];
#endif
        if(m>64){cerr<<"m too large for 64-bit quota mask\n";exit(2);}
        {
            int c1pts[64],k=0;
            for(int t=0;t<n;t++)if(C1.test(t))c1pts[k++]=t;
            qm.assign(n,0);
            for(int v=0;v<n;v++)
                for(int i=0;i<k;i++)
                    if(cov[v].test(c1pts[i]))qm[v]|=1ULL<<i;
        }
        for(int g=1;g<=n;g++) for(int v=1;v<=n;v++){
            int r=(int)(1LL*g*v%p); int c=min(r,p-r);
            trans[g-1][v-1]=(uint16_t)c;
            if(c==1)invfold[v-1]=(uint16_t)g;
        }
        out.reserve(1<<20);
#ifdef BGK_INCREMENTAL_GAIN
        gain_kernel=std::make_unique<BgkIncrementalGain<MAXN>>(*this);
#endif
    }

    string pack(const vector<uint16_t>& v)const{
        if(wb==8){
            string z;z.reserve(v.size());
            for(uint16_t c:v)z.push_back((char)(c&255));
            return z;
        }
        string z(((size_t)wb*v.size()+7)/8,'\0');
        size_t bit=0;
        for(uint16_t c:v)
            for(int b=wb-1;b>=0;b--,bit++)
                if((c>>b)&1) z[bit>>3]|=(char)(0x80>>(bit&7));
        return z;
    }
    int unpack(const string& s,int idx)const{
        if(wb==8)return (unsigned char)s[(size_t)idx];
        int c=0;size_t bit=(size_t)wb*idx;
        for(int b=0;b<wb;b++,bit++)
            c=(c<<1)|(((unsigned char)s[bit>>3]>>(7-(bit&7)))&1);
        return c;
    }

    string canonical_k(const vector<int>& support,const vector<int>& mult,int KK){
        canon_calls++;
        vector<uint16_t> best(KK,0xFFFF),cur;cur.reserve(KK);
        for(int g=0;g<n;g++){
            cur.clear();
            for(size_t j=0;j<support.size();j++){
                uint16_t c=trans[g][support[j]-1];
                for(int z=0;z<mult[j];z++)cur.push_back(c);
            }
            sort(cur.begin(),cur.end());if(cur<best)best=cur;
        }
        return pack(best);
    }
    string canonical(const vector<int>& support,const vector<int>& mult){ return canonical_k(support,mult,coordK); }

    void compositions_rec(const vector<int>& support,int idx,int remaining,vector<int>& mult,
                          bool force_mult1_one, function<void(const vector<int>&)> emit){
        int s=(int)support.size();
        if(idx==s-1){
            if(remaining>=1){mult[idx]=remaining;emit(mult);}return;
        }
        int lo=1, hi=remaining-(s-idx-1);
        if(force_mult1_one && idx==0){lo=hi=1;}
        for(int x=lo;x<=hi;x++){
            mult[idx]=x;
            compositions_rec(support,idx+1,remaining-x,mult,force_mult1_one,emit);
        }
    }

    bool final_irredundant(const vector<int>& sel){
        vector<int> cnt(n,0);
        for(int v:sel) for(int t=0;t<n;t++) if(cov[v-1].test(t)) cnt[t]++;
        for(int v:sel){bool priv=false;for(int t=0;t<n;t++)if(cov[v-1].test(t)&&cnt[t]==1){priv=true;break;}if(!priv)return false;}
        return true;
    }
    bool partial_irredundant(const vector<int>& sel){return final_irredundant(sel);}

    void emit_direct_support(const vector<int>& sel){
        leaves++;
        vector<int> mult(sel.size());
        compositions_rec(sel,0,coordK,mult,false,[&](const vector<int>& mm){out.insert(canonical(sel,mm));});
    }
    void emit_km1_support(const vector<int>& sel, unordered_set<string>& coversKm1){
        leaves++;
        vector<int> mult(sel.size());
        compositions_rec(sel,0,KM1,mult,true,[&](const vector<int>& mm){
            coversKm1.insert(canonical_k(sel,mm,KM1));
        });
    }

    inline int gain_of(int v,const Bits&unc
#if VEC
        ,__m512i uncvec
#endif
    )const{
#if VEC
        __m512i cv=_mm512_load_si512((const void*)cov8[v].data());
        return (int)_mm512_reduce_add_epi64(_mm512_popcnt_epi64(_mm512_and_si512(cv,uncvec)));
#else
        return (cov[v]&unc).count();
#endif
    }

    int rarest_uncovered(const Bits&unc,const Bits&avail){
        int bt=-1,bc=INT_MAX;
#if VEC
        __m512i av=_mm512_maskz_loadu_epi64(VMASK,(const void*)avail.w);
#endif
        for(int j=0;j<NW;j++){
            uint64_t w=unc.w[j];
            while(w){
                int t=64*j+__builtin_ctzll(w); w&=w-1;
#if VEC
                __m512i cv=_mm512_load_si512((const void*)cand8[t].data());
                int c=(int)_mm512_reduce_add_epi64(_mm512_popcnt_epi64(_mm512_and_si512(cv,av)));
#else
                int c=(cand[t]&avail).count();
#endif
                if(c<bc){bc=c;bt=t;if(c<=1)return bt;}
            }
        }
        return bt;
    }

    bool gain_bound(const Bits&unc,int need,const Bits&avail,int slots,int&bsum){
        if(need>slots*m)return false;
        if(slots<=4){
            int top[4]={},ng=0;
            for(int j=0;j<NW;j++){
                uint64_t w=avail.w[j];
                while(w){
                    int v=64*j+__builtin_ctzll(w);w&=w-1;ng++;
                    int g=(cov[v]&unc).count();
                    for(int k=0;k<slots;k++)if(g>top[k]){
                        for(int z=slots-1;z>k;z--)top[z]=top[z-1];
                        top[k]=g;break;
                    }
                }
            }
            if(ng<slots)return false;
            int sum=0;for(int k=0;k<slots;k++)sum+=top[k];
            bsum=sum-top[slots-1];
            return sum>=need;
        }
        int gains[MAXN]; int ng=0; long long total=0;
        for(int j=0;j<NW;j++){
            uint64_t w=avail.w[j];
            while(w){
                int v=64*j+__builtin_ctzll(w); w&=w-1;
                int g=(cov[v]&unc).count(); gains[ng++]=g; total+=g;
            }
        }
        if(ng<slots)return false;
        if(total<need)return false;
        if(ng>slots)nth_element(gains,gains+slots,gains+ng,greater<int>());
        long long sum=0;int mn=INT_MAX;
        for(int i=0;i<slots;i++){sum+=gains[i];mn=min(mn,gains[i]);}
        bsum=(int)(sum-mn);
        return sum>=need;
    }

    int ordered_choices(const Bits&choices,const Bits&unc,int*outv,int*outg=nullptr){
        pair<int,int> buf[MAXM]; int k=0;
        for(int j=0;j<NW;j++){
            uint64_t w=choices.w[j];
            while(w){
                int v=64*j+__builtin_ctzll(w); w&=w-1;
                buf[k++]={-(int)(cov[v]&unc).count(),v};
            }
        }
        sort(buf,buf+k);
        for(int i=0;i<k;i++){outv[i]=buf[i].second;if(outg)outg[i]=-buf[i].first;}
        return k;
    }

    // exact gains over avail; proves the node bound and feeds child ordering
    bool gain_bound_precomputed(const Bits&unc,int need,const Bits&avail,int slots,
                                int&bsum,uint8_t*gains,int*hist
#ifdef BGK_INCREMENTAL_GAIN
                                ,const uint8_t*parent_gains=nullptr,const Bits*removed=nullptr
#endif
                                )const{
        if(need>slots*m)return false;
#ifdef BGK_INCREMENTAL_GAIN
        gain_kernel->compute(unc,need,gains,parent_gains,removed);
#ifdef BGK_AUDIT_GAINS
        if(parent_gains&&removed&&removed->count()<need)++gain_delta_calls;
        for(int v=0;v<n;++v){
            if(gains[v]!=(cov[v]&unc).count())throw runtime_error("cached gain mismatch");
            ++gain_values_checked;
        }
        for(int v=n;v<gain_kernel->stride();++v)
            if(gains[v])throw runtime_error("gain padding mismatch");
#endif
#endif
#ifdef BGK_THRESHOLD_WITNESS
#ifndef BGK_INCREMENTAL_GAIN
#error "BGK_THRESHOLD_WITNESS requires BGK_INCREMENTAL_GAIN."
#endif
        // A cover with 'slots' choices needs at least one available gain
        // >= ceil(need/slots). No maximum value or sorted list is needed.
        if(need>0){
            const __m256i threshold=_mm256_set1_epi8((char)((need+slots-1)/slots));
            bool possible=false;
            for(int base=0;base<gain_kernel->stride();base+=32){
                const __m256i values=_mm256_loadu_si256((const __m256i*)(gains+base));
                const unsigned eligible=(unsigned)_mm256_movemask_epi8(
                    _mm256_cmpeq_epi8(_mm256_max_epu8(values,threshold),values));
                const uint32_t available=(uint32_t)(avail.w[base/64]>>(base%64));
                if(eligible&available){possible=true;break;}
            }
            if(!possible)return false;
        }
#endif
        memset(hist,0,(MAXM+1)*sizeof(*hist));
        int bestT[4]={},probe[4],np=0,ng=0;
#if VEC
        __m512i uncvec=_mm512_maskz_loadu_epi64(VMASK,(const void*)unc.w);
#endif
#ifndef BGK_DEFERRED_PROBES
        if(slots==2||slots==3){
            for(int j=0;j<NW&&np<4;j++){
                uint64_t w=unc.w[j];
                while(w&&np<4){probe[np++]=64*j+__builtin_ctzll(w);w&=w-1;}
            }
        }
#endif
        for(int j=0;j<NW;j++){
            uint64_t w=avail.w[j];
            while(w){
                int v=64*j+__builtin_ctzll(w);w&=w-1;++ng;
#ifdef BGK_INCREMENTAL_GAIN
                int g=gains[v];
#elif VEC
                int g=gain_of(v,unc,uncvec);
#else
                int g=(cov[v]&unc).count();
#endif
                gains[v]=(uint8_t)g;++hist[g];
#ifndef BGK_DEFERRED_PROBES
                for(int q=0;q<np;q++)if(cov[v].test(probe[q]))
                    bestT[q]=max(bestT[q],g);
#endif
            }
        }
        if(ng<slots)return false;
        int top[3]={},left=slots,sum=0,last=0,ntop=0;
        for(int g=m;g>=0&&left;g--){
            int take=min(left,hist[g]);sum+=take*g;left-=take;last=g;
            for(int z=0;z<take&&ntop<3;z++)top[ntop++]=g;
        }
        bsum=sum-last;
        if(sum<need)return false;
#ifdef BGK_DEFERRED_PROBES
        // Same test, evaluated after the global bound: max gain among
        // candidates covering t >= need - sum of the largest slots-1 gains.
        // cand[t] is the exact transpose of cov. Only an existence test is
        // needed, so stop at the first qualifying candidate. Empty maxima
        // have value zero, exactly as in the original bestT initialization.
        if(slots==2||slots==3){
            int threshold=need-top[0]-(slots==3?top[1]:0);
            if(threshold>0){
                int checked=0;
                for(int j=0;j<NW&&checked<4;++j){
                    uint64_t w=unc.w[j];
                    while(w&&checked<4){
                        int t=64*j+__builtin_ctzll(w);w&=w-1;++checked;
                        bool found=false;
                        for(int h=0;h<NW&&!found;++h){
                            uint64_t eligible=cand[t].w[h]&avail.w[h];
                            while(eligible){
                                int v=64*h+__builtin_ctzll(eligible);eligible&=eligible-1;
                                if(gains[v]>=threshold){found=true;break;}
                            }
                        }
#ifdef BGK_AUDIT_PROBES
                        int reference=0;
                        for(int v=0;v<n;++v)if(avail.test(v)&&cov[v].test(t))
                            reference=max(reference,(int)gains[v]);
                        if(found!=(reference>=threshold))throw runtime_error("deferred probe mismatch");
#endif
                        if(!found)return false;
                    }
                }
            }
        }
#else
        if(slots==2)
            for(int q=0;q<np;q++)if(top[0]+bestT[q]<need)return false;
        if(slots==3)
            for(int q=0;q<np;q++)if(top[0]+top[1]+bestT[q]<need)return false;
#endif
        return true;
    }

    int ordered_choices_gains(const Bits&choices,const uint8_t*gains,int*outv,int*outg=nullptr)const{
        uint16_t val[MAXM];uint8_t gv[MAXM];int count[MAXM+1]={},k=0;
        for(int j=0;j<NW;j++){
            uint64_t w=choices.w[j];
            while(w){
                int v=64*j+__builtin_ctzll(w);w&=w-1;int g=gains[v];
                val[k]=(uint16_t)v;gv[k]=(uint8_t)g;++count[g];++k;
            }
        }
        int next[MAXM+1],pos=0;
        for(int g=m;g>=0;g--){next[g]=pos;pos+=count[g];}
        for(int i=0;i<k;i++){
            int g=gv[i],at=next[g]++;
            outv[at]=val[i];if(outg)outg[at]=g;
        }
        return k;
    }

    struct TopKState {int sum=INT_MIN/2,cut=-1,used_cut=0;};
    TopKState init_topk(const int*hist,int k)const{
        TopKState s;int left=k,sum=0;
        for(int g=m;g>=0&&left;g--){
            int take=min(left,hist[g]);
            if(take){sum+=take*g;left-=take;s.cut=g;s.used_cut=take;}
        }
        if(!left)s.sum=sum;
        return s;
    }
    void remove_topk(int*hist,TopKState&s,int g)const{
        int before=hist[g];--hist[g];
        if(s.cut<0||g<s.cut)return;
        if(g==s.cut&&before>s.used_cut)return;
        if(g>s.cut&&hist[s.cut]>=s.used_cut+1){
            s.sum+=s.cut-g;++s.used_cut;return;
        }
        int next=s.cut-1;while(next>=0&&!hist[next])--next;
        if(next<0){s.sum=INT_MIN/2;s.cut=-1;s.used_cut=0;return;}
        s.sum+=next-g;s.cut=next;s.used_cut=1;
    }

    // general support enumeration (any support size), fixed class 1
    template<class Emit>
    void dfs_support(vector<int>&sel,Bits covered,Bits avail,int slots,Emit emit){
        nodes++;
        if(quota_limit>=0){
            uint64_t u=0;for(size_t j=1;j<sel.size();j++)u|=qm[sel[j]-1];
            if(__builtin_popcountll(u)>quota_limit)return;
        }
        if(irredundant && !partial_irredundant(sel))return;
        if(slots==0){if(covers_all(covered,ALL))emit(sel);return;}
        if(avail.count()<slots)return;
        Bits unc=andnot(ALL,covered);int need=unc.count();
        if(need==0){
            vector<int> av;for(int v=0;v<n;v++)if(avail.test(v))av.push_back(v+1);
            function<void(int,int,Bits,Bits)> fill=[&](int pos,int need,Bits covd,Bits avbits){
                if(need==0){emit(sel);return;}
                if((int)av.size()-pos<need)return;
                for(int j=pos;j<=(int)av.size()-need;j++){
                    int v=av[j];sel.push_back(v);Bits nextavail=avbits;nextavail.reset(v-1);
                    if(quota_limit<0 || ([&](){uint64_t u=0;for(size_t k=1;k<sel.size();k++)u|=qm[sel[k]-1];return __builtin_popcountll(u)<=quota_limit;})()){
                        if(!irredundant || partial_irredundant(sel))fill(j+1,need-1,covd|cov[v-1],nextavail);
                    }
                    sel.pop_back();
                }
            };
            fill(0,slots,covered,avail);return;
        }
        int bsum;
        if(!gain_bound(unc,need,avail,slots,bsum))return;
        int bt=rarest_uncovered(unc,avail);if(bt<0)return;
        Bits choices=cand[bt]&avail;
        int ch[MAXM],gg[MAXM];int nch=ordered_choices(choices,unc,ch,gg);
        Bits rem=avail;
        for(int ci=0;ci<nch;ci++){
            if(gg[ci]+bsum<need)break;
            int vi=ch[ci];
            rem.reset(vi);
            sel.push_back(vi+1);
            dfs_support(sel,covered|cov[vi],rem,slots-1,emit);
            sel.pop_back();
        }
    }

    vector<int> root_choices(int &bt_out){
        vector<int>sel{1};Bits covered=C1,avail=ALL;avail.reset(0);
        Bits unc=andnot(ALL,covered);
        int bt=rarest_uncovered(unc,avail);bt_out=bt;
        Bits choices=cand[bt]&avail;vector<int>ch;for(int v=0;v<n;v++)if(choices.test(v))ch.push_back(v);
        sort(ch.begin(),ch.end(),[&](int a,int b){int ga=(cov[a]&unc).count(),gb=(cov[b]&unc).count();if(ga!=gb)return ga>gb;return a<b;});
        return ch;
    }

    // v6.6 bitmask irredundance state: covered2 = points covered >= 2 times,
    // priv[j] = current private points of the j-th selected class.  A class
    // only LOSES private points during an add (a point moves cnt 1->2 exactly
    // once), so "some old class ends with zero private points" is detected
    // exactly as with per-point counters.
    struct IrState {
        Bits covered2;
        Bits priv[K];
    };
    void init_irstate(IrState&st){ st.covered2=Bits(); st.priv[0]=C1; }
    bool add_irred(IrState &st,const Bits&covered,int pos,int v){
        const Bits&cv=cov[v-1];
        Bits newpriv=andnot(cv,covered);
        Bits newly2=andnot(cv&covered,st.covered2);
        st.covered2=st.covered2|newly2;
        st.priv[pos]=newpriv;
        if(newly2.any())
            for(int j=0;j<pos;j++){
                Bits np=andnot(st.priv[j],newly2);
                st.priv[j]=np;
                if(!np.any())return false;
            }
        return newpriv.any();
    }
    struct IrUndo { Bits covered2; Bits priv[K]; int n=0; };
    bool add_irred_logged(IrState &st,const Bits&covered,int pos,int v,IrUndo &u){
        u.covered2=st.covered2;u.n=pos;
        for(int j=0;j<pos;j++)u.priv[j]=st.priv[j];
        return add_irred(st,covered,pos,v);
    }
    void undo_irred(IrState &st,const IrUndo &u){
        st.covered2=u.covered2;
        for(int j=0;j<u.n;j++)st.priv[j]=u.priv[j];
    }

#ifdef EXP_MINPRIV
    // Experiment: "min-private rooting".  Every unit orbit of covers has a
    // rooting (translate containing class 1) in which class 1 has the minimum
    // private count.  Then priv(1) <= #priv/KK <= (n-|covered2|)/KK, i.e. the
    // other classes must double-cover at least D = m - floor((n-|covered2|)/KK)
    // - |C1 & covered2| more points of class 1.  Bound: the top-`slots`
    // values of |cov[v] & priv[0]| over avail must sum to >= D.  This pass
    // only MEASURES how many nodes that bound would cut (cut_state marks the
    // subtree) and checks completeness of the min-priv leaves at the orbit level.
    long long nodes_cut=0,leaves_all=0,leaves_minpriv=0;
    long long nodes_by_slots[K+1]={},pruned_gain[K+1]={},dual_would_prune[K+1]={};
    bool cut_state=false;IrState*cur_state=nullptr;
    unordered_set<string> orb_all,orb_minpriv;
    string abs_canon(const vector<int>&S)const{
        uint64_t best[NW];bool have=false;
        for(int j=0;j<K;j++){
            int g=invfold[S[(size_t)j]-1];
            uint64_t bits[NW]={};
            for(int z=0;z<K;z++){int c=(g==1)?S[(size_t)z]:trans[g-1][S[(size_t)z]-1];int v=c-1;bits[v>>6]|=1ULL<<(v&63);}
            if(!have){memcpy(best,bits,sizeof bits);have=true;continue;}
            for(int k=0;k<NW;k++){uint64_t d=bits[k]^best[k];if(d){if(bits[k]&d&-d)memcpy(best,bits,sizeof bits);break;}}
        }
        return string((const char*)best,sizeof best);
    }
    bool minpriv_cut(const Bits&avail,int slots,const vector<int>&sel,const IrState&state){
        int KK=(int)sel.size()+slots;
        int c2=state.covered2.count();int c1cov=m-state.priv[0].count();
        int P=(n-c2)/KK;int D=m-P-c1cov;
        if(D<=0)return false;
        int top[K]={};
        for(int j=0;j<NW;j++){uint64_t w=avail.w[j];while(w){int v=64*j+__builtin_ctzll(w);w&=w-1;
            int h=(cov[v]&state.priv[0]).count();
            for(int k=0;k<slots;k++)if(h>top[k]){for(int z=slots-1;z>k;z--)top[z]=top[z-1];top[k]=h;break;}}}
        int s=0;for(int k=0;k<slots;k++)s+=top[k];
        return s<D;
    }
    void minpriv_leaf(const vector<int>&S){
        leaves_all++;const IrState&st=*cur_state;
        int p1=st.priv[0].count();bool ismin=true;
        for(int j=1;j<(int)S.size();j++)if(st.priv[j].count()<p1){ismin=false;break;}
        string c=abs_canon(S);orb_all.insert(c);
        if(ismin){leaves_minpriv++;orb_minpriv.insert(c);}
    }
    template<class Emit>
    void dfs_irred(vector<int>&sel,Bits covered,Bits avail,int slots,uint64_t qmask,IrState&state,Emit emit
#ifdef BGK_INCREMENTAL_GAIN
                   ,const uint8_t*parent_gains=nullptr,const Bits*removed=nullptr
#endif
                   ){
        bool saved=cut_state;cur_state=&state;
        if(cut_state)nodes_cut++;               // subtree below a dual-bound cut
        // (min-priv D-cut measured earlier: 68 of 246,819 nodes at p=83 -> dead idea)
        dfs_irred_impl(sel,covered,avail,slots,qmask,state,emit);
        cut_state=saved;cur_state=&state;
    }
    template<class Emit>
    void dfs_irred_impl(vector<int>&sel,Bits covered,Bits avail,int slots,uint64_t qmask,IrState&state,Emit emit){
#else
    template<class Emit>
    void dfs_irred(vector<int>&sel,Bits covered,Bits avail,int slots,uint64_t qmask,IrState&state,Emit emit
#ifdef BGK_INCREMENTAL_GAIN
                   ,const uint8_t*parent_gains=nullptr,const Bits*removed=nullptr
#endif
                   ){
#endif
        nodes++;
#ifdef EXP_SYM
        sym_record(sel);
#endif
#ifdef EXP_MINPRIV
        nodes_by_slots[slots]++;
#endif
        if(__builtin_popcountll(qmask)>quota_limit)return;
        if(slots==0){if(covers_all(covered,ALL))emit(sel);return;}
        Bits unc=andnot(ALL,covered);int need=unc.count();
        if(need<slots)return;
        if(avail.count()<slots)return;
#ifdef BGK_INCREMENTAL_GAIN
        alignas(32) uint8_t gains[(MAXN+31)&~31];
#else
        uint8_t gains[MAXN];
#endif
        int bsum,gain_hist[MAXM+1];
        if(!gain_bound_precomputed(unc,need,avail,slots,bsum,gains,gain_hist
#ifdef BGK_INCREMENTAL_GAIN
                                  ,parent_gains,removed
#endif
                                  )){
#ifdef EXP_MINPRIV
            pruned_gain[slots]++;
#endif
            return;
        }
#ifdef EXP_MINPRIV
        {   // dual bound y_t = 1/gmax(t): LB on classes still needed
            static double inv[MAXM+1];if(inv[1]==0){for(int g=1;g<=MAXM;g++)inv[g]=1.0/g;}
            uint8_t gmax[MAXN];memset(gmax,0,sizeof gmax);
            for(int j=0;j<NW;j++){uint64_t w=avail.w[j];while(w){int v=64*j+__builtin_ctzll(w);w&=w-1;uint8_t g=gains[v];if(!g)continue;
                Bits cu=cov[v]&unc;for(int jj=0;jj<NW;jj++){uint64_t ww=cu.w[jj];while(ww){int t=64*jj+__builtin_ctzll(ww);ww&=ww-1;if(g>gmax[t])gmax[t]=g;}}}}
            double lb=0;for(int j=0;j<NW;j++){uint64_t w=unc.w[j];while(w){int t=64*j+__builtin_ctzll(w);w&=w-1;lb+=inv[gmax[t]];}}
            if(lb>slots+1e-9){dual_would_prune[slots]++;if(!cut_state){cut_state=true;nodes_cut++;}}
        }
#endif
        if(slots==2){
            int b1=bsum;
            int bt=rarest_uncovered(unc,avail);if(bt<0)return;
            Bits choices=cand[bt]&avail;
            int ch[MAXM],gg[MAXM];int nch=ordered_choices_gains(choices,gains,ch,gg);
            Bits rem=avail;TopKState later=init_topk(gain_hist,1);
            for(int ci=0;ci<nch;ci++){
                if(gg[ci]+b1<need)break;
                int vi=ch[ci];
                remove_topk(gain_hist,later,gains[vi]);
                if(gg[ci]+later.sum<need)break;
                rem.reset(vi); IrUndo u1;
                bool ok1=add_irred_logged(state,covered,(int)sel.size(),vi+1,u1);
                if(!ok1){undo_irred(state,u1);continue;}
                uint64_t q1=qmask|qm[vi];
                if(__builtin_popcountll(q1)>quota_limit){undo_irred(state,u1);continue;}
                sel.push_back(vi+1);
                Bits cov1=covered|cov[vi], unc1=andnot(ALL,cov1);
                if(unc1.any()){
                    Bits completions=rem;
                    for(int j=0;j<NW;j++){
                        uint64_t w=unc1.w[j];
                        while(w){
                            int t=64*j+__builtin_ctzll(w); w&=w-1;
                            completions=completions&cand[t];
                            if(!completions.any())break;
                        }
                        if(!completions.any())break;
                    }
                    if(completions.any()){
                        int bt1=rarest_uncovered(unc1,rem);
                        Bits choices1=cand[bt1]&rem;
                        int ch1[MAXM];int nch1=ordered_choices(choices1,unc1,ch1);
                        for(int cj=0;cj<nch1;cj++){
                            int wi=ch1[cj];if(!completions.test(wi))continue;
                            IrUndo u2;
                            bool ok2=add_irred_logged(state,cov1,(int)sel.size(),wi+1,u2);
                            if(ok2){
                                uint64_t q2=q1|qm[wi];
                                if(__builtin_popcountll(q2)<=quota_limit){
                                    sel.push_back(wi+1);emit(sel);sel.pop_back();
                                }
                            }
                            undo_irred(state,u2);
                        }
                    }
                }
                sel.pop_back();undo_irred(state,u1);
            }
            return;
        }
        int bt=rarest_uncovered(unc,avail);if(bt<0)return;
        Bits choices=cand[bt]&avail;
        int ch[MAXM],gg[MAXM];int nch=ordered_choices_gains(choices,gains,ch,gg);
        Bits rem=avail;TopKState later=init_topk(gain_hist,slots-1);
        for(int ci=0;ci<nch;ci++){
            if(gg[ci]+bsum<need)break;
            int vi=ch[ci];
            remove_topk(gain_hist,later,gains[vi]);
            if(gg[ci]+later.sum<need)break;
            rem.reset(vi);IrUndo u;
            bool ok=add_irred_logged(state,covered,(int)sel.size(),vi+1,u);
            if(ok){
                uint64_t nq=qmask|qm[vi];
                if(__builtin_popcountll(nq)<=quota_limit){
                    sel.push_back(vi+1);
#ifdef BGK_INCREMENTAL_GAIN
                    Bits newly_covered=unc&cov[vi];
                    dfs_irred(sel,covered|cov[vi],rem,slots-1,nq,state,emit,gains,&newly_covered);
#else
                    dfs_irred(sel,covered|cov[vi],rem,slots-1,nq,state,emit);
#endif
                    sel.pop_back();
                }
            }
            undo_irred(state,u);
        }
    }

    // ---- K-1 covers (reducible base) ----
    void set_quota_km1(){ int q=(max(0,2*n-KM1*m)+KM1-1)/KM1;quota_limit=m-q;irredundant=false; }
    void set_quota_k(){ coordK=K;int q=(max(0,2*n-K*m)+K-1)/K;quota_limit=m-q;irredundant=false; }

    unordered_set<string> enumerate_km1root(int ridx){
        unordered_set<string> covers;covers.reserve(1<<18);
        set_quota_km1();
        int bt=-1;auto ch=root_choices(bt);if(ridx<0||ridx>=(int)ch.size())throw runtime_error("bad root index");
        Bits avail=ALL;avail.reset(0);Bits rem=avail;
        for(int j=0;j<=ridx;j++)rem.reset(ch[j]);
        int vi=ch[ridx];vector<int>sel{1};
        IrState st;init_irstate(st);
        if(!add_irred(st,C1,(int)sel.size(),vi+1))return covers;
        sel.push_back(vi+1);
        dfs_irred(sel,C1|cov[vi],rem,KM1-2,qm[vi],st,[&](const vector<int>&S){emit_km1_support(S,covers);});
        return covers;
    }
    unordered_set<string> enumerate_km1covers_upto(int smax){
        unordered_set<string> covers;covers.reserve(1<<20);
        set_quota_km1();
        for(int s=1;s<=smax;s++){
            vector<int>sel{1};Bits avail=ALL;avail.reset(0);
            dfs_support(sel,C1,avail,s-1,[&](const vector<int>&S){emit_km1_support(S,covers);});
        }
        return covers;
    }
    unordered_set<string> enumerate_km1covers(){
        unordered_set<string> covers;covers.reserve(1<<20);
        set_quota_km1();
        for(int s=1;s<=KM1;s++){
            vector<int>sel{1};Bits avail=ALL;avail.reset(0);
            dfs_support(sel,C1,avail,s-1,[&](const vector<int>&S){emit_km1_support(S,covers);});
        }
        return covers;
    }
    void reducible_from_km1(const unordered_set<string>&covers){
        coordK=K;
        vector<int> support,mult;
        for(auto &key:covers){
            vector<int> base(KM1);
            for(int i=0;i<KM1;i++)base[i]=unpack(key,i);
            for(int v=1;v<=n;v++){
                vector<int> full=base;full.push_back(v);sort(full.begin(),full.end());
                support.clear();mult.clear();
                for(int c:full){if(support.empty()||support.back()!=c){support.push_back(c);mult.push_back(1);}else mult.back()++;}
                out.insert(canonical(support,mult));
            }
        }
    }

    // ---- irredundant K-covers ----
    unordered_set<string> enumerate_irroot(int ridx){
        unordered_set<string> z;z.reserve(1<<20);set_quota_k();
        int bt=-1;auto ch=root_choices(bt);if(ridx<0||ridx>=(int)ch.size())throw runtime_error("bad irroot index");
        Bits avail=ALL;avail.reset(0);Bits rem=avail;for(int j=0;j<=ridx;j++)rem.reset(ch[j]);
        int vi=ch[ridx];vector<int>sel{1};
        IrState st;init_irstate(st);
        if(!add_irred(st,C1,(int)sel.size(),vi+1))return z;
        sel.push_back(vi+1);
        dfs_irred(sel,C1|cov[vi],rem,K-2,qm[vi],st,[&](const vector<int>&S){vector<int>mm(S.size(),1);z.insert(canonical(S,mm));leaves++;});
        return z;
    }

    string hexkey(const string& s)const{
        static const char*h="0123456789abcdef";
        int k=(int)(8*s.size()/(size_t)wb);
        string z;z.reserve((size_t)hexw*k);
        for(int i=0;i<k;i++){
            int c=unpack(s,i);
            if(hexw==4){z.push_back(h[(c>>12)&15]);z.push_back(h[(c>>8)&15]);}
            z.push_back(h[(c>>4)&15]);z.push_back(h[c&15]);
        }
        return z;
    }

    struct IrPrefix2 {
        vector<int> sel;
        Bits covered, avail;
        uint64_t other=0;
        IrState st;
        vector<int> choices;
        int bt=-1;
    };
    IrPrefix2 make_irprefix2(int ridx){
        set_quota_k();
        int bt0=-1;auto ch0=root_choices(bt0);if(ridx<0||ridx>=(int)ch0.size())throw runtime_error("bad ir root");
        Bits avail=ALL;avail.reset(0);for(int j=0;j<=ridx;j++)avail.reset(ch0[j]);
        int vi=ch0[ridx];vector<int>sel{1};
        IrState st;init_irstate(st);
        if(!add_irred(st,C1,(int)sel.size(),vi+1))throw runtime_error("dead ir root");
        sel.push_back(vi+1);Bits covered=C1|cov[vi];
        Bits unc=andnot(ALL,covered);
        int bt=rarest_uncovered(unc,avail);Bits choices=cand[bt]&avail;vector<int>ch;
        for(int v=0;v<n;v++)if(choices.test(v))ch.push_back(v);
        sort(ch.begin(),ch.end(),[&](int a,int b){int ga=(cov[a]&unc).count(),gb=(cov[b]&unc).count();if(ga!=gb)return ga>gb;return a<b;});
        return {sel,covered,avail,qm[vi],st,ch,bt};
    }

    void write_rawhex(ofstream&f,const vector<int>&S)const{
        uint16_t row[K];
        for(int i=0;i<K;i++)row[i]=(uint16_t)S[(size_t)i];
        sort(row,row+K);
        static const char*h="0123456789abcdef";
        char outrow[K*4+1];char*q=outrow;
        for(uint16_t c:row){
            if(hexw==4){*q++=h[(c>>12)&15];*q++=h[(c>>8)&15];}
            *q++=h[(c>>4)&15];*q++=h[c&15];
        }
        *q++='\n';f.write(outrow,q-outrow);
    }

    // orbit-canonicity relative to the enumerated set (class 1 in S, quota ok)
    bool is_canonical(const vector<int>&S){
        uint64_t id[NW]={};
        for(int j=0;j<K;j++){
            int v=S[(size_t)j]-1;id[v>>6]|=1ULL<<(v&63);
        }
        for(int j=0;j<K;j++){
            int g=invfold[S[(size_t)j]-1];
            if(g==1)continue;
            uint16_t image[K];uint64_t image_bits[NW]={};
            for(int z=0;z<K;z++){
                int c=trans[g-1][S[(size_t)z]-1];image[z]=(uint16_t)c;
                int v=c-1;image_bits[v>>6]|=1ULL<<(v&63);
            }
            bool smaller=false;
            for(int k=0;k<NW;k++){
                uint64_t d=image_bits[k]^id[k];
                if(d){smaller=(image_bits[k]&d&-d)!=0;break;}
            }
            if(!smaller)continue;
            uint64_t no=0;bool one_seen=false;
            for(int z=0;z<K;z++){
                if(image[z]==1&&!one_seen){one_seen=true;continue;}
                no|=qm[image[z]-1];
            }
            if(__builtin_popcountll(no)<=quota_limit)return false;
        }
        return true;
    }

    uint64_t enumerate_irsubroot_raw(int ridx,int sidx,ofstream &f,bool canonOnly,uint64_t &scanned){
        auto pre=make_irprefix2(ridx);
        if(sidx<0||sidx>=(int)pre.choices.size())throw runtime_error("bad ir subroot");
        Bits rem=pre.avail;for(int j=0;j<=sidx;j++)rem.reset(pre.choices[j]);
        int vi=pre.choices[sidx];IrState ns=pre.st;if(!add_irred(ns,pre.covered,(int)pre.sel.size(),vi+1))return 0;
        uint64_t nq=pre.other|qm[vi];if(__builtin_popcountll(nq)>quota_limit)return 0;
        pre.sel.push_back(vi+1);uint64_t cnt=0;
        dfs_irred(pre.sel,pre.covered|cov[vi],rem,K-3,nq,ns,[&](const vector<int>&S){
            scanned++;leaves++;
            if(canonOnly&&!is_canonical(S))return;
            write_rawhex(f,S);cnt++;});
        return cnt;
    }

    void enumerate_irredundantK(){
        set_quota_k();
        vector<int>sel{1};Bits avail=ALL;avail.reset(0);
        IrState st;init_irstate(st);
        dfs_irred(sel,C1,avail,K-1,0ULL,st,[&](const vector<int>&S){vector<int>mm(S.size(),1);out.insert(canonical(S,mm));leaves++;});
    }

    size_t run_decomposed(){
        nodes=leaves=canon_calls=0;out.clear();
        auto c=enumerate_km1covers();
        cerr<<"p="<<p<<" canonical"<<KM1<<"="<<c.size()<<" nodes_after="<<nodes<<"\n";
        reducible_from_km1(c);
        size_t red=out.size();
        cerr<<"p="<<p<<" reducible"<<K<<"="<<red<<"\n";
        enumerate_irredundantK();
        cerr<<"p="<<p<<" union"<<K<<"="<<out.size()<<" nodes_total="<<nodes<<" leaves="<<leaves<<" canon="<<canon_calls<<"\n";
        return out.size();
    }
    size_t run_direct(){
        nodes=leaves=canon_calls=0;out.clear();quota_limit=-1;irredundant=false;coordK=K;
        for(int s=1;s<=K;s++){
            vector<int>sel{1};Bits avail=ALL;avail.reset(0);
            dfs_support(sel,C1,avail,s-1,[&](const vector<int>&S){emit_direct_support(S);});
        }
        cerr<<"p="<<p<<" direct"<<K<<"="<<out.size()<<" nodes="<<nodes<<" leaves="<<leaves<<" canon="<<canon_calls<<"\n";
        return out.size();
    }
};

int main(int argc,char**argv){
    if(argc<3){cerr<<"usage: prog p mode(irroots|irsubroots r|irsubrootrawc r s out|irsubrootraw r s out|km1roots|km1root r out|km1all out|both|direct|decomp)\n";return 2;}
    int p=atoi(argv[1]);string mode=argv[2];
    auto t0=chrono::steady_clock::now();
    Generator G(p);
    auto secs=[&](){return chrono::duration<double>(chrono::steady_clock::now()-t0).count();};
    if(mode=="info"){cout<<"K="<<K<<" p="<<p<<" n="<<G.n<<" m="<<G.m<<" VEC="<<VEC<<" NW="<<NW<<"\n";}
    else if(mode=="direct")cout<<G.run_direct()<<"\n";
    else if(mode=="km1roots"){
        int bt=-1;auto ch=G.root_choices(bt);cout<<"bt="<<bt+1<<" roots="<<ch.size()<<"\n";for(int i=0;i<(int)ch.size();i++)cout<<i<<" "<<ch[i]+1<<"\n";
    }
    else if(mode=="km1root"){
        if(argc<5){cerr<<"need root-index output-file\n";return 2;}int r=atoi(argv[3]);auto z=G.enumerate_km1root(r);vector<string>v(z.begin(),z.end());sort(v.begin(),v.end());ofstream f(argv[4]);for(auto&s:v)f<<G.hexkey(s)<<"\n";cout<<"root="<<r<<" canonical"<<KM1<<"="<<z.size()<<" nodes="<<G.nodes<<" leaves="<<G.leaves<<" secs="<<secs()<<"\n";
    }
    else if(mode=="irroots"){
        int bt=-1;auto ch=G.root_choices(bt);cout<<"bt="<<bt+1<<" roots="<<ch.size()<<"\n";for(int i=0;i<(int)ch.size();i++)cout<<i<<" "<<ch[i]+1<<"\n";
    }
    else if(mode=="irroot"){
        if(argc<5){cerr<<"need root-index output-file\n";return 2;}int r=atoi(argv[3]);auto z=G.enumerate_irroot(r);vector<string>v(z.begin(),z.end());sort(v.begin(),v.end());ofstream f(argv[4]);for(auto&s:v)f<<G.hexkey(s)<<"\n";cout<<"root="<<r<<" canonical_irred"<<K<<"="<<z.size()<<" nodes="<<G.nodes<<" leaves="<<G.leaves<<" secs="<<secs()<<"\n";
#ifdef EXP_SYM
        cout<<"depth: visits distinct_orbits visits/distinct\n";
        for(size_t d=0;d<G.sym_visits.size();d++)if(G.sym_visits[d]){cout<<"  "<<d<<": "<<G.sym_visits[d]<<" "<<G.sym_sets[d].size();if(G.sym_sets[d].size())cout<<" "<<(double)G.sym_visits[d]/G.sym_sets[d].size();cout<<"\n";}
#endif
    }
#ifdef EXP_SYM
    else if(mode=="irallsym"){
        int bt=-1;auto ch=G.root_choices(bt);uint64_t tot=0;
        for(int r=0;r<(int)ch.size();r++){auto z=G.enumerate_irroot(r);tot+=z.size();cout<<"root="<<r<<" canon="<<z.size()<<" nodes="<<G.nodes<<" leaves="<<G.leaves<<"\n";}
        cout<<"ALL canonical_irred"<<K<<"="<<tot<<" nodes="<<G.nodes<<" leaves="<<G.leaves<<" secs="<<secs()<<"\n";
        cout<<"depth: visits distinct_orbits visits/distinct (across all roots)\n";
        for(size_t d=0;d<G.sym_visits.size();d++)if(G.sym_visits[d]){cout<<"  "<<d<<": "<<G.sym_visits[d]<<" "<<G.sym_sets[d].size();if(G.sym_sets[d].size())cout<<" "<<(double)G.sym_visits[d]/G.sym_sets[d].size();cout<<"\n";}
    }
#endif
    else if(mode=="irsubroots"){
        if(argc<4){cerr<<"need root-index\n";return 2;}int r=atoi(argv[3]);auto pre=G.make_irprefix2(r);cout<<"root="<<r<<" bt="<<pre.bt+1<<" subroots="<<pre.choices.size()<<"\n";for(int i=0;i<(int)pre.choices.size();i++)cout<<i<<" "<<pre.choices[i]+1<<"\n";
    }
    else if(mode=="irsubrootraw"||mode=="irsubrootrawc"){
        if(argc<6){cerr<<"need root-index sub-index output-file\n";return 2;}int r=atoi(argv[3]),ss=atoi(argv[4]);ofstream f(argv[5]);uint64_t scanned=0;auto cnt=G.enumerate_irsubroot_raw(r,ss,f,mode=="irsubrootrawc",scanned);f.close();cout<<"root="<<r<<" sub="<<ss<<" raw_irred"<<K<<"="<<cnt<<" rawscan="<<scanned<<" nodes="<<G.nodes<<" leaves="<<G.leaves<<" secs="<<secs()<<"\n";
    }
    else if(mode=="km1low"){
        if(argc<5){cerr<<"need smax output-file\n";return 2;}int smax=atoi(argv[3]);auto c=G.enumerate_km1covers_upto(smax);vector<string>v(c.begin(),c.end());sort(v.begin(),v.end());ofstream f(argv[4]);for(auto&x:v)f<<G.hexkey(x)<<"\n";   // hex rows: same format as km1root, consumable by the extension cascade
        cout<<"km1low smax="<<smax<<" canonical="<<c.size()<<" nodes="<<G.nodes<<" leaves="<<G.leaves<<" secs="<<secs()<<"\n";
    }
    else if(mode=="km1all"){
        if(argc<4){cerr<<"need output-file\n";return 2;}
        auto c=G.enumerate_km1covers();
        vector<string>v(c.begin(),c.end());sort(v.begin(),v.end());
        ofstream f(argv[3]);for(auto&s:v)f<<G.hexkey(s)<<"\n";
        cout<<"canonical"<<KM1<<"="<<c.size()<<" nodes="<<G.nodes<<" leaves="<<G.leaves<<" secs="<<secs()<<"\n";
    }
#ifdef EXP_MINPRIV
    else if(mode=="minprivtest"){
        // optional 4th arg: root index (default: all roots)
        int only=(argc>=4)?atoi(argv[3]):-1;
        int bt=-1;auto ch=G.root_choices(bt);
        auto run=[&](bool quota){
            G.nodes=0;G.leaves=0;G.nodes_cut=0;G.leaves_all=0;G.leaves_minpriv=0;G.orb_all.clear();G.orb_minpriv.clear();
            for(int r=0;r<(int)ch.size();r++){
                if(only>=0&&r!=only)continue;
                G.set_quota_k();if(!quota)G.quota_limit=G.m;
                Bits avail=G.ALL;avail.reset(0);Bits rem=avail;for(int j=0;j<=r;j++)rem.reset(ch[j]);
                int vi=ch[r];vector<int>sel{1};
                Generator::IrState st;G.init_irstate(st);
                if(!G.add_irred(st,G.C1,(int)sel.size(),vi+1))continue;
                sel.push_back(vi+1);
                G.dfs_irred(sel,G.C1|G.cov[vi],rem,K-2,G.qm[vi],st,[&](const vector<int>&S){G.leaves++;G.minpriv_leaf(S);});
            }
            cout<<(quota?"QUOTA-ON ":"QUOTA-OFF")<<" nodes="<<G.nodes<<" nodes_in_dual_cut_subtrees="<<G.nodes_cut
                <<" leaves="<<G.leaves<<" leaves_minpriv="<<G.leaves_minpriv
                <<" orbits_all="<<G.orb_all.size()<<" orbits_minpriv="<<G.orb_minpriv.size()<<" secs="<<secs()<<"\n";
            cout<<"slots: nodes / pruned_by_gain / dual_prunes_extra\n";
            for(int s=K;s>=0;s--)if(G.nodes_by_slots[s])cout<<"  "<<s<<": "<<G.nodes_by_slots[s]<<" / "<<G.pruned_gain[s]<<" / "<<G.dual_would_prune[s]<<"\n";
            memset(G.nodes_by_slots,0,sizeof G.nodes_by_slots);memset(G.pruned_gain,0,sizeof G.pruned_gain);memset(G.dual_would_prune,0,sizeof G.dual_would_prune);
        };
        run(true);if(only<0)run(false);
    }
#endif
    else if(mode=="decomp")cout<<G.run_decomposed()<<"\n";
    else if(mode=="both"){
        size_t a=G.run_direct();auto A=G.out;
        size_t b=G.run_decomposed();
        size_t miss=0,extra=0;for(auto&s:A)if(!G.out.count(s))miss++;for(auto&s:G.out)if(!A.count(s))extra++;
        cout<<"direct="<<a<<" decomp="<<b<<" missing="<<miss<<" extra="<<extra<<" secs="<<secs()<<"\n";
    } else return 2;
}
