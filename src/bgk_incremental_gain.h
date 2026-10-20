#ifndef BGK_INCREMENTAL_GAIN_H
#define BGK_INCREMENTAL_GAIN_H
#include <cstdint>
#include <stdexcept>
#include <vector>
#include <immintrin.h>
#ifndef __AVX2__
#error "BGK_INCREMENTAL_GAIN requires AVX2; compile both comparison binaries with the same target flags."
#endif

// Exact incidence sums only: no filtering, ordering, or tree changes.
// U_child = U_parent minus R, R = U_parent intersect C(chosen).
// Thus g_child(v) = g_parent(v) - |R intersect C(v)| for EVERY class v.
template<int Capacity>
class BgkIncrementalGain {
    int n_, stride_;
    std::vector<uint8_t> rows_;
    template<int Tile>
    void tile(const uint8_t*const*selected,int count,int offset,uint8_t*out
#ifdef BGK_FUSED_DELTA
              ,const uint8_t*parent
#endif
              )const {
        __m256i acc[Tile];
        for(int j=0;j<Tile;++j)acc[j]=_mm256_setzero_si256();
        for(int i=0;i<count;++i){
            const uint8_t*row=selected[i]+offset;
            #pragma GCC unroll 8
            for(int j=0;j<Tile;++j)acc[j]=_mm256_add_epi8(acc[j],
                _mm256_loadu_si256(reinterpret_cast<const __m256i*>(row+32*j)));
        }
        for(int j=0;j<Tile;++j){
#ifdef BGK_FUSED_DELTA
            if(parent)acc[j]=_mm256_sub_epi8(_mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(parent+offset+32*j)),acc[j]);
#endif
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(out+offset+32*j),acc[j]);
        }
    }
    template<class Bits>
    void sum(const Bits&points,uint8_t*out
#ifdef BGK_FUSED_DELTA
             ,const uint8_t*parent=nullptr
#endif
             )const {
        const uint8_t* selected[Capacity];int count=0;
        for(int base=0;base<n_;base+=64){
            uint64_t w=points.w[base/64];
            if(n_-base<64)w&=(uint64_t{1}<<(n_-base))-1;
            while(w){
                int t=base+__builtin_ctzll(w);w&=w-1;
                selected[count++]=rows_.data()+std::size_t(t)*stride_;
            }
        }
        int offset=0;
#ifdef BGK_FUSED_DELTA
        for(;offset+128<=stride_;offset+=128)tile<4>(selected,count,offset,out,parent);
        if(offset+64<=stride_){tile<2>(selected,count,offset,out,parent);offset+=64;}
        if(offset<stride_)tile<1>(selected,count,offset,out,parent);
#else
        for(;offset+128<=stride_;offset+=128)tile<4>(selected,count,offset,out);
        if(offset+64<=stride_){tile<2>(selected,count,offset,out);offset+=64;}
        if(offset<stride_)tile<1>(selected,count,offset,out);
#endif
    }
public:
    template<class Generator>
    explicit BgkIncrementalGain(const Generator&g)
      : n_(g.n),stride_((g.n+31)&~31),rows_(std::size_t(g.n)*stride_,0) {
        static_assert(Capacity>0,"positive capacity");
        if(n_<=0||n_>Capacity||g.m>64)throw std::runtime_error("gain kernel capacity");
        for(int v=0;v<n_;++v){
            unsigned column=0;
            for(int t=0;t<n_;++t)if(g.cov[v].test(t)){
                rows_[std::size_t(t)*stride_+v]=1;++column;
            }
            if(column!=unsigned(g.m)||column>64)throw std::runtime_error("gain incidence column");
        }
    }
    int stride()const{return stride_;}
    template<class Bits>
    void compute(const Bits&unc,int need,uint8_t*out,
                 const uint8_t*parent,const Bits*removed)const {
        // The parent stores ALL columns, even for unavailable candidates.
        // Histogram/sibling removals never modify that array. Both parent
        // and removed are caller-owned and live through this synchronous call.
        if(parent&&removed&&removed->count()<need){
#ifdef BGK_FUSED_DELTA
            sum(*removed,out,parent);
#else
            sum(*removed,out);
            for(int j=0;j<stride_;j+=32){
                __m256i old=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(parent+j));
                __m256i delta=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(out+j));
                // R subset U_parent => 0 <= delta <= old <= m <= 64.
                _mm256_storeu_si256(reinterpret_cast<__m256i*>(out+j),_mm256_sub_epi8(old,delta));
            }
#endif
        }else sum(unc,out); // Fresh root or a cheaper direct recomputation.
    }
};
#endif
