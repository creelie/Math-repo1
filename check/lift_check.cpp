// lift_check.cpp - independent re-implementation of the binary lift cascade
// (written from Definition 2.1 of Sungkawichai-Trakulthongchai, not from the
// Allikvere sources). For a base tuple v in Z_p^K it counts, at levels
// L = 2, 4, 8, 16, 32, the lifts w in Z_{Lp}^K (w = v mod p) that are
// (K,p,L)-improper: no time T in Z_{Lp} has (K+1)*dist(T*w_i, Lp) >= Lp for
// every i, and fewer than K-1 coordinates are even (the gcd clause for the
// only prime 2 dividing L). Every time T = 1..Lp-1 is checked, both parities.
// A lift proper at level L stays proper at 2L (its witness T gives 2T), so
// level-2L candidates are the children w + a*Lp of improper level-L lifts.
// Usage: lift_check P v1 ... vK      (one base)
//        lift_check P -             (bases from stdin, one per line)
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#ifndef K
#define K 15
#endif
typedef long long ll;
static const int NA = 1 << K, NWD = NA / 64;
static int P;
static uint64_t BITMASK[K][2][NWD];   // BITMASK[i][b] = set of a with a_i == b

static inline bool bad(ll T, ll w, ll N) {
    ll y = (T % N) * (((w % N) + N) % N) % N;
    ll d = y < N - y ? y : N - y;
    return (ll)(K + 1) * d < N;
}

// improper children of parent w (improper at level M): w_i + a_i*M mod 2M
static void children(const std::vector<ll>& w, ll M, std::vector<int>& out) {
    ll N = 2 * M;
    static uint64_t imp[NWD], cov[NWD];
    for (int j = 0; j < NWD; j++) imp[j] = ~0ULL;
    for (ll T = 1; T < N; T++) {
        bool all = false, any = false;
        for (int j = 0; j < NWD; j++) cov[j] = 0;
        for (int i = 0; i < K; i++) {
            bool b0 = bad(T, w[i], N), b1 = bad(T, w[i] + M, N);
            if (b0 && b1) { all = true; break; }
            if (!b0 && !b1) continue;
            const uint64_t* m = BITMASK[i][b1 ? 1 : 0];
            for (int j = 0; j < NWD; j++) cov[j] |= m[j];
            any = true;
        }
        if (all) continue;                 // T covered for every lift
        if (!any) { for (int j = 0; j < NWD; j++) imp[j] = 0; break; }   // T is a witness for every lift
        uint64_t nz = 0;
        for (int j = 0; j < NWD; j++) { imp[j] &= cov[j]; nz |= imp[j]; }
        if (!nz) break;                    // every lift already has a witness
    }
    out.clear();
    for (int a = 0; a < NA; a++) if ((imp[a >> 6] >> (a & 63)) & 1) {
        int even = 0;
        for (int i = 0; i < K; i++) { ll x = w[i] + (((a >> i) & 1) ? M : 0); if (x % 2 == 0) even++; }
        if (even >= K - 1) continue;       // gcd clause: proper
        out.push_back(a);
    }
}

static ll cnt[6];
static void rec(const std::vector<ll>& w, ll M, int lev) {   // w improper at level M = 2^lev * P
    if (lev == 5) return;
    std::vector<int> ch; children(w, M, ch);
    cnt[lev + 1] += (ll)ch.size();
    std::vector<ll> x(K);
    for (int a : ch) {
        for (int i = 0; i < K; i++) x[i] = w[i] + (((a >> i) & 1) ? M : 0);
        rec(x, 2 * M, lev + 1);
    }
}

static void run(const std::vector<ll>& v) {
    for (int l = 0; l < 6; l++) cnt[l] = 0;
    bool cover = true;
    for (ll T = 1; T < P && cover; T++) {
        bool c = false;
        for (int i = 0; i < K; i++) if (bad(T, v[i], P)) { c = true; break; }
        if (!c) cover = false;
    }
    for (int i = 0; i < K; i++) printf("%lld ", v[i]);
    printf("| level1 %s;", cover ? "improper" : "PROPER");
    if (cover) rec(v, P, 0);
    for (int l = 1; l <= 5; l++) printf(" l%d=%lld", 1 << l, cnt[l]);
    printf("%s\n", cnt[5] ? " PERSISTENT" : "");
    fflush(stdout);
}

int main(int argc, char** argv) {
    for (int i = 0; i < K; i++) for (int a = 0; a < NA; a++) {
        int b = (a >> i) & 1;
        BITMASK[i][b][a >> 6] |= 1ULL << (a & 63);
    }
    if (argc < 3) { fprintf(stderr, "usage: lift_check P v1..v%d | lift_check P -\n", K); return 2; }
    P = atoi(argv[1]);
    std::vector<ll> v(K);
    if (argv[2][0] == '-') {
        while (true) {
            int ok = 1;
            for (int i = 0; i < K && ok; i++) ok = (scanf("%lld", &v[i]) == 1);
            if (!ok) break;
            run(v);
        }
        return 0;
    }
    if (argc != K + 2) { fprintf(stderr, "need %d values\n", K); return 2; }
    for (int i = 0; i < K; i++) v[i] = atoll(argv[2 + i]);
    run(v);
}
