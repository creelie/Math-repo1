// lift_check_sym.cpp - independent recount of the binary lift cascade to any
// binary level, with the shift lemma, written from the definitions (it shares
// no code with src/cascade_filter_k.cpp).
//
// For a level-one cover v of Z_p by K speed classes it walks the lifts
// w in Z_{Lp}^K (w = v mod p) for L = 2, 4, ..., 2^MAXLEV. A lift is settled
// when some time T in 1..Lp-1 (both parities, checked one by one) has
// (K+1)*dist(T*w_i, Lp) >= Lp for every i, when at least K-1 coordinates are
// even (gcd clause), or when for some d = 2^j dividing L, j >= 2, the
// coordinates not divisible by d forbid fewer than d of the shifts q/d in
// total (the number forbidden by w_i is the largest number of the d points
// q*w_i/d mod 1 in an open arc of length 2/(K+1), found by brute force).
// Settled lifts stay settled at higher levels, so only the children of
// unsettled lifts are examined.
//
// Coordinates with equal values are interchangeable. Among the children
// w + a*M of a node, one per orbit of the node's stabiliser is expanded: the
// one whose bits are non-decreasing (zeros first) along every group of equal
// coordinates, with weight prod C(g, ones in the group). (The cascade uses the
// opposite representative, ones first.) The printed counts are the exact
// numbers of labelled unsettled lifts.
//
// Usage: lift_check_sym MAXLEV P v1 ... vK
//        lift_check_sym MAXLEV P -        (bases from stdin, one per line)
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#ifndef K
#define K 15
#endif
typedef long long ll;
typedef unsigned __int128 u128;
static const int NA = 1 << K, NWD = NA / 64;
static int P, MAXLEV;
static uint64_t BITMASK[K][2][NWD];   // BITMASK[i][b] = set of a with a_i == b
static ll BIN[K + 1][K + 1];
static u128 cnt[32], shiftcnt;

static inline bool bad(ll T, ll w, ll N) {
    ll y = (T % N) * (((w % N) + N) % N) % N;
    ll d = y < N - y ? y : N - y;
    return (ll)(K + 1) * d < N;
}

static int max_in_arc(ll x, int d) {
    std::vector<int> r(d);
    for (int q = 0; q < d; q++) r[q] = (int)(((ll)q * (x % d)) % d);
    int best = 0;
    for (int s = 0; s < d; s++) {
        int c = 0;
        for (int q = 0; q < d; q++) { int diff = ((r[q] - r[s]) % d + d) % d; if ((ll)(K + 1) * diff < 2LL * d) c++; }
        if (c > best) best = c;
    }
    return best;
}

static bool shift_settled(const std::vector<ll>& x, ll L) {
    for (ll d = 4; d <= L && L % d == 0; d *= 2) {
        int forb = 0;
        for (int i = 0; i < K; i++) if (x[i] % d != 0) forb += max_in_arc(x[i], (int)d);
        if (forb < d) return true;
    }
    return false;
}

// unsettled-by-witness children a of the node w at level M (bit i adds M to w_i)
static void witness_free_children(const std::vector<ll>& w, ll M, std::vector<int>& out) {
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
        if (all) continue;
        if (!any) { for (int j = 0; j < NWD; j++) imp[j] = 0; break; }
        uint64_t nz = 0;
        for (int j = 0; j < NWD; j++) { imp[j] &= cov[j]; nz |= imp[j]; }
        if (!nz) break;
    }
    out.clear();
    for (int a = 0; a < NA; a++) if ((imp[a >> 6] >> (a & 63)) & 1) out.push_back(a);
}

static void rec(const std::vector<ll>& w, ll M, int lev, u128 weight) {   // w unsettled at level M = 2^lev * P
    if (lev == MAXLEV) return;
    // groups of equal coordinates, in position order
    std::vector<int> grp(K, -1);
    int ng = 0;
    for (int i = 0; i < K; i++) {
        if (grp[i] >= 0) continue;
        for (int j = i; j < K; j++) if (w[j] == w[i]) grp[j] = ng;
        ng++;
    }
    std::vector<int> ch;
    witness_free_children(w, M, ch);
    std::vector<ll> x(K);
    for (int a : ch) {
        // orbit representative: along each group the bits read 0...01...1
        bool canon = true;
        std::vector<int> seenOne(ng, 0), size(ng, 0), ones(ng, 0);
        for (int i = 0; i < K && canon; i++) {
            int g = grp[i], b = (a >> i) & 1;
            if (b) { seenOne[g] = 1; ones[g]++; }
            else if (seenOne[g]) canon = false;
            size[g]++;
        }
        if (!canon) continue;
        ll wt = 1;
        for (int g = 0; g < ng; g++) wt *= BIN[size[g]][ones[g]];
        int even = 0;
        for (int i = 0; i < K; i++) { x[i] = w[i] + (((a >> i) & 1) ? M : 0); if (x[i] % 2 == 0) even++; }
        if (even >= K - 1) continue;
        if (shift_settled(x, 2 * M / P)) { shiftcnt += weight * (u128)wt; continue; }
        cnt[lev + 1] += weight * (u128)wt;
        rec(x, 2 * M, lev + 1, weight * (u128)wt);
    }
}

static void pr(u128 x) {
    char b[48]; int n = 0;
    do { b[n++] = (char)('0' + (int)(x % 10)); x /= 10; } while (x);
    while (n) putchar(b[--n]);
}

static void run(const std::vector<ll>& v) {
    for (int l = 0; l < 32; l++) cnt[l] = 0;
    shiftcnt = 0;
    bool cover = true;
    for (ll T = 1; T < P && cover; T++) {
        bool c = false;
        for (int i = 0; i < K; i++) if (bad(T, v[i], P)) { c = true; break; }
        if (!c) cover = false;
    }
    for (int i = 0; i < K; i++) printf("%lld ", v[i]);
    printf("| level1 %s;", cover ? "improper" : "PROPER");
    if (cover) rec(v, P, 0, 1);
    int last = 0;
    for (int l = 1; l <= MAXLEV; l++) { printf(" l%d=", 1 << l); pr(cnt[l]); last = l; if (cnt[l] == 0) break; }
    printf(" shift="); pr(shiftcnt);
    printf("%s\n", (cover && cnt[last] != 0) ? " PERSISTENT" : "");
    fflush(stdout);
}

int main(int argc, char** argv) {
    for (int i = 0; i < K; i++) for (int a = 0; a < NA; a++) {
        int b = (a >> i) & 1;
        BITMASK[i][b][a >> 6] |= 1ULL << (a & 63);
    }
    for (int n = 0; n <= K; n++) { BIN[n][0] = 1; for (int k = 1; k <= n; k++) BIN[n][k] = BIN[n - 1][k - 1] + (k < n ? BIN[n - 1][k] : 0); }
    if (argc < 4) { fprintf(stderr, "usage: lift_check_sym MAXLEV P v1..v%d | lift_check_sym MAXLEV P -\n", K); return 2; }
    MAXLEV = atoi(argv[1]);
    P = atoi(argv[2]);
    if (MAXLEV < 1 || MAXLEV > 12 || P < 3 || P % 2 == 0) { fprintf(stderr, "bad MAXLEV or P\n"); return 2; }
    std::vector<ll> v(K);
    if (argv[3][0] == '-') {
        while (true) {
            int ok = 1;
            for (int i = 0; i < K && ok; i++) ok = (scanf("%lld", &v[i]) == 1);
            if (!ok) break;
            run(v);
        }
        return 0;
    }
    if (argc != K + 3) { fprintf(stderr, "need %d values\n", K); return 2; }
    for (int i = 0; i < K; i++) v[i] = atoll(argv[3 + i]);
    run(v);
}
