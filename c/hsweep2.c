/* Independent re-implementation of the hull-unit sweep (written without reusing hsweep.c).
 *
 * For a level lam, let S'(lam) = {(m,k) in Z^2 : m,k >= 1, m^2 + k^2 <= lam} and let
 * H(lam) be the number of primitive lattice units on the non-axis-parallel edges of the
 * upper-right boundary chain of conv S'(lam).  The paper proves: if H >= 2 r - 3, where
 * r(lam) = #{(m,k) in Z_{>=1}^2 : m^2 + k^2 = lam}, then the extended Courant property
 * holds at lam.  This program checks that inequality at every eigenvalue lam <= X.
 *
 * Differences from hsweep.c, on purpose:
 *   - r(lam) comes from the factorisation of lam (smallest-prime-factor sieve) and the
 *     formula r = prod_{p = 1 mod 4} (e_p + 1) - [lam is a square], valid when every
 *     p = 3 mod 4 has even exponent (r = 0 otherwise); hsweep.c counts pairs directly;
 *   - column heights by a monotone pointer, not by a square root per column;
 *   - the hull is built from right to left (m = M down to 1).
 * Output: per-decade minima of the slack ceil(H/2) - (r - 1), the number of failures, and
 * checksums sum(H), sum(r) over all eigenvalues, for comparison with other programs.
 *
 * Build: cc -O2 -fopenmp -o hsweep2 hsweep2.c -lm
 * Usage: ./hsweep2 X        or     ./hsweep2 -l lam1 lam2 ...
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

typedef long long ll;

static ll gcd(ll a, ll b) { if (a < 0) a = -a; if (b < 0) b = -b; while (b) { ll t = a % b; a = b; b = t; } return a; }

static ll isqrt_floor(ll n) {           /* integer square root by Newton, no floating point */
    if (n < 2) return n < 0 ? -1 : n;
    ll x = n, y = (x + 1) / 2;
    while (y < x) { x = y; y = (x + n / x) / 2; }
    return x;
}

/* H(lam); buf must hold isqrt(lam)+2 entries of each array */
static ll hull_units(ll lam, ll *px, ll *py) {
    ll M = isqrt_floor(lam - 1);        /* largest m with m^2 + 1 <= lam */
    if (M < 1) return 0;
    int n = 0;
    ll k = 1;
    /* walk m from M down to 1; column height c(m) grows as m decreases */
    for (ll m = M; m >= 1; m--) {
        while (m * m + (k + 1) * (k + 1) <= lam) k++;
        /* point (m, k); keep the chain strictly convex (upper-right hull, right to left) */
        while (n >= 2) {
            ll ax = px[n - 1] - px[n - 2], ay = py[n - 1] - py[n - 2];
            ll bx = m - px[n - 2],        by = k - py[n - 2];
            /* going leftwards along the upper hull we must turn clockwise... for a chain
               traversed with decreasing x, an upper hull turns left (cross > 0) */
            if (ax * by - ay * bx <= 0) n--; else break;
        }
        px[n] = m; py[n] = k; n++;
    }
    ll H = 0;
    for (int i = 1; i < n; i++) {
        ll dx = px[i] - px[i - 1], dy = py[i] - py[i - 1];
        if (dx != 0 && dy != 0) H += gcd(dx, dy);
    }
    return H;
}

static uint32_t *spf;                    /* smallest prime factor */

static int mult_r(ll lam) {
    ll n = lam; ll prod = 1;
    while (n > 1) {
        ll p = spf[n]; int e = 0;
        while (n % p == 0) { n /= p; e++; }
        if (p % 4 == 3 && (e & 1)) return 0;
        if (p % 4 == 1) prod *= (e + 1);
    }
    ll s = isqrt_floor(lam);
    return (int)(prod - (s * s == lam ? 1 : 0));
}

static void sieve(ll X) {
    spf = calloc(X + 1, sizeof(uint32_t));
    for (ll i = 2; i <= X; i++) if (!spf[i]) for (ll j = i; j <= X; j += i) if (!spf[j]) spf[j] = (uint32_t)i;
}

int main(int argc, char **argv) {
    if (argc > 2 && !strcmp(argv[1], "-l")) {
        ll mx = 0; for (int i = 2; i < argc; i++) { ll v = atoll(argv[i]); if (v > mx) mx = v; }
        sieve(mx);
        ll *px = malloc(sizeof(ll) * (isqrt_floor(mx) + 4)), *py = malloc(sizeof(ll) * (isqrt_floor(mx) + 4));
        for (int i = 2; i < argc; i++) {
            ll lam = atoll(argv[i]); ll H = hull_units(lam, px, py); int r = mult_r(lam);
            printf("%lld H=%lld r=%d slack=%lld\n", lam, H, r, (H + 1) / 2 - (r - 1));
        }
        return 0;
    }
    ll X = atoll(argv[1]);
    sieve(X);
    ll minslack[20], minlam[20], nlev[20], fails = 0, sumH = 0, sumr = 0;
    for (int d = 0; d < 20; d++) { minslack[d] = (ll)1 << 60; minlam[d] = 0; nlev[d] = 0; }
    #pragma omp parallel
    {
        ll *px = malloc(sizeof(ll) * (isqrt_floor(X) + 4)), *py = malloc(sizeof(ll) * (isqrt_floor(X) + 4));
        ll ms[20], ml[20], nl[20], f = 0, sH = 0, sr = 0;
        for (int d = 0; d < 20; d++) { ms[d] = (ll)1 << 60; ml[d] = 0; nl[d] = 0; }
        #pragma omp for schedule(dynamic, 8192)
        for (ll lam = 2; lam <= X; lam++) {
            int r = mult_r(lam);
            if (r == 0) continue;
            ll H = hull_units(lam, px, py);
            int d = 0; for (ll t = lam; t >= 10; t /= 10) d++;
            ll s = (H + 1) / 2 - (r - 1);
            nl[d]++; sH += H; sr += r;
            if (s < ms[d] || (s == ms[d] && lam < ml[d])) { ms[d] = s; ml[d] = lam; }
            if (s < 0) { f++;
                #pragma omp critical
                printf("FAIL lam=%lld r=%d H=%lld\n", lam, r, H);
            }
        }
        #pragma omp critical
        {
            for (int d = 0; d < 20; d++) {
                if (ms[d] < minslack[d] || (ms[d] == minslack[d] && ml[d] < minlam[d])) { minslack[d] = ms[d]; minlam[d] = ml[d]; }
                nlev[d] += nl[d];
            }
            fails += f; sumH += sH; sumr += sr;
        }
    }
    printf("X=%lld fails=%lld sumH=%lld sumr=%lld\n", X, fails, sumH, sumr);
    for (int d = 0; d < 20; d++) if (nlev[d])
        printf("[1e%d,1e%d): levels=%lld min slack=%lld at %lld\n", d, d + 1, nlev[d], minslack[d], minlam[d]);
    return 0;
}
