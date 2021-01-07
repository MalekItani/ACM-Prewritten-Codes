// PRIMES
ll A = 100013;
ll B = 1000000007;
ll C = 989999993;
ll D = 999999929;

// SIEVE + prime factor decomposition
bool prime[1000005];
int spf[1000005];
vector<int> sieve(int n){
	vector<int> v;
    memset(prime, true, sizeof(prime));
    spf[0] = 1;
    spf[1] = 1;
    for (int p=2; p*p<=n; p++){
        if (prime[p] == true){
            for (int i=p*p; i<=n; i += p){
                if(prime[i]) spf[i] = p;
                prime[i] = false;
            }
        }
    }

    for (int p=2; p<=n; p++) if (prime[p]){
      v.push_back(p);
      spf[p] = p;
    }
	return v;
}

map<int, int> pfd(int x){
    map<int, int> res;
    while(x > 1){
        res[spf[x]]++;
        x/=spf[x];
    }
    return res;
}

// EULER TOTIENT FUNCTION:
// tot(p) = Number of integers k such that gcd(k, p) = 1
// Notes: ( a ^ n ) mod p = ( a ^ (tot(p) + (n mod tot(p)) ) ) mod p
// tot(tot( ... tot(p))) = 1 in ~ O(logn) steps
ll tot(ll n) {
    ll result = n;
    for (ll i = 2; i * i <= n; i++) {
        if(n % i == 0) {
            while(n % i == 0)
                n /= i;
            result -= result / i;
        }
    }
    if(n > 1)
        result -= result / n;
    return result;
}


// MODULAR INVERSE
ll modInverse(ll a, ll m){
    ll m0 = m;
    ll y = 0, x = 1;
    if (m == 1) return 0;
    while (a > 1){
        ll q = a / m; ll t = m;
        m = a % m, a = t;
        t = y; y = x - q * y; x = t;
    }
    if (x < 0) x += m0;
    return x;
}

// MODULAR INT - https://codeforces.com/contest/1336/submission/76882437

#include <cstdint>
template <std::uint_fast64_t Modulus> class modint {
    using u64 = std::uint_fast64_t;

  public:
    u64 a;
    constexpr modint(const u64 x = 0) noexcept : a(x % Modulus) {}
    constexpr u64 &value() noexcept { return a; }
    constexpr const u64 &value() const noexcept { return a; }
    constexpr modint operator+(const modint rhs) const noexcept { return modint(*this) += rhs; }
    constexpr modint operator-(const modint rhs) const noexcept { return modint(*this) -= rhs; }
    constexpr modint operator*(const modint rhs) const noexcept { return modint(*this) *= rhs; }
    constexpr modint operator/(const modint rhs) const noexcept { return modint(*this) /= rhs; }
    constexpr modint &operator+=(const modint rhs) noexcept {
        a += rhs.a;
        if(a >= Modulus) {
            a -= Modulus;
        }
        return *this;
    }
    constexpr modint &operator-=(const modint rhs) noexcept {
        if(a < rhs.a) {
            a += Modulus;
        }
        a -= rhs.a;
        return *this;
    }
    constexpr modint &operator*=(const modint rhs) noexcept {
        a = a * rhs.a % Modulus;
        return *this;
    }
    constexpr modint &operator/=(modint rhs) noexcept {
        u64 exp = Modulus - 2;
        while(exp) {
            if(exp % 2) {
                *this *= rhs;
            }
            rhs *= rhs;
            exp /= 2;
        }
        return *this;
    }
};
using mint = modint<MOD>;


// LCM MODULO mod
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ START ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
#define MAX 10000003
ll prime[MAX];


// Function to return a^n
ll power(ll a, ll n){
    if (n == 0) return 1;
    ll p = power(a, n / 2) % mod;
    p = (p * p) % mod;
    if (n & 1) p = (p * a) % mod;
    return p;
}

// Function to find the smallest prime factors
// of numbers upto MAX
void sieve(){
    prime[0] = prime[1] = 1;
    for (ll i = 2; i < MAX; i++) {
        if (prime[i] == 0) {
            for (ll j = i * 2; j < MAX; j += i) if (prime[j] == 0) prime[j] = i;
            prime[i] = i;
        }
    }
}

// Function to return the LCM modulo M
ll lcmModuloM(vi ar){
  unordered_map<ll, ll> max_map;
  if(sz(ar) == 0) return 1;
    for (ll i = 0; i < sz(ar); i++) {
        ll num = ar[i];
        unordered_map<ll, ll> temp;
        while (num > 1) {
            ll factor = prime[num];
            temp[factor]++;
            num /= factor;
        }
        for (auto it : temp) max_map[it.first] = max(max_map[it.first], it.second);
    }
    ll ans = 1;
    for (auto it : max_map) ans = (ans * power(it.first, it.second)) % mod;
    return ans;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ END ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
