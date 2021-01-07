#include <bits/stdc++.h>
#define ll long long
#define ld long double
#define pi pair<ll, ll>
#define mp make_pair
#define pb push_back
#define vi vector<ll>
#define pd pair<double, double>
#define vp vector<pi>
#define f(i, n) for(int i = 0; i < n; i++)
#define fo(i, a, n) for(int i = a; i < n; i++)
#define mod (1000*1000*1000+7)
#define sz(x) int((x).size())
#define all(x) x.begin(),x.end()
#define DBG(v) cerr << #v << " = " << (v) << endl;
ll INF = (1ll<<62);

using namespace std;

// INCLUSION-EXCLUSION
// Finds the number of numbers from 1 to m (inclusive) that are
// divisible by any of the primes in a list "a".
ll find_mults(vi &a, ll m){
    ll n = a.size();
    ll odd = 0, even = 0;
    ll counter, j, p = 1;
    ll pow_set_size = (1 << n);
    for (counter = 1; counter < pow_set_size; counter++){
        p = 1;
        for (j = 0; j < n; j++) if (counter & (1 << j)) p *= a[j];
        if (__builtin_popcount(counter) & 1) odd += ( m / p );
        else even += ( m / p );
    }
    return odd - even;
}

// SLIDING WINDOW MINIMUM/MAXIMUM
template<class T, class R>
void sw(T arr[], R res[], int n, int window_size){
	deque<int> dq;
    int x = 0;
	for(int i = 0; i < n; ++i) {
		while(dq.size() && arr[dq.back()] > arr[i]) dq.pop_back();
		dq.push_back(i);
		if(i >= window_size-1) {
            res[x++] = arr[dq.front()];
			if(dq.front() <= i-(window_size-1)) dq.pop_front();
		}
	}
}

// Hack-proof hashing
struct custom_hash {
    static ll splitmix64(ll x) {
        x += 0x9e3779b97f4a7c15;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
        x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
        return x ^ (x >> 31);
    }

    ll operator()(ll x) const {
        static const ll FIXED_RANDOM = chrono::steady_clock::now().time_since_epoch().count();
        return splitmix64(x + FIXED_RANDOM);
    }
};

// Ternary Search: MUST TAKE MIN AT THE END
ll i = 0, j = 1e10;
while(i+2 < j){
	ll m1 = i + (j-i)/3;
	ll m2 = j - (j-i)/3;
	if(check(m1) < check(m2)){
		j = m2;
	}else{
		i = m1;
	}
}
ll ans = min({check(i), check(i+1), check(i+2)});


// 4D prefix sum (0 based)
ans[i+1][j+1][k+1][l+1] = ans[i][j+1][k+1][l+1]
			+ans[i+1][j][k+1][l+1]
			+ans[i+1][j+1][k][l+1]
			+ans[i+1][j+1][k+1][l]
			-ans[i][j][k+1][l+1]
			-ans[i][j+1][k][l+1]
			-ans[i][j+1][k+1][l]
			-ans[i+1][j][k][l+1]
			-ans[i+1][j][k+1][l]
			-ans[i+1][j+1][k][l]
			+ans[i+1][j][k][l]
			+ans[i][j+1][k][l]
			+ans[i][j][k+1][l]
			+ans[i][j][k][l+1]
			-ans[i][j][k][l]
			+pg2[i+1][j+1][k+1][l+1];

cout << ans[c+1][d+1][c+1][d+1]
	-ans[a][d+1][c+1][d+1]
	-ans[c+1][b][c+1][d+1]
	-ans[c+1][d+1][a][d+1]
	-ans[c+1][d+1][c+1][b]
	+ans[a][b][c+1][d+1]
	+ans[a][d+1][a][d+1]
	+ans[a][d+1][c+1][b]
	+ans[c+1][b][a][d+1]
	+ans[c+1][b][c+1][b]
	+ans[c+1][d+1][a][b]
	-ans[c+1][b][a][b]
	-ans[a][d+1][a][b]
	-ans[a][b][c+1][b]
	-ans[a][b][a][d+1]
	+ans[a][b][a][b] << endl;



// BURGERS - CHARBEL
ll solve_burgers(vi nums){
	sort(nums.begin(), nums.end());
	ll res = 0;

	if (nums[0] <= nums[2] - nums[1]) {
		return nums[0] + nums[1];
	}

	res += nums[2] - nums[1];
	nums[0] -= (nums[2] - nums[1]);
	nums[1] -= (nums[0]+1)/2;
	res += nums[0] + nums[1];

	return res;
}
