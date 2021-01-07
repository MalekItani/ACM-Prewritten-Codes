// CAREFUL MAKE SURE TO CHANGE N IF MATRIX SIZE > 100 !!!!!
// ALSO BE CAREFUL OF INT/LONG LONG
namespace mat{
	int n; // Actual matrix size. Use mat::n = n in main() or something
	const int N = 101; // Max matrix size DONT TOUCH outside this namespace
	const int SIZE = N * N;
	struct mat{
		int v[SIZE];
		mat(int x=0){
			fill(v, v + SIZE, x);
		}
		int &at(int i, int j){
			return v[i * n + j];
		}
	};

	mat mult(mat &a, mat &b){
		mat c(0); // THIS NEEDS TO BE THE IDENTITY ELEMENT (i.e. 0 for normal mat mult, -INF for sum of max)
		f(i, n) f(j, n) f(k, n) {
			c.at(i, k) += (a.at(i, j) * 1ll * b.at(j, k))%mod;
			c.at(i, k) %= mod;
		}
		return c;
	}

	vector<int> apply(mat &a, vector<int> &b){
		vector<int> c(n, 0);
		f(i, n) f(j, n) {
			c[j] += (a.at(i, j) * 1ll * b[i])%mod;
			c[j] %= mod;
		}
		return c;
	}

	mat xp(mat &x, int k){
		mat ans; f(i, n) ans.at(i, i) = 1;
		while(k){
			if(k & 1) ans = mult(ans, x);
			x = mult(x, x);
			k >>= 1;
		}
		return ans;
	}
}
