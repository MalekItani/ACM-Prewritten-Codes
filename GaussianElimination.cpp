#define real long double

const long double EPS = 1e-12;
inline bool IsZero(real x){ return x>=-EPS && x<=EPS; }

int GaussDouble(vector<vector<real>> &a, vector<real> &b, vector<real> &x) {
	int m = a.size(), n = a[0].size(), k, r;
	vi q; // which column in taken care of on the k-th step
	for (k = 0; k < min(m, n); k++) {
		int i, j;
		for (j = k; j < n; j++)
			for (i = k; i < m; i++)
				if (!IsZero(a[i][j])) goto found;
		break;
		found:
		if (j != k) f(t, m) swap(a[t][j], a[t][k]);
		q.pb(j);
		if (i != k) {
			swap(a[i], a[k]);
			swap(b[i], b[k]);
		}
		fo(j,k+1,m) if (!IsZero(a[j][k])) {
			real l = (a[j][k] / a[k][k]);
			fo(i,k,n) a[j][i] = a[j][i] - (l * a[k][i]);
			b[j] = b[j] - (l * b[k]);
		}
	}
	r = k; // r == rank of matrix a
	x.clear(); x.resize(n, 0);
	fo(k,r,m) if (!IsZero(b[k])) return -1; // no solution
	for(int k = r-1; k >= 0; k--) {
		real s = b[k];
		fo(j,k+1,r) s = s - (a[k][j] * x[j]);
		x[k] = s / a[k][k];
	}
	for(int k = r-1; k>= 0; k--) swap(x[k], x[q[k]]); // determinant = product of a[i][i] over i
	return n - r; // the dimension of the space of solutions
}
