// SEGMENT TREE - ANTOINE

template<class V, class U> struct ST {
	int n;
	vector<V> st;
	vector<U> lazy;
	vector<int> leaves;

	ST( int n) : n(n) {
		leaves.resize(n);
		init(1, 0, n - 1);
		lazy.resize(st.size());
	}

	void init( int si,  int lo,  int hi) { // init st[] and leaves[]
		if (lo == hi) {
			if(si >= (int)st.size()) st.resize(si + 1);
			st[si] = V(); // optional if V(V(), V()) = V() ---- | ---- cin >> ??
			leaves[lo] = si;
		} else {
			 int mid = (lo + hi) >> 1;
			init(si << 1, lo, mid);
			init(si << 1 | 1, mid + 1, hi);
			st[si] = V(st[si << 1], st[si << 1 | 1]); // optional if T(T(), T()) = T()
		}
	}

	void updateLazy( int si,  int lo,  int hi) {
		lazy[si].apply(st[si], lo, hi);
		if (lo != hi) {
			lazy[si << 1] = U(lazy[si << 1], lazy[si]);
			lazy[si << 1 | 1] = U(lazy[si << 1 | 1], lazy[si]);
		}
		lazy[si] = U();
	}

	V query( int l,  int r) { return (l <= r && l < n && r >= 0) ? query(l, r, 1, 0, n - 1) : V(); }
	V query( int l,  int r,  int si,  int lo,  int hi) {
		updateLazy(si, lo, hi);
		if (l <= lo && hi <= r) return st[si];

		 int mid = (lo + hi) >> 1;
		if (r <= mid) return query(l, r, si << 1, lo, mid);
		if (mid < l) return query(l, r, si << 1 | 1, mid + 1, hi);
		return V(query(l, r, si << 1, lo, mid), query(l, r, si << 1 | 1, mid + 1, hi));
	}


	void update( int l,  int r,  U u) { if (l <= r) update(l, r, u, 1, 0, n - 1); }
	void update( int l,  int r,  U &u,  int si,  int lo,  int hi) {
		if (l <= lo && hi <= r) {
			lazy[si] = U(lazy[si], u);
			updateLazy(si, lo, hi);
		} else {
			updateLazy(si, lo, hi);
			if (hi < l || r < lo) return;
			 int mid = (lo + hi) >> 1;
			update(l, r, u, si << 1, lo, mid);
			update(l, r, u, si << 1 | 1, mid + 1, hi);
			st[si] = V(st[si << 1], st[si << 1 | 1]);
		}
	}

	// ------------------------------------------| ________ |---------------------------------------------------------
	// ------------------------------------------| NON-LAZY |---------------------------------------------------------

	void update( int i, U u) { // single index update, [!!!] NON-LAZY SEGMENT TREE ONLY !!!
		int si = leaves[i];
		u.apply(st[si]);
		for (si >>= 1; si; si >>= 1)
			st[si] = V(st[si << 1], st[si << 1 | 1]);
	}

	// ------------------------------------------| ________ |---------------------------------------------------------
	// ------------------------------------------| OPTIONAL |---------------------------------------------------------

	// remove and return all non-zero values in a range, assuming all values >= 0
	vector<int> remAll( int l,  int r) { vector<int> res; remAll(l, r, 1, 0, n - 1, res); return res; }
	void remAll( int l,  int r,  int si,  int lo,  int hi, vector<int> &res) {
		updateLazy(si, lo, hi);
		if (!st[si] || hi < l || r < lo) return;
		if (lo == hi) {
			res.push_back(lo);
			st[si] = 0;
		}
		else {
			 int mid = (lo + hi) >> 1;
			remAll(l, r, si << 1, lo, mid, res);
			remAll(l, r, si << 1 | 1, mid + 1, hi, res);
			st[si] = V(st[si << 1], st[si << 1 | 1]);
		}
	}

	void updateAllLazy() { updateAllLazy(1, 0, n - 1); }
	void updateAllLazy( int si,  int lo,  int hi) {
		updateLazy(si, lo, hi);
		if (lo != hi) {
			 int mid = (lo + hi) >> 1;
			updateAllLazy(si << 1, lo, mid);
			updateAllLazy(si << 1 | 1, mid + 1, hi);
		}
	}

	vector<V> retrieve() {
		updateAllLazy();
		vector<V> res(n);
		for (int i = 0; i < n; ++i)
			res[i] = st[leaves[i]];
		return res;
	}
};

struct StVal { // define
	ll v = 0; // optional
	StVal() {}
	StVal( ll v) : v(v) {} // optional
	StVal( StVal &v1,  StVal &v2) { v = v1 + v2; }
	operator ll()  { return v; } // optional
};

struct StUpdate { // define
	ll v = 0; // optional
	StUpdate() {}
	StUpdate( ll v) : v(v) {} // optional
	StUpdate( StUpdate &u1,  StUpdate &u2) { v = u1 + u2; }
	operator ll()  { return v; } // optional
	void apply(StVal &v,  int lo,  int hi) { // omit last 2 arguments if non-lazy segment tree
		v.v += this->v * (hi - lo + 1);
	}
};

typedef ST<StVal, StUpdate> MyST;
