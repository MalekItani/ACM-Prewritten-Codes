// ORDER STATISTIC TREE

// find_by_order(k): Returns an iterator to the k-th element
// erase(k): Erases elements with value k
// order_of_key(k): Returns the "index" of some value k. If it is not in the set, then returns 0.
#include <ext/pb_ds/assoc_container.hpp> // Common file
#include <ext/pb_ds/tree_policy.hpp> // Including tree_order_statistics_node_update

using namespace __gnu_pbds;
typedef tree<pair<int, int>, null_type, less<pair<int, int>>, rb_tree_tag, tree_order_statistics_node_update> ordered_set;
// typedef tree<key, value, comparator, base structure, update to use when nodes added/removed> ordered_set;


// SEGMENT TREE - ANTOINE

struct MT {
	int n;
	vector<ordered_set> st;
	vector<int> leaves;

	MT(const int n, vector<int> &in) : n(n) {
		leaves.resize(n);
		init(1, 0, n - 1, in);
	}

	void init(const int si, const int lo, const int hi, vector<int> &val) { // init st[] and leaves[]
		if (lo == hi) {
			if(si >= (int)st.size()) st.resize(si + 1);
			leaves[lo] = si;
			st[si].insert(mp(val[lo], lo));
		} else {
			const int mid = (lo + hi) >> 1;
			init(si << 1, lo, mid, val);
			init(si << 1 | 1, mid + 1, hi, val);

			vector<pair<int, int>> v;
			ordered_set &a = st[si << 1];
			ordered_set &b = st[si << 1 | 1];
			merge(a.begin(), a.end(), b.begin(), b.end(), back_inserter(v));
			st[si] = ordered_set(v.begin(), v.end());
		}
	}

	int query(const int &l, const int &r, int &p, int &u, int t) { return (l <= r && l < n && r >= 0) ? query(l, r, 1, 0, n - 1, p, u, t) : 0; }
	int query(const int &l, const int &r, const int &si, const int &lo, const int &hi, int &p, int &u, int t) {
		if (l <= lo && hi <= r){
			if(t == 0){
				pair<int, int> k1 = mp(p, 1e8);
				pair<int, int> k2 = mp(u, 1e8);

				int a = sz(st[si]) - st[si].order_of_key(k1);
				int b = sz(st[si]) - st[si].order_of_key(k2);

				return b - a;
			}else{
				pair<int, int> k1 = mp(p, -1e8);
				pair<int, int> k2 = mp(u, -1e8);

				int a = st[si].order_of_key(k1);
				int b = st[si].order_of_key(k2);

				return b - a;
			}
		}
		const int mid = (lo + hi) >> 1;
		if (r <= mid) return query(l, r, si << 1, lo, mid, p, u, t);
		if (mid < l) return query(l, r, si << 1 | 1, mid + 1, hi,p, u, t);
		return query(l, r, si << 1, lo, mid, p, u, t) + query(l, r, si << 1 | 1, mid + 1, hi, p, u, t);
	}

	void update(const int i, int u) { // single index update, *replaces ith element with u*
		int si = leaves[i];
		pair<int, int> old = *st[si].begin();
		pair<int, int> nw{u, i};
		st[si].erase(old);
		st[si].insert(nw);
		for (si >>= 1; si; si >>= 1){
			st[si].erase(old);
			st[si].insert(nw);
		}
	}
};
