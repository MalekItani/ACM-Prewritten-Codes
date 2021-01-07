// SPARSE TABLE / STATIC RANGE QUERIES - Homemade
// Usage: Provide the custom func you want in the constructor, along with the array.
// Works with any type T so long as func takes two Ts and returns another T.
// To query between [L, R], call q(L, R)
template<class C> struct rq{
	vector<vector<C>> v;
	function<C(C, C)> func;
	rq(){}
	rq(vector<C> &arr, function<C(C, C)> _f){
		func = _f;
		int n = sz(arr), k = 32 - __builtin_clz(n);
		v.resize(n, vector<C>(k+1));
		for (int i = 0; i < n; i++) v[i][0] = arr[i];
		for(int j = 1; j < k; j++ ) for(int i = 0; i <= n-(1<<j); i++) v[i][j] = func(v[i][j-1],v[i+(1<<(j-1))][j-1]);
	}
	C q(int l, int r){
		int j = 31 - __builtin_clz(r - l + 1);
		return func(v[l][j], v[r - (1 << j) + 1][j]);
	}
};

// CONVEX HULL TRICK
const ll INF = 1LL<<61; bool cmpA;
struct Line { // ax+b. xl = the x-coord of intersection with previous line in set
    ll a, b;
    mutable long double xl;
    bool operator < (const Line &l) const { if (cmpA) return a < l.a; else return xl
    < l.xl; }
};
struct DynamicHull : multiset<Line> {
    iterator prev (iterator y) { return --y; } // skip if using c++0x/c++11
    iterator next (iterator y) { return ++y; } // skip if using c++0x/c++11
    bool bad (iterator y) {
        iterator z = next(y), x;
        if (y == begin()) {
        if (z == end()) return 0;
        return y->a == z->a && y->b <= z->b;
        }
        x = prev(y);
        if (z == end()) return y->a == x->a && y->b <= x->b;
        return 1.0L * (x->b - y->b) * (z->a - y->a) >= 1.0L * (y->b - z->b) * (y->a -
        x->a);
    }
    void add (ll a, ll b) { // add a line ax+b to the convex hull
        cmpA = 1; iterator y = insert((Line){a,b,-INF});
        if (bad(y)) { erase(y); return; }
        while (next(y) != end() && bad(next(y))) erase(next(y));
        while (y != begin() && bad(prev(y))) erase(prev(y));
        if (next(y) != end()) next(y)->xl = 1.0L * (y->b - next(y)->b) / (next(y)->a -
        y->a);
        if (y != begin()) y->xl = 1.0L * (y->b - prev(y)->b) / (prev(y)->a - y->a);
    }
    ll eval (ll x) { // find max ax+b, over all lines added, for given x
        if (empty()) return -INF;
        cmpA = 0; iterator it = prev(lower_bound((Line){0,0,1.0L*x}));
        return it->a * x + it->b;
    }
};


// ORDER STATISTIC TREE

// find_by_order(k): Returns an iterator to the k-th element
// erase(k): Erases elements with value k
// order_of_key(k): Returns the "index" of some value k. If it is not in the set, then returns 0.
#include <ext/pb_ds/assoc_container.hpp> // Common file
#include <ext/pb_ds/tree_policy.hpp> // Including tree_order_statistics_node_update

using namespace __gnu_pbds;
typedef tree<int, null_type, less<int>, rb_tree_tag, tree_order_statistics_node_update> ordered_set;
// typedef tree<key, value, comparator, base structure, update to use when nodes added/removed> ordered_set;

// FASTER HASHTABLES
#include <ext/pb_ds/assoc_container.hpp>
using namespace __gnu_pbds;
gp_hash_table<int, int> table;


// AHO-CORASICK
// usage: first add all words (addWord(s)), then calcLink(), then make searches
// using searchAll or searchFirst
struct mkmp {
	struct node {
		map<char,node*> son, cache;
		node *lnk, *wo; // pi function; link to max prefsufa which is a pattern
		int el; // number of pattern (-1 if not a pattern)
		node() : el(-1) {}
	};
	vi len;
	node root;
	node* mv (node *w, char l) {
		node*& r = w->cache[l];
		return r ? r : r = (w->son.count(l) ? w->son[l] : w == &root ? w : mv(w->lnk,l)); // watch out if there are stack limit – deep recursion!
	}
	// adds a pattern and assigns a number to it; returns this number
	int addWord (string &s) {
		int l = sz(s);
		node *p = &root;
		for (int i  = 0; i < sz(s); ++i) {
			auto e = p->son.find(s[i]);
			p = (e == p->son.end()) ? p->son[s[i]] = new node : e->second;
		}
		if (p->el == -1) { p->el = sz(len); len.pb(l); }
		return p->el;
	}

	void calcLink () {
	 vector<node*> l;
	 node *w;
	 root.lnk = root.wo = 0;
	 for(auto it = root.son.begin(); it != root.son.end(); it++) {
		 l.pb(it->second); it->second->lnk = &root;
	 }
	 f(x, sz(l)) {
		 l[x]->wo = (l[x]->lnk->el != -1) ? l[x]->lnk : l[x]->lnk->wo;
		 for(auto it = l[x]->son.begin(); it != l[x]->son.end(); it++) {
		 	l.pb(it->second); w = l[x]->lnk; w = mv(w,it->first); it->second->lnk = w;
		 }
	 }
}
 // finds all occurences of every pattern, result: {<position, number of pattern>}
 vector<pi> searchAll(string &s) {
	 vector<pi> ans;
	 node *p = &root, *r;
	 for (int x = 0; s[x]; ++x) {
		 p = mv(p, s[x]);
		 for (r = p; r; r = r->wo)
		 	if (r->el != -1) ans.pb(mp(x - len[r->el] + 1, r->el));
	 }
	 return ans;
 }
// finds at most one (the first) occurence of every pattern,
// result: {<position, number of pattern>} (watch out: it breaks ->el pointers)
 vector<pi> searchFirst(string &s) {
	 vector<pi> ans;
	 node *p = &root, *r, *t;
	 for (int x = 0; s[x]; ++x) {
	 	r = p = mv(p, s[x]);
		while (r) {
			if (r->el != -1) ans.pb(mp(x - len[r->el] + 1, r->el));
			r->el = -1; t = r; r = r->wo; t->wo = 0;
			}
	 }
	 return ans;
 }
 // finds at every position the longest pattern occuring there;
 // result: {<position, number of pattern>}
 vector<pi> searchLongest(string &s) {
	 vector<pi> ans;
	 node *p = &root, *r;
	 for (int x = 0; s[x]; ++x) {
	 	r = p = mv(p, s[x]);
	 	if (r->el == -1) {
	 		r = r->wo;
	 		if (!r || r->el == -1) continue;
 		}
		ans.pb(mp(x - len[r->el] + 1, r->el));
	}
	return ans;
	}
};

// SUFFIX AUTOMATON

struct SuffixAutomaton {
  // ROOT IS 0
  vector<map<char,int>> edges; // edges[i]  : the labeled edges from node i
  vector<int> link;            // link[i]   : the parent of i
  vector<int> length;          // length[i] : the length of the longest string in the ith class
  int last;                    // the index of the equivalence class of the whole string
  vector<int> terminals;

  SuffixAutomaton(string s) {
    // add the initial node
    edges.push_back(map<char,int>());
    link.push_back(-1);
    length.push_back(0);
    last = 0;

    for(int i=0;i<s.size();i++) {
      // construct r
      edges.push_back(map<char,int>());
      length.push_back(i+1);
      link.push_back(0);
      int r = edges.size() - 1;

      // add edges to r and find p with link to q
      int p = last;
      while(p >= 0 && edges[p].find(s[i]) == edges[p].end()) {
        edges[p][s[i]] = r;
        p = link[p];
      }
      if(p != -1) {
        int q = edges[p][s[i]];
        if(length[p] + 1 == length[q]) {
          // we do not have to split q, just set the correct suffix link
          link[r] = q;
        } else {
          // we have to split, add q'
          edges.push_back(edges[q]); // copy edges of q
          length.push_back(length[p] + 1);
          link.push_back(link[q]); // copy parent of q
          int qq = edges.size()-1;
          // add qq as the new parent of q and r
          link[q] = qq;
          link[r] = qq;
          // move short classes pointing to q to point to q'
          while(p >= 0 && edges[p][s[i]] == q) {
            edges[p][s[i]] = qq;
            p = link[p];
          }
        }
      }
      last = r;
    }

    int p = last;
    while(p > 0) {
      terminals.push_back(p);
      p = link[p];
    }
  }
  bool has_as_next(int i, char c){
    return edges[i].find(c) != edges[i].end();
  }
  int next(int i, char c){
    return edges[i][c];
  }
};


// SUFFIX ARRAY - https://codeforces.com/blog/entry/4025

namespace SA{
	const int MAXN = 1 << 21;
	string S;
	int N, gap;
	int sa[MAXN], pos[MAXN], tmp[MAXN], lcp[MAXN];
	// sa[i] = If you sort all suffixes, the suffix starting at sa[i] comes in the ith position
	bool sufCmp(int i, int j){
		if (pos[i] != pos[j])
			return pos[i] < pos[j];
		i += gap;
		j += gap;
		return (i < N && j < N) ? pos[i] < pos[j] : i > j;
	}

	void buildSA(){
		N = sz(S);
    for (int i = 0; i < (int)(N); ++i){
      sa[i] = i, pos[i] = S[i];
    }
		for (gap = 1;; gap *= 2){
			sort(sa, sa + N, sufCmp);
      for (int i = 0; i < N-1; ++i){
        tmp[i + 1] = tmp[i] + sufCmp(sa[i], sa[i + 1]);
      }
      for (int i = 0; i < (int)(N); ++i){
        pos[sa[i]] = tmp[i];
      }

			if (tmp[N - 1] == N - 1) break;
		}
	}

	void buildLCP(){
		for (int i = 0, k = 0; i < N; ++i) if (pos[i] != N - 1){
			for (int j = sa[pos[i] + 1]; S[i + k] == S[j + k];) ++k;
			lcp[pos[i]] = k;
			if (k)--k;
		}
	}
}
