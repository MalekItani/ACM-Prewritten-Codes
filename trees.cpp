// LCA (UNWEIGHTED) - Homemade
// O(1) per query

struct rq{
	vector<vector<int>> v;
	rq(){}
	rq(vector<int> &arr){
		int n = sz(arr), k = 32 - __builtin_clz(n);
		v.resize(n, vector<int>(k+1));
		for (int i = 0; i < n; i++) v[i][0] = arr[i];
		for(int j = 1; j < k; j++ ) for(int i = 0; i <= n-(1<<j); i++) v[i][j] = min(v[i][j-1],v[i+(1<<(j-1))][j-1]);
	}
	int q(int l, int r){
		int j = 31 - __builtin_clz(r - l + 1);
		return min(v[l][j], v[r - (1 << j) + 1][j]);
	}
};

struct LCA{
  vector<int> height, id, trav, foc, inv;
  rq mq;
  int i = 0;
  LCA(){}
  LCA(vvi &adj, int root = 1){
    int n = sz(adj);
    height.resize(n, 0);
    inv = id = foc = height;
    dfs(root, -1, 0, adj);
    mq = rq(trav);
  }
  void dfs(int v, int pV, int h, vvi &g){
    height[v] = h, id[v] = i, inv[i] = v, foc[v] = sz(trav);
    trav.pb(i++);
    for(auto &k : g[v]){
      if(k == pV) continue;
      dfs(k, v, h+1, g);
      trav.pb(id[v]);
    }
  }
  int lca(int u, int v){
    int x = foc[u], y = foc[v];
    if(x > y) swap(x, y);
    return inv[mq.q(x, y)];
  }
  int dist(int u, int v){
    int x = lca(u, v);
    return abs(height[x] - height[u]) + abs(height[x] - height[v]);
  }
};

// LCA (WEIGHTED) - Homemade
// O(1) per query

struct rq{
	vector<vector<int>> v;
	rq(){}
	rq(vector<int> &arr){
		int n = sz(arr), k = 32 - __builtin_clz(n);
		v.resize(n, vector<int>(k+1));
		for (int i = 0; i < n; i++) v[i][0] = arr[i];
		for(int j = 1; j < k; j++ ) for(int i = 0; i <= n-(1<<j); i++) v[i][j] = min(v[i][j-1],v[i+(1<<(j-1))][j-1]);
	}
	int q(int l, int r){
		int j = 31 - __builtin_clz(r - l + 1);
		return min(v[l][j], v[r - (1 << j) + 1][j]);
	}
};

#define vvp vector<vp>

struct LCA{
  vector<int> id, trav, foc, inv;
  vi height;
  rq mq;
  int i = 0;
  LCA(){}
  LCA(vvp &adj, int root = 1){
    int n = sz(adj);
    height.resize(n);
    id.resize(n, 0);
    inv = foc = id;
    dfs(root, -1, 0, adj);
    mq = rq(trav);
  }
  void dfs(int v, int pV, int h, vvp &g){
    height[v] = h, id[v] = i, inv[i] = v, foc[v] = sz(trav);
    trav.pb(i++);
    for(auto &k : g[v]){
      if(k.first == pV) continue;
      dfs(k.first, v, h+k.second, g);
      trav.pb(id[v]);
    }
  }
  int lca(int u, int v){
    int x = foc[u], y = foc[v];
    if(x > y) swap(x, y);
    return inv[mq.q(x, y)];
  }
  ll dist(int u, int v){
    int x = lca(u, v);
    return abs(height[x] - height[u]) + abs(height[x] - height[v]);
  }
};


// Binary Lifting (Unweighted, for weighted, use vvp and modify dfs as usual) - Homemade
// O(logn) per query
// up(x, k): Go k steps up from node x.

struct BL{
	vector<vector<int>> lift; int root;
//	vector<vector<ll>> eval; // OPTIONAL
	BL(){}
	BL(vvi &g, int r = 1) : root(r){
		lift.resize(sz(g) + 1);
//		eval.resize(sz(g) + 1); // OPTIONAL
		dfs(g, r, -1);
	}
	void dfs(vvi &g, int v, int pv){
		if(pv != -1){
			lift[v].pb(pv);
//			eval[v].pb(k); // OPTIONAL (k is value from v to pv like edge weight)
			for(int i = 1; i < 21; i++){
				if(sz(lift[lift[v][i-1]]) < i) break;
				lift[v].pb(lift[lift[v][i-1]][i-1]);
//				 eval[v].pb(min(eval[v][i-1], eval[lift[v][i-1]][i-1])); // (OPTIONAL);
			}
		}
		for(auto &k : g[v]){
			if(k == pv) continue;
			dfs(g, k, v);
		}
	}
	int up(int x, int k){
		if(k == 0) return x;
		int y = 31 - __builtin_clz(k);
		assert(sz(lift[x]) > y); // Error if there aren't that many nodes from x to root. Change to "return root;" if needed.
		return up(lift[x][y], k - (1 << y));
	}

//	Useful to query eval from x to k steps above, use with LCA to query along paths (here it is min)
//	int q(int x, int k){
//		if(k == 0) return 1e9;
//		int y = 31 - __builtin_clz(k);
//		assert(sz(lift[x]) > y); // Error if there aren't that many nodes from x to root. Change to "return root;" if needed.
//		return min(eval[x][y], up(lift[x][y], k - (1 << y)));
//	}
};

// SegTree on Tree (Ex: dfs order)

struct treeST{
	ST<StVal, StUpdate> stree;
	vector<int> st, en, order;
	treeST(){}
	treeST(vvi &g, int root = 1){
		st.resize(sz(g), 0);
		en = st;
		int t = 0;
		dfs(g, root, -1, t);
		stree = ST<StVal, StUpdate>(sz(g));
	}
	void dfs(vvi &g, int v, int pv, int &t){
		order.pb(v);
		st[v] = t;
		for(auto &k : g[v]){
			if(k == pv) continue;
			dfs(g, k, v, ++t);
		}
		en[v] = t;
	}
	StVal query(int v){ return stree.query(st[v], en[v]);} // Subtree query
	void update(int v, StUpdate u){ return stree.update(st[v], en[v], u);} // Subtree update
};




// Number of rooted forests of x nodes with y different trees is:
// f(x, y) = y * x ^ (x-y-1)



