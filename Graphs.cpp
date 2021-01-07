// DFS
void dfs(ll v, ll pV){
  vi &adj = g[v];
  for(auto k : adj){
    if(k == pV) continue;
    dfs(k, v);
  }
}

// BFS
void bfs(int s, int n){
  vi vis(n+1, 0);
  vi dis(n+1, 0);
  queue<int> q;
  q.push(s);
  vis[s] = 1;
  while(!q.empty()){
    s = q.front(); q.pop();
    vi &adj = g[s];
    for(auto k : adj){
      if(vis[k]) continue;
      dis[k] = dis[s] + 1;
      vis[k] = 1;
      q.push(k);
    }
  }
}


// DIJKSTRA
vi dijk(vvp &g, int s){
	int n = sz(g);
	vi dis(n, INF), vis(n, 0);
	priority_queue<pair<ll, int> > pq;
	pq.push({0, s});
	dis[s] = 0;
	while(!pq.empty()){
		pair<ll, ll> p = pq.top(); pq.pop();
		if(vis[p.second]) continue;
		vis[p.second] = 1;
		for(auto &k : g[p.second]){
			if(vis[k.second]) continue;
			dis[k.second] = min(dis[k.second], dis[p.second] + k.first);
			pq.push({-dis[k.second], k.second});
		}
	}
	return dis;
}

// PRIM
// input: (n,adj) – connected graph on vertices 1..n
// output: cost of MST.
set<int> s; // This will contain the ids of the edges in the MST.
ll prim (vvvi &g, int st) {
	int n = sz(g);
	vis.resize(n,0);
	ll cost = 0;
	priority_queue<pair<ll,vector<ll>> > q;
	q.push(mp(0,vector<ll>({-1,st,-1}))); // minus cost, from, to, edgeid
	while (!q.empty()) {
		ll v = q.top().second[1], d = -q.top().first, from = q.top().second[0];
		ll id = q.top().second[2];
		q.pop();
		if (vis[v]) continue;
		vis[v] = 1;
		cost += d;

		if (from != -1) s.insert(id);

		for(auto &k : g[v]) if (!vis[k[1]] && ex != k[2]) q.push(mp(-k[0],vector<ll>({v, k[1], k[2]})));
	}
	return cost;
}


// DFS + TOPOLOGICAL SORT
vector<int> *g;
stack<int> order;
int vis[3005];
bool dfs(int node){
    vis[node] = 1;
    bool b = 0;
    for(int i = 0; i < (int) g[node].size(); i++){
        if(vis[g[node][i]] == 0){
            b = b || dfs(g[node][i]);
        }else if(vis[node] == 1){
            b = 1; break;
        }
    }
    vis[node] = 2;
    order.push(node);
    return b;
}


// LONGEST PATH IN DAG
vi *g;
int vis[1005];
int dis[1005];
int ans = 0;

void dfs(int Node){
    vis[Node] = 1;
    for(int i = 0; i < g[Node].size(); ++ i){
        if(!vis[g[Node][i]]) dfs(g[Node][i]);
        dis[Node] = max(dis[Node], dis[ g[Node][i] ] + 1);
    }
    dis[Node] = max(dis[Node], 1);
    ans = max(ans, dis[Node]);
}


// SHORTEST CYCLE (**DIRECTED** graph)
// Returns vector of ints representing nodes in shortest cycle, if no cycle then the vector is empty
// n: Number of nodes
// g: Adjacency list
vi shortestCycle(int n, vi *g){
	vi sc;
	fo(i, 1, n+1){
		vi p(n+1, 0);
		queue<int> q;
		q.push(i);
		p[i] = -1;
		vi c;
		while(!q.empty()){
			int t = q.front();
			q.pop();
			vi &adj = g[t];
			for(auto k : adj){
				if(p[k] == 0){
					p[k] = t;
					q.push(k);
				}else if(k == i){
					c.pb(k);
					while(t != i){
						c.pb(t);
						t = p[t];
					}
					if(sz(sc) == 0 || sz(c) < sz(sc)){
						sc = c;
					}
					c.clear();
					while(!q.empty()) q.pop();
					break;
				}
			}
		}
	}
	return sc;
}

// SHORTEST CYCLE
// Returns length of the shortest cycle, if no cycle then returns 1e9 (1-BASED) (**UNDIRECTED**)
// n: Number of nodes
// g: Adjacency list
int shortestCycle(int n, vi *g){
	int sc = 1e9;
	fo(i, 1, n+1){
		vi dis(n+1, 1e9);
		vi p(n+1, 0);
		queue<int> q;
		q.push(i);
		dis[i] = 0;
		while(!q.empty()){
			int s;
			s = q.front(); q.pop();
			vi &adj = g[s];
			for(auto k : adj){
				if(dis[k] < 1e9) {
					if(k != p[s]){
						sc = min((ll) sc, dis[k] + dis[s] + 1);
					}
					continue;
				}
				dis[k] = dis[s] + 1;
				p[k] = s;
				q.push(k);
			}
		}
	}
	return sc;
}
