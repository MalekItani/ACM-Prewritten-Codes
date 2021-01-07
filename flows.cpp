#include <bits/stdc++.h>
#define ll long long
#define ld long double
#define pi pair<ll, ll>
#define mp make_pair
#define pb push_back
#define vi vector<ll>
#define vvi vector<vi>
#define pd pair<double, double>
#define vp vector<pi>
#define f(i, n) for(int i = 0; i < n; i++)
#define fo(i, a, n) for(int i = a; i < n; i++)
#define mod (1000*1000*1000+7)
#define sz(x) int((x).size())
#define all(x) x.begin(),x.end()
#define DBG(v) cerr << #v << " = " << (v) << endl;
#define tc(i) int test_casesq; cin >> test_casesq; f(i, test_casesq)
ll INF = (1ll<<61);
using namespace std;

// MAX FLOW - DINIC
struct Edge {
    ll u, v;
    ll cap, flow;
    ll id;
    Edge() {}
    Edge(ll u, ll v, ll cap, ll id = 0) : u(u), v(v), cap(cap), flow(0), id(id) {}
};

int cnt = 1;
struct Dinic {
    // 0 based
    ll N;
    vector<Edge> E;
    vector<vector<ll>> g;
    vector<ll> d, pt;
    Dinic(ll N) : N(N), E(0), g(N), d(N), pt(N) {}
    void AddEdge(ll u, ll v, ll capUV, ll capVU = 0, ll id = 0) {
        if (u != v) {
            E.emplace_back(Edge(u, v, capUV, cnt));
            g[u].emplace_back(E.size() - 1);
            E.emplace_back(Edge(v, u, capVU, -cnt));
            g[v].emplace_back(E.size() - 1);
            cnt++;
        }
    }
    bool BFS(ll S, ll T) {
        queue<ll> q({ S });
        fill(d.begin(), d.end(), N + 1);
        d[S] = 0;
        while (!q.empty()) {
            ll u = q.front(); q.pop();
            if (u == T) break;
            for (ll k : g[u]) {
                Edge& e = E[k];
                if (e.flow < e.cap && d[e.v] > d[e.u] + 1) {
                    d[e.v] = d[e.u] + 1;
                    q.emplace(e.v);
                }
            }
        }
        return d[T] != N + 1;
    }
    ll DFS(ll u, ll T, ll flow = -1) {
        if (u == T || flow == 0) return flow;
        for (ll& i = pt[u]; i < g[u].size(); ++i) {
            Edge& e = E[g[u][i]];
            Edge& oe = E[g[u][i] ^ 1];
            if (d[e.v] == d[e.u] + 1) {
                ll amt = e.cap - e.flow;
                if (flow != -1 && amt > flow) amt = flow;
                if (ll pushed = DFS(e.v, T, amt)) {
                    e.flow += pushed;
                    oe.flow -= pushed;
                    return pushed;
                }
            }
        }
        return 0;
    }
    ll MaxFlow(ll S, ll T) {
        ll total = 0;
        while (BFS(S, T)) {
            fill(pt.begin(), pt.end(), 0);
            while (ll flow = DFS(S, T))
                total += flow;
        }
        return total;
    }
};

// GETTING ALL FLOW PATHS IN TERMS OF EDGES
vi vis;
ll dfs(ll v, Dinic &d, ll t, vector<Edge> &path){
  if(v == t){
    return INF;
  }
  // DBG(v)
  vis[v] = 1;
  vector<ll> &a = d.g[v];
  ll mflow = 0;
  for(auto k : a){
    Edge &e = d.E[k];
    if(vis[e.v]) continue;
    if(e.flow > 0){
      path.pb(e);
      mflow = min(e.flow, dfs(e.v, d, t, path));
      break;
    }
  }
  return mflow;
}

void update(ll v, Dinic &d, ll flow, vector<Edge> &path, ll i){
  if(i == sz(path)) return ;
  vector<ll> &a = d.g[v];
  vis[v] = 0;
  bool flag = 0;
  for(int j = 0; j < int(a.size()); j++){
    int k = a[j];
    Edge &e = d.E[k];
    if(e.id == path[i].id){
      assert(e.flow >= flow);
      e.flow -= flow;
      update(e.v, d, flow, path, i+1);
      if(e.flow == 0) {
        swap(a[j], a.back());
        a.pop_back();
      }
      flag = 1;
    }
    else if(e.flow == 0){
      swap(a[j], a.back());
      j--;
      a.pop_back();
    }
    if(flag) break;
  }
}

// RETURNS LIST OF PATHS WHERE EACH PATH IS GIVEN BY ITS EDGES.
// N is the number of nodes in the graph, S is source node, T is target node.
vector<vector<Edge> > get_decomposition(int n, Dinic &din, int s, int t){
  vector<Edge> p;
  vector<vector<Edge> > paths;
  ll res;
  vis.clear();
  vis.resize(n+1, 0);
  while((res = dfs(s, din, t, p)) ){
    assert(res > 0);
    paths.pb(p);
    update(s, din, res, p, 0);
    p.clear();
  }
  return paths;
}
