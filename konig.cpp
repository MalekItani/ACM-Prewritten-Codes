// KONIG THEOREM - FIND VERTEX COVER GIVEN MAXIMUM BIPARTITE MATCHING
// FOR MAXIMUM INDEPENDENT SET TAKE COMPELEMENT OF VERTEX COVER

struct Edge{
  int u, v;
  Edge(){}
  Edge(int u, int v):u(u), v(v){}
};

// Find vertex cover given graph and matching edges (does not modify G)
// G: Bipartite graph with n by m nodes, numbered 1 ... n on left and n+1 ... n+m on right.
// E: Edges that form the matching.
// n: Size of first partition of bipartite graph
// m: Size of second partition of bipartite graph
// Returns pair of vector of ints: vertices included on left side numbered from 1 ... n and
// vertices included on right side numbered from 1 ... m
pair<vector<int>, vector<int> > konig(vector<set<int> > &G, vector<Edge> &E, int n, int m){

  vector<set<int> > Gp = G;
  vector<int> l, r;
  vector<int> inc(n+m+1, 1);

  for(int i = 0; i < (int) E.size(); i++){
    int u = E[i].u;
    int v = E[i].v;
    if(u > v) swap(u, v);
    Gp[u].erase(v);
    inc[u] = 0;
  }

  for(int u = 1; u <= n; u++){
    for(auto v : Gp[u]){
      Gp[v].erase(u);
    }
  }

  for(int i = 0; i < n + m + 1; i++) vis[i] = 0;

  for(int i = 1; i <= n; i++){
    if(inc[i]){
      dfs(Gp, i);
    }
  }

  for(int i = 1; i <= n; i++){
    if(!vis[i]) l.pb(i);
  }

  for(int i = n+1; i <= n+m; i++) if(vis[i]) r.pb(i - n);

  return make_pair(l, r);
}
