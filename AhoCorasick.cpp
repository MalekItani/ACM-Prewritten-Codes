ll nnodes = 0;

struct node{
    node *bck, *wo; // Direct fail link, fail link to longest terminal node.
    char c; // Character for this node.
    int id;
    // Put additional member variables here

    set<int> s; // Set of word ids that end at this node.
    map<char, node*> next; // Set of next characters for current node.
    node(char c):c(c), id(nnodes++){}
    node* add(char c, node *root, bool i, int k){
        if(!next.count(c)){
            next[c] = new node(c);
            next[c]->bck = root;
            next[c]->wo = root;
        }
        // Update **NEXT** members here:
        // For example, if nodes have member called "d", and you
        // want to set d to 0, write next[c]->d = 0;

        if(i) next[c]->s.insert(k);

        return next[c];
    }
    bool has(char c){
        return next.count(c);
    }
    node* to(char c){
        if(has(c)) return next[c];
        if(bck->id == 0 && !bck->has(c)) return bck;
        return bck->to(c);
    }
};

// Aho-Corasick Trie (Trie + fail links) - Homemade
// Usage: add() strings of dictionary and then call calcLinks();
// Move to the next character using next()
struct Trie{
    node *root, *curr;
    int nw = 0; // Number of words
    vi visw, visn, len; // Optional
    Trie(){
        nnodes = 0;
        root = new node('.');
        root->bck = root;
        root->wo = root;
        curr = root;
    }
    void add(string &s){
        // Adds a string to the trie.
        curr = root;
        len.pb(sz(s));
        f(i, sz(s)) curr = curr->add(s[i], root, i==sz(s) - 1, nw);
        nw++;
    }
    void reset(){
        // Places current pointer back to root
        curr = root;
    }
    node * next(char c){
        // Moves to the next node with the specified char (possibly using fail links).
        // If there is no such node, goes to root.
        return curr = curr->to(c);
    }

    // Optional if using basic trie.
    void calcLinks(){
        // Calculates fail links.
        vector<node*> l;
        for(auto &n : root->next){
            l.pb(n.second);
            n.second->bck = root;
        }
        f(x, sz(l)) {
            l[x]->wo = sz(l[x]->bck->s) ? l[x]->bck : l[x]->bck->wo;
            for(auto &k : l[x]->next){
                l.pb(k.second);
                curr = l[x]->bck;
                curr = next(k.first);
                k.second->bck = curr;
            }
        }
    }

    // OPTIONAL (Need calcLinks())

    // findExisting: Given a string, returns a vector of bools where v[i] = 1 if the i-th string
    // can be found in s as a substring. Useful for offline substring queries.
    vi findExisting(string &s){
        visw.resize(nw, 0); visn.resize(nnodes+1, 0);
        visn[0] = 1;
        reset();
        f(i, sz(s)){
            node *r = next(s[i]);
            while(!visn[r->id]){
                for(auto &k : r->s) visw[k] = 1;
                visn[r->id] = 1;
                r = r->wo;
            }
            visn[curr->id] = 1;
        }
        return visw;
    }

    // findAll: Given a string, returns a vector of pairs (pos, type)
    // where pos is the position of the string and type is the index
    // (in insertion order)
    vp findAll(string &s){
        reset();
        vp res; node *r;
        f(i, sz(s)){
            r = next(s[i]);
            for(auto p = r; p->id ; p=p->wo)
            	for(auto &k : p->s) res.pb({i-len[k]+1, k});
        }
        return res;
    }
    // findAll: Given a string, returns a pair (pos, type) first occurence
    // of a word from the dictionary. It is found at position pos and it is
    // the type-th word (in insertion order) // UNTESTED
    pi findFirst(string &s){
        reset();
        node *r;
        f(i, sz(s)){
            r = next(s[i]);
            if(sz(r->s)) return {i-len[*r->s.begin()]+1, *r->s.begin()};
        }
        return {-1,-1};
    }

};
