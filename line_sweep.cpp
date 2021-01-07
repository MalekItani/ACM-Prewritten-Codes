// GIVEN A SET OF SEGMENTS, COMPUTE ALL INTERSECTION POINTS
// Runs in O(nlogn) if checking for points, O(# intersection points + nlogn) when finding them


#define ld long double

ld EPS = 1e-15;

struct pt{
  ld x, y;
  pt() {}
  pt(ld x, ld y) : x(x), y(y) {}
  pt(const pt &p) : x(p.x), y(p.y)    {}
  pt operator + (const pt &p)  const { return pt(x+p.x, y+p.y); }
  pt operator - (const pt &p)  const { return pt(x-p.x, y-p.y); }
  pt operator * (ld c)     const { return pt(x*c,   y*c  ); }
  pt operator / (ld c)     const { return pt(x/c,   y/c  ); }
  bool operator == (const pt &p)     const { return fabs(x - p.x) < EPS && fabs(y - p.y) < EPS; }
};

double dot(pt p, pt q)     { return p.x*q.x+p.y*q.y; }
double dist2(pt p, pt q)   { return dot(p-q,p-q); }
double cross(pt p, pt q)   { return p.x*q.y-p.y*q.x; }
ostream &operator<<(ostream &os, const pt &p) {
  return os << "(" << p.x << "," << p.y << ")";
}

struct seg {
    pt p, q;
    int id;
    seg(pt& x, pt& y, int id) : p(x), q(y), id(id){}
    ld get_y (ld x) const {
        if (abs (p.x - q.x) < EPS)  return p.y;
        return p.y + (q.y - p.y) * (x - p.x) / (q.x - p.x);
    }
    bool operator<(const seg & b) const {
        ld x = max (min (p.x, q.x), min (b.p.x, b.q.x));
        return get_y(x) < b.get_y(x) - EPS;
    }
};


inline bool intersect1d (ld l1, ld r1, ld l2, ld r2) {
    if (l1 > r1)  swap (l1, r1);
    if (l2 > r2)  swap (l2, r2);
    return max (l1, l2) <= min (r1, r2) + EPS;
}

inline int vec (const pt & a, const pt & b, const pt & c) {
    ld s = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
    return abs(s)<EPS ? 0 : s>0 ? +1 : -1;
}


bool intersect (const seg & a, const seg & b) {
   return intersect1d (a.p.x, a.q.x, b.p.x, b.q.x)
          && intersect1d (a.p.y, a.q.y, b.p.y, b.q.y)
          && vec (a.p, a.q, b.p) * vec (a.p, a.q, b.q) <= 0
          && vec (b.p, b.q, a.p) * vec (b.p, b.q, a.q) <= 0;
}

pt ComputeLineIntersection(pt a, pt b, pt c, pt d) {
  b=b-a; d=c-d; c=c-a;
  return a + b*cross(c, d)/cross(b, d);
}

pt find_intp(const seg & a, const seg & b) {
  return ComputeLineIntersection(a.p, a.q, b.p, b.q);
}

struct line_sweep {
    struct event {
        ld x;
        int tp, id;
        event() {}
        event(ld x, int tp, int id) : x(x), tp(tp), id(id) {}
        bool operator<(const event &e) const {
            if (abs(x - e.x) > EPS) return x < e.x;
            return tp > e.tp;
        }
    };
    set<seg> s;
    vector<set<seg>::iterator> where;
    inline set<seg>::iterator prev(set<seg>::iterator it) {
        return it == s.begin() ? s.end() : --it;
    }
    inline set<seg>::iterator next(set<seg>::iterator it) {
        return ++it;
    }
    vector<pt> solve(const vector<seg> &a, bool find_points=0) {
        vector<pt> res;
        int n = (int) a.size();
        vector<event> e;
        for (int i = 0; i < n; ++i) {
            e.push_back(event(min(a[i].p.x, a[i].q.x), +1, i));
            e.push_back(event(max(a[i].p.x, a[i].q.x), -1, i));
        }
        sort(e.begin(), e.end());
        s.clear();
        where.resize(a.size());
        for (size_t i = 0; i < e.size(); ++i) {
            int id = e[i].id;
            if (e[i].tp == +1) {
                auto nxt = s.lower_bound(a[id]), prv = prev(nxt);
                if (nxt != s.end() && intersect(*nxt, a[id])){
                  pt p = find_intp(*nxt, a[id]);
                  res.pb(p);
                  if(!find_points) return res;
                }
                if (prv != s.end() && intersect(*prv, a[id])){
                  pt p = find_intp(*prv, a[id]);
                  res.pb(p);
                  if(!find_points) return res;
                }
                where[id] = s.insert(nxt, a[id]);
            } else {
                if (where[id] == s.end()) continue;
                auto nxt = next(where[id]), prv = prev(where[id]);
                if (nxt != s.end() && prv != s.end() && intersect(*nxt, *prv)){
                  pt p = find_intp(*nxt, *prv);
                  res.pb(p);
                  if(!find_points) return res;
                }
                s.erase(where[id]);
            }
        }
        return res;
    }
};
