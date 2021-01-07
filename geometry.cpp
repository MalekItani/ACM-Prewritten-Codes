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

using namespace std;


#define real ld

ld EPS = 1e-12;
ld PI = acos(-1.0);

struct pt {
  real x, y;
  pt() {}
  pt(real x, real y) : x(x), y(y) {}
  pt(const pt &p) : x(p.x), y(p.y)    {}
  pt operator + (const pt &p)  const { return pt(x+p.x, y+p.y); }
  pt operator - (const pt &p)  const { return pt(x-p.x, y-p.y); }
  pt operator * (real c)     const { return pt(x*c,   y*c  ); }
  pt operator / (real c)     const { return pt(x/c,   y/c  ); }
  real norm(){return sqrt(x * x + y * y);}
  bool operator == (const pt &p)     const { return fabs(x - p.x) < EPS && fabs(y - p.y) < EPS; }
  bool operator<(pt &p) { return tie(x, y) < tie(p.x, p.y); }
};

real dot(pt p, pt q)     { return p.x*q.x+p.y*q.y; }
real dist2(pt p, pt q)   { return dot(p-q,p-q); }
real cross(pt p, pt q)   { return p.x*q.y-p.y*q.x; }
ostream &operator<<(ostream &os, const pt &p) {
  return os << "(" << p.x << "," << p.y << ")";
}

// rotate a point CCW or CW around the origin
pt RotateCCW90(pt p)   { return pt(-p.y,p.x); }
pt RotateCW90(pt p)    { return pt(p.y,-p.x); }
pt RotateCCW(pt p, real t) {
  return pt(p.x*cos(t)-p.y*sin(t), p.x*sin(t)+p.y*cos(t));
}


// Basic Line, use when you care about evaluating points given some X.
struct Line{
  // ax + by + c = 0
  real a, b, c;
  pt p1, p2;
  Line(){}
  Line(pt p1, pt p2):p1(p1), p2(p2){
    a = p2.y - p1.y;
    b = p1.x - p2.x;
    c = -a * p1.x - b*p1.y;
  }
  Line(real slope, real y_int){
    a = -slope;
    b = 1;
    c = -y_int;
  }
  bool contains(pt p){
    return fabs(a * p.x + b * p.y + c) < EPS;
  }

  // USE ONLY WHEN REAL == (LONG) DOUBLE
  pt eval(real x){
    return pt(x, (-c - a * x)/b);
  }
  real slope(){
    if(b == 0) return INF;
    return -a/b;
  }
  bool operator<(Line &L){
    return slope() < L.slope();
  }
};

// Parametrized Line (PLine)
// Use for lines in higher dimensions or when you care about checking if points belong to a line.
struct PLine{
  // L: A + Bt, t in R
  pt A, B;
  pt N; // Normal
  pt _p2; // (OPTIONAL) Secondary point
  PLine(){}
  PLine(pt p1, pt p2){
    A = p1;
    B = p2 - p1;
    // OPTIONAL
    N = pt(-B.y, B.x);
    _p2 = p2;
  }
  bool contains(pt p){
    return fabs(cross(p - A, B-A)) < EPS;
  }
  // USE ONLY WHEN REAL == (LONG) DOUBLE
  // Compute point given parameter.
  pt eval(real t){
    return A + B * t;
  }
  // Compute distance from line to point
  real dist(pt p){
    return fabs(dot(p-A, N)/N.norm());
  }
  bool operator<(PLine &L){
    return B < L.B;
  }
  pt intersection(PLine &L){
    real D = cross(B, L.B);
    real Dv = cross(L.A - A, B);
    return L.eval(Dv/D);
  }
};

struct seg : PLine{
  real mag;
  seg(){}
  seg(pt p1, pt p2):PLine(p1, p2){
    mag = (p2 - p1).norm();
  }
};

// Checks if a point p is inside or on the border of a CONVEX polygon. O(logn)
// ASSUMES POLYGON IS CW.
bool PointInConvexPol(vector<pt>& l, pt p) {
  auto vec = [](pt a,pt b){
    return a.x*b.y-a.y*b.x;
  };
  int a = 1, b = sz(l)-1, c;
  if (vec(l[a] - l[0], l[b] - l[0]) > 0) swap(a,b);
  if (vec(l[a] - l[0], p - l[0]) > 0 || vec(l[b] - l[0], p - l[b]) < 0) return 0;
  while(abs(a-b) > 1) {
    c = (a+b)/2;
    if (vec(l[c] - l[0], p - l[0]) > 0) b = c; else a = c;
  }
  return vec(l[b] - l[a], p - l[a]) <= 0;
}


// returns the 2D convex hull in counterclockwise direction.
// consecutive vertices are not collinear.
// also works with pii instead of pt, after just changing x->fi, y->se
// and pt->pii in function right. Then you don't have to define pt, == or <
vector<pt> convexhull(vector<pt> ab) {
  auto right = [](pt a, pt b, pt c){return cross(b-a, c-a) <= EPS;};
  sort(ab.begin(),ab.end());
  ab.erase(unique(ab.begin(),ab.end()),ab.end());
  int l=sz(ab), i, j, k;
  vector<pt> res(l+1);
  if (l<3) return ab;
  j=0;
  for(i=0;i<l;i++) {
    while (j-2>=0 && right(res[j-2],res[j-1],ab[i])) j--;
    res[j++]=ab[i];
  }
  k=j;
  for(i=l-2;i>=0;i--) {
  while (j-1>=k && right(res[j-2],res[j-1],ab[i])) j--;
  res[j++]=ab[i];
  }
  if (res[j-1]==res[0]) j--;
  return vector<pt>(res.begin(),res.begin()+j);
}

// Distance between angle a and angle b. Where a, b in [0, 2*Pi]
// Directed means the angle has to go from a to b.
real distAnAn(real a, real b, bool directed = 0){
  real p = fmod(2 * PI + fmod(b - a , 2 * PI), 2 * PI);
  if(directed) p = p > PI ? 2*PI - p : p;
  return p;
}

// Returns true if phi is inside the rays theta = a & theta = b (In this direction)
bool anBetween(real phi, real a, real b){
  if (a < b) return a <= phi && phi <= b;
  return a <= phi || phi <= b;
}


// project point c onto line through a and b
// assuming a != b
pt ProjectPointLine(pt a, pt b, pt c) {
  return a + (b-a)*dot(c-a, b-a)/dot(b-a, b-a);
}

// project point c onto line segment through a and b
pt ProjectPointSegment(pt a, pt b, pt c) {
  double r = dot(b-a,b-a);
  if (fabs(r) < EPS) return a;
  r = dot(c-a, b-a)/r;
  if (r < 0) return a;
  if (r > 1) return b;
  return a + (b-a)*r;
}

// compute distance from c to segment between a and b
double DistancePointSegment(pt a, pt b, pt c) {
  return sqrt(dist2(c, ProjectPointSegment(a, b, c)));
}

// compute distance between point (x,y,z) and plane ax+by+cz=d
double DistancePointPlane(double x, double y, double z,
                          double a, double b, double c, double d)
{
  return fabs(a*x+b*y+c*z-d)/sqrt(a*a+b*b+c*c);
}

// determine if lines from a to b and c to d are parallel or collinear
bool LinesParallel(pt a, pt b, pt c, pt d) {
  return fabs(cross(b-a, c-d)) < EPS;
}

bool LinesCollinear(pt a, pt b, pt c, pt d) {
  return LinesParallel(a, b, c, d)
      && fabs(cross(a-b, a-c)) < EPS
      && fabs(cross(c-d, c-a)) < EPS;
}

// determine if line segment from a to b intersects with
// line segment from c to d
bool SegmentsIntersect(pt a, pt b, pt c, pt d) {
  if (LinesCollinear(a, b, c, d)) {
    if (dist2(a, c) < EPS || dist2(a, d) < EPS ||
      dist2(b, c) < EPS || dist2(b, d) < EPS) return true;
    if (dot(c-a, c-b) > 0 && dot(d-a, d-b) > 0 && dot(c-b, d-b) > 0)
      return false;
    return true;
  }
  if (cross(d-a, b-a) * cross(c-a, b-a) > 0) return false;
  if (cross(a-c, d-c) * cross(b-c, d-c) > 0) return false;
  return true;
}

// compute intersection of line passing through a and b
// with line passing through c and d, assuming that unique
// intersection exists; for segment intersection, check if
// segments intersect first
pt ComputeLineIntersection(pt a, pt b, pt c, pt d) {
  b=b-a; d=c-d; c=c-a;
  assert(dot(b, b) > EPS && dot(d, d) > EPS);
  return a + b*cross(c, d)/cross(b, d);
}

// compute center of circle given three points
pt ComputeCircleCenter(pt a, pt b, pt c) {
  b=(a+b)/2;
  c=(a+c)/2;
  return ComputeLineIntersection(b, b+RotateCW90(a-b), c, c+RotateCW90(a-c));
}

// determine if point is in a possibly non-convex polygon (by William
// Randolph Franklin); returns 1 for strictly interior points, 0 for
// strictly exterior points, and 0 or 1 for the remaining points.
// Note that it is possible to convert this into an *exact* test using
// integer arithmetic by taking care of the division appropriately
// (making sure to deal with signs properly) and then by writing exact
// tests for checking point on polygon boundary
bool PointInPolygon(const vector<pt> &p, pt q) {
  bool c = 0;
  for (int i = 0; i < p.size(); i++){
    int j = (i+1)%p.size();
    if ((p[i].y <= q.y && q.y < p[j].y ||
      p[j].y <= q.y && q.y < p[i].y) &&
      q.x < p[i].x + (p[j].x - p[i].x) * (q.y - p[i].y) / (p[j].y - p[i].y))
      c = !c;
  }
  return c;
}

// determine if point is on the boundary of a polygon
bool PointOnPolygon(const vector<pt> &p, pt q) {
  for (int i = 0; i < p.size(); i++)
    if (dist2(ProjectPointSegment(p[i], p[(i+1)%p.size()], q), q) < EPS)
      return true;
    return false;
}

// compute intersection of line through points a and b with
// circle centered at c with radius r > 0
vector<pt> CircleLineIntersection(pt a, pt b, pt c, double r) {
  vector<pt> ret;
  b = b-a;
  a = a-c;
  double A = dot(b, b);
  double B = dot(a, b);
  double C = dot(a, a) - r*r;
  double D = B*B - A*C;
  if (D < -EPS) return ret;
  ret.push_back(c+a+b*(-B+sqrt(D+EPS))/A);
  if (D > EPS)
    ret.push_back(c+a+b*(-B-sqrt(D))/A);
  return ret;
}


// compute intersection of circle centered at a with radius r
// with circle centered at b with radius R
vector<pt> CircleCircleIntersection(pt a, pt b, real r, real R) {
  vector<pt> ret;
  real d = sqrt(dist2(a, b));
  if (d > r+R || d+min(r, R) < max(r, R)) return ret;
  real x = (d*d-R*R+r*r)/(2*d);
  real y = sqrt(r*r-x*x);
  pt v = (b-a)/d;
  ret.push_back(a+v*x + RotateCCW90(v)*y);
  if (y > 0)
    ret.push_back(a+v*x - RotateCCW90(v)*y);
  return ret;
}

// This code computes the area or centroid of a (possibly nonconvex)
// polygon, assuming that the coordinates are listed in a clockwise or
// counterclockwise fashion.  Note that the centroid is often known as
// the "center of gravity" or "center of mass".
double ComputeSignedArea(const vector<pt> &p) {
  double area = 0;
  for(int i = 0; i < p.size(); i++) {
    int j = (i+1) % p.size();
    area += p[i].x*p[j].y - p[j].x*p[i].y;
  }
  return area / 2.0;
}

double ComputeArea(const vector<pt> &p) {
  return fabs(ComputeSignedArea(p));
}

pt ComputeCentroid(const vector<pt> &p) {
  pt c(0,0);
  double scale = 6.0 * ComputeSignedArea(p);
  for (int i = 0; i < p.size(); i++){
    int j = (i+1) % p.size();
    c = c + (p[i]+p[j])*(p[i].x*p[j].y - p[j].x*p[i].y);
  }
  return c / scale;
}

// tests whether or not a given polygon (in CW or CCW order) is simple
bool IsSimple(const vector<pt> &p) {
  for (int i = 0; i < p.size(); i++) {
    for (int k = i+1; k < p.size(); k++) {
      int j = (i+1) % p.size();
      int l = (k+1) % p.size();
      if (i == l || j == k) continue;
      if (SegmentsIntersect(p[i], p[j], p[k], p[l]))
        return false;
    }
  }
  return true;
}

int main() {
  cout << fixed << setprecision(14);

  // Line test:
  Line L(pt(0,0), pt(2, 2));

  assert(L.contains(pt(1, 1)));
  assert(L.contains(pt(1e9, 1e9)));
  assert(L.contains(pt(-1e5, -1e5)));

  cout << L.slope() << endl;

  Line L2(pt(-2, 7), pt(6, 3));

  assert(L2.contains(pt(0, 6)));
  assert(L2.contains(pt(17900, -8944)));

  cout << L2.slope() << endl;

  // Test when real == long double
  Line L3(7.45, -17.2);

  cout << L3.slope() << endl;

  cout << L3.eval(1) << endl;
  cout << L3.eval(1e5) << endl;
  cout << L3.eval(1e9) << endl;

  cout << L3.eval(-1) << endl;
  cout << L3.eval(-1e5) << endl;
  cout << L3.eval(-1e9) << endl;


  // PLine test:
  PLine PL(pt(0,0), pt(2, 2));

  PLine PL2(pt(0,0.5), pt(-1,1.5));

  assert(PL.contains(pt(1, 1)));
  assert(PL.contains(pt(1e9, 1e9)));
  assert(PL.contains(pt(-1e5, -1e5)));

  cout << PL.intersection(PL2) << endl;

  cout << ComputeLineIntersection(PL.A, PL._p2, PL2.A, PL2._p2) << endl;

  auto disll = [](pt p, pt a, pt b){
    pt x = ProjectPointLine(a, b, p);
    return (p-x).norm();
  };

  cout << disll(pt{12, 13.5}, PL.A, PL._p2) << endl;
  cout << PL.dist(pt{12, 13.5}) << endl;


  // Point in poly test
  vector<pt> poly = {pt(0,0), pt(0,2), pt(2,2), pt(2, 0)};
  assert(point_in_convex(poly, pt(0,0)));
  assert(point_in_convex(poly, pt(1,1)));
  assert(point_in_convex(poly, pt(1,0)));
  assert(!point_in_convex(poly, pt(2.1,1)));
  assert(point_in_convex(poly, pt(2,2)));
  assert(point_in_convex(poly, pt(0.1,0.4)));

  poly = {pt(0,0), pt(-1,1), pt(-5,10), pt(3, 10), pt(4, 5), pt(2, 1)};
  assert(point_in_convex(poly, pt(1.6,0.8)));
  assert(point_in_convex(poly, pt(1, 0.5)));
  assert(point_in_convex(poly, pt(-0.6, 0.6)));
  assert(point_in_convex(poly, pt(1, 0.5)));
  assert(!point_in_convex(poly, pt(1.2, 0.5)));
  assert(point_in_convex(poly, pt(0.4, 0.6)));
  assert(point_in_convex(poly, pt(2, 10)));
  assert(point_in_convex(poly, pt(-2, 6)));
  assert(!point_in_convex(poly, pt(2, 10.1)));
  assert(!point_in_convex(poly, pt(4, 1)));

  // expected: (-5,2)
  cerr << RotateCCW90(pt(2,5)) << endl;

  // expected: (5,-2)
  cerr << RotateCW90(pt(2,5)) << endl;

  // expected: (-5,2)
  cerr << RotateCCW(pt(2,5),M_PI/2) << endl;

  // expected: (5,2)
  cerr << ProjectPointLine(pt(-5,-2), pt(10,4), pt(3,7)) << endl;

  // expected: (5,2) (7.5,3) (2.5,1)
  cerr << ProjectPointSegment(pt(-5,-2), pt(10,4), pt(3,7)) << " "
       << ProjectPointSegment(pt(7.5,3), pt(10,4), pt(3,7)) << " "
       << ProjectPointSegment(pt(-5,-2), pt(2.5,1), pt(3,7)) << endl;

  // expected: 6.78903
  cerr << DistancePointPlane(4,-4,3,2,-2,5,-8) << endl;

  // expected: 1 0 1
  cerr << LinesParallel(pt(1,1), pt(3,5), pt(2,1), pt(4,5)) << " "
       << LinesParallel(pt(1,1), pt(3,5), pt(2,0), pt(4,5)) << " "
       << LinesParallel(pt(1,1), pt(3,5), pt(5,9), pt(7,13)) << endl;

  // expected: 0 0 1
  cerr << LinesCollinear(pt(1,1), pt(3,5), pt(2,1), pt(4,5)) << " "
       << LinesCollinear(pt(1,1), pt(3,5), pt(2,0), pt(4,5)) << " "
       << LinesCollinear(pt(1,1), pt(3,5), pt(5,9), pt(7,13)) << endl;

  // expected: 1 1 1 0
  cerr << SegmentsIntersect(pt(0,0), pt(2,4), pt(3,1), pt(-1,3)) << " "
       << SegmentsIntersect(pt(0,0), pt(2,4), pt(4,3), pt(0,5)) << " "
       << SegmentsIntersect(pt(0,0), pt(2,4), pt(2,-1), pt(-2,1)) << " "
       << SegmentsIntersect(pt(0,0), pt(2,4), pt(5,5), pt(1,7)) << endl;

  // expected: (1,2)
  cerr << ComputeLineIntersection(pt(0,0), pt(2,4), pt(3,1), pt(-1,3)) << endl;

  // expected: (1,1)
  cerr << ComputeCircleCenter(pt(-3,4), pt(6,1), pt(4,5)) << endl;

  vector<pt> v;
  v.push_back(pt(0,0));
  v.push_back(pt(5,0));
  v.push_back(pt(5,5));
  v.push_back(pt(0,5));

  // expected: 1 1 1 0 0
  cerr << PointInPolygon(v, pt(2,2)) << " "
       << PointInPolygon(v, pt(2,0)) << " "
       << PointInPolygon(v, pt(0,2)) << " "
       << PointInPolygon(v, pt(5,2)) << " "
       << PointInPolygon(v, pt(2,5)) << endl;

  // expected: 0 1 1 1 1
  cerr << PointOnPolygon(v, pt(2,2)) << " "
       << PointOnPolygon(v, pt(2,0)) << " "
       << PointOnPolygon(v, pt(0,2)) << " "
       << PointOnPolygon(v, pt(5,2)) << " "
       << PointOnPolygon(v, pt(2,5)) << endl;

  // expected: (1,6)
  //           (5,4) (4,5)
  //           blank line
  //           (4,5) (5,4)
  //           blank line
  //           (4,5) (5,4)
  vector<pt> u = CircleLineIntersection(pt(0,6), pt(2,6), pt(1,1), 5);
  for (int i = 0; i < u.size(); i++) cerr << u[i] << " "; cerr << endl;
  u = CircleLineIntersection(pt(0,9), pt(9,0), pt(1,1), 5);
  for (int i = 0; i < u.size(); i++) cerr << u[i] << " "; cerr << endl;
  u = CircleCircleIntersection(pt(1,1), pt(10,10), 5, 5);
  for (int i = 0; i < u.size(); i++) cerr << u[i] << " "; cerr << endl;
  u = CircleCircleIntersection(pt(1,1), pt(8,8), 5, 5);
  for (int i = 0; i < u.size(); i++) cerr << u[i] << " "; cerr << endl;
  u = CircleCircleIntersection(pt(1,1), pt(4.5,4.5), 10, sqrt(2.0)/2.0);
  for (int i = 0; i < u.size(); i++) cerr << u[i] << " "; cerr << endl;
  u = CircleCircleIntersection(pt(1,1), pt(4.5,4.5), 5, sqrt(2.0)/2.0);
  for (int i = 0; i < u.size(); i++) cerr << u[i] << " "; cerr << endl;

  // area shoureal be 5.0
  // centroid shoureal be (1.1666666, 1.166666)
  pt pa[] = { pt(0,0), pt(5,0), pt(1,1), pt(0,5) };
  vector<pt> p(pa, pa+4);
  pt c = ComputeCentroid(p);
  cerr << "Area: " << ComputeArea(p) << endl;
  cerr << "Centroid: " << c << endl;

  return 0;
}
