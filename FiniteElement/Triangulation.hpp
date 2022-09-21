#ifndef Point_hpp
#define Point_hpp

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <numeric>
#include <math.h>
#include <valarray>
#include <exception>


using namespace std;

/*------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------*/
inline double getAverage(vector<double> const& v) {
    if (v.empty()) {
        return 0;
    }
    return accumulate(v.begin(), v.end(), 0.0) / v.size();
}

inline valarray<double> vector2valarray(vector<double> const& v) {
    valarray<double> c(v.size());
    copy(begin(v), end(v), begin(c));
    return c;
}

inline vector<size_t> argsort(vector<double> const& v) {
    vector<pair<double,size_t>> vect;
    for (size_t i(0); i<v.size(); i++) {
            vect.push_back(make_pair(v[i],i));
    }
    sort(vect.begin(), vect.end());
    
    vector<size_t> sortIndex;
    for (size_t i(0); i<v.size(); i++) {
        sortIndex.push_back(vect[i].second);
    }
    return sortIndex;
}

inline bool in(size_t const& d, vector<size_t> const& list) {
    for (size_t i(0); i<list.size(); i++) {
        if (d == list[i]) {
            return true;
        }
    }
    return false;
}

inline double rad2deg(double const& radians) {
    double degrees = radians*180/M_PI;
    return degrees;
}

inline double deg2rad(double const& degrees) {
    double radians = degrees*M_PI/180;
    return radians;
}

inline vector<pair<double,double>> readGeometry(string const& pathFilename) {
    ifstream file(pathFilename);
    string line;
    double x,y;
    vector<pair<double,double>> geometry;
    while (getline(file,line)) {
        stringstream line_stream(line);
        if (line_stream >> x >> y) {
            geometry.push_back(make_pair(x,y));
        }
    }
    return geometry;
}

/*------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------*/
class Point
{
public:
    Point();
    Point(double const& x, double const& y);
    Point(pair<double,double> const& point);
    
    virtual bool isEqual(Point *p2) const;
    virtual double getX() const;
    virtual double getY() const;
    virtual void addX(double d);
    virtual void addY(double d);
    virtual double distance(Point *p2) const;
    virtual bool collinear(Point *p2, Point *p3) const;
    virtual bool isBewteen(Point *p1, Point *p2) const;
    
    virtual void affiche() const;
    virtual ~Point();

protected:
    double m_x;
    double m_y;
};

class Circle
{
public:
    Circle();
    Circle(Point const& center, double const& radius);
    
    virtual bool inside(Point *p);
    
    virtual Point* getCenter() const;
    virtual double getRadius() const;
    
    virtual void affiche() const;
    virtual ~Circle();
    
protected:
    Point *m_center;
    double m_radius;
    
};

class listPoints
{
public:
    listPoints();
    listPoints(Point const& p);
    listPoints(vector<Point> const& list);
    listPoints(vector<Point*> list);
    listPoints(vector<pair<double,double>> const& list);
    
    virtual int listSize() const;
    virtual pair<bool,size_t> inList(Point *p) const;
    virtual void add(Point *p);
    
    virtual vector<double> xList() const;
    virtual vector<double> yList() const;
    virtual Point geometricCenter() const;
    virtual listPoints translatedGeometricCenter() const;
    virtual listPoints sortIncreasingTrigonometricOrder() const;

    virtual vector<Point*> getList() const;
    virtual void affiche() const;
    virtual ~listPoints();
    
protected:
    vector<Point*> m_list;
};

class Edge
{
public:
    Edge();
    Edge(Point const& p1, Point const& p2);
    Edge(Point *p1, Point *p2);
    virtual bool isEqual(Edge *e) const;
    virtual Point midpoint() const;
    virtual double length() const;
    
    virtual pair<bool,size_t> inList(vector<Edge*> list) const;
    virtual Circle circumcircle() const;
    
    virtual pair<size_t,vector<Edge*>> nextEdge(vector<Edge*> list) const;
    
    virtual bool belongs(Point *p) const;
    virtual bool subsegment(Edge *e) const;
    
    virtual Point* getP1() const;
    virtual Point* getP2() const;
    virtual void affiche() const;
    virtual ~Edge();
    
protected:
    Point *m_p1;
    Point *m_p2;
};

class Line
{
public:
    Line();
    Line(double a, double b);
    
    virtual pair<bool, Point> intersect(Line *l) const;
    virtual pair<bool, Point> intersection(Edge *e) const;
    virtual ~Line();
    
protected:
    double m_a;
    double m_b;
};

class Polygon
{
public:
    Polygon();
    Polygon(vector<Edge*> polygon);
    Polygon(vector<Point> list);
    Polygon(vector<Point*> list);
    Polygon(listPoints list);
    Polygon(listPoints *list);
    Polygon(Circle *c, int const& N);
    
    virtual listPoints poly2list() const;
    virtual int rayCastingNumber(Point *p) const;
    virtual bool belongs(Point *p) const;
    virtual bool circumTriangle(vector<Point*> triangle) const;
    virtual bool inside(Point *p) const;
    virtual vector<Edge*> tooLongEdges(double const& criterion) const;
    
    virtual void prepare(double const& criterion);
    virtual vector<Point> superTriangle() const;
    
    virtual bool subsegment(Edge *e) const;
    
    virtual vector<Edge*> getPolygon() const;
    virtual void affiche() const;
    virtual ~Polygon();
    
protected:
    vector<Edge*> m_polygon;
};

class Triangle
{
public:
    Triangle();
    Triangle(vector<Point> const& triangle);
    Triangle(vector<Point*> triangle);
    
    virtual pair<bool,Circle> circumcircle() const;
    virtual pair<Edge*, double> shortestEdge() const;
    virtual pair<Edge*, double> longestEdge() const;
    virtual bool commonEdge(Polygon *p) const;
    virtual bool inscribed(Polygon *p) const;
    virtual bool inside(Polygon *p) const;
    
    virtual double minAngle() const;
    virtual pair<bool,vector<Edge*>> hasSubsegments(Polygon *geometry) const;
    virtual bool quality(double const& minAngle) const;
    
    virtual vector<Point*> getTriangle() const;
    virtual void affiche() const;
    virtual ~Triangle();
    
protected:
    vector<Point*> m_triangle;
};

class Geometries
{
public:
    Geometries();
    Geometries(vector<Polygon*> geometries);
    Geometries(vector<listPoints*> geometries);
    
    virtual bool outside(Point *p) const;
    virtual bool commonEdge(Triangle *t) const;
    
    virtual void prepare(double const& criterion);
    
    virtual vector<Polygon*> getGeometries() const;
    virtual ~Geometries();
    
protected:
    vector<Polygon*> m_geometries;
};

class Triangulation
{
public:
    Triangulation();
    Triangulation(vector<Triangle*> triangulation);
    Triangulation(Triangle const& t);
    
    virtual vector<Triangle*> getTriangulation() const;
    virtual listPoints pointsList() const;
    virtual size_t local2global(listPoints const& points, size_t const& triangle, size_t const& sommet) const;
    virtual vector<vector<size_t>> connectivityTable() const;
    
    virtual pair<Edge*, double> shortestEdge() const;
    virtual pair<Edge*, double> longestEdge() const;
    virtual vector<Edge*> tooLongEdges(double const& citerion) const;
    
    virtual vector<size_t> cavityIndex(Point *p) const;
    virtual vector<size_t> triangleIndexBelongs(Edge *e) const;
    virtual pair<Polygon*,vector<size_t>> cavity(Point *p) const;
    virtual void add(Point *p);
    
    virtual void domainSuperTriangulation(Polygon *domain);
    virtual void removeSuperTriangle(Polygon *domain);
    virtual void removeOutTriangle(Polygon *domain);
    virtual void removeInTriangle(Polygon *geometry);
    
    virtual vector<Edge*> missingSubsegments(Edge *segment) const;
    virtual bool inTriangulation(Edge *segment) const;
    virtual vector<Edge*> missingGeometrySegment(Polygon *geometry) const;
    virtual void constrainedTriangulation(Geometries *geometries);
    
    virtual void domainTriangulation(Polygon *domain);
    
    virtual bool encroached(Edge *e) const;
    virtual vector<Edge*> encroachedSegments(Polygon *geometry) const;
    virtual pair<bool, Edge*> encroachedUpon(Point *p, Triangle *triangle) const;
    virtual bool encroachedCorrection(Polygon *geometry);
    
    virtual pair<vector<size_t>, size_t> skinnyTriangle(double const& minAngle) const;
    virtual bool qualityCorrection(double const& minAngle, Geometries *geometries);
    
    virtual bool meshSizeCorrection(double const& h);
    
    virtual void triangulation(Geometries *geometries, double const& minAngle, double const& h);
    
    virtual void affiche() const;
    virtual void saveVTU(string const& path) const;
    virtual ~Triangulation();
    
protected:
    vector<Triangle*> m_triangulation;
};

#endif /* Point_hpp */
