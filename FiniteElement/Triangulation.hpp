#ifndef Point_hpp
#define Point_hpp

#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <numeric>
#include <math.h>
#include <valarray>

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

/*------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------*/
class Point
{
public:
    Point(double x, double y);
    virtual bool isEqual(Point *p2) const;
    virtual double getX() const;
    virtual double getY() const;
    virtual double distance(Point const& p2) const;
    virtual bool collinear(Point const& p2, Point const& p3) const;
    virtual bool isBewteen(Point const& p1, Point const& p2) const;
    virtual void affiche() const;
    virtual ~Point();

protected:
    double m_x;
    double m_y;
};

class listPoints
{
public:
    listPoints();
    listPoints(Point p);
    listPoints(vector<Point> list);
    
    virtual bool inList(Point *p) const;
    virtual int inListIndex(Point *p) const;
    
    virtual vector<double> xList() const;
    virtual vector<double> yList() const;
    virtual Point geometricCenter() const;
    virtual listPoints translatedGeometricCenter() const;
    virtual listPoints sortIncreasingTrigonometricOrder() const;
    
    virtual void affiche() const;
    virtual ~listPoints();
    
protected:
    vector<Point*> m_list;
};

class Edge
{
public:
    Edge(Point p1, Point p2);
    virtual bool isEqual(Edge const& e) const;
    virtual void affiche() const;
    virtual ~Edge();
    
protected:
    Point *m_p1;
    Point *m_p2;
};

class Triangle
{
public:
    Triangle();
    Triangle(vector<Point> triangle);
    virtual void affiche() const;
    virtual ~Triangle();
    
protected:
    vector<Point*> m_triangle;
};

class Circle
{
public:
    Circle(Point center, double radius);
    virtual void affiche() const;
    virtual ~Circle();
    
protected:
    Point *m_center;
    double m_radius;
    
};

class Triangulation
{
public:
    Triangulation();
    Triangulation(Triangle t);
    virtual void affiche() const;
    virtual ~Triangulation();
    
protected:
    vector<Triangle*> m_triangulation;
};

#endif /* Point_hpp */
