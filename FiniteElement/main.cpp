#include <iostream>
#include "Triangulation.hpp"

using namespace std;

int main()
{
    /*
    Point p1(0,0), p2(0,1), p3(0.5,0.25), p4(0.5,1), p5(1,0.25), p6(1,1), p7(1.75,0);
    listPoints l({p1, p2, p3, p4, p5, p6, p7});
    Polygon *geometry = new Polygon(l);
    Geometries *geometries = new Geometries({geometry});
    
    Polygon *oven = new Polygon({Point(-1,-1), Point(1,-1), Point(1,1), Point(-1,1)});
    Polygon *piece = new Polygon({Point(-0.5,-0.2), Point(0.5,-0.2), Point(0.5,0.2), Point(-0.5,0.2)});
    */
    
    /*
    int N = 21;
    Polygon *circle1 = new Polygon(new Circle(Point(-0.75,-0.75),0.05), N);
    Polygon *circle2 = new Polygon(new Circle(Point(0.75,-0.75),0.05), N);
    Polygon *circle3 = new Polygon(new Circle(Point(0.75,0.75),0.05), N);
    Polygon *circle4 = new Polygon(new Circle(Point(-0.75,0.75),0.05), N);
    Geometries *geometries = new Geometries({oven, piece, circle1, circle2, circle3, circle4});
    */
    
    Polygon *domain = new Polygon({Point(-0.2,-0.25), Point(1.2,-0.25), Point(1.2,0.25), Point(-0.2,0.25)});
    Polygon *naca = new Polygon(listPoints(readGeometry("/Users/fp/Desktop/Ugo/Projets/C++/FiniteElement/Mesh/NACA6409.txt")));
    Geometries *geometries = new Geometries({domain, naca});
    
    double minAngle = 25, meshSize = 100;
    Triangulation delaunayTriangulation;
    delaunayTriangulation.triangulation(geometries, minAngle, meshSize);
    
    string const path = "/Users/fp/Desktop/Ugo/Projets/C++/FiniteElement/Mesh/";
    delaunayTriangulation.saveVTU(path);
    
    return 0;
}
