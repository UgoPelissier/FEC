#include <iostream>
#include "Triangulation.hpp"

using namespace std;

int main()
{
    Point p1(0,0), p2(0,1), p3(0.5,0.25), p4(0.5,1), p5(1,0.25), p6(1,1), p7(1.75,0);
    
    listPoints l({p1, p2, p3, p4, p5, p6, p7});
    
    Polygon *geometry = new Polygon(l);
    
    Geometries *geometries = new Geometries({geometry});
    
    Triangulation delaunayTriangulation;
    delaunayTriangulation.triangulation(geometries, 15, 0.5);
    
    string const path = "/Users/fp/Desktop/Ugo/Projets/C++/FiniteElement/Mesh/";
    delaunayTriangulation.saveVTU(path);
    
    return 0;
}
