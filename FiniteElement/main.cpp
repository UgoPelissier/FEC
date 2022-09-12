#include <iostream>
#include "Triangulation.hpp"

using namespace std;

int main()
{
    Point p1(0,0), p2(1,1), p3(1,0), p4(0,1);
    
    listPoints l({p1, p2, p3, p4});
    
    l.sortIncreasingTrigonometricOrder().affiche();
    
    return 0;
}
