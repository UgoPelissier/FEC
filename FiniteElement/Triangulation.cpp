#include "Triangulation.hpp"

using namespace std;

/*------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------*/
Point::Point(double x, double y)
    :m_x(x), m_y(y)
{}

bool Point::isEqual(Point *p2) const
{
    if (abs(m_x - p2->m_x)<1e-10 && abs(m_y - p2->m_y)<1e-10)
        return true;
    else
        return false;
}

double Point::getX() const
{
    return m_x;
}

double Point::getY() const
{
    return m_y;
}

double Point::distance(Point const& p2) const
{
    double d = sqrt(pow((m_x-p2.m_x),2)+pow((m_y-p2.m_y),2));
    return d;
}

bool Point::collinear(Point const& p2, Point const& p3) const
{
    double a = m_x * (p2.m_y - p3.m_y) + p2.m_x * (p3.m_y - m_y) + p3.m_x * (m_y - p2.m_y);
    if (abs(a)<1e-10)
        return true;
    else
        return false;
}

bool Point::isBewteen(Point const& p1, Point const& p2) const
{
    if (abs( (this->distance(p1) + this->distance(p2)) - p1.distance(p2) ) < 1e-10) {
        return true;
    }
    return false;
}

void Point::affiche() const
{
    cout << "[" << m_x << " " << m_y << "]";
}

Point::~Point()
{}

/*------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------*/
listPoints::listPoints()
    :m_list(0)
{}

listPoints::listPoints(Point p)
{
    m_list.push_back(new Point(p));
}

listPoints::listPoints(vector<Point> list)
{
    for (size_t i(0); i<list.size(); i++) {
        m_list.push_back(new Point(list[i]));
    }
}

bool listPoints::inList(Point *p) const
{
    for (size_t i(0); i<m_list.size(); i++) {
        if (m_list[i]->isEqual(p)) {
            return true;
        }
    }
    return false;
}

int listPoints::inListIndex(Point *p) const
{
    int index(-1);
    for (size_t i(0); i<m_list.size(); i++) {
        if (m_list[i]->isEqual(p)) {
            index = (int)i;
        }
    }
    return index;
}

vector<double> listPoints::xList() const
{
    vector<double> X;
    for (size_t i(0); i<m_list.size(); i++) {
        X.push_back(m_list[i]->getX());
    }
    return X;
}

vector<double> listPoints::yList() const
{
    vector<double> Y;
    for (size_t i(0); i<m_list.size(); i++) {
        Y.push_back(m_list[i]->getY());
    }
    return Y;
}

Point listPoints::geometricCenter() const
{
    vector<double> X = this->xList();
    vector<double> Y = this->yList();
    
    double x = getAverage(X);
    double y = getAverage(Y);
    
    return Point(x,y);
}

listPoints listPoints::translatedGeometricCenter() const
{
    Point GC = this->geometricCenter();
    listPoints translated;
    
    for (size_t i(0); i<m_list.size(); i++) {
        Point p(m_list[i]->getX()-GC.getX(),m_list[i]->getY()-GC.getY());
        translated.m_list.push_back(new Point(p));
    }
    
    return translated;
}

listPoints listPoints::sortIncreasingTrigonometricOrder() const
{
    listPoints translated = this->translatedGeometricCenter();
    vector<double> angles;
    for (size_t i(0); i<translated.m_list.size(); i++) {
        angles.push_back(atan2(translated.m_list[i]->getY(), translated.m_list[i]->getX()));
    }
    
    vector<size_t> index = argsort(angles);
    listPoints sorted;
    for (size_t i(0); i<index.size(); i++) {
        Point *p = m_list[index[i]];
        sorted.m_list.push_back(new Point(p->getX(),p->getY()));
    }
    return sorted;
}

void listPoints::affiche() const
{
    for (size_t i(0); i<m_list.size(); i++) {
        m_list[i]->affiche();
        cout << endl;
    }
}

listPoints::~listPoints()
{
    for(size_t i(0); i<m_list.size(); ++i)
    {
        delete m_list[i];  //On libère la i-ème case mémoire allouée
        m_list[i] = 0;  //On met le pointeur à 0 pour éviter les soucis
    }
}

/*------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------*/
Edge::Edge(Point p1, Point p2)
{
    m_p1 = new Point(p1);
    m_p2 = new Point(p2);
}

bool Edge::isEqual(Edge const& e) const
{
    if (((e.m_p1->isEqual(m_p1)) and (m_p2->isEqual(e.m_p2))) or ((m_p1->isEqual(e.m_p2)) and (m_p2->isEqual(e.m_p1)))) {
        return true;
    }
    return false;
}

void Edge::affiche() const
{
    m_p1->affiche();
    cout << " ";
    m_p2->affiche();
    cout << endl;
}

Edge::~Edge()
{
    delete m_p1;
    delete m_p2;
}

/*------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------*/
Triangle::Triangle()
    :m_triangle(0)
{}

Triangle::Triangle(vector<Point> triangle)
{
    for (size_t i(0); i<triangle.size(); i++) {
        m_triangle.push_back(new Point(triangle[i]));
    }
}

void Triangle::affiche() const
{
    for (size_t i(0); i<m_triangle.size(); i++) {
        m_triangle[i]->affiche();
        cout << " ";
    }
    cout << endl;
}

Triangle::~Triangle()
{
    for(int i(0); i<m_triangle.size(); ++i)
    {
        delete m_triangle[i];  //On libère la i-ème case mémoire allouée
        m_triangle[i] = 0;  //On met le pointeur à 0 pour éviter les soucis
    }
}

/*------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------*/
Circle::Circle(Point center, double radius)
{
    m_center = new Point(center);
    m_radius = radius;
}

void Circle::affiche() const
{
    cout << "[";
    m_center->affiche();
    cout << ", " << m_radius << "]" << endl;
}

Circle::~Circle()
{
    delete m_center;
}

/*------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------*/
Triangulation::Triangulation()
    :m_triangulation(0)
{}

Triangulation::Triangulation(Triangle t)
{
    m_triangulation.push_back(new Triangle(t));
}

void Triangulation::affiche() const
{
    for(size_t i(0); i<m_triangulation.size(); ++i)
    {
        m_triangulation[i]->affiche();
    }
}
    
Triangulation::~Triangulation()
{
    for(size_t i(0); i<m_triangulation.size(); ++i)
    {
        delete m_triangulation[i];  //On libère la i-ème case mémoire allouée
        m_triangulation[i] = 0;  //On met le pointeur à 0 pour éviter les soucis
    }
}
