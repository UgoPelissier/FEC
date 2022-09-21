#include "Triangulation.hpp"

using namespace std;

/*------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------*/
Point::Point()
    :m_x(0), m_y(0)
{}

Point::Point(double const& x, double const& y)
    :m_x(x), m_y(y)
{}

Point::Point(pair<double,double> const& point)
    :m_x(point.first), m_y(point.second)
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

void Point::addX(double d)
{
    m_x+=d;
}

void Point::addY(double d)
{
    m_y+=d;
}

double Point::distance(Point *p2) const
{
    return sqrt(pow((m_x-p2->m_x),2)+pow((m_y-p2->m_y),2));
}

bool Point::collinear(Point *p2, Point *p3) const
{
    double a = m_x * (p2->m_y - p3->m_y) + p2->m_x * (p3->m_y - m_y) + p3->m_x * (m_y - p2->m_y);
    if (abs(a)<1e-10)
        return true;
    else
        return false;
}

bool Point::isBewteen(Point *p1, Point *p2) const
{
    if (abs( (this->distance(p1) + this->distance(p2)) - p1->distance(p2) ) < 1e-10) {
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
Circle::Circle()
    :m_center(0), m_radius(0)
{}

Circle::Circle(Point const& center, double const& radius)
{
    m_center = new Point(center);
    m_radius = radius;
}

bool Circle::inside(Point *p)
{
    if (p->distance(m_center)<m_radius+1e-10) {
        return true;
    } else {
        return false;
    }
}

Point* Circle::getCenter() const
{
    return m_center;
}

double Circle::getRadius() const
{
    return m_radius;
}

void Circle::affiche() const
{
    cout << "[";
    m_center->affiche();
    cout << ", " << m_radius << "]" << endl;
}

Circle::~Circle()
{
    // delete m_center;
}

/*------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------*/
listPoints::listPoints()
    :m_list(0)
{}

listPoints::listPoints(Point const& p)
{
    m_list.push_back(new Point(p));
}

listPoints::listPoints(vector<Point> const& list)
{
    for (size_t i(0); i<list.size(); i++) {
        m_list.push_back(new Point(list[i]));
    }
}

listPoints::listPoints(vector<Point*> list)
{
    for (size_t i(0); i<list.size(); i++) {
        m_list.push_back(list[i]);
    }
}

listPoints::listPoints(vector<pair<double,double>> const& list)
{
    for (size_t i(0); i<list.size(); i++) {
        m_list.push_back(new Point(list[i]));
    }
}

int listPoints::listSize() const
{
    return (int)m_list.size();
}

pair<bool,size_t> listPoints::inList(Point *p) const
{
    for (size_t i(0); i<m_list.size(); i++) {
        if (m_list[i]->isEqual(p)) {
            return make_pair(true, i);
        }
    }
    return make_pair(false, 0);
}

void listPoints::add(Point *p)
{
    m_list.push_back(p);
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

vector<Point*> listPoints::getList() const
{
    return m_list;
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
        // delete m_list[i];  //On libère la i-ème case mémoire allouée
        // m_list[i] = 0;  //On met le pointeur à 0 pour éviter les soucis
    }
}

/*------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------*/
Edge::Edge()
    : m_p1(0), m_p2(0)
{}

Edge::Edge(Point const& p1, Point const& p2)
{
    m_p1 = new Point(p1);
    m_p2 = new Point(p2);
}

Edge::Edge(Point *p1, Point *p2)
    : m_p1(p1), m_p2(p2)
{}

bool Edge::isEqual(Edge *e) const
{
    if (((e->m_p1->isEqual(m_p1)) and (m_p2->isEqual(e->m_p2))) or ((m_p1->isEqual(e->m_p2)) and (m_p2->isEqual(e->m_p1)))) {
        return true;
    }
    return false;
}

Point Edge::midpoint() const
{
    double xp, yp;
    xp = (m_p1->getX()+m_p2->getX())/2;
    yp = (m_p1->getY()+m_p2->getY())/2;
    return Point(xp,yp);
}

double Edge::length() const
{
    return m_p1->distance(m_p2);
}

pair<bool,size_t> Edge::inList(vector<Edge*> list) const
{
    for (size_t i(0); i<list.size(); i++) {
        if (this->isEqual(list[i])) {
            return make_pair(true,i);
        }
    }
    return make_pair(false,0);
}

Circle Edge::circumcircle() const
{
    return Circle(this->midpoint(), this->length()/2);
}

pair<size_t,vector<Edge*>> Edge::nextEdge(vector<Edge*> list) const
{
    for (size_t i(0); i<list.size(); i++) {
        if (not this->isEqual(list[i])) {
            if (m_p2->isEqual(list[i]->getP1())) {
                return make_pair(i,list);
            }
            else if (m_p2->isEqual(list[i]->getP2())) {
                list[i] = new Edge(list[i]->getP2(), list[i]->getP1());
                return make_pair(i,list);
            }
        }
    }
    return make_pair(size_t(-1), list);
}

bool Edge::belongs(Point *p) const
{
    if (p->isBewteen(m_p1, m_p2)) {
        return true;
    }
    else {
        return false;
    }
}

bool Edge::subsegment(Edge *e) const
{
    if (not this->isEqual(e)) {
        if ((e->belongs(m_p1)) and (e->belongs(m_p2))) {
            return true;
        }
    }
    return false;
}

Point* Edge::getP1() const
{
    return m_p1;
}

Point* Edge::getP2() const
{
    return m_p2;
}

void Edge::affiche() const
{
    cout << "[";
    m_p1->affiche();
    cout << " ";
    m_p2->affiche();
    cout << "]";
    cout << endl;
}

Edge::~Edge()
{
    // delete m_p1;
    // delete m_p2;
}

/*------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------*/
Line::Line()
    :m_a(0), m_b(0)
{}

Line::Line(double a, double b)
    :m_a(a), m_b(b)
{}

pair<bool, Point> Line::intersect(Line *l) const
{
    Point p(0,0);
    if (abs(m_a-l->m_a)<1e-10) {
        if (abs(m_b-l->m_b)<1e-10)
        {
            return make_pair(true, Point(0,m_b));
        }
        else {
            return make_pair(false, Point());
        }
    }
    else {
        double x = (m_b-l->m_b)/(l->m_a-m_a);
        double y = m_a*x+m_b;
        return make_pair(true, Point(x,y));
    }
}

pair<bool, Point> Line::intersection(Edge *e) const
{
    Point p(e->getP1()->getX(),m_a*e->getP1()->getX()+m_b);
    if (abs(e->getP2()->getX()-e->getP1()->getX())<1e-10) {
        if (p.isBewteen(e->getP1(), e->getP2()))
            return make_pair(true, p);
        else
            return make_pair(false, Point());
    }
    else {
        double c = (e->getP2()->getY()-e->getP1()->getY())/(e->getP2()->getX()-e->getP1()->getX());
        double d = e->getP1()->getY() - c*e->getP1()->getX();
        if (this->intersect(new Line(c,d)).first) {
            if (abs(m_a-c)<1e-10) {
                if (abs((m_a*e->getP1()->getX()+m_b)-e->getP1()->getY())<1e-10) {
                    return make_pair(true, Point(e->getP1()->getX(),e->getP1()->getY()));
                }
                else
                    return make_pair(false, Point());
            }
            else {
                if (this->intersect(new Line(c,d)).second.isBewteen(e->getP1(), e->getP2())) {
                    return make_pair(true,this->intersect(new Line(c,d)).second);
                } else {
                    return make_pair(false, Point());
                }
            }
        }
        else
            return make_pair(false, Point());
    }
}

Line::~Line()
{}

/*------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------*/
Polygon::Polygon()
    :m_polygon(0)
{}

Polygon::Polygon(vector<Edge*> polygon)
    :m_polygon(polygon)
{}

Polygon::Polygon(vector<Point> list)
{
    if (not list[0].isEqual(new Point(list[list.size()-1]))) {
        list.push_back(list[0]);
    }
    for (size_t i(0); i<list.size()-1; i++) {
        m_polygon.push_back(new Edge(list[i], list[i+1]));
    }
}

Polygon::Polygon(vector<Point*> list)
{
    if (not list[0]->isEqual(list[list.size()-1])) {
        list.push_back(list[0]);
    }
    for (size_t i(0); i<list.size()-1; i++) {
        m_polygon.push_back(new Edge(list[i], list[i+1]));
    }
}

Polygon::Polygon(listPoints list)
{
    if (not list.getList()[0]->isEqual(list.getList()[list.listSize()-1])) {
        list.add(list.getList()[0]);
    }
    for (size_t i(0); i<list.listSize()-1; i++) {
        m_polygon.push_back(new Edge(list.getList()[i], list.getList()[i+1]));
    }
}

Polygon::Polygon(listPoints *list)
{
    if (not list->getList()[0]->isEqual(list->getList()[list->listSize()-1])) {
        list->add(list->getList()[0]);
    }
    for (size_t i(0); i<list->listSize()-1; i++) {
        m_polygon.push_back(new Edge(list->getList()[i], list->getList()[i+1]));
    }
}

Polygon::Polygon(Circle *c, int const& N)
{
    double r = c->getRadius();
    Point *center = c->getCenter();
    double x,y, theta;
    Point *p1, *p2;
    for (size_t i(0); i<(N-1); i++) {
        theta = (2*M_PI*i)/(N-1);
        x = r*cos(theta)+center->getX();
        y = r*sin(theta)+center->getY();
        p1 = new Point(x,y);
        theta = (2*M_PI*(i+1))/(N-1);
        x = r*cos(theta)+center->getX();
        y = r*sin(theta)+center->getY();
        p2 = new Point(x,y);
        m_polygon.push_back(new Edge(p1,p2));
    }
}

listPoints Polygon::poly2list() const
{
    listPoints list;
    for (size_t i(0); i<m_polygon.size(); i++) {
        list.add(m_polygon[i]->getP1());
        if (i==m_polygon.size()) {
            list.add(m_polygon[i]->getP2());
        }
    }
    return list;
}

int Polygon::rayCastingNumber(Point *p) const
{
    int left(0), right(0);
    Line l(0,p->getY());
    listPoints previous;
    for (size_t i(0); i<m_polygon.size(); i++) {
        bool intersect;
        Point intersectPoint;
        tie(intersect, intersectPoint) = l.intersection(m_polygon[i]);
        if (intersect and (not previous.inList(new Point(intersectPoint)).first)) {
            previous.add(new Point(intersectPoint));
            if (intersectPoint.getX()<p->getX()) {
                left++;
            } else {
                right++;
            }
        }
    }
    return max(left, right);
}

bool Polygon::belongs(Point *p) const
{
    for (size_t i(0); i<m_polygon.size(); i++) {
        if (p->isBewteen(m_polygon[i]->getP1(), m_polygon[i]->getP2())) {
            return true;
        }
    }
    return false;
}

bool Polygon::circumTriangle(vector<Point*> triangle) const
{
    int compteur = 0;
    for (size_t i(0); i<triangle.size(); i++) {
        if (this->belongs(triangle[i])) {
            compteur++;
        }
    }
    if (compteur == triangle.size()) {
        return true;
    } else {
        return false;
    }
}

bool Polygon::inside(Point *p) const
{
    if (this->belongs(p)) {
        return true;
    } else {
        int n = this->rayCastingNumber(p);
        if ((n%2)==0) {
            return false;
        } else {
            return true;
        }
    }
}

vector<Edge*> Polygon::tooLongEdges(double const& criterion) const
{
    vector<Edge*> L;
    for (size_t i(0); i<m_polygon.size(); i++) {
        double d = m_polygon[i]->length();
        if (d>criterion) {
            if (not m_polygon[i]->inList(L).first) {
                L.push_back(m_polygon[i]);
            }
        }
    }
    return L;
}

void Polygon::prepare(double const& criterion)
{
    vector<Edge*> list = this->tooLongEdges(criterion);
    while (not list.empty()) {
        vector<Edge*> newGeometry;
        for (size_t i(0); i<m_polygon.size(); i++) {
            if (not m_polygon[i]->inList(list).first) {
                newGeometry.push_back(m_polygon[i]);
            } else {
                Point *p = new Point(m_polygon[i]->midpoint());
                newGeometry.push_back(new Edge(m_polygon[i]->getP1(), p));
                newGeometry.push_back(new Edge(p, m_polygon[i]->getP2()));
            }
        }
        m_polygon = newGeometry;
        list = this->tooLongEdges(criterion);
    }
}

vector<Point> Polygon::superTriangle() const {
    double eps(0.1);
    vector<double> x = this->poly2list().xList(), y = this->poly2list().yList();
    auto xit = std::minmax_element(x.begin(), x.end()), yit = std::minmax_element(y.begin(), y.end());
    double xmin = *xit.first, xmax = *xit.second, ymin = *yit.first, ymax = *yit.second;
    Point A(xmin-2*eps, ymin-eps), B(xmin+2*(xmax-xmin)+3*eps,ymin-eps), C(xmin-2*eps,ymin+2*(ymax-ymin)+3*eps);
    return {A,B,C};
}

bool Polygon::subsegment(Edge *e) const
{
    for (size_t i(0); i<m_polygon.size(); i++) {
        if (e->subsegment(m_polygon[i])) {
            return true;
        }
    }
    return false;
}

vector<Edge*> Polygon::getPolygon() const
{
    return m_polygon;
}

void Polygon::affiche() const
{
    for (size_t i(0); i<m_polygon.size(); i++) {
        m_polygon[i]->affiche();
    }
}

Polygon::~Polygon()
{
    for(size_t i(0); i<m_polygon.size(); ++i)
    {
        // delete m_polygon[i];  //On libère la i-ème case mémoire allouée
        // m_polygon[i] = 0;  //On met le pointeur à 0 pour éviter les soucis
    }
}

/*------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------*/
Triangle::Triangle()
    :m_triangle(0)
{}

Triangle::Triangle(vector<Point> const& triangle)
{
    for (size_t i(0); i<triangle.size(); i++) {
        m_triangle.push_back(new Point(triangle[i]));
    }
}

Triangle::Triangle(vector<Point*> triangle)
    :m_triangle(triangle)
{}

pair<bool,Circle> Triangle::circumcircle() const
{
    double a, b, c, x1, y1, x2, y2, x3, y3, xp, yp, d, radius;
    x1 = m_triangle[0]->getX();
    y1 = m_triangle[0]->getY();
    x2 = m_triangle[1]->getX();
    y2 = m_triangle[1]->getY();
    x3 = m_triangle[2]->getX();
    y3 = m_triangle[2]->getY();
    a = sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));
    b = sqrt((x3 - x1)*(x3 - x1) + (y3 - y1)*(y3 - y1));
    c = sqrt((x3 - x2)*(x3 - x2) + (y3 - y2)*(y3 - y2));
    if (((a+b+c)*(b+c-a)*(c+a-b)*(a+b-c))<=0) {
        return make_pair(false, Circle());
    } else {
        radius = (a*b*c) / (sqrt((a+b+c)*(b+c-a)*(c+a-b)*(a+b-c)));
        d = 2*(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2));
        xp = ((x1*x1 + y1*y1)*(y2-y3) + (x2*x2 + y2*y2)*(y3-y1) + (x3*x3 + y3*y3)*(y1-y2))/d;
        yp = ((x1*x1 + y1*y1)*(x3-x2) + (x2*x2 + y2*y2)*(x1-x3) + (x3*x3 + y3*y3)*(x2-x1))/d;
        return make_pair(true,Circle(Point(xp,yp), radius));
    }
}

pair<Edge*, double> Triangle::shortestEdge() const
{
    Polygon triangle(m_triangle);
    Edge *shortest = triangle.getPolygon()[0];
    double d = triangle.getPolygon()[0]->length();
    for (size_t i(0); i<triangle.getPolygon().size(); i++) {
        if ((triangle.getPolygon()[i]->length())<d) {
            d = triangle.getPolygon()[i]->length();
            shortest = triangle.getPolygon()[i];
        }
    }
    return make_pair(shortest, d);
}

pair<Edge*, double> Triangle::longestEdge() const
{
    Polygon triangle(m_triangle);
    Edge *longest = triangle.getPolygon()[0];
    double d = triangle.getPolygon()[0]->length();
    for (size_t i(0); i<triangle.getPolygon().size(); i++) {
        if ((triangle.getPolygon()[i]->length())>d) {
            d = triangle.getPolygon()[i]->length();
            longest = triangle.getPolygon()[i];
        }
    }
    return make_pair(longest, d);
}

bool Triangle::commonEdge(Polygon *p) const
{
    Polygon triangle(m_triangle);
    for (size_t i(0); i<triangle.getPolygon().size(); i++) {
        if (triangle.getPolygon()[i]->inList(p->getPolygon()).first) {
            return true;
        }
    }
    return false;
}

bool Triangle::inscribed(Polygon *p) const
{
    if (p->belongs(m_triangle[0]) and p->belongs(m_triangle[1]) and p->belongs(m_triangle[2])) {
        return true;
    }
    return false;
}

bool Triangle::inside(Polygon *p) const
{
    if (p->inside(m_triangle[0]) and p->inside(m_triangle[1]) and p->inside(m_triangle[2])) {
        return true;
    }
    return false;
}

double Triangle::minAngle() const
{
    double radius = this->circumcircle().second.getRadius();
    double d = this->shortestEdge().second;
    double theta = asin(d/(2*radius));
    return theta;
}

pair<bool,vector<Edge*>> Triangle::hasSubsegments(Polygon *geometry) const
{
    Polygon triangle(m_triangle);
    vector<Edge*> list;
    bool has = false;
    for (size_t i(0); i<triangle.getPolygon().size(); i++)  {
        if (geometry->subsegment(triangle.getPolygon()[i])) {
            list.push_back(triangle.getPolygon()[i]);
            has = true;
        }
    }
    return make_pair(has, list);
}

bool Triangle::quality(double const& minAngle) const
{
    double theta = rad2deg(this->minAngle());
    if (theta < minAngle) {
        return false;
    } else {
        return true;
    }
}

vector<Point*> Triangle::getTriangle() const
{
    return m_triangle;
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
    for(int i(0); i<m_triangle.size(); ++i) {
        // delete m_triangle[i];
        // m_triangle[i] = 0;
    }
}

/*------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------*/
Geometries::Geometries()
    :m_geometries(0)
{}

Geometries::Geometries(vector<Polygon*> geometries)
    :m_geometries(geometries)
{}

Geometries::Geometries(vector<listPoints*> geometries)
{
    for (size_t i(0); i<geometries.size(); i++) {
        m_geometries.push_back(new Polygon(geometries[i]));
    }
}

bool Geometries::outside(Point *p) const
{
    for(int i(0); i<m_geometries.size(); ++i) {
        if (m_geometries[i]->inside(p)) {
            return false;
        }
    }
    return true;
}

bool Geometries::commonEdge(Triangle *t) const
{
    for(int i(0); i<m_geometries.size(); ++i) {
        if (t->commonEdge(m_geometries[i])) {
            return true;
        }
    }
    return false;
}

void Geometries::prepare(double const& criterion)
{
    for (size_t i(0); i<m_geometries.size(); i++) {
        m_geometries[i]->prepare(criterion);
    }
}

vector<Polygon*> Geometries::getGeometries() const
{
    return m_geometries;
}

Geometries::~Geometries()
{
    for(int i(0); i<m_geometries.size(); ++i) {
        // delete m_geometries[i];
        // m_geometries[i] = 0;
    }
}

/*------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------*/
Triangulation::Triangulation()
    :m_triangulation(0)
{}

Triangulation::Triangulation(vector<Triangle*> triangulation)
    :m_triangulation(triangulation)
{}

Triangulation::Triangulation(Triangle const& t)
{
    m_triangulation.push_back(new Triangle(t));
}

vector<Triangle*> Triangulation::getTriangulation() const
{
    return m_triangulation;
}

listPoints Triangulation::pointsList() const
{
    listPoints list;
    for (size_t i(0); i<m_triangulation.size(); i++) {
        Triangle *t = m_triangulation[i];
        for (size_t j(0); j<t->getTriangle().size(); j++) {
            if (not list.inList(t->getTriangle()[j]).first) {
                list.add(t->getTriangle()[j]);
            }
        }
    }
    return list;
}

size_t Triangulation::local2global(listPoints const& points, size_t const& triangle, size_t const& sommet) const
{
    size_t index;
    Point *p = m_triangulation[triangle]->getTriangle()[sommet];
    index = points.inList(p).second;
    return index;
}

vector<vector<size_t>> Triangulation::connectivityTable() const
{
    listPoints points = this->pointsList();
    vector<vector<size_t>> cT;
    for (size_t i(0); i<m_triangulation.size(); i++) {
        cT.push_back({local2global(points,i,0),local2global(points,i,1),local2global(points,i,2)});
    }
    return cT;
}

pair<Edge*, double> Triangulation::shortestEdge() const
{
    pair<Edge*, double> shortest = m_triangulation[0]->shortestEdge();
    for (size_t i(0); i<m_triangulation.size(); i++) {
        if (m_triangulation[i]->shortestEdge().second<shortest.second) {
            shortest = m_triangulation[i]->shortestEdge();
        }
    }
    return shortest;
}

pair<Edge*, double> Triangulation::longestEdge() const
{
    pair<Edge*, double> longest = m_triangulation[0]->longestEdge();
    for (size_t i(0); i<m_triangulation.size(); i++) {
        if (m_triangulation[i]->longestEdge().second>longest.second) {
            longest = m_triangulation[i]->longestEdge();
        }
    }
    return longest;
}

vector<Edge*> Triangulation::tooLongEdges(double const& criterion) const
{
    vector<Edge*> L;
    for (size_t i(0); i<m_triangulation.size(); i++) {
        Triangle t(m_triangulation[i]->getTriangle());
        for (size_t j(0); j<t.getTriangle().size()-1; j++) {
            Edge *e = new Edge(t.getTriangle()[j], t.getTriangle()[j+1]);
            double d = e->length();
            if (d>criterion) {
                if (not e->inList(L).first) {
                    L.push_back(e);
                }
            }
        }
    }
    return L;
}

vector<size_t> Triangulation::cavityIndex(Point *p) const
{
    vector<size_t> index;
    for (size_t i(0); i<m_triangulation.size(); i++) {
        pair<bool, Circle> circum = m_triangulation[i]->circumcircle();
        if (circum.first == true) {
            if (circum.second.inside(p)) {
                index.push_back(i);
            }
        }
    }
    return index;
}

vector<size_t> Triangulation::triangleIndexBelongs(Edge *e) const
{
    vector<size_t> index;
    for (size_t i(0); i<m_triangulation.size(); i++) {
        Polygon triangle(m_triangulation[i]->getTriangle());
        if (e->inList(triangle.getPolygon()).first) {
            index.push_back(i);
        }
    }
    return index;
}

pair<Polygon*,vector<size_t>> Triangulation::cavity(Point *p) const
{
    vector<Edge*> cav;
    vector<size_t> index = cavityIndex(p);
    vector<size_t> start;
    for (size_t i(0); i<index.size(); i++) {
        Polygon triangle(m_triangulation[index[i]]->getTriangle());
        for (size_t j(0); j<triangle.getPolygon().size(); j++) {
            Edge *e = triangle.getPolygon()[j];
            vector<size_t> triangleIndex = this->triangleIndexBelongs(e);
            if (triangleIndex.size() == 1) {
                cav.push_back(e);
                if (triangle.inside(p) or triangle.belongs(p)) {
                    start.push_back(cav.size()-1);
                }
            }
            else if (triangleIndex.size() > 1) {
                size_t indexOpposedTriangle(0);
                for (size_t k(0); k<triangleIndex.size(); k++) {
                    if (index[i] != triangleIndex[k]) {
                        indexOpposedTriangle = triangleIndex[k];
                    }
                }
                if (not in(indexOpposedTriangle,index)) {
                    cav.push_back(e);
                    if (triangle.inside(p) or triangle.belongs(p)) {
                        start.push_back(cav.size()-1);
                    }
                }
            }
        }
    }
    if (start.empty()) {
        start.push_back(0);
    }
    if (cav.empty()) {
        return make_pair(new Polygon(cav), index);
    }
    bool pointInCavity(false);
    vector<size_t> sortedIndex;
    vector<Edge*> convexCavity;
    Polygon poly;
    size_t i(0), j(0);
    while (not pointInCavity) {
        sortedIndex.clear();
        convexCavity.clear();
        sortedIndex.push_back(start[j]);
        tie(i,cav) = cav[start[j]]->nextEdge(cav);
        sortedIndex.push_back(i);
        while(i!=start[j]) {
            tie(i, cav) = cav[i]->nextEdge(cav);
            sortedIndex.push_back(i);
        }
        sortedIndex.pop_back();
        for (size_t k(0); k<sortedIndex.size(); k++) {
            convexCavity.push_back(cav[sortedIndex[k]]);
        }
        poly = Polygon(convexCavity);
        pointInCavity = Polygon(convexCavity).inside(p);
        i = 0;
        j++;
    }
    Polygon *convCav = new Polygon(convexCavity);
    vector<size_t> cavIndex;
    for (size_t i(0); i<index.size(); i++) {
        if (m_triangulation[index[i]]->commonEdge(convCav) or m_triangulation[index[i]]->inscribed(convCav) or m_triangulation[index[i]]->inside(convCav)) {
            cavIndex.push_back(index[i]);
        }
    }
    return make_pair(convCav, cavIndex);
}

void Triangulation::add(Point *p)
{
    vector<Triangle*> newTriangulation;
    Polygon *cav;
    vector<size_t> cavIndex;
    tie(cav,cavIndex) = this->cavity(p);
    for (size_t i(0); i<m_triangulation.size(); i++) {
        if (not in(i, cavIndex)) {
            newTriangulation.push_back(m_triangulation[i]);
        }
    }
    for (size_t i(0); i<cav->getPolygon().size(); i++) {
        if (not p->collinear(cav->getPolygon()[i]->getP1(), cav->getPolygon()[i]->getP2())) {
            vector<Point*> newTriangle = {cav->getPolygon()[i]->getP1(), p, cav->getPolygon()[i]->getP2()};
            newTriangulation.push_back(new Triangle(newTriangle));
        }
    }
    m_triangulation = newTriangulation;
}

void Triangulation::domainSuperTriangulation(Polygon *domain)
{
    m_triangulation = {new Triangle(domain->superTriangle())};
    listPoints list = domain->poly2list();
    for (size_t i(0); i<list.getList().size(); i++) {
        this->add(list.getList()[i]);
    }
}

void Triangulation::removeSuperTriangle(Polygon *domain)
{
    vector<Triangle*> newTriangulation;
    for (size_t i(0); i<m_triangulation.size(); i++) {
        if (domain->circumTriangle(m_triangulation[i]->getTriangle())) {
            newTriangulation.push_back(m_triangulation[i]);
        }
    }
    m_triangulation = newTriangulation;
}

void Triangulation::removeOutTriangle(Polygon *domain)
{
    vector<Triangle*> newTriangulation;
    for (size_t i(0); i<m_triangulation.size(); i++) {
        int c(0);
        Polygon triangle = Polygon(m_triangulation[i]->getTriangle());
        for (size_t j(0); j<triangle.getPolygon().size(); j++) {
            Edge *e = triangle.getPolygon()[j];
            Point *mid = new Point(e->midpoint());
            mid->addX(1e-9);
            mid->addY(1e-9);
            if (not domain->inside(mid)) {
                c++;
            }
        }
        if (c == 0) {
            newTriangulation.push_back(m_triangulation[i]);
        }
    }
    m_triangulation = newTriangulation;
}

void Triangulation::removeInTriangle(Polygon *geometry)
{
    vector<Triangle*> newTriangulation;
    for (size_t i(0); i<m_triangulation.size(); i++) {
        if (not m_triangulation[i]->inside(geometry)) {
            newTriangulation.push_back(m_triangulation[i]);
        }
    }
    m_triangulation = newTriangulation;
}

vector<Edge*> Triangulation::missingSubsegments(Edge *segment) const
{
    vector<Edge*> missings;
    vector<Edge*> subsegments;
    listPoints points;
    for (size_t i(0); i<m_triangulation.size(); i++) {
        Polygon triangle = Polygon(m_triangulation[i]->getTriangle());
        for (size_t j(0); j<triangle.getPolygon().size(); j++) {
            Edge *e = triangle.getPolygon()[j];
            if (e->subsegment(segment)) {
                if (not e->inList(subsegments).first) {
                    subsegments.push_back(e);
                }
                if (not points.inList(e->getP1()).first) {
                    points.add(e->getP1());
                }
                if (not points.inList(e->getP2()).first) {
                    points.add(e->getP2());
                }
            }
            if (segment->belongs(e->getP1())) {
                if (not points.inList(e->getP1()).first) {
                    points.add(e->getP1());
                }
            }
        }
    }
    if (points.getList().empty()) {
        missings.push_back(segment);
        return missings;
    }
    vector<double> distance;
    for (size_t i(0); i<points.listSize(); i++) {
        distance.push_back(segment->getP1()->distance(points.getList()[i]));
    }
    vector<size_t> index = argsort(distance);
    listPoints sorted;
    for (size_t i(0); i<index.size(); i++) {
        sorted.add(points.getList()[index[i]]);
    }
    for (size_t i(0); i<sorted.listSize()-1; i++) {
        Edge *e = new Edge(sorted.getList()[i], sorted.getList()[i+1]);
        if ( (not e->inList(subsegments).first) and (not e->inList(missings).first) ) {
            missings.push_back(e);
        }
    }
    return missings;
}

bool Triangulation::inTriangulation(Edge *segment) const
{
    for (size_t i(0); i<m_triangulation.size(); i++) {
        Polygon triangle = Polygon(m_triangulation[i]->getTriangle());
        for (size_t j(0); j<triangle.getPolygon().size(); j++) {
            Edge *e = triangle.getPolygon()[j];
            if (e->isEqual(segment)) {
                return true;
            }
        }
    }
    if (this->missingSubsegments(segment).empty()) {
        return true;
    }
    return false;
}

vector<Edge*> Triangulation::missingGeometrySegment(Polygon *geometry) const
{
    vector<Edge*> missings;
    for (size_t i(0); i<geometry->getPolygon().size(); i++) {
        if (not this->inTriangulation(geometry->getPolygon()[i])) {
            missings.push_back(geometry->getPolygon()[i]);
        }
    }
    return missings;
}

void Triangulation::constrainedTriangulation(Geometries *geometries)
{
    for (size_t j(0); j<geometries->getGeometries().size(); j++) {
        Polygon *geometry = geometries->getGeometries()[j];
        vector<Edge*> missings = this->missingGeometrySegment(geometry);
        while (not missings.empty()) {
            for (size_t i(0); i<missings.size(); i++) {
                Edge *segment = missings[i];
                vector<Edge*> missingsSubsegments = this->missingSubsegments(segment);
                while (not missingsSubsegments.empty()) {
                    for (size_t j(0); j<missingsSubsegments.size(); j++) {
                        Edge *subsegment = missingsSubsegments[j];
                        Point *p = new Point(subsegment->midpoint());
                        this->add(p);
                    }
                    missingsSubsegments = this->missingSubsegments(segment);
                }
            }
            missings = this->missingGeometrySegment(geometry);
        }
    }
    
}

void Triangulation::domainTriangulation(Polygon *domain)
{
    this->domainSuperTriangulation(domain);
    this->removeSuperTriangle(domain);
    this->removeOutTriangle(domain);
    this->constrainedTriangulation(new Geometries({domain}));
}

bool Triangulation::encroached(Edge *e) const
{
    Circle c = e->circumcircle();
    for(size_t i(0); i<m_triangulation.size(); ++i) {
        Triangle *t = m_triangulation[i];
        for (size_t j(0); j<t->getTriangle().size(); j++) {
            Point *p = t->getTriangle()[j];
            if ( (not p->isEqual(e->getP1())) and (not p->isEqual(e->getP2())) and c.inside(p) ) {
                return true;
            }
        }
    }
    return false;
}

vector<Edge*> Triangulation::encroachedSegments(Polygon *geometry) const
{
    vector<Edge*> segments;
    for(size_t i(0); i<m_triangulation.size(); ++i) {
        Polygon t = Polygon(m_triangulation[i]->getTriangle());
        for (size_t j(0); j<t.getPolygon().size(); j++) {
            Edge *e = t.getPolygon()[j];
            if ( (this->encroached(e)) and ( (e->inList(geometry->getPolygon()).first) or (geometry->subsegment(e)) ) ) {
                if (not e->inList(segments).first) {
                    segments.push_back(e);
                }
            }
        }
    }
    return segments;
}

pair<bool, Edge*> Triangulation::encroachedUpon(Point *p, Triangle *triangle) const
{
    Polygon tri = Polygon(triangle->getTriangle());
    for(size_t i(0); i<m_triangulation.size(); ++i) {
        Polygon t = Polygon(m_triangulation[i]->getTriangle());
        for (size_t j(0); j<t.getPolygon().size(); j++) {
            Edge *e = t.getPolygon()[j];
            Circle c = e->circumcircle();
            if ( (c.inside(p)) and (not e->inList(tri.getPolygon()).first) ) {
                return make_pair(true, e);
            }
        }
    }
    return make_pair(false, new Edge(0,0));
}

bool Triangulation::encroachedCorrection(Polygon *geometry)
{
    vector<Edge*> segments = this->encroachedSegments(geometry);
    if (segments.empty()) {
        return true;
    }
    else {
        while (not segments.empty()) {
            this->add(new Point(segments[0]->midpoint()));
            segments = this->encroachedSegments(geometry);
        }
    }
    return false;
}

pair<vector<size_t>, size_t> Triangulation::skinnyTriangle(double const& minAngle) const
{
    vector<size_t> index;
    size_t skinny(0);
    double theta = m_triangulation[0]->minAngle();
    for(size_t i(0); i<m_triangulation.size(); ++i) {
        if (not m_triangulation[i]->quality(minAngle)) {
            index.push_back(i);
            if (m_triangulation[i]->minAngle()<theta) {
                theta = m_triangulation[i]->minAngle();
                skinny = i;
            }
        }
    }
    return make_pair(index, skinny);
}

bool Triangulation::qualityCorrection(double const& minAngle, Geometries *geometries)
{
    vector<size_t> index;
    size_t skinny;
    tie(index, skinny) = this->skinnyTriangle(minAngle);
    if (index.empty()) {
        return true;
    }
    while (not index.empty()) {
        Triangle *triangle = m_triangulation[skinny];
        Point *p = triangle->circumcircle().second.getCenter();
        bool encroach;
        Edge *e;
        tie(encroach, e) = this->encroachedUpon(p, triangle);
        if (encroach) {
            this->add(new Point(e->midpoint()));
        }
        else {
            if (geometries->outside(p)) {
                bool correction = false;;
                Polygon t = Polygon(triangle->getTriangle());
                for (size_t i(0); i<t.getPolygon().size(); i++) {
                    Edge *e = t.getPolygon()[i];
                    for (size_t j(0); j<geometries->getGeometries().size(); j++) {
                        Polygon *geometry = geometries->getGeometries()[j];
                        if (geometry->subsegment(e)) {
                            this->add(new Point(e->midpoint()));
                            correction = true;
                        }
                    }
                }
                if (not correction) {
                    for (size_t i(0); i<t.getPolygon().size(); i++) {
                        Edge *e = t.getPolygon()[i];
                        this->add(new Point(e->midpoint()));
                    }
                }
            }
            else {
                this->add(p);
            }
        }
        tie(index, skinny) = this->skinnyTriangle(minAngle);
    }
    return false;
}

bool Triangulation::meshSizeCorrection(double const& h)
{
    Edge *e = this->longestEdge().first;
    vector<Edge*> list = this->tooLongEdges(h);
    if (list.empty()) {
        return true;
    }
    while (not list.empty()) {
        this->add(new Point(e->midpoint()));
        e = this->longestEdge().first;
        list = this->tooLongEdges(h);
    }
    return false;
}

void Triangulation::triangulation(Geometries *geometries, double const& minAngle, double const& h)
{
    Polygon *domain = geometries->getGeometries()[0];
    this->domainTriangulation(domain);
    for (size_t i(1); i<geometries->getGeometries().size(); i++) {
        Polygon *geometry = geometries->getGeometries()[i];
        for (size_t j(0); j<geometry->getPolygon().size(); j++) {
            this->add(geometry->getPolygon()[j]->getP1());
        }
    }
    bool state = (this->qualityCorrection(minAngle, geometries) and this->meshSizeCorrection(h));
    while (not state) {
        for (size_t i(0); i<geometries->getGeometries().size(); i++) {
            state = (this->qualityCorrection(minAngle, geometries) and this->meshSizeCorrection(h));
        }
    }
    this->constrainedTriangulation(geometries);
    for (size_t i(1); i<geometries->getGeometries().size(); i++) {
        Polygon *geometry = geometries->getGeometries()[i];
        this->removeInTriangle(geometry);
    }
}

void Triangulation::affiche() const
{
    for(size_t i(0); i<m_triangulation.size(); ++i)
    {
        m_triangulation[i]->affiche();
    }
}

void Triangulation::saveVTU(string const& path) const
{
    vector<Point*> points = this->pointsList().getList();
    vector<vector<size_t>> cT = this->connectivityTable();
    string const name = "mesh.vtu";
    string const filename = path + name;
    ofstream f(filename.c_str());
    f << "<VTKFile type='UnstructuredGrid' version='0.1' byte_order='LittleEndian'>" << endl;
    f << "  <UnstructuredGrid>" << endl;
    f << "    <Piece NumberOfPoints='" << points.size() << "' "<< " NumberOfCells='" << m_triangulation.size() << "'>" << endl;
    f << "      <Points>" << endl;
    f << "        <DataArray type = 'Float32' NumberOfComponents='3' format='ascii'>" << endl;
    for (size_t i(0); i<points.size(); i++) {
        f << points[i]->getX() << " " << points[i]->getY() << " " << 0.0 << endl;
    }
    f << "        </DataArray>" << endl;
    f << "      </Points>" << endl;
    f << "      <Cells>" << endl;
    f << "        <DataArray type='Int32' Name='connectivity' format='ascii'>" << endl;
    for (size_t i(0); i<m_triangulation.size(); i++) {
        f << cT[i][0] << " " << cT[i][1] << " " << cT[i][2] << endl;
    }
    f << "        </DataArray>" << endl;
    f << "        <DataArray type='Int32' Name='offsets' format='ascii'>" << endl;
    for (size_t i(0); i<m_triangulation.size(); i++) {
        f << 3*int(i+1) << endl;
    }
    f << "        </DataArray>" << endl;
    f << "        <DataArray type='UInt8' Name='types' format='ascii'>" << endl;
    for (size_t i(0); i<m_triangulation.size(); i++) {
        f << 5 << endl;
    }
    f << "        </DataArray>" << endl;
    f << "      </Cells>" << endl;
    f << "<PointData>" << endl;
    f << "  <DataArray type = 'Float32' Name='Erreur' format='ascii'>" << endl;
    for (size_t i(0); i<points.size(); i++) {
        f << 0 << endl;
    }
    f << "  </DataArray>" << endl;
    f << "</PointData>" << endl;
    f << "<CellData>" << endl;
    f << "</CellData>" << endl;
    f << "<UserData>" << endl;
    f << "  <DataArray type = 'Float32' Name='CompteurTemps' format='ascii'>" << endl;
    f << 0 << endl;
    f << "  </DataArray>" << endl;
    f << "  <DataArray type = 'Float32' Name='Temps' format='ascii'>" << endl;
    f << 0 << endl;
    f << "  </DataArray>" << endl;
    f << "</UserData>" << endl;
    f << "    </Piece>" << endl;
    f << "  </UnstructuredGrid>" << endl;
    f << "</VTKFile>" << endl;
    f.close();
}
    
Triangulation::~Triangulation()
{
    for(size_t i(0); i<m_triangulation.size(); ++i)
    {
        // delete m_triangulation[i];  //On libère la i-ème case mémoire allouée
        // m_triangulation[i] = 0;  //On met le pointeur à 0 pour éviter les soucis
    }
}
