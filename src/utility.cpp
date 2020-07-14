//
// Created by Yucheng Shi on 7/8/20.
//

#include <memory>
#include "StaticBTESolver/utility.h"

double dot_prod(const Point& p1, const Point& p2) {
    return p1.x * p2.x + p1.y * p2.y + p1.z * p2.z;
}

double dot_prod(const std::shared_ptr<Point>& p1, const std::shared_ptr<Point>& p2) {
    return p1->x * p2->x + p1->y * p2->y + p1->z * p2->z;
}


Point cross_prod(const Point& p1, const Point& p2) {
    return Point(p1.y * p2.z - p2.y * p1.z, p2.x * p1.z - p1.x * p2.z, p1.x * p2.y - p2.x * p1.y);
}

std::shared_ptr<Point> cross_prod(const std::shared_ptr<Point>& p1, const std::shared_ptr<Point>& p2) {
    auto ptr = std::make_shared<Point>(p1->y * p2->z - p2->y * p1->z, p2->x * p1->z - p1->x * p2->z, p1->x * p2->y - p2->x * p1->y);
    return ptr;
}

double getLength(const Point& p1, const Point& p2) {
    return p1.distance(p2);
}

double getLength(const std::shared_ptr<Point>& p1, const std::shared_ptr<Point>& p2) {
    double x = p1->x - p2->x;
    double y = p1->y - p2->y;
    double z = p1->z - p2->z;
    return sqrt(x * x + y * y + z * z);
}


double getArea(const Point& p1, const Point& p2, const Point& p3) {
    return 0.5 * cross_prod(p3 - p1, p2 - p1).length();
}

double getArea(const std::shared_ptr<Point>& p1, const std::shared_ptr<Point>& p2, const std::shared_ptr<Point>& p3) {
    auto p31 = std::make_shared<Point>(p3->x - p1->x, p3->y - p1->y, p3->z - p1->z);
    auto p32 = std::make_shared<Point>(p3->x - p2->x, p3->y - p2->y, p3->z - p2->z);
    return cross_prod(p31, p32)->length();
}

double getVolume(const Point& p1, const Point& p2, const Point& p3, const Point& p4) {
    Point p1p(p1 - p4), p2p(p2 - p4), p3p(p3 - p4);
    double det = p1p.x * (p2p.y * p3p.z - p2p.z * p3p.y) - p2p.x * (p1p.y * p3p.z - p1p.z * p3p.y) + p3p.x * (p1p.y * p2p.z - p1p.z * p2p.y);
    return 1.0 / 6 * abs(det);
}

double getVolume(const std::shared_ptr<Point>& p1, const std::shared_ptr<Point>& p2, const std::shared_ptr<Point>& p3, const std::shared_ptr<Point>& p4) {
    auto p1p = std::make_shared<Point>(p1->x - p4->x, p1->y - p4->y, p1->z - p4->z);
    auto p2p = std::make_shared<Point>(p2->x - p4->x, p2->y - p4->y, p2->z - p4->z);
    auto p3p = std::make_shared<Point>(p3->x - p4->x, p3->y - p4->y, p3->z - p4->z);
    double det = p1p->x * (p2p->y * p3p->z - p2p->z * p3p->y) - p2p->x * (p1p->y * p3p->z - p1p->z * p3p->y) + p3p->x * (p1p->y * p2p->z - p1p->z * p2p->y);
    return 1.0 / 6 * abs(det);
}