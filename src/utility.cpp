//
// Created by Yucheng Shi on 7/8/20.
//

#include "StaticBTESolver/utility.h"
#include <cstdlib>
#include <cstdio>
#include <cstring>

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
    return 0.5 * cross_prod(p31, p32)->length();
}

double getVolume(const Point& p1, const Point& p2, const Point& p3, const Point& p4) {
    Point p1p(p1 - p4), p2p(p2 - p4), p3p(p3 - p4);
    double det = p1p.x * (p2p.y * p3p.z - p2p.z * p3p.y) - p2p.x * (p1p.y * p3p.z - p1p.z * p3p.y) + p3p.x * (p1p.y * p2p.z - p1p.z * p2p.y);
    return 1.0 / 6 * std::abs(det);
}

double getVolume(const std::shared_ptr<Point>& p1, const std::shared_ptr<Point>& p2, const std::shared_ptr<Point>& p3, const std::shared_ptr<Point>& p4) {
    auto p1p = std::make_shared<Point>(p1->x - p4->x, p1->y - p4->y, p1->z - p4->z);
    auto p2p = std::make_shared<Point>(p2->x - p4->x, p2->y - p4->y, p2->z - p4->z);
    auto p3p = std::make_shared<Point>(p3->x - p4->x, p3->y - p4->y, p3->z - p4->z);
    double det = p1p->x * (p2p->y * p3p->z - p2p->z * p3p->y) - p2p->x * (p1p->y * p3p->z - p1p->z * p3p->y) + p3p->x * (p1p->y * p2p->z - p1p->z * p2p->y);
    return 1.0 / 6 * std::abs(det);
}

double margin(std::vector<double>& a, std::vector<double>& b) {
    double res = 0;
    for (int i = 0; i < a.size(); i++) {
        res = std::max(std::abs(a[i] - b[i]), res);
    }
    return res;
}

// TODO: need to rewrite
std::vector<double> GaussIntegrationPoints(double a, double b, int N) {
    N = N - 1;
    int N1 = N + 1;
    int N2 = N + 2;
    std::vector<double> gauss,y,y0;
    auto *Lp = new double[N1];
    auto **L = new double*[N1];
    for (int i = 0; i < N1; i++) {
        L[i]=new double[N2];
    }
    for (int i=0; i < N1; i++) {
        L[i][0] = 1;
    }
    for (int i=0; i < N1; i++) {
        y.push_back(cos((2 * (N - i) + 1) * PI / (2 * N + 2)));
    }
    for (int i = 0; i < N1; i++) {
        y0.push_back(2);
    }
    while (margin(y, y0) > eps) {
        for (int i=0; i<N1; i++) {
            L[i][1]=y[i];
        }
        for (int j = 2; j < N1 + 1; j++) {
            for (int i=0; i < N1; i++) {
                L[i][j]=((2*j-1)*y[i]*L[i][j-1]-(j-1)*L[i][j-2])/j;
            }
        }
        for (int i = 0; i < N1; i++) {
            Lp[i]=(N2)*( L[i][N1-1]-y[i]*L[i][N2-1])/(1-y[i]*y[i]);
            y0[i]=y[i];
            y[i]=y0[i]-L[i][N2-1]/Lp[i];
        }
    }
    for (int i = 0; i < N1; i++) {
        gauss.push_back((a*(1-y[i])+b*(1+y[i]))/2);
        gauss.push_back((b-a)/((1-y[i]*y[i])*Lp[i]*Lp[i])*N2/N1*N2/N1);
    }
    for (int i = 0; i < N1; i++) {
        delete [] L[i];
    }
    delete [] L;
    delete [] Lp;
    return gauss;
}

int parseLine(char* line){
    // This assumes that a digit will be found and the line ends in " Kb".
    int i = strlen(line);
    const char* p = line;
    while (*p <'0' || *p > '9') p++;
    line[i-3] = '\0';
    i = atoi(p);
    return i;
}

size_t get_host_memory() {
        FILE* file = fopen("/proc/self/status", "r");
        int result = -1;
        char line[128];
        while (fgets(line, 128, file) != nullptr){
            if (strncmp(line, "VmRSS:", 6) == 0){
                result = parseLine(line);
                break;
            }
        }
        fclose(file);
        return result;
}