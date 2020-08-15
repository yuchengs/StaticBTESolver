//
// Created by Yucheng Shi on 7/6/20.
//

#ifndef BTESOLVER_UTILITY_H
#define BTESOLVER_UTILITY_H

#include <climits>
#include <vector>
#include <string>
#include <cmath>
#include <memory>
#include <iomanip>
#include <algorithm>
#include <iostream>

namespace staticbtesolver {
#define PI M_PI
#define eps pow(2,-52)

    struct Point {
        double x, y, z;

        explicit Point(double xx = 0, double yy = 0, double zz = 0) : x{xx}, y{yy}, z{zz} {}

        Point(const Point &pt) = default;

        Point &operator=(const Point &rhs) = default;

        Point operator+(const Point &rhs) const {
            return Point(x + rhs.x, y + rhs.y, z + rhs.z);
        }

        Point operator-(const Point &rhs) const {
            return Point(x - rhs.x, y - rhs.y, z - rhs.z);
        }

        Point operator/(const double rhs) const {
            return Point(x / rhs, y / rhs, z / rhs);
        }

        double operator*(const Point &rhs) const {
            return x * rhs.x + y * rhs.y + z * rhs.z;
        }

        double length() const {
            return sqrt(x * x + y * y + z * z);
        }

        double distance(const Point &p) const {
            return sqrt((x - p.x) * (x - p.x) + (y - p.y) * (y - p.y) + (z - p.z) * (z - p.z));
        }

        friend std::ostream &operator<<(std::ostream &os, const Point &pt) {
            os << std::fixed << std::setprecision(18) << "(" << pt.x << ", " << pt.y << " " << pt.z << ")";
            return os;
        }
    };

    double dot_prod(const Point &p1, const Point &p2);

    double dot_prod(const std::shared_ptr<Point> &p1, const std::shared_ptr<Point> &p2);

    Point cross_prod(const Point &p1, const Point &p2);

    std::shared_ptr<Point> cross_prod(const std::shared_ptr<Point> &p1, const std::shared_ptr<Point> &p2);

    double getLength(const Point &p1, const Point &p2);

    double getLength(const std::shared_ptr<Point> &p1, const std::shared_ptr<Point> &p2);

    double getArea(const Point &p1, const Point &p2, const Point &p3);

    double
    getArea(const std::shared_ptr<Point> &p1, const std::shared_ptr<Point> &p2, const std::shared_ptr<Point> &p3);

    double getVolume(const Point &p1, const Point &p2, const Point &p3, const Point &p4);

    double
    getVolume(const std::shared_ptr<Point> &p1, const std::shared_ptr<Point> &p2, const std::shared_ptr<Point> &p3,
              const std::shared_ptr<Point> &p4);

    template<typename A>
    using vector2D = std::vector<std::vector<A>>;

    template<typename A>
    using vector3D = std::vector<vector2D<A>>;

    template<typename A>
    using vector4D = std::vector<vector3D<A>>;

    struct Segment {
        int index[2]{};
        int entity_index;

        explicit Segment(int a = 0, int b = 0, int e = 0) {
            index[0] = a;
            index[1] = b;
            entity_index = e;
        }
    };

    struct Triangle {
        int index[3]{};
        int entity_index;

        explicit Triangle(int a = 0, int b = 0, int c = 0, int e = 0) {
            index[0] = a;
            index[1] = b;
            index[2] = c;
            entity_index = e;
        }
    };

    struct Tetrahedron {
        int index[4]{};

        explicit Tetrahedron(int a = 0, int b = 0, int c = 0, int d = 0) {
            index[0] = a;
            index[1] = b;
            index[2] = c;
            index[3] = d;
        }
    };


    struct BoundaryCondition {
        int index;
        int type;
        double temperature;

        explicit BoundaryCondition(int index = 0, int type = 0, int temp = 0) : index(index), type(type),
                                                                                temperature(temp) {}

        friend std::ostream &operator<<(std::ostream &os, const BoundaryCondition &bc) {
            os << "Index: " << bc.index
               << ", Condition: " << bc.type
               << ", Temperature: " << bc.temperature;
            return os;
        }
    };

    struct Band {
        double group_velocity;
        double relaxation_time;
        double Ctot;
        double Lr;

        explicit Band(double vg = 0, double rt = 0, double ct = 0, double Lr = 0) : group_velocity(vg),
                                                                                    relaxation_time(rt), Ctot(ct),
                                                                                    Lr(Lr) {}

        friend std::ostream &operator<<(std::ostream &os, const Band &band) {
            os << "Group velocity: " << band.group_velocity
               << ", relaxation time: " << band.relaxation_time
               << ", Ctot: " << band.Ctot
               << ", Lr: " << band.Lr;
            return os;
        }
    };

    class ContinuousArray {
    public:
        double *data{};
        int dim1{};
        int dim2{};

        ContinuousArray() = default;

        ContinuousArray(int dim1, int dim2) {
            init(dim1, dim2);
        }

        ContinuousArray(const ContinuousArray &rhs) {
            this->dim1 = rhs.dim1;
            this->dim2 = rhs.dim2;
            data = new double[this->dim1 * this->dim2];
            std::copy(rhs.data, rhs.data + rhs.dim1 * rhs.dim2, data);
        }

        ContinuousArray &operator=(const ContinuousArray &rhs) {
            if (this != &rhs) {
                delete[] this->data;
                this->dim1 = rhs.dim1;
                this->dim2 = rhs.dim2;
                data = new double[this->dim1 * this->dim2];
                std::copy(rhs.data, rhs.data + rhs.dim1 * rhs.dim2, data);
            }
            return *this;
        }

        void init(int d1, int d2) {
            delete[] data;
            data = new double[d1 * d2]();
            this->dim1 = d1;
            this->dim2 = d2;
        }

        double *get_ptr(int i, int j) const {
            return data + (i * this->dim2 + j);
        }

        double get(int i, int j) const {
            return *get_ptr(i, j);
        }

        void clear() const {
            std::fill_n(this->data, this->dim1 * this->dim2, 0.0);
        }

        void print() const {
            for (int i = 0; i < dim1; i++) {
                for (int j = 0; j < dim2; j++) {
                    std::cout << data[i * dim2 + j] << "  ";
                }
                std::cout << std::endl;
            }
        }

        ~ContinuousArray() {
            delete[] data;
        }
    };

    std::vector<double> GaussIntegrationPoints(double a, double b, int N);

    size_t get_host_memory();
}
#endif //BTESOLVER_UTILITY_H
