////////////////////////////////////////////////////////////////////////////////
#include <algorithm>
#include <complex>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>

#include <Eigen/Dense>
// Shortcut to avoid  everywhere, DO NOT USE IN .h
using namespace Eigen;
////////////////////////////////////////////////////////////////////////////////

const std::string root_path = DATA_DIR;

// Computes the determinant of the matrix whose columns are the vector u and v
double inline det(const Vector2d &u, const Vector2d &v)
{
    // Well known formula
    return (u.x() * v.y()) - (u.y() * v.x());
}

// Return true iff [a,b] intersects [c,d]
bool intersect_segment(const Vector2d &a, const Vector2d &b, const Vector2d &c, const Vector2d &d)
{
    // parametric coordinates of the intersection
    double t = det(a - c, c - d) / det(a - b, c - d);

    double u = det(a - b, a - c) / det(a - b, c - d);

    // t, u in [0, 1] for both segments, so they intersect if they are both in [0, 1]
    return (0.0 <= t) && (t <= 1.0) && (0.0 <= u) && (u <= 1.0);
}

////////////////////////////////////////////////////////////////////////////////

bool is_inside(const std::vector<Vector2d> &poly, const Vector2d &query)
{
    // 1. Compute bounding box and set coordinate of a point outside the polygon
    Vector2d top_left(0, 0);
    Vector2d bot_right(0, 0);

    for (const auto &point : poly)
    {
        if (point.x() < top_left.x())
            top_left.x() = point.x();
        else if (point.x() > bot_right.x())
            bot_right.x() = point.x();

        if (point.y() > top_left.y())
            top_left.y() = point.y();
        else if (point.y() < bot_right.y())
            bot_right.y() = point.y();
    }
    

    Vector2d outside(top_left.x() - 1, top_left.y() + 1);

    // 2. Cast a ray from the query point to the 'outside' point, count number of intersections
    bool inside = false;
    for (int i = 0; i < poly.size(); i++)
    {
        if (intersect_segment(outside, query, poly[i], poly[(i + 1) % poly.size()]))
            inside = !inside;
    }

    return inside;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Vector2d> load_xyz(const std::string &filename)
{
    std::vector<Vector2d> points;
    std::ifstream in(filename);
    // TODO
    return points;
}

void save_xyz(const std::string &filename, const std::vector<Vector2d> &points)
{
    // TODO
}

std::vector<Vector2d> load_obj(const std::string &filename)
{
    std::ifstream in(filename);
    std::vector<Vector2d> points;
    std::vector<Vector2d> poly;
    char key;
    while (in >> key)
    {
        if (key == 'v')
        {
            double x, y, z;
            in >> x >> y >> z;
            points.push_back(Vector2d(x, y));
        }
        else if (key == 'f')
        {
            std::string line;
            std::getline(in, line);
            std::istringstream ss(line);
            int id;
            while (ss >> id)
            {
                poly.push_back(points[id - 1]);
            }
        }
    }
    return poly;
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
    const std::string points_path = root_path + "/points.xyz";
    const std::string poly_path = root_path + "/polygon.obj";

    std::vector<Vector2d> points = load_xyz(points_path);

    ////////////////////////////////////////////////////////////////////////////////
    //Point in polygon
    std::vector<Vector2d> poly = load_obj(poly_path);
    std::vector<Vector2d> result;
    for (size_t i = 0; i < points.size(); ++i)
    {
        if (is_inside(poly, points[i]))
        {
            result.push_back(points[i]);
        }
    }
    save_xyz("output.xyz", result);

    return 0;
}
