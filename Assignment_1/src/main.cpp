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
    return (u.x() * v.y()) - (v.x() * u.y());
}

// Return true iff [a,b] intersects [c,d]
bool intersect_segment(const Vector2d &a, const Vector2d &b, const Vector2d &c, const Vector2d &d)
{
    // NOTE: Colinear lines are assumed to not intersect
    // However, in practice it is extremely rare for lines to be both colinear and overlapping

    // Check first line between second line ends
    double d1 = det(b - a, c - a);
    double d2 = det(b - a, d - a);

    if (d1 * d2 >= 0.0)
        return false;

    // Check second line between first line ends
    double p1 = det(d - c, a - c);
    double p2 = det(d - c, b - c);

    return (p1 * p2 < 0.0);
}

////////////////////////////////////////////////////////////////////////////////

bool is_inside(const std::vector<Vector2d> &poly, const Vector2d &query)
{
    // 1. Compute bounding box and set coordinate of a point outside the polygon
    Vector2d outside(0, 0);

    // Get the max x and y
    for (const auto &point : poly)
    {
        if (point.x() > outside.x())
            outside.x() = point.x();

        if (point.y() > outside.y())
            outside.y() = point.y();
    }

    // Ensure it's outside by moving away by 1 in each direction
    outside += Vector2d(1, 1);

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

    /* WARNING: Very little error handling
     * Does NOT check for valid numerical strings
     * Does NOT check for correct file length
     * I am VERY lazy and I WELCOME security issues with OPEN arms
     */
    if (!in.is_open())
    {
        std::cout << "Unable to open " << filename << std::endl;
        return points;
    }

    std::string line;

    // First line is num of points
    std::getline(in, line);
    int num_points = std::stoi(line);

    // Reduces the amount of resizing
    points.reserve(num_points);

    // Read until all points read or no more lines (hopefully not)
    for (int i = 0; i < num_points && std::getline(in, line); i++)
    {
        std::istringstream iss(line);
        double x, y;

        // Get the doubles from the line, make a new Vector2d, and put it into points
        if (iss >> x >> y)
        {
            Vector2d new_vector(x, y);
            points.push_back(new_vector);
        }
    }

    in.close();

    return points;
}

void save_xyz(const std::string &filename, const std::vector<Vector2d> &points)
{
    std::ofstream out(filename);

    if (!out.is_open())
    {
        std::cout << "Unable to open " << filename << std::endl;
        return;
    }

    // Write the number of points 
    out << points.size() << "\n";

    // Loop through every point and write it to the file
    for (const auto &point : points)
    {
        double x = point.x();
        double y = point.y();

        out << x << " " << y << " 0\n";
    }

    // Add an extra line
    out << std::endl;

    out.close();
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
