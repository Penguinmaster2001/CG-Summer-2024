
////////////////////////////////////////////////////////////////////////////////


// C++ include
#include <iostream>
#include <string>
#include <vector>
#include <queue>
#include <limits>
#include <fstream>
#include <algorithm>
#include <numeric>


// I didn't use __DBL_EPSILON__ or limits because they were too small and caused issues
#define EPSILON 1.0E-12

// Utilities for the Assignment
#include "utils.h"


// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"


// Shortcut to avoid Eigen:: everywhere, DO NOT USE IN .h
using namespace Eigen;




////////////////////////////////////////////////////////////////////////////////
// Class to store tree
////////////////////////////////////////////////////////////////////////////////


class AABBTree
{
public:
    class Node
    {
    public:
        AlignedBox3d bbox;
        int parent;   // Index of the parent node (-1 for root)
        int left;     // Index of the left child (-1 for a leaf)
        int right;    // Index of the right child (-1 for a leaf)
        int triangle; // Index of the node triangle (-1 for internal nodes)
    };

    std::vector<Node> nodes;
    int root;

    AABBTree() = default;                           // Default empty constructor
    AABBTree(const MatrixXd &V, const MatrixXi &F); // Build a BVH from an existing mesh

private:
    // builds the bvh recursively
    int build_recursive(const MatrixXd &V, const MatrixXi &F, const MatrixXd &centroids,
        int from, int to, int parent, std::vector<int> &triangles);
};




////////////////////////////////////////////////////////////////////////////////
// Scene setup, global variables
////////////////////////////////////////////////////////////////////////////////


const std::string data_dir = DATA_DIR;
const std::string filename("raytrace.png");
const std::string mesh_filename(data_dir + "dragon.off");


//Camera settings
const double focal_length = 5;
const double field_of_view = 0.7854; //45 degrees
const bool is_perspective = true;
const Vector3d camera_position(0, 0, 2);



// Maximum number of recursive calls
// Reflection and refraction use this independently
// So the true maximum is like 2^max_bounce or something idk
const int max_bounce = 4;


// Triangle Mesh
MatrixXd vertices; // n x 3 matrix (n points)
MatrixXi facets;   // m x 3 matrix (m triangles)
AABBTree bvh;
const bool use_bvh_tree = true;

// Objects
std::vector<Vector3d> sphere_centers;
std::vector<double> sphere_radii;
std::vector<Matrix3d> parallelograms;


//Material for the object, same material for all objects
const Vector4d obj_ambient_color(0.0, 0.5, 0.0, 0);
const Vector4d obj_diffuse_color(0.5, 0.5, 0.5, 0);
const Vector4d obj_specular_color(0.2, 0.2, 0.2, 0);
const double obj_specular_exponent = 256.0;
const Vector4d obj_reflection_color(0.7, 0.7, 0.7, 0);
const Vector4d obj_refraction_color(0.7, 0.7, 0.7, 0);
// Diamond
const double obj_refractive_index = 2.42;


// Precomputed (or otherwise) gradient vectors at each grid node
const int grid_size = 20;
std::vector<std::vector<Vector2d>> grid;


//Lights
std::vector<Vector3d> light_positions;
std::vector<Vector4d> light_colors;


//Ambient light
const Vector4d ambient_light(0.2, 0.2, 0.2, 0);



//Fills the different arrays
void setup_scene()
{
    //Loads file
    std::ifstream in(mesh_filename);
    std::string token;
    in >> token;
    int nv, nf, ne;
    in >> nv >> nf >> ne;
    vertices.resize(nv, 3);
    facets.resize(nf, 3);
    for (int i = 0; i < nv; ++i)
    {
        in >> vertices(i, 0) >> vertices(i, 1) >> vertices(i, 2);
    }
    for (int i = 0; i < nf; ++i)
    {
        int s;
        in >> s >> facets(i, 0) >> facets(i, 1) >> facets(i, 2);
        assert(s == 3);
    }

    //setup tree
    bvh = AABBTree(vertices, facets);

    //Lights
    light_positions.emplace_back(8, 8, 0);
    light_colors.emplace_back(16, 16, 16, 0);

    light_positions.emplace_back(6, -8, 0);
    light_colors.emplace_back(16, 16, 16, 0);

    light_positions.emplace_back(4, 8, 0);
    light_colors.emplace_back(16, 16, 16, 0);

    light_positions.emplace_back(2, -8, 0);
    light_colors.emplace_back(16, 16, 16, 0);

    light_positions.emplace_back(0, 8, 0);
    light_colors.emplace_back(16, 16, 16, 0);

    light_positions.emplace_back(-2, -8, 0);
    light_colors.emplace_back(16, 16, 16, 0);

    light_positions.emplace_back(-4, 8, 0);
    light_colors.emplace_back(16, 16, 16, 0);

    // I made this one bigger because there's this really random and weird bug that causes
    // y1 to be 21 in perlin() and I don't want to figure out why
    grid.resize(grid_size + 2);
    for (int i = 0; i <= grid_size + 1; ++i)
    {
        grid[i].resize(grid_size + 2);
        for (int j = 0; j <= grid_size + 1; ++j)
            grid[i][j] = Vector2d::Random().normalized();
    }

    //Spheres
    sphere_centers.emplace_back(10, 0, 1);
    sphere_radii.emplace_back(1);

    sphere_centers.emplace_back(7, 0.05, -1);
    sphere_radii.emplace_back(1);

    sphere_centers.emplace_back(4, 0.1, 1);
    sphere_radii.emplace_back(1);

    sphere_centers.emplace_back(1, 0.2, -1);
    sphere_radii.emplace_back(1);

    sphere_centers.emplace_back(-2, 0.4, 1);
    sphere_radii.emplace_back(1);

    sphere_centers.emplace_back(-5, 0.8, -1);
    sphere_radii.emplace_back(1);

    sphere_centers.emplace_back(-8, 1.6, 1);
    sphere_radii.emplace_back(1);

    //parallelograms
    parallelograms.emplace_back();
    parallelograms.back() << -100, 100, -100,
        -1.25, 0, -1.2,
        -100, -100, 100;
}





////////////////////////////////////////////////////////////////////////////////
// Perlin noise code
////////////////////////////////////////////////////////////////////////////////

// Function to linearly interpolate between a0 and a1
// Weight w should be in the range [0.0, 1.0]
double lerp(double a0, double a1, double w)
{
    assert(w >= 0);
    assert(w <= 1);
    
    // Linear
    // return a0 + w * (a1 - a0);

    // Cubic
    return (a1 - a0) * (3.0 - w * 2.0) * w * w + a0;
}



// Computes the dot product of the distance and gradient vectors.
double dotGridGradient(int ix, int iy, double x, double y)
{
    // Compute the distance vector
    float dx = x - (float) ix;
    float dy = y - (float) iy;

    // Compute and return the dot-product
    return (dx * grid[iy][ix][0]) + (dy * grid[iy][ix][1]);
}



// Compute Perlin noise at coordinates x, y
double perlin(double x, double y)
{
    // Determine grid cell coordinates x0, y0
    int x0 = int(x);
    int x1 = x0 + 1;
    int y0 = int(y);
    int y1 = y0 + 1;

    // Determine interpolation weights
    double sx = x - x0;
    double sy = y - y0;

    // Interpolate between grid point gradients
    double n0 = dotGridGradient(x0, y0, x, y);
    double n1 = dotGridGradient(x1, y0, x, y);

    double ix0 = lerp(n0, n1, sx);

    n0 = dotGridGradient(x0, y1, x, y);
    n1 = dotGridGradient(x1, y1, x, y);

    double ix1 = lerp(n0, n1, sx);
    double value = lerp(ix0, ix1, sy);

    return value;
}



Vector4d procedural_texture(const double tu, const double tv)
{
    assert(tu >= 0);
    assert(tv >= 0);

    assert(tu <= 1);
    assert(tv <= 1);

    // Uncomment these lines once you implement the perlin noise
    const double color = (perlin(tu * grid_size, tv * grid_size) + 1) / 2;

    return Vector4d(0, color, 0, 0);

    // Example of checkerboard texture
    // const double color = (int(tu * grid_size) + int(tv * grid_size)) % 2 == 0 ? 0 : 1;
    // return Vector4d(0, color, 0, 0);
}





////////////////////////////////////////////////////////////////////////////////
// BVH Code
////////////////////////////////////////////////////////////////////////////////


AlignedBox3d bbox_from_triangle(const Vector3d &a, const Vector3d &b, const Vector3d &c)
{
    AlignedBox3d box;
    box.extend(a);
    box.extend(b);
    box.extend(c);
    return box;
}



AABBTree::AABBTree(const MatrixXd &V, const MatrixXi &F)
{
    // Compute the centroids of all the triangles in the input mesh
    MatrixXd centroids(F.rows(), V.cols());
    centroids.setZero();
    for (int i = 0; i < F.rows(); ++i)
    {
        for (int k = 0; k < F.cols(); ++k)
        {
            centroids.row(i) += V.row(F(i, k));
        }

        centroids.row(i) /= F.cols();
    }

    //Vector containing the list of triangle indices
    std::vector<int> triangles(F.rows());

    std::iota(triangles.begin(), triangles.end(), 0);

    root = build_recursive(V, F, centroids, 0, triangles.size(), -1, triangles);
}



// V are the vertices n x 3
// F are the facets m x 3
// centroids are the centers of the facets m x 3
// from is the index of the first tri possibly in this box
// to is the index after the last tri possibly in this box
// triangles is the order of the triangles along the largest axis
int AABBTree::build_recursive(const MatrixXd &V, const MatrixXi &F, const MatrixXd &centroids,
    int from, int to, int parent, std::vector<int> &triangles)
{
    // Scene is empty, so is the aabb tree
    if (to - from <= 0)
    {
        return -1;
    }


    AABBTree::Node new_node;

    int index = nodes.size();


    // If there is only 1 triangle left, then we are at a leaf
    if (to - from == 1)
    {
        int vert_ind0 = F(triangles[from], 0);
        int vert_ind1 = F(triangles[from], 1);
        int vert_ind2 = F(triangles[from], 2);

        new_node.bbox = bbox_from_triangle(V.row(vert_ind0), V.row(vert_ind1), V.row(vert_ind2));
        new_node.parent = parent;
        new_node.left = -1;
        new_node.right = -1;
        new_node.triangle = triangles[from];

        nodes.push_back(new_node);

        return index;
    }


    AlignedBox3d centroid_box(3);

    // Use AlignedBox3d to find the box around the current centroids
    for (int i = from; i < to; ++i)
    {
        int vert_ind0 = F(triangles[i], 0);
        int vert_ind1 = F(triangles[i], 1);
        int vert_ind2 = F(triangles[i], 2);

        centroid_box.extend(bbox_from_triangle(V.row(vert_ind0), V.row(vert_ind1), V.row(vert_ind2)));
    }


    // Diagonal of the box
    Vector3d extent = centroid_box.diagonal();


    // Find the largest dimension
    int longest_dim = extent.x() > extent.y() ? (extent.x() > extent.z() ? 0 : 2) : (extent.y() > extent.z() ? 1 : 2);



    // Sort centroids along the longest dimension
    std::sort(triangles.begin() + from, triangles.begin() + to, [&](int f1, int f2)
        {
            // return true if triangle f1 comes before triangle f2
            return centroids(f1, longest_dim) < centroids(f2, longest_dim);
        });


    int split = (from + to) / 2;

    // Create a new internal node and do a recursive call to build the left and right part of the tree

    // Trust me there's a reason I did it this way
    nodes.push_back(new_node);

    new_node.bbox = centroid_box;
    new_node.parent = parent;
    new_node.left = build_recursive(V, F, centroids, from, split, index, triangles);
    new_node.right = build_recursive(V, F, centroids, split, to, index, triangles);
    new_node.triangle = -1;

    // Yeah this is bad
    nodes[index] = new_node;


    // Finally return the correct index
    return index;
}




////////////////////////////////////////////////////////////////////////////////
// Intersection code
////////////////////////////////////////////////////////////////////////////////


//Compute the intersection between a ray and a sphere, return -1 if no intersection
double ray_sphere_intersection(const Vector3d &ray_origin, const Vector3d &ray_direction, int index, Vector3d &p, Vector3d &N)
{
    //return t or -1 if no intersection

    // From Assignment 2

    const Vector3d sphere_center = sphere_centers[index];
    const double sphere_radius = sphere_radii[index];

    // Moving the origin
    Vector3d c = ray_origin - sphere_center;

    double ray_dir_dot_c = ray_direction.dot(c);
    double ray_dir_sqr = ray_direction.dot(ray_direction);

    double discriminant = (ray_dir_dot_c * ray_dir_dot_c) - (ray_dir_sqr * (c.dot(c) - (sphere_radius * sphere_radius)));
    if (discriminant < 0.0) return -1;

    double t = (-ray_dir_dot_c - std::sqrt(discriminant)) / ray_dir_sqr;
    if (t < 0.0) return -1;

    p = ray_origin + (t * ray_direction);
    N = (p - sphere_center) / sphere_radius;

    return t;
}


double det_of_3_col_vectors(Vector3d a, Vector3d b, Vector3d c)
{
    return a.x() * (b.y() * c.z() - c.y() * b.z())
         + b.x() * (c.y() * a.z() - a.y() * c.z())
         + c.x() * (a.y() * b.z() - b.y() * a.z());
}


double ray_triangle_intersection(const Vector3d &ray_origin, const Vector3d &ray_direction,
    const Vector3d &a, const Vector3d &b, const Vector3d &c, Vector3d &p, Vector3d &N)
{
    // We need to solve a linear system of 3 equations
    // The equations are parametrization of lines
    // Eqn 1 is the ray in terms of t
    // Eqn 2 is the first edge of the triangle in terms of s
    // Eqn 3 is the second edge of the triangle in terms of r

    // The solution to these are the parametrization variables
    // The ray intersects if t > 0 and 0 < s, r < 1, s + r < 1

    // This acts as a change of coords to set vert_a at the origin
    // The names are awful, I know
    Vector3d u = a - b;
    Vector3d v = a - c;
    Vector3d e = a - ray_origin;
    Vector3d d = ray_direction;

    // This is a common denominator from using Cramer's rule
    double det = det_of_3_col_vectors(u, v, d);

    // Ray is parallel
    if (fabs(det) < __DBL_EPSILON__) return -1;

    double inv_det = 1.0 / det;

    // Magic (Cramer's rule)
    double t = det_of_3_col_vectors(u, v, e) * inv_det;
    if (t < 0.0) return -1;

    double s = det_of_3_col_vectors(u, e, d) * inv_det;
    if (s < 0.0 || 1.0 < s) return -1;

    double r = det_of_3_col_vectors(e, v, d) * inv_det;
    if (r < 0.0 || 1.0 < r + s) return -1;

    p = ray_origin + (t * ray_direction);
    N = v.cross(u).normalized();
    return t;
}



//Compute the intersection between a ray and a parallelogram, return -1 if no intersection
// From Assignment 2
double ray_parallelogram_intersection(const Vector3d &ray_origin, const Vector3d &ray_direction, int index, Vector3d &p, Vector3d &N)
{
    //return t or -1 if no intersection

    const Vector3d pgram_origin = parallelograms[index].col(0);
    const Vector3d A = parallelograms[index].col(1);
    const Vector3d B = parallelograms[index].col(2);
    const Vector3d pgram_u = A - pgram_origin;
    const Vector3d pgram_v = B - pgram_origin;


    // To fix the confusion around the definition of the parallelogram
    const Vector3d o = pgram_origin;
    const Vector3d u = pgram_origin + pgram_u;
    const Vector3d v = pgram_origin + pgram_v;


    // A parallelogram is just two triangles
    // First tri:
    // o -> u -> v
    double t = ray_triangle_intersection(ray_origin, ray_direction, o, u, v, p, N);

    if (t >= 0) return t;

    // Second tri:
    // v -> o -> v + u - o
    t = ray_triangle_intersection(ray_origin, ray_direction, v, u, v + u - o, p, N);

    return t;
}



// Compute whether the ray intersects the given box.
bool ray_box_intersection(const Vector3d &ray_origin, const Vector3d &ray_direction, const AlignedBox3d &box)
{
    // we are not testing with the real surface here anyway

   double tmin = -std::numeric_limits<double>::infinity();
    double tmax = std::numeric_limits<double>::infinity();

    for (int i = 0; i < 3; ++i)
    {
        double invD = 1.0 / ray_direction[i];
        double t0 = (box.min()[i] - ray_origin[i]) * invD;
        double t1 = (box.max()[i] - ray_origin[i]) * invD;

        if (invD < 0.0) std::swap(t0, t1);

        tmin = std::max(tmin, t0);
        tmax = std::min(tmax, t1);

        if (tmax < tmin) return false;
    }

    return true;
}



//Finds the closest intersecting object returns its index
//In case of intersection it writes into p and N (intersection point and normals)
bool find_nearest_object(const Vector3d &ray_origin, const Vector3d &ray_direction, Vector3d &p, Vector3d &N)
{
    Vector3d tmp_p, tmp_N;

    // Method (1): Traverse every triangle and return the closest hit.
    bool intersect = false;
    double closest_t = std::numeric_limits<double>::max(); //closest t is "+ infinity"

    if (!use_bvh_tree)
    {
        // C++ is truly an amazing language in the way it automatically obfuscates itself
        for (int facet_ind = 0; facet_ind < facets.rows(); ++facet_ind)
        {
            int vert_ind0 = facets(facet_ind, 0);
            int vert_ind1 = facets(facet_ind, 1);
            int vert_ind2 = facets(facet_ind, 2);

            //returns t and writes on tmp_p and tmp_N
            const double t = ray_triangle_intersection(ray_origin, ray_direction,
                vertices.row(vert_ind0), vertices.row(vert_ind1), vertices.row(vert_ind2), tmp_p, tmp_N);
            
            //We have intersection
            if (t >= 0)
            {
                //The point is before our current closest t
                if (t < closest_t)
                {
                    intersect = true;
                    closest_t = t;
                    p = tmp_p;
                    N = tmp_N;
                }
            }
        }
    }
    else
    {
        // Method (2): Traverse the BVH tree and test the intersection with a
        // triangles at the leaf nodes that intersects the input ray.

        // Level order traversal

        // The tree is empty
        if (bvh.root < 0) return false;

        std::queue<AABBTree::Node> current_nodes;

        current_nodes.push(bvh.nodes[bvh.root]);
        
        while (!current_nodes.empty())
        {
            AABBTree::Node current_node = current_nodes.front();

            // Found the triangle thingy
            if (current_node.triangle >= 0)
            {
                int vert_ind0 = facets(current_node.triangle, 0);
                int vert_ind1 = facets(current_node.triangle, 1);
                int vert_ind2 = facets(current_node.triangle, 2);

                //returns t and writes on tmp_p and tmp_N
                const double t = ray_triangle_intersection(ray_origin, ray_direction,
                    vertices.row(vert_ind0), vertices.row(vert_ind1), vertices.row(vert_ind2), tmp_p, tmp_N);
                
                //We have intersection
                if (t >= 0)
                {
                    //The point is before our current closest t
                    if (t < closest_t)
                    {
                        intersect = true;
                        closest_t = t;
                        p = tmp_p;
                        N = tmp_N;
                    }
                }
            }

            current_nodes.pop();

            // Enqueue left child
            if (current_node.left >= 0)
            {
                AABBTree::Node left_node = bvh.nodes[current_node.left];

                if (ray_box_intersection(ray_origin, ray_direction, left_node.bbox))
                {
                    current_nodes.push(left_node);
                }
            }

            // Enqueue right child
            if (current_node.right >= 0)
            {
                AABBTree::Node right_node = bvh.nodes[current_node.right];

                if (ray_box_intersection(ray_origin, ray_direction, right_node.bbox))
                {
                    current_nodes.push(right_node);
                }
            }
        }
    }



    for (int i = 0; i < sphere_centers.size(); ++i)
    {
        //returns t and writes on tmp_p and tmp_N
        const double t = ray_sphere_intersection(ray_origin, ray_direction, i, tmp_p, tmp_N);
        //We have intersection
        if (t >= 0)
        {
            //The point is before our current closest t
            if (t < closest_t)
            {
                intersect = true;

                closest_t = t;
                p = tmp_p;
                N = tmp_N;
            }
        }
    }


    for (int i = 0; i < parallelograms.size(); ++i)
    {
        //returns t and writes on tmp_p and tmp_N
        const double t = ray_parallelogram_intersection(ray_origin, ray_direction, i, tmp_p, tmp_N);
        //We have intersection
        if (t >= 0)
        {
            //The point is before our current closest t
            if (t < closest_t)
            {
                intersect = true;

                closest_t = t;
                p = tmp_p;
                N = tmp_N;
            }
        }
    }


    return intersect;
}




////////////////////////////////////////////////////////////////////////////////
// Raytracer code
////////////////////////////////////////////////////////////////////////////////


//Checks if the light is visible
bool is_light_visible(const Vector3d &ray_origin, const Vector3d &ray_direction, const Vector3d &light_position)
{
    // Determine if the light is visible here
    // Use find_nearest_object
    Vector3d p, N;

    // -1 means no object was hit so the light is visible.
    return (find_nearest_object(ray_origin, ray_direction, p, N) < 0);
}



// Returns true and sets the direction if refracting
// false otherwise
bool calc_refraction(const Vector3d &ray_direction, const Vector3d &N, Vector3d &refraction_direction)
{
    double d_dot_n = ray_direction.dot(N);

    double det = 1.0 - ((obj_refractive_index * obj_refractive_index) * (1.0 - (d_dot_n * d_dot_n)));

    // Total internal reflection
    if (det < 0.0) return false;

    refraction_direction = (obj_refractive_index * (ray_direction - (d_dot_n * N)))
        - det * N;

    return true;
}


Vector4d shoot_ray(const Vector3d &ray_origin, const Vector3d &ray_direction, int bounce)
{
    //Intersection point and normal, these are output of find_nearest_object
    Vector3d p, N;

    const bool nearest_object = find_nearest_object(ray_origin, ray_direction, p, N);

    if (!nearest_object)
    {
        // Return a transparent color
        return Vector4d(0, 0, 0, 0);
    }

    N = N.normalized();
    Vector3d offset_p = p + (EPSILON * N);
    Vector3d internal_offset_p = p - (EPSILON * N);

    // Ambient light contribution
    const Vector4d ambient_color = obj_ambient_color.array() * ambient_light.array();

    // Punctual lights contribution (direct lighting)
    Vector4d lights_color(0, 0, 0, 0);
    for (int i = 0; i < light_positions.size(); ++i)
    {
        const Vector3d &light_position = light_positions[i];
        const Vector4d &light_color = light_colors[i];

        // Shadow ray direction
        const Vector3d Li = (light_position - p).normalized();

        // Shoot a shadow ray to determine if the light should affect the intersection point and call is_light_visible
        if (!is_light_visible(offset_p, Li, light_position)) continue;


        Vector4d diff_color = obj_diffuse_color;

        //Compute UV coordinates for the point on the sphere
        // const double x = p(0) - sphere_centers[nearest_object][0];
        // const double y = p(1) - sphere_centers[nearest_object][1];
        // const double z = p(2) - sphere_centers[nearest_object][2];
        // double tu = acos(z / sphere_radii[nearest_object]) / 3.1415;
        // double tv = (3.1415 + atan2(y, x)) / (2 * 3.1415);
        // tu = std::min(tu, 1.0);
        // tu = std::max(tu, 0.0);

        // tv = std::min(tv, 1.0);
        // tv = std::max(tv, 0.0);

        // diff_color = procedural_texture(tu, tv);

        // Diffuse contribution
        const Vector4d diffuse = diff_color * std::max(Li.dot(N), 0.0);

        // Specular contribution
        const Vector3d Hi = (Li - ray_direction).normalized();
        const Vector4d specular = obj_specular_color * std::pow(std::max(N.dot(Hi), 0.0), obj_specular_exponent);
        // Vector3d specular(0, 0, 0);

        // Attenuate lights according to the squared distance to the lights
        const Vector3d D = light_position - p;
        lights_color += (diffuse + specular).cwiseProduct(light_color) / D.squaredNorm();
    }

    Vector4d refl_color = obj_reflection_color;


    Vector3d reflection_direction = ray_direction - ((2.0 * ray_direction.dot(N)) * N);

    Vector3d refraction_direction;

    // Compute the color of the reflected ray and add its contribution to the current point color.
    // use refl_color
    Vector4d reflection_color = bounce > 0
        ? refl_color.cwiseProduct(shoot_ray(offset_p, reflection_direction, bounce - 1))
        : Vector4d(0, 0, 0, 0);

    // Compute the color of the refracted ray and add its contribution to the current point color.
    // Make sure to check for total internal reflection before shooting a new ray.
    Vector4d refraction_color = (bounce > 0 && calc_refraction(ray_direction, N, refraction_direction))
        ? obj_refraction_color.cwiseProduct(shoot_ray(internal_offset_p, refraction_direction, bounce - 1))
        : Vector4d(0, 0, 0, 0);

    // Rendering equation
    Vector4d C = ambient_color + lights_color + reflection_color + refraction_color;

    //Set alpha to 1
    C(3) = 1;

    return C;
}




////////////////////////////////////////////////////////////////////////////////


void raytrace_scene()
{
    std::cout << "Simple ray tracer." << std::endl;

    const int w = 640;
    const int h = 480;
    MatrixXd R = MatrixXd::Zero(w, h);
    MatrixXd G = MatrixXd::Zero(w, h);
    MatrixXd B = MatrixXd::Zero(w, h);
    MatrixXd A = MatrixXd::Zero(w, h); // Store the alpha mask

    // The camera always points in the direction -z
    // The sensor grid is at a distance 'focal_length' from the camera center,
    // and covers an viewing angle given by 'field_of_view'.
    const double aspect_ratio = double(w) / double(h);
    const double image_y = tan(field_of_view / 2.0) * focal_length; // Compute the correct pixels size
    const double image_x = image_y * aspect_ratio; // Compute the correct pixels size

    // The pixel grid through which we shoot rays is at a distance 'focal_length'
    const Vector3d image_origin(-image_x, image_y, camera_position[2] - focal_length);
    const Vector3d x_displacement(2.0 / w * image_x, 0, 0);
    const Vector3d y_displacement(0, -2.0 / h * image_y, 0);

    for (unsigned i = 0; i < w; ++i)
    {
        for (unsigned j = 0; j < h; ++j)
        {
            const Vector3d pixel_center = image_origin + (i + 0.5) * x_displacement + (j + 0.5) * y_displacement;

            // Prepare the ray
            Vector3d ray_origin;
            Vector3d ray_direction;

            if (is_perspective)
            {
                // Perspective camera
                ray_origin = camera_position;
                ray_direction = (pixel_center - camera_position).normalized();
            }
            else
            {
                // Orthographic camera
                ray_origin = pixel_center;
                ray_direction = Vector3d(0, 0, -1);
            }

            const Vector4d C = shoot_ray(ray_origin, ray_direction, max_bounce);
            R(i, j) = C(0);
            G(i, j) = C(1);
            B(i, j) = C(2);
            A(i, j) = C(3);
        }

        std::cout << (100.0 * double(i) / double(w)) << "\% done (" << i << " / " << w << " rows)\n";
    }

    // Save to png
    write_matrix_to_png(R, G, B, A, filename);

    // For completeness
    std::cout << "100\% done!" << std::endl;
}




////////////////////////////////////////////////////////////////////////////////


int main(int argc, char *argv[])
{
    setup_scene();

    raytrace_scene();
    return 0;
}
