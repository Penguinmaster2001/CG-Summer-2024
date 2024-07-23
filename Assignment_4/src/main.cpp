
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
const double focal_length = 2;
const double field_of_view = 0.7854; //45 degrees
const bool is_perspective = true;
const Vector3d camera_position(0, 0, 2);


// Triangle Mesh
MatrixXd vertices; // n x 3 matrix (n points)
MatrixXi facets;   // m x 3 matrix (m triangles)
AABBTree bvh;
bool use_bvh_tree = true;


//Material for the object, same material for all objects
const Vector4d obj_ambient_color(0.0, 0.5, 0.0, 0);
const Vector4d obj_diffuse_color(0.5, 0.5, 0.5, 0);
const Vector4d obj_specular_color(0.2, 0.2, 0.2, 0);
const double obj_specular_exponent = 256.0;
const Vector4d obj_reflection_color(0.7, 0.7, 0.7, 0);


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

        return intersect;
    }

    // TODO
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

    return intersect;
}




////////////////////////////////////////////////////////////////////////////////
// Raytracer code
////////////////////////////////////////////////////////////////////////////////


Vector4d shoot_ray(const Vector3d &ray_origin, const Vector3d &ray_direction)
{
    //Intersection point and normal, these are output of find_nearest_object
    Vector3d p, N;

    const bool nearest_object = find_nearest_object(ray_origin, ray_direction, p, N);

    if (!nearest_object)
    {
        // Return a transparent color
        return Vector4d(0, 0, 0, 0);
    }

    // Ambient light contribution
    const Vector4d ambient_color = obj_ambient_color.array() * ambient_light.array();

    // Punctual lights contribution (direct lighting)
    Vector4d lights_color(0, 0, 0, 0);
    for (int i = 0; i < light_positions.size(); ++i)
    {
        const Vector3d &light_position = light_positions[i];
        const Vector4d &light_color = light_colors[i];

        Vector4d diff_color = obj_diffuse_color;

        // Diffuse contribution
        const Vector3d Li = (light_position - p).normalized();
        const Vector4d diffuse = diff_color * std::max(Li.dot(N), 0.0);

        // Specular contribution
        const Vector3d Hi = (Li - ray_direction).normalized();
        const Vector4d specular = obj_specular_color * std::pow(std::max(N.dot(Hi), 0.0), obj_specular_exponent);
        // Vector3d specular(0, 0, 0);

        // Attenuate lights according to the squared distance to the lights
        const Vector3d D = light_position - p;
        lights_color += (diffuse + specular).cwiseProduct(light_color) / D.squaredNorm();
    }

    // Rendering equation
    Vector4d C = ambient_color + lights_color;

    //Set alpha to 1
    C(3) = 1;

    return C;
}




////////////////////////////////////////////////////////////////////////////////


void raytrace_scene()
{
    std::cout << "Simple ray tracer." << std::endl;

    int w = 640;
    int h = 480;
    MatrixXd R = MatrixXd::Zero(w, h);
    MatrixXd G = MatrixXd::Zero(w, h);
    MatrixXd B = MatrixXd::Zero(w, h);
    MatrixXd A = MatrixXd::Zero(w, h); // Store the alpha mask

    // The camera always points in the direction -z
    // The sensor grid is at a distance 'focal_length' from the camera center,
    // and covers an viewing angle given by 'field_of_view'.
    double aspect_ratio = double(w) / double(h);
    double image_y = tan(field_of_view / 2.0) * focal_length; // Compute the correct pixels size
    double image_x = image_y * aspect_ratio; // Compute the correct pixels size

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

            const Vector4d C = shoot_ray(ray_origin, ray_direction);
            R(i, j) = C(0);
            G(i, j) = C(1);
            B(i, j) = C(2);
            A(i, j) = C(3);
        }
    }

    // Save to png
    write_matrix_to_png(R, G, B, A, filename);
}




////////////////////////////////////////////////////////////////////////////////


int main(int argc, char *argv[])
{
    setup_scene();

    raytrace_scene();
    return 0;
}
