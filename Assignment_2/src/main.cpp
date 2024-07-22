// C++ include
#include <iostream>
#include <string>
#include <vector>

// Utilities for the Assignment
#include "utils.h"

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"

// Shortcut to avoid Eigen:: everywhere, DO NOT USE IN .h
using namespace Eigen;

void raytrace_sphere()
{
    std::cout << "Simple ray tracer, one sphere with orthographic projection" << std::endl;

    const std::string filename("sphere_orthographic.png");
    MatrixXd C = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is orthographic, pointing in the direction -z and covering the
    // unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / C.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / C.rows(), 0);

    // Single light source
    const Vector3d light_position(-1, 1, 1);

    for (unsigned i = 0; i < C.cols(); ++i)
    {
        for (unsigned j = 0; j < C.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            // Prepare the ray
            const Vector3d ray_origin = pixel_center;
            const Vector3d ray_direction = camera_view_direction;

            // Intersect with the sphere
            // NOTE: this is a special case of a sphere centered in the origin and for orthographic rays aligned with the z axis
            Vector2d ray_on_xy(ray_origin(0), ray_origin(1));
            const double sphere_radius = 0.9;

            if (ray_on_xy.norm() < sphere_radius)
            {
                // The ray hit the sphere, compute the exact intersection point
                Vector3d ray_intersection(
                    ray_on_xy(0), ray_on_xy(1),
                    sqrt(sphere_radius * sphere_radius - ray_on_xy.squaredNorm()));

                // Compute normal at the intersection point
                Vector3d ray_normal = ray_intersection.normalized();

                // Simple diffuse model
                C(i, j) = (light_position - ray_intersection).normalized().transpose() * ray_normal;

                // Clamp to zero
                C(i, j) = std::max(C(i, j), 0.);

                // Disable the alpha mask for this pixel
                A(i, j) = 1;
            }
        }
    }

    // Save to png
    write_matrix_to_png(C, C, C, A, filename);
}



double det_of_3_col_vectors(Vector3d a, Vector3d b,Vector3d c)
{
    return a.x() * (b.y() * c.z() - c.y() * b.z())
         + b.x() * (c.y() * a.z() - a.y() * c.z())
         + c.x() * (a.y() * b.z() - b.y() * a.z());
}



// Returns if the ray with direction and origin intersects the tri defined by 3 verts
// False if the ray is parallel
// The intersection is stored in intersection_point if it exists
bool intersect_triangle(Vector3d ray_orig, Vector3d ray_dir, Vector3d vert_a, Vector3d vert_b, Vector3d vert_c,
    std::shared_ptr<Vector3d> intersection_point, std::shared_ptr<Vector3d> normal)
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
    Vector3d u = vert_a - vert_b;
    Vector3d v = vert_a - vert_c;
    Vector3d e = vert_a - ray_orig;
    Vector3d d = ray_dir;

    // This is a common denominator from using Cramer's rule
    double c = det_of_3_col_vectors(u, v, d);

    // Keep everything positive to reduce branching
    // I have no idea if this is actually faster
    // however it is more readable if we use a generous definition of readable
    // int sign = c > 0.0 ? 1 : -1;
    // c *= sign;
    
    // Ray is parallel
    // if (c < __DBL_EPSILON__) return false;

    // Magic (Cramer's rule)
    double t = det_of_3_col_vectors(u, v, e) / c;
    if (t < 0.0) return false;

    double s = det_of_3_col_vectors(u, e, d) / c;
    if (s < 0.0 || 1.0 < s) return false;

    double r = det_of_3_col_vectors(e, v, d) / c;
    if (r < 0.0 || 1.0 < r + s) return false;

    *intersection_point = ray_orig + (t * ray_dir);
    *normal = v.cross(u).normalized();
    return true;
}



void raytrace_parallelogram()
{
    std::cout << "Simple ray tracer, one parallelogram with orthographic projection" << std::endl;

    const std::string filename("plane_orthographic.png");
    MatrixXd C = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / C.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / C.rows(), 0);

    // Parameters of the parallelogram (position of the lower-left corner + two sides)
    const Vector3d pgram_origin(-0.5, -0.5, 0);
    const Vector3d pgram_u(0, 0.7, -10);
    const Vector3d pgram_v(1, 0.4, 0);

    // Single light source
    const Vector3d light_position(-1, 1, 1);

    for (unsigned i = 0; i < C.cols(); ++i)
    {
        for (unsigned j = 0; j < C.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            // Prepare the ray
            const Vector3d ray_origin = pixel_center;
            const Vector3d ray_direction = camera_view_direction;

            // Check if the ray intersects with the parallelogram
            std::shared_ptr<Vector3d> ray_intersection = std::make_shared<Vector3d>(0, 0, 0);
            std::shared_ptr<Vector3d> ray_normal = std::make_shared<Vector3d>(0, 0, 0);

            // To fix the confusion around the definition of the parallelogram
            const Vector3d o = pgram_origin;
            const Vector3d u = pgram_origin + pgram_u;
            const Vector3d v = pgram_origin + pgram_v;

            // A parallelogram is just two triangles
            // First tri:
            // o -> u -> v
            if (intersect_triangle(ray_origin, ray_direction, o, u, v, ray_intersection, ray_normal)
            
            // Second tri:
            // v -> o -> v + u - o
                || intersect_triangle(ray_origin, ray_direction, v, u, v + u - o, ray_intersection, ray_normal))
            {
                // Simple diffuse model
                C(i, j) = (light_position - (*ray_intersection)).normalized().dot(*ray_normal);

                // Clamp to zero
                C(i, j) = std::max(C(i, j), 0.);

                // Disable the alpha mask for this pixel
                A(i, j) = 1;
            }
        }
    }

    // Save to png
    write_matrix_to_png(C, C, C, A, filename);
}



void raytrace_perspective()
{
    std::cout << "Simple ray tracer, one parallelogram with perspective projection" << std::endl;

    const std::string filename("plane_perspective.png");
    MatrixXd C = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is perspective, pointing in the direction -z and covering the unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / C.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / C.rows(), 0);

    // TODO: Parameters of the parallelogram (position of the lower-left corner + two sides)
    const Vector3d pgram_origin(-0.5, -0.5, 0);
    const Vector3d pgram_u(0, 0.7, -10);
    const Vector3d pgram_v(1, 0.4, 0);

    // Single light source
    const Vector3d light_position(-1, 1, 1);

    for (unsigned i = 0; i < C.cols(); ++i)
    {
        for (unsigned j = 0; j < C.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            // Prepare the ray (origin point and direction)
            const Vector3d ray_origin = pixel_center;
            const Vector3d ray_direction = pixel_center - camera_origin;

            // Check if the ray intersects the parallelogram
            std::shared_ptr<Vector3d> ray_intersection = std::make_shared<Vector3d>(0, 0, 0);
            std::shared_ptr<Vector3d> ray_normal = std::make_shared<Vector3d>(0, 0, 0);

            // To fix the confusion around the definition of the parallelogram
            const Vector3d o = pgram_origin;
            const Vector3d u = pgram_origin + pgram_u;
            const Vector3d v = pgram_origin + pgram_v;

            // A parallelogram is just two triangles
            // First tri:
            // o -> u -> v
            if (intersect_triangle(ray_origin, ray_direction, o, u, v, ray_intersection, ray_normal)
            
            // Second tri:
            // v -> o -> v + u - o
                || intersect_triangle(ray_origin, ray_direction, v, u, v + u - o, ray_intersection, ray_normal))
            {
                // Simple diffuse model
                C(i, j) = (light_position - (*ray_intersection)).normalized().dot(*ray_normal);

                // Clamp to zero
                C(i, j) = std::max(C(i, j), 0.);

                // Disable the alpha mask for this pixel
                A(i, j) = 1;
            }
        }
    }

    // Save to png
    write_matrix_to_png(C, C, C, A, filename);
}




bool intersect_sphere(Vector3d ray_orig, Vector3d ray_dir, Vector3d center, double radius,
    std::shared_ptr<Vector3d> intersection_point, std::shared_ptr<Vector3d> normal)
{
    // Moving the origin
    Vector3d c = ray_orig - center;

    double discriminant = std::pow(ray_dir.dot(c), 2.0) - (ray_dir.dot(ray_dir) * (c.dot(c) - (radius * radius)));
    if (discriminant < 0.0) return false;

    double t = (-ray_dir.dot(c) - std::sqrt(discriminant)) / ray_dir.dot(ray_dir);
    if (t < 0.0) return false;
    
    Vector3d point = ray_orig + (t * ray_dir);

    (*intersection_point) = point;
    (*normal) = (point - center) / radius;
    return true;
}



void raytrace_shading()
{
    std::cout << "Simple ray tracer, one sphere with different shading" << std::endl;

    const std::string filename("shading.png");
    MatrixXd R = MatrixXd::Zero(800, 800); // Store red
    MatrixXd G = MatrixXd::Zero(800, 800); // Store green
    MatrixXd B = MatrixXd::Zero(800, 800); // Store blue
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is perspective, pointing in the direction -z and covering the unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / A.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / A.rows(), 0);

    //Sphere setup
    const Vector3d sphere_center(0, 0, 0);
    const double sphere_radius = 0.9;

    //material params
    const Vector3d diffuse_color(1, 0, 1);
    const double specular_exponent = 100;
    const Vector3d specular_color(0.0, 0, 1);

    // Single light source
    const Vector3d light_position(-1, 1, 1);
    double ambient = 0.1;

    for (unsigned i = 0; i < A.cols(); ++i)
    {
        for (unsigned j = 0; j < A.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            // Prepare the ray (origin point and direction)
            const Vector3d ray_origin = pixel_center;
            const Vector3d ray_direction = pixel_center - camera_origin;

            // Intersect with the sphere
            std::shared_ptr<Vector3d> ray_intersection = std::make_shared<Vector3d>(0, 0, 0);
            std::shared_ptr<Vector3d> ray_normal = std::make_shared<Vector3d>(0, 0, 0);

            if (intersect_sphere(ray_origin, ray_direction, sphere_center, sphere_radius,
                ray_intersection, ray_normal))
            {
                // TODO: Add shading parameter here
                const Vector3d view_direction = -ray_direction.normalized();
                const Vector3d light_direction = (light_position - (*ray_intersection)).normalized();
                const Vector3d half_angle = (view_direction + light_direction).normalized();

                const Vector3d diffuse = diffuse_color * light_direction.dot(*ray_normal);
                const Vector3d specular = specular_color * std::pow(half_angle.dot(*ray_normal), specular_exponent);

                // Simple diffuse model
                R(i, j) = ambient + std::max(diffuse.x(), 0.0) + std::max(specular.x(), 0.0);
                G(i, j) = ambient + std::max(diffuse.y(), 0.0) + std::max(specular.y(), 0.0);
                B(i, j) = ambient + std::max(diffuse.z(), 0.0) + std::max(specular.z(), 0.0);

                // Disable the alpha mask for this pixel
                A(i, j) = 1;
            }
        }
    }

    // Save to png
    write_matrix_to_png(R, G, B, A, filename);
}



int main()
{
    raytrace_sphere();
    raytrace_parallelogram();
    raytrace_perspective();
    raytrace_shading();

    return 0;
}
