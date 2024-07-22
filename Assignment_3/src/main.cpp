
////////////////////////////////////////////////////////////////////////////////
// C++ include
#include <iostream>
#include <string>
#include <vector>
#include <limits>


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
// Scene setup, global variables
////////////////////////////////////////////////////////////////////////////////

const std::string filename("raytrace.png");


//Camera settings
const double focal_length = 10;
const double field_of_view = 0.7854; //45 degrees
const double image_z = 5;
// const bool is_perspective = true;
// Doing it this way made it noticeably faster. The value can't be changed after compiling anyway
#define perspective
const Vector3d camera_position(0, 0, 5);

// The dof effect is really slow, and unfocused (the focal plane doesn't really change)
const double camera_aperture = 0.0; //0.05;
const int dof_samples = 1; //200;



// Maximum number of recursive calls
// Reflection and refraction use this independently
// So the true maximum is like 2^max_bounce or something idk
// Compile with optimizations on
const int max_bounce = 10;


// Objects
std::vector<Vector3d> sphere_centers;
std::vector<double> sphere_radii;
std::vector<Matrix3d> parallelograms;


//Material for the object, same material for all objects
const Vector4d obj_ambient_color(0.5, 0.1, 0.1, 0);
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




// Return the color of the ray
Vector4d shoot_ray(const Vector3d &ray_origin, const Vector3d &ray_direction, int max_bounce);





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



// Comment :-/
// From Assignment 2
double ray_triangle_intersection(Vector3d ray_orig, Vector3d ray_dir, Vector3d vert_a, Vector3d vert_b, Vector3d vert_c,
                        Vector3d &intersection_point, Vector3d &normal)
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

    // Ray is parallel
    if (fabs(c) < __DBL_EPSILON__) return -1;

    double inv_c = 1.0 / c;

    // Magic (Cramer's rule)
    double t = det_of_3_col_vectors(u, v, e) * inv_c;
    if (t < 0.0) return -1;

    double s = det_of_3_col_vectors(u, e, d) * inv_c;
    if (s < 0.0 || 1.0 < s) return -1;

    double r = det_of_3_col_vectors(e, v, d) * inv_c;
    if (r < 0.0 || 1.0 < r + s) return -1;

    intersection_point = ray_orig + (t * ray_dir);
    normal = v.cross(u).normalized();
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



//Finds the closest intersecting object returns its index
//In case of intersection it writes into p and N (intersection point and normals)
int find_nearest_object(const Vector3d &ray_origin, const Vector3d &ray_direction, Vector3d &p, Vector3d &N)
{
    // Find the object in the scene that intersects the ray first
    // we store the index and the 'closest_t' to their expected values
    int closest_index = -1;
    double closest_t = std::numeric_limits<double>::max(); //closest t is "+ infinity"

    Vector3d tmp_p, tmp_N;
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
                closest_index = i;
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
                closest_index = sphere_centers.size() + i;
                closest_t = t;
                p = tmp_p;
                N = tmp_N;
            }
        }
    }

    return closest_index;
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



Vector4d shoot_ray(const Vector3d &ray_origin, const Vector3d &ray_direction, int max_bounce)
{
    //Intersection point and normal, these are output of find_nearest_object
    Vector3d p, N;

    const int nearest_object = find_nearest_object(ray_origin, ray_direction, p, N);

    if (nearest_object < 0)
    {
        // Return a transparent color
        return Vector4d(0, 0, 0, 0);
    }

    N = N.normalized();
    Vector3d offset_p = p + (EPSILON * N);
    Vector3d internal_offset_p = p - (EPSILON * N);

    // Ambient light contribution
    const Vector4d ambient_color = obj_ambient_color.cwiseProduct(ambient_light);

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

        if (nearest_object == 4)
        {
            //Compute UV coordinates for the point on the sphere
            const double x = p(0) - sphere_centers[nearest_object][0];
            const double y = p(1) - sphere_centers[nearest_object][1];
            const double z = p(2) - sphere_centers[nearest_object][2];
            double tu = acos(z / sphere_radii[nearest_object]) / 3.1415;
            double tv = (3.1415 + atan2(y, x)) / (2 * 3.1415);
            tu = std::min(tu, 1.0);
            tu = std::max(tu, 0.0);

            tv = std::min(tv, 1.0);
            tv = std::max(tv, 0.0);

            diff_color = procedural_texture(tu, tv);
        }

        // Add shading parameters

        const Vector3d half_angle = (Li - ray_direction.normalized()).normalized();

        // Diffuse contribution
        const Vector4d diffuse = std::max(Li.dot(N), 0.0) * diff_color;

        // Specular contribution, use obj_specular_color
        const Vector4d specular = std::pow(half_angle.dot(N), obj_specular_exponent) * obj_specular_color;

        // Attenuate lights according to the squared distance to the lights
        const Vector3d D = light_position - p;
        lights_color += (diffuse + specular).cwiseProduct(light_color) / D.squaredNorm();
    }

    Vector4d refl_color = obj_reflection_color;
    if (nearest_object == 4)
    {
        refl_color = Vector4d(0.5, 0.5, 0.5, 0);
    }


    Vector3d reflection_direction = ray_direction - ((2.0 * ray_direction.dot(N)) * N);

    Vector3d refraction_direction;

    // Compute the color of the reflected ray and add its contribution to the current point color.
    // use refl_color
    Vector4d reflection_color = max_bounce > 0
        ? refl_color.cwiseProduct(shoot_ray(offset_p, reflection_direction, max_bounce - 1))
        : Vector4d(0, 0, 0, 0);

    // Compute the color of the refracted ray and add its contribution to the current point color.
    // Make sure to check for total internal reflection before shooting a new ray.
    Vector4d refraction_color = (max_bounce > 0 && calc_refraction(ray_direction, N, refraction_direction))
        ? obj_refraction_color.cwiseProduct(shoot_ray(internal_offset_p, refraction_direction, max_bounce - 1))
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

    int w = 1600;
    int h = 800;
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

    // The pixel grid through which we shoot rays is at a distance 'focal_length'.
    // Relative to camera.
    const Vector3d image_origin(-image_x, image_y, -image_z);
    const Vector3d x_displacement(2.0 / w * image_x, 0, 0);
    const Vector3d y_displacement(0, -2.0 / h * image_y, 0);

    #ifdef perspective
    const double inv_dof_samples = 1.0 / double(dof_samples);

    // C++ dares to maximize complexity where every other language cowers behind simplicity
    std::random_device rd;
    // The classic 32 bit Mersenne Twister! Of course! Classic!
    // Everyone has their favorite 623-Dimensionally Equidistributed Uniform Pseudorandom Number
    // Although it may be time for MT19937_64 to take its place as we are working with doubles
    // It's truly vital that we take full advantage of the super astronomical period of 2^19937 − 1
    // After all this program needs to be run less than 2.6971405e5994 times before it repeats a dof pattern (a thought that keeps me up at night)
    std::mt19937_64 eng(rd());
    std::uniform_real_distribution<> distr(0, camera_aperture); 
    #endif

    for (unsigned x = 0; x < w; ++x)
    {
        for (unsigned y = 0; y < h; ++y)
        {

            #ifdef perspective
            // TODO: Implement depth of field.
            // This is very much hacked in
            Vector4d C(0, 0, 0, 0);
            for (unsigned dof_sample = 0; dof_sample < dof_samples; ++dof_sample)
            {
                // Relative to camera.
                const Vector3d pixel_center = image_origin
                    + (x + 0.5) * x_displacement
                    + (y + 0.5) * y_displacement;
                
                double random_x = distr(eng) - (camera_aperture * 0.5);
                double random_y = distr(eng) - (camera_aperture * 0.5);

                Vector3d aperture_displacement = (random_x * x_displacement.normalized())
                    + (random_y * y_displacement.normalized());


                // Prepare the ray.
                // Relative to world.
                Vector3d ray_origin = camera_position + aperture_displacement;

                // Vector3d initial_ray_direction = (pixel_center - camera_position).normalized();
                // Vector3d focal_point = camera_position + (focal_distance * initial_ray_direction);
                // Vector3d ray_direction = focal_point + aperture_displacement - ray_origin;
                Vector3d ray_direction = pixel_center - camera_position;

                C += inv_dof_samples * shoot_ray(ray_origin, ray_direction, max_bounce);
            }
            #else
            // Relative to camera.
            const Vector3d pixel_center = image_origin + (x + 0.5) * x_displacement + (y + 0.5) * y_displacement;

            // Prepare the ray.
            // Relative to world.
            Vector3d ray_origin;
            Vector3d ray_direction;

            // Orthographic camera
            ray_origin = camera_position + Vector3d(pixel_center[0], pixel_center[1], 0);
            ray_direction = Vector3d(0, 0, -1);

            const Vector4d C = shoot_ray(ray_origin, ray_direction, max_bounce);
            #endif

            R(x, y) = C(0);
            G(x, y) = C(1);
            B(x, y) = C(2);
            A(x, y) = C(3);
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
