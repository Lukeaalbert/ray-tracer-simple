/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: Luke Albert
 * *************************
*/

#ifdef WIN32
  #include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
  #include <GL/gl.h>
  #include <GL/glut.h>
#elif defined(__APPLE__)
  #include <OpenGL/gl.h>
  #include <GLUT/glut.h>
#endif

#include <cmath>
#include <vector>
#include <cfloat>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef WIN32
  #define strcasecmp _stricmp
#endif

#include <imageIO.h>

#include "vec_math.hpp" // my own custom include

#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100

char * filename = NULL;

// The different display modes.
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY;

// While solving the homework, it is useful to make the below values smaller for debugging purposes.
// The still images that you need to submit with the homework should be at the below resolution (640x480).
// However, for your own purposes, after you have solved the homework, you can increase those values to obtain higher-resolution images.
#define WIDTH 640
#define HEIGHT 480

// The field of view of the camera, in degrees.
#define fov 60.0

// Buffer to store the image.
unsigned char buffer[HEIGHT][WIDTH][3];

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

struct Triangle
{
  Vertex v[3];

void project_2d(const double triangle_normal[3],
  const double project_point_3d[3], double a_projected[2],
  double b_projected[2], double c_projected[2],
  double p_projected[2])
  {
    // pick projection axis based on the largest component of normal
    double normal_x = fabs(triangle_normal[0]);
    double normal_y = fabs(triangle_normal[1]);
    double normal_z = fabs(triangle_normal[2]);

    // kill x
    if (normal_x > normal_y && normal_x > normal_z) {
        a_projected[0] = v[0].position[1];
        a_projected[1] = v[0].position[2];
        b_projected[0] = v[1].position[1];
        b_projected[1] = v[1].position[2];
        c_projected[0] = v[2].position[1];
        c_projected[1] = v[2].position[2];
        p_projected[0] = project_point_3d[1];
        p_projected[1] = project_point_3d[2];
    }
    // kill y
    else if (normal_y > normal_z) {
        a_projected[0] = v[0].position[0];
        a_projected[1] = v[0].position[2];
        b_projected[0] = v[1].position[0];
        b_projected[1] = v[1].position[2];
        c_projected[0] = v[2].position[0];
        c_projected[1] = v[2].position[2];
        p_projected[0] = project_point_3d[0];
        p_projected[1] = project_point_3d[2];
    }
    // kill z
    else {
        a_projected[0] = v[0].position[0];
        a_projected[1] = v[0].position[1];
        b_projected[0] = v[1].position[0];
        b_projected[1] = v[1].position[1];
        c_projected[0] = v[2].position[0];
        c_projected[1] = v[2].position[1];
        p_projected[0] = project_point_3d[0];
        p_projected[1] = project_point_3d[1];
    }
  }

  void barycentric(const double A[2], const double B[2], const double C[2],
    const double P[2], double& alpha, double& beta, double& gamma) {
    double v0[2], v1[2], v2[2];
    v0[0] = B[0] - A[0]; 
    v0[1] = B[1] - A[1];
    v1[0] = C[0] - A[0]; 
    v1[1] = C[1] - A[1];
    v2[0] = P[0] - A[0]; 
    v2[1] = P[1] - A[1];

    double d00 = v0[0] * v0[0] + v0[1] * v0[1];
    double d01 = v0[0] * v1[0] + v0[1] * v1[1];
    double d11 = v1[0] * v1[0] + v1[1] * v1[1];
    double d20 = v2[0] * v0[0] + v2[1] * v0[1];
    double d21 = v2[0] * v1[0] + v2[1] * v1[1];

    double denom = d00 * d11 - d01 * d01;
    beta  = (d11 * d20 - d01 * d21) / denom;
    gamma = (d00 * d21 - d01 * d20) / denom;
    alpha = 1.0 - beta - gamma;
  }

  // t is the distance from the ray origin to intersection point
  bool intersect(const double ray_origin[3], const double ray_direction[3], double& t)
  {
    // triangle normal
    double n[3];
    get_triangle_normal(n);

    // ray is parallel to plane formed by triangle vertices
    if (dot(n, ray_direction) == 0) {
      return false;
    }

    // solve for t = −[ (n dot p0 + −n dot A) / (n dot direction) ]
    {
      // numerator
      double numerator = -dot(n, v[0].position) + dot(n, ray_origin);
      // denominator
      double denominator = dot(n, ray_direction);
      t = -(numerator/denominator);
    }

    // intersection is behind ray
    if (t <= 0) {
      return false;
    }

    // at this point, we have a valid ray PLANE intersection
    // for the plane formed by the 3 points of the triangle.
    // now we have to use barycentric coordinates to make sure
    // that the intersetion point is actually inside the triangle.

    // finding intersection point
    double intersection_point[3];
    double ray_direction_scaled[3];
    deep_copy(ray_direction_scaled, ray_direction);
    scale(ray_direction_scaled, t);
    add(ray_origin, ray_direction_scaled, intersection_point);

    // we want to flatten an axis to preform barycentric coordinate
    // calculations. we have to do this dynamically though, to make sure
    // the triangle isn't perpendicular to the axis. we'll try the z-axis
    // first and waterfall to the y-axis as needed.

    // project 2-d
    double A[2], B[2], C[2], P[2];
    project_2d(n, intersection_point, A, B, C, P);
    double alpha, beta, gamma;
    barycentric(A, B, C, P, alpha, beta, gamma);

    // p inside triangle iff 0 ≤ α, β, γ ≤ 1, α + β + γ = 1 
    return (0 <= alpha && alpha <= 1)
      && (0 <= beta && beta <= 1)
      && (0 <= gamma && gamma <= 1)
      && (abs(alpha + beta + gamma - 1) < 1e-6);
  }

  void get_triangle_normal(double normal[3]) {
    // n = (B−A)×(C−A)
    // where A, B, and C are triangle vertices
    double b_sub_a[3];
    subtract(v[1].position, v[0].position, b_sub_a);
    double c_sub_a[3];
    subtract(v[2].position, v[0].position, c_sub_a);
    cross(b_sub_a, c_sub_a, normal);
    normalize(normal);
  }
};

struct Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;

  // t is the distance from the ray origin to intersection point
  bool intersect(const double ray_origin[3], const double ray_direction[3], double& t)
  {
    // vector from ray origin to sphere center is just negative position vector.
    double origin_to_sphere_center[3];
    subtract(ray_origin, position, origin_to_sphere_center);

    // quadratic equation: t^2 + bt + c = 0
    // assumes a = 1 (direction is normalized)

    // b = 2 * (ray_direction dot origin_to_sphere_center)
    double b = 2.0 * dot(ray_direction, origin_to_sphere_center);

    // c = (origin_to_sphere_center)^2 - radius^2
    double c = dot(origin_to_sphere_center, origin_to_sphere_center) - (radius * radius);

    // discriminant
    double discriminant = b * b - 4.0 * c;

    // ignore when intersection if discriminant is negative
    if (discriminant < 0.0) return false;

    double sqrt_discriminant = sqrt(discriminant);

    // two solutions: t0 and t1. prioritize the closest one.
    double t0 = (-b - sqrt_discriminant) / 2.0;
    double t1 = (-b + sqrt_discriminant) / 2.0;
    if (t0 > 0.0) {
      t = t0;
      return true;
    }
    else if (t1 > 0.0) {
      t = t1;
      return true;
    }

    return false;
  }

  void get_normal_at_intersection_point(const double intersection_point[3], double normal[3])
  {
    normal[0] = (intersection_point[0] - position[0]) / radius;
    normal[1] = (intersection_point[1] - position[1]) / radius;
    normal[2] = (intersection_point[2] - position[2]) / radius;
  }

};

struct Light
{
  double position[3];
  double color[3];
};

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;

// Populated in compute_eye_to_image_plane_rays().
// Every entry has an x, y, z coordinate.
double unit_rays[WIDTH][HEIGHT][3];

enum ObjectType {TRIANGLE, SPHERE, LIGHT, UNDEFINED};

struct TracedRay {
  // unit ray casted to yield intersection
  double unit_ray[3];
  // ray object intersection data
  bool intersects_object;
  double intersection_point[3];
  ObjectType object_type;
  double object_idx;
  double t;
  // object light intersection data
  std::vector<bool> is_illuminated_by_light; // whether object is illuminated by light at light idx i
  TracedRay() : intersects_object(false),
    object_type(UNDEFINED),
    object_idx(-1.0),
    t(DBL_MAX)
  {
    unit_ray[0] = 0.0;
    unit_ray[1] = 0.0;
    unit_ray[2] = 0.0;
    intersection_point[0] = 0.0;
    intersection_point[1] = 0.0;
    intersection_point[2] = 0.0;
  }
};

TracedRay traced_rays[WIDTH][HEIGHT];

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

// populates unit_rays global variable with unit rays from eye to image plane
void compute_unit_rays() {
  // Calculate image plane attributes. This should be done
  // before calculating the actual rays from the origin to the
  // image plane (obviously).
  double fov_radians = fov * M_PI / 180.0f;
  double aspect_ratio = ((double)WIDTH/(double)HEIGHT);
  double image_plane_height = 2.0f * tan(fov_radians/2);
  double image_plane_width = 2.0f * (aspect_ratio * tan(fov_radians/2));
  // Find rays from origin to the image plane.
  for (double x=0.0f; x<WIDTH; x += 1.0f)
  {
    // -1 to 1
    double u = ( (x + 0.5) / WIDTH ) * 2 - 1;
    // x coordinate for point on image plane 
    double px = u * (image_plane_width/2);
    for (double y=0.0f; y<HEIGHT; y += 1.0f)
    {
      // -1 to 1
      double v = ((y + 0.5) / HEIGHT) * 2 - 1;
      // y coordinate for point on image plane 
      double py = v * (image_plane_height/2);
      // y coordinate for point on image plane (constant -1)
      double pz = -1.0f;
      // construct and normalize unit ray
      double unit_ray[3] = {px, py, pz};
      normalize(unit_ray);
      // insert in unit_rays array
      unit_rays[(int)x][(int)y][0] = unit_ray[0];
      unit_rays[(int)x][(int)y][1] = unit_ray[1];
      unit_rays[(int)x][(int)y][2] = unit_ray[2];
    }
  }
}

// casts all unit rays and finds the cloest object they intersects.
// populates 'traced_rays' global var.
void cast_primary_rays() {
  double origin[3] = {0.0, 0.0, 0.0};
  // iterate through all rays
  for (unsigned int i = 0; i < WIDTH; ++i) {
    for (int j = 0; j < HEIGHT; ++j) {
      // get current ray
      double curr_ray[3];
      deep_copy(curr_ray, unit_rays[i][j]);
      // variables used for tracking the closest intersection
      ObjectType closest_object_type = UNDEFINED; 
      double closest_t = DBL_MAX;
      int closest_obj_idx = -1;
      int idx = 0; 
      // looping condition (while there are still shapes to explore)
      while (idx < num_spheres || idx < num_triangles || idx < num_lights) {
        // sphere intersection
        if (idx < num_spheres) {
          double t;
          if (spheres[idx].intersect(origin, curr_ray, t)) {
            // intersection is valid and less is new closest intersection
            if (t < closest_t) {
              closest_t = t;
              closest_object_type = SPHERE;
              closest_obj_idx = idx;
            }
          }
        }
        // triangle intersection
        if (idx < num_triangles) {
          double t;
          if (triangles[idx].intersect(origin, curr_ray, t)) {
            if (t < closest_t) {
              closest_t = t;
              closest_object_type = TRIANGLE;
              closest_obj_idx = idx;
            }
          }
        }
        // light intersection
        if (idx < num_lights) {
          // TODO
        }
        // increment idx (always!)
        idx++;
      }
      TracedRay traced_ray;
      // valid ray object intersection
      if (closest_t != DBL_MAX) {
        traced_ray.intersects_object = true;
        traced_ray.object_type = closest_object_type;
        traced_ray.object_idx = closest_obj_idx;
        traced_ray.t = closest_t;
        traced_ray.intersection_point[0] = closest_t * curr_ray[0];
        traced_ray.intersection_point[1] = closest_t * curr_ray[1];
        traced_ray.intersection_point[2] = closest_t * curr_ray[2];
        traced_ray.unit_ray[0] = curr_ray[0];
        traced_ray.unit_ray[1] = curr_ray[1];
        traced_ray.unit_ray[2] = curr_ray[2];
      }
      traced_rays[i][j] = traced_ray;
    }
  }
}

// populates 'is_illuminated_by_light' fields in 'traced_rays'.
void cast_shadow_rays() {
  for (unsigned int x = 0; x < WIDTH; ++x) {
    for (unsigned int y = 0; y < HEIGHT; ++y) {
      // skip if no intersection
      if (!traced_rays[x][y].intersects_object) {
        continue;
      }
      // resize illumination vector
      traced_rays[x][y].is_illuminated_by_light.resize(num_lights, false);
      // iterate through lights
      for (int l = 0; l < num_lights; ++l) {
        Light light = lights[l];
        
        // direction from intersection to light
        double intersection_to_curr_light[3];
        subtract(light.position, traced_rays[x][y].intersection_point, intersection_to_curr_light);
        double distance_to_light = magnitude(intersection_to_curr_light);
        normalize(intersection_to_curr_light);
        
        // to check if path to light is blocked
        bool intersection_found = false;
        
        // spheres checks
        for (int idx = 0; idx < num_spheres; idx++) {
          if (idx == traced_rays[x][y].object_idx
            && traced_rays[x][y].object_type == SPHERE) continue;
          double t;
          if (spheres[idx].intersect(traced_rays[x][y].intersection_point, intersection_to_curr_light, t)) {
            if (t > 0 && t < distance_to_light) {
              intersection_found = true;
              break;
            }
          }
        }
        
        // triangles checks
        if (!intersection_found) {
          for (int idx = 0; idx < num_triangles; idx++) {
            if (idx == traced_rays[x][y].object_idx
              && traced_rays[x][y].object_type == TRIANGLE) continue;
              double t;
              if (triangles[idx].intersect(traced_rays[x][y].intersection_point, intersection_to_curr_light, t)) {
                if (t > 0 && t < distance_to_light) {
                intersection_found = true;
                break;
              }
            }
          }
        }
        
        // if no blocking object found, point is illuminated
        if (!intersection_found) {
          traced_rays[x][y].is_illuminated_by_light[l] = true;
        }
      }
    }
  }
}

void compute_phong_color(TracedRay& ray, unsigned char& r, unsigned char& g, unsigned char& b) {
  // return background color if no intersection
  if (!ray.intersects_object) {
    r = 0;
    g = 0;
    b = 0;
    return;
  }
  
  // get surface properties (dynamically, based on object type) at intersection point
  double normal[3];
  double kd[3];
  double ks[3];
  double shininess;
  // sphere surface properties case
  if (ray.object_type == SPHERE) {
    int idx = ray.object_idx;
    spheres[idx].get_normal_at_intersection_point(ray.intersection_point, normal);
    
    // diffuse
    kd[0] = spheres[idx].color_diffuse[0];
    kd[1] = spheres[idx].color_diffuse[1];
    kd[2] = spheres[idx].color_diffuse[2];
    // specular
    ks[0] = spheres[idx].color_specular[0];
    ks[1] = spheres[idx].color_specular[1];
    ks[2] = spheres[idx].color_specular[2];
    
    shininess = spheres[idx].shininess;
  }
  // triangle surface properties case
  else if (ray.object_type == TRIANGLE) {
    int idx = ray.object_idx;
    Triangle& tri = triangles[idx];

    tri.get_triangle_normal(normal);

    // get barycentric coordinates
    double A[2], B[2], C[2], P[2];
    tri.project_2d(normal, ray.intersection_point, A, B, C, P);
    double alpha, beta, gamma;
    tri.barycentric(A, B, C, P, alpha, beta, gamma);

    // diffuse interpolation
    kd[0] = alpha * tri.v[0].color_diffuse[0] + beta
      * tri.v[1].color_diffuse[0] + gamma * tri.v[2].color_diffuse[0];
    kd[1] = alpha * tri.v[0].color_diffuse[1] + beta
      * tri.v[1].color_diffuse[1] + gamma * tri.v[2].color_diffuse[1];
    kd[2] = alpha * tri.v[0].color_diffuse[2] + beta
      * tri.v[1].color_diffuse[2] + gamma * tri.v[2].color_diffuse[2];

    // specular interpolation
    ks[0] = alpha * tri.v[0].color_specular[0] + beta
      * tri.v[1].color_specular[0] + gamma * tri.v[2].color_specular[0];
    ks[1] = alpha * tri.v[0].color_specular[1] + beta
      * tri.v[1].color_specular[1] + gamma * tri.v[2].color_specular[1];
    ks[2] = alpha * tri.v[0].color_specular[2] + beta
      * tri.v[1].color_specular[2] + gamma * tri.v[2].color_specular[2];

    // shininess interpolation
    shininess = alpha * tri.v[0].shininess + beta
      * tri.v[1].shininess + gamma * tri.v[2].shininess;
  }

  // final color (populated by phong shading)
  double final_color[3];
  
  // ambient component. ALWAYS ADDED
  final_color[0] = ambient_light[0] * kd[0];
  final_color[1] = ambient_light[1] * kd[1];
  final_color[2] = ambient_light[2] * kd[2];
  
  // v = view direction vector from intersection point to origin (0,0,0)
  double v[3];
  v[0] = -ray.intersection_point[0];
  v[1] = -ray.intersection_point[1];
  v[2] = -ray.intersection_point[2];
  normalize(v);
  
  // add contributions from each light source that illuminates this point
  for (int l = 0; l < num_lights; l++) {
    if (!ray.is_illuminated_by_light[l]) {
      continue;
    }
    
    Light& light = lights[l];
    
    // compute light direction ld
    double ld[3];
    subtract(light.position, ray.intersection_point, ld);
    normalize(ld);
    
    // compute normal dot ld for diffuse
    // (i.e., diffuse term = normal_dot_ld)
    double normal_dot_ld = dot(normal, ld);
    // clamp to 0 if negative
    if (normal_dot_ld < 0.0) normal_dot_ld = 0.0;
    
    // reflection direction
    double reflection[3];
    reflection[0] = 2.0 * normal_dot_ld * normal[0] - ld[0];
    reflection[1] = 2.0 * normal_dot_ld * normal[1] - ld[1];
    reflection[2] = 2.0 * normal_dot_ld * normal[2] - ld[2];
    normalize(reflection);
    
    // (reflection dot v) for specular term
    double reflection_dot_v = dot(reflection, v);
    // again, clamp if negative
    if (reflection_dot_v < 0.0) reflection_dot_v = 0.0;
    
    // specular term: (reflection_dot_v)^shininess
    double specular_term = pow(reflection_dot_v, shininess);
    
    // Finally, apply the phong formula
    final_color[0] += light.color[0] * (kd[0] * normal_dot_ld + ks[0] * specular_term);
    final_color[1] += light.color[1] * (kd[1] * normal_dot_ld + ks[1] * specular_term);
    final_color[2] += light.color[2] * (kd[2] * normal_dot_ld + ks[2] * specular_term);
  }
  
  // clamp final color from 0.0 min to 1.0 max
  if (final_color[0] > 1.0) final_color[0] = 1.0;
  if (final_color[1] > 1.0) final_color[1] = 1.0;
  if (final_color[2] > 1.0) final_color[2] = 1.0;
  if (final_color[0] < 0.0) final_color[0] = 0.0;
  if (final_color[1] < 0.0) final_color[1] = 0.0;
  if (final_color[2] < 0.0) final_color[2] = 0.0;
  
  // convert to unsigned char (what is used in plot_pixel function
  // that uses these RBG vals)
  r = (unsigned char)(final_color[0] * 255);
  g = (unsigned char)(final_color[1] * 255);
  b = (unsigned char)(final_color[2] * 255);
}

void draw_scene()
{
  for(unsigned int x = 0; x < WIDTH; x++)
  {
    glPointSize(2.0);  
    glBegin(GL_POINTS);
    for(unsigned int y = 0; y < HEIGHT; y++)
    {
      unsigned char r, g, b;
      compute_phong_color(traced_rays[x][y], r, g, b);
      plot_pixel(x, y, r, g, b);
    }
    glEnd();
    glFlush();
  }
  printf("Ray tracing completed.\n"); 
  fflush(stdout);
}

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  buffer[y][x][0] = r;
  buffer[y][x][1] = g;
  buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
    plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
  printf("Saving JPEG file: %s\n", filename);

  ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
  if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
    printf("Error in saving\n");
  else 
    printf("File saved successfully\n");
}

void parse_check(const char *expected, char *found)
{
  if(strcasecmp(expected,found))
  {
    printf("Expected '%s ' found '%s '\n", expected, found);
    printf("Parsing error; abnormal program abortion.\n");
    exit(0);
  }
}

void parse_doubles(FILE* file, const char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE *file, double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE *file, double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE * file = fopen(argv,"r");
  if (!file)
  {
    printf("Unable to open input file %s. Program exiting.\n", argv);
    exit(0);
  }

  int number_of_objects;
  char type[50];
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i", &number_of_objects);

  printf("number of objects: %i\n",number_of_objects);

  parse_doubles(file,"amb:",ambient_light);

  for(int i=0; i<number_of_objects; i++)
  {
    fscanf(file,"%s\n",type);
    printf("%s\n",type);
    if(strcasecmp(type,"triangle")==0)
    {
      printf("found triangle\n");
      for(int j=0;j < 3;j++)
      {
        parse_doubles(file,"pos:",t.v[j].position);
        parse_doubles(file,"nor:",t.v[j].normal);
        parse_doubles(file,"dif:",t.v[j].color_diffuse);
        parse_doubles(file,"spe:",t.v[j].color_specular);
        parse_shi(file,&t.v[j].shininess);
      }

      if(num_triangles == MAX_TRIANGLES)
      {
        printf("too many triangles, you should increase MAX_TRIANGLES!\n");
        exit(0);
      }
      triangles[num_triangles++] = t;
    }
    else if(strcasecmp(type,"sphere")==0)
    {
      printf("found sphere\n");

      parse_doubles(file,"pos:",s.position);
      parse_rad(file,&s.radius);
      parse_doubles(file,"dif:",s.color_diffuse);
      parse_doubles(file,"spe:",s.color_specular);
      parse_shi(file,&s.shininess);

      if(num_spheres == MAX_SPHERES)
      {
        printf("too many spheres, you should increase MAX_SPHERES!\n");
        exit(0);
      }
      spheres[num_spheres++] = s;
    }
    else if(strcasecmp(type,"light")==0)
    {
      printf("found light\n");
      parse_doubles(file,"pos:",l.position);
      parse_doubles(file,"col:",l.color);

      if(num_lights == MAX_LIGHTS)
      {
        printf("too many lights, you should increase MAX_LIGHTS!\n");
        exit(0);
      }
      lights[num_lights++] = l;
    }
    else
    {
      printf("unknown type in scene description:\n%s\n",type);
      exit(0);
    }
  }
  return 0;
}

void display()
{
}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  // Hack to make it only draw once.
  static int once=0;
  if(!once)
  {
    // Construct the rays from eye to image plane
    compute_unit_rays();
    cast_primary_rays();
    cast_shadow_rays();
    draw_scene();
    if(mode == MODE_JPEG)
      save_jpg();
  }
  once=1;
}

int main(int argc, char ** argv)
{
  if ((argc < 2) || (argc > 3))
  {  
    printf ("Usage: %s <input scenefile> [output jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
  {
    mode = MODE_JPEG;
    filename = argv[2];
  }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  #ifdef __APPLE__
    // This is needed on recent Mac OS X versions to correctly display the window.
    glutReshapeWindow(WIDTH - 1, HEIGHT - 1);
  #endif
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}

