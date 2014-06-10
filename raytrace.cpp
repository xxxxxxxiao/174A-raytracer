//
// raytrace.cpp
// Thomas Lutton
// UCLA CS 174A 
// Spring 2014
// Section 1B
//

#define _CRT_SECURE_NO_WARNINGS
#include "matm.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

int g_width;
int g_height;

struct Ray
{
    vec4 origin;
    vec4 dir;
};

struct Sphere {
	string m_name;
	vec4 m_pos;
	vec4 m_color;
	vec3 m_scale;
	float m_k_a, m_k_d, m_k_s, m_k_r;
	float m_n;
};

struct Light {
	string m_name;
	vec4 m_pos;
	vec4 m_color;
};

vector<vec4> g_colors;
vector<Sphere> g_spheres;
vector<Light> g_lights;

vec4 g_ambient_color;
vec4 g_back_color;

float g_left;
float g_right;
float g_top;
float g_bottom;
float g_near;

string g_output;


// -------------------------------------------------------------------
// Input file parsing

vec4 toVec4(const string& s1, const string& s2, const string& s3)
{
    stringstream ss(s1 + " " + s2 + " " + s3);
    vec4 result;
    ss >> result.x >> result.y >> result.z;
    result.w = 1.0f;
    return result;
}

// to be used when parsing m_scale for spheres. Allows easy usage of mat4 Scale from matm.h
vec3 toVec3(const string& s1, const string& s2, const string& s3)
{
	stringstream ss(s1 + " " + s2 + " " + s3);
	vec3 result;
	ss >> result.x >> result.y >> result.z;
	return result;
}

float toFloat(const string& s)
{
    stringstream ss(s);
    float f;
    ss >> f;
    return f;
}



void parseLine(const vector<string>& vs)
{
    if (vs[0] == "RES") {
        g_width = (int)toFloat(vs[1]);
        g_height = (int)toFloat(vs[2]);
        g_colors.resize(g_width * g_height);
    } else if(vs[0] == "NEAR") {
		g_near = toFloat(vs[1]);
	} else if(vs[0] == "LEFT") {
		g_left = toFloat(vs[1]);
	} else if(vs[0] == "RIGHT") {
		g_right = toFloat(vs[1]);
	} else if(vs[0] == "BOTTOM") {
		g_bottom = toFloat(vs[1]);
	} else if(vs[0] == "TOP") {
		g_top = toFloat(vs[1]);
	} else if(vs[0] == "SPHERE") {
		Sphere sphere;
		sphere.m_name  = vs[1];
		sphere.m_pos   = toVec4(vs[2], vs[3], vs[4]);
		sphere.m_scale = toVec3(vs[5], vs[6], vs[7]);
		sphere.m_color = toVec4(vs[8], vs[9], vs[10]);
		sphere.m_k_a   = toFloat(vs[11]);
		sphere.m_k_d   = toFloat(vs[12]);
		sphere.m_k_s   = toFloat(vs[13]);
		sphere.m_k_r   = toFloat(vs[14]);
		sphere.m_n     = toFloat(vs[15]);
		g_spheres.push_back(sphere);
	} else if(vs[0] == "LIGHT") {
		Light light;
		light.m_name  = vs[1];
		light.m_pos   = toVec4(vs[2], vs[3], vs[4]);
		light.m_color = toVec4(vs[5], vs[6], vs[7]);
		g_lights.push_back(light);
	} else if(vs[0] == "BACK") {
		g_back_color = toVec4(vs[1], vs[2], vs[3]);
	} else if(vs[0] == "AMBIENT") {
		g_ambient_color = toVec4(vs[1], vs[2], vs[3]);
	} else if(vs[0] == "OUTPUT") {
		g_output = vs[1];
	}
}

void loadFile(const char* filename)
{
    ifstream is(filename);
    if (is.fail())
    {
        cout << "Could not open file " << filename << endl;
        exit(1);
    }
    string s;
    vector<string> vs;
    while(!is.eof())
    {
        vs.clear();
        getline(is, s);
        istringstream iss(s);
        while (!iss.eof())
        {
            string sub;
            iss >> sub;
            vs.push_back(sub);
        }
        parseLine(vs);
    }
}


// -------------------------------------------------------------------
// Utilities

void setColor(int ix, int iy, const vec4& color)
{
    int iy2 = g_height - iy - 1; // Invert iy coordinate.
    g_colors[iy2 * g_width + ix] = color;
}


// -------------------------------------------------------------------
// Intersection routine

vec4 intersect(const Ray& ray, int& sphere_index, int& intersect_type)
{
	vec4 return_point(0.0f, 0.0f, 0.0f, 1.0f);
	float t_h = FLT_MAX;
	for(int i = 0; i < g_spheres.size(); i++)
	{
		// Get matrix to invert scale of sphere (transform sphere back to unit sphere)
		mat4 invertedScale;
		InvertMatrix(Scale(g_spheres[i].m_scale), invertedScale); //use Scale function from matm.h to create mat4 scale matrix

		// Get vector of sphere pos w.r.t origin
		vec4 ots = g_spheres[i].m_pos - ray.origin;

		// |C|^2 * t^2 + 2(S dot C) + ( |S|^2 - 1 )
		// Use quadratic formula to solve for t1 and t2
		vec4 S = invertedScale * ots;
		vec4 C = invertedScale * ray.dir;

		float a = dot(C,C);
		float b = dot(S,C);
		float c = dot(S,S) - 1;

		// if discriminant >= 0, 1 or 2 real solutions
		float discriminant = b*b - a*c;
		if(discriminant >= 0)
		{
			float r = sqrt(discriminant) / a;

			float t1 = b/a - r;
			float t2 = b/a + r;

			// t1 > near plane and less than t_h (meaning closest intersection point)
			if(t1 > 1.0f && t1 < t_h)
			{
				t_h = t1;
				return_point = ray.origin + t_h * ray.dir;
				if(t2 > 1.0f) // sphere is not intersected by near plane
					intersect_type = 1;
				else // sphere intersected by near plane
					intersect_type = 2;
				sphere_index = i;
			}

			if(t2 > 1.0f && t2 < t_h)
			{
				t_h = t2;
				return_point = ray.origin + t_h * ray.dir;
				if(t1 > 1.0f)
					intersect_type = 1;
				else
					intersect_type = 2;
				sphere_index = i;
			}
		}
	}

	return return_point;
}


// -------------------------------------------------------------------
// Ray tracing

vec4 trace(const Ray& ray, int recursive_level)
{
	if(recursive_level > 3)
		return g_back_color;

	// intersect_type == 0, no intersection
	// intersect_type == 1, regular sphere
	// intersect_type == 2, hollow/middle of sphere (part of sphere intersects image plane)
	int sphere_index = -1, intersect_type = 0;
    vec4 point = intersect(ray, sphere_index, intersect_type);

	if(intersect_type == 0) // if no intersection
		return g_back_color;

	// Compute ambient component of pixel color
	vec4 pixel_color = g_spheres[sphere_index].m_k_a * g_spheres[sphere_index].m_color * g_ambient_color;

	//Get normal at the point using the gradient method
	vec4 normal_at_point = point - g_spheres[sphere_index].m_pos;
	mat4 invertedScale;
	InvertMatrix(Scale(g_spheres[sphere_index].m_scale), invertedScale);
	normal_at_point.w = 0;
	normal_at_point = normalize(2 * invertedScale * invertedScale * normal_at_point);


	vec4 diffuse_color(0.0f, 0.0f, 0.0f, 0.0f);
	vec4 specular_color(0.0f, 0.0f, 0.0f, 0.0f);
	
	// Compute v (vector from point to eye) for use with specular component
	vec4 v = ray.origin - point;
	v = normalize(v);

	for(int i=0; i < g_lights.size(); i++)
	{
		Ray shadow_ray;
		shadow_ray.origin = point;
		shadow_ray.dir = normalize(g_lights[i].m_pos - point);

		int si, it = 0;
		intersect(shadow_ray, si, it); //Look for intersections with other spheres between point and light[i]
		if(it == 0)  // No intersection with sphere
		{
			// Compute diffuse component of pixel color
			float ndotl = dot(normal_at_point, shadow_ray.dir);
			if(ndotl < 0) // no contribution to diffuse light
				continue;

			diffuse_color += ndotl * g_lights[i].m_color * g_spheres[sphere_index].m_color * g_spheres[sphere_index].m_k_d;

			// If sphere is hollow/intersecting near plane no specular component
			if(intersect_type == 2)
				continue;

			// Compute specular component of pixel_color
			vec4 r = (2 * dot(normal_at_point, shadow_ray.dir) * normal_at_point) - shadow_ray.dir;
			r = normalize(r);

			specular_color += powf(dot(r, v), g_spheres[sphere_index].m_n) * g_lights[i].m_color * g_spheres[sphere_index].m_k_s;
		}
	}
	
	pixel_color += diffuse_color + specular_color;

	if(intersect_type == 2) // We're actually inside the sphere. Don't reflect
		return pixel_color;

	// Compute reflected ray to be recursively passed to trace
	Ray reflected;
	reflected.origin = point;
	reflected.dir = normalize(ray.dir - 2.0f * dot(normal_at_point, ray.dir) * normal_at_point);

	vec4 reflection_color = trace(reflected, recursive_level + 1);

	// if reflection_color returned from recursion != background color, then it has meaningful contribution to color
	if(reflection_color.x != g_back_color.x || reflection_color.y != g_back_color.y || reflection_color.z != g_back_color.z)
		pixel_color += reflection_color * g_spheres[sphere_index].m_k_r;

	return pixel_color;
}

vec4 getDir(int ix, int iy)
{
    // This should return the direction from the origin
    // to pixel (ix, iy), normalized.
	float x = g_left + (g_right - g_left) * ( (float)ix / ((float)g_width - 1));
	float y = g_bottom + (g_top - g_bottom) * ( (float)iy / ((float)g_height - 1));

    return normalize(vec4(x, y, -g_near, 0.0f) );
}

void renderPixel(int ix, int iy)
{
    Ray ray;
    ray.origin = vec4(0.0f, 0.0f, 0.0f, 1.0f);
    ray.dir = getDir(ix, iy);
    vec4 color = trace(ray, 0);
    setColor(ix, iy, color);
}

void render()
{
    for (int iy = 0; iy < g_height; iy++)
        for (int ix = 0; ix < g_width; ix++)
            renderPixel(ix, iy);
}


// -------------------------------------------------------------------
// PPM saving

void savePPM(int Width, int Height, const char* fname, unsigned char* pixels) 
{
    FILE *fp;
    const int maxVal=255;

    printf("Saving image %s: %d x %d\n", fname, Width, Height);
    fp = fopen(fname,"wb");
    if (!fp) {
        printf("Unable to open file '%s'\n", fname);
        return;
    }
    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", Width, Height);
    fprintf(fp, "%d\n", maxVal);

    for(int j = 0; j < Height; j++) {
        fwrite(&pixels[j*Width*3], 3, Width, fp);
    }

    fclose(fp);
}

void saveFile()
{
    // Convert color components from floats to unsigned chars.
    // clamp values if out of range.
    unsigned char* buf = new unsigned char[g_width * g_height * 3];
    for (int y = 0; y < g_height; y++)
        for (int x = 0; x < g_width; x++)
            for (int i = 0; i < 3; i++) {
				float g = ((float*)g_colors[y*g_width+x])[i];
				unsigned char uc;
				if(g > 1.0f) // clamp g and scale by 255
					uc = (unsigned char)(255.9f);
				else if(g < 0.0f)
					uc = (unsigned char)(0.0f);
				else
					uc = (unsigned char)(g * 255.9f);
				buf[y*g_width*3+x*3+i] = uc;
			}
    savePPM(g_width, g_height, g_output.c_str(), buf);
    delete[] buf;
}


// -------------------------------------------------------------------
// Main

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        cout << "Usage: template-rt <input_file.txt>" << endl;
        exit(1);
    }
    loadFile(argv[1]);
    render();
    saveFile();
	return 0;
}

