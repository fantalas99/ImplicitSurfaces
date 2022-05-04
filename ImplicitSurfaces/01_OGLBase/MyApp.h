#pragma once

// MarchingCubes
#include "Dependencies/implicit_marching_cubes-master/mc.h"

// Mesh
#include "Dependencies/TriMesh/trimesh.h"
#include "HalfEdgeMesh.h"


// C++ includes
#include <memory>


// GLEW
#include <GL/glew.h>

// SDL
#include <SDL.h>
#include <SDL_opengl.h>

// GLM
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform2.hpp>
#include "includes/gCamera.h"
#include "includes/ProgramObject.h"
#include "includes/BufferObject.h"
#include "includes/VertexArrayObject.h"
#include "includes/TextureObject.h"

// mesh
#include "includes/ObjParser_OGL3.h"

class HalfEdgeMesh;
class Function_T;

class ImplicitSurfaceApp
{
public:
	ImplicitSurfaceApp();
	~ImplicitSurfaceApp();

	bool Init();
	void Clean();

	void Update();
	void Render();

	void KeyboardDown(SDL_KeyboardEvent&);
	void KeyboardUp(SDL_KeyboardEvent&);
	void MouseMove(SDL_MouseMotionEvent&);
	void MouseDown(SDL_MouseButtonEvent&);
	void MouseUp(SDL_MouseButtonEvent&);
	void MouseWheel(SDL_MouseWheelEvent&);
	void Resize(int, int);

protected:

	ProgramObject		m_program;			// mesh shader

	VertexArrayObject	m_Vao;			// VAO
	IndexBuffer			m_Indices;		// index buffer
	ArrayBuffer			m_VertexBuffer;	// VBO

	gCamera				m_camera;

	int TriangeCount = 0;

	Geometry::TriMesh mesh;
	std::vector<glm::vec3> auto_normals;
	std::vector<glm::vec3> estimated_normals;
	std::vector<float> auto_curvatures;
	std::vector<float> estimated_curvatures;

	struct Vertex
	{
		glm::vec3 p;
		glm::vec3 n;
		glm::vec2 t;
		float c;
		float c_est;
	};

	typedef std::array<size_t, 3> Triangle;

	std::vector<Function_T> Functions;

	glm::vec3 AutoDiffNormal(Geometry::Vector3D v);
	float AutoDiffCurvature(Geometry::Vector3D v);

	glm::vec3 EstimateNormal(int vertex_index);
	float EstimateCurvature(int vertex_index);

	std::vector<Triangle> OrderTriangles(std::vector<Triangle> tris);

	double Det_2x2(double a1, double a2, double b1, double b2);

	void InitFunctions();
	void InitShaders();
	void LoadMeshIntoBuffer(Geometry::TriMesh m);
	void DebugNormals();
	void DebugEdges();
	void GenerateMesh(int function_index);

	std::vector< std::array<size_t, 3>> GetAllTriangles(int vertex_id);

	glm::vec3 ConvertVector(Geometry::Vector3D v) {
	
		return glm::vec3(v[0], v[1], v[2]);

	}
	float GetTriangleArea(std::array<size_t, 3> tri);
	float GetTriangleAngle(Triangle t);
	glm::vec3 GetTriangleNormal(Triangle t);
	bool IsValidTriangle(Triangle t);

	bool UseNormalScale = true;
	bool ImguiDebugEdges = false;
	bool ImguiDebugNormals = false;
	bool VisualizeCurvature = true;
	bool UseAutoNormals = true;
	bool UseAutoCurvature = false;
	bool UseMeanCurvature = true;
	int ActiveFunctionIndex = 0;
	int DebugVertexCount = 0;
	std::vector<std::string> FunctionNames = 
	{"Surface One ","Surface Two", "Surface Three", "Surface Four",
		"Surface One HD","Surface Two HD", "Surface Three HD", "Surface Four HD" };


};

