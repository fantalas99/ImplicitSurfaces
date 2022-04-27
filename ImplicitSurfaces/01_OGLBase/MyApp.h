#pragma once

// MarchingCubes
#include "Dependencies/implicit_marching_cubes-master/mc.h"

// OpenMesh
#include <OpenMesh/Core/Mesh/Traits.hh>
#include <OpenMesh/Core/Mesh/DefaultTriMesh.hh>
#include <OpenMesh/Core/Mesh/TriMeshT.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/Handles.hh>
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
	ProgramObject		m_programSkybox;	// skybox shader

	VertexArrayObject	m_CubeVao;			// VAO
	IndexBuffer			m_CubeIndices;		// index buffer
	ArrayBuffer			m_CubeVertexBuffer;	// VBO

	//VertexArrayObject		m_DebugVao;			// VAO
	//IndexBuffer			m_DebugIndices;		// index buffer
	//ArrayBuffer			m_DebugVertexBuffer;// VBO

	gCamera				m_camera;

	Texture2D			m_Texture_Red;
	Texture2D			m_Texture_Green;

	int TriangeCount;

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
	};

	std::unique_ptr<Mesh> m_mesh;

	typedef std::array<size_t, 3> Triangle;

	trimesh::trimesh_t tmesh;

	HalfEdgeMesh hmesh;

	std::vector<Function_T> Functions;

	glm::vec3 AutoDiffNormal(Geometry::Vector3D v);
	float AutoDiffCurvature(Geometry::Vector3D v);

	glm::vec3 EstimateNormal(int vertex_index);
	float EstimateCurvature(int vertex_index);

	double Det_2x2(double a1, double a2, double b1, double b2);

	void InitFunctions();
	void DrawCoordinateSystem();
	void InitShaders();
	void LoadMeshIntoBuffer(Geometry::TriMesh m);
	void DrawLine(glm::vec3 p_1, glm::vec3 p_2, float w);
	void DebugNormals();
	void DebugEdges();
	void GenerateMesh(int function_index);

	std::vector< std::array<size_t, 3>> GetAllTriangles(int vertex_id, Geometry::TriMesh mesh);

	glm::vec3 ConvertVector(Geometry::Vector3D v) {
	
		return glm::vec3(v[0], v[1], v[2]);

	}
	float GetTriangleArea(std::array<size_t, 3> tri, Geometry::TriMesh mesh);
	float GetTriangleAngle(Triangle t);

	bool UseNormalScale = true;
	bool ImguiDebugEdges = false;
	bool ImguiDebugNormals = false;
	bool VisualizeCurvature = true;
	bool UseAutoNormals = true;
	bool UseAutoCurvature = false;
	bool UseMeanCurvature = true;
	int ActiveFunctionIndex = 0;
	int DebugVertexCount = 0;
	std::vector<std::string> FunctionNames = {"Surface One","Surface Two", "Surface Three", "Surface Four"};


	struct MyTraits : public OpenMesh::DefaultTraits {
		using Point = OpenMesh::Vec3d; // the default would be Vec3f
		using Normal = OpenMesh::Vec3d;
		VertexTraits{
		  double mean;              // approximated mean curvature
		};
	};

};

