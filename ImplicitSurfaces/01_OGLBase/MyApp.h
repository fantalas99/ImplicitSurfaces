#pragma once

// MarchingCubes
#include "Dependencies/implicit_marching_cubes-master/mc.h"

// OpenMesh
//#include <OpenMesh/Core/Mesh/Traits.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
//#include <OpenMesh/Core/Mesh/Handles.hh>

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

class Function_T;

class CMyApp
{
public:
	CMyApp();
	~CMyApp();

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

	gCamera				m_camera;

	Texture2D			m_Texture_Red;
	Texture2D			m_Texture_Green;

	int TriangeCount;

	Geometry::TriMesh mesh;
	std::vector<glm::vec3> auto_normals;
	std::vector<glm::vec3> estimated_normals;

	struct Vertex
	{
		glm::vec3 p;
		glm::vec3 n;
		glm::vec2 t;
	};

	std::unique_ptr<Mesh> m_mesh;

	std::vector<Function_T> Functions;

	glm::vec3 AutoDiffNormal(Geometry::Vector3D v);

	void InitFunctions();
	void DrawCoordinateSystem();
	void InitShaders();
	void LoadMeshIntoBuffer(Geometry::TriMesh m);
	void DrawLine(glm::vec3 p_1, glm::vec3 p_2, float w);
	void DebugNormals();
	void DebugEdges();

	bool ImGUIDebugNormals = false;
	bool UseAutoNormals = false;

	struct MyTraits : public OpenMesh::DefaultTraits {
		using Point = OpenMesh::Vec3d; // the default would be Vec3f
		using Normal = OpenMesh::Vec3d;
		VertexTraits{
		  double mean;              // approximated mean curvature
		};
	};

	//typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> open_mesh;

};

