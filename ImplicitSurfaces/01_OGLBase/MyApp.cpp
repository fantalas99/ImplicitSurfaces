﻿#include "MyApp.h"

#include <math.h>
#include <vector>
#include <autodiff/forward/dual.hpp>
#include "Function.h"

#include <array>
#include <list>
#include <tuple>
#include <imgui/imgui.h>
#include "includes/GLUtils.hpp"

CMyApp::CMyApp(void)
{
	m_camera.SetView(glm::vec3(5, 5, 5), glm::vec3(0, 0, 0), glm::vec3(0, 1, 0));
	m_mesh = nullptr;
}

CMyApp::~CMyApp(void)
{
}

void CMyApp::LoadMeshIntoBuffer(Geometry::TriMesh m)
{

	std::vector<Vertex>vertices;

	TriangeCount = m.triangles().size() * 3;

	std::vector<int> indices(TriangeCount);

	using Triangle = std::array<size_t, 3>;

	using open_mesh = OpenMesh::TriMesh_ArrayKernelT<MyTraits>;

	open_mesh omesh;

	std::vector<open_mesh::VertexHandle> openmesh_vertices, openmesh_triangles;

	for (int i = 0; i < m.points().size(); ++i) {

		Geometry::Vector3D v = m.points()[i];

		glm::vec3 n = AutoDiffNormal(v);

		openmesh_vertices.push_back(omesh.add_vertex(MyTraits::Point(v[0], v[1], v[2])));

		auto_normals.push_back(n);



	}


	for (Triangle t : m.triangles())
	{
		indices.push_back(t[0]);
		indices.push_back(t[1]);
		indices.push_back(t[2]);
		openmesh_triangles.push_back(openmesh_vertices[t[0]]);
		openmesh_triangles.push_back(openmesh_vertices[t[1]]);
		openmesh_triangles.push_back(openmesh_vertices[t[2]]);
		omesh.add_face(openmesh_triangles);

	}
	for (open_mesh::VertexHandle v : openmesh_vertices) {

		MyTraits::Normal est_n;
		omesh.calc_vertex_normal_correct(v, est_n);
		estimated_normals.push_back(glm::vec3(est_n[0], est_n[1], est_n[2]));

	}

	for (int i = 0; i < m.points().size(); ++i) {
	
		Geometry::Vector3D v = m.points()[i];
		glm::vec3 n = (UseAutoNormals ? auto_normals[i] : estimated_normals[i]);
		vertices.push_back(Vertex{ glm::vec3(v[0],v[1],v[2]), n , glm::vec2()});

	
	}

	

	m_CubeVertexBuffer.BufferData(vertices);

	// és a primitíveket alkotó csúcspontok indexei (az előző tömbökből) - triangle list-el való kirajzolásra felkészülve
	m_CubeIndices.BufferData(indices);

	// geometria VAO-ban való regisztrálása
	m_CubeVao.Init(
		{
			// 0-ás attribútum "lényegében" glm::vec3-ak sorozata és az adatok az m_CubeVertexBuffer GPU pufferben vannak
			{ CreateAttribute<		0,						// attribútum: 0
									glm::vec3,				// CPU oldali adattípus amit a 0-ás attribútum meghatározására használtunk <- az eljárás a glm::vec3-ból kikövetkezteti, hogy 3 darab float-ból áll a 0-ás attribútum
									0,						// offset: az attribútum tároló elejétől vett offset-je, byte-ban
									sizeof(Vertex)			// stride: a következő csúcspont ezen attribútuma hány byte-ra van az aktuálistól
								>, m_CubeVertexBuffer },
			{ CreateAttribute<1, glm::vec3, (sizeof(glm::vec3)), sizeof(Vertex)>, m_CubeVertexBuffer },
			{ CreateAttribute<2, glm::vec2, (2 * sizeof(glm::vec3)), sizeof(Vertex)>, m_CubeVertexBuffer },
		},
		m_CubeIndices
	);
}

glm::vec3 CMyApp::AutoDiffNormal(Geometry::Vector3D v)
{
	using namespace autodiff;

	dual x = v[0];
	dual y = v[1];
	dual z = v[2];

	double dx = derivative(Functions[0].func_dual, wrt(x), at(x, y, z));
	double dy = derivative(Functions[0].func_dual, wrt(y), at(x, y, z));
	double dz = derivative(Functions[0].func_dual, wrt(z), at(x, y, z));

	//std::cout << dx << " | " << dy << " | " << dz << std::endl;

	glm::vec3 normal(dx, dy, dz);
	normal = glm::normalize(normal);

	return normal;

}

void CMyApp::InitFunctions()
{

	//std::function<autodiff::dual(Geometry::Vector3D)> d = 
	using namespace autodiff;
	Functions.push_back(
		Function_T(
			[](Geometry::Vector3D p) {
				double x, y, z;
				x = p[0]; y = p[1]; z = p[2];
				return  x * x + y * y + z * z - 1;

			},
			[](dual x, dual y, dual z) -> autodiff::dual {


				return x * x + y * y + z * z - 1;

			}
		));

}

void CMyApp::DebugNormals() {

	for (int i = 0; i < mesh.points().size(); ++i) {
	
		Geometry::Vector3D vec = mesh.points()[i];
		Vertex vertex;
		vertex.p = glm::vec3(vec[0], vec[1], vec[2]);
		DrawLine(vertex.p, vertex.p + 0.1f * auto_normals[i], 5);
	
	}

}

void CMyApp::DebugEdges() {


}

void CMyApp::DrawLine(glm::vec3 p_1, glm::vec3 p_2, float w = 2) {

	//glLineWidth(w);
	//glBegin(GL_LINES);
	glVertex3f(p_1.x, p_1.y, p_1.z);
	glVertex3f(p_2.x, p_2.y, p_2.z);
	//glEnd();
}

void CMyApp::DrawCoordinateSystem()
{
	// draw some lines
	glColor3f(1.0, 0.0, 0.0); // red x
	glBegin(GL_LINES);
	// x aix

	glVertex3f(-4.0, 0.0f, 0.0f);
	glVertex3f(4.0, 0.0f, 0.0f);

	glVertex3f(4.0, 0.0f, 0.0f);
	glVertex3f(3.0, 1.0f, 0.0f);

	glVertex3f(4.0, 0.0f, 0.0f);
	glVertex3f(3.0, -1.0f, 0.0f);
	glEnd();

	// y 
	glColor3f(0.0, 1.0, 0.0); // green y
	glBegin(GL_LINES);
	glVertex3f(0.0, -4.0f, 0.0f);
	glVertex3f(0.0, 4.0f, 0.0f);

	glVertex3f(0.0, 4.0f, 0.0f);
	glVertex3f(1.0, 3.0f, 0.0f);

	glVertex3f(0.0, 4.0f, 0.0f);
	glVertex3f(-1.0, 3.0f, 0.0f);
	glEnd();

	// z 
	glColor3f(0.0, 0.0, 1.0); // blue z
	glBegin(GL_LINES);
	glVertex3f(0.0, 0.0f, -4.0f);
	glVertex3f(0.0, 0.0f, 4.0f);


	glVertex3f(0.0, 0.0f, 4.0f);
	glVertex3f(0.0, 1.0f, 3.0f);

	glVertex3f(0.0, 0.0f, 4.0f);
	glVertex3f(0.0, -1.0f, 3.0f);
	glEnd();
}

void CMyApp::InitShaders()
{
	// a shadereket tároló program létrehozása az OpenGL-hez hasonló módon:
	m_program.AttachShaders({
		{ GL_VERTEX_SHADER, "myVert.vert"},
		{ GL_FRAGMENT_SHADER, "myFrag.frag"}
	});

	// attributomok osszerendelese a VAO es shader kozt
	m_program.BindAttribLocations({
		{ 0, "vs_in_pos" },				// VAO 0-as csatorna menjen a vs_in_pos-ba
		{ 1, "vs_in_norm" },			// VAO 1-es csatorna menjen a vs_in_norm-ba
		{ 2, "vs_in_tex" },				// VAO 2-es csatorna menjen a vs_in_tex-be
	});

	m_program.LinkProgram();

	// shader program rövid létrehozása, egyetlen függvényhívással a fenti három:
	m_programSkybox.Init(
		{
			{ GL_VERTEX_SHADER, "skybox.vert" },
			{ GL_FRAGMENT_SHADER, "skybox.frag" }
		},
		{
			{ 0, "vs_in_pos" },				// VAO 0-as csatorna menjen a vs_in_pos-ba
		}
	);
}

bool CMyApp::Init()
{
	// törlési szín legyen kékes
	glClearColor(0.125f, 0.25f, 0.5f, 1.0f);

	//glEnable(GL_CULL_FACE); 
	// kapcsoljuk be a hatrafele nezo lapok eldobasat
	glEnable(GL_DEPTH_TEST); // mélységi teszt bekapcsolása (takarás)

	InitShaders();
	InitFunctions();

	// egyéb  betöltése
	m_Texture_Red.FromFile("assets/red.jpg");
	m_Texture_Green.FromFile("assets/green.jpg");

	//m_mesh = std::unique_ptr<Mesh>(ObjParser::parse("assets/Suzanne.obj"));
	//m_mesh->initBuffers();
	
	m_camera.SetProj(glm::radians(60.0f), 640.0f / 480.0f, 0.01f, 1000.0f);

	std::function<double(Geometry::Vector3D)> sphere_surface = [](Geometry::Vector3D p) {return sin(p[0] * p[2]) - p[1]; };

	//auto bounding box calculation
	std::array<Geometry::Vector3D, 2> boundingBox = { Geometry::Vector3D(-5,-5,-5), Geometry::Vector3D(5,5,5) };

	mesh = IMC::marching_cubes(Functions[0].func , 0, boundingBox, { 100, 100, 100 }, true);

	LoadMeshIntoBuffer(mesh);

	

	return true;
}

void CMyApp::Clean()
{
}

void CMyApp::Update()
{
	static Uint32 last_time = SDL_GetTicks();
	float delta_time = (SDL_GetTicks() - last_time) / 1000.0f;

	m_camera.Update(delta_time);

	last_time = SDL_GetTicks();
}

void CMyApp::Render()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glm::mat4 viewProj = m_camera.GetViewProj();
	glm::mat4 MeshWorld = glm::mat4(1.0f);

	//DrawCoordinateSystem();

	m_program.Use();
	m_program.SetTexture("texImage", 0, m_Texture_Red);
	m_program.SetUniform("MVP", viewProj * MeshWorld);
	m_program.SetUniform("world", MeshWorld);
	m_program.SetUniform("worldIT", glm::inverse(glm::transpose(MeshWorld)));
	m_program.SetUniform("viewPos", m_camera.GetEye());
	//m_mesh->draw();

	m_CubeVao.Bind();
	glm::mat4 cubeWorld;

	
	m_program.SetTexture("texImage", 0, m_Texture_Green);

	{
		glDrawElements(GL_TRIANGLES, TriangeCount * 3, GL_UNSIGNED_INT, nullptr);
	}

	if (ImGUIDebugNormals) {

		glLineWidth(5);
		glBegin(GL_LINES);
		DebugNormals();
		glEnd();

	}

	m_program.Unuse();


	ImGui::ShowTestWindow();
	ImGui::Checkbox("Show Normals",&ImGUIDebugNormals);

}

void CMyApp::KeyboardDown(SDL_KeyboardEvent& key)
{
	m_camera.KeyboardDown(key);
}

void CMyApp::KeyboardUp(SDL_KeyboardEvent& key)
{
	m_camera.KeyboardUp(key);
}

void CMyApp::MouseMove(SDL_MouseMotionEvent& mouse)
{
	m_camera.MouseMove(mouse);
}

void CMyApp::MouseDown(SDL_MouseButtonEvent& mouse)
{
}

void CMyApp::MouseUp(SDL_MouseButtonEvent& mouse)
{
}

void CMyApp::MouseWheel(SDL_MouseWheelEvent& wheel)
{
}

void CMyApp::Resize(int _w, int _h)
{
	glViewport(0, 0, _w, _h );

	m_camera.Resize(_w, _h);
}
