﻿#include "MyApp.h"

#include <math.h>
#include <vector>
#include <autodiff/forward/dual.hpp>
#include "Function.h"
#include "Dependencies/TriMesh/trimesh.cpp"

#include <array>
#include <list>
#include <tuple>
#include <imgui/imgui.h>
#include "includes/GLUtils.hpp"


ImplicitSurfaceApp::ImplicitSurfaceApp(void)
{
	m_camera.SetView(glm::vec3(5, 5, 5), glm::vec3(0, 0, 0), glm::vec3(0, 1, 0));
	m_mesh = nullptr;
}

ImplicitSurfaceApp::~ImplicitSurfaceApp(void)
{
}

void ImplicitSurfaceApp::LoadMeshIntoBuffer(Geometry::TriMesh m) 
{

	std::vector<Vertex> vertices;

	auto_normals.clear();
	auto_curvatures.clear();

	TriangeCount = m.triangles().size() * 3;

	std::vector<int> indices(TriangeCount);

	using Triangle = std::array<size_t, 3>;

	int tricount = m.triangles().size();
	std::vector<trimesh::triangle_t> tris;
	std::vector<trimesh::edge_t> edges;

	for (Triangle t : m.triangles())
	{
		indices.push_back(t[0]);
		indices.push_back(t[1]);
		indices.push_back(t[2]);
		trimesh::triangle_t tri;
		tri.v[0] = t[0];
		tri.v[1] = t[1];
		tri.v[2] = t[2];
		tris.push_back(tri);
	
	}

	trimesh::unordered_edges_from_triangles(tricount, &tris[0], edges);

	tmesh.build(m.points().size(), tricount, &tris[0], edges.size(), &edges[0]);

	for (int i = 0; i < m.points().size(); ++i) {

		Geometry::Vector3D v = m.points()[i];

		glm::vec3 n = AutoDiffNormal(v);

		auto_normals.push_back(n);
		estimated_normals.push_back(EstimateNormal(i));
		auto_curvatures.push_back(AutoDiffCurvature(v));
		estimated_curvatures.push_back(EstimateCurvature(i));

	}


	float max_curv = -100000;
	float min_curv = 100000;

	for (float f : auto_curvatures) {

		if (f > max_curv) {

			max_curv = f;

		}
		if (f < min_curv) {
		
			min_curv = f;

		}

	}

	for (int i = 0; i < m.points().size(); ++i) {
	
		Geometry::Vector3D v = m.points()[i];
		glm::vec3 n = (UseAutoNormals ? auto_normals[i] : estimated_normals[i]);
		float curv = UseAutoCurvature ? auto_curvatures[i] : estimated_curvatures[i];
		float curv_norm = 0;
		if (UseNormalScale) {

			max_curv = 1;
			min_curv = -1;
			if (curv > 1) { curv = 1; }
			if (curv < -1) { curv = -1; }
		}
		if (curv >= 0) {
		
			curv_norm = 0.5 + (curv / (max_curv * 2));

		}
		else {

			curv_norm = 0.5 - (curv / (min_curv * 2));

		}

		vertices.push_back(Vertex{ glm::vec3(v[0],v[1],v[2]), n , glm::vec2(), curv_norm});

	}

	m_CubeVertexBuffer.Clean();
	m_CubeVertexBuffer.BufferData(vertices);

	m_CubeIndices.Clean();
	m_CubeIndices.BufferData(indices);

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
			{ CreateAttribute<3, float, (2 * sizeof(glm::vec3)) + sizeof(glm::vec2), sizeof(Vertex)>, m_CubeVertexBuffer },

		},
		m_CubeIndices
	);
}

glm::vec3 ImplicitSurfaceApp::AutoDiffNormal(Geometry::Vector3D v) {

	using namespace autodiff;

	dual x = v[0];
	dual y = v[1];
	dual z = v[2];

	double dx = derivative(Functions[0].func_dual, wrt(x), at(x, y, z));
	double dy = derivative(Functions[0].func_dual, wrt(y), at(x, y, z));
	double dz = derivative(Functions[0].func_dual, wrt(z), at(x, y, z));

	glm::vec3 normal(dx, dy, dz);
	normal = glm::normalize(normal);

	return normal;

}

float ImplicitSurfaceApp::AutoDiffCurvature(Geometry::Vector3D v)
{

	using namespace autodiff;
	
	dual2nd x = v[0];
	dual2nd y = v[1];
	dual2nd z = v[2];
	double dx = derivative(Functions[ActiveFunctionIndex].func_dual_2nd, wrt(x), at(x, y, z));
	double dy = derivative(Functions[ActiveFunctionIndex].func_dual_2nd, wrt(y), at(x, y, z));
	double dz = derivative(Functions[ActiveFunctionIndex].func_dual_2nd, wrt(z), at(x, y, z));

	glm::vec3 grad = glm::vec3(dx, dy, dz);
	float grad_length = glm::length(grad);

	glm::mat3x3 hess;

	hess[0][0] = derivative<2>(Functions[ActiveFunctionIndex].func_dual_2nd, wrt(x, x), at(x, y, z));
	hess[0][1] = derivative<2>(Functions[ActiveFunctionIndex].func_dual_2nd, wrt(x, y), at(x, y, z));
	hess[0][2] = derivative<2>(Functions[ActiveFunctionIndex].func_dual_2nd, wrt(x, z), at(x, y, z));
	hess[1][0] = derivative<2>(Functions[ActiveFunctionIndex].func_dual_2nd, wrt(y, x), at(x, y, z));
	hess[1][1] = derivative<2>(Functions[ActiveFunctionIndex].func_dual_2nd, wrt(y, y), at(x, y, z));
	hess[1][2] = derivative<2>(Functions[ActiveFunctionIndex].func_dual_2nd, wrt(y, z), at(x, y, z));
	hess[2][0] = derivative<2>(Functions[ActiveFunctionIndex].func_dual_2nd, wrt(z, x), at(x, y, z));
	hess[2][1] = derivative<2>(Functions[ActiveFunctionIndex].func_dual_2nd, wrt(z, y), at(x, y, z));
	hess[2][2] = derivative<2>(Functions[ActiveFunctionIndex].func_dual_2nd, wrt(z, z), at(x, y, z));

	if (UseMeanCurvature) {

	
		float trace = hess[0][0] + hess[1][1] + hess[2][2];

		glm::vec3 hess_x_grad = hess * grad;
		float gradt_x_hess_x_grad = dx * hess_x_grad[0] + dy * hess_x_grad[1] + dz * hess_x_grad[2];

		return (grad_length * grad_length * trace - gradt_x_hess_x_grad) / (grad_length * grad_length * grad_length * 2);

	}
	else {
	
	
		glm::mat3x3 adj;

		adj[0][0] = Det_2x2(hess[1][1],hess[1][2],hess[2][1],hess[2][2]);
		adj[0][1] = -1 * Det_2x2(hess[0][1], hess[0][2], hess[2][1], hess[2][2]);
		adj[0][2] = Det_2x2(hess[0][1], hess[0][2], hess[1][1], hess[1][2]);
		adj[1][0] = -1 * Det_2x2(hess[1][0], hess[1][2], hess[2][0], hess[2][2]);
		adj[1][1] = Det_2x2(hess[0][0], hess[0][2], hess[2][0], hess[2][2]);
		adj[1][2] = -1 * Det_2x2(hess[0][0], hess[0][2], hess[1][0], hess[1][2]);
		adj[2][0] = Det_2x2(hess[1][0], hess[1][1], hess[2][0], hess[2][1]);
		adj[2][1] = -1 * Det_2x2(hess[0][0], hess[0][1], hess[2][0], hess[2][1]);
		adj[2][2] = Det_2x2(hess[0][0], hess[0][1], hess[1][0], hess[1][1]);

		glm::vec3 adj_x_grad = adj * grad;
		float gradt_x_adj_x_grad = dx * adj_x_grad[0] + dy * adj_x_grad[1] + dz * adj_x_grad[2];

		return (gradt_x_adj_x_grad) / (pow(grad_length, 4));
	}
}

glm::vec3 ImplicitSurfaceApp::EstimateNormal(int vertex_index)
{

	std::vector<trimesh::index_t> neighbours = tmesh.vertex_vertex_neighbors(vertex_index);
	
	return glm::vec3(0, 0, 1);

}

float ImplicitSurfaceApp::EstimateCurvature(int vertex_index)
{
	return 0.0f;
}

double ImplicitSurfaceApp::Det_2x2(double a1, double a2, double b1, double b2)
{
	return a1* b2 - a2 * b1;
}

void ImplicitSurfaceApp::InitFunctions()
{

	//std::function<autodiff::dual(Geometry::Vector3D)> d = 
	using namespace autodiff;
	Functions.push_back(
		Function_T(
			[](Geometry::Vector3D p) {
				double x, y, z;
				x = p[0]; y = p[1]; z = p[2];
				return  x * x + y * y + z - 1;

			},
			[](dual x, dual y, dual z) -> autodiff::dual {

				
				return x * x + y * y + z - 1;

			},
			[](dual2nd x, dual2nd y, dual2nd z) -> autodiff::dual2nd {


				return x * x + y * y + z - 1;

			}
		));
	Functions.push_back(
		Function_T(
			[](Geometry::Vector3D p) {
				double x, y, z;
				x = p[0]; y = p[1]; z = p[2];
				return  sin(x * z) + y - 1;

			},
			[](dual x, dual y, dual z) -> autodiff::dual {


				return sin(x * z) + y - 1;

			},
				[](dual2nd x, dual2nd y, dual2nd z) -> autodiff::dual2nd {


				return sin(x * z) + y - 1;

			}
			));
	Functions.push_back(
		Function_T(
			[](Geometry::Vector3D p) {
		double x, y, z;
		x = p[0]; y = p[1]; z = p[2];
		return 12.5990051961515 * pow(x, 2) * pow(y, 2) * pow(z, 2) + 10.0000000000000 * pow(x, 2) * pow(y, 2) + 2.34314575050762 * pow(x, 2) * pow(z, 2) + 10.0000000000000 * pow(y, 2) * pow(z, 2) + (x - 0.500000000000000) * pow(x, 2) + (2 * y - 0.400000000000000) * pow(y, 2) + (z - 0.500000000000000) * pow(z, 2);

	},
		[](dual x, dual y, dual z) -> autodiff::dual {


		return 12.5990051961515 * pow(x, 2) * pow(y, 2) * pow(z, 2) + 10.0000000000000 * pow(x, 2) * pow(y, 2) + 2.34314575050762 * pow(x, 2) * pow(z, 2) + 10.0000000000000 * pow(y, 2) * pow(z, 2) + (x - 0.500000000000000) * pow(x, 2) + (2 * y - 0.400000000000000) * pow(y, 2) + (z - 0.500000000000000) * pow(z, 2);

	},
		[](dual2nd x, dual2nd y, dual2nd z) -> autodiff::dual2nd {


		return 12.5990051961515 * pow(x, 2) * pow(y, 2) * pow(z, 2) + 10.0000000000000 * pow(x, 2) * pow(y, 2) + 2.34314575050762 * pow(x, 2) * pow(z, 2) + 10.0000000000000 * pow(y, 2) * pow(z, 2) + (x - 0.500000000000000) * pow(x, 2) + (2 * y - 0.400000000000000) * pow(y, 2) + (z - 0.500000000000000) * pow(z, 2);

	}
	));

}

void ImplicitSurfaceApp::DebugNormals() {
	
	glBegin(GL_LINES);

	for (int i = 0; i < mesh.points().size(); ++i) {
	
		Geometry::Vector3D vec = mesh.points()[i];
		
		glm::vec3 v(vec[0], vec[1], vec[2]);
		glm::vec3 v_n = v + (auto_normals[i] * 0.1f);
		glVertex3f(v[0], v[1], v[2]);
		glVertex3f(v_n[0], v_n[1], v_n[2]);

	}

	glEnd();

}

void ImplicitSurfaceApp::DebugEdges() {

	glBegin(GL_LINES);

	for (const Geometry::TriMesh::Triangle & t : mesh.triangles()) {
	
		Geometry::Vector3D v_1 = mesh.points()[t[0]];
		Geometry::Vector3D v_2 = mesh.points()[t[1]];
		Geometry::Vector3D v_3 = mesh.points()[t[2]];

		glVertex3f(v_1[0], v_1[1], v_1[2]);
		glVertex3f(v_2[0], v_2[1], v_2[2]);

		glVertex3f(v_2[0], v_2[1], v_2[2]);
		glVertex3f(v_3[0], v_3[1], v_3[2]);

		glVertex3f(v_3[0], v_3[1], v_3[2]);
		glVertex3f(v_1[0], v_1[1], v_1[2]);
	}

	glEnd();
}

void ImplicitSurfaceApp::GenerateMesh(int function_index)
{

	std::array<Geometry::Vector3D, 2> boundingBox = { Geometry::Vector3D(-5,-5,-5), Geometry::Vector3D(5,5,5) };

	mesh = IMC::marching_cubes(Functions[function_index].func, 0, boundingBox, { 100, 100, 100 }, true);

	LoadMeshIntoBuffer(mesh);

}

void ImplicitSurfaceApp::DrawLine(glm::vec3 p_1, glm::vec3 p_2, float w = 2) {

	//glLineWidth(w);
	//glBegin(GL_LINES);
	glVertex3f(p_1.x, p_1.y, p_1.z);
	glVertex3f(p_2.x, p_2.y, p_2.z);
	//glEnd();
}

void ImplicitSurfaceApp::DrawCoordinateSystem()
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

void ImplicitSurfaceApp::InitShaders()
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
		{ 3, "vs_in_curv" },			// VAO 3-es csatorna menjen a vs_in_curv-be
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

bool ImplicitSurfaceApp::Init()
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

	GenerateMesh(0);

	return true;
}

void ImplicitSurfaceApp::Clean()
{
}

void ImplicitSurfaceApp::Update()
{
	static Uint32 last_time = SDL_GetTicks();
	float delta_time = (SDL_GetTicks() - last_time) / 1000.0f;

	m_camera.Update(delta_time);

	last_time = SDL_GetTicks();
}

void ImplicitSurfaceApp::Render()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glm::mat4 viewProj = m_camera.GetViewProj();
	glm::mat4 MeshWorld = glm::mat4(1.0f);

	m_program.Use();
	m_program.SetTexture("texImage", 0, m_Texture_Red);
	m_program.SetUniform("MVP", viewProj * MeshWorld);
	m_program.SetUniform("world", MeshWorld);
	m_program.SetUniform("showCurv", VisualizeCurvature ? 1 : 0);
	m_program.SetUniform("worldIT", glm::inverse(glm::transpose(MeshWorld)));
	m_program.SetUniform("viewPos", m_camera.GetEye());

	m_CubeVao.Bind();

	glm::mat4 cubeWorld;

	
	m_program.SetTexture("texImage", 0, m_Texture_Green);

	{
		glDrawElements(GL_TRIANGLES, TriangeCount * 3, GL_UNSIGNED_INT, nullptr);
	}
	if (ImguiDebugNormals) {

		DebugNormals();

	}
	if (ImguiDebugEdges) {

		DebugEdges();

	}
	m_program.Unuse();



	ImGui::ShowTestWindow();
	//ImGui::Checkbox("Debug Normals",&ImguiDebugNormals);
	//ImGui::Checkbox("Debug Edges", &ImguiDebugEdges);
	ImGui::Checkbox("Visualize Curvature", &VisualizeCurvature);

	if (ImGui::Checkbox("Use Estimated Curvature", &UseAutoCurvature)) {
	
		LoadMeshIntoBuffer(mesh);

	}

	if (ImGui::Checkbox("Use Mean Curvature", &UseMeanCurvature)) {
	
		LoadMeshIntoBuffer(mesh);
	
	}
	if (ImGui::ListBox("Surfaces", &ActiveFunctionIndex, FunctionNames)) {

		std::cout << ActiveFunctionIndex << std::endl;
		GenerateMesh(ActiveFunctionIndex);
	
	}
	
}

void ImplicitSurfaceApp::KeyboardDown(SDL_KeyboardEvent& key)
{
	m_camera.KeyboardDown(key);
}

void ImplicitSurfaceApp::KeyboardUp(SDL_KeyboardEvent& key)
{
	m_camera.KeyboardUp(key);
}

void ImplicitSurfaceApp::MouseMove(SDL_MouseMotionEvent& mouse)
{
	m_camera.MouseMove(mouse);
}

void ImplicitSurfaceApp::MouseDown(SDL_MouseButtonEvent& mouse)
{
}

void ImplicitSurfaceApp::MouseUp(SDL_MouseButtonEvent& mouse)
{
}

void ImplicitSurfaceApp::MouseWheel(SDL_MouseWheelEvent& wheel)
{
}

void ImplicitSurfaceApp::Resize(int _w, int _h)
{
	glViewport(0, 0, _w, _h );

	m_camera.Resize(_w, _h);
}
