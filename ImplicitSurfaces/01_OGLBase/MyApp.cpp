#include "MyApp.h"

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
}

ImplicitSurfaceApp::~ImplicitSurfaceApp(void)
{
}

void ImplicitSurfaceApp::LoadMeshIntoBuffer(Geometry::TriMesh m) 
{

	std::vector<Vertex> vertices;
	if (vertex_trimap != nullptr) {
		delete[] vertex_trimap;
	}
	vertex_trimap = new std::vector<int>[m.points().size()];
	vertex_trimap->reserve(m.points().size());

	for (int i = 0; i < m.points().size(); ++i) {
	
		vertex_trimap[i] = std::vector<int>();

	}

	auto_normals.clear();
	auto_curvatures.clear();
	estimated_curvatures.clear();
	estimated_normals.clear();
	mesh_triangles.clear();

	TriangeCount = m.triangles().size() * 3;

	std::vector<int> indices(TriangeCount);

	int tricount = m.triangles().size();

	for (const Triangle & t : m.triangles()) {
	
		mesh_triangles.push_back(t);

	}

	for (int i = 0; i < mesh_triangles.size(); ++i)
	{

		Triangle t = mesh_triangles[i];
		indices.push_back(t[0]);
		indices.push_back(t[1]);
		indices.push_back(t[2]);
		vertex_trimap[t[0]].push_back(i);
		vertex_trimap[t[1]].push_back(i);
		vertex_trimap[t[2]].push_back(i);

	}
	
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

	for (int i = 0; i < m.points().size(); ++i) {
	
		Geometry::Vector3D v = m.points()[i];
		glm::vec3 n = (UseAutoNormals ? auto_normals[i] : estimated_normals[i]);
		float curv =  auto_curvatures[i];
		float curv_est = estimated_curvatures[i];
		float curv_norm = 0;
		float curv_norm_est = 0;
		if (UseNormalScale) {

			max_curv = 1;
			min_curv = -1;

			//0.9-et használ 1 helyett numerikus hibák elkerülésésre
			if (curv > 0.9) { curv = 0.9; }
			if (curv < -0.9) { curv = -0.9; }
			if (curv_est > 0.9) { curv_est = 0.9; }
			if (curv_est < -0.9) { curv_est = -0.9; }

		}
		if (curv >= 0) {
		
			curv_norm = 0.5 + (curv / (max_curv * 2));
			curv_norm_est = 0.5 + (curv_est / (max_curv * 2));
		}
		else {

			curv_norm = 0.5 - (curv / (min_curv * 2));
			curv_norm_est = 0.5 - (curv_est / (min_curv * 2));


		}
		if (curv_norm_est > 1 || curv_norm_est < 0) {
		
			std::cout << "ERROR CURVATURE NOT NORMALIZED " << curv_norm_est << std::endl;
		
		}

		vertices.push_back(Vertex{ glm::vec3(v[0],v[1],v[2]), n , glm::vec2(), curv_norm, curv_norm_est});

	}

	m_VertexBuffer.Clean();
	m_VertexBuffer.BufferData(vertices);

	m_Indices.Clean();
	m_Indices.BufferData(indices);

	m_Vao.Init(
		{
			// 0-ás attribútum "lényegében" glm::vec3-ak sorozata és az adatok az m_CubeVertexBuffer GPU pufferben vannak
			{ CreateAttribute<		0,						// attribútum: 0
									glm::vec3,				// CPU oldali adattípus amit a 0-ás attribútum meghatározására használtunk <- az eljárás a glm::vec3-ból kikövetkezteti, hogy 3 darab float-ból áll a 0-ás attribútum
									0,						// offset: az attribútum tároló elejétől vett offset-je, byte-ban
									sizeof(Vertex)			// stride: a következő csúcspont ezen attribútuma hány byte-ra van az aktuálistól
								>, m_VertexBuffer },
			{ CreateAttribute<1, glm::vec3, (sizeof(glm::vec3)), sizeof(Vertex)>, m_VertexBuffer },
			{ CreateAttribute<2, glm::vec2, (2 * sizeof(glm::vec3)), sizeof(Vertex)>, m_VertexBuffer },
			{ CreateAttribute<3, float, (2 * sizeof(glm::vec3)) + sizeof(glm::vec2), sizeof(Vertex)>, m_VertexBuffer },
			{ CreateAttribute<4, float, (2 * sizeof(glm::vec3)) + sizeof(glm::vec2) + sizeof(float), sizeof(Vertex)>, m_VertexBuffer },
		},
		m_Indices
	);
}

glm::vec3 ImplicitSurfaceApp::AutoDiffNormal(Geometry::Vector3D v) 
{

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

	std::vector<Triangle> tris = GetAllTriangles(vertex_index);
	glm::vec3 aggr_normal;
	float area = 0;
	for (Triangle t : tris) {
	
		area += GetTriangleArea(t);
		aggr_normal += GetTriangleNormal(t);

	}

	return (0.5f * area) * aggr_normal;

}

float ImplicitSurfaceApp::EstimateCurvature(int vertex_index)
{
	if(UseMeanCurvature) {
	
		std::vector<Triangle> tris = GetAllTriangles(vertex_index);
		tris = ReorientTriangles(tris, vertex_index);
		tris = OrderTriangles(tris);

		float numerator = 0;
		float area = 0;
		
		for (int i = 0; i < tris.size(); ++i) {

			if (IsValidTriangle(tris[i])) {

				glm::vec3 side = ConvertVector(mesh.points()[tris[i][0]] - mesh.points()[tris[i][1]]);
				glm::vec3 normal = GetTriangleNormal(tris[i]);
				glm::vec3 prev_normal = GetTriangleNormal(tris[i == 0 ? tris.size() - 1 : i - 1]);

				float angle;
				float f = dot(normal, prev_normal) / (length(normal) * length(prev_normal));

				if (f < 0.999 && f > -0.999) 
				{
					angle = acos(f);
				}
				else if (f >= 0.99) 
				{
					angle = 0;
				}
				else 
				{
					angle = M_PI;
				}
				float side_length = length(side);

				numerator += (angle * side_length);
				area += GetTriangleArea(tris[i]);

			}

		}

		return (3 * numerator) / (4 * area);
	
	}
	else {

		float area = 0;
		float angle = 2 * M_PI;

		std::vector<Triangle> tris = GetAllTriangles(vertex_index);
		tris = ReorientTriangles(tris, vertex_index);

		for (int i = 0; i < tris.size(); ++i){

			Triangle t = tris[i];
			area += GetTriangleArea(t);
			angle -= GetTriangleAngle(t);

		}

		float curv = angle * (3 / (area));

		return curv;

	}
}

std::vector<std::array<size_t, 3>> ImplicitSurfaceApp::OrderTriangles(std::vector<Triangle> tris)
{
	if (tris.size() > 1) {
		std::vector<Triangle> out_tris;
		Triangle& current_tri = tris[0];
		out_tris.push_back(current_tri);

		for (int i = 1; i < tris.size(); ++i) {

			for (const Triangle& t : tris) {

				if ((t[2] == current_tri[1] || t[1] == current_tri[1] || t[2] == current_tri[2] || t[1] == current_tri[2]) && !std::count(out_tris.begin(), out_tris.end(), t))
				{

					current_tri = t;
					out_tris.push_back(current_tri);
					break;

				}
			}
		}

		return out_tris;
	}
	else{

		return tris;

	}
}

std::vector<std::array<size_t, 3>> ImplicitSurfaceApp::ReorientTriangles(std::vector<Triangle> tris, int center_vertex_id)
{

	std::vector<Triangle> out_tris;

	for (int i = 0; i < tris.size(); ++i) {

		Triangle t = tris[i];
		for (int k = 0; k < 3; k++) {
			if (t[k] == center_vertex_id) {
				out_tris.push_back(Triangle{ t[k], t[(k+1) % 3], t[(k+2) % 3] });
			}
		}
	}

	return out_tris;
}

double ImplicitSurfaceApp::Det_2x2(double a1, double a2, double b1, double b2)
{
	return a1* b2 - a2 * b1;
}

void ImplicitSurfaceApp::InitFunctions()
{

	using namespace autodiff;

	//HD
	Functions.push_back(
		Function_T(
			[](Geometry::Vector3D p) {
				double x, y, z;
				x = p[0]; y = p[1]; z = p[2];
				return  x * x + y * y + z * z - 1;

			},
			[](dual x, dual y, dual z) -> autodiff::dual {


				return x * x + y * y + z * z - 1;

			},
				[](dual2nd x, dual2nd y, dual2nd z) -> autodiff::dual2nd {


				return x * x + y * y + z * z - 1;

			}, 2, 50
				)
	);

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

			}, 3, 50
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

			}, 4, 40
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

			}, 3, 50
				));

	Functions.push_back(
		Function_T(
			[](Geometry::Vector3D p) {
				double x, y, z;
				x = p[0]; y = p[1]; z = p[2];
				 return 1.0 * pow(x, 2) * pow(y, 2) * pow(z, 2) * (-30.400000000000006 + 16.0 * M_SQRT2) * pow(x - 1, 2) * pow(z - 1, 2) - pow(x, 2) * pow(y, 2) * pow(z, 2) * (-8.0 + 4 * M_SQRT2) * pow(x - 1, 2) - pow(x, 2) * pow(y, 2) * pow(z, 2) * (-8.0 + 4 * M_SQRT2) * pow(z - 1, 2) - 0.59999999999999964 * pow(x, 2) * pow(y, 2) * pow(x - 1, 2) * pow(z - 1, 2) + pow(x, 2) * pow(y, 2) * (y - 0.5) * pow(z - 1, 2) + pow(x, 2) * pow(z, 2) * pow(x - 1, 2) * (x - 0.5) + pow(x, 2) * pow(z, 2) * (-8.0 + 4 * M_SQRT2) * pow(x - 1, 2) * pow(z - 1, 2) + pow(x, 2) * pow(z, 2) * pow(z - 1, 2) * (z - 0.5) + pow(y, 2) * pow(z, 2) * pow(x - 1, 2) * (y - 0.5) - 0.59999999999999964 * pow(y, 2) * pow(z, 2) * pow(x - 1, 2) * pow(z - 1, 2) + (1.0 / 2.0) * pow(y, 2) * pow(x - 1, 2) * pow(z - 1, 2) * (x + 2 * y + z - 1.6000000000000001);
			},
			[](dual x, dual y, dual z) -> autodiff::dual {


				 return 1 * pow(x, 2) * pow(y, 2) * pow(z, 2) * (-30.400000000000006 + 16.0 * M_SQRT2) * pow(x - 1, 2) * pow(z - 1, 2) - pow(x, 2) * pow(y, 2) * pow(z, 2) * (-8.0 + 4 * M_SQRT2) * pow(x - 1, 2) - pow(x, 2) * pow(y, 2) * pow(z, 2) * (-8.0 + 4 * M_SQRT2) * pow(z - 1, 2) - 0.59999999999999964 * pow(x, 2) * pow(y, 2) * pow(x - 1, 2) * pow(z - 1, 2) + pow(x, 2) * pow(y, 2) * (y - 0.5) * pow(z - 1, 2) + pow(x, 2) * pow(z, 2) * pow(x - 1, 2) * (x - 0.5) + pow(x, 2) * pow(z, 2) * (-8.0 + 4 * M_SQRT2) * pow(x - 1, 2) * pow(z - 1, 2) + pow(x, 2) * pow(z, 2) * pow(z - 1, 2) * (z - 0.5) + pow(y, 2) * pow(z, 2) * pow(x - 1, 2) * (y - 0.5) - 0.59999999999999964 * pow(y, 2) * pow(z, 2) * pow(x - 1, 2) * pow(z - 1, 2) + (1.0 / 2.0) * pow(y, 2) * pow(x - 1, 2) * pow(z - 1, 2) * (x + 2 * y + z - 1.6000000000000001);
			},
				[](dual2nd x, dual2nd y, dual2nd z) -> autodiff::dual2nd {


				 return 1 * pow(x, 2) * pow(y, 2) * pow(z, 2) * (-30.400000000000006 + 16.0 * M_SQRT2) * pow(x - 1, 2) * pow(z - 1, 2) - pow(x, 2) * pow(y, 2) * pow(z, 2) * (-8.0 + 4 * M_SQRT2) * pow(x - 1, 2) - pow(x, 2) * pow(y, 2) * pow(z, 2) * (-8.0 + 4 * M_SQRT2) * pow(z - 1, 2) - 0.59999999999999964 * pow(x, 2) * pow(y, 2) * pow(x - 1, 2) * pow(z - 1, 2) + pow(x, 2) * pow(y, 2) * (y - 0.5) * pow(z - 1, 2) + pow(x, 2) * pow(z, 2) * pow(x - 1, 2) * (x - 0.5) + pow(x, 2) * pow(z, 2) * (-8.0 + 4 * M_SQRT2) * pow(x - 1, 2) * pow(z - 1, 2) + pow(x, 2) * pow(z, 2) * pow(z - 1, 2) * (z - 0.5) + pow(y, 2) * pow(z, 2) * pow(x - 1, 2) * (y - 0.5) - 0.59999999999999964 * pow(y, 2) * pow(z, 2) * pow(x - 1, 2) * pow(z - 1, 2) + (1.0 / 2.0) * pow(y, 2) * pow(x - 1, 2) * pow(z - 1, 2) * (x + 2 * y + z - 1.6000000000000001);
			}, 1, 50
				));

	//return 1.0*pow(x, 2)*pow(y, 2)*pow(z, 2)*(-30.400000000000006 + 16.0*M_SQRT2)*pow(x - 1, 2)*pow(z - 1, 2) - pow(x, 2)*pow(y, 2)*pow(z, 2)*(-8.0 + 4*M_SQRT2)*pow(x - 1, 2) - pow(x, 2)*pow(y, 2)*pow(z, 2)*(-8.0 + 4*M_SQRT2)*pow(z - 1, 2) - 0.59999999999999964*pow(x, 2)*pow(y, 2)*pow(x - 1, 2)*pow(z - 1, 2) + pow(x, 2)*pow(y, 2)*(y - 0.5)*pow(z - 1, 2) + pow(x, 2)*pow(z, 2)*pow(x - 1, 2)*(x - 0.5) + pow(x, 2)*pow(z, 2)*(-8.0 + 4*M_SQRT2)*pow(x - 1, 2)*pow(z - 1, 2) + pow(x, 2)*pow(z, 2)*pow(z - 1, 2)*(z - 0.5) + pow(y, 2)*pow(z, 2)*pow(x - 1, 2)*(y - 0.5) - 0.59999999999999964*pow(y, 2)*pow(z, 2)*pow(x - 1, 2)*pow(z - 1, 2) + (1.0/2.0)*pow(y, 2)*pow(x - 1, 2)*pow(z - 1, 2)*(x + 2*y + z - 1.6000000000000001);

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
	int acc = Functions[function_index].accuracy;
	mesh = IMC::marching_cubes(Functions[function_index].func, 0, Functions[function_index].boundingBox, { acc, acc, acc }, true);
	//delete[] vertex_trimap;
	LoadMeshIntoBuffer(mesh);

}

std::vector<std::array<size_t, 3>> ImplicitSurfaceApp::GetAllTriangles(int vertex_id)
{

	std::vector<std::array<size_t, 3>> out;
	for (int i : vertex_trimap[vertex_id]) {
		if(IsValidTriangle(mesh_triangles[i]))
		out.push_back(mesh_triangles[i]);

	}
	return out;

}

float ImplicitSurfaceApp::GetTriangleArea(std::array<size_t, 3> tri)
{
	{
		std::vector<glm::vec3> pos;

		for (size_t i : tri) {

			pos.push_back(ConvertVector(mesh.points()[i]));

		}

		return 0.5f * length((cross(pos[0] - pos[1], pos[1] - pos[2])));

	}
}

float ImplicitSurfaceApp::GetTriangleAngle(Triangle t)
{

	glm::vec3 side1 = ConvertVector(mesh.points()[t[0]]) - ConvertVector(mesh.points()[t[1]]);
	glm::vec3 side2 = ConvertVector(mesh.points()[t[0]]) - ConvertVector(mesh.points()[t[2]]);

	float sides = (glm::length(side1) * glm::length(side2));

	float x = (glm::dot(side1, side2) / sides);

	if (x > 1.f) {
		return 0;
	}
	if(x < -1){
		return M_PI;
	}

	return glm::acos(x);
}

glm::vec3 ImplicitSurfaceApp::GetTriangleNormal(Triangle t)
{
	return cross(ConvertVector(mesh.points()[t[0]] - mesh.points()[t[1]]), ConvertVector(mesh.points()[t[0]] - mesh.points()[t[2]]));
}

bool ImplicitSurfaceApp::IsValidTriangle(Triangle t)
{
	return t[0] != t[1] && t[0] != t[2] && t[1] != t[2];
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
		{ 4, "vs_in_curv_est" },		// VAO 4-es csatorna menjen a vs_in_curv_est-be
	});

	m_program.LinkProgram();

}

bool ImplicitSurfaceApp::Init()
{
	// törlési szín legyen kékes
	glClearColor(0.125f, 0.25f, 0.5f, 1.0f);
	glEnable(GL_DEPTH_TEST); // mélységi teszt bekapcsolása (takarás)

	InitShaders();
	InitFunctions();
	
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
	m_program.SetUniform("MVP", viewProj * MeshWorld);
	m_program.SetUniform("world", MeshWorld);
	m_program.SetUniform("showCurv", VisualizeCurvature ? 1 : 0);
	m_program.SetUniform("useEstimatedCurv", UseAutoCurvature ? 0 : 1);
	m_program.SetUniform("worldIT", glm::inverse(glm::transpose(MeshWorld)));
	m_program.SetUniform("viewPos", m_camera.GetEye());

	m_Vao.Bind();

	glm::mat4 cubeWorld;


	{
		glDrawElements(GL_TRIANGLES, TriangeCount * 3, GL_UNSIGNED_INT, nullptr);
	}
	if (ImguiDebugNormals) {

		//DebugNormals();

	}
	if (ImguiDebugEdges) {

		//DebugEdges();

	}
	m_program.Unuse();



	ImGui::ShowTestWindow();
	//ImGui::Checkbox("Debug Normals",&ImguiDebugNormals);
	//ImGui::Checkbox("Debug Edges", &ImguiDebugEdges);
	ImGui::Checkbox("Visualize Curvature", &VisualizeCurvature);

	if (ImGui::Checkbox("Use Differentiated Curvature", &UseAutoCurvature)) {
	
		//LoadMeshIntoBuffer(mesh);

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
