#pragma once
#include <vector>
#include <glm/glm.hpp>
#include <iostream>

class HalfEdgeMesh
{

public:

	struct half_edge {

		int vertex_id;
		int opposite_edge_id;
		int next_edge_id;


		half_edge(int vtx_id, int op_id, int nx_id) : vertex_id(vtx_id), opposite_edge_id(op_id), next_edge_id(nx_id)
		{}

	};

	struct Vertex {

		glm::vec3 p;

		Vertex(glm::vec3 pos) : p(pos) {}

	};

	struct edge_with_endpoints {

		int edge_index;
		int endpoint_1_index;
		int endpoint_2_index;


	};

	std::vector<half_edge> edges;
	std::vector<Vertex> vertices;


	std::vector<Vertex> get_vertex_neighbours(int vertex_id) {
	
		std::vector<Vertex> out;
		std::vector<int> v_edges = get_edges_of_vertex(vertex_id);

		for (int i : v_edges) {
		
			out.push_back(vertices[ edges[i].vertex_id  ]);
		
		}

		return out;

	}

	std::vector<std::vector<glm::vec3>> get_all_neighbouring_triangles(int vertex_id) {
	
		std::vector<std::vector<glm::vec3>> out;

		std::vector<int> v_edges = get_edges_of_vertex(vertex_id);

		for (int i : v_edges) {

			int l = i;
			std::vector<glm::vec3> face;

			face.push_back(vertices[edges[l].vertex_id].p);
			l = edges[l].next_edge_id;

			while (l != i) {

				if (l < edges.size()) {

					face.push_back(vertices[edges[l].vertex_id].p);
					l = edges[l].next_edge_id;

				}
				else {
				
					std::cout << "EDGE INDEX ERROR INDEXING " << l << " FROM A VECTOR OF SIZE " << edges.size() << std::endl;
					break;
				}
			}

			out.push_back(face);

		}

		return out;
	}

	float get_triangle_area(std::vector<glm::vec3> tri) {
	
		return glm::length(glm::cross(tri[0],tri[1])) * 0.5f;
	
	}

	glm::vec3 get_triangle_normal(std::vector<glm::vec3> tri) {
	
		return glm::normalize(glm::cross(tri[0], tri[1]));
	
	}

	float get_angle_of_triangle(std::vector<glm::vec3> tri) {
	
		return glm::acos(glm::dot(tri[0], tri[2]) / ((glm::length(tri[0]) * glm::length(tri[2]))));
	
	}

	void add_vertex(Vertex v) {

			vertices.push_back(v);

		}

	void add_edge(int vtx_id, int op_id, int nx_id) {

			edges.push_back(half_edge(vtx_id, op_id, nx_id));

		}

	int get_end_vertex_id(int edge_id) {

			return get_next_edge(edge_id).vertex_id;
			//edges[edges[edge_id].opposite_edge_id].vertex_id;

		}

	half_edge get_next_edge(int edge_id) {

			return edges[edges[edge_id].next_edge_id];

		}

	std::vector<int> enumerate_all_faces() {

			std::vector<int> output;

			for (int i = 0; i < edges.size(); ++i) {

				std::vector<int> edges_in_face = get_edges_of_face(i);

				int minimum = 10000;
				for (int k = 0; k < edges_in_face.size(); ++k) {

					if (edges_in_face[k] < minimum) {

						minimum = edges_in_face[k];

					}
				}




				if (std::find(output.begin(), output.end(), minimum) == output.end()) {

					//if(minimum >= 0)
					output.push_back(minimum);

				}


			}

			return output;

		}

	std::vector<int> enumerate_all_full_edges() {

			std::vector<int> output;

			for (int i = 0; i < edges.size(); ++i) {

				if (std::find(output.begin(), output.end(), edges[i].opposite_edge_id) == output.end()) {

					/*if (i < output[index]) {
						output[index] = i;
					}*/
					output.push_back(i);
				}

			}

			return output;
		}

	glm::vec3 get_face_point(int edge_index) {

			float n = 0;
			glm::vec3 sum = glm::vec3(0, 0, 0);

			int loop_edge = edge_index;
			//if (edge_index >= edges.size()) { return glm::vec3(10, 5, 0); }
			do {

				n++;
				sum = sum + vertices[edges[loop_edge].vertex_id].p;
				loop_edge = edges[loop_edge].next_edge_id;

			} while (loop_edge != edge_index);

			return sum / n;

		}

	std::vector<int> get_edges_of_vertex(int vertex_index) {

			std::vector<int> output;

			for (int i = 0; i < edges.size(); ++i) {

				if (edges[i].vertex_id == vertex_index) {

					output.push_back(i);
				}
			}
			return output;
		}

	glm::vec3 get_half_point_of_edge(int edge_index) {


			return (vertices[edges[edge_index].vertex_id].p +
				vertices[
					edges[
						edges[edge_index].opposite_edge_id
					].vertex_id
				].p)
				/ 2.f;

		}

	std::vector<int> get_edges_of_face(int face_edge_id) {

			std::vector<int> edges_in_face;
			if (face_edge_id >= edges.size()) { return std::vector<int>(); }
			int loop_edge = face_edge_id;
			do {

				edges_in_face.push_back(loop_edge);
				loop_edge = edges[loop_edge].next_edge_id;

			} while (loop_edge != face_edge_id);

			return edges_in_face;

		}

	glm::vec3 get_vertices_average_position() {

			glm::vec3 sum = glm::vec3(0, 0, 0);

			for (Vertex v : vertices) {

				sum += v.p;

			}

			float n = (float)vertices.size();

			return sum / n;

		}

	void recalculate_opposites() {

			std::vector<edge_with_endpoints> full_edges_duplicate;

			for (int i = 0; i < edges.size(); ++i) {

				full_edges_duplicate.push_back(edge_with_endpoints{ i, edges[i].vertex_id, get_next_edge(i).vertex_id });

			}

			for (int i = 0; i < full_edges_duplicate.size(); ++i) {

				int endpoint_1_1 = full_edges_duplicate[i].endpoint_1_index;
				int endpoint_2_1 = full_edges_duplicate[i].endpoint_2_index;

				for (int j = 0; j < full_edges_duplicate.size() && j != i; ++j) {

					int endpoint_1_2 = full_edges_duplicate[j].endpoint_1_index;
					int endpoint_2_2 = full_edges_duplicate[j].endpoint_2_index;

					if ((endpoint_1_1 == endpoint_2_2 && endpoint_1_2 == endpoint_2_1) || (endpoint_1_1 == endpoint_1_2 && endpoint_2_2 == endpoint_2_1)) {

						//match
						if (true) {
							edges[full_edges_duplicate[i].edge_index].opposite_edge_id = full_edges_duplicate[j].edge_index;
							edges[full_edges_duplicate[j].edge_index].opposite_edge_id = full_edges_duplicate[i].edge_index;
						}
					}

				}

			}

		}

};

