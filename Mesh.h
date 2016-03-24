#ifndef _MESH_H_
#define _MESH_H_

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>

template<typename T> void Remove(std::vector<T> & c, T t) { auto it = std::find(begin(c), end(c), t); c.erase(it); }

typedef unsigned int uint;

class Vertex;
class Face;

class Vertex {
public:
	float x; //3d coordinates
	float y;
	float z;

	int id; //vertex index in list
	std::vector<Vertex *> neighbour; //adjacent vertex
	std::vector<Face *> face; //adjacent faces
	float cost; //cost to collapse edge
	Vertex * collapse; //vertex for collapse

	Vertex(float x, float y, float z, int id); //vertex constructor
	~Vertex(); //vertex destructor

	void removeIfNonNeighbour(Vertex *n);
};

class Face {
public:
	Vertex* vertex[3]; //3 points of face
	float nx; // normal coordinates
	float ny;
	float nz;

	Face(Vertex *v0, Vertex *v1, Vertex *v2); //face constructor
	~Face(); //face destructor

	void computeNormal();
	void replaceVertex(Vertex *vold, Vertex *vnew); //replace with new vertex
	bool hasVertex(Vertex *v);
};

typedef struct {
	float x;
	float y;
	float z;
} VertexH;

typedef struct {
	uint a, b, c;
} FaceH;

//data structure of Halfedge
typedef struct {
	std::vector<VertexH>  oppositeHalfEdge; //twin half edge
	std::vector<VertexH>  nextHalfEdge; //next halfedge
	std::vector<VertexH>  prevHalfEdge; //previous halfedge
	VertexH *vertex_h; //vertex pointing to 
	FaceH *face_h; //adjacent face
} HEEdge;

typedef struct {
	VertexH *vertex_h;
	HEEdge *half_edge;
} HEVertex;

typedef struct {
	HEEdge *half_edge;
} HEFace;

class Mesh {
public:
	Mesh() {};
	Mesh(const char*);
	//load a Mesh from .mesh file
	void loadMF(const char*);
	//write a Mesh to .mesh file (no header)
	void writeMF(const char*);
	//simplify a mesh
	void simplifyMesh(const char* input, const char* output, int faceCnt);

	//compute edge collapse cost
	float computeEdgeCollapseCost(Vertex *u, Vertex *v);
	//compute edge cost at vertex
	void computeEdgeCostAtVertex(Vertex *v);
	//collapse edge
	void collapse(Vertex *u, Vertex *v);

	//return vertex count
	int Vcnt();
	//return face count
	int Fcnt();
};
#endif