#include "Mesh.h"

std::vector<Vertex *> V;
std::vector<Face *> F;

Vertex::Vertex(float x, float y, float z, int id) {
	this->x = x;
	this->y = y;
	this->z = z;
	this->id = id;
	V.push_back(this);
}

Vertex::~Vertex() {
	while (neighbour.size()) {
		Remove(neighbour[0]->neighbour, this);
		Remove(neighbour, neighbour[0]);
	}
	Remove(V, this);
}


void Vertex::removeIfNonNeighbour(Vertex *n) {
	int hasNeighbour = std::count(begin(neighbour), end(neighbour), n);
	if (!hasNeighbour) return;
	for (int i = 0; i<face.size(); i++) {
		if (face[i]->hasVertex(n)) return;
	}
	Remove(neighbour, n);
}


Face::Face(Vertex *v0, Vertex *v1, Vertex *v2) {
	vertex[0] = v0;
	vertex[1] = v1;
	vertex[2] = v2;
	computeNormal();
	F.push_back(this);
	for (int i = 0; i < 3; i++) {
		vertex[i]->face.push_back(this);
		for (int j = 0; j < 3; j++) {
			if (i != j) {
				int hasVertexJ = std::count(begin(vertex[i]->neighbour), end(vertex[i]->neighbour), vertex[j]);
				if (!hasVertexJ) vertex[i]->neighbour.push_back(vertex[j]);
			}
		}
	}
}


Face::~Face() {
	Remove(F, this);
	for (int i = 0; i < 3; i++) {
		if (vertex[i]) Remove(vertex[i]->face, this);
	}

	// remove a vertex from another vertex if the former vertex is not a neighbour of latter
	for (int i = 0; i < 3; i++) {
		int j = (i + 1) % 3;
		if (!vertex[i] || !vertex[j])
			continue;
		vertex[i]->removeIfNonNeighbour(vertex[j]);
		vertex[j]->removeIfNonNeighbour(vertex[i]);
	}
}

//Does the 3 points of the face contain the vertex
bool Face::hasVertex(Vertex *v) {
	return (v == vertex[0] || v == vertex[1] || v == vertex[2]);
}


Mesh::Mesh(const char* filename) {
	loadMF(filename);
}

void Mesh::loadMF(const char* filename) {
	if (Vcnt()>0) V.clear();
	if (Fcnt()>0) F.clear();
	std::ifstream infile;
	infile.open(filename, std::ios::in);
	std::string strbuff;
	int id = 1;
	while (std::getline(infile, strbuff)) {
		std::stringstream ss;
		ss << strbuff;
		char type;
		ss >> type;
		if (type == 'v') {
			float vx, vy, vz;
			ss >> vx >> vy >> vz;
			Vertex *v = new Vertex(vx, vy, vz, id);
			id++;
		}
		else if (type == 'f') {
			int a, b, c;
			ss >> a >> b >> c;
			Face *f = new Face(V[a - 1], V[b - 1], V[c - 1]);
		}
	}
	infile.close();
}

void Mesh::writeMF(const char* filename) {
	std::ofstream outfile;
	outfile.open(filename, std::ios::out);
	std::string strbuff;
	for (uint i = 0; i<Vcnt(); i++) {
		outfile << "v " << V[i]->x << " " << V[i]->y << " " << V[i]->z << std::endl;
	}
	for (uint i = 0; i<Fcnt(); i++) {
		outfile << "f " << F[i]->vertex[0]->id << " " << F[i]->vertex[1]->id << " " << F[i]->vertex[2]->id << std::endl;
	}
	outfile.close();
}

void Mesh::simplifyMesh(const char* input, const char* output, int faceCnt) {
	// you may assume inputs are always valid
	loadMF(input);
	std::cout << "Original face count: " << Fcnt() << std::endl;
	for (int i = 0; i < Vcnt(); i++) {
		computeEdgeCostAtVertex(V[i]);
	}

	// while face is more than desired amount
	while (Fcnt() > faceCnt) {
		Vertex* minimum = V[0];
		for (int i = 0; i < Vcnt(); i++) {
			if (V[i]->cost < minimum->cost) {
				minimum = V[i];
			}
		}
		collapse(minimum, minimum->collapse);
	}

	//update ids
	for (int i = 0; i < Vcnt(); i++) {
		V[i]->id = i + 1;
	}

	std::cout << "Finished collapsing" << std::endl;
	writeMF(output);
}

int Mesh::Vcnt() {
	return V.size();
}

int Mesh::Fcnt() {
	return F.size();
}

float Mesh::computeEdgeCollapseCost(Vertex *u, Vertex *v) {

	float edgeLength = sqrtf((v->x - u->x)*(v->x - u->x) + (v->y - u->y)*(v->y - u->y) + (v->z - u->z)*(v->z - u->z));
	float curvature = 0;

	// faces that have edge uv
	std::vector<Face*> sides;
	for (int i = 0; i < u->face.size(); i++) {
		if (u->face[i]->hasVertex(v)) {
			sides.push_back(u->face[i]);
		}
	}

	//determine curvature term
	for (int i = 0; i<u->face.size(); i++) {
		float minCurv = 1;
		for (int j = 0; j<sides.size(); j++) {
			float dotProd = ((u->face[i]->nx * sides[j]->nx) + (u->face[i]->ny * sides[j]->ny) + (u->face[i]->nz * sides[j]->nz));
			minCurv = std::min(minCurv, (1 - dotProd) / 2.0f);
		}
		curvature = std::max(curvature, minCurv);
	}
	return edgeLength * curvature;
}


void Mesh::computeEdgeCostAtVertex(Vertex *v) {
	if (v->neighbour.size() == 0) {
		v->collapse = NULL;
		v->cost = -0.01f;
		return;
	}
	v->cost = 1000000;
	v->collapse = NULL;

	// search all neighbouring edges for "least cost" edge
	for (int i = 0; i < v->neighbour.size(); i++) {
		float c = computeEdgeCollapseCost(v, v->neighbour[i]);
		if (c < v->cost) {
			v->collapse = v->neighbour[i];
			v->cost = c;
		}
	}
}

void Mesh::collapse(Vertex *u, Vertex *v) {
	// Collapse the edge uv by moving vertex u onto v
	if (!v) { //u is a vector all by itself so can just delete
		delete u;
		return;
	}

	std::vector<Vertex*> tmp; //make temp a vector of all the neighbours of u
	for (int i = 0; i < u->neighbour.size(); i++) {
		tmp.push_back(u->neighbour[i]);
	}

	//delete faces on uv
	for (int i = u->face.size() - 1; i >= 0; i--) {
		if (u->face[i]->hasVertex(v)) {
			delete(u->face[i]);
		}
	}

	// update remaining faces to have v instead of u
	for (int i = u->face.size() - 1; i >= 0; i--) {
		u->face[i]->replaceVertex(u, v);
	}

	delete u;

	// recompute the edge collapse costs in neighborhood
	for (int i = 0; i<tmp.size(); i++) {
		computeEdgeCostAtVertex(tmp[i]);
	}
}


void Face::computeNormal() {
	float ax = vertex[0]->x - vertex[1]->x;
	float ay = vertex[0]->y - vertex[1]->y;
	float az = vertex[0]->z - vertex[1]->z;

	float bx = vertex[2]->x - vertex[1]->x;
	float by = vertex[2]->y - vertex[1]->y;
	float bz = vertex[2]->z - vertex[1]->z;

	//cross pdt
	nx = (ay * bz) - (az * by);
	ny = (az * bx) - (ax * bz);
	nz = (ax * by) - (ay * bx);

	float length = sqrt((nx*nx) + (ny*ny) + (nz*nz));
	if (length > 0.0f) {
		nx = nx / length;
		ny = ny / length;
		nz = nz / length;
	}

}

void Face::replaceVertex(Vertex *oldVertex, Vertex *newVertex) {
	if (oldVertex == vertex[0]) {
		vertex[0] = newVertex;
	}
	else if (oldVertex == vertex[1]) {
		vertex[1] = newVertex;
	}
	else {
		vertex[2] = newVertex;
	}

	Remove(oldVertex->face, this);
	newVertex->face.push_back(this);

	for (int i = 0; i < 3; i++) {
		oldVertex->removeIfNonNeighbour(vertex[i]);
		vertex[i]->removeIfNonNeighbour(oldVertex);
	}

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			if (i != j) {
				int hasVertexJ = std::count(begin(vertex[i]->neighbour), end(vertex[i]->neighbour), vertex[j]);
				if (!hasVertexJ) vertex[i]->neighbour.push_back(vertex[j]);
			}
		}
	}

	computeNormal();

}