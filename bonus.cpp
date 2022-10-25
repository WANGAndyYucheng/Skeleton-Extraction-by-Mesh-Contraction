#include "deformer.h"
#include <iostream>

Deformer::Deformer() : mMesh(nullptr),
                       mCholeskySolver(nullptr) {
}

Deformer::~Deformer() {
	clear();
}

void Deformer::clear() {
	if (mCholeskySolver) {
		delete mCholeskySolver;
	}
	mCholeskySolver = nullptr;
	mRoiList.clear();
}

void Deformer::setMesh(Mesh* mesh) {
	mMesh = mesh;
	clear();
	// Record the handle vertices
	for (Vertex* vert : mMesh->vertices()) {
		if (vert->flag() > 0 || vert->isBoundary()) {
			mRoiList.push_back(vert);
		}
	}
	// Build system matrix for deformation
	buildSystemMat();
}


void Deformer::buildSystemMat() {
	/*====== Programming Assignment 2 ======*/

	/**********************************************/
	/*          Insert your code here.            */
	/**********************************************/
	/*
	/* Build the matrix of the linear system for 
	/* deformation and do factorization, in order
	/* to reuse and speed up in Deformer::deform().
	/* Handle vertices are maked by Vertex::flag() > 0
	/* Movements of the specified handle are already
	/* recorded in Vertex::position()
	/**********************************************/

	// get the number of vertices and constraints
    int n_vertices = mMesh->vertices().size();
    int n_constraints = mRoiList.size();

	if(n_inter == 1){
		ori_area.resize(n_vertices+200,0);
		curr_area.resize(n_vertices+200,0);
	}

	// create init mA and mB
    mA = Eigen::SparseMatrix<double>(n_vertices + n_constraints, n_vertices);
    mB = Eigen::MatrixX3d(n_vertices + n_constraints, 3);
	
    for (Vertex *vert : mMesh->vertices()) {
		// get the index
        int v_i = vert->index();

		// get adjacent vertices
        OneRingVertex one_ring(vert);
        Vertex *curr = nullptr;
        std::vector<Vertex *> adj_vertices;
        while (curr = one_ring.nextVertex()) {
            adj_vertices.push_back(curr);
        }
		int n = adj_vertices.size();

		// get the weights
        std::vector<double> weights(n);
        double weightSum = 0;
		double areasum = 0;
        for (int k = 0; k < n; ++k) {
			// same as assignment one, get the cot value and normalized weight
            const Eigen::Vector3f &prev = adj_vertices[k]->position();
            const Eigen::Vector3f &curr = adj_vertices[(k + 1) % n]->position();
            const Eigen::Vector3f &next = adj_vertices[(k + 2) % n]->position();
            double cot1 = triangleCot(vert->position(), prev, curr);
            double cot2 = triangleCot(vert->position(), next, curr);
            double weight = cot1 + cot2;
            weights[(k + 1) % n] = weight;
            weightSum += weight;
			areasum += triangleArea(vert->position(), next, curr);
        }
		if(n_inter == 1){
			ori_area[v_i] = areasum;
		}
		curr_area[v_i] = areasum;
		//calculate delta_i and generate mA and mB
        for (int k = 0; k < n; ++k) {
            int j = adj_vertices[k]->index();
            mA.insert(v_i, j) = weights[k];
        }
        mA.insert(v_i, v_i) = -1 * weightSum;
        mB.row(v_i) = Eigen::MatrixXd::Zero(1,3);
    }
	
	mA = mA * n_inter;
	if(n_inter == 1){
		for (int k = 0; k < n_constraints; ++k) {
			int i = mRoiList[k]->index();
			Wl0 += ori_area[i];
		}
		Wl0 = sqrt(Wl0) / 1000;
	}
	std:: cout << Wl0 << std::endl;
	// for the rest of the mA
    for (int k = 0; k < n_constraints; ++k) {
        int i = mRoiList[k]->index();
		if (curr_area[i] != 0)
        	mA.insert(n_vertices + k, i) = Wl0 * sqrt(ori_area[i]/curr_area[i]);
		else 
			mA.insert(n_vertices + k, i) = Wl0;
    }

	n_inter = n_inter * 2;
    Eigen::SparseMatrix< double > systemMat = mA.transpose() * mA;

	/*====== Programming Assignment 2 ======*/

	// Please refer to the following link for the usage of sparse linear system solvers in Eigen
	// https://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html

	// Do factorization
	if (systemMat.nonZeros() > 0) {
		mCholeskySolver = new Eigen::SimplicialLDLT< Eigen::SparseMatrix< double > >();
		mCholeskySolver->compute(systemMat);
		if (mCholeskySolver->info() != Eigen::Success) {
			// Decomposition failed
			std::cout << "Sparse decomposition failed\n";
		} else {
			std::cout << "Sparse decomposition succeeded\n";
		}
	}
}

void Deformer::deform() {
	if (mCholeskySolver == nullptr) {
		return;
	}

	/*====== Programming Assignment 2 ======*/

	/**********************************************/
	/*          Insert your code here.            */
	/**********************************************/
	/*
	/* This is the place where the editing techniques 
	/* take place.
	/* Solve for the new vertex positions after the 
	/* specified handles move using the factorized
	/* matrix from Deformer::buildSystemMat(), i.e.,
	/* mCholeskySolver defined in deformer.h
	/**********************************************/

	// Please refer to the following link for the usage of sparse linear system solvers in Eigen
	// https://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html

	// get the number of vertices and constraints
	int n_vertices = mMesh->vertices().size();
    int n_constraints = mRoiList.size();
	
	// for the rest of the mB
    for (int k = 0; k < n_constraints; ++k) {
		int i = mRoiList[k]->index();
		if (curr_area[i] != 0)
        	mB.row(n_vertices + k) = (mRoiList[k]->position().cast<double>()) * Wl0 * sqrt(ori_area[i]/curr_area[i]);
		else 
			mB.row(n_vertices + k) = (mRoiList[k]->position().cast<double>()) * Wl0;
    }
	// calculate new position	
    Eigen::MatrixX3f newPositions(n_vertices, 3);
    for (int d = 0; d < 3; ++d) {
		// rhs matrix = mA.T * mB
        Eigen::VectorXd rhs = mA.transpose() * mB.col(d);
        Eigen::VectorXd x = mCholeskySolver->solve(rhs);
		// rechange the ans to float
        newPositions.col(d) = x.cast<float>();
    }

	// set new position for all the vertices
    for (Vertex *vert : mMesh->vertices()) {
        int i = vert->index();
		if (curr_area[i] <= 0.7)
			vert->setColor(VCOLOR_RED);
        Eigen::Vector3f newPosition = newPositions.row(i);
        vert->setPosition(newPosition);
    }
	mMesh->setVertexColorDirty(true);
	/*====== Programming Assignment 2 ======*/
}
