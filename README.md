# Skeleton Extraction by Mesh Contraction

![1.png](https://s2.loli.net/2022/12/15/ScfGwdQUuhM9zIm.png)
![2.png](https://s2.loli.net/2022/12/15/SckHQXUhjBGiZen.png)

## Methodology

![m.png](https://s2.loli.net/2022/12/15/fC4RkeZh5QDJG87.png)

## Core code for Deformer::deform()
```c++

// get the number of vertices and constraints
int n_vertices = mMesh->vertices().size();
int n_constraints = mRoiList.size();
	
// for the rest of the mB
for (int k = 0; k < n_constraints; ++k) {
  int i = mRoiList[k]->index();
  if (curr_area[i] != 0) mB.row(n_vertices + k) = (mRoiList[k]->position().cast<double>()) * Wl0 * sqrt(ori_area[i]/curr_area[i]);
  else mB.row(n_vertices + k) = (mRoiList[k]->position().cast<double>()) * Wl0;
}
  
// calculate new position	
Eigen::MatrixX3f newPositions(n_vertices, 3);
for (int d = 0; d < 3; ++d) {
  // rhs matrix = mA.T * mB
  Eigen::VectorXd rhs = mA.transpose() * mB.col(d);
  Eigen::VectorXd x = mCholeskySolver->solve(rhs);
  // change the ans to float
  newPositions.col(d) = x.cast<float>();
}

// set new position for all the vertices
for (Vertex *vert : mMesh->vertices()) {
  int i = vert->index();
  if (curr_area[i] <= 0.7) vert->setColor(VCOLOR_RED);
  Eigen::Vector3f newPosition = newPositions.row(i);
  vert->setPosition(newPosition);
}
	
mMesh->setVertexColorDirty(true);
  ```
