import numpy as np

def unit_vector(data, axis=None):
    """Return ndarray normalized by length, i.e. Euclidean norm, along axis.
    """
    data = np.array(data, dtype=np.float64, copy=True)
    if data.ndim == 1:
        data /= math.sqrt(np.dot(data, data))
        return data
        
    length = np.atleast_1d(np.sum(data * data, axis))
    np.sqrt(length, length)
    if axis is not None:
        length = np.expand_dims(length, axis)
    data /= length
    return data
    
if __name__ == "__main__":
    basis_from = np.array([[10, 0, 0], [0, 100, 0], [0, 0, 1000]])
    print(f"normalized vector: {unit_vector(basis_from, axis = 1)}")

def load_rlt(force, moment, origin_from, basis_from, origin_to, basis_to):
  """
  Transforms a force vector and calculates its moment in another coordinate system.

  Args:
    force: A 3D numpy array representing the force vector in the source system.
    moment: A 3D numpy array representing the moment vector in the source system.
    origin_from: A 3D numpy array representing the origin of the source system.
    basis_from: A 3x3 numpy array representing the basis vectors of the source system (from global)
    origin_to: A 3D numpy array representing the origin of the target system.
    basis_to: A 3x3 numpy array representing the basis vectors of the target system(from global)

  Returns:
    A tuple containing:
      transformed_force: A 3D numpy array representing the transformed force vector.
      transformed_moment: A 3D numpy array representing the moment in the target system.
  """

  # Check if basis vectors are orthonormal (normalized and orthogonal)
  for i in range(3):
    for j in range(i + 1, 3):
      dot_product = np.dot(basis_from[:, i], basis_from[:, j])
      if abs(dot_product) > 1e-6:
        raise ValueError("Basis vectors in basis_from must be orthonormal.")
      dot_product = np.dot(basis_to[:, i], basis_to[:, j])
      if abs(dot_product) > 1e-6:
        raise ValueError("Basis vectors in basis_to must be orthonormal.")
          
  # normalized basis (axis)        
  basis_from = unit_vector(basis_from, axis=1)
  basis_to = unit_vector(basis_to, axis=1)
    
  # Transform force and moment vector
  transformed_force  = np.dot(basis_to.T, np.dot(basis_from, force))
  transformed_moment = np.dot(basis_to.T, np.dot(basis_from, moment))
    
  # Calculate position vector relative to origin of target system
  position_global = np.dot(basis_to.T, (origin_from - origin_to))
  # Calculate transformed moment components
  transformed_force_moment = np.cross(position_global, transformed_force)
  transformed_total_force = transformed_force_moment + transformed_moment
  return transformed_force, transformed_total_force


if __name__ == "__main__":
    force = np.array([10, 5, 2])
    moment = np.array([10, 5, 2])
    
    origin_from = np.array([0, 0, 0])
    basis_from = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    origin_to = np.array([1, 2, 3])
    basis_to = np.array([[0.707, -0.707, 0], [0.707, 0.707, 0], [0, 0, 1]])
    
    transformed_force, transformed_moment = load_rlt(force,moment, origin_from, basis_from, origin_to, basis_to)
    print(f"force RLT = {transformed_force}")
    print(f"Moment RLT= {transformed_moment}")
