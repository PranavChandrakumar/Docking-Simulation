import numpy as np
import sympy as sp
from scipy.optimize import newton
from scipy.linalg import norm
import Parameters
mu = Parameters.Planet_mu

def normalize_vector(vector):
    return vector / np.linalg.norm(vector)

def are_vectors_equal(vec1, vec2, tolerance=1e-7):
    return np.allclose(vec1, vec2, atol=tolerance)


def InclinationChangePosition(SpacecraftPos, TargetPos, SpacecraftVel, TargetVel, dt_array):
    # Convert inputs to numpy arrays
    SpacecraftPos, TargetPos = np.array(SpacecraftPos), np.array(TargetPos)
    
    # Normalize vectors in both arrays
    normalized_spacecraft = np.array([normalize_vector(pos) for pos in SpacecraftPos])
    normalized_target = np.array([normalize_vector(pos) for pos in TargetPos])
    
    # Compare normalized vectors
    equal_vectors = []
    for i, sc_vec in enumerate(normalized_spacecraft):
        for j, targ_vec in enumerate(normalized_target):
            if are_vectors_equal(sc_vec, targ_vec):
                equal_vectors.append((i, j))
    
    if equal_vectors:
        print(f"Found {len(equal_vectors)} pairs of equal normalized vectors:")
        for sc_idx, targ_idx in equal_vectors:
            print(f"Spacecraft index: {sc_idx}, Target index: {targ_idx}")
    else:
        print("No equal normalized vectors found.")

    # Placeholder for deltaPos calculation
    deltaPos = np.zeros_like(SpacecraftPos[0])  # Assuming deltaPos is a single vector

    return deltaPos