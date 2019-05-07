# Calculate different SGS stress models in 3D

# import modules
import numpy as np
import spherical

# Everything is in 3D
dims = 3

def remove_trace(field):
    """
    Make field traceless.

    Parameters
    ----------
    field : ndarray
        Array to be made traceless. Shape must be (dims,dims,:)

    Returns
    -------
    field : ndarray
        Same field, but without the trace.
    """ 
    trace = np.sum([field[i,i] for i in range(dims)], axis=0)
    for i in range(dims):
        field[i,i] -= (1/dims)*trace
    return field

def get_tau_dns(filt_velos, filt_prods):
    """
    Get the exact SGS stresses. Result is already traceless.

    Parameters
    ----------
    filt_velos : ndarray
        Filtered velocity fields.
    filt_prods : ndarray
        Filtered product fields.

    Returns
    -------
    tau : exact SGS stress, already traceless.
    """
    tau_dns = np.array([[filt_prods[i,j] - filt_velos[i]*filt_velos[j]
                         for i in range(dims)]
                         for j in range(dims)])
    return remove_trace(tau_dns)

def tensor_dot(A,B):
    """
    Calculate the dot product between tensors A and B.
    
    Parameters
    ----------
    A, B : ndarray
        The two fields to be correlated.

    Returns
    -------
    result : float
    """
    aux  =    A[0,0]*B[0,0] + A[1,1]*B[1,1] + A[2,2]*B[2,2]
    aux += 2*(A[0,1]*B[0,1] + A[0,2]*B[0,2] + A[1,2]*B[1,2])
    return aux


def get_tau_smag(strain,delta,c_s=0.16):
    """
    Get the smagorinsky SGS stresses. Result is already traceless.

    Parameters
    ----------
    strain : ndarray
        Filtered strain rates.

    Returns
    -------
    tau : Smagorinsky stress tensor, result is already traceless.
    """
    # Calculate eddy viscosity
    char_strain = np.sqrt(2*tensor_dot(strain,strain))
    nu_smag = (c_s*delta)**2*char_strain

    # Calcualte Smagorinsky stress tensor
    tau_smag = -2*nu_smag*strain
    return remove_trace(tau_smag)

def correlations(A,B):
    """
    Calculate the correlations between tensors A and B, minus their respective
    means. The result is normalized.
    
    Parameters
    ----------
    A, B : ndarray
        The two fields to be correlated.

    Returns
    -------
    result : float
        Correlation coefficient of fields A and B.
    """
    prod = np.mean(A,axis=(-1,-2))*np.mean(B,axis=(-1,-2))
    prod = np.array([[A[i,j]*B[i,j]-prod[i,j] for i in range(dims)]
                                              for j in range(dims)])
    nume = np.mean(prod)

    prod = np.mean(A,axis=(-1,-2))*np.mean(A,axis=(-1,-2))
    prod = np.array([[A[i,j]*A[i,j]-prod[i,j] for i in range(dims)]
                                              for j in range(dims)])
    s1   = np.sqrt(np.mean(prod))

    prod = np.mean(B,axis=(-1,-2))*np.mean(B,axis=(-1,-2))
    prod = np.array([[B[i,j]*B[i,j]-prod[i,j] for i in range(dims)]
                                              for j in range(dims)])
    s2   = np.sqrt(np.mean(prod))

    return nume/(s1*s2)
