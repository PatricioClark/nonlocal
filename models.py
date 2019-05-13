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
        The two fields to be dotted.

    Returns
    -------
    result : float
    """
    aux  = np.array([[A[i,j]*B[i,j] for i in range(dims)]
                                    for j in range(dims)])
    return np.sum(aux, axis=(0,1))


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

def correlations(A,B,normalized=False):
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

    # Check dims
    sa = np.shape(A)
    if sa!=np.shape(B):
        raise TypeError('Arrays must have the same shape!')
    
    # Single component or multi? 1, 2 or 3D?
    if   sa[0]==dims and sa[1]==dims:
        idxs     = [(i,j) for i in range(dims) for j in range(dims)]
        sum_axis = tuple([-1*(i+1) for i in range(len(sa)-2)])
    elif sa[0]==dims and sa[1]>dims:
        idxs = range(dims)
        sum_axis = tuple([-1*(i+1) for i in range(len(sa)-1)])
    elif sa[0]>dims and sa[1]>dims:
        idxs = [()]
        sum_axis = tuple([-1*(i+1) for i in range(len(sa))])

    prod = np.mean(A,axis=sum_axis)*np.mean(B,axis=sum_axis)
    prod = np.array([A[ip]*B[ip]-prod[ip] for ip in idxs])
    nume = np.mean(prod)

    if normalized:
        prod = np.mean(A,axis=sum_axis)*np.mean(A,axis=sum_axis)
        prod = np.array([A[ip]*A[ip]-prod[ip] for ip in idxs])
        s1   = np.sqrt(np.mean(prod))

        prod = np.mean(B,axis=sum_axis)*np.mean(B,axis=sum_axis)
        prod = np.array([B[ip]*B[ip]-prod[ip] for ip in idxs])
        s2   = np.sqrt(np.mean(prod))

        return nume/(s1*s2)
    else:
        return nume
