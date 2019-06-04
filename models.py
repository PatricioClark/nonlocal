# Calculate different SGS stress models in 3D

# import modules
import numpy as np
import spherical

# Everything is in 3D
dims = 3

def remove_trace(field, print_trace=False):
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
    if print_trace: print(np.sum(trace))
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

def correlations(A,B,normalized=False,fluctuations=True):
    """
    Calculate the correlations between tensors A and B. A and B can be zero,
    one, or two-ranked tensors (each rank is of order dims), and one, two, or
    three-dimensional fields.
    
    Parameters
    ----------
    A, B : ndarray
        The two fields to be correlated.
    normalized : bool, optional
        If True, normalize the output. Default is False.
    fluctuations : bool, optional
        If True, calculate the correlations between the fluctuations, i.e.
        remove the mean of A and B. Default is True.

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
    if   len(sa)==1:
        idxs     = [()]
        sum_axis = 0
    elif sa[0]==dims and sa[1]==dims:
        idxs     = [(i,j) for i in range(dims) for j in range(dims)]
        sum_axis = tuple([-1*(i+1) for i in range(len(sa)-2)])
    elif sa[0]==dims and sa[1]>dims:
        idxs = range(dims)
        sum_axis = tuple([-1*(i+1) for i in range(len(sa)-1)])
    elif sa[0]>dims and sa[1]>dims:
        idxs = [()]
        sum_axis = tuple([-1*(i+1) for i in range(len(sa))])

    prod = np.mean(A,axis=sum_axis)*np.mean(B,axis=sum_axis)
    if not fluctuations: prod = np.zeros((dims,dims))
    prod = np.array([A[ip]*B[ip]-prod[ip] for ip in idxs])
    nume = np.mean(prod)

    if normalized:
        prod = np.mean(A,axis=sum_axis)*np.mean(A,axis=sum_axis)
        if not fluctuations: prod = np.zeros((dims,dims))
        prod = np.array([A[ip]*A[ip]-prod[ip] for ip in idxs])
        s1   = np.sqrt(np.mean(prod))

        prod = np.mean(B,axis=sum_axis)*np.mean(B,axis=sum_axis)
        if not fluctuations: prod = np.zeros((dims,dims))
        prod = np.array([B[ip]*B[ip]-prod[ip] for ip in idxs])
        s2   = np.sqrt(np.mean(prod))

        return nume/(s1*s2)
    else:
        return nume

def two_point_corr(A,B,axes=2,one_directional=False,fluctuations=True):
    """
    Calculate the two-point correlations between two-rank tensors A and B along
    a given direction.

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

    L = sa[axes]//2
    A = np.swapaxes(A,axes,-1)
    B = np.swapaxes(B,axes,-1)

    corrs = [correlations(A[...,:L],B[...,:L],fluctuations=fluctuations)]
    for r in range(1,L):
        local = slice(None,L)
        prime = slice(r,r+L)
        aux1 = correlations(A[...,local],B[...,prime],fluctuations=fluctuations)
        aux2 = correlations(A[...,prime],B[...,local],fluctuations=fluctuations)
        if one_directional:
            corrs.append(0.5*(aux1+aux2))
        else:
            corrs = [aux2]+corrs+[aux1]
    return np.array(corrs)

