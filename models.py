# Calculate different SGS stress models in 3D

# import modules
import numpy as np

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

def get_tau_nonl(strain_nonl,strain_local,delta,alpha):
    """
    Get the nonlocal SGS stresses. Result is already traceless.

    Parameters
    ----------
    strain_nonl : ndarray
        Filtered non-local strain rates of order alpha.
    strain_local : ndarray
        Filtered local strain rates.
    delta : float
        Filter size
    alpha : float
        Non-local order

    Returns
    -------
    tau : Non-local stress tensor, result is already traceless.
    """
    # Calculate eddy viscosity
    C_K = 1.5
    char_strain = np.sqrt(2*tensor_dot(abs(strain_nonl),abs(strain_local)))
    nu_nonl  = ((alpha+1./3)/(2*C_K))**(3./2)
    nu_nonl *= (delta/np.pi)**((3.*alpha+1)/2)
    nu_nonl *= char_strain

    # Calcualte Non-local stress tensor
    tau_nonl = -2*nu_nonl*strain_nonl
    return remove_trace(tau_nonl)

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

def two_point_corr(A, B,
                   axes=2,
                   one_directional=False,
                   fluctuations=True,
                   homogeneous=True):
    """
    Calculate the two-point correlations between two-rank tensors A and B along
    a given direction.

    Parameters
    ----------
    A, B : ndarray
        The two fields to be correlated.
    axes : int, optional
        Direction into which the correlations are to be calculated. Default is 2.
    one_directional : bool, optional
        If True, it keeps both -r and +r directions. Otherwise it averages.
        If flow is not homogenous the branches actually mean A(x)*B(x+r) and
        A(x+r)*B(x), respectively.  Default is False.
    fluctuations : bool, optional
        If True, calculate the correlations between the fluctuations, i.e.
        remove the mean of A and B. Default is True.
    homogeneous : bool, optional 
        If True assume the fields are statistically homogenous and average on
        the direction of the correlations. If False it calculates the
        correlations between the starting x0 and x0+r only. Default is True.

    Returns
    -------
    result : float
        Correlation coefficient of fields A and B.
    """
    # Check dims
    sa = np.shape(A)
    if sa!=np.shape(B):
        raise TypeError('Arrays must have the same shape!')

    # Determine length
    Lc = sa[axes]//2
    if homogeneous:
        La = Lc
    else:
        La = 1

    # Change axes
    A = np.swapaxes(A,axes,-1)
    B = np.swapaxes(B,axes,-1)

    corrs = [correlations(A[...,:La],B[...,:La],fluctuations=fluctuations)]
    for r in range(1,Lc):
        local = slice(None,La)
        prime = slice(r,r+La)
        aux1 = correlations(A[...,local],B[...,prime],fluctuations=fluctuations)
        aux2 = correlations(A[...,prime],B[...,local],fluctuations=fluctuations)
        if one_directional:
            corrs.append(0.5*(aux1+aux2))
        else:
            corrs = [aux2]+corrs+[aux1]
    return np.array(corrs)

def symmetrize(corr):
    '''Makes correlation function symmetric/even'''
    return 0.5*(corr+corr[::-1])
