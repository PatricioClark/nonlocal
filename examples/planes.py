# Calculate the nonlocal SGS model at different heights and different orders
# all for fixed R=5\Delta

# import modules
from   pylab import *
import scipy.signal
import spherical
import models

# Get info
dest   = 'data/'
dims, dx, y_domain, dz, delta, bound_x, bound_z, RD = load(dest+'params.npy',
                                                           allow_pickle=True)
# y0_data = 125
R  = 5*delta
dr = dz*3

# Send filtered strain field
spherical.set_domain(dx,y_domain,dz)

# Create domain
x_domain = arange(spherical.Nx)*dx
z_domain = arange(spherical.Nz)*dz
xs       = x_domain[bound_x*RD:-bound_x*RD]
zs       = z_domain[bound_z*RD:-bound_z*RD]

for t0 in range(1,73):
    strain = load(dest+'filt_strain.{}.npy'.format(t0))

    # Iterate over heights
    for ii, y0 in enumerate([120,70,40,23],start=1):
    # for ii, y0 in enumerate([256-y0_data],start=5):
        ys = y_domain[y0]
        print(ys)
        xx, yy, zz = meshgrid(xs, ys, zs, indexing='ij')
        nonloc     = zeros((dims,dims,*shape(xx)))
        result     = zeros(shape(xx))

        # iterate over alpha
        for jj, alpha in enumerate(arange(0.1,1,0.1),start=1):

            # Loop over dimensions
            for i in range(dims):
                for j in range(dims):
                    if j<i: continue
                    # Send field
                    spherical.set_field(strain[i,j])

                    # Get nonlocal strain
                    nonloc[i,j,:] = spherical.IntegrateInVolume(result,
                                                                dr,
                                                                R,
                                                                alpha,
                                                                xx,
                                                                yy,
                                                                zz)

            # Calculate the nonlocal model
            nonloc[1,0] = nonloc[0,1]
            nonloc[2,0] = nonloc[0,2]
            nonloc[2,1] = nonloc[1,2]
            tau_nl = models.get_tau_smag(nonloc,delta*R**(alpha-1))
            save(dest+'tau_nl.{}.{}.{}.npy'.format(t0, ii, jj), tau_nl)
    print(t0)
