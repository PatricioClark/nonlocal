# Calculate the nonlocal SGS model at different heights and different orders
# all for fixed R=5\Delta

# import modules
from   pylab import *
import scipy.signal
import spherical
import models

# Get info
dest   = 'data/'
dims, dx, y_domain_bot, dz, delta, bound_x, bound_z, RD = load(dest+'params.npy',
                                                           allow_pickle=True)
_, _, y_domain_top, _, _, _, _, _ = load(dest+'params_top.npy',
                                         allow_pickle=True)
R  = 5*delta
dr = dz*3

# Create domain and sample at LES resolution
x_domain = arange(spherical.Nx)*dx
z_domain = arange(spherical.Nz)*dz
xs       = x_domain[bound_x*RD:-bound_x*RD]
xs       = xs[::bound_x//2]
zs       = z_domain[bound_z*RD:-bound_z*RD]
zs       = zs[::bound_z//2]

for t0 in range(1,109):
    strain = load(dest+'filt_strain.{}.npy'.format(t0))

    # Send filtered strain field
    y_domain = y_domain_bot
    spherical.set_domain(dx,y_domain,dz)

    # Iterate over heights
    for ii, y0 in enumerate([120,70,40,23,256],start=1):

        # Changed to top planes?
        if ii>=5:
            strain = load(dest+'filt_strain_top.{}.npy'.format(t0))
            y0_data = 125
            y0 = y0-y0_data
            y_domain = y_domain_top
            spherical.set_domain(dx,y_domain,dz)

        # Get local strain
        local  = strain[:,:,bound_x*RD:-bound_x*RD,y0:y0+1,bound_z*RD:-bound_z*RD]
        local  = local[:,:,::bound_x//2,0:1,::bound_z//2]

        ys = y_domain[y0]
        xx, yy, zz = meshgrid(xs, ys, zs, indexing='ij')
        nonloc     = zeros((dims,dims,*shape(xx)))
        result     = zeros(shape(xx))

        # iterate over alpha
        for jj, alpha in zip([2,5,7],[0.2,0.5,0.7]):

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

            # Enforce symmetry
            nonloc[1,0] = nonloc[0,1]
            nonloc[2,0] = nonloc[0,2]
            nonloc[2,1] = nonloc[1,2]
            
            # Calculate and save the nonlocal stress model
            tau_nl = models.get_tau_nonl(nonloc,local,delta,alpha)[:,:,:,0,:]
            save(dest+'tau_nl.{}.{}.{}.npy'.format(t0, ii, jj), tau_nl)
    print(t0)
