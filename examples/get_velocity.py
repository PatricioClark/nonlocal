# Get velocity field from the JHUTDB

# The velocity fields are filtered and then the exact SGS stresses, the
# filtered viscous stresses, and the Smagorinsky SGS stresses are calculated.

# The filtered velocity fields saved are wider than the stresses because the
# boundaries are needed in the nonlocal integration. The stresses have the same
# horizontal dimensions as those that will be calcualted nonlocally.

# import modules
import numpy as np
import pyJHTDB
import scipy.signal
import models
import sys

# Initialize and add token
auth_token  = "edu.jhu.pato-ca56ca00" 
lJHTDB = pyJHTDB.libJHTDB()
lJHTDB.initialize()
lJHTDB.add_token(auth_token)

# Dimensions
dims        = 3
Nx          = 20
Ny          = 20
Nz          = 20
visc_length = 1.0006e-3
dt          = 0.0065
dx          = 8*np.pi/2048
dz          = 3*np.pi/1536 
dx_plus     = dx/visc_length
dz_plus     = dz/visc_length
y_points    = np.loadtxt('y_points.txt')
y_plus      = y_points/visc_length + 1/visc_length
y0          = 0
nu          = 5e-5
RD          = 5
fixed_frame = True

# Filter
delta_plus = 40
delta      = delta_plus*visc_length
delta_x    = round(delta_plus/dx_plus)
delta_z    = round(delta_plus/dz_plus)
bound_x    = 4 # approx delta/dx
bound_z    = 8 # approx delta/dz
cut_dims   = [Nx+2*bound_x, Ny, Nz+2*bound_z]
x_domain   = np.arange(Nx+2*bound_x)*dx
y_domain   = y_points[y0:y0+Ny]
z_domain   = np.arange(Nz+2*bound_z)*dz
filt_box   = np.full((delta_x,1,delta_z),fill_value=1.0/(delta_x*1*delta_z))
def les_filter(field):
    return scipy.signal.convolve(field, filt_box, mode='same')

t0 = int(sys.argv[1])
if t0<=40:
    tms = 5
    ti  = 0
    # Time
    if   t0%tms==1:
        tidx = 0
    elif t0%tms==2:
        tidx = 2000
    elif t0%tms==3:
        tidx = 3900
    elif t0%tms==4:
        tidx = 1000
    elif t0%tms==0:
        tidx = 3000
else:
    tms = 4
    ti  = 40
    # Time
    if   t0%tms==1:
        tidx = 500
    elif t0%tms==2:
        tidx = 1500
    elif t0%tms==3:
        tidx = 2500
    elif t0%tms==0:
        tidx = 3500

# Location
if   t0<=tms+ti:
    x0 = 0
    z0 = 0
elif t0<=tms*2+ti:
    x0 = 1000
    z0 = 0
elif t0<=tms*3+ti:
    x0 = 0
    z0 = 600
elif t0<=tms*4+ti:
    x0 = 1000
    z0 = 600
elif t0<=tms*5+ti:
    x0 = 500
    z0 = 0
elif t0<=tms*6+ti:
    x0 = 500
    z0 = 600
elif t0<=tms*7+ti:
    x0 = 1500
    z0 = 0
elif t0<=tms*8+ti:
    x0 = 1500
    z0 = 600
print(t0,tidx,x0,z0)

# Get velocity field (unfiltered)
if fixed_frame:
    points = np.zeros((cut_dims[2],cut_dims[1],cut_dims[0],3), np.float32)
    points[...,0] = x_domain[np.newaxis, np.newaxis, :] + x0*dx
    points[...,1] = y_domain[np.newaxis, :, np.newaxis]
    points[...,2] = z_domain[:, np.newaxis, np.newaxis] + z0*dz
    velocity = lJHTDB.getData(
        tidx*dt,
        points,
        data_set = 'channel',
        getFunction='getVelocity')
else:
    velocity = lJHTDB.getCutout(
        data_set = 'channel',
        field='u',
        time=tidx,
        start = np.array([x0, y0, z0], dtype = np.int),
        size  = np.array(cut_dims,  dtype = np.int),
        step  = np.array([1, 1, 1], dtype = np.int),
        filter_width = 1)

# Make the shape of velocity equal to (dims,Nx,Ny,Nz)
velocity = np.transpose(velocity)

# Get products
products = np.array([[velocity[i]*velocity[j] for i in range(dims)]
                                              for j in range(dims)])

# Get filtered velocities
filt_velos = np.array([les_filter(velocity[i])    for i in range(dims)])
filt_prods = np.array([[les_filter(products[i,j]) for i in range(dims)]
                                                  for j in range(dims)])

# Remove boundaries affected by the filtering
filt_velos = filt_velos[:,bound_x:-bound_x,:,bound_z:-bound_z]
filt_prods = filt_prods[:,:,bound_x:-bound_x,:,bound_z:-bound_z]

# Get exact SGS stress tensor
tau_dns = models.get_tau_dns(filt_velos, filt_prods)
tau_dns = tau_dns[:,:,bound_x*RD:-bound_x*RD,:,bound_z*RD:-bound_z*RD]

# Calculate filtered gradient
filt_grad = np.array([[np.gradient(filt_velos[i],difs,axis=j,edge_order=2)
                        for i in range(dims)]
                        for (difs,j) in zip((dx,y_domain,dz),range(dims))])

# Get viscous stress tensor
tau_visc = 2*nu*np.array([[filt_grad[i,0]*filt_grad[j,0] +
                           filt_grad[i,1]*filt_grad[j,1] + 
                           filt_grad[i,2]*filt_grad[j,2] for i in range(dims)]
                                                         for j in range(dims)])
tau_visc = tau_visc[:,:,bound_x*RD:-bound_x*RD,:,bound_z*RD:-bound_z*RD]

# Get strain rate
strain = 0.5*(filt_grad + np.swapaxes(filt_grad,0,1))

# Get Smagorinsky SGS stress
tau_smag = models.get_tau_smag(strain, delta)
tau_smag = tau_smag[:,:,bound_x*RD:-bound_x*RD,:,bound_z*RD:-bound_z*RD]

# Save fields
dest = 'data/'
np.save(dest+'params', [dims,
                        dx,
                        y_domain,
                        dz,
                        delta,
                        bound_x,
                        bound_z,
                        RD], allow_pickle=True)
# np.save(dest+'filt_strain_top.{}'.format(t0), strain)
# for ii, yi in enumerate([256,215,171],start=5):
np.save(dest+'filt_strain.{}'.format(t0), strain)
for ii, yi in enumerate([120,70,40,23],start=1):
    np.save(dest+'filt_strain.{}.{}'.format(t0,ii), strain[:,:,:,yi-y0,:])
    np.save(dest+'tau_dns.{}.{}'.format(t0,ii),     tau_dns[:,:,:,yi-y0,:])
    np.save(dest+'tau_visc.{}.{}'.format(t0,ii),    tau_visc[:,:,:,yi-y0,:])
    np.save(dest+'tau_smag.{}.{}'.format(t0,ii),    tau_smag[:,:,:,yi-y0,:])
