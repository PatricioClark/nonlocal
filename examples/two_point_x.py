# Calculate the correlations between the different models and tau_dns

# import modules
from   pylab     import *
from   dom       import implot
from   itertools import cycle
from   scipy.optimize import curve_fit
from   dom import idx_nearest
import spherical
import models

# LaTex
rc('font',**{'size':16})
rc('text',usetex=True)
rc('text.latex', preamble=r'\usepackage{bm}')

# Data dir
ddir = 'data/'
dims, dx, y_domain, dz, delta, bound_x, bound_z, RD = load(ddir+'params.npy',
                                                           allow_pickle=True)
R = delta*RD
visc_length = 1.0006e-3
y_plus = y_domain/visc_length + 1/visc_length
mm = cycle(['^','o'])
heights = {120:1,70:2,40:3,23:4,256:5}
alfas = list(arange(0.1,1,0.1))
num_res = 108

# Direction
dr  = 2*dx
cps = 99
dps = 2*cps+1
dd  = 2
xxss = arange(-49,50)*dr/visc_length
dnxs = arange(-99,100)*dx/visc_length
correction = 0.5**1.5

# Different heights
for y0 in [23,40,70,120,256]:
    ii = heights[y0]
    try:
        corrs_dns  = load(ddir+'corrs_dns_x_{}.npy'.format(ii))
        corrs_smag = load(ddir+'corrs_smag_x_{}.npy'.format(ii))
        corrs_nonl = load(ddir+'corrs_nonl_x_{}.npy'.format(ii))
    except:
        corrs_dns  = zeros((num_res,dps))
        corrs_smag = zeros((num_res,cps))
        corrs_nonl = zeros((3,num_res,cps))
        for t0 in range(1,num_res+1):

            # Get fields
            strain   = load(ddir+'filt_strain.{}.{}.npy'.format(t0,ii))
            strain   = strain[:,:,bound_x*RD:-bound_x*RD,bound_z*RD:-bound_z*RD]
            tau_dns  = load(ddir+'tau_dns.{}.{}.npy'.format(t0,ii))
            tau_smag = load(ddir+'tau_smag.{}.{}.npy'.format(t0,ii))
            strmod   = strain[:,:,::bound_x//2,::bound_z//2]
            tau_smag = tau_smag[:,:,::bound_x//2,::bound_z//2]
            corrs_dns[t0-1]  = models.two_point_corr(tau_dns,  strain, axes=dd)
            corrs_smag[t0-1] = models.two_point_corr(tau_smag, strmod, axes=dd)
            for ia,alpha in enumerate([2,5,7]):
                aux = load(ddir+'tau_nl.{}.{}.{}.npy'.format(t0,ii,alpha))
                corrs_nonl[ia,t0-1] = models.two_point_corr(aux, strmod, axes=dd)
        save(ddir+'corrs_dns_x_{}.npy'.format(ii),  corrs_dns)
        save(ddir+'corrs_smag_x_{}.npy'.format(ii), corrs_smag)
        save(ddir+'corrs_nonl_x_{}.npy'.format(ii), corrs_nonl)

    # Average over realizations and normalize
    corrs_dns  = -corrs_dns.mean(axis=0)
    corrs_smag = -corrs_smag.mean(axis=0)
    corrs_nonl = -corrs_nonl.mean(axis=1)
    corrs_dns  /= corrs_dns[dps//2]
    corrs_smag /= corrs_smag[cps//2]
    for ia in range(3):
        corrs_nonl[ia] /= corrs_nonl[ia][cps//2]

    # Make symmetric/even
    corrs_dns = models.symmetrize(corrs_dns)
    corrs_smag = models.symmetrize(corrs_smag)
    for ia in range(3):
        corrs_nonl[ia] = models.symmetrize(corrs_nonl[ia])

    figure(y0)
    clf()
    plot(dnxs, corrs_dns,  label='DNS', ls='--')
    plot(xxss, corrs_smag, label=r'$\alpha=1$')
    plot(xxss, corrs_nonl[2], label=r'$\alpha=0.7$')
    plot(xxss, corrs_nonl[1], label=r'$\alpha=0.5$')
    plot(xxss, corrs_nonl[0], label=r'$\alpha=0.2$')
    # save('r_models', xxss)
    # save('r_dns', dnxs)
    # save('corr_dns_{}'.format(ii), corrs_dns)
    # save('corr_smag_{}'.format(ii), corrs_smag)
    # save('corr_nonl_0.7_{}'.format(ii), corrs_nonl[2])
    # save('corr_nonl_0.5_{}'.format(ii), corrs_nonl[1])
    # save('corr_nonl_0.2_{}'.format(ii), corrs_nonl[0])
    xlabel('$r^+$')
    ylabel(r'$C_x(r)/C_x(0)$')
    xlim(-1000,1000)
    if y0!=23:
        yscale('log')
        ylim(bottom=1e-3)
    legend()
    if y0==256:
        title(r'$y^+=1000$')
    else:
        title(r'$y^+={:.1f}$'.format(y_plus[y0]))
    savefig('ms/fig2px{}b'.format(ii))

draw()
show()
