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
num_res = 72

# Direction
dr  = 2*dx
cps = 99
dd  = 2
xxss = arange(-49,50)*dr/visc_length

# Different heights
for y0 in [120,70,40,23,256]:
    if y0==256:
        num_res = 24
    ii = heights[y0]
    try:
        corrs_dns  = load('corrs_dns_x_{}.npy'.format(ii))
        corrs_smag = load('corrs_smag_x_{}.npy'.format(ii))
        corrs_nonl = load('corrs_nonl_x_{}.npy'.format(ii))
    except:
        corrs_dns  = zeros((num_res,cps))
        corrs_smag = zeros((num_res,cps))
        corrs_nonl = zeros((9,num_res,cps))
        for t0 in range(1,num_res+1):

            # Get fields
            strain   = load(ddir+'filt_strain.{}.{}.npy'.format(t0,ii))
            strain   = strain[:,:,bound_x*RD:-bound_x*RD,bound_z*RD:-bound_z*RD]
            tau_dns  = load(ddir+'tau_dns.{}.{}.npy'.format(t0,ii))
            tau_smag = load(ddir+'tau_smag.{}.{}.npy'.format(t0,ii))
            strain   = strain[:,:,::2,::4]
            tau_dns  = tau_dns[:,:,::2,::4]
            tau_smag = tau_smag[:,:,::2,::4]
            corrs_dns[t0-1]  = models.two_point_corr(tau_dns,  strain, axes=dd)
            corrs_smag[t0-1] = models.two_point_corr(tau_smag, strain, axes=dd)
            for ia in range(9):
                aux = load(ddir+'tau_nl.{}.{}.{}.npy'.format(t0,ii,ia+1))[:,:,:,0,:]
                aux = aux[:,:,::2,::4]
                corrs_nonl[ia,t0-1] = models.two_point_corr(aux, strain, axes=dd)
            save('corrs_dns_x_{}.npy'.format(ii),  corrs_dns)
            save('corrs_smag_x_{}.npy'.format(ii), corrs_smag)
            save('corrs_nonl_x_{}.npy'.format(ii), corrs_nonl)

    # Average over realizations and normalize
    corrs_dns  = -corrs_dns.mean(axis=0)
    corrs_smag = -corrs_smag.mean(axis=0)
    corrs_nonl = -corrs_nonl.mean(axis=1)
    corrs_dns  /= corrs_dns[cps//2]
    corrs_smag /= corrs_smag[cps//2]
    for ia in range(9):
        corrs_nonl[ia] /= corrs_nonl[ia][cps//2]

    figure(y0)
    clf()
    plot(xxss, corrs_dns,  label='DNS', ls='--')
    plot(xxss, corrs_smag, label=r'$\alpha=1$')
    plot(xxss, corrs_nonl[6], label=r'$\alpha=0.7$')
    plot(xxss, corrs_nonl[4], label=r'$\alpha=0.5$')
    plot(xxss, corrs_nonl[1], label=r'$\alpha=0.2$')
    xlabel('$r^+$')
    ylabel(r'$C(r)/C(0)$')
    xlim(-1000,1000)
    if y0!=23:
        yscale('log')
        ylim(bottom=1e-3)
    legend()
    if y0==256:
        title(r'$y^+=1000$')
    else:
        title(r'$y^+={:.1f}$'.format(y_plus[y0]))
    savefig('ms/fig2px{}a'.format(ii))

    xlim(-300,300)
    if y0!=23:
        ylim(bottom=1e-2)
    if y0==120:
        gca().get_legend().remove()
        axvspan(-200,200,color='k',alpha=0.1)
        text(0,1e-2,'Non-local kernel volume\n$R=5\Delta$',horizontalalignment='center')
    savefig('ms/fig2px{}b'.format(ii))

draw()
show()
