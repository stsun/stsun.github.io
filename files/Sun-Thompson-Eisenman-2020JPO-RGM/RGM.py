'''
This is a 1.5-layer model for the overturning circulation.
'''
from numba import jit
from config import *

def initialize(nIter0):
    '''
    Here, let's initialize our model
    '''
    if nIter0  == 0:
        # Initialize h
        h   = np.zeros((ny, nx)) + 500.0 * msk
        u   = np.zeros((ny, nx))
        v   = np.zeros((ny, nx))
    else:
        with open("output/pickup.%10.10d.bin"%nIter0, "rb") as fid:
            D = pickle.load(fid)
        h  = D["h"]
        u  = D["u"]
        v  = D["v"]
    return h,u,v


# I tried to parallel the following code using the default prange, but
# it appears not to help. So I will stick to the current form.
@jit(nopython=True)
def obtain_step_uv(u, v, h, ctime):
    '''
    This will return du/dt, dv/dt, dh/dt
    '''
    # Calculate P1 and P2
    dudt = np.zeros((ny, nx))
    dvdt = np.zeros((ny, nx))
    dhdt = np.zeros((ny, nx))
    P   = g * h
    for i in range(1, ny-1):
        for j in range(0, nx):
            ip1   = i+1
            jp1   = j+1
            if jp1 == nx:
                jp1 = 0
            im1   = i-1
            jm1   = j-1

            #--- U:
            if msk[i,j] * msk[i,jp1] >0:     # No Normal flow
                vm = 0.25 * (v[i,j]+v[im1,j]+v[i,jp1]+v[im1,jp1])
                AU = -0.5*(u[i,j]+abs(u[i,j]))*(u[i,j]  - u[i,jm1])/dx[i,j] \
                     -0.5*(u[i,j]-abs(u[i,j]))*(u[i,jp1]- u[i,j]  )/dx[i,j]
                AV = -0.5*(vm+abs(vm))*(u[i,j]-u[im1,j])/dy[i,j]  \
                     -0.5*(vm-abs(vm))*(u[ip1,j]-u[i,j])/dy[i,j]
                CF = f[i,j]*vm
                PG = -(P[i,jp1]-P[i,j])/dx[i,j]
                AH = Ah*(u[i,jp1]+u[i,jm1]-2.0*u[i,j])/dx[i,j]**2 + \
                     Ah*(u[ip1,j]+u[im1,j]-2.0*u[i,j])/dy[i,j]**2
                RR = -rb * u[i,j]
                TAU= taux[i,j]/rho0/max(minH, 0.5*(h[i,j]+h[i,jp1]))

                dudt[i,j] = AU + AV + CF + PG + AH + RR + TAU

            # -- V
            if msk[i,j] * msk[ip1,j] > 0:
                um = 0.25 * (u[i,jm1] + u[i,j] + u[ip1, jm1] + u[ip1,j])
                AU = -0.5*(um+abs(um))*(v[i,j]-v[i,jm1])/dx[i,j] \
                     -0.5*(um-abs(um))*(v[i,jp1]-v[i,j])/dx[i,j]
                AV = -0.5*(v[i,j]+abs(v[i,j]))*(v[i,j]-v[im1,j])/dy[i,j] \
                     -0.5*(v[i,j]-abs(v[i,j]))*(v[ip1,j]-v[i,j])/dy[i,j]
                CF = -f[i,j]*um
                PG = -(P[ip1,j]-P[i,j])/dy[i,j]
                AH = Ah*(v[i,jp1]+v[i,jm1]-2.0*v[i,j])/dx[i,j]**2 + \
                     Ah*(v[ip1,j]+v[im1,j]-2.0*v[i,j])/dy[i,j]**2
                RR = -rb * v[i,j]

                dvdt[i,j] = AU + AV + CF + PG + AH + RR

            # -- H
            if msk[i,j] > 0:
                AU = -(u[i,j]  *dy[i,j]   *(h[i,jp1]+h[i,j])-
                       u[i,jm1]*dy[i,jm1] *(h[i,jm1]+h[i,j]))/2.0
                AV = -(v[i,j]  *dxV[i,j]  *(h[ip1,j]+h[i,j]) -
                       v[im1,j]*dxV[im1,j]*(h[im1,j]+h[i,j]))/2.0
                # AU = -(u[i,j]  * dy[i,j]   * (h[i,jp1]-h[i,j]) +
                #        u[i,jm1]* dy[i,jm1] * (h[i,j]-h[i,jm1]))/2.0
                # AV = -(v[i,j]  * dxV[i,j]   * (h[ip1,j]-h[i,j]) +
                #        v[i,jm1]* dxV[im1,j] * (h[i,j]-h[im1,j]))/2.0
                GMU= Kgm*(msk[i,jp1]*dy[i,j]   *(h[i,jp1]-h[i,j])/dx[i,j] -
                          msk[i,jm1]*dy[i,jm1] *(h[i,j]-h[i,jm1])/dx[i,j])
                GMV= Kgm*(msk[ip1,j]*dxV[i,j]  *(h[ip1,j]-h[i,j])/dy[i,j] -
                          msk[im1,j]*dxV[im1,j]*(h[i,j]-h[im1,j])/dy[i,j])

                if lat[i] > -45 and lat[i]<67:
                    DIA= kappa/max(minH, h[i,j]) * dxdy[i,j]
                else:
                    DIA= 0

                if i < 10: # In the southern margin, a relaxation is applied
                    relax = (minH-h[i,j])/(tau_relax*i)/86400*dxdy[i,j]
                else:
                    relax = 0
                    
                # use relaxation in the Southernmost grid points
                if h[i, j]<minH: 
                    relax += (minH-h[i, j])/strong_relax*dxdy[i,j]

                if T_NADW == 0:
                    WW_NADW = -wNADW[i,j] * dxdy[i,j]
                else:
                    # T_NADW is in years
                    #WW_NADW = -(wNADW[i,j] +
                    #            A_NADW *np.sin(2*PI*time/T_NADW/86400/360))*dxdy[i,j]
                    WW_NADW = -wNADW[i,j]*dxdy[i,j]* \
                              (1+A_NADW/nadw*np.sin(2*PI*ctime/T_NADW/86400/360))
                

                dhdt[i,j] = (AU + AV + GMU + GMV + DIA + relax + WW_NADW)/dxdy[i,j]
                
    return dudt,dvdt,dhdt
                
@jit(nopython=True)                
def obtain_transport_depth(V, H):
    '''
    Here we obtain the transport at 30S
    '''
    kk  = 27    # 30S
    ATL = 0
    PAC = 0

    mk  = msk[kk,:] * msk[kk+1,:]
    VH  = V[kk,:]*dxV[kk,:]*(H[kk,:]+H[kk+1,:])/2.0 * mk
    VG  = Kgm * mk * (H[kk+1,:]-H[kk,:])/dy[kk,0]*dx[kk,0]
    ATL = np.sum(atl[kk+1,:]*(VH-VG))
    PAC = np.sum(pac[kk+1,:]*(VH-VG))
    # for j in range(0, nx):
    #     VH  = V[kk,j]*dxV[kk,j]*(H[kk,j]+H[kk+1,j])/2.0
    #     VG  = Kgm * msk[kk,j]*msk[kk+1,j]*(H[kk+1,j]-H[kk,j])/dy[kk,j] * dx[kk,j]
    #     if atl[kk,j] == 1:
    #         ATL += VH - VG
    #     else:
    #         PAC += VH - VG


    # Calculate the mean depth of ATL and PAC
    h_atl = np.sum(H*dxdy*atl)/A_ATL
    h_pac = np.sum(H*dxdy*pac)/A_PAC
    # Now append it to the files:
    #fid.write("%10.2f %10.2f %10.2f %10.2f \n"%(ATL/1e6, PAC/1e6, h_atl, h_pac))
    #fid.flush()
    return ATL/1e6, PAC/1e6, h_atl, h_pac

@jit(nopython=True)
def use_ab3(uold, ddt):
    unew = uold + dt*(23*ddt[-1,:,:]-16*ddt[-2,:,:]+5*ddt[-3,:,:])/12
    return unew
                
            
                
                    
