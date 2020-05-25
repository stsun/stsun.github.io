#!/central/software/python/3.7.0/bin/python3

# I will use the third-order Adams-Bashforth method to advance the
# momentum and thickness equation

import time
from glob import glob
from RGM import *
from subprocess import call

if restart:
    list_pickup = glob("output/pickup.*.bin")
    list_pickup.sort()
    nIter0      = int(list_pickup[-1].split(".")[-2])
else:
    nIter0      = 0

# Let's initialize the simulations
H  = np.zeros((2, ny, nx))
U  = np.zeros((2, ny, nx))
V  = np.zeros((2, ny, nx))
H[0,:,:],U[0,:,:],V[0,:,:] = initialize(nIter0)

if new_diag:
    fid_diag = open("diagnosis.txt", "w+")
else:
    fid_diag = open("diagnosis.txt", "a+")

nSteps  = int(nyear*360*86400/dt)

# Initialize dudt,dvdt, and dhdt
DUDT = np.zeros((3, ny, nx))
DVDT = np.zeros((3, ny, nx))
DHDT = np.zeros((3, ny, nx))

old  = 0
new  = 1

# Now start the iteration: 1st step
for i,kstep in enumerate(range(nIter0, nIter0+nSteps)):
    ctime = i*dt + PerturbStart
    if i < 2:
        DUDT[i,:,:],DVDT[i,:,:],DHDT[i,:,:] = obtain_step_uv(U[old,:,:], V[old,:,:],
                                                             H[old,:,:], ctime)
        U[new,:,:] = U[old, :, :] + dt * DUDT[i, :, :]
        V[new,:,:] = V[old, :, :] + dt * DVDT[i, :, :]
        H[new,:,:] = H[old, :, :] + dt * DHDT[i, :, :]
    else:
        DUDT[-1,:,:],DVDT[-1,:,:],DHDT[-1,:,:] = obtain_step_uv(U[old,:,:], V[old,:,:],
                                                                H[old,:,:], ctime)
        U[new,:,:] = use_ab3(U[old,:,:], DUDT)
        V[new,:,:] = use_ab3(V[old,:,:], DVDT)
        H[new,:,:] = use_ab3(H[old,:,:], DHDT)

        DUDT = np.roll(DUDT, 2, axis=0)
        DVDT = np.roll(DVDT, 2, axis=0)
        DHDT = np.roll(DHDT, 2, axis=0)


    
    # Now shift U, V, and H
    U = np.roll(U, 1, axis=0)
    V = np.roll(V, 1, axis=0)
    H = np.roll(H, 1, axis=0)
    
    # ATL,PAC  = obtain_transport(V[new_step, :,:], H[new_step,:,:])
    # print("%8.8d: Maximum of h: %f; Mimimum of h: %f ; ATL: %f; PAC: %f"%(kstep, np.nanmax(h),
    #                                                                       np.nanmin(h), ATL, PAC))
    h  = H[old, :, :] + 0.0
    h[h==0] = np.nan
    if np.nanmax(h) > 1e4:
        # send an email to my email account if it fails
        fid = open("email.txt", "w+")
        fid.writelines("Subject: %s failed \n"%runname)
        fid.writelines("\n")
        fid.writelines("%s has failed. Do something!!!"%runname)
        fid.close()
        call("sendmail stsun1989@gmail.com < email.txt", shell=True)
        print("You are dead!!!")
        exit()

    if np.mod(ctime, monFreq) < dt:
        print("%8.8d: max h as %f and min of h as %f"%(kstep, np.nanmax(h), np.nanmin(h)))

    if np.mod(ctime, diagFreq) < dt:
        # diagnose the mean depth, transport, and write to files
        ATL,PAC,h_atl,h_pac = obtain_transport_depth( V[old,:,:],
                                                      H[old,:,:])
        fid_diag.write("%10.5f %10.5f %10.5f %10.5f \n"%(ATL, PAC, h_atl, h_pac))
        fid_diag.flush()

        
    if np.mod(ctime, pickupFreq) < dt:
        # dump pickup files to restart
        D      = {}
        D["h"] = H[old,:,:]
        D["u"] = U[old,:,:]
        D["v"] = V[old,:,:]
        
        with open("output/pickup.%10.10d.bin"%kstep, "wb") as fid:
            pickle.dump(D, fid)

    # The following code will cause the memory usage to
    # increase. Don't turn this on for long-term runs.
    
    # if np.mod(time, figFreq) == 0:
    #     fig = plt.figure(figsize=(12, 6))
    #     ax1 = fig.add_subplot(121)
    #     ax2 = fig.add_subplot(122)
    #     #lev = np.arange(0, 3000, 300)

    #     #h[h==0] = np.nan
    #     ax1.contourf(lon, lat, h, levels=np.arange(50,500, 10), extend="both", cmap="jet")
    #     ax2.contourf(lon, lat, U[2,:,:]) #, levels=lev, extend="both")

    #     #for ax in [ax1, ax2]:
    #     #    ax.patch.set_color("0.9")
    #     plt.savefig("Figures/Depth-%6.6d.png"%i, bbox_inches="tight")

fid.close()
# fid = open("email.txt", "w+")
# fid.writelines("Subject: %s is finished \n"%runname)
# fid.writelines("\n")
# fid.writelines("%s is finished."%runname)
# fid.close()
# call("sendmail stsun1989@gmail.com < email.txt", shell=True)
print("Congratulations!!!")

if resubmit > 0:
    # If resubmitting, I will continue the time counting 
    call('sed -i "/PerturbStart/c\\PerturbStart = %f" config.py'%(ctime),
         shell=True)
    # modify resubmit in config.py
    call('sed -i "/resubmit/c\\resubmit = %d" config.py'%(resubmit-1),
         shell=True)
    # make sure restart is true
    call('sed -i "/restart   =/c\\restart   =True" config.py', shell=True)
    call("sbatch submit.sh", shell=True)

    
