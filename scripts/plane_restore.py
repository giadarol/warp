"""PlaneRestore class that is used to restore particle and field data
at a specified z plane, data saved by PlaneSave.
The two simulations are linked together.
"""

__all__ = ['PlaneRestore']

from warp import *
import cPickle

class PlaneRestore:
  """
Saves the particle data and phi to a file just after the
field solve is finished. The positions and velocities are staggered.
It is automatically called every time step after the field solve
Input:
  - filename=runid.plane: filename where data is stored
  - zplane: location where simulations are linked together. Units of meters
            relative to the lab frame. Defaults to w3d.zmmin.
  - js: species which are saved. Defaults to all species. Can be single
        integer or a list of integers.
  - l_restore_phi=1: flag for restoring phi or not.
  - lrestoreparticles=1: flag for restoring particles
  - starttime=None: If specified, the time at which to start the simulation.
                    This can be used to skip part of the saved data, or to
                    start at an earlier time before saved data is available.
  - verbose=False: When True, prints messages about what is happening
  """

  def __init__(self,filename,zplane=None,js=None,
               l_restore_phi=1,lrestoreparticles=1,starttime=None,verbose=False):

    # --- Save some input values
    self.filename = filename
    self.zplane = zplane
    self.js = js
    self.l_restore_phi = l_restore_phi
    self.lrestoreparticles = lrestoreparticles
    self.starttime = starttime
    self.verbose = verbose

    # --- Install the routines that do the work.
    installuserinjection(self.restoreparticles)
    installbeforefs(self.restoreplane_bfs)
    installafterfs(self.restoreplane_afs)

    # --- Initialize self.f to None to flag that the initialization
    # --- has not yet happened.
    self.f = None

  def read(self):
    if self.f is None:
      self.f = open(self.filename,'r')

    return cPickle.load(self.f)

  def readnextstep(self):
    self.data = {}
    while True:
      name,val = self.read()
      self.data[name] = val
      if name.startswith('phiplane'): break

    # --- name will be 'phiplane%09d'%it
    self.data['it'] = int(name[8:])

    if self.verbose: print "PlaneRestore: read in data from step %d"%self.data['it']

  def initrestoreplane(self):

    if self.zplane is None: self.zplane = w3d.zmmin

    # --- Read in the initial data
    initdata = {}
    while True:
      name,val = self.read()
      initdata[name] = val
      if name == 'solvergeom': break

    self.readnextstep()

    self.zshift = self.zplane - initdata['zplane']

    self.lsavephi = initdata['lsavephi']
    self.lsaveparticles = initdata['lsaveparticles']
    if not self.lsavephi: self.l_restore_phi = false
    if not self.lsaveparticles: self.lrestoreparticles = false

    # --- get time level of first plane and subtract 1
    self.it_restore = -1

    # --- get time step, tmin, tmax
    self.dt = initdata['dt']
    self.tmin = initdata['tmin']
    top.time = self.tmin
    self.time_restore = self.tmin

    if self.lrestoreparticles:
      # --- initializes list of species
      if self.js is None:
        self.jslist = range(top.ns)
      else:
        try:
          list(self.js)
          self.jslist = self.js
        except TypeError:
          self.jslist= [self.js]

      # --- restore particle charge, mass, weight
      for js in self.jslist:
        top.pgroup.sq[js] = initdata['sq_%d'%js]
        top.pgroup.sm[js] = initdata['sm_%d'%js]
        top.pgroup.sw[js] = initdata['sw_%d'%js]

      # --- make sure that pid will be allocated
      #top.npid = initdata['npid']
      #setuppgroup(top.pgroup)

    if self.l_restore_phi:
      # --- restore solver geometry of the saved data
      try:
        self.solvergeom = initdata['solvergeom']
      except:
        self.solvergeom = w3d.XYZgeom

      # set up indices which specify transverse extent of saved and restored phi
      # _r for restored phi array, _s for saved phi array
      # '0' is minimum index, 'm' is maximum index

      self.sym_plane = initdata['sym_plane']
      ixa_plane = initdata['ixa_plane']
      iya_plane = initdata['iya_plane']
      self.nx_plane = initdata['nx_plane']
      self.ny_plane = initdata['ny_plane']
      self.xmmin = initdata['xmmin']
      self.xmmax = initdata['xmmax']
      self.ymmin = initdata['ymmin']
      self.ymmax = initdata['ymmax']

      #self.nx0_r = max(0,int(floor((self.xmmin - w3d.xmmin)/w3d.dx)))
      #self.ny0_r = max(0,int(floor((self.ymmin - w3d.ymmin)/w3d.dy)))
      #self.nxm_r = min(w3d.nx,int(floor((self.xmmax - w3d.xmmin)/w3d.dx)))
      #self.nym_r = min(w3d.ny,int(floor((self.ymmax - w3d.ymmin)/w3d.dy)))

      self.nx0_r = max(0, 0 - ixa_plane + w3d.ix_axis)
      self.ny0_r = max(0, 0 - iya_plane + w3d.iy_axis)
      self.nxm_r = min(w3d.nx, self.nx_plane - ixa_plane + w3d.ix_axis)
      self.nym_r = min(w3d.ny, self.ny_plane - iya_plane + w3d.iy_axis)
      self.nx0_s = self.nx0_r - w3d.ix_axis + ixa_plane
      self.ny0_s = self.ny0_r - w3d.iy_axis + iya_plane
      self.nxm_s = self.nxm_r - w3d.ix_axis + ixa_plane
      self.nym_s = self.nym_r - w3d.iy_axis + iya_plane

      # --- deal with symmetries
      # --- if saved is 2 or 4 fold symmetric and restored isn't,
      # --- lower half of restored is filled with inverted saved phi
      if ((self.sym_plane == 2 and (not w3d.l2symtry and not w3d.l4symtry)) or
          (self.sym_plane == 4 and (not w3d.l2symtry and not w3d.l4symtry))):
        self.ny0_r2 = max(0, - self.ny_plane - iya_plane + w3d.iy_axis)
        self.nym_r2 = min(w3d.ny, 0 - iya_plane + w3d.iy_axis)
        self.ny0_s2 = - self.ny0_r + w3d.iy_axis + iya_plane
        self.nym_s2 =   self.nym_r - w3d.iy_axis + iya_plane
      if ((self.sym_plane == 4 and (not w3d.l2symtry and not w3d.l4symtry)) or
          (self.sym_plane == 4 and (    w3d.l2symtry and not w3d.l4symtry))):
        self.nx0_r2 = max(0, - self.nx_plane - ixa_plane + w3d.ix_axis)
        self.nxm_r2 = min(w3d.nx, 0 - ixa_plane + w3d.ix_axis)
        self.nx0_s2 = self.nxm_r - w3d.ix_axis + ixa_plane
        self.nxm_s2 = - self.nx0_r + w3d.ix_axis + ixa_plane

    # --- Reset the time to the start time if specified.
    if self.starttime is not None:
      top.time = self.starttime
      while self.time_restore+self.dt <= top.time:
        # --- increment the timelevel of the plane
        self.it_restore += 1
        self.time_restore += self.dt
      # --- Setup phi at the start time
      self.restoreplane_bfs()

    if self.verbose:
      print "PlaneRestore: initial data"
      print "  File",self.filename
      print "  Restoring phi",self.lsavephi
      print "  Restoring particles",self.lsaveparticles
      print "  Start time",self.tmin
      print "  Time step",self.dt

  ###########################################################################
  def disable_plane_restore(self):
    # for some reason, does not work!
    uninstalluserinjection(self.restoreparticles)
    uninstallbeforefs(self.restoreplane_bfs)
    uninstallafterfs(self.restoreplane_afs)

  def jumptotime(self,time):
    """Jump to the specified time and set the phi boundary condition.
No particles are loaded."""

    # --- Do the initialization if it hasn't been done yet.
    if self.f is None:
      self.initrestoreplane()

    # --- Set the time to the desired time.
    top.time = time

    while self.time_restore+self.dt <= top.time:
      # --- increment the timelevel of the plane
      self.it_restore += 1
      self.time_restore += self.dt

    # --- restore phi only if between grid bounds
    if (self.zplane < w3d.zmmin+top.zbeam or
        self.zplane+top.zbeam >= w3d.zmmax): return

    # --- calculate grid location of new_plane
    iz = nint((self.zplane - top.zbeam - w3d.zmmin)/w3d.dz)

    # --- load saved phi into the phi array
    self.restore_phi(iz,self.it_restore)

  ###########################################################################
  def restoreparticles(self):
    "Restore the particles"
    # --- Do the initialization if it hasn't been done yet.
    if self.f is None:
      self.initrestoreplane()

    # --- restore only if between grid bounds
    if (self.zplane < w3d.zmmin+top.zbeam or
        self.zplane+top.zbeam >= w3d.zmmax): return

    # --- Loop over restored data, restoring the data up to the current
    # --- time level of the simulation. This allows the stored data dt
    # --- to be different than the current simulation dt.
    while self.time_restore <= top.time:

      # --- increment the timelevel of the plane
      self.it_restore += 1
      self.time_restore += self.dt

      # --- load particles for each species
      for js in self.jslist:
        self.restoreparticlespecies(js,self.it_restore)

  def restoreparticlespecies(self,js=0,it=0):
    if not self.lrestoreparticles: return

    while self.data['it'] < it:
      self.readnextstep()

    # --- Apparently, no data was written for this step, so do nothing
    if self.data['it'] > it:
      if self.verbose: print "PlaneRestore: no data for step",it
      return

    # --- put restored data into particle arrays, adjusting the z location
    suffix = '%09d_%d'%(it,js)

    # --- Check if data was written for this step.
    if 'xp'+suffix not in self.data:
      if self.verbose: print "PlaneRestore: no particle data for step",it
      return

    xx = self.data['xp'+suffix]
    yy = self.data['yp'+suffix]
    zz = self.data['zp'+suffix] + self.zshift
    ux = self.data['uxp'+suffix]
    uy = self.data['uyp'+suffix]
    uz = self.data['uzp'+suffix]
    gi = self.data['gaminv'+suffix]
    id = self.data['pid'+suffix]

    # --- Do some fudging to get the shape of pid correct. This is not
    # --- perfect, since it will may munge the data in pid if things
    # --- are arranged differently.
    if id.shape[1] < top.npid:
      newid = fzeros((id.shape[0],top.npid),'d')
      newid[:,:id.shape[1]] = id
      id = newid
    elif id.shape[1] > top.npid:
      id = id[:,:top.npid]

    # --- Check if particles are being added out of bounds
    if min(zz) < top.zpmin+top.zbeam or max(zz) > top.zpmax+top.zbeam:
      print "PlaneRestore: restored particles are out of bounds."
      print "\nThe extent of the simulation is %f to %f"%(top.zpmin+top.zbeam,top.zpmax+top.zbeam)
      print "The extent of the restored particles is %f to %f\n"%(min(zz),max(zz))
      raise Exception("PlaneRestore: restored particles are out of bounds.")

    # --- Note that all processors read in the data, but only particles
    # --- within the processors domain are added.
    if self.verbose: print "PlaneRestore: Restoring %d particles on step %d"%(len(xx),it)
    addparticles(xx,yy,zz,ux,uy,uz,gi,id,
                 js=js,
                 lallindomain=false,
                 lmomentum=true,
                 resetrho=false)

  ###########################################################################
  # --- restore the next plane of data
  def restoreplane_bfs(self):

    # --- Do the initialization if it hasn't been done yet.
    if self.f is None:
      self.initrestoreplane()

    # --- restore only if between grid bounds
    if (self.zplane < w3d.zmmin+top.zbeam or
        self.zplane+top.zbeam >= w3d.zmmax): return

    # --- calculate grid location of new_plane
    iz = nint((self.zplane - top.zbeam - w3d.zmmin)/w3d.dz)

    # --- load saved phi into the phi array
    self.restore_phi(iz,self.it_restore)

    if self.verbose: print "PlaneRestore: Restoring phi on step %d"%(len(xx),it)

  def restoreplane_afs(self):
    # --- this routine resets the potential at the plane iz=-1 after the
    # --- field solve if this is needed

    # --- restore only if between grid bounds
    if (self.zplane < w3d.zmmin+top.zbeam or
        self.zplane+top.zbeam >= w3d.zmmax): return

    # --- calculate grid location of new_plane
    iz = nint((self.zplane - top.zbeam - w3d.zmmin)/w3d.dz)

    # --- reset phi at plane iz=-1 if zplane is at iz=0
    if (iz == 0):
      self.restore_phi(iz,self.it_restore)

  #######################################################################
  def restore_phi(self,iz,it):
    # --- return if flag indicates phi not to be restored
    if self.l_restore_phi is 0: return

    while self.data['it'] < it:
      self.readnextstep()

    # --- Apparently, no data was written for this step, so do nothing
    if self.data['it'] > it: return

    # --- Read in the phi data if it is available.
    if 'phiplane%09d'%it not in self.data: return

    savedphi = self.data['phiplane%09d'%it]

    solver = getregisteredsolver()
    if solver is None: solver = w3d

    if self.solvergeom == solver.solvergeom and solver.solvergeom == w3d.XYZgeom:
      self.restore_phi_3d_to_3d(iz-top.izfsslave[me],it,savedphi,solver.phi,
                                solver)
      if top.izpslave[me] != top.izfsslave[me]:
        # --- This is not really correct, since phip will have a different
        # --- shape and phi, so the dimensions should be passed in too.
        self.restore_phi_3d_to_3d(iz-top.izpslave[me],it,savedphi,solver.phip,
                                  solver)
    elif self.solvergeom == solver.solvergeom and solver.solvergeom == w3d.RZgeom:
      self.restore_phi_rz_to_rz(iz-top.izfsslave[me],it,savedphi,solver.phi,
                                solver)
    elif self.solvergeom == w3d.RZgeom and solver.solvergeom == w3d.XYZgeom:
      self.restore_phi_rz_to_3d(iz,it,savedphi,solver.phi)
      if top.izpslave[me] != top.izfsslave[me]:
        self.restore_phi_rz_to_3d(iz,it,savedphi,solver.phip)

  #######################################################################
  # This routine copies the saved phi plane into the current phi array
  # making use of different numbers of grid cells and differing symmetries.
  # Both saved and restored phi are 3-D.
  def restore_phi_3d_to_3d(self,iz,it,savedphi,phi,solver):
    if iz < 0 or iz > solver.nzlocal: return
    for i in range(2):
      grid2grid(phi[self.nx0_r:self.nxm_r+1,self.ny0_r:self.nym_r+1,iz+i],
                solver.nx,solver.ny,
                solver.xmmin,solver.xmmax,solver.ymmin,solver.ymmax,
                savedphi[...,i],self.nx_plane,self.ny_plane,
                self.xmmin,self.xmmax,self.ymmin,self.ymmax)

    if ((self.sym_plane == 2 and (not solver.l2symtry and not solver.l4symtry)) or
        (self.sym_plane == 4 and (not solver.l2symtry and not solver.l4symtry))):
    #     phi(self.nx0_r:self.nxm_r,self.ny0_r2:self.nym_r2,iz-1:iz)=
    #      phi_plane(nx0_s:nxm_s,nym_s2:ny0_s2:-1,)
      for i in range(2):
        grid2grid(phi[self.nx0_r:self.nxm_r,self.ny0_r2:self.nym_r2,iz+i],
                  solver.nx,solver.ny,
                  solver.xmmin,solver.xmmax,solver.ymmin,solver.ymmax,
                  savedphi[...,i],self.nx_plane,self.ny_plane,
                  self.xmmin,self.xmmax,self.ymmin,self.ymmax)

    if ((self.sym_plane == 4 and ( solver.l2symtry and not solver.l4symtry))):
    #  phi(nx0_r2:nxm_r2,ny0_r:nym_r,iz-1:iz)=
    #    phi_plane(nx0_s2:nxm_s2:-1,ny0_s:nym_s,)
      for i in range(2):
        grid2grid(phi[self.nx0_r:self.nxm_r,self.ny0_r2:self.nym_r2,iz+i],
                  solver.nx,solver.ny,
                  solver.xmmin,solver.xmmax,solver.ymmin,solver.ymmax,
                  savedphi[...,i],self.nx_plane,self.ny_plane,
                  self.xmmin,self.xmmax,self.ymmin,self.ymmax)

    if ((self.sym_plane == 4 and (not solver.l2symtry and not solver.l4symtry))):
    #  phi(nx0_r2:nxm_r2,ny0_r2:nym_r,iz-1:iz)=
    #    phi_plane(nx0_s2:nxm_s2:-1,ny0_s2:nym_s,)
      for i in range(2):
        grid2grid(phi[self.nx0_r:self.nxm_r,self.ny0_r2:self.nym_r2,iz+i],
                  solver.nx,solver.ny,
                  solver.xmmin,solver.xmmax,solver.ymmin,solver.ymmax,
                  savedphi[...,i],self.nx_plane,self.ny_plane,
                  self.xmmin,self.xmmax,self.ymmin,self.ymmax)

  def restore_phi_rz_to_rz(self,iz,it,savedphi,phi,solver):
    if iz < 0 or iz > solver.nzlocal: return
    # --- For now, this assumes that the arrays are the same shape.
    phi[1:-1,0,iz] = savedphi[:,0]
    phi[1:-1,0,iz+1] = savedphi[:,1]

  #######################################################################
  # This routine copies the saved phi plane into the current phi array
  # making use of different numbers of grid cells and differing symmetries.
  # Saved phi is rz and restored phi is 3-D.
  def restore_phi_rz_to_3d(self,iz,it,savedphi,phi):
    if iz < 0 or iz > w3d.nzlocal: return
    try:
      self._rz_to_3d_inited
    except:
      self._rz_to_3d_inited = 1
      xmmin = w3d.xmmin + w3d.dx*self.nx0_r
      nx = self.nxm_r - self.nx0_r
      ymmin = w3d.ymmin + w3d.dy*self.ny0_r
      ny = self.nym_r - self.ny0_r
      xmesh,ymesh = getmesh2d(xmmin,w3d.dx,nx,ymmin,w3d.dy,ny)
      rmesh = sqrt(xmesh**2 + ymesh**2)
      dr = (self.xmmax - self.xmmin)/self.nx_plane
      self.irmesh = int(rmesh/dr)
      self.wrmesh =     rmesh/dr  - self.irmesh
      self.wrmesh = where(self.irmesh >= self.nx_plane,1,self.wrmesh)
      self.irmesh = where(self.irmesh >= self.nx_plane,self.nx_plane-1,self.irmesh)

    i1 = self.nx0_r
    i2 = self.nxm_r+1
    j1 = self.ny0_r
    j2 = self.nym_r+1
    for i in range(2):
      phi[i1:i2,j1:j2,iz+i] = (
          take(savedphi[:,0,i],self.irmesh  )*(1.-self.wrmesh) +
          take(savedphi[:,0,i],self.irmesh+1)*self.wrmesh)

#============================================================================
#============================================================================
#============================================================================
#============================================================================
class PlaneRestoreOriginal:
  """
OBSOLETE
Use only when reading in old PlaneSave datafiles that used PyPDB format.
Saves the particle data and phi to a file just after the
field solve is finished. The positions and velocities are staggered.
It is automatically called every time step after the field solve
Input:
  - filename=runid.plane: filename where data is stored
  - zplane: location where simulations are linked together. Units of meters
            relative to the lab frame. Defaults to w3d.zmmin.
  - js: species which are saved. Defaults to all species. Can be single
        integer or a list of integers.
  - l_restore_phi=1: flag for restoring phi or not.
  - lrestoreparticles=1: flag for restoring particles
  - starttime=None: If specified, the time at which to start the simulation.
                    This can be used to skip part of the saved data, or to
                    start at an earlier time before saved data is available.
  """

  def __init__(self,filename,zplane=None,js=None,
               l_restore_phi=1,lrestoreparticles=1,starttime=None):

    # --- Save some input values
    self.filename = filename
    self.zplane = zplane
    self.js = js
    self.l_restore_phi = l_restore_phi
    self.lrestoreparticles = lrestoreparticles
    self.starttime = starttime

    # --- Install the routines that do the work.
    installuserinjection(self.restoreparticles)
    installbeforefs(self.restoreplane_bfs)
    installafterfs(self.restoreplane_afs)

    # --- Initialize self.f to None to flag that the initialization
    # --- has not yet happened.
    self.f = None

  def initrestoreplane(self):

    ##############################
    # --- Initialization stuff

    if self.zplane is None: self.zplane = w3d.zmmin

    # --- open the file which holds the data
    self.f = PR.PR(self.filename)
    self.zshift = self.zplane - self.f.zplane

    self.lsavephi = self.f.lsavephi
    self.lsaveparticles = self.f.lsaveparticles
    if not self.lsavephi: self.l_restore_phi = false
    if not self.lsaveparticles: self.lrestoreparticles = false

    # --- get time level of first plane and subtract 1
    self.it_restore = 0

    # --- get time step, tmin, tmax
    self.dt = self.f.dt
    self.tmin = self.f.tmin
    top.time = self.tmin
    self.tmax = self.f.tmax
    self.time_restore = self.tmin

    if self.lrestoreparticles:
      # --- initializes list of species
      if self.js is None:
        self.jslist = range(top.ns)
      else:
        try:
          list(self.js)
          self.jslist = self.js
        except TypeError:
          self.jslist= [self.js]

      # --- restore particle charge, mass, weight
      for js in self.jslist:
        top.pgroup.sq[js] = self.f.read('sq_%d'%js)
        top.pgroup.sm[js] = self.f.read('sm_%d'%js)
        top.pgroup.sw[js] = self.f.read('sw_%d'%js)

      # --- make sure that pid will be allocated
      #top.npid = self.f.npid
      #setuppgroup(top.pgroup)

    if self.l_restore_phi:
      # --- restore solver geometry of the saved data
      try:
        self.solvergeom = self.f.solvergeom
      except:
        self.solvergeom = w3d.XYZgeom

      # set up indices which specify transverse extent of saved and restored phi
      # _r for restored phi array, _s for saved phi array
      # '0' is minimum index, 'm' is maximum index

      #self.nx0_r = max(0,int(floor((self.f.xmmin - w3d.xmmin)/w3d.dx)))
      #self.ny0_r = max(0,int(floor((self.f.ymmin - w3d.ymmin)/w3d.dy)))
      #self.nxm_r = min(w3d.nx,int(floor((self.f.xmmax - w3d.xmmin)/w3d.dx)))
      #self.nym_r = min(w3d.ny,int(floor((self.f.ymmax - w3d.ymmin)/w3d.dy)))

      self.nx0_r = max(0, 0 - self.f.ixa_plane + w3d.ix_axis)
      self.ny0_r = max(0, 0 - self.f.iya_plane + w3d.iy_axis)
      self.nxm_r = min(w3d.nx, self.f.nx_plane - self.f.ixa_plane + w3d.ix_axis)
      self.nym_r = min(w3d.ny, self.f.ny_plane - self.f.iya_plane + w3d.iy_axis)
      self.nx0_s = self.nx0_r - w3d.ix_axis + self.f.ixa_plane
      self.ny0_s = self.ny0_r - w3d.iy_axis + self.f.iya_plane
      self.nxm_s = self.nxm_r - w3d.ix_axis + self.f.ixa_plane
      self.nym_s = self.nym_r - w3d.iy_axis + self.f.iya_plane

      # --- deal with symmetries
      # --- if saved is 2 or 4 fold symmetric and restored isn't,
      # --- lower half of restored is filled with inverted saved phi
      if ((self.f.sym_plane == 2 and (not w3d.l2symtry and not w3d.l4symtry)) or
          (self.f.sym_plane == 4 and (not w3d.l2symtry and not w3d.l4symtry))):
        self.ny0_r2 = max(0, - self.f.ny_plane - self.f.iya_plane + w3d.iy_axis)
        self.nym_r2 = min(w3d.ny, 0 - self.f.iya_plane + w3d.iy_axis)
        self.ny0_s2 = - self.ny0_r + w3d.iy_axis + self.f.iya_plane
        self.nym_s2 =   self.nym_r - w3d.iy_axis + self.f.iya_plane
      if ((self.f.sym_plane == 4 and (not w3d.l2symtry and not w3d.l4symtry)) or
          (self.f.sym_plane == 4 and (    w3d.l2symtry and not w3d.l4symtry))):
        self.nx0_r2 = max(0, - self.f.nx_plane - self.f.ixa_plane + w3d.ix_axis)
        self.nxm_r2 = min(w3d.nx, 0 - self.f.ixa_plane + w3d.ix_axis)
        self.nx0_s2 = self.nxm_r - w3d.ix_axis + self.f.ixa_plane
        self.nxm_s2 = - self.nx0_r + w3d.ix_axis + self.f.ixa_plane

    # --- Reset the time to the start time if specified.
    if self.starttime is not None:
      top.time = self.starttime
      while self.time_restore <= top.time:
        # --- increment the timelevel of the plane
        self.it_restore += 1
        self.time_restore += self.dt
      # --- Setup phi at the start time
      self.restoreplane_bfs()

  ###########################################################################
  def disable_plane_restore(self):
    # for some reason, does not work!
    uninstalluserinjection(self.restoreparticles)
    uninstallbeforefs(self.restoreplane_bfs)
    uninstallafterfs(self.restoreplane_afs)

  def jumptotime(self,time):
    """Jump to the specified time and set the phi boundary condition.
No particles are loaded."""

    # --- Do the initialization if it hasn't been done yet.
    if self.f is None:
      self.initrestoreplane()

    # --- Set the time to the desired time.
    top.time = time

    while self.time_restore <= top.time:
      # --- increment the timelevel of the plane
      self.it_restore += 1
      self.time_restore += self.dt

    # --- restore phi only if between grid bounds
    if (self.zplane < w3d.zmmin+top.zbeam or
        self.zplane+top.zbeam >= w3d.zmmax): return

    # --- calculate grid location of new_plane
    iz = nint((self.zplane - top.zbeam - w3d.zmmin)/w3d.dz)

    # --- load saved phi into the phi array
    self.restore_phi(iz,self.it_restore)


  ###########################################################################
  def restoreparticles(self):
    "Restore the particles"
    # --- Do the initialization if it hasn't been done yet.
    if self.f is None:
      self.initrestoreplane()

    # --- restore only if between grid bounds
    if (self.zplane < w3d.zmmin+top.zbeam or
        self.zplane+top.zbeam >= w3d.zmmax): return

    # --- Loop over restored data, restoring the data up to the current
    # --- time level of the simulation. This allows the stored data dt
    # --- to be different than the current simulation dt.
    while self.time_restore <= top.time:
      # --- increment the timelevel of the plane
      self.it_restore += 1
      self.time_restore += self.dt
      it = self.it_restore

      # --- load particles for each species
      for js in self.jslist:
        self.restoreparticlespecies(js,it)

  def restoreparticlespecies(self,js=0,it=0):
    if not self.lrestoreparticles: return

    # --- put restored data into particle arrays, adjusting the z location
    # --- Check if data was written for this step.
    suffix = '%09d_%d'%(it,js)
    oldsuffix = '%d_%d'%(it,js)
    if ('xp'+suffix in self.f.inquire_names()
         or 'xp'+oldsuffix in self.f.inquire_names()):
      try:
        xx = self.f.read('xp'+suffix)
      except:
        suffix = oldsuffix
        xx = self.f.read('xp'+suffix)
      yy = self.f.read('yp'+suffix)
      zz = self.f.read('zp'+suffix)+self.zshift
      ux = self.f.read('uxp'+suffix)
      uy = self.f.read('uyp'+suffix)
      uz = self.f.read('uzp'+suffix)
      gi = self.f.read('gaminv'+suffix)
      id = self.f.read('pid'+suffix)
      # --- Do some fudging to get the shape of pid correct. This is not
      # --- perfect, since it will may munge the data in pid if things
      # --- are arranged differently.
      if id.shape[1] < top.npid:
        newid = fzeros((id.shape[0],top.npid),'d')
        newid[:,:id.shape[1]] = id
        id = newid
      elif id.shape[1] > top.npid:
        id = id[:,:top.npid]
      addparticles(xx,yy,zz,ux,uy,uz,gi,id,
                   js=js,
                   lallindomain=false,
                   lmomentum=true,
                   resetrho=false)

  ###########################################################################
  # --- restore the next plane of data
  def restoreplane_bfs(self):

    # --- Do the initialization if it hasn't been done yet.
    if self.f is None:
      self.initrestoreplane()

    # --- restore only if between grid bounds
    if (self.zplane < w3d.zmmin+top.zbeam or
        self.zplane+top.zbeam >= w3d.zmmax): return

    # --- calculate grid location of new_plane
    iz = nint((self.zplane - top.zbeam - w3d.zmmin)/w3d.dz)

    # --- load saved phi into the phi array
    self.restore_phi(iz,self.it_restore)

  def restoreplane_afs(self):
    # --- this routine resets the potential at the plane iz=-1 after the
    # --- field solve if this is needed

    # --- restore only if between grid bounds
    if (self.zplane < w3d.zmmin+top.zbeam or
        self.zplane+top.zbeam >= w3d.zmmax): return

    # --- calculate grid location of new_plane
    iz = nint((self.zplane - top.zbeam - w3d.zmmin)/w3d.dz)

    # --- reset phi at plane iz=-1 if zplane is at iz=0
    if (iz == 0):
      self.restore_phi(iz,self.it_restore)

  #######################################################################
  def restore_phi(self,iz,it):
    # --- return if flag indicates phi not to be restored
    if self.l_restore_phi is 0: return

    # --- Read in the phi data if it is available. The second read is looking
    # --- for the old suffix format.
    savedphi = None
    if 'phiplane%09d'%it in self.f.inquire_names():
      savedphi = self.f.read('phiplane%09d'%it)
    elif 'phiplane%d'%it in self.f.inquire_names():
      savedphi = self.f.read('phiplane%d'%it)
    if savedphi is None: return

    solver = getregisteredsolver()
    if solver is None: solver = w3d

    if self.solvergeom == solver.solvergeom and solver.solvergeom == w3d.XYZgeom:
      self.restore_phi_3d_to_3d(iz-top.izfsslave[me],it,savedphi,solver.phi,
                                solver)
      if top.izpslave[me] != top.izfsslave[me]:
        # --- This is not really correct, since phip will have a different
        # --- shape and phi, so the dimensions should be passed in too.
        self.restore_phi_3d_to_3d(iz-top.izpslave[me],it,savedphi,solver.phip,
                                  solver)
    elif self.solvergeom == solver.solvergeom and solver.solvergeom == w3d.RZgeom:
      self.restore_phi_rz_to_rz(iz-top.izfsslave[me],it,savedphi,solver.phi,
                                solver)
    elif self.solvergeom == w3d.RZgeom and solver.solvergeom == w3d.XYZgeom:
      self.restore_phi_rz_to_3d(iz,it,savedphi,solver.phi)
      if top.izpslave[me] != top.izfsslave[me]:
        self.restore_phi_rz_to_3d(iz,it,savedphi,solver.phip)

  #######################################################################
  # This routine copies the saved phi plane into the current phi array
  # making use of different numbers of grid cells and differing symmetries.
  # Both saved and restored phi are 3-D.
  def restore_phi_3d_to_3d(self,iz,it,savedphi,phi,solver):
    if iz < 0 or iz > solver.nzlocal: return
    for i in range(2):
      grid2grid(phi[self.nx0_r:self.nxm_r+1,self.ny0_r:self.nym_r+1,iz+i],
                solver.nx,solver.ny,
                solver.xmmin,solver.xmmax,solver.ymmin,solver.ymmax,
                savedphi[...,i],self.f.nx_plane,self.f.ny_plane,
                self.f.xmmin,self.f.xmmax,self.f.ymmin,self.f.ymmax)

    if ((self.f.sym_plane == 2 and (not solver.l2symtry and not solver.l4symtry)) or
        (self.f.sym_plane == 4 and (not solver.l2symtry and not solver.l4symtry))):
    #     phi(self.nx0_r:self.nxm_r,self.ny0_r2:self.nym_r2,iz-1:iz)=
    #      self.f.phi_plane(nx0_s:nxm_s,nym_s2:ny0_s2:-1,)
      for i in range(2):
        grid2grid(phi[self.nx0_r:self.nxm_r,self.ny0_r2:self.nym_r2,iz+i],
                  solver.nx,solver.ny,
                  solver.xmmin,solver.xmmax,solver.ymmin,solver.ymmax,
                  savedphi[...,i],self.f.nx_plane,self.f.ny_plane,
                  self.f.xmmin,self.f.xmmax,self.f.ymmin,self.f.ymmax)

    if ((self.f.sym_plane == 4 and ( solver.l2symtry and not solver.l4symtry))):
    #  phi(nx0_r2:nxm_r2,ny0_r:nym_r,iz-1:iz)=
    #    self.f.phi_plane(nx0_s2:nxm_s2:-1,ny0_s:nym_s,)
      for i in range(2):
        grid2grid(phi[self.nx0_r:self.nxm_r,self.ny0_r2:self.nym_r2,iz+i],
                  solver.nx,solver.ny,
                  solver.xmmin,solver.xmmax,solver.ymmin,solver.ymmax,
                  savedphi[...,i],self.f.nx_plane,self.f.ny_plane,
                  self.f.xmmin,self.f.xmmax,self.f.ymmin,self.f.ymmax)

    if ((self.f.sym_plane == 4 and (not solver.l2symtry and not solver.l4symtry))):
    #  phi(nx0_r2:nxm_r2,ny0_r2:nym_r,iz-1:iz)=
    #    self.f.phi_plane(nx0_s2:nxm_s2:-1,ny0_s2:nym_s,)
      for i in range(2):
        grid2grid(phi[self.nx0_r:self.nxm_r,self.ny0_r2:self.nym_r2,iz+i],
                  solver.nx,solver.ny,
                  solver.xmmin,solver.xmmax,solver.ymmin,solver.ymmax,
                  savedphi[...,i],self.f.nx_plane,self.f.ny_plane,
                  self.f.xmmin,self.f.xmmax,self.f.ymmin,self.f.ymmax)

  def restore_phi_rz_to_rz(self,iz,it,savedphi,phi,solver):
    if iz < 0 or iz > solver.nzlocal: return
    # --- For now, this assumes that the arrays are the same shape.
    phi[1:-1,0,iz] = savedphi[:,0]
    phi[1:-1,0,iz+1] = savedphi[:,1]

  #######################################################################
  # This routine copies the saved phi plane into the current phi array
  # making use of different numbers of grid cells and differing symmetries.
  # Saved phi is rz and restored phi is 3-D.
  def restore_phi_rz_to_3d(self,iz,it,savedphi,phi):
    if iz < 0 or iz > w3d.nzlocal: return
    try:
      self._rz_to_3d_inited
    except:
      self._rz_to_3d_inited = 1
      xmmin = w3d.xmmin + w3d.dx*self.nx0_r
      nx = self.nxm_r - self.nx0_r
      ymmin = w3d.ymmin + w3d.dy*self.ny0_r
      ny = self.nym_r - self.ny0_r
      xmesh,ymesh = getmesh2d(xmmin,w3d.dx,nx,ymmin,w3d.dy,ny)
      rmesh = sqrt(xmesh**2 + ymesh**2)
      dr = (self.f.xmmax - self.f.xmmin)/self.f.nx
      self.irmesh = int(rmesh/dr)
      self.wrmesh =     rmesh/dr  - self.irmesh
      self.wrmesh = where(self.irmesh >= self.f.nx,1,self.wrmesh)
      self.irmesh = where(self.irmesh >= self.f.nx,self.f.nx-1,self.irmesh)

    i1 = self.nx0_r
    i2 = self.nxm_r+1
    j1 = self.ny0_r
    j2 = self.nym_r+1
    for i in range(2):
      phi[i1:i2,j1:j2,iz+i] = (
          take(savedphi[:,0,i],self.irmesh  )*(1.-self.wrmesh) +
          take(savedphi[:,0,i],self.irmesh+1)*self.wrmesh)

