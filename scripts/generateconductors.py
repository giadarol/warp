"""
This module contains classes for generating the conductor data from a
combination of simple geometrical elements.
The following elements are defined:

Plane(z0=0.,zsign=1,theta=0.,phi=0.,...)
Box(xsize,ysize,zsize,...)

Cylinders:
 Cylinder(radius,length,theta=0.,phi=0.,...)
 ZCylinder(radius,length,...)
 ZCylinderOut(radius,length,...)
 ZCylinderElliptic(ellipticity,radius,length,...)
 ZCylinderEllipticOut(ellipticity,radius,length,...)
 ZRoundedCylinder(radius,length,radius2,...)
 ZRoundedCylinderOut(radius,length,radius2,...)
 XCylinder(radius,length,...)
 XCylinderOut(radius,length,...)
 XCylinderElliptic(ellipticity,radius,length,...)
 XCylinderEllipticOut(ellipticity,radius,length,...)
 YCylinder(radius,length,...)
 YCylinderOut(radius,length,...)
 YCylinderElliptic(ellipticity,radius,length,...)
 YCylinderEllipticOut(ellipticity,radius,length,...)
 ZAnnulus(rmin,rmax,length,...)
 ZAnnulusElliptic(ellipticity,rmin,rmax,length,...)

Cones:
 Cone(r_zmin,r_zmax,length,theta=0.,phi=0.,...)
 ZCone(r_zmin,r_zmax,length,...)
 ZConeOut(r_zmin,r_zmax,length,...)
 ConeSlope(slope,intercept,length,theta=0.,phi=0.,...)
 ZConeSlope(slope,intercept,length,...)
 ZConeOutSlope(slope,intercept,length,...)

Other:
 Sphere(radius,...)
 ZElliptoid(ellipticity,radius,...)
 ZTorus(r1,r2,...)
 Beamletplate(za,zb,z0,thickness,...)

Surfaces of revolution:
 ZSrfrvOut(rofzfunc,zmin,zmax,rmax,...)
 ZSrfrvIn(rofzfunc,zmin,zmax,rmin,...)
 ZSrfrvInOut(rminofz,rmaxofz,zmin,zmax,...)
 ZSrfrvEllipticOut(ellipticity,rofzfunc,zmin,zmax,rmax,...)
 ZSrfrvEllipticIn(ellipticity,rofzfunc,zmin,zmax,rmin,...)
 ZSrfrvEllipticInOut(ellipticity,rminofz,rmaxofz,zmin,zmax,...)
 XSrfrvOut(rofzfunc,zmin,zmax,rmax,...)
 XSrfrvIn(rofzfunc,zmin,zmax,rmin,...)
 XSrfrvInOut(rminofz,rmaxofz,zmin,zmax,...)
 YSrfrvOut(rofzfunc,zmin,zmax,rmax,...)
 YSrfrvIn(rofzfunc,zmin,zmax,rmin,...)
 YSrfrvInOut(rminofz,rmaxofz,zmin,zmax,...)

Note that all take the following additional arguments:
voltage=0.,xcent=0.,ycent=0.,zcent=0.,condid=1,
name=None,material='SS',laccuimagecharge=0,neumann=0

installconductors(a,...): generates the data needed for the fieldsolve
                          See its documentation for the additional arguments.

All of the conductor objects have the following methods:
  distance(xx,yy,zz): calculates the shortest distance between each of the
                      points and the conductor. It returns an instance
                      of a Distance object whichs contains the distance in
                      an attribute named distance.
  isinside(xx,yy,zz): determines whether or not each of the points is inside
                      of the conductor. It returns an instance
                      of an Isinside object whichs contains the flag in
                      an attribute named isinside.
  intercept(xx,yy,zz,vx,vy,vz): calculates the location were a particle with
                                the given velocities most recently intersected
                                a conductor. It returns an instance of an
                                Intercept object whichs contains the data in
                                the attributes xi, yi, zi, and angles itheta
                                and iphi of the surface normal.

A set of classes for generating and manipulating surfaces of revolution
using lines and arcs primitives (SRFRVLA) are available:
SRFRVLAfromfile(filename,voltages,condids,zshifts=None,rshifts=None,install=1)
SRFRVLAsystem(SRFRVLAconds,install=1)
SRFRVLAcond(name,parts,voltage,condid,install=1)
SRFRVLApart(name,data)
SRFRVLA_circle(name,c,r)
SRFRVLA_rectangle(name,c,l,h)
SRFRVLA_rnd_rectangle(name,c,l,h,r)
"""

# The following classes are also defined but should not be directly used.
# Grid
# Assembly
# AssemblyNot
# AssemblyAnd
# AssemblyPlus
# AssemblyMinus
# Delta

from warp import *
import operator
import Opyndx
from string import *
from appendablearray import *

generateconductorsversion = "$Id: generateconductors.py,v 1.183 2008/04/03 18:25:34 dave Exp $"
def generateconductors_doc():
  import generateconductors
  print generateconductors.__doc__

##############################################################################
installedconductors = []
def installconductors(a,xmin=None,xmax=None,ymin=None,ymax=None,
                        zmin=None,zmax=None,dfill=2.,
                        zbeam=None,
                        nx=None,ny=None,nzlocal=None,nz=None,
                        xmmin=None,xmmax=None,ymmin=None,ymmax=None,
                        zmmin=None,zmmax=None,zscale=1.,l2symtry=None,l4symtry=None,
                        installrz=None,gridmode=1,solvergeom=None,
                        conductors=None,gridrz=None,mgmaxlevels=None,
                        my_index=None,nslaves=None,izfsslave=None,nzfsslave=None):
  """
Installs the given conductors.
  - a: the assembly of conductors
  - xmin,xmax,ymin,ymax,zmin,zmax: extent of conductors. Defaults to the
    mesh size. These can be set for optimization, to avoid looking
    for conductors where there are none. Also, they can be used crop a
    conductor
  - dfill=2.: points at a depth in the conductor greater than dfill
              are skipped.
  - zbeam=top.zbeam: location of the beam frame
  - nx,ny,nz: Number of grid cells in the mesh. Defaults to values from w3d
  - xmmin,xmmax,ymmin,ymmax,zmmin,zmmax: extent of mesh. Defaults to values
                                         from w3d
  - zscale=1.: scale factor on dz. This is used when the relativistic scaling
              is done for the longitudinal dimension
  - l2symtry,l4symtry: assumed transverse symmetries. Defaults to values
                       from w3d
  """
  if conductors is None and gridrz is None:
    # --- If conductors was not specified, first check if mesh-refinement
    # --- or other special solver is being used.
    solver = getregisteredsolver()
    import __main__
    if solver is not None:
      solver.installconductor(a,dfill=dfill)
      return
    elif __main__.__dict__.has_key("AMRtree"):
      __main__.__dict__["AMRtree"].installconductor(a,dfill=dfill)
      return

  # --- Use whatever conductors object was specified, or
  # --- if no special solver is being used, use f3d.conductors.
  if conductors is None: conductors = f3d.conductors

  # --- Set the installrz argument if needed.
  if installrz is None:
    installrz = (frz.getpyobject('basegrid') is not None)

  # First, create a grid object
  g = Grid(xmin,xmax,ymin,ymax,zmin,zmax,zbeam,nx,ny,nzlocal,nz,
           xmmin,xmmax,ymmin,ymmax,zmmin,zmmax,zscale,l2symtry,l4symtry,
           installrz,gridrz,
           mgmaxlevels=mgmaxlevels,
           my_index=my_index,nslaves=nslaves,
           izslave=izfsslave,nzslave=nzfsslave)
  # Generate the conductor data
  g.getdata(a,dfill)
  # Then install it
  g.installdata(installrz,gridmode,solvergeom,conductors,gridrz)
  installedconductors.append(a)

def uninstallconductors(a):
  "Removes the conductors from the list of installed conductors"
  installedconductors.remove(a)

##############################################################################
##############################################################################
##############################################################################
listofallconductors = []
class Assembly(Opyndx.Visualizable):
  """
Class to hold assemblies of conductors.  Base class of all conductors.
Should never be directly created by the user.
 - v=0.: voltage on conductor
 - x,y,z=0.,0.,0: center of conductor
 - condid=1: conductor identification number, can be 'next' in which case
             a unique ID is chosen
 - kwlist=[]: list of string names of variable describing conductor
 - generatorf=None: function which generates the distances between the points
                    and the conductors along the axis.
 - generatord=None: function which generates the smallest distance between the
                    points and the conductor surface.
 - name=None: conductor name (string)
 - material='SS': conductor material
 - laccuimagecharge=0: Flags accumulation of image charges
  """

  voltage = 0.
  xcent = 0.
  ycent = 0.
  zcent = 0.
  nextcondid = 1

  __inputs__ = {'name':None,'material':'SS','laccuimagecharge':0,'neumann':0}

  def __init__(self,v=0.,x=0.,y=0.,z=0.,condid=1,kwlist=[],
                    generatorf=None,generatord=None,generatori=None,kw={}):
    listofallconductors.append(self)
    self.voltage = v
    self.xcent = x
    self.ycent = y
    self.zcent = z
    if condid == 'next':
      self.condid = Assembly.nextcondid
      Assembly.nextcondid += 1
    else:
      self.condid = condid
    self.kwlist = kwlist
    self.generatorf = generatorf
    self.generatord = generatord
    self.generatori = generatori

    while 'kw' in kw:
      kwtemp = kw['kw']
      del kw['kw']
      kw.update(kwtemp)

    for name,default in Assembly.__inputs__.items():
      self.__dict__[name] = kw.get(name,default)
      if name in kw: del kw[name]
    assert len(kw) == 0,"Invalid keyword arguments "+str(kw.keys())

    self.lostparticles_data  = AppendableArray(typecode='d',unitshape=[4])
    self.emitparticles_data  = AppendableArray(typecode='d',unitshape=[4])
    self.imageparticles_data = AppendableArray(typecode='d',unitshape=[2])
    self.lostparticles_angles    = {}
    self.lostparticles_energies  = {}
    self.lostparticles_minenergy = {}
    self.lostparticles_maxenergy = {}
    self.emitparticles_angles    = {}
    self.emitparticles_energies  = {}
    self.emitparticles_minenergy = {}
    self.emitparticles_maxenergy = {}

    self.accuimagechargeenabled = 0
    if self.laccuimagecharge: self.enable_accuimagecharge()

  def getkwlist(self):
    kwlist = []
    for k in self.kwlist:
      kwlist.append(self.__dict__[k])
    kwlist.append(self.__dict__['xcent'])
    kwlist.append(self.__dict__['ycent'])
    kwlist.append(self.__dict__['zcent'])
    return kwlist

  def getextent(self):
    return self.extent

  def createextent(self,mins,maxs):
    self.extent = ConductorExtent(
               [self.xcent+mins[0],self.ycent+mins[1],self.zcent+mins[2]],
               [self.xcent+maxs[0],self.ycent+maxs[1],self.zcent+maxs[2]])

  def griddistance(self,ix,iy,iz,xx,yy,zz):
    result = Delta(ix,iy,iz,xx,yy,zz,voltage=self.voltage,condid=self.condid,
                   generator=self.generatorf,neumann=self.neumann,
                   kwlist=self.getkwlist())
    return result

  def distance(self,xx,yy,zz):
    result = Distance(xx,yy,zz,generator=self.generatord,
                      kwlist=self.getkwlist())
    return result

  def isinside(self,xx,yy,zz,aura=0.):
    result = IsInside(xx,yy,zz,generator=self.generatord,
                      condid=self.condid,aura=aura,kwlist=self.getkwlist())
    return result

  def intercept(self,xx,yy,zz,vx,vy,vz):
    result = Intercept(xx,yy,zz,vx,vy,vz,generator=self.generatori,
                       condid=self.condid,conductor=self,
                       kwlist=self.getkwlist())
    return result

  def draw(self,**kw):
    nowarn = kw.get('nowarn',0)
    if not nowarn:
      print 'draw method not implemented for '+self.__class__.__name__

  def drawzx(self,**kw):
    self.draw(**kw)

  def drawzr(self,**kw):
    kw.setdefault('fullplane',0)
    self.draw(**kw)

  def plotdata(self,r,z,color='fg',filled=None,fullplane=1):
    z = array(z) + self.zcent
    r = array(r)
    if color is not None:
      plg(self.xcent+r,z,color=color)
      if fullplane:
        plg(self.xcent-array(r),z,color=color)
    if filled is not None:
      if filled == 'condid': filled = self.condid
      c = array([filled]).astype(ubyte)
      plfp(c,self.xcent+r,z,[len(r)])
      if fullplane:
        plfp(c,self.xcent-array(r),z,[len(r)])

  def get_current_history(self,js=None,l_lost=1,l_emit=1,l_image=1,tmin=None,tmax=None,nt=100):
    """
  Returns conductor current history:
    - js=None    : select species to consider (default None means that contribution from all species are added)
    - l_lost=1   : logical flag setting wether lost particles are taken into account
    - l_emit=1   : logical flag setting wether emitted particles are taken into account
    - l_image=1  : logical flag setting wether image "particles" are taken into account
    - t_min=None : min time
    - t_max=None : max time
    - nt=100     : nb of cells
    """
    tminl=0.
    tmaxl=top.time
    tmine=0.
    tmaxe=top.time
    tmini=0.
    tmaxi=top.time
    # collect lost particles data
    nl = 0
    datal = self.lostparticles_data.data()
    if l_lost and len(datal) > 0:
      ql = datal[:,1].copy()
      tl = datal[:,0].copy()
      if tmin is None:tminl=min(tl)
      if tmax is None:tmaxl=max(tl)
      if js is not None:
        jl = datal[:,3].copy()
        ql = compress(jl==js,ql)
        tl = compress(jl==js,tl)
      nl = shape(ql)[0]
    # collect emitted particles data
    ne = 0
    datae = self.emitparticles_data.data()
    if l_emit and len(datae) > 0:
      qe = -datae[:,1].copy()
      te =  datae[:,0].copy()
      if tmin is None:tmine=min(te)
      if tmax is None:tmaxe=max(te)
      if js is not None:
        je = datae[:,3].copy()
        qe = compress(je==js,qe)
        te = compress(je==js,te)
      ne = shape(qe)[0]
    # collect accumulated image data
    ni = 0
    datai = self.imageparticles_data.data()
    if l_image and len(datai) > 2:
      qi = 0.5*(datai[2:,1]-datai[:-2,1])
      ti = datai[1:-1,0].copy()
      if tmin is None:tmini=min(ti)
      if tmax is None:tmaxi=max(ti)
      ni = shape(qi)[0]
    # setup time min/max and arrays
    if tmin is None:tmin=min(tminl,tmine,tmini)
    if tmax is None:tmax=max(tmaxl,tmaxe,tmaxi)
    qt = zeros(nt+1,Float)
    qtmp = zeros(nt+1,Float)
    dt = (tmax-tmin)/nt
    # accumulate data
    if nl>0:deposgrid1d(1,nl,tl,ql,nt,qt,qtmp,tmin,tmax)
    if ne>0:deposgrid1d(1,ne,te,qe,nt,qt,qtmp,tmin,tmax)
    if ni>0:deposgrid1d(1,ni,ti,qi,nt,qt,qtmp,tmin,tmax)
    return arange(tmin,tmax+0.1*dt,dt),qt/dt

  def plot_current_history(self,js=None,l_lost=1,l_emit=1,l_image=1,tmin=None,tmax=None,nt=100,color=black,width=1,type='solid'):
    """
  Plots conductor current history:
    - js=None      : select species to consider (default None means that contribution from all species are added)
    - l_lost=1     : logical flag setting wether lost particles are taken into account
    - l_emit=1     : logical flag setting wether emitted particles are taken into account
    - l_image=1    : logical flag setting wether image "particles" are taken into account
    - t_min=None   : min time
    - t_max=None   : max time
    - nt=100       : nb of cells
    - color=black  : line color
    - widht=1      : line width
    - type='solid' : line type
    """
    if me<>0:return
    time,current=self.get_current_history(js=js,l_lost=l_lost,l_emit=l_emit,l_image=l_image,tmin=tmin,tmax=tmax,nt=nt)
    plg(current,time,color=color,width=width,type=type)
    ptitles('Current history at '+self.name,'time (s)','I (A)')

  def enable_accuimagecharge(self):
    if self.accuimagechargeenabled: return
    self.accuimagechargeenabled = 1
    if not isinstalledafterfs(self.accuimagecharge):
      installafterfs(self.accuimagecharge)

  def disable_accuimagecharge(self):
    if not self.accuimagechargeenabled: return
    self.accuimagechargeenabled = 0
    if isinstalledafterfs(self.accuimagecharge):
      uninstallafterfs(self.accuimagecharge)

  def accuimagecharge(self,doplot=false,l_verbose=false):
    """
    """
    # get extents
    mins = self.getextent().mins
    maxs = self.getextent().maxs
    g = getregisteredsolver()
    if g is None and (w3d.solvergeom in [w3d.RZgeom,w3d.XYgeom,w3d.XZgeom]):
      g = frz.basegrid
      # --- Note that frz.basegrid.phi has fortran ordering and the .shape
      # --- attribute can only be changed on C ordered arrays. The transpose
      # --- converts from fortran to C ordering so shape can be applied.
      # --- The modified array is then tranposed back.
      phit = transpose(frz.basegrid.phi)
      phit.shape = [phit.shape[0],1,phit.shape[1]]
      phi = transpose(phit)[1:-1,:,:]
      rhot = transpose(frz.basegrid.rho)
      rhot.shape = [rhot.shape[0],1,rhot.shape[1]]
      rho = transpose(rhot)
      nx = frz.basegrid.nr
      ny = 0
      nz = frz.basegrid.nz
      if lparallel: nzlocal = frz.basegrid.nzpar
      else:         nzlocal = nz
      dx = frz.basegrid.dr
      dy = 1.
      dz = frz.basegrid.dz
      xmmin = frz.basegrid.rmin
      xmmax = frz.basegrid.rmax
      ymmin = mins[1]
      ymmax = mins[1]
      zmmin = frz.basegrid.zmin
      zmmax = frz.basegrid.zmax
      izfsslave = top.izfsslave
      # --- This needs to be fixed to handle mesh refinement properly
      l2symtry = w3d.l2symtry
      l4symtry = w3d.l4symtry
    else:
      if g is None:
        g = w3d
        interior = f3d.conductors.interior
      else:
        try:
          interior = g.getconductorobject().interior
        except AttributeError:
          interior = g.conductors.interior
      phi = g.phi[1:-1,1:-1,1:-1]
      rho = g.rho
      nx = g.nx
      ny = g.ny
      nz = g.nz
      nzlocal = g.nzlocal
      dx = g.dx
      dy = g.dy
      dz = g.dz
      xmmin = g.xmmin
      xmmax = g.xmmax
      ymmin = g.ymmin
      ymmax = g.ymmax
      zmmin = g.zmmin
      zmmax = g.zmmax
      izfsslave = g.izfsslave
      l2symtry = g.l2symtry
      l4symtry = g.l4symtry

    # compute mins and maxs
    xmin = max(xmmin,mins[0])
    ymin = max(ymmin,mins[1])
    zmin = max(zmmin,mins[2])
    xmax = min(xmmax,maxs[0])
    ymax = min(ymmax,maxs[1])
    zmax = min(zmmax,maxs[2])

    # get box boundaries at nodes locations
    ixmin = max(0,  int((xmin-xmmin)/dx))
    iymin = max(0,  int((ymin-ymmin)/dy))
    izmin = max(0,  int((zmin-zmmin)/dz))
    izminlocal = max(0,  int((zmin-zmmin)/dz) - izfsslave[me])
    ixmax = min(nx, int((xmax-xmmin)/dx)+1)
    iymax = min(ny, int((ymax-ymmin)/dy)+1)
    izmax = min(nz, int((zmax-zmmin)/dz)+1)
    izmaxlocal = min(nzlocal, int((zmax-zmmin)/dz)+1 - izfsslave[me])

    # --- When in parallel, there is some overlap of the rho and phi arrays.
    # --- This statement increases izminlocal so that the overlap region
    # --- does not get multiply counted.
    if me>0:izminlocal=max(izminlocal,2)

    # --- The "+1" is needed because of the guard cells in z for phi.
    izminp = izminlocal+1
    izmaxp = izmaxlocal+1

    if w3d.solvergeom in [w3d.RZgeom]:
      # --- The rfac is to take into account the r dtheta term in the integrals.
      rfac = 2.*pi*(iota(ixmin,ixmax)*dx + xmmin)
      rfac.shape = (rfac.shape[0],1)
      rmin = 2.*pi*((ixmin - 0.5)*dx + xmmin)
      rmax = 2.*pi*((ixmax + 0.5)*dx + xmmin)
    else:
      rfac = ones([ixmax-ixmin+1,iymax-iymin+1],'d')
      rmin = 1.
      rmax = 1.

    # accumulate charge due to integral form of Gauss Law
    q = 0.
    qc = 0.
    if ixmax >= ixmin and iymax >= iymin and izmaxp >= izminp:

      # compute total charge inside volume
      qc = sum(sum(sum(rho[ixmin:ixmax+1,iymin:iymax+1,izminlocal:izmaxlocal+1]*
                       rfac[:,:,NewAxis])))*dx*dy*dz

# --- This block of code is needed if the rho in conductor interiors is
# --- not zeroed out. Note that the setting of interior above is also needed.
      # --- subtract off any charge inside of the conductors (This charge is
      # --- included in the sum over rho, but it does not affect the
      # --- potential and so is not represented in sum of Enormal and so is not
      # --- accounted for properly. It needs to be explicitly subtracted off
      # --- since it should not be included as image charge.)
      qinterior=zeros(1,'d')
      if (w3d.solvergeom in [w3d.RZgeom,w3d.XYgeom,w3d.XZgeom] and
          getregisteredsolver() is None):
        cond_sumrhointerior2d(qinterior,g,nx,nzlocal,rho[:,0,:],
                              ixmin,ixmax,izminlocal,izmaxlocal,dx,xmmin)
      else:
        subcond_sumrhointerior(qinterior,interior,nx,ny,nzlocal,rho,
                               ixmin,ixmax,iymin,iymax,izminlocal,izmaxlocal)
      qc = qc - qinterior[0]*dx*dy*dz

      # --- Sum the normal E field on the surface of the volume
      # --- When in parallel, z planes are skipped if they are not at the edge of
      # --- the region of interest.
      if 0<=izmaxlocal<nzlocal and izmax == izmaxlocal+izfsslave[me]:
        q += sum(sum(( phi[ixmin:ixmax+1, iymin:iymax+1, izmaxp+1       ]
                      -phi[ixmin:ixmax+1, iymin:iymax+1, izmaxp         ])*rfac))*dx*dy/dz
      if 0<izminlocal<=nzlocal and izmin == izminlocal+izfsslave[me]:
        q += sum(sum(( phi[ixmin:ixmax+1, iymin:iymax+1, izminp-1       ]
                      -phi[ixmin:ixmax+1, iymin:iymax+1, izminp         ])*rfac))*dx*dy/dz
      if 0<=ixmax<nx:
        q += sum(sum(( phi[ixmax+1,       iymin:iymax+1, izminp:izmaxp+1]
                      -phi[ixmax,         iymin:iymax+1, izminp:izmaxp+1])*rmax))*dz*dy/dx
      if 0<ixmin<=nx:
        q += sum(sum(( phi[ixmin-1,       iymin:iymax+1, izminp:izmaxp+1]
                      -phi[ixmin,         iymin:iymax+1, izminp:izmaxp+1])*rmin))*dz*dy/dx
      if 0<=iymax<ny:
        q += sum(sum(( phi[ixmin:ixmax+1, iymax+1,       izminp:izmaxp+1]
                      -phi[ixmin:ixmax+1, iymax,         izminp:izmaxp+1])))*dx*dz/dy
      if 0<iymin<=ny:
        q += sum(sum(( phi[ixmin:ixmax+1, iymin-1,       izminp:izmaxp+1]
                      -phi[ixmin:ixmax+1, iymin,         izminp:izmaxp+1])))*dx*dz/dy

      # correct for symmetries (this will never be done for RZ so rfac is not needed)
      if l4symtry:
       q=q*4.
       qc=qc*4.
      elif l2symtry:
       q=q*2.
       qc=qc*2.
      if l2symtry or l4symtry:
       if iymin==0:
         if 0<=izmaxlocal< nzlocal and izmax == izmaxlocal+izfsslave[me]: q -= 2.*sum(phi[ixmin:ixmax+1,iymin,izmaxp+1] -phi[ixmin:ixmax+1,iymin,izmaxp] )*dx*dy/dz
         if 0< izminlocal<=nzlocal and izmin == izminlocal+izfsslave[me]: q -= 2.*sum(phi[ixmin:ixmax+1,iymin,izminp-1] -phi[ixmin:ixmax+1,iymin,izminp] )*dx*dy/dz
         if 0<=ixmax< nx: q -= 2.*sum(phi[ixmax+1,iymin,izminp:izmaxp+1]-phi[ixmax,iymin,izminp:izmaxp+1])*dz*dy/dx
         if 0< ixmin<=nx: q -= 2.*sum(phi[ixmin-1,iymin,izminp:izmaxp+1]-phi[ixmin,iymin,izminp:izmaxp+1])*dz*dy/dx
         qc -= 2.*sum(sum(rho[ixmin:ixmax+1,iymin,izminlocal:izmaxlocal+1]))*dx*dy*dz
      if l4symtry:
       if ixmin==0:
         if 0<=izmaxlocal< nzlocal and izmax == izmaxlocal+izfsslave[me]: q -= 2.*sum(phi[ixmin,iymin:iymax+1,izmaxp+1] -phi[ixmin,iymin:iymax+1,izmaxp] )*dx*dy/dz
         if 0< izminlocal<=nzlocal and izmin == izminlocal+izfsslave[me]: q -= 2.*sum(phi[ixmin,iymin:iymax+1,izminp-1] -phi[ixmin,iymin:iymax+1,izminp] )*dx*dy/dz
         if 0<=iymax< ny: q -= 2.*sum(phi[ixmin,iymax+1,izminp:izmaxp+1]-phi[ixmin,iymax,izminp:izmaxp+1])*dx*dz/dy
         if 0< iymin<=ny: q -= 2.*sum(phi[ixmin,iymin-1,izminp:izmaxp+1]-phi[ixmin,iymin,izminp:izmaxp+1])*dx*dz/dy
         qc -= 2.*sum(sum(rho[ixmin,iymin:iymax+1,izminlocal:izmaxlocal+1]))*dx*dy*dz
       if ixmin==0 and iymin==0:
         if 0<=izmaxlocal< nzlocal: q += (phi[ixmin,iymin,izmaxp+1]-phi[ixmin,iymin,izmaxp])*dx*dy/dz
         if 0< izminlocal<=nzlocal: q += (phi[ixmin,iymin,izminp-1]-phi[ixmin,iymin,izminp])*dx*dy/dz
         qc += sum(rho[ixmin,iymin,izminlocal:izmaxlocal+1])*dx*dy*dz

    # --- Gather up the charges from all of the parallel processors.
    q  = parallelsum(q)
    qc = parallelsum(qc)

    self.imageparticles_data.append(array([top.time,q*eps0+qc]))
    if l_verbose:print self.name,q*eps0,qc
    if doplot:
      window(0)
      pldj([zmin,zmin,zmin,zmax],[xmin,xmin,xmax,xmin],[zmax,zmin,zmax,zmax],[xmin,xmax,xmax,xmax],color=red,width=3)
      window(1)
      pldj([zmin,zmin,zmin,zmax],[ymin,ymin,ymax,ymin],[zmax,zmin,zmax,zmax],[ymin,ymax,ymax,ymax],color=red,width=3)
      window(0)
      zmin=zmmin+izminlocal*dz
      zmax=zmmin+izmaxlocal*dz
      xmin=xmmin+ixmin*dx
      xmax=xmmin+ixmax*dx
      ymin=ymmin+iymin*dy
      ymax=ymmin+iymax*dy
      pldj([zmin,zmin,zmin,zmax],[xmin,xmin,xmax,xmin],[zmax,zmin,zmax,zmax],[xmin,xmax,xmax,xmax],color=blue,width=3)
      window(1)
      pldj([zmin,zmin,zmin,zmax],[ymin,ymin,ymax,ymin],[zmax,zmin,zmax,zmax],[ymin,ymax,ymax,ymax],color=blue,width=3)
      window(0)

   # Operations which return an Assembly expression.
  def __mul__(self,right):
    if right is None: return self
    return AssemblyAnd(self,right)
  def __add__(self,right):
    if right is None: return self
    return AssemblyPlus(self,right)
  def __sub__(self,right):
    if right is None: return self
    return AssemblyMinus(self,right)
  def __neg__(self):
    return AssemblyNot(self)
  def __pos__(self):
    return self

  def __rmul__(self,left):
    return self*left
  def __radd__(self,left):
    return self+left
  def __rsub__(self,left):
    return (-self)+left


class AssemblyNot(Assembly):
  """
AssemblyNot class.  Represents 'not' of assemblies.
  """
  def __init__(self,l):
    Assembly.__init__(self,0.,l.xcent,l.ycent,l.zcent,l.condid,
                      kw={'material':l.material,'name':l.name,'neumann':l.neumann})
    self.left = l
  def getextent(self):
    return (-self.left.getextent())
  def griddistance(self,ix,iy,iz,xx,yy,zz):
    return (-(self.left.griddistance(ix,iy,iz,xx,yy,zz)))
  def distance(self,xx,yy,zz):
    return (-(self.left.distance(xx,yy,zz)))
  def isinside(self,xx,yy,zz,aura=0.):
    return (-(self.left.isinside(xx,yy,zz,aura)))
  def intercept(self,xx,yy,zz,vx,vy,vz):
    return (-(self.left.intercept(xx,yy,zz,vx,vy,vz)))
  def createdxobject(self,kwdict={},**kw):
    kw.update(kwdict)
    self.dxobject = self.left.getdxobject(kwdict=kw)
  def draw(self,**kw):
    self.left.draw(**kw)


class AssemblyAnd(Assembly):
  """
AssemblyAnd class.  Represents 'and' of assemblies.
  """
  def __init__(self,l,r):
    Assembly.__init__(self,0.,l.xcent,l.ycent,l.zcent,l.condid,
                      kw={'material':l.material,'name':(l.name or r.name),'neumann':l.neumann})
    self.left = l
    self.right = r
  def getextent(self):
    return (self.left.getextent()*self.right.getextent())
  def griddistance(self,ix,iy,iz,xx,yy,zz):
    return (self.left.griddistance(ix,iy,iz,xx,yy,zz) *
            self.right.griddistance(ix,iy,iz,xx,yy,zz))
  def distance(self,xx,yy,zz):
    return (self.left.distance(xx,yy,zz) *
            self.right.distance(xx,yy,zz))
  def isinside(self,xx,yy,zz,aura=0.):
    return (self.left.isinside(xx,yy,zz,aura) *
            self.right.isinside(xx,yy,zz,aura))
  def intercept(self,xx,yy,zz,vx,vy,vz):
    return (self.left.intercept(xx,yy,zz,vx,vy,vz) *
            self.right.intercept(xx,yy,zz,vx,vy,vz))
  def createdxobject(self,kwdict={},**kw):
    kw.update(kwdict)
    l = self.left.getdxobject(kwdict=kw)
    r = self.right.getdxobject(kwdict=kw)
    self.dxobject = Opyndx.DXCollection(l,r)
  def draw(self,**kw):
    self.left.draw(**kw)
    self.right.draw(**kw)


class AssemblyPlus(Assembly):
  """
AssemblyPlus class.  Represents 'or' of assemblies.
  """
  def __init__(self,l,r):
    Assembly.__init__(self,0.,l.xcent,l.ycent,l.zcent,l.condid,
                      kw={'material':l.material,'name':(l.name or r.name),'neumann':l.neumann})
    self.left = l
    self.right = r
  def getextent(self):
    return (self.left.getextent()+self.right.getextent())
  def griddistance(self,ix,iy,iz,xx,yy,zz):
    return (self.left.griddistance(ix,iy,iz,xx,yy,zz) +
            self.right.griddistance(ix,iy,iz,xx,yy,zz))
  def distance(self,xx,yy,zz):
    return (self.left.distance(xx,yy,zz) +
            self.right.distance(xx,yy,zz))
  def isinside(self,xx,yy,zz,aura=0.):
    return (self.left.isinside(xx,yy,zz,aura) +
            self.right.isinside(xx,yy,zz,aura))
  def intercept(self,xx,yy,zz,vx,vy,vz):
    return (self.left.intercept(xx,yy,zz,vx,vy,vz) +
            self.right.intercept(xx,yy,zz,vx,vy,vz))
  def createdxobject(self,kwdict={},**kw):
    kw.update(kwdict)
    l = self.left.getdxobject(kwdict=kw)
    r = self.right.getdxobject(kwdict=kw)
    self.dxobject = Opyndx.DXCollection(l,r)
  def draw(self,**kw):
    self.left.draw(**kw)
    self.right.draw(**kw)


class AssemblyMinus(Assembly):
  """
AssemblyMinus class.
  """
  def __init__(self,l,r):
    Assembly.__init__(self,0.,l.xcent,l.ycent,l.zcent,l.condid,
                      kw={'material':l.material,'name':(l.name or r.name),'neumann':l.neumann})
    self.left = l
    self.right = r
  def getextent(self):
    return (self.left.getextent()-self.right.getextent())
  def griddistance(self,ix,iy,iz,xx,yy,zz):
    return (self.left.griddistance(ix,iy,iz,xx,yy,zz) -
            self.right.griddistance(ix,iy,iz,xx,yy,zz))
  def distance(self,xx,yy,zz):
    return (self.left.distance(xx,yy,zz) -
            self.right.distance(xx,yy,zz))
  def isinside(self,xx,yy,zz,aura=0.):
    return (self.left.isinside(xx,yy,zz,aura) -
            self.right.isinside(xx,yy,zz,aura))
  def intercept(self,xx,yy,zz,vx,vy,vz):
    return (self.left.intercept(xx,yy,zz,vx,vy,vz) -
            self.right.intercept(xx,yy,zz,vx,vy,vz))
  def createdxobject(self,kwdict={},**kw):
    kw.update(kwdict)
    l = self.left.getdxobject(kwdict=kw)
    r = self.right.getdxobject(kwdict=kw)
    self.dxobject = Opyndx.DXCollection(l,r)
  def draw(self,**kw):
    self.left.draw(**kw)
    self.right.draw(**kw)

#============================================================================
class EllipticAssembly(Assembly):
  """
Elliptic assembly
  """
  def __init__(self,ellipticity,v=0.,x=0.,y=0.,z=0.,condid=1,kwlist=[],
                    generatorf=None,generatord=None,generatori=None,**kw):
    Assembly.__init__(self,v,x,y,z,condid,kwlist,
                           self.ellipseconductorf,self.ellipseconductord,
                           self.ellipseintercept,
                           kw=kw)
    self.ellipticity = ellipticity
    self.circlegeneratorf = generatorf
    self.circlegeneratord = generatord
    self.circlegeneratori = generatori
    self.extent.toellipse(self.ellipticity)

  def ellipseconductorf(self,*argtuple):
    arglist = list(argtuple)
    if self.ellipticity != 1.:
      y = arglist[-9]
      delmy = arglist[-5]
      delpy = arglist[-4]
      arglist[-9] = y/self.ellipticity

    apply(self.circlegeneratorf,arglist)

    if self.ellipticity != 1.:
      delmy[:] = delmy*self.ellipticity
      delpy[:] = delpy*self.ellipticity

  def ellipseconductord(self,*argtuple):
    arglist = list(argtuple)
    if self.ellipticity != 1.:
      x = arglist[-4]
      y = arglist[-3]
      distance = arglist[-1]
      arglist[-3] = y/self.ellipticity

    apply(self.circlegeneratord,arglist)

    if self.ellipticity != 1.:
      tt = arctan2(y,x)
      dx = distance*cos(tt)
      dy = distance*sin(tt)*self.ellipticity
      distance[:] = sqrt(dx**2 + dy**2)*sign(distance)

  def ellipseintercept(self,*argtuple):
    arglist = list(argtuple)
    if self.ellipticity != 1.:
      y = arglist[-10]
      xi = arglist[-5]
      yi = arglist[-4]
      iphi = arglist[-1]
      arglist[-10] = y/self.ellipticity

    apply(self.circlegeneratori,arglist)

    if self.ellipticity != 1.:
      yi[:] = yi*self.ellipticity
      iphi[:] = arctan2(sin(iphi),self.ellipticity*cos(iphi))

  def __getstate__(self):
    """
An explicit getstate is needed in order for this to be picklable. The
generator attributes are functions and cannot be pickled so they must be
removed from the dictionary when pickling.
    """
    dict = self.__dict__.copy()
    del dict['generatorf']
    del dict['generatord']
    del dict['generatori']
    return dict

  def __setstate__(self,dict):
    """
This explicit setstate restores the generator attributes.
    """
    self.__dict__.update(dict)
    self.generatorf = self.ellipseconductorf
    self.generatord = self.ellipseconductord
    self.generatori = self.ellipseintercept

#============================================================================
class XAssembly(Assembly):
  """
Assembly aligned along X axis
  """
  def __init__(self,v=0.,x=0.,y=0.,z=0.,condid=1,kwlist=[],
                    generatorf=None,generatord=None,generatori=None,**kw):
    Assembly.__init__(self,v,x,y,z,condid,kwlist,
                           self.xconductorf,self.xconductord,
                           self.xintercept,
                           kw=kw)
    self.zgeneratorf = generatorf
    self.zgeneratord = generatord
    self.zgeneratori = generatori
    self.extent.toX(self.xcent,self.ycent,self.zcent)

  def xconductorf(self,*argtuple):
    arglist = list(argtuple)
    # --- permutate the object center
    arglist[-14] = argtuple[-13]
    arglist[-13] = argtuple[-12]
    arglist[-12] = argtuple[-14]
    # --- permutate coordinates
    arglist[-10] = argtuple[- 9]
    arglist[- 9] = argtuple[- 8]
    arglist[- 8] = argtuple[-10]
    # --- permutate deltas
    arglist[-7] = argtuple[-5]
    arglist[-6] = argtuple[-4]
    arglist[-5] = argtuple[-3]
    arglist[-4] = argtuple[-2]
    arglist[-3] = argtuple[-7]
    arglist[-2] = argtuple[-6]
    apply(self.zgeneratorf,arglist)

  def xconductord(self,*argtuple):
    arglist = list(argtuple)
    # --- permutate the object center
    arglist[-8] = argtuple[-7]
    arglist[-7] = argtuple[-6]
    arglist[-6] = argtuple[-8]
    # --- permutate coordinates
    arglist[-4] = argtuple[-3]
    arglist[-3] = argtuple[-2]
    arglist[-2] = argtuple[-4]
    apply(self.zgeneratord,arglist)

  def xintercept(self,*argtuple):
    arglist = list(argtuple)
    # --- permutate the object center
    arglist[-15] = argtuple[-14]
    arglist[-14] = argtuple[-13]
    arglist[-13] = argtuple[-15]
    # --- permutate coordinates
    arglist[-11] = argtuple[-10]
    arglist[-10] = argtuple[- 9]
    arglist[- 9] = argtuple[-11]
    # --- permutate velocities
    arglist[-8] = argtuple[-7]
    arglist[-7] = argtuple[-6]
    arglist[-6] = argtuple[-8]
    # --- permutate intercept coordinates
    arglist[-5] = argtuple[-4]
    arglist[-4] = argtuple[-3]
    arglist[-3] = argtuple[-5]
    # --- Create temps for the surface normal angles
    arglist[-2] = zeros(shape(argtuple[-2]),'d')
    arglist[-1] = zeros(shape(argtuple[-1]),'d')
    apply(self.zgeneratori,arglist)
    # --- Undo the surface normals
    ttheta = arglist[-2]
    tphi = arglist[-1]
    if 0:  # original coding from DPG
      itheta = arctan2(sqrt(cos(ttheta)**2 + (cos(tphi)*sin(ttheta))**2),
                       sin(tphi)*sin(ttheta))
      iphi = arctan2(cos(tphi)*sin(ttheta),cos(ttheta))
    else: # alternative from JLV
      y=cos(tphi)*sin(ttheta)
      z=sin(tphi)*sin(ttheta)
      x=cos(ttheta)
      r=sqrt(x*x+y*y)
      itheta = arctan2(r,z)
      iphi   = arctan2(y,x)
    argtuple[-2][:] = itheta
    argtuple[-1][:] = iphi

  def __getstate__(self):
    """
An explicit getstate is needed in order for this to be picklable. The
generator attributes are functions and cannot be pickled so they must be
removed from the dictionary when pickling.
    """
    dict = self.__dict__.copy()
    del dict['generatorf']
    del dict['generatord']
    del dict['generatori']
    return dict

  def __setstate__(self,dict):
    """
This explicit setstate restores the generator attributes.
    """
    self.__dict__.update(dict)
    self.generatorf = self.xconductorf
    self.generatord = self.xconductord
    self.generatori = self.xintercept

#============================================================================
class YAssembly(Assembly):
  """
Assembly aligned along Y axis
  """
  def __init__(self,v=0.,x=0.,y=0.,z=0.,condid=1,kwlist=[],
                    generatorf=None,generatord=None,generatori=None,**kw):
    Assembly.__init__(self,v,x,y,z,condid,kwlist,
                           self.yconductorf,self.yconductord,
                           self.yintercept,
                           kw=kw)
    self.zgeneratorf = generatorf
    self.zgeneratord = generatord
    self.zgeneratori = generatori
    self.extent.toY(self.xcent,self.ycent,self.zcent)

  def yconductorf(self,*argtuple):
    arglist = list(argtuple)
    # --- permutate the object center
    arglist[-14] = argtuple[-12]
    arglist[-13] = argtuple[-14]
    arglist[-12] = argtuple[-13]
    # --- permutate coordinates
    arglist[-10] = argtuple[- 8]
    arglist[- 9] = argtuple[-10]
    arglist[- 8] = argtuple[- 9]
    # --- permutate deltas
    arglist[-7] = argtuple[-3]
    arglist[-6] = argtuple[-2]
    arglist[-5] = argtuple[-7]
    arglist[-4] = argtuple[-6]
    arglist[-3] = argtuple[-5]
    arglist[-2] = argtuple[-4]
    apply(self.zgeneratorf,arglist)

  def yconductord(self,*argtuple):
    arglist = list(argtuple)
    # --- permutate the object center
    arglist[-8] = argtuple[-6]
    arglist[-7] = argtuple[-8]
    arglist[-6] = argtuple[-7]
    # --- permutate coordinates
    arglist[-4] = argtuple[-2]
    arglist[-3] = argtuple[-4]
    arglist[-2] = argtuple[-3]
    apply(self.zgeneratord,arglist)

  def yintercept(self,*argtuple):
#    l=len(argtuple)
#    for i in range(l):
#      print -i-1,argtuple[-i-1]
    arglist = list(argtuple)
    # --- permutate the object center
    arglist[-15] = argtuple[-13]
    arglist[-14] = argtuple[-15]
    arglist[-13] = argtuple[-14]
    # --- permutate coordinates
    arglist[-11] = argtuple[- 9]
    arglist[-10] = argtuple[-11]
    arglist[- 9] = argtuple[-10]
    # --- permutate velocities
    arglist[-8] = argtuple[-6]
    arglist[-7] = argtuple[-8]
    arglist[-6] = argtuple[-7]
    # --- permutate intercept coordinates
    arglist[-5] = argtuple[-3]
    arglist[-4] = argtuple[-5]
    arglist[-3] = argtuple[-4]
    # --- Create temps for the surface normal angles
    arglist[-2] = zeros(shape(argtuple[-2]),'d')
    arglist[-1] = zeros(shape(argtuple[-1]),'d')
    apply(self.zgeneratori,arglist)
    # --- Undo the surface normals
    ttheta = arglist[-2]
    tphi = arglist[-1]
    if 0:  # original coding from DPG
      itheta = arctan2(sqrt((sin(tphi)*sin(ttheta))**2 + cos(ttheta)**2),
                     cos(tphi)*sin(ttheta))
      iphi = arctan2(cos(ttheta),sin(tphi)*sin(ttheta))
    else: # alternative from JLV
      z=cos(tphi)*sin(ttheta)
      x=sin(tphi)*sin(ttheta)
      y=cos(ttheta)
      r=sqrt(x*x+y*y)
      itheta = arctan2(r,z)
      iphi   = arctan2(y,x)
    argtuple[-2][:] = itheta
    argtuple[-1][:] = iphi

  def __getstate__(self):
    """
An explicit getstate is needed in order for this to be picklable. The
generator attributes are functions and cannot be pickled so they must be
removed from the dictionary when pickling.
    """
    dict = self.__dict__.copy()
    del dict['generatorf']
    del dict['generatord']
    del dict['generatori']
    return dict

  def __setstate__(self,dict):
    """
This explicit setstate restores the generator attributes.
    """
    self.__dict__.update(dict)
    self.generatorf = self.yconductorf
    self.generatord = self.yconductord
    self.generatori = self.yintercept

##############################################################################
class ConductorExtent:
  """
Class to hold the extent of a conductor. This is somewhat overkill for a
class, but it does provide a nice way of putting this into one spot.
  """
  def __init__(self,mins,maxs):
    self.mins = copy.copy(mins)
    self.maxs = copy.copy(maxs)
  def toellipse(self,ellipticity):
    self.mins[1] = self.mins[1]*ellipticity
    self.maxs[1] = self.maxs[1]*ellipticity
  def toX(self,xcent,ycent,zcent):
    cent = array([xcent,ycent,zcent])
    minscopy = self.mins - cent
    maxscopy = self.maxs - cent
    self.mins[0] = minscopy[2] + xcent
    self.maxs[0] = maxscopy[2] + xcent
    self.mins[1] = minscopy[0] + ycent
    self.maxs[1] = maxscopy[0] + ycent
    self.mins[2] = minscopy[1] + zcent
    self.maxs[2] = maxscopy[1] + zcent
  def toY(self,xcent,ycent,zcent):
    cent = array([xcent,ycent,zcent])
    minscopy = self.mins - cent
    maxscopy = self.maxs - cent
    self.mins[0] = minscopy[1] + xcent
    self.maxs[0] = maxscopy[1] + xcent
    self.mins[1] = minscopy[2] + ycent
    self.maxs[1] = maxscopy[2] + ycent
    self.mins[2] = minscopy[0] + zcent
    self.maxs[2] = maxscopy[0] + zcent

  def __neg__(self):
    "This one is doesn't help much"
    return ConductorExtent([-largepos,-largepos,-largepos],
                           [+largepos,+largepos,+largepos])
  def __add__(self,right):
    return ConductorExtent(minimum(self.mins,right.mins),
                           maximum(self.maxs,right.maxs))
  def __mul__(self,right):
    return ConductorExtent(maximum(self.mins,right.mins),
                           minimum(self.maxs,right.maxs))
  def __sub__(self,right):
    return ConductorExtent(self.mins,self.maxs)

##############################################################################

class Delta:
  """
Class to hold the set of distances in each of the six directions.
Distances have the sign of the outward normal surface vector, i.e.
distances to outside the surface are positive, inside negative.
  """

  def __init__(self,ix=None,iy=None,iz=None,xx=None,yy=None,zz=None,
                    dels=None,vs=None,ns=None,parity=None,
                    voltage=0.,condid=1,generator=None,neumann=0,
                    kwlist=[]):
    self.datalist = []
    self.neumann = neumann
    if ix is None:
      self.ndata = 0
      self.ix = None
      self.iy = None
      self.iz = None
      self.xx = None
      self.yy = None
      self.zz = None
      self.dels = None
      self.vs = None
      self.ns = None
      self.parity = None
      self.mglevel = None
    elif generator is not None:
      self.ndata = len(ix)
      self.ix = ix
      self.iy = iy
      self.iz = iz
      self.xx = xx
      self.yy = yy
      self.zz = zz
      self.dels = zeros((6,self.ndata),'d')
      fuzz = 1.e-13
      apply(generator,kwlist + [self.ndata,self.xx,self.yy,self.zz,
                                self.dels[0,:],self.dels[1,:],
                                self.dels[2,:],self.dels[3,:],
                                self.dels[4,:],self.dels[5,:]] + [fuzz])
      self.setvoltages(voltage)
      self.setcondids(condid)
      self.setlevels(0)
    else:
      self.ndata = len(ix)
      self.ix = ix
      self.iy = iy
      self.iz = iz
      self.xx = xx
      self.yy = yy
      self.zz = zz
      self.dels = dels
      self.vs = vs
      self.ns = int(ns)
      self.parity = parity
      self.setlevels(0)
    self.fuzzsign = -1

    self.append(self)

  def setvoltages(self,voltage):
    "Routine to set appropriate voltages."
    # --- The voltage can be a number, a class instance, or a function.
    # --- If it is an instance, the class must have a method getvolt
    # --- that takes the time as an argument. If a function, it must
    # --- take one argument, the time.
    if type(voltage) == FunctionType:
      v = voltage(top.time)
    elif type(voltage) == InstanceType:
      v = voltage.getvolt(top.time)
    else:
      v = voltage
    self.vs = v + zeros((6,self.ndata),'d')

  def setcondids(self,condid):
    "Routine to setcondid condids."
    self.ns = int(condid) + zeros((6,self.ndata),'l')

  def setlevels(self,level):
    self.mglevel = level + zeros(self.ndata,'l')

  def normalize(self,dx,dy,dz):
    """
Normalizes the data with respect to the grid cell sizes.
dx,dy,dz: the grid cell sizes
    """
    self.dels[:,:] = self.dels/array([dx,dx,dy,dy,dz,dz])[:,NewAxis]
    if self.neumann:
      # --- For points that are within fuzz of 0 or 1, force them to be 0 or 1.
      # --- For Neumann boundaries, dels=0 is a valid value and this deals
      # --- with roundoff problems since a conductor that is aligned with the
      # --- grid should be producing subgrid points (with dels=0) and not
      # --- interior points, which would happen if dels=-fuzz. (The Neumann
      # --- method ignores all interior points.)
      # --- For points near 1, force them to one since they may be duplicating
      # --- neighboring points with dels=0. Neumann ignores points with
      # --- dels==1.
      fuzz = 1.e-13
      self.fixneumannzeros(0,1,fuzz)
      self.fixneumannzeros(1,0,fuzz)
      self.fixneumannzeros(2,3,fuzz)
      self.fixneumannzeros(3,2,fuzz)
      self.fixneumannzeros(4,5,fuzz)
      self.fixneumannzeros(5,4,fuzz)
      self.dels[:,:] = where((1.-fuzz < self.dels[:,:])&(self.dels[:,:] < 1.),
                             1.,self.dels[:,:])
      # --- This deals with points that straddle a surface, where some
      # --- directions are inside and others outside. This ensures that all
      # --- directions are inside in these cases by changing dels that are
      # --- outside so that they are deep inside. This still doesn't seem to
      # --- be enough to fix the problems with convergence.
      delsmin = minimum.reduce(self.dels)
      delsmax = maximum.reduce(self.dels)
      ii = compress((delsmin < 0.)&(delsmax > 0.),range(self.dels.shape[1]))
      for i in ii:
        self.dels[:,i] = where(self.dels[:,i] > 0.,-2.,self.dels[:,i])
      # --- This deals with special case points. It can sometimes happen that
      # --- when a point is on a surface, one of the dels can be zero, but others
      # --- can be between -1 and 0. This is an attempt to fix those points.
      # --- It forces the dels which are between -1 and 0 to be deep inside the
      # --- conductor, i.e., -2. If this is not done, then the normal Dirichlet
      # --- subgrid algorithm would be applied, leading to serious errors.
      # --- However, this fix doesn't seen to fix the whole problem since the
      # --- code still has problems converging when there are points with dels=0.
      delsmin = minimum.reduce(abs(self.dels))
      for i in range(6):
        ccc = where((-1+fuzz<self.dels[i,:])&(self.dels[i,:]<0.),-2.,self.dels[i,:])
        self.dels[i,:] = where(delsmin==0.,ccc,self.dels[i,:])

  def fixneumannzeros(self,i0,i1,fuzz):
    # --- If a point is between -fuzz and 0, make sure it gets put on the
    # --- conductor surface (i.e. set del = 0). Note that the direction
    # --- of the point must be switched (between plus and minus) since
    # --- the point is being moved from inside the surface to on the
    # --- surface.
    self.dels[i0,:] = where((-fuzz < self.dels[i1,:])&(self.dels[i1,:] < 0.),
                            0.,self.dels[i0,:])
    self.dels[i1,:] = where((-fuzz < self.dels[i1,:])&(self.dels[i1,:] < 0.),
                            -2.,self.dels[i1,:])

  def setparity(self,dfill,fuzzsign):
    """
Set parity. For points inside, this is set to -1. For points near the surface,
this is set to the parity of ix+iy+iz. Otherwise defaults to large integer.
This assumes that the data has already been normalized with respect to the
grid cell sizes.
    """
    # --- Using the inplace add is slightly faster since it doesn't have to
    # --- allocate a new array.
    self.parity = zeros(self.ndata,'l')
    add(self.parity,999,self.parity)
    self.fuzzsign = fuzzsign
    if self.neumann: fuzz0 = 0.
    else:            fuzz0 = +1.e-9
    fuzz1 = 1.e-9
    if self.neumann: dfill = 0.
    # --- A compiled routine is called for optimization
    setconductorparity(self.ndata,self.ix,self.iy,self.iz,
                       self.dels,self.parity,fuzz0,fuzz1,fuzzsign,dfill)

  def clean(self):
    """
Removes the data which is far from any conductors. Assumes that setparity
has already been called.
    """
    if with_numpy:
      ii = nonzero(self.parity < 2)[0]
    else:
      ii = nonzero(self.parity < 2)
    self.ix    = take(self.ix,ii)
    self.iy    = take(self.iy,ii)
    self.iz    = take(self.iz,ii)
    self.xx    = take(self.xx,ii)
    self.yy    = take(self.yy,ii)
    self.zz    = take(self.zz,ii)
    self.dels  = take(self.dels,ii,1)
    self.vs    = take(self.vs,ii,1)
    self.ns    = take(self.ns,ii,1)
    self.parity= take(self.parity,ii)
    self.ndata = len(self.ix)

  def append(self,d):
    self.datalist.append(d)

  def install(self,installrz=1,solvergeom=None,conductors=None,grid=None):
    """
Installs the data into the WARP database
    """
    # --- If no conductors object was passed in, use the default one
    # --- from the f3d package.
    if conductors is None: conductors = f3d.conductors

    conductors.fuzzsign = self.fuzzsign
    if self.neumann: delssign = -1
    else:            delssign = +1

    # --- If the RZ solver is being used and the data is to be installed,
    # --- then clear out an existing conductor data in the database first.
    # --- If this is not done, then data from conductor
    # --- objects will be copied in the RZ database multiple time if
    # --- multiple objects are installed separately.
    # --- At the end of the install, the install_conductor_rz routine
    # --- will copy data back from the RZ database to the database.
    # --- This is slightly inefficient since for each object installed, all
    # --- of the accumulated data will be copied back into the database.
    # --- The way around that is to make the call to install_conductors_rz
    # --- only after all of the objects have been installed.
    if solvergeom is None: solvergeom = w3d.solvergeom
    if(installrz and
       (solvergeom in [w3d.RZgeom,w3d.XZgeom,w3d.XYgeom])):
      conductors.interior.n = 0
      conductors.evensubgrid.n = 0
      conductors.oddsubgrid.n = 0

    # --- Count how much data there is
    ncnew = 0
    nenew = 0
    nonew = 0
    for data in self.datalist:
      if data.parity is None: continue
      ncnew += sum(where(data.parity == -1,1,0))
      nenew += sum(where(data.parity == 0,1,0))
      nonew += sum(where(data.parity == 1,1,0))

    ntot = ncnew + nenew + nonew

    nc = conductors.interior.n
    ne = conductors.evensubgrid.n
    no = conductors.oddsubgrid.n

    if ncnew + nc > conductors.interior.nmax:
      conductors.interior.nmax = ncnew + nc
    if nenew + ne > conductors.evensubgrid.nmax:
      conductors.evensubgrid.nmax = nenew + ne
    if nonew + no > conductors.oddsubgrid.nmax:
      conductors.oddsubgrid.nmax = nonew + no

    conductors.gchange("*")

    conductors.interior.n = nc + ncnew
    conductors.evensubgrid.n = ne + nenew
    conductors.oddsubgrid.n = no + nonew

    # --- Install all of the conductor data into the database.
    for data in self.datalist:
      if data.parity is None: continue

      ncnew = sum(where(data.parity == -1,1,0))
      if ncnew > 0:
        if with_numpy:
          ii = nonzero(data.parity == -1)[0]
        else:
          ii = nonzero(data.parity == -1)
        conductors.interior.indx[0,nc:nc+ncnew] = take(data.ix,ii)
        conductors.interior.indx[1,nc:nc+ncnew] = take(data.iy,ii)
        conductors.interior.indx[2,nc:nc+ncnew] = take(data.iz,ii)
        conductors.interior.volt[nc:nc+ncnew] = take(data.vs[0,:],ii)
        conductors.interior.numb[nc:nc+ncnew] = take(data.ns[0,:],ii)
        conductors.interior.ilevel[nc:nc+ncnew] = take(data.mglevel,ii)
        nc = nc + ncnew

      nenew = sum(where(data.parity == 0,1,0))
      if nenew > 0:
        if with_numpy:
          ii = nonzero(data.parity == 0)[0]
        else:
          ii = nonzero(data.parity == 0)
        conductors.evensubgrid.indx[0,ne:ne+nenew] = take(data.ix,ii)
        conductors.evensubgrid.indx[1,ne:ne+nenew] = take(data.iy,ii)
        conductors.evensubgrid.indx[2,ne:ne+nenew] = take(data.iz,ii)
        conductors.evensubgrid.dels[:,ne:ne+nenew] = take(data.dels,ii,1)*delssign
        conductors.evensubgrid.volt[:,ne:ne+nenew] = take(data.vs,ii,1)
        conductors.evensubgrid.numb[:,ne:ne+nenew] = take(data.ns,ii,1)
        conductors.evensubgrid.ilevel[ne:ne+nenew] = take(data.mglevel,ii)
        ne = ne + nenew


      nonew = sum(where(data.parity == 1,1,0))
      if nonew > 0:
        if with_numpy:
          ii = nonzero(data.parity == 1)[0]
        else:
          ii = nonzero(data.parity == 1)
        conductors.oddsubgrid.indx[0,no:no+nonew] = take(data.ix,ii)
        conductors.oddsubgrid.indx[1,no:no+nonew] = take(data.iy,ii)
        conductors.oddsubgrid.indx[2,no:no+nonew] = take(data.iz,ii)
        conductors.oddsubgrid.dels[:,no:no+nonew] = take(data.dels,ii,1)*delssign
        conductors.oddsubgrid.volt[:,no:no+nonew] = take(data.vs,ii,1)
        conductors.oddsubgrid.numb[:,no:no+nonew] = take(data.ns,ii,1)
        conductors.oddsubgrid.ilevel[no:no+nonew] = take(data.mglevel,ii)
        no = no + nonew

    # --- If the RZ solver is being used, then copy the data into that
    # --- database. This also copies all of the accumulated data back into
    # --- the database to allow for plotting and diagnostics.
    if ntot > 0 and installrz:
      if solvergeom in [w3d.RZgeom,w3d.XZgeom,w3d.XYgeom]:
        if grid is None:
          frz.install_conductors_rz(conductors,frz.basegrid)
        else:
          frz.install_conductors_rz(conductors,grid)

  def __neg__(self):
    "Delta not operator."
    return Delta(self.ix,self.iy,self.iz,self.xx,self.yy,self.zz,
                 -self.dels,self.vs,self.ns,neumann=self.neumann)

  def __mul__(self,right):
    "'and' operator, returns maximum of distances to surfaces."
    assert self.neumann == right.neumann,\
      "Neumann objects cannot be mixed with Dirichlet objects"
    c = less(self.dels,right.dels)
    return Delta(self.ix,self.iy,self.iz,self.xx,self.yy,self.zz,
                 choose(c,(self.dels,right.dels)),
                 choose(c,(self.vs  ,right.vs)),
                 choose(c,(self.ns  ,right.ns)),neumann=right.neumann)

  def __add__(self,right):
    "'or' operator, returns minimum of distances to surfaces."
    assert self.neumann == right.neumann,\
      "Neumann objects cannot be mixed with Dirichlet objects"
    c = greater(self.dels,right.dels)
    return Delta(self.ix,self.iy,self.iz,self.xx,self.yy,self.zz,
                 choose(c,(self.dels,right.dels)),
                 choose(c,(self.vs  ,right.vs)),
                 choose(c,(self.ns  ,right.ns)),neumann=right.neumann)

  def __sub__(self,right):
    "'or' operator, returns minimum of distances to surfaces."
    assert self.neumann == right.neumann,\
      "Neumann objects cannot be mixed with Dirichlet objects"
    rdels = -right.dels
    c = less(self.dels,rdels)
    result = Delta(self.ix,self.iy,self.iz,self.xx,self.yy,self.zz,
                   choose(c,(self.dels,rdels)),
                   choose(c,(self.vs  ,right.vs)),
                   choose(c,(self.ns  ,right.ns)),neumann=right.neumann)
    # --- This is a kludgy fix for problems with subtracting elements.
    # --- If the subtractee has surfaces in common with the subtractor,
    # --- the above algorithm leaves a zero thickness shell there.
    # --- This fix effectively removes those points in common so the shell
    # --- does not appear.
    result.dels = where(abs(self.dels - right.dels)<1.e-10,largepos,result.dels)
    return result

  def __str__(self):
    "Prints out delta"
    return repr(self.dels)+" "+repr(self.vs)+" "+repr(self.ns)

##############################################################################

class Distance:
  """
Class to hold the distance between points and a conductor
Distances have the sign of the outward normal surface vector, i.e.
distances outside the surface are positive, inside negative.
The attribute 'distance' holds the calculated distance.
  """

  def __init__(self,xx=None,yy=None,zz=None,
                    distance=None,generator=None,kwlist=[]):
    if generator is not None:
      self.ndata = len(xx)
      self.xx = xx
      self.yy = yy
      self.zz = zz
      self.distance = zeros(self.ndata,'d')
      apply(generator,kwlist + [self.ndata,self.xx,self.yy,self.zz,
                                self.distance[:]])
    else:
      self.ndata = len(xx)
      self.xx = xx
      self.yy = yy
      self.zz = zz
      self.distance = distance

  def __neg__(self):
    "Distance not operator."
    return Distance(self.xx,self.yy,self.zz, -self.distance)

  def __mul__(self,right):
    "'and' operator, returns maximum of distances to surfaces."
    c = less(self.distance,right.distance)
    return Distance(self.xx,self.yy,self.zz,
                    choose(c,(self.distance,right.distance)))

  def __add__(self,right):
    "'or' operator, returns minimum of distances to surfaces."
    c = greater(self.distance,right.distance)
    dd = Distance(self.xx,self.yy,self.zz,
                  choose(c,(self.distance,right.distance)))
    dd.distance = where((self.distance < 0.) & (right.distance >= 0.),
                        self.distance,dd.distance)
    dd.distance = where((self.distance >= 0.) & (right.distance < 0.),
                        right.distance,dd.distance)
    return dd

  def __sub__(self,right):
    "'or' operator, returns minimum of distances to surfaces."
    # --- Warning - while this work for many cases, there is no
    # --- gaurantee of robustness! It should only work when the right
    # --- hand side is a cylinder.
    rdistance = -right.distance
    dd = Distance(self.xx,self.yy,self.zz,self.distance)
    dd.distance = where((rdistance >= 0.) & (self.distance >= 0.),
                        sqrt(rdistance**2+self.distance**2),
                        dd.distance)
    dd.distance = where((rdistance >= 0.) & (self.distance <= 0.),
                        rdistance,dd.distance)
    dd.distance = where((rdistance < 0.) & (self.distance <= 0.),
                        maximum(rdistance,self.distance),dd.distance)
#   # --- This is an alternate solution that gaurantees that distance will
#   # --- have the correct sign.
#   c = greater(abs(self.distance),abs(right.distance))
#   dd = Distance(self.xx,self.yy,self.zz,
#                 choose(c,(self.distance,right.distance)))
#   dd.distance = where((self.distance < 0.) & (right.distance > 0.),
#                       maximum(-right.distance,self.distance),dd.distance)
#   dd.distance = where((self.distance < 0.) & (right.distance <= 0.),
#                       -right.distance,dd.distance)
#   dd.distance = where((self.distance >= 0.) & (right.distance < 0.),
#                       -right.distance,dd.distance)
#   dd.distance = where((self.distance >= 0.) & (right.distance > 0.),
#                       self.distance,dd.distance)
#   dd.distance = where((self.distance == 0.) & (right.distance == 0.),
#                       1.,dd.distance)
    return dd


  def __str__(self):
    "Prints out delta"
    return repr(self.distance)

##############################################################################

class IsInside:
  """
Class to hold flag whether or not a point is inside a conductor.
The attribute 'distance' holds the calculated distance.
The attribute 'isinside' holds the flag specifying whether a point is in or out
  """

  def __init__(self,xx=None,yy=None,zz=None,
                    isinside=None,generator=None,condid=1,aura=0.,kwlist=[]):
    self.condid = condid
    if generator is not None:
      self.ndata = len(xx)
      self.xx = xx
      self.yy = yy
      self.zz = zz
      self.condid = condid
      self.aura = aura
      distance = zeros(self.ndata,'d')
      apply(generator,kwlist + [self.ndata,self.xx,self.yy,self.zz,
                                distance[:]])
      self.isinside = where(distance <= aura,condid,0.)
    else:
      self.ndata = len(xx)
      self.xx = xx
      self.yy = yy
      self.zz = zz
      self.aura = aura
      self.isinside = isinside*self.condid

  def __neg__(self):
    "IsInside not operator."
    return IsInside(self.xx,self.yy,self.zz,
                    logical_not(self.isinside),
                    condid=self.condid,aura=self.aura)

  def __mul__(self,right):
    "'and' operator, returns logical and of isinsides."
    return IsInside(self.xx,self.yy,self.zz,
                    logical_and(self.isinside,right.isinside),
                    condid=self.condid,aura=self.aura)

  def __add__(self,right):
    "'or' operator, returns logical or of isinsides."
    return IsInside(self.xx,self.yy,self.zz,
                    logical_or(self.isinside,right.isinside),
                    condid=self.condid,aura=self.aura)

  def __sub__(self,right):
    "'or' operator, returns logical or of isinsides."
    return IsInside(self.xx,self.yy,self.zz,
                    logical_and(self.isinside,logical_not(right.isinside)),
                    condid=self.condid,aura=self.aura)

  def __str__(self):
    "Prints out delta"
    return repr(self.isinside)

##############################################################################

class Intercept:
  """
Class to hold information about where a trajectory intercepted a conductor.
The attributes xi, yi, and zi, hold the calculated intersection location.
The attribute angle holds the angle relative to the normal of the surface
at the intersection point.
  """

  def __init__(self,xx=None,yy=None,zz=None,vx=None,vy=None,vz=None,
                    xi=None,yi=None,zi=None,itheta=None,iphi=None,
                    generator=None,condid=1,conductor=None,kwlist=[]):
    self.condid = condid
    self.conductor = conductor
    self.ndata = len(xx)
    self.xx = xx
    self.yy = yy
    self.zz = zz
    self.vx = vx
    self.vy = vy
    self.vz = vz
    if generator is not None:
      self.xi = zeros(self.ndata,'d')
      self.yi = zeros(self.ndata,'d')
      self.zi = zeros(self.ndata,'d')
      self.itheta = zeros(self.ndata,'d')
      self.iphi = zeros(self.ndata,'d')
      apply(generator,kwlist + [self.ndata,self.xx,self.yy,self.zz,
                                self.vx,self.vy,self.vz,
                                self.xi,self.yi,self.zi,self.itheta,self.iphi])
    else:
      self.xi = xi
      self.yi = yi
      self.zi = zi
      self.itheta = itheta
      self.iphi = iphi

  def __neg__(self):
    "Intercept not operator."
    return Intercept(self.xx,self.yy,self.zz,self.vx,self.vy,self.vz,
                     self.xi,self.yi,self.zi,self.itheta+pi,self.iphi,
                     conductor=self.conductor,condid=self.condid)

  def magsq(self):
    return (self.xx-self.xi)**2 + (self.yy-self.yi)**2 + (self.zz-self.zi)**2

  def binaryop(self,right,cond,addpi):
    """All binary operations are the same, and depend only on the conductor
obtained by applying the same operation on the self and right conductors.
    """
    # --- Fuzz value - points less than this from the surface are
    # --- considered to be on the surface.
    surffuzz = 1.e-9
    # --- Get distances from combined conductor to the intercept points
    # --- of each part. Extract absolute value of the distance.
    selfdist = cond.distance(self.xi,self.yi,self.zi)
    rightdist = cond.distance(right.xi,right.yi,right.zi)
    si = abs(selfdist.distance)
    ri = abs(rightdist.distance)
    # --- Get distances between original point and intercept points.
    ds = self.magsq()
    dr = right.magsq()
    # --- Choose intercept points whose distance to the combined conductor
    # --- is the smallest (i.e. get points that lie on its surface).
    cc = si < ri
    # --- If the intercepts points from both conductors lie on the surface,
    # --- Choose the one closest to the original point.
    cc = where((si < surffuzz) & (ri < surffuzz),ds < dr,cc)
    # --- Pick the intercept point which satisfies the above criteria.
    xi = where(cc,self.xi,right.xi)
    yi = where(cc,self.yi,right.yi)
    zi = where(cc,self.zi,right.zi)
    itheta = where(cc,self.itheta,right.itheta+addpi)
    iphi = where(cc,self.iphi,right.iphi)
    # --- Check for cases where neither point lies on the surface. This is
    # --- when the distances to the surface are greater than the fuzz value.
    dd = (minimum(si,ri) < surffuzz)
    xi = where(dd,xi,largepos)
    yi = where(dd,yi,largepos)
    zi = where(dd,zi,largepos)
    itheta = where(dd,itheta,0.)
    iphi = where(dd,iphi,0.)
    return Intercept(self.xx,self.yy,self.zz,self.vx,self.vy,self.vz,
                     xi,yi,zi,itheta,iphi,conductor=cond,condid=self.condid)

  def __mul__(self,right):
    "'and' operator, returns logical and."
    cond = self.conductor*right.conductor
    return self.binaryop(right,cond,0.)

  def __add__(self,right):
    "'or' operator, returns logical or."
    cond = self.conductor + right.conductor
    return self.binaryop(right,cond,0.)

  def __sub__(self,right):
    "'or' operator, returns logical or."
    cond = self.conductor - right.conductor
    return self.binaryop(right,cond,pi)

  def __str__(self):
    "Prints out delta"
    return (repr(self.xi)+' '+repr(self.yi)+' '+repr(self.zi)+' '+
            repr(self.itheta)+' '+repr(self.iphi))

##############################################################################
##############################################################################
##############################################################################

# This is extremely kludgey. The z grid cell size is saved in this variable
# by the getdata function in Grid and then picked up in the Srfrv routines
# to pass into the fortran. There is not simple way of passing this
# information into the conductor object.
_griddzkludge = [0.]

##############################################################################
class Grid:
  """
Class holding the grid info.
Constructor arguments:
  - xmin,xmax,ymin,ymax,zmin,zmax: extent of conductors. Defaults to the
    mesh size. These only need to be set for optimization, to avoid looking
    for conductors where there are none. They can also be used to crop a
    conductor. Note that these are relative to the lab frame.
  - zbeam=top.zbeam: location of grid frame relative to lab frame
  - nx,ny,nz: Number of grid cells in the mesh. Defaults to values from w3d
  - xmmin,xmmax,ymmin,ymmax,zmmin,zmmax: extent of mesh. Defaults to values
                                         from w3d
  - zscale=1.: scale factor on dz. This is used when the relativistic scaling
              is done for the longitudinal dimension. This is only applied
              when getting the coarsening levels.
  - l2symtry,l4symtry: assumed transverse symmetries. Defaults to values
                       from w3d
  - gridrz: RZ grid block to consider
  - my_index,nslaves,izslave,nzslave: information for parallelization
  - solver=w3d: object from which to get the grid size information
Call getdata(a,dfill) to generate the conductor data. 'a' is a geometry object.
Call installdata(installrz,gridmode) to install the data into the WARP database.
  """

  def __init__(self,xmin=None,xmax=None,ymin=None,ymax=None,
                    zmin=None,zmax=None,zbeam=None,
                    nx=None,ny=None,nzlocal=None,nz=None,
                    xmmin=None,xmmax=None,ymmin=None,ymmax=None,
                    zmmin=None,zmmax=None,zscale=1.,
                    l2symtry=None,l4symtry=None,
                    installrz=None,gridrz=None,
                    my_index=None,nslaves=None,izslave=None,nzslave=None,
                    solver=None,mgmaxlevels=None):
    """
Creates a grid object which can generate conductor data.
    """
    _default = lambda x,d: (x,d)[x is None]
    self.zbeam = zbeam
    if self.zbeam is None: zbeam = top.zbeam
    else:                  zbeam = self.zbeam

    if solver is None:
      solver = w3d
      solvertop = top
    else:
      solvertop = solver

    self.nx = _default(nx,solver.nx)
    self.ny = _default(ny,solver.ny)
    self.nz = _default(nz,solver.nz)
    self.nzlocal = _default(nzlocal,solver.nzlocal)
    if self.nzlocal == 0: self.nzlocal = self.nz
    if nz is not None and nzlocal is None and not lparallel:
      self.nzlocal = self.nz
    self.xmmin = _default(xmmin,solver.xmmin)
    self.ymmin = _default(ymmin,solver.ymmin)
    self.zmmin = _default(zmmin,solver.zmmin)
    self.xmmax = _default(xmmax,solver.xmmax)
    self.ymmax = _default(ymmax,solver.ymmax)
    self.zmmax = _default(zmmax,solver.zmmax)
    self.l2symtry = _default(l2symtry,solver.l2symtry)
    self.l4symtry = _default(l4symtry,solver.l4symtry)
    self.my_index = _default(my_index,solvertop.my_index)
    self.nslaves = _default(nslaves,solvertop.nslaves)
    self.izslave = _default(izslave,solvertop.izfsslave)
    self.nzslave = _default(nzslave,solvertop.nzfsslave)

    self.xmin = _default(xmin,self.xmmin)
    self.xmax = _default(xmax,self.xmmax)
    self.ymin = _default(ymin,self.ymmin)
    self.ymax = _default(ymax,self.ymmax)
    self.zmin = _default(zmin,self.zmmin+zbeam)
    self.zmax = _default(zmax,self.zmmax+zbeam)

    self.mgmaxlevels = mgmaxlevels

    # --- Check for symmetries
    if self.l2symtry:
      if self.ymin  < 0.: self.ymin  = 0.
      if self.ymmin < 0.: self.ymmin = 0.
    elif self.l4symtry:
      if self.xmin  < 0.: self.xmin  = 0.
      if self.xmmin < 0.: self.xmmin = 0.
      if self.ymin  < 0.: self.ymin  = 0.
      if self.ymmin < 0.: self.ymmin = 0.

    # --- Note that for the parallel version, the values of zmmin and zmmax
    # --- will be wrong if this is done before the generate, during which time
    # --- the decomposition is done.

    # --- Calculate dx, dy, and dz in case this is called before
    # --- the generate.
    self.dx = (self.xmmax - self.xmmin)/self.nx
    if self.ny > 0: self.dy = (self.ymmax - self.ymmin)/self.ny
    else:           self.dy = self.dx
    # --- z is different since it is not affected by transverse symmetries
    # --- but is affected by parallel decomposition.
    if self.nz > 0: self.dz = (self.zmmax - self.zmmin)/self.nz
    else:           self.dz = (self.zmmax - self.zmmin)
    #if w3d.solvergeom==w3d.XYgeom:self.dz=1.

    # --- Check if frz.basegrid is allocated if installrz is set.
    # --- If not, then turn off installrz.
    if installrz is None:
      installrz = (frz.getpyobject('basegrid') is not None)

    if w3d.solvergeom not in [w3d.RZgeom,w3d.XZgeom,w3d.XYgeom] or not installrz:
      conductors = ConductorType()
      if self.ny > 0: ny = self.ny
      else:           ny = self.nx
      getmglevels(self.nx,ny,self.nzlocal,self.nz,self.dx,self.dy,self.dz*zscale,
                  conductors,
                  self.my_index,self.nslaves,self.izslave,self.nzslave)
      self.mglevels = conductors.levels
      self.mgleveliz = conductors.leveliz[:self.mglevels].copy()
      self.mglevelnz = conductors.levelnz[:self.mglevels].copy()
      self.mglevellx = conductors.levellx[:self.mglevels].copy()
      self.mglevelly = conductors.levelly[:self.mglevels].copy()
      self.mglevellz = conductors.levellz[:self.mglevels].copy()
    else:
      if gridrz is None:gridrz=frz.basegrid
      setmglevels_rz(gridrz)
      self.mglevels = f3d.mglevels
      self.mgleveliz = f3d.mglevelsiz[:f3d.mglevels].copy()
      self.mglevelnz = f3d.mglevelsnz[:f3d.mglevels].copy()
      self.mglevellx = f3d.mglevelslx[:f3d.mglevels].copy()
      self.mglevelly = f3d.mglevelsly[:f3d.mglevels].copy()
      self.mglevellz = f3d.mglevelslz[:f3d.mglevels].copy()

    if self.mgmaxlevels is not None:
      self.mglevels = self.mgmaxlevels

    # --- Create empty lists of conductors
    self.dlist = []
    self.dlistinstalled = []

  def getmeshsize(self,mglevel=0):
    dx = self.dx*self.mglevellx[mglevel]
    dy = self.dy*self.mglevelly[mglevel]
    dz = self.dz*self.mglevellz[mglevel]
    nx = nint(self.nx/self.mglevellx[mglevel])
    ny = nint(self.ny/self.mglevelly[mglevel])
    iz = self.mgleveliz[mglevel]
    nzlocal = self.mglevelnz[mglevel]
    return dx,dy,dz,nx,ny,nzlocal,iz

  def getmesh(self,mglevel=0,extent=None):
    dx,dy,dz,nx,ny,nzlocal,iz = self.getmeshsize(mglevel)
    _griddzkludge[0] = dz

    if self.zbeam is None: zbeam = top.zbeam
    else:                  zbeam = self.zbeam

    xmin,ymin = self.xmin,self.ymin
    xmax,ymax = self.xmax,self.ymax
    zmin = max(self.zmin,self.zmmin+iz*dz+zbeam)
    zmax = min(self.zmax,self.zmmin+(iz+nzlocal)*dz+zbeam)
    if extent is not None:
      xmin,ymin,zmin = maximum(array(extent.mins),array([xmin,ymin,zmin]))
      xmax,ymax,zmax = minimum(array(extent.maxs),array([xmax,ymax,zmax]))

    # --- The conductor extent is completely outside the grid
    if xmin-dx > xmax or ymin-dy > ymax or zmin-dz > zmax:
      return [],[],[],[],[],[],0.,0.,0.,0.,0,0,0,[],0

    zmmin = self.zmmin + iz*dz

    xmesh = self.xmmin + dx*arange(nx+1)
    ymesh = self.ymmin + dy*arange(ny+1)
    zmesh =      zmmin + dz*arange(nzlocal+1) + zbeam
    xmesh = compress(logical_and(xmin-dx <= xmesh,xmesh <= xmax+dx),xmesh)
    ymesh = compress(logical_and(ymin-dy <= ymesh,ymesh <= ymax+dy),ymesh)
    zmesh = compress(logical_and(zmin-dz <= zmesh,zmesh <= zmax+dz),zmesh)
    x = ravel(xmesh[:,NewAxis]*ones(len(ymesh)))
    y = ravel(ymesh*ones(len(xmesh))[:,NewAxis])
    z = zeros(len(xmesh)*len(ymesh),'d')
    ix = nint((x - self.xmmin)/dx)
    iy = nint((y - self.ymmin)/dy)
    iz = zeros(len(xmesh)*len(ymesh),'l')
    return ix,iy,iz,x,y,z,zmmin,dx,dy,dz,nx,ny,nzlocal,zmesh,zbeam

  def checkoverlap(self,mglevel,extent):
    if extent is None: return 1

    dx,dy,dz,nx,ny,nzlocal,iz = self.getmeshsize(mglevel)

    xmin,ymin = self.xmin,self.ymin
    xmax,ymax = self.xmax,self.ymax
    if lparallel:
      if self.zbeam is None: zbeam = top.zbeam
      else:                  zbeam = self.zbeam
      zmin = max(self.zmin,self.zmmin+iz*dz+zbeam)
      zmax = min(self.zmax,self.zmmin+(iz+nzlocal)*dz+zbeam)
    else:
      zmin = self.zmin
      zmax = self.zmax

    xmin,ymin,zmin = maximum(array(extent.mins),array([xmin,ymin,zmin]))
    xmax,ymax,zmax = minimum(array(extent.maxs),array([xmax,ymax,zmax]))

    if zmin-dz > zmax or xmin-dx > xmax or ymin-dy > ymax: return 0
    return 1

  def getdata(self,a,dfill=2.,fuzzsign=-1):
    """
Given an Assembly, accumulate the appropriate data to represent that
Assembly on this grid.
 - a: the assembly or a list of assemblies
 - dfill=2.: points at a depth in the conductor greater than dfill
             are skipped.
    """

    # --- If 'a' is a list, then recursively call this routine for each
    # --- element of 'a'. Note that this will be recursive if some elements
    # --- of 'a' are themselves lists.
    if type(a) == ListType:
      for c in a: self.getdata(c,dfill=dfill,fuzzsign=fuzzsign)
      return

    # --- Check if total conductor overlaps with the grid. If it doesn't,
    # --- then there is no need to check any individual pieces.
    if not lparallel:
      aextent = a.getextent()
      if not self.checkoverlap(0,aextent): return
    else:
      # --- Note that for the parallel version, each level must be checked
      # --- since the levels can have a difference z extent.
      doesoverlap = 0
      aextent = a.getextent()
      for i in range(self.mglevels):
        if self.checkoverlap(i,aextent): doesoverlap = 1
      if not doesoverlap: return

    # --- If 'a' is an AssemblyPlus, save time by generating the conductor
    # --- data for each part separately. Time is saved since only data within
    # --- the extent of each part is checked. Note that this will be recursive
    # --- one of the parts of 'a' are themselves an AssemblyPlus.
    if a.__class__ == AssemblyPlus:
      self.getdata(a.left,dfill=dfill,fuzzsign=fuzzsign)
      self.getdata(a.right,dfill=dfill,fuzzsign=fuzzsign)
      return

    starttime = wtime()
    timeit = 0
    if timeit: tt1 = wtime()
    if timeit: tt2 = zeros(10,'d')
    aextent = a.getextent()
    # --- Leave dall empty. It is created if needed below.
    dall = None
    if timeit: tt2[8] = tt2[8] + wtime() - tt1
    for i in range(self.mglevels):
      if timeit: tt1 = wtime()
      ix,iy,iz,x,y,z,zmmin,dx,dy,dz,nx,ny,nzlocal,zmesh,zbeam=self.getmesh(i,aextent)
      if timeit: tt2[0] = tt2[0] + wtime() - tt1
      if len(x) == 0: continue
      for zz in zmesh:
        if timeit: tt1 = wtime()
        z[:] = zz
        iz[:] = nint((zz - zmmin - zbeam)/dz)
        if timeit: tt2[1] = tt2[1] + wtime() - tt1
        if timeit: tt1 = wtime()
        d = a.griddistance(ix,iy,iz,x,y,z)
        if timeit: tt2[2] = tt2[2] + wtime() - tt1
        if timeit: tt1 = wtime()
        d.normalize(dx,dy,dz)
        if timeit: tt2[3] = tt2[3] + wtime() - tt1
        if timeit: tt1 = wtime()
        d.setparity(dfill,fuzzsign)
        if timeit: tt2[4] = tt2[4] + wtime() - tt1
        if timeit: tt1 = wtime()
        d.clean()
        if timeit: tt2[5] = tt2[5] + wtime() - tt1
        if timeit: tt1 = wtime()
        d.setlevels(i)
        if timeit: tt2[6] = tt2[6] + wtime() - tt1
        if timeit: tt1 = wtime()
        if dall is None:
          # --- Only create the Delta instance if it is actually needed.
          dall = Delta(neumann=a.neumann)
          self.dlist.append(dall)
        dall.append(d)
        if timeit: tt2[7] = tt2[7] + wtime() - tt1
    endtime = wtime()
    self.generatetime = endtime - starttime
    if timeit: tt2[9] = endtime - starttime
    if timeit: print tt2

  def installdata(self,installrz=1,gridmode=1,solvergeom=None,
                  conductors=f3d.conductors,gridrz=None):
    """
Installs the conductor data into the fortran database
    """
    conductors.levels = self.mglevels
    conductors.leveliz[:self.mglevels] = self.mgleveliz[:self.mglevels]
    conductors.levelnz[:self.mglevels] = self.mglevelnz[:self.mglevels]
    conductors.levellx[:self.mglevels] = self.mglevellx[:self.mglevels]
    conductors.levelly[:self.mglevels] = self.mglevelly[:self.mglevels]
    conductors.levellz[:self.mglevels] = self.mglevellz[:self.mglevels]
    for d in self.dlist:
      d.install(installrz,solvergeom,conductors,gridrz)
    self.dlistinstalled += self.dlist
    self.dlist = []
    if gridmode is not None:
      f3d.gridmode = gridmode

  def getdistances(self,a,mglevel=0):
    """
Given an Assembly, accumulate the distances between the assembly and the
grid points.
 - a: the assembly
 - mglevel=0: coarsening level to use
    """
    starttime = wtime()
    tt2 = zeros(4,'d')
    tt1 = wtime()
    ix,iy,iz,x,y,z,zmmin,dx,dy,dz,nx,ny,nzlocal,zmesh,zbeam = self.getmesh(mglevel)
    try:
      self.distances[0,0,0]
    except AttributeError:
      self.distances = fzeros((1+nx,1+ny,1+nzlocal),'d')
    ix1 = min(ix)
    ix2 = max(ix)
    iy1 = min(iy)
    iy2 = max(iy)
    tt2[0] = tt2[0] + wtime() - tt1
    if len(x) == 0: return
    for zz in zmesh:
      tt1 = wtime()
      z[:] = zz
      iz[:] = nint((zz - zmmin - zbeam)/dz)
      tt2[1] = tt2[1] + wtime() - tt1
      tt1 = wtime()
      d = a.distance(x,y,z)
      tt2[2] = tt2[2] + wtime() - tt1
      tt1 = wtime()
      dd = d.distance
      dd.shape = (ix2-ix1+1,iy2-iy1+1)
      self.distances[ix1:ix2+1,iy1:iy2+1,iz[0]] = dd
      tt2[3] = tt2[3] + wtime() - tt1
    endtime = wtime()
    self.generatetime = endtime - starttime
    #print tt2

  def resetisinside(self,mglevel=0):
    """
Clears out any data in isinside by recreating the array.
    """
    ix,iy,iz,x,y,z,zmmin,dx,dy,dz,nx,ny,nzlocal,zmesh,zbeam = self.getmesh(mglevel)
    self.isinside = fzeros((1+nx,1+ny,1+nzlocal),'d')

  def removeisinside(self,a,nooverlap=0):
    """
Removes an assembly from the isinside data. This assumes that that assembly
does not overlap any others. The flag nooverlap must be set as a reminder.
    """
    if nooverlap:
      self.isinside = where(self.isinside==a.condid,0,self.isinside)
    else:
      print "removeisinside only works when the assembly does not overlap any others"
      print 'Set the nooverlap flag to true is this is the case.'
      raise ''

  def getisinside(self,a,mglevel=0,aura=0.):
    """
Given an Assembly, set flag for each grid point whether it is inside the
assembly.
 - a: the assembly
 - mglevel=0: coarsening level to use
    """
    starttime = wtime()
    tt2 = zeros(4,'d')
    tt1 = wtime()
    ix,iy,iz,x,y,z,zmmin,dx,dy,dz,nx,ny,nzlocal,zmesh,zbeam = self.getmesh(mglevel)
    try:
      self.isinside[0,0,0]
    except AttributeError:
      self.resetisinside(mglevel)
    ix1 = min(ix)
    ix2 = max(ix)
    iy1 = min(iy)
    iy2 = max(iy)
    tt2[0] = tt2[0] + wtime() - tt1
    if len(x) == 0: return
    for zz in zmesh:
      tt1 = wtime()
      z[:] = zz #####
      iz[:] = nint((zz - zmmin - zbeam)/dz)  #####
      tt2[1] = tt2[1] + wtime() - tt1
      tt1 = wtime()
      d = a.isinside(x,y,z,aura)  #####
      tt2[2] = tt2[2] + wtime() - tt1
      tt1 = wtime()
      dd = d.isinside
      dd.shape = (ix2-ix1+1,iy2-iy1+1)
      self.isinside[ix1:ix2+1,iy1:iy2+1,iz[0]] = where(dd>0,dd,
                                   self.isinside[ix1:ix2+1,iy1:iy2+1,iz[0]])
      tt2[3] = tt2[3] + wtime() - tt1
    endtime = wtime()
    self.generatetime = endtime - starttime
    #print tt2

##############################################################################
##############################################################################
##############################################################################

#============================================================================
class Plane(Assembly):
  """
Plane class
  - z0=0: locate of plane relative to zcent
  - zsign=1: when positive, conductor is in the z>0 side
  - theta=0,phi=0: normal of surface defining plane relative to z-axis
    theta is angle in z-x plane
    phi is angle in z-y plane
  - voltage=0: box voltage
  - xcent=0.,ycent=0.,zcent=0.: center of box
  - condid=1: conductor id of box, must be integer, or can be 'next' in which
              case a unique ID is chosen
  """
  def __init__(self,z0=0.,zsign=1.,theta=0.,phi=0.,
                    voltage=0.,xcent=0.,ycent=0.,zcent=0.,
                    condid=1,**kw):
    kwlist=['z0','zsign','theta','phi']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                           planeconductorf,planeconductord,planeintercept,
                           kw=kw)
    self.z0 = z0
    self.zsign = zsign
    self.theta = theta
    self.phi = phi
    if self.theta == 0. and self.phi == 0.: z1,z2 = z0,z0
    else:                                   z1,z2 = -largepos,+largepos
    self.createextent([-largepos,-largepos,z1],[+largepos,+largepos,z2])

#============================================================================
class Box(Assembly):
  """
Box class
  - xsize,ysize,zsize: box size
  - voltage=0: box voltage
  - xcent=0.,ycent=0.,zcent=0.: center of box
  - condid=1: conductor id of box, must be integer, or can be 'next' in which
              case a unique ID is chosen
  """
  def __init__(self,xsize,ysize,zsize,voltage=0.,xcent=0.,ycent=0.,zcent=0.,
                    condid=1,**kw):
    kwlist=['xsize','ysize','zsize']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                           boxconductorf,boxconductord,boxintercept,
                           kw=kw)
    self.xsize = xsize
    self.ysize = ysize
    self.zsize = zsize
    self.createextent([-self.xsize/2.,-self.ysize/2.,-self.zsize/2.],
                      [+self.xsize/2.,+self.ysize/2.,+self.zsize/2.])

  def draw(self,color='fg',filled=None,fullplane=1,**kw):
    r = self.xsize/2.*array([-1,+1,+1,-1,-1])
    z = self.zsize/2.*array([-1,-1,+1,+1,-1])
    self.plotdata(r,z,color=color,filled=filled,fullplane=fullplane)

  def createdxobject(self,kwdict={},**kw):
    kw.update(kwdict)
    v = Opyndx.VisualBox(self.xsize,self.ysize,self.zsize,
                         self.xcent,self.ycent,self.zcent,
                         kwdict=kw)
    self.dxobject = v

#============================================================================
class Cylinder(Assembly):
  """
Cylinder class
  - radius,length: cylinder size
  - theta=0,phi=0: angle of cylinder relative to z-axis
    theta is angle in z-x plane
    phi is angle in z-y plane
  - voltage=0: cylinder voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cylinder
  - condid=1: conductor id of cylinder, must be integer, or can be 'next' in
              which case a unique ID is chosen
  """
  def __init__(self,radius,length,theta=0.,phi=0.,
                    voltage=0.,xcent=0.,ycent=0.,zcent=0.,condid=1,**kw):
    kwlist = ['radius','length','theta','phi']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                           cylinderconductorf,cylinderconductord,
                           cylinderintercept,
                           kw=kw)
    self.radius = radius
    self.length = length
    self.theta  = theta
    self.phi    = phi

    # --- This is the easiest thing to do without thinking.
    ll = sqrt(self.radius**2 + (self.length/2.)**2)
    self.createextent([-ll,-ll,-ll],[+ll,+ll,+ll])

  def draw(self,color='fg',filled=None,fullplane=1,**kw):
    """
Plots the r versus z
 - color='fg': color of outline, set to None to not plot the outline
 - filled=None: when set to an integer, fills the outline with the color
                specified from the current palette. Should be between 0 and 199.
 - fullplane=1: when true, plot the top and bottom, i.e. r vs z, and -r vs z.
 - rmin=0.: inner range in r to include in plot
    """
    # --- This is kind of a hack, but this routine doesn't make much sense
    # --- for an arbitrarily rotated cylinder.
    r = array([self.radius,self.radius,-self.radius,-self.radius,self.radius])
    z = self.length*array([-0.5,0.5,0.5,-0.5,-0.5])

    ct = cos(self.theta)
    st = sin(self.theta)
    cp = cos(self.phi)
    sp = sin(self.phi)
    xp = +r*ct - 0*st*sp + z*st*cp
    yp =       + 0*cp    + z*sp
    zp = -r*st - 0*ct*sp + z*ct*cp

    self.plotdata(xp,zp,color=color,filled=filled,fullplane=0)

  def createdxobject(self,kwdict={},**kw):
    kw.update(kwdict)
    v = Opyndx.VisualRevolution(
                       zzmin=-self.length/2.,zzmax=+self.length/2.,
                       rendzmin=0.,rendzmax=0.,
                       xoff=self.xcent,yoff=self.ycent,zoff=self.zcent,
                       theta=self.theta,phi=self.phi,
                       rofzdata=[self.radius,self.radius],
                       zdata=[-self.length/2.,+self.length/2.],
                       largepos=largepos,
                       kwdict=kw)
    self.dxobject = v

#============================================================================
class Cylinders(Assembly):
  """
Cylinders class for a list of cylinders
  - radius,length: cylinder size
  - theta=0,phi=0: angle of cylinder relative to z-axis
    theta is angle in z-x plane
    phi is angle in z-y plane
  - voltage=0: cylinder voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cylinder
  - condid=1: conductor id of cylinder, must be integer, or can be 'next' in
              which case a unique ID is chosen
  """
  def __init__(self,radius,length,theta=0.,phi=0.,
                    voltage=0.,xcent=0.,ycent=0.,zcent=0.,condid=1,**kw):
    kwlist = ['ncylinders','radius','length','theta','phi']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                           cylindersconductorf,cylindersconductord,
                           cylindersintercept,
                           kw=kw)
    self.ncylinders = 0
    self.radius = radius
    self.length = length
    self.theta  = theta
    self.phi    = phi
    # --- Find the first input argument that has a length > 1. This is done
    # --- since some of the inputs may be arrays of length 1, which is OK,
    # --- but does not mean that there is only 1 cone. Some other input
    # --- may be a longer array. getkwlist is used so that the xcent etc.
    # --- are included.
    kwlist = self.getkwlist()
    for k in kwlist:
      try:
        self.ncylinders = len(k)
        if self.ncylinders > 1: break
      except TypeError:
        pass

    assert self.ncylinders > 0,"At least on of the input arguments must be a list!"
    self.radius = self.radius*ones(self.ncylinders)
    self.length = self.length*ones(self.ncylinders)
    self.theta  = self.theta*ones(self.ncylinders)
    self.phi    = self.phi*ones(self.ncylinders)
    self.xcent  = self.xcent*ones(self.ncylinders)
    self.ycent  = self.ycent*ones(self.ncylinders)
    self.zcent  = self.zcent*ones(self.ncylinders)

    rmax = sqrt(self.radius**2 + (self.length/2.)**2)
    xmin = min(self.xcent-rmax)
    ymin = min(self.ycent-rmax)
    zmin = min(self.zcent-rmax)
    xmax = max(self.xcent+rmax)
    ymax = max(self.ycent+rmax)
    zmax = max(self.zcent+rmax)
    self.extent = ConductorExtent([xmin,ymin,zmin],[xmax,ymax,zmax])

#============================================================================
class ZCylinder(Assembly):
  """
Cylinder aligned with z-axis
  - radius,length: cylinder size
  - voltage=0: cylinder voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cylinder
  - condid=1: conductor id of cylinder, must be integer, or can be 'next' in
              which case a unique ID is chosen
  - zlower,zupper: Optionally, the lower and upper extent of the cylinder
                   can be specified instead of length and zcent. If both are
                   given, then length and zcent are ignored. If only one is
                   given, it is ignored.
  """
  def __init__(self,radius,length=None,voltage=0.,xcent=0.,ycent=0.,zcent=0.,
                    condid=1,zlower=None,zupper=None,**kw):
    if zlower is not None and zupper is not None:
      length = zupper - zlower
      zcent = 0.5*(zlower + zupper)
      self.zlower = zlower
      self.zupper = zupper
    assert length is not None,\
      "ZCylinder: either length or both zlower and zupper must be specified"
    kwlist = ['radius','length']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                           zcylinderconductorf,zcylinderconductord,
                           zcylinderintercept,
                           kw=kw)
    self.radius = radius
    self.length = length
    self.createextent([-self.radius,-self.radius,-self.length/2.],
                      [+self.radius,+self.radius,+self.length/2.])

  def createdxobject(self,kwdict={},**kw):
    kw.update(kwdict)
    v = Opyndx.VisualRevolution(
                       zzmin=-self.length/2.,zzmax=+self.length/2.,
                       rendzmin=0.,rendzmax=0.,
                       xoff=self.xcent,yoff=self.ycent,zoff=self.zcent,
                       rofzdata=[self.radius,self.radius],
                       zdata=[-self.length/2.,+self.length/2.],
                       largepos=largepos,
                       kwdict=kw)
    self.dxobject = v

  def draw(self,color='fg',filled=None,fullplane=1,**kw):
    """
Plots the r versus z
 - color='fg': color of outline, set to None to not plot the outline
 - filled=None: when set to an integer, fills the outline with the color
                specified from the current palette. Should be between 0 and 199.
 - fullplane=1: when true, plot the top and bottom, i.e. r vs z, and -r vs z.
 - rmin=0.: inner range in r to include in plot
    """
    rmin = kw.get('rmin',None)
    if rmin is None: rmin = 0.
    r = [self.radius,self.radius,rmin,rmin,self.radius]
    z = self.length*array([-0.5,0.5,0.5,-0.5,-0.5])
    self.plotdata(r,z,color=color,filled=filled,fullplane=fullplane)

#============================================================================
class ZRoundedCylinder(Assembly):
  """
Cylinder with rounded corners aligned with z-axis
  - radius,length: cylinder size
  - radius2: radius of rounded corners
  - voltage=0: cylinder voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cylinder
  - condid=1: conductor id of cylinder, must be integer, or can be 'next' in
              which case a unique ID is chosen
  """
  def __init__(self,radius,length,radius2,voltage=0.,
                    xcent=0.,ycent=0.,zcent=0.,
                    condid=1,**kw):
    kwlist = ['radius','length','radius2']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      zroundedcylinderconductorf,zroundedcylinderconductord,
                      zroundedcylinderintercept,
                      kw=kw)
    self.radius = radius
    self.length = length
    self.radius2 = radius2
    self.createextent([-self.radius,-self.radius,-self.length/2.],
                      [+self.radius,+self.radius,+self.length/2.])

  def createdxobject(self,kwdict={},**kw):
    kw.update(kwdict)
    rr = [self.radius-self.radius2,self.radius,
          self.radius,self.radius-self.radius2]
    zz = [-self.length/2.,-self.length/2.+self.radius2,
          +self.length/2.-self.radius2,+self.length/2.]
    rad = [self.radius2,None,self.radius2]
    zc = [None,None,None]
    rc = [None,None,None]
    if zz[1] > zz[2]:
      zz[1] = 0.5*(zz[1] + zz[2])
      zz[2] = zz[1]
      rr[1:3] = sqrt(self.radius2**2 - (self.radius2 - self.length/2.)**2)
    Srfrv.checkarcs(Srfrv(),zz,rr,rad,zc,rc)

    v = Opyndx.VisualRevolution(
                       zzmin=-self.length/2.,zzmax=+self.length/2.,
                       rendzmin=0.,rendzmax=0.,
                       xoff=self.xcent,yoff=self.ycent,zoff=self.zcent,
                       rofzdata=rr,zdata=zz,raddata=rad,zcdata=zc,rcdata=rc,
                       largepos=largepos,
                       kwdict=kw)
    self.dxobject = v

#============================================================================
class ZCylinderOut(Assembly):
  """
Outside of a cylinder aligned with z-axis
  - radius,length: cylinder size
  - voltage=0: cylinder voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cylinder
  - condid=1: conductor id of cylinder, must be integer, or can be 'next' in
              which case a unique ID is chosen
  - zlower,zupper: Optionally, the lower and upper extent of the cylinder
                   can be specified instead of length and zcent. If both are
                   given, then length and zcent are ignored. If only one is
                   given, it is ignored.
  """
  def __init__(self,radius,length=None,voltage=0.,xcent=0.,ycent=0.,zcent=0.,
                    condid=1,zlower=None,zupper=None,**kw):
    if zlower is not None and zupper is not None:
      length = zupper - zlower
      zcent = 0.5*(zlower + zupper)
      self.zlower = zlower
      self.zupper = zupper
    assert length is not None,\
      "ZCylinderOut: either length or both zlower and zupper must be specified"
    kwlist = ['radius','length']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      zcylinderoutconductorf,zcylinderoutconductord,
                      zcylinderoutintercept,
                      kw=kw)
    self.radius = radius
    self.length = length
    self.createextent([-largepos,-largepos,-self.length/2.],
                      [+largepos,+largepos,+self.length/2.])

  def createdxobject(self,rend=1.,kwdict={},**kw):
    kw.update(kwdict)
    v = Opyndx.VisualRevolution(
                       zzmin=-self.length/2.,zzmax=+self.length/2.,
                       rendzmin=rend,rendzmax=rend,
                       xoff=self.xcent,yoff=self.ycent,zoff=self.zcent,
                       rofzdata=[self.radius,self.radius],
                       zdata=[-self.length/2.,+self.length/2.],
                       raddata=[largepos],zcdata=[largepos],rcdata=[largepos],
                       normalsign=-1,
                       largepos=largepos,
                       kwdict=kw)
    self.dxobject = v

  def draw(self,color='fg',filled=None,fullplane=1,**kw):
    """
Plots the r versus z
 - color='fg': color of outline, set to None to not plot the outline
 - filled=None: when set to an integer, fills the outline with the color
                specified from the current palette. Should be between 0 and 199.
 - fullplane=1: when true, plot the top and bottom, i.e. r vs z, and -r vs z.
 - rmax=w3d.xmmax: outer range in r to include in plot
    """
    rmax = kw.get('rmax',None)
    if rmax is None: rmax = w3d.xmmax
    r = [self.radius,self.radius,rmax,rmax,self.radius]
    z = self.length*array([-0.5,0.5,0.5,-0.5,-0.5])
    self.plotdata(r,z,color=color,filled=filled,fullplane=fullplane)

#============================================================================
class ZRoundedCylinderOut(Assembly):
  """
Outside of a cylinder with rounded corners aligned with z-axis
  - radius,length: cylinder size
  - radius2: radius of rounded corners
  - voltage=0: cylinder voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cylinder
  - condid=1: conductor id of cylinder, must be integer, or can be 'next' in
              which case a unique ID is chosen
  """
  def __init__(self,radius,length,radius2,voltage=0.,
                    xcent=0.,ycent=0.,zcent=0.,
                    condid=1,**kw):
    kwlist = ['radius','length','radius2']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      zroundedcylinderoutconductorf,
                      zroundedcylinderoutconductord,
                      zroundedcylinderoutintercept,
                      kw=kw)
    self.radius = radius
    self.length = length
    self.radius2 = radius2
    self.createextent([-largepos,-largepos,-self.length/2.],
                      [+largepos,+largepos,+self.length/2.])

  def createdxobject(self,rend=1.,kwdict={},**kw):
    kw.update(kwdict)

    rr = [self.radius+self.radius2,self.radius,
          self.radius,self.radius+self.radius2]
    zz = [-self.length/2.,-self.length/2.+self.radius2,
          +self.length/2.-self.radius2,+self.length/2.]
    rad = [-self.radius2,None,-self.radius2]
    zc = [None,None,None]
    rc = [None,None,None]
    if zz[1] > zz[2]:
      zz[1] = 0.5*(zz[1] + zz[2])
      zz[2] = zz[1]
      rfixed = (self.radius + self.radius2 -
                 sqrt(self.radius2**2 - (self.radius2 - self.length/2.)**2))
      rr[1] = rfixed
      rr[2] = rfixed
    Srfrv.checkarcs(Srfrv(),zz,rr,rad,zc,rc)

    v = Opyndx.VisualRevolution(
                       zzmin=-self.length/2.,zzmax=+self.length/2.,
                       rendzmin=rend,rendzmax=rend,
                       xoff=self.xcent,yoff=self.ycent,zoff=self.zcent,
                       rofzdata=rr,zdata=zz,raddata=rad,zcdata=zc,rcdata=rc,
                       normalsign=-1,
                       largepos=largepos,
                       kwdict=kw)
    self.dxobject = v

#============================================================================
class XCylinder(ZCylinder,XAssembly):
  """
Cylinder aligned with X-axis
  - radius,length: cylinder size
  - voltage=0: cylinder voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cylinder
  - condid=1: conductor id of cylinder, must be integer, or can be 'next' in
              which case a unique ID is chosen
  - xlower,xupper: Optionally, the lower and upper extent of the cylinder
                   can be specified instead of length and xcent. If both are
                   given, then length and xcent are ignored. If only one is
                   given, it is ignored.
  """
  def __init__(self,radius,length=None,voltage=0.,xcent=0.,ycent=0.,zcent=0.,
                    condid=1,xlower=None,xupper=None,**kw):
    if xlower is not None and xupper is not None:
      length = xupper - xlower
      xcent = 0.5*(xlower + xupper)
      self.xlower = xlower
      self.xupper = xupper
    assert length is not None,\
      "XCylinder: either length or both xlower and xupper must be specified"
    ZCylinder.__init__(self,radius,length,
                            voltage,xcent,ycent,zcent,condid=1)
    XAssembly.__init__(self,voltage,xcent,ycent,zcent,condid,self.kwlist,
                            self.generatorf,self.generatord,self.generatori,
                            kw=kw)

#============================================================================
class XCylinderOut(ZCylinderOut,XAssembly):
  """
Cylinder aligned with X-axis
  - radius,length: cylinder size
  - voltage=0: cylinder voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cylinder
  - condid=1: conductor id of cylinder, must be integer, or can be 'next' in
              which case a unique ID is chosen
  - xlower,xupper: Optionally, the lower and upper extent of the cylinder
                   can be specified instead of length and xcent. If both are
                   given, then length and xcent are ignored. If only one is
                   given, it is ignored.
  """
  def __init__(self,radius,length=None,voltage=0.,xcent=0.,ycent=0.,zcent=0.,
                    condid=1,xlower=None,xupper=None,**kw):
    if xlower is not None and xupper is not None:
      length = xupper - xlower
      xcent = 0.5*(xlower + xupper)
      self.xlower = xlower
      self.xupper = xupper
    assert length is not None,\
      "XCylinder: either length or both xlower and xupper must be specified"
    ZCylinderOut.__init__(self,radius,length,
                               voltage,xcent,ycent,zcent,condid=1)
    XAssembly.__init__(self,voltage,xcent,ycent,zcent,condid,self.kwlist,
                            self.generatorf,self.generatord,self.generatori,
                            kw=kw)

#============================================================================
class YCylinder(ZCylinder,YAssembly):
  """
Cylinder aligned with Y-axis
  - radius,length: cylinder size
  - voltage=0: cylinder voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cylinder
  - condid=1: conductor id of cylinder, must be integer, or can be 'next' in
              which case a unique ID is chosen
  - ylower,yupper: Optionally, the lower and upper extent of the cylinder
                   can be specified instead of length and ycent. If both are
                   given, then length and ycent are ignored. If only one is
                   given, it is ignored.
  """
  def __init__(self,radius,length=None,voltage=0.,xcent=0.,ycent=0.,zcent=0.,
                    condid=1,ylower=None,yupper=None,**kw):
    if ylower is not None and yupper is not None:
      length = yupper - ylower
      ycent = 0.5*(ylower + yupper)
      self.ylower = ylower
      self.yupper = yupper
    assert length is not None,\
      "YCylinder: either length or both ylower and yupper must be specified"
    ZCylinder.__init__(self,radius,length,
                            voltage,xcent,ycent,zcent,condid=1)
    YAssembly.__init__(self,voltage,xcent,ycent,zcent,condid,self.kwlist,
                            self.generatorf,self.generatord,self.generatori,
                            kw=kw)

#============================================================================
class YCylinderOut(ZCylinderOut,YAssembly):
  """
Cylinder aligned with Y-axis
  - radius,length: cylinder size
  - voltage=0: cylinder voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cylinder
  - condid=1: conductor id of cylinder, must be integer, or can be 'next' in
              which case a unique ID is chosen
  - ylower,yupper: Optionally, the lower and upper extent of the cylinder
                   can be specified instead of length and ycent. If both are
                   given, then length and ycent are ignored. If only one is
                   given, it is ignored.
  """
  def __init__(self,radius,length=None,voltage=0.,xcent=0.,ycent=0.,zcent=0.,
                    condid=1,ylower=None,yupper=None,**kw):
    if ylower is not None and yupper is not None:
      length = yupper - ylower
      ycent = 0.5*(ylower + yupper)
      self.ylower = ylower
      self.yupper = yupper
    assert length is not None,\
      "YCylinder: either length or both ylower and yupper must be specified"
    ZCylinderOut.__init__(self,radius,length,
                               voltage,xcent,ycent,zcent,condid=1)
    YAssembly.__init__(self,voltage,xcent,ycent,zcent,condid,self.kwlist,
                            self.generatorf,self.generatord,self.generatori,
                            kw=kw)

#============================================================================
class ZCylinderElliptic(ZCylinder,EllipticAssembly):
  """
Elliptical cylinder aligned with z-axis
  - ellipticity: ratio of y radius to x radius
  - radius,length: cylinder size
  - voltage=0: cylinder voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cylinder
  - condid=1: conductor id of cylinder, must be integer, or can be 'next' in
              which case a unique ID is chosen
  - zlower,zupper: Optionally, the lower and upper extent of the cylinder
                   can be specified instead of length and zcent. If both are
                   given, then length and zcent are ignored. If only one is
                   given, it is ignored.
  """
  def __init__(self,ellipticity,radius,length=None,
                    voltage=0.,xcent=0.,ycent=0.,zcent=0.,condid=1,
                    zlower=None,zupper=None,**kw):
    ZCylinder.__init__(self,radius,length,
                            voltage,xcent,ycent,zcent,condid,
                            zlower=zlower,zupper=zupper)
    EllipticAssembly.__init__(self,ellipticity,
                              voltage,xcent,ycent,zcent,condid,self.kwlist,
                              self.generatorf,self.generatord,self.generatori,
                              kw=kw)

  def createdxobject(self,kwdict={},**kw):
    kw.update(kwdict)
    e = self.ellipticity
    v = Opyndx.VisualRevolutionEllipse(
                       zzmin=-self.length/2.,zzmax=+self.length/2.,
                       rxendzmin=0.,rxendzmax=0.,
                       ryendzmin=0.,ryendzmax=0.,
                       xoff=self.xcent,yoff=self.ycent,zoff=self.zcent,
                       rxofzdata=[self.radius,self.radius],
                       ryofzdata=[e*self.radius,e*self.radius],
                       zdata=[-self.length/2.,+self.length/2.],
                       largepos=largepos,
                       kwdict=kw)
    self.dxobject = v

#============================================================================
class ZCylinderEllipticOut(ZCylinderOut,EllipticAssembly):
  """
Outside an elliptical cylinder aligned with z-axis
  - ellipticity: ratio of y radius to x radius
  - radius,length: cylinder size
  - voltage=0: cylinder voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cylinder
  - condid=1: conductor id of cylinder, must be integer, or can be 'next' in
              which case a unique ID is chosen
  - zlower,zupper: Optionally, the lower and upper extent of the cylinder
                   can be specified instead of length and zcent. If both are
                   given, then length and zcent are ignored. If only one is
                   given, it is ignored.
  """
  def __init__(self,ellipticity,radius,length=None,
                    voltage=0.,xcent=0.,ycent=0.,zcent=0.,condid=1,
                    zlower=None,zupper=None,**kw):
    ZCylinderOut.__init__(self,radius,length,
                               voltage,xcent,ycent,zcent,condid,
                               zlower=zlower,zupper=zupper)
    EllipticAssembly.__init__(self,ellipticity,
                              voltage,xcent,ycent,zcent,condid,self.kwlist,
                              self.generatorf,self.generatord,self.generatori,
                              kw=kw)

  def createdxobject(self,rend=1.,kwdict={},**kw):
    kw.update(kwdict)
    e = self.ellipticity
    v = Opyndx.VisualRevolutionEllipse(
                       zzmin=-self.length/2.,zzmax=+self.length/2.,
                       rxendzmin=rend,rxendzmax=rend,
                       ryendzmin=e*rend,ryendzmax=e*rend,
                       xoff=self.xcent,yoff=self.ycent,zoff=self.zcent,
                       rxofzdata=[self.radius,self.radius],
                       ryofzdata=[e*self.radius,e*self.radius],
                       zdata=[-self.length/2.,+self.length/2.],
                       largepos=largepos,
                       kwdict=kw)
    self.dxobject = v

#============================================================================
class XCylinderElliptic(ZCylinder,EllipticAssembly,XAssembly):
  """
Elliptical cylinder aligned with x-axis
  - ellipticity: ratio of z radius to x radius
  - radius,length: cylinder size
  - voltage=0: cylinder voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cylinder
  - condid=1: conductor id of cylinder, must be integer, or can be 'next' in
              which case a unique ID is chosen
  - xlower,xupper: Optionally, the lower and upper extent of the cylinder
                   can be specified instead of length and xcent. If both are
                   given, then length and xcent are ignored. If only one is
                   given, it is ignored.
  """
  def __init__(self,ellipticity,radius,length=None,
                    voltage=0.,xcent=0.,ycent=0.,zcent=0.,condid=1,
                    xlower=None,xupper=None,**kw):
    if xlower is not None and xupper is not None:
      length = xupper - xlower
      xcent = 0.5*(xlower + xupper)
      self.xlower = xlower
      self.xupper = xupper
    assert length is not None,\
      "XCylinder: either length or both xlower and xupper must be specified"
    ZCylinder.__init__(self,radius,length,
                            voltage,xcent,ycent,zcent,condid)
    EllipticAssembly.__init__(self,ellipticity,
                              voltage,xcent,ycent,zcent,condid,self.kwlist,
                              self.generatorf,self.generatord,self.generatori)
    XAssembly.__init__(self,
                       voltage,xcent,ycent,zcent,condid,self.kwlist,
                       self.generatorf,self.generatord,self.generatori,
                       kw=kw)

#============================================================================
class XCylinderEllipticOut(ZCylinderOut,EllipticAssembly,XAssembly):
  """
Outside of an elliptical cylinder aligned with x-axis
  - ellipticity: ratio of z radius to x radius
  - radius,length: cylinder size
  - voltage=0: cylinder voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cylinder
  - condid=1: conductor id of cylinder, must be integer, or can be 'next' in
              which case a unique ID is chosen
  - xlower,xupper: Optionally, the lower and upper extent of the cylinder
                   can be specified instead of length and xcent. If both are
                   given, then length and xcent are ignored. If only one is
                   given, it is ignored.
  """
  def __init__(self,ellipticity,radius,length=None,
                    voltage=0.,xcent=0.,ycent=0.,zcent=0.,condid=1,
                    xlower=None,xupper=None,**kw):
    if xlower is not None and xupper is not None:
      length = xupper - xlower
      xcent = 0.5*(xlower + xupper)
      self.xlower = xlower
      self.xupper = xupper
    assert length is not None,\
      "XCylinder: either length or both xlower and xupper must be specified"
    ZCylinderOut.__init__(self,radius,length,
                               voltage,xcent,ycent,zcent,condid)
    EllipticAssembly.__init__(self,ellipticity,
                              voltage,xcent,ycent,zcent,condid,self.kwlist,
                              self.generatorf,self.generatord,self.generatori)
    XAssembly.__init__(self,
                       voltage,xcent,ycent,zcent,condid,self.kwlist,
                       self.generatorf,self.generatord,self.generatori,
                       kw=kw)

#============================================================================
class YCylinderElliptic(ZCylinder,EllipticAssembly,YAssembly):
  """
Elliptical cylinder aligned with y-axis
  - ellipticity: ratio of x radius to z radius
  - radius,length: cylinder size
  - voltage=0: cylinder voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cylinder
  - condid=1: conductor id of cylinder, must be integer, or can be 'next' in
              which case a unique ID is chosen
  - ylower,yupper: Optionally, the lower and upper extent of the cylinder
                   can be specified instead of length and ycent. If both are
                   given, then length and ycent are ignored. If only one is
                   given, it is ignored.
  """
  def __init__(self,ellipticity,radius,length=None,
                    voltage=0.,xcent=0.,ycent=0.,zcent=0.,condid=1,ylower=None,yupper=None,**kw):
    if ylower is not None and yupper is not None:
      length = yupper - ylower
      ycent = 0.5*(ylower + yupper)
      self.ylower = ylower
      self.yupper = yupper
    assert length is not None,\
      "YCylinder: either length or both ylower and yupper must be specified"
    ZCylinder.__init__(self,radius,length,
                            voltage,xcent,ycent,zcent,condid)
    EllipticAssembly.__init__(self,ellipticity,
                              voltage,xcent,ycent,zcent,condid,self.kwlist,
                              self.generatorf,self.generatord,self.generatori)
    YAssembly.__init__(self,
                       voltage,xcent,ycent,zcent,condid,self.kwlist,
                       self.generatorf,self.generatord,self.generatori,
                       kw=kw)

#============================================================================
class YCylinderEllipticOut(ZCylinderOut,EllipticAssembly,YAssembly):
  """
Outside of an elliptical cylinder aligned with y-axis
  - ellipticity: ratio of x radius to z radius
  - radius,length: cylinder size
  - voltage=0: cylinder voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cylinder
  - condid=1: conductor id of cylinder, must be integer, or can be 'next' in
              which case a unique ID is chosen
  - ylower,yupper: Optionally, the lower and upper extent of the cylinder
                   can be specified instead of length and ycent. If both are
                   given, then length and ycent are ignored. If only one is
                   given, it is ignored.
  """
  def __init__(self,ellipticity,radius,length=None,
                    voltage=0.,xcent=0.,ycent=0.,zcent=0.,condid=1,ylower=None,yupper=None,**kw):
    if ylower is not None and yupper is not None:
      length = yupper - ylower
      ycent = 0.5*(ylower + yupper)
      self.ylower = ylower
      self.yupper = yupper
    assert length is not None,\
      "YCylinder: either length or both ylower and yupper must be specified"
    ZCylinderOut.__init__(self,radius,length,
                               voltage,xcent,ycent,zcent,condid)
    EllipticAssembly.__init__(self,ellipticity,
                              voltage,xcent,ycent,zcent,condid,self.kwlist,
                              self.generatorf,self.generatord,self.generatori)
    YAssembly.__init__(self,
                       voltage,xcent,ycent,zcent,condid,self.kwlist,
                       self.generatorf,self.generatord,self.generatori,
                       kw=kw)

#============================================================================
class Sphere(Assembly):
  """
Sphere
  - radius: radius
  - voltage=0: sphere voltage
  - xcent=0.,ycent=0.,zcent=0.: center of sphere
  - condid=1: conductor id of sphere, must be integer, or can be 'next' in
              which case a unique ID is chosen
  """
  def __init__(self,radius,voltage=0.,xcent=0.,ycent=0.,zcent=0.,
                    condid=1,**kw):
    kwlist = ['radius']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      sphereconductorf,sphereconductord,sphereintercept,
                      kw=kw)
    self.radius = radius
    self.createextent([-self.radius,-self.radius,-self.radius],
                      [+self.radius,+self.radius,+self.radius])

  def createdxobject(self,kwdict={},**kw):
    kw.update(kwdict)
    v = Opyndx.VisualRevolution(
                       zzmin=-self.radius,zzmax=+self.radius,
                       rendzmin=0.,rendzmax=0.,
                       xoff=self.xcent,yoff=self.ycent,zoff=self.zcent,
                       rofzdata=[0.,0.],
                       zdata=[-self.radius,self.radius],
                       raddata=[self.radius],zcdata=[0.],rcdata=[0.],
                       largepos=largepos,
                       kwdict=kw)
    self.dxobject = v

  def draw(self,color='fg',filled=None,fullplane=1,**kw):
    narcpoints = kw.get('narcpoints',64)
    theta = span(0.,2.*pi,narcpoints+1)
    r = self.radius*cos(theta)
    z = self.radius*sin(theta)
    self.plotdata(r,z,color=color,filled=filled,fullplane=fullplane)

#============================================================================
class ZElliptoid(Sphere,EllipticAssembly):
  """
Elliptoidal sphere
  - ellipticity: ratio of y radius to x radius
  - radius: radius
  - voltage=0: sphere voltage
  - xcent=0.,ycent=0.,zcent=0.: center of sphere
  - condid=1: conductor id of sphere, must be integer, or can be 'next' in
              which case a unique ID is chosen
  """
  def __init__(self,ellipticity,radius,
                    voltage=0.,xcent=0.,ycent=0.,zcent=0.,condid=1,**kw):
    Sphere.__init__(self,radius,
                         voltage,xcent,ycent,zcent,condid)
    EllipticAssembly.__init__(self,ellipticity,
                              voltage,xcent,ycent,zcent,condid,self.kwlist,
                              self.generatorf,self.generatord,self.generatori,
                              kw=kw)

#============================================================================
class Cone(Assembly):
  """
Cone
  - r_zmin: radius at z min
  - r_zmax: radius at z max
  - length: length
  - theta=0,phi=0: angle of cylinder relative to z-axis
    theta is angle in z-x plane
    phi is angle in z-y plane
  - voltage=0: cone voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cone
  - condid=1: conductor id of cone, must be integer, or can be 'next' in
              which case a unique ID is chosen
  """
  def __init__(self,r_zmin,r_zmax,length,theta,phi,voltage=0.,
                    xcent=0.,ycent=0.,zcent=0.,condid=1,**kw):
    kwlist = ['r_zmin','r_zmax','length','theta','phi']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      coneconductorf,coneconductord,coneintercept,
                      kw=kw)
    self.r_zmin = r_zmin
    self.r_zmax = r_zmax
    self.theta = theta
    self.phi = phi
    self.length = length

    rmax = max(sqrt(self.r_zmin**2+(self.length/2.)**2),
               sqrt(self.r_zmax**2+(self.length/2.)**2))
    self.createextent([-rmax,-rmax,-self.length/2.],
                      [+rmax,+rmax,+self.length/2.])

#============================================================================
class ConeSlope(Assembly):
  """
Cone
  - slope: ratio of radius at zmax minus radius at zmin over length
  - intercept: location where line defining cone crosses the axis, relative
               to zcent
  - r_zmax: radius at z max
  - length: length
  - theta=0,phi=0: angle of cylinder relative to z-axis
    theta is angle in z-x plane
    phi is angle in z-y plane
  - voltage=0: cone voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cone
  - condid=1: conductor id of cone, must be integer, or can be 'next' in
              which case a unique ID is chosen
  """
  def __init__(self,slope,length,theta,phi,voltage=0.,
                    xcent=0.,ycent=0.,zcent=0.,condid=1,**kw):
    kwlist = ['r_zmin','r_zmax','length','theta','phi']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      coneconductorf,coneconductord,coneintercept,
                      kw=kw)
    self.slope = slope
    self.intercept = intercept
    self.theta = theta
    self.phi = phi
    self.length = length

    self.r_zmin = self.slope*(-self.length/2. - self.intercept)
    self.r_zmax = self.slope*(+self.length/2. - self.intercept)
    rmax = max(sqrt(self.r_zmin**2+(self.length/2.)**2),
               sqrt(self.r_zmax**2+(self.length/2.)**2))
    self.createextent([-rmax,-rmax,-self.length/2.],
                      [+rmax,+rmax,+self.length/2.])

  def getkwlist(self):
    self.r_zmin = self.slope*(-self.length/2. - self.intercept)
    self.r_zmax = self.slope*(+self.length/2. - self.intercept)
    return Assembly.getkwlist(self)

#============================================================================
class Cones(Assembly):
  """
Cones
  - r_zmin: radius at z min
  - r_zmax: radius at z max
  - length: length
  - theta=0,phi=0: angle of cylinder relative to z-axis
    theta is angle in z-x plane
    phi is angle in z-y plane
  - voltage=0: cone voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cone
  - condid=1: conductor id of cone, must be integer, or can be 'next' in
              which case a unique ID is chosen
  """
  def __init__(self,r_zmin,r_zmax,length,theta,phi,voltage=0.,
                    xcent=0.,ycent=0.,zcent=0.,condid=1,**kw):
    kwlist = ['ncones','r_zmin','r_zmax','length','theta','phi']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      conesconductorf,conesconductord,conesintercept,
                      kw=kw)
    self.ncones = 0
    self.r_zmin = r_zmin
    self.r_zmax = r_zmax
    self.length = length
    self.theta = theta
    self.phi = phi
    # --- Find the first input argument that has a length > 1. This is done
    # --- since some of the inputs may be arrays of length 1, which is OK,
    # --- but does not mean that there is only 1 cone. Some other input
    # --- may be a longer array. getkwlist is used so that the xcent etc.
    # --- are included.
    kwlist = self.getkwlist()
    for k in kwlist:
      try:
        self.ncones = len(k)
        if self.ncones > 1: break
      except TypeError:
        pass

    assert self.ncones > 0,"At least on of the input arguments must be a list!"
    self.r_zmin = self.r_zmin*ones(self.ncones)
    self.r_zmax = self.r_zmax*ones(self.ncones)
    self.length = self.length*ones(self.ncones)
    self.theta  = self.theta*ones(self.ncones)
    self.phi    = self.phi*ones(self.ncones)
    self.xcent  = self.xcent*ones(self.ncones)
    self.ycent  = self.ycent*ones(self.ncones)
    self.zcent  = self.zcent*ones(self.ncones)

    rmax = sqrt(maximum(abs(self.r_zmin),abs(self.r_zmax))**2 +
                        (self.length/2.)**2)
    xmin = min(self.xcent-rmax)
    ymin = min(self.ycent-rmax)
    zmin = min(self.zcent-rmax)
    xmax = max(self.xcent+rmax)
    ymax = max(self.ycent+rmax)
    zmax = max(self.zcent+rmax)
    self.extent = ConductorExtent([xmin,ymin,zmin],[xmax,ymax,zmax])

#============================================================================
class ZTorus(Assembly):
  """
Torus
  - r1: toroidal radius
  - r2: poloidal radius
  - voltage=0: torus voltage
  - xcent=0.,ycent=0.,zcent=0.: center of torus
  - condid=1: conductor id of torus, must be integer, or can be 'next' in
              which case a unique ID is chosen
  """
  def __init__(self,r1,r2,voltage=0.,xcent=0.,ycent=0.,zcent=0.,condid=1,**kw):
    kwlist = ['r1','r2']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      ztorusconductorf,ztorusconductord,ztorusintercept,
                      kw=kw)
    self.r1 = r1
    self.r2 = r2

    rmax = self.r1 + self.r2
    self.createextent([-rmax,-rmax,-self.r2],
                      [+rmax,+rmax,+self.r2])

  def createdxobject(self,kwdict={},**kw):
    kw.update(kwdict)
    v = Opyndx.VisualRevolution(
                       zzmin=-self.r2,zzmax=+self.r2,
                       rendzmin=self.r1,rendzmax=self.r1,
                       xoff=self.xcent,yoff=self.ycent,zoff=self.zcent,
                       rofzdata=[self.r1,self.r1,self.r1],
                       zdata=[-self.r2,self.r2,-self.r2],
                       raddata=[self.r2,self.r2],
                       zcdata=[0.,0.],rcdata=[self.r1,self.r1],
                       largepos=largepos,
                       kwdict=kw)
    self.dxobject = v

  def draw(self,color='fg',filled=None,fullplane=1,**kw):
    narcpoints = kw.get('narcpoints',64)
    theta = span(0.,2.*pi,narcpoints+1)
    r = self.r2*cos(theta) + self.r1
    z = self.r2*sin(theta)
    self.plotdata(r,z,color=color,filled=filled,fullplane=fullplane)

#============================================================================
class Beamletplate(Assembly):
  """
Plate from beamlet pre-accelerator
  - za: location of spherical center in the x-plane
  - zb: location of spherical center in the y-plane
  - z0: location of the center of the plate on the z-axis
  - thickness: thickness of the plate
  - voltage=0: beamlet plate voltage
  - xcent=0.,ycent=0.,zcent=0.: center of beamlet plate
  - condid=1: conductor id of beamlet plate, must be integer, or can be 'next'
              in which case a unique ID is chosen
  """
  def __init__(self,za,zb,z0,thickness,voltage=0.,
               xcent=0.,ycent=0.,zcent=0.,condid=1,**kw):
    kwlist = ['za','zb','z0','thickness']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      beamletplateconductorf,beamletplateconductord,
                      beamletplateintercept,
                      kw=kw)
    self.za = za
    self.zb = zb
    self.z0 = z0
    self.thickness = thickness

    # --- Give a somewhat thoughtful result.
    if za < 1.e10: zza = za - sqrt((za-z0)**2 - w3d.xmmax**2 - w3d.ymmax**2)
    else:          zza = z0
    if zb < 1.e10: zzb = zb - sqrt((zb-z0)**2 - w3d.xmmax**2 - w3d.ymmax**2)
    else:          zzb = z0
    zmin = z0 - thickness
    zmax = max(zza,zzb) + 5*thickness

    self.createextent([-largepos,-largepos,zmin],
                      [+largepos,+largepos,zmax])

  def createdxobject(self,xmin=None,xmax=None,ymin=None,ymax=None,
                nx=None,ny=None,nz=None,
                xmmin=None,xmmax=None,ymmin=None,ymmax=None,
                zmmin=None,zmmax=None,l2symtry=None,l4symtry=None):
    _default = lambda x,d: (x,d)[x is None]
    xmmin = _default(xmmin,w3d.xmmin)
    xmmax = _default(xmmax,w3d.xmmax)
    ymmin = _default(ymmin,w3d.ymmin)
    ymmax = _default(ymmax,w3d.ymmax)
    zmmin = _default(zmmin,w3d.zmmin)
    zmmax = _default(zmmax,w3d.zmmax)
    nx = _default(nx,w3d.nx)
    ny = _default(ny,w3d.ny)
    nz = _default(nz,w3d.nz)
    l2symtry = _default(l2symtry,w3d.l2symtry)
    l4symtry = _default(l4symtry,w3d.l4symtry)

    xmin = _default(xmin,xmmin)
    xmax = _default(xmax,xmmax)
    ymin = _default(ymin,ymmin)
    ymax = _default(ymax,ymmax)
    nx = min(nx,20)
    ny = min(ny,20)

    # --- Check for symmetries
    if l2symtry:
      ymin = 0.
      ymmin = 0.
    elif l4symtry:
      xmin = 0.
      xmmin = 0.
      ymin = 0.
      ymmin = 0.

    # --- Calculate dx, dy, and dz in case this is called before
    # --- the generate.
    dx = (xmmax - xmmin)/nx
    dy = (ymmax - ymmin)/ny
    dz = (zmmax - zmmin)/nz
    if ny == 0: dy = dx

    xmesh = xmmin + dx*arange(nx+1)
    ymesh = ymmin + dy*arange(ny+1)
    xmesh = compress(logical_and(xmin-dx <= xmesh,xmesh <= xmax+dx),xmesh)
    ymesh = compress(logical_and(ymin-dy <= ymesh,ymesh <= ymax+dy),ymesh)
    x = ravel(xmesh[:,NewAxis]*ones(len(ymesh)))
    y = ravel(ymesh*ones(len(xmesh))[:,NewAxis])
    ix = nint((x - xmmin)/dx)
    iy = nint((y - ymmin)/dy)
    if len(x) == 0: return

    xx = x
    yy = y
    xx.shape = (len(xmesh),len(ymesh))
    yy.shape = (len(xmesh),len(ymesh))

    # --- Outer face
    z = self.z0*ones(len(xmesh)*len(ymesh),'d') - self.thickness
    iz = nint((z - zmmin)/dz)
    d = self.griddistance(ix,iy,iz,x,y,z)
    zl = z[0] + d.dels[5,:]
    zl.shape = (len(xmesh),len(ymesh))
    ml = Opyndx.VisualMesh(xx,yy,zl,twoSided=true)

    # --- Inner face
    z = self.z0*ones(len(xmesh)*len(ymesh),'d') + 0.5*self.za
    iz = nint((z - zmmin)/dz)
    d = self.griddistance(ix,iy,iz,x,y,z)
    zr = z[0] - d.dels[4,:]
    zr.shape = (len(xmesh),len(ymesh))
    mr = Opyndx.VisualMesh(xx,yy,zr,twoSided=true)

    # --- Four sides between faces
    xside = xx[:,0]*ones(2)[:,NewAxis]
    yside = yy[:,0]*ones(2)[:,NewAxis]
    zside = array([zl[:,0],zr[:,0]])
    ms1 = Opyndx.VisualMesh(xside,yside,zside,twoSided=true)

    xside = xx[:,-1]*ones(2)[:,NewAxis]
    yside = yy[:,-1]*ones(2)[:,NewAxis]
    zside = array([zl[:,-1],zr[:,-1]])
    ms1 = Opyndx.VisualMesh(xside,yside,zside,twoSided=true)

    xside = xx[0,:]*ones(2)[:,NewAxis]
    yside = yy[0,:]*ones(2)[:,NewAxis]
    zside = array([zl[0,:],zr[0,:]])
    ms1 = Opyndx.VisualMesh(xside,yside,zside,twoSided=true)

    xside = xx[-1,:]*ones(2)[:,NewAxis]
    yside = yy[-1,:]*ones(2)[:,NewAxis]
    zside = array([zl[-1,:],zr[-1,:]])
    ms1 = Opyndx.VisualMesh(xside,yside,zside,twoSided=true)

#============================================================================
#============================================================================
class Srfrv:
  """
Generic surface of revolution. Contains utility routines to check the
data and make sure it is consistent.
  """
  def checkarcs(self,zz,rr,rad,zc,rc):
    """Utility routine for surface of revoluation routines.
  Checks consistency of input and calculates any parameters not given.
    """
    for i in range(len(zz)-1):
      if ((rad[i] == None or rad[i] == largepos) and
          (zc[i] == None or zc[i] == largepos) and
          (rc[i] == None or rc[i] == largepos)):
        # --- When there is a straight line, then set the radius to a large
        # --- number (used as a flag in the code).
        rad[i] = largepos
        zc[i] = largepos
        rc[i] = largepos
      elif ((zc[i] == None or zc[i] == largepos) or
            (rc[i] == None or rc[i] == largepos)):
        # --- Given a radius and the two endpoints, the center of the
        # --- circle can be found.
        assert 4*rad[i]**2 > ((zz[i] - zz[i+1])**2 + (rr[i] - rr[i+1])**2),\
             "Radius of circle must be larger than the distance between points %d,%d\n%e < %e"%\
             (i,i+1,2*rad[i],sqrt((zz[i] - zz[i+1])**2 + (rr[i] - rr[i+1])**2))
        zm = 0.5*(zz[i] + zz[i+1])
        rm = 0.5*(rr[i] + rr[i+1])
        dbm = sqrt((zm - zz[i+1])**2 + (rm - rr[i+1])**2)
        dcm = sqrt(rad[i]**2 - dbm**2)
        angle1 = arcsin((rm - rr[i+1])/dbm)
        if rad[i] < 0:
          zc[i] = zm + dcm*sin(angle1)
          rc[i] = rm + dcm*cos(angle1)
        else:
          zc[i] = zm - dcm*sin(angle1)
          rc[i] = rm - dcm*cos(angle1)
      elif (rad[i] == None or rad[i] == largepos):
        # --- Given the center, the radius can be found. With the two end
        # --- points of the arc given, the data is redundant, so check
        # --- to be sure it is consistent.
        rad[i] = sqrt((zz[i] - zc[i])**2 + (rr[i] - rc[i])**2)
        rad2 = sqrt((zz[i+1] - zc[i])**2 + (rr[i+1] - rc[i])**2)
        assert (abs(rad[i] - rad2)/rad[i] < 1.e-2),\
           "Points %d and %d are not at the same radius relative to the arc center. The radii are %e and %e"%(i,i+1,rad[i],rad2)
        # --- Make sure the radius has the correct sign.
        if rc[i] > rr[i] or rc[i] > rr[i+1]: rad[i] = -rad[i]

  def setdatadefaults(self,data,ndata,default):
    if data is None:
      data = ndata*[default]
    else:
      assert len(data) == ndata,\
             "Some of the input surface data is not the correct length"
    return data

  def getplotdata(self,lrfunc,rfunc,
                  npoints,rdata,zdata,raddata,rcdata,zcdata,narcpoints):
    r = []
    z = []
    if lrfunc:
      # --- Get the data from the rofz function
      if npoints is None:
        solver = getregisteredsolver()
        if solver is None: dz = w3d.dz
        else:              dz = solver.dz
        npoints = max(100,nint((self.zmax-self.zmin)/dz))
      dz = (self.zmax - self.zmin)/npoints
      for i in range(npoints+1):
        f3d.srfrv_z = self.zmin + i*dz
        rfunc()
        r.append(f3d.srfrv_r)
        z.append(f3d.srfrv_z)
    else:
      # --- Get the data from the rofz table
      for i in range(len(rdata)-1):
        if raddata[i] == largepos:
          r.append(rdata[i])
          z.append(zdata[i])
        else:
          zz = span(zdata[i],zdata[i+1],narcpoints)
          if raddata[i] > 0.:
            rr = rcdata[i] + sqrt(maximum(0,raddata[i]**2 - (zz-zcdata[i])**2))
          else:
            rr = rcdata[i] - sqrt(maximum(0,raddata[i]**2 - (zz-zcdata[i])**2))
          r = r + list(rr)
          z = z + list(zz)
      r.append(rdata[-1])
      z.append(zdata[-1])

      # --- zcent needs to be added in when chopping the data at zmin and zmax,
      # --- but then subtracted out since the base plotting routines adds
      # --- in zcent.
      z = array(z) + self.zcent
      z = minimum(self.zmax,maximum(self.zmin,array(z)))
      z = list(z - self.zcent)

    return r,z

#============================================================================
# --- These handle the callbacks when a python function is used to
# --- describe the surface of revolution.
def rofzfunc():
  try:
    f3d.srfrv_r = rofzfunc.rofzfunc(f3d.srfrv_z)
  except TypeError:
    rofzfunc.rofzfunc()
def rminofz():
  try:
    f3d.srfrv_r = rminofz.rminofz(f3d.srfrv_z)
  except TypeError:
    rminofz.rminofz()
def rmaxofz():
  try:
    f3d.srfrv_r = rmaxofz.rmaxofz(f3d.srfrv_z)
  except TypeError:
    rmaxofz.rmaxofz()

#============================================================================
class ZSrfrvOut(Srfrv,Assembly):
  """
Outside of a surface of revolution
  - rofzfunc=None: python function describing surface
  - zmin=None,zmax=None: z-extent of the surface, will be obtained from
                         any tablized data if not given
  - rmax=largepos: max radius of the surface
  - voltage=0: conductor voltage
  - xcent=0.,ycent=0.,zcent=0.: center of conductor
  - condid=1: conductor id of conductor, must be integer, or can be 'next' in
              which case a unique ID is chosen
  - rofzdata=None: optional tablized data of radius of surface
  - zdata=None: optional tablized data of z locations of rofzdata
      raddata[i] is radius for segment from zdata[i] to zdata[i+1]
  - zcdata=None: z center of circle or curved segment
  - rcdata=None: r center of circle or curved segment
  - raddata=None: optional radius of curvature of segments
      The centers of the circles will be calculated automatically if
      not supplied, or if the centers are supplied, the radii will be
      calculated automatically.
      The length of the radii and centers lists is one less than the length
      of the list of r and z data.
    Note that if tablized data is given, the first argument is ignored.
Methods:
  - draw: draws the object's r versus z
  - createdxobject: creates (internally) the object for visualization using
                    OpenDX. This can be used to specify options on how the
                    image is made before passing the object to DXImage.
  - getdxobject: returns a object for visualization. This can be used to
                 specify options on how the image is made. The returned object
                 is then passed to DXImage
  """
  def __init__(self,rofzfunc=None,zmin=None,zmax=None,rmax=largepos,
                    voltage=0.,xcent=0.,ycent=0.,zcent=0.,condid=1,
                    rofzdata=None,zdata=None,raddata=None,
                    zcdata=None,rcdata=None,**kw):
    kwlist = ['lrofzfunc','zmin','zmax','rmax','griddz']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      zsrfrvoutconductorf,zsrfrvoutconductord,
                      zsrfrvoutintercept,
                      kw=kw)
    self.rofzfunc = rofzfunc
    self.rmax = rmax

    # --- Deal with tablized data.
    # --- Make sure the input is consistent
    if operator.isSequenceType(rofzdata):
      self.lrofzfunc = false
      self.zdata = zdata
      self.rofzdata = self.setdatadefaults(rofzdata,len(zdata),rmax)
      self.raddata = self.setdatadefaults(raddata,len(zdata)-1,None)
      self.zcdata = self.setdatadefaults(zcdata,len(zdata)-1,None)
      self.rcdata = self.setdatadefaults(rcdata,len(zdata)-1,None)
      self.checkarcs(self.zdata,self.rofzdata,self.raddata,
                     self.zcdata,self.rcdata)
      if zmin is None: zmin = self.zdata[0]
      if zmax is None: zmax = self.zdata[-1]
    else:
      assert type(self.rofzfunc) in [FunctionType,StringType],\
             'The rofzfunc is not properly specified'
      self.lrofzfunc = true
      if isinstance(self.rofzfunc,StringType):
        # --- Check if the rofzfunc is in main. Complain if it is not.
        import __main__
        self.rofzfunc = __main__.__dict__[self.rofzfunc]
      self.rofzdata = None
      self.zdata = None
      self.raddata = None
      self.rcdata = None
      self.zcdata = None

    assert zmin is not None,'zmin must be specified'
    assert zmax is not None,'zmin must be specified'
    self.zmin = zmin
    self.zmax = zmax

    self.createextent([-self.rmax,-self.rmax,self.zmin],
                      [+self.rmax,+self.rmax,self.zmax])

  def getkwlist(self):
    self.griddz = _griddzkludge[0]
    # --- If data arrays are specified, then put the data in the right place
    if self.lrofzfunc:
      rofzfunc.rofzfunc = self.rofzfunc
    else:
      f3d.npnts_sr = len(self.zdata)
      f3d.z_sr = self.zdata
      f3d.r_sr = self.rofzdata
      f3d.rad_sr = self.raddata
      f3d.zc_sr = self.zcdata
      f3d.rc_sr = self.rcdata

    return Assembly.getkwlist(self)

  def draw(self,color='fg',filled=None,fullplane=1,nzpoints=None,**kw):
    """
Plots the r versus z
 - color='fg': color of outline, set to None to not plot the outline
 - filled=None: when set to an integer, fills the outline with the color
                specified from the current palette. Should be between 0 and 199.
 - fullplane=1: when true, plot the top and bottom, i.e. r vs z, and -r vs z.
 - rmax=None: when given, overrides the instance's value of rmax, useful in
              cases when the instance's value of rmax is largepos.
 - narcpoints=40: number of points to draw along any circular arcs
    """
    narcpoints = kw.get('narcpoints',40)
    rmax = kw.get('rmax',None)
    rofzfunc.rofzfunc = self.rofzfunc
    r,z = self.getplotdata(self.lrofzfunc,rofzfunc,
                           nzpoints,
                           self.rofzdata,self.zdata,self.raddata,
                           self.rcdata,self.zcdata,narcpoints)
    if rmax is None: rmax = self.rmax
    r = [rmax] + r + [rmax]
    z = [self.zmin] + z + [self.zmax]
    self.plotdata(r,z,color=color,filled=filled,fullplane=fullplane)

  def createdxobject(self,rmax=None,kwdict={},**kw):
    """
Creates internally the object to be used for visualization.
 - rmax=None: when given, overrides the instance's value of rmax, useful in
              cases when the instance's value of rmax is largepos.
For other options, see documentation of Opyndx.VisualRevolution.
    """
    kw.update(kwdict)
    if rmax is None: rmax = self.rmax
    v = Opyndx.VisualRevolution(self.rofzfunc,self.zmin,self.zmax,
                       rendzmin=rmax,rendzmax=rmax,
                       xoff=self.xcent,yoff=self.ycent,zoff=self.zcent,
                       rofzdata=self.rofzdata,zdata=self.zdata,
                       raddata=self.raddata,zcdata=self.zcdata,
                       rcdata=self.rcdata,
                       normalsign=-1,
                       largepos=largepos,
                       kwdict=kw)
    self.dxobject = v

#============================================================================
class ZSrfrvIn(Srfrv,Assembly):
  """
Inside of a surface of revolution
  - rofzfunc=None: python function describing surface
  - zmin=None,zmax=None: z-extent of the surface, will be obtained from
                         any tablized data if not given
  - rmin=0: min radius of the surface
  - voltage=0: conductor voltage
  - xcent=0.,ycent=0.,zcent=0.: center of conductor
  - condid=1: conductor id of conductor, must be integer, or can be 'next' in
              which case a unique ID is chosen
  - rofzdata=None: optional tablized data of radius of surface
  - zdata=None: optional tablized data of z locations of rofzdata
  - raddata=None: optional radius of curvature of segments
      raddata[i] is radius for segment from zdata[i] to zdata[i+1]
  - zcdata=None: z center of circle or curved segment
  - rcdata=None: r center of circle or curved segment
      The centers of the circles will be calculated automatically if
      not supplied, or if the centers are supplied, the radii will be
      calculated automatically.
      The length of the radii and centers lists is one less than the length
      of the list of r and z data.
    Note that if tablized data is given, the first argument is ignored.
Methods:
  - draw: draws the object's r versus z
  - createdxobject: creates (internally) the object for visualization using
                    OpenDX. This can be used to specify options on how the
                    image is made before passing the object to DXImage.
  - getdxobject: returns a object for visualization. This can be used to
                 specify options on how the image is made. The returned object
                 is then passed to DXImage
  """
  def __init__(self,rofzfunc=None,zmin=None,zmax=None,rmin=0,
                    voltage=0.,xcent=0.,ycent=0.,zcent=0.,condid=1,
                    rofzdata=None,zdata=None,raddata=None,
                    zcdata=None,rcdata=None,**kw):
    kwlist = ['lrofzfunc','zmin','zmax','rmin','griddz']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      zsrfrvinconductorf,zsrfrvinconductord,
                      zsrfrvinintercept,
                      kw=kw)
    self.rofzfunc = rofzfunc
    self.rmin = rmin

    # --- Deal with tablized data.
    # --- Make sure the input is consistent
    if operator.isSequenceType(rofzdata):
      self.lrofzfunc = false
      self.zdata = zdata
      self.rofzdata = self.setdatadefaults(rofzdata,len(zdata),rmin)
      self.raddata = self.setdatadefaults(raddata,len(zdata)-1,None)
      self.zcdata = self.setdatadefaults(zcdata,len(zdata)-1,None)
      self.rcdata = self.setdatadefaults(rcdata,len(zdata)-1,None)
      self.checkarcs(self.zdata,self.rofzdata,self.raddata,
                     self.zcdata,self.rcdata)
      if zmin is None: zmin = self.zdata[0]
      if zmax is None: zmax = self.zdata[-1]
    else:
      assert type(self.rofzfunc) in [FunctionType,StringType],\
             'The rofzfunc is not properly specified'
      self.lrofzfunc = true
      if isinstance(self.rofzfunc,StringType):
        # --- Check if the rofzfunc is in main. Complain if it is not.
        import __main__
        self.rofzfunc = __main__.__dict__[self.rofzfunc]
      self.rofzdata = None
      self.zdata = None
      self.raddata = None
      self.rcdata = None
      self.zcdata = None

    assert zmin is not None,'zmin must be specified'
    assert zmax is not None,'zmin must be specified'
    self.zmin = zmin
    self.zmax = zmax

    if self.lrofzfunc: rmax = largepos
    else:              rmax = max(self.rofzdata)
    self.createextent([-rmax,-rmax,self.zmin],
                      [+rmax,+rmax,self.zmax])

  def getkwlist(self):
    self.griddz = _griddzkludge[0]
    # --- If data arrays are specified, then put the data in the right place
    if self.lrofzfunc:
      rofzfunc.rofzfunc = self.rofzfunc
    else:
      f3d.npnts_sr = len(self.zdata)
      f3d.z_sr = self.zdata
      f3d.r_sr = self.rofzdata
      f3d.rad_sr = self.raddata
      f3d.zc_sr = self.zcdata
      f3d.rc_sr = self.rcdata

    return Assembly.getkwlist(self)

  def draw(self,color='fg',filled=None,fullplane=1,nzpoints=None,**kw):
    """
Plots the r versus z
 - color='fg': color of outline, set to None to not plot the outline
 - filled=None: when set to an integer, fills the outline with the color
                specified from the current palette. Should be between 0 and 199.
 - fullplane=1: when true, plot the top and bottom, i.e. r vs z, and -r vs z.
 - rmin: when given, overrides the instance's values of rmin - this is not
         really useful
 - narcpoints=40: number of points to draw along any circular arcs
    """
    narcpoints = kw.get('narcpoints',40)
    rmin = kw.get('rmin',None)
    rofzfunc.rofzfunc = self.rofzfunc
    r,z = self.getplotdata(self.lrofzfunc,rofzfunc,
                           nzpoints,
                           self.rofzdata,self.zdata,self.raddata,
                           self.rcdata,self.zcdata,narcpoints)
    if rmin is None: rmin = self.rmin
    r = [rmin] + r + [rmin]
    z = [self.zmin] + z + [self.zmax]
    self.plotdata(r,z,color=color,filled=filled,fullplane=fullplane)

  def createdxobject(self,kwdict={},**kw):
    """
Creates internally the object to be used for visualization.
For options, see documentation of Opyndx.VisualRevolution.
    """
    kw.update(kwdict)
    v = Opyndx.VisualRevolution(self.rofzfunc,self.zmin,self.zmax,
                       rendzmin=self.rmin,rendzmax=self.rmin,
                       xoff=self.xcent,yoff=self.ycent,zoff=self.zcent,
                       rofzdata=self.rofzdata,zdata=self.zdata,
                       raddata=self.raddata,zcdata=self.zcdata,
                       rcdata=self.rcdata,
                       largepos=largepos,
                       kwdict=kw)
    self.dxobject = v

#============================================================================
class ZSrfrvInOut(Srfrv,Assembly):
  """
Between surfaces of revolution
  - rminofz=None,rmaxofz=None: python functions describing surfaces
  - zmin=None,zmax=None: z-extent of the surface, will be obtained from
                         any tablized if not given
  - voltage=0: conductor voltage
  - xcent=0.,ycent=0.,zcent=0.: center of conductor
  - condid=1: conductor id of conductor, must be integer, or can be 'next' in
              which case a unique ID is chosen
  - rminofzdata,rmaxofzdata=None: optional tablized data of radii of surface
  - zmindata,zmaxdata=None: optional tablized data of z locations of r data
  - radmindata,radmaxdata=None: optional radius of curvature of segments
      radmindata[i] is radius for segment from zmindata[i] to zmindata[i+1]
  - zcmindata,zcmaxdata=None: z center of circle or curved segment
  - rcmindata,rcmaxdata=None: r center of circle or curved segment
      The centers of the circles will be calculated automatically if
      not supplied, or if the centers are supplied, the radii will be
      calculated automatically.
      The length of the radii and centers lists is one less than the length
      of the list of r and z data.
    Note that if tablized data is given, the first two arguments are ignored.
Methods:
  - draw: draws the object's r versus z
  - createdxobject: creates (internally) the object for visualization using
                    OpenDX. This can be used to specify options on how the
                    image is made before passing the object to DXImage.
  - getdxobject: returns a object for visualization. This can be used to
                 specify options on how the image is made. The returned object
                 is then passed to DXImage
  """
  def __init__(self,rminofz=None,rmaxofz=None,zmin=None,zmax=None,
                    voltage=0.,xcent=0.,ycent=0.,zcent=0.,condid=1,
                    rminofzdata=None,zmindata=None,radmindata=None,
                    rcmindata=None,zcmindata=None,
                    rmaxofzdata=None,zmaxdata=None,radmaxdata=None,
                    rcmaxdata=None,zcmaxdata=None,**kw):
    kwlist = ['lrminofz','lrmaxofz','zmin','zmax','griddz']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      zsrfrvinoutconductorf,zsrfrvinoutconductord,
                      zsrfrvinoutintercept,
                      kw=kw)
    self.rminofz = rminofz
    self.rmaxofz = rmaxofz

    # --- Deal with tablized data.
    # --- Making sure the input is consistent
    if operator.isSequenceType(zmindata):
      self.lrminofz = false
      self.zmindata = zmindata
      self.rminofzdata = self.setdatadefaults(rminofzdata,len(zmindata),0.)
      self.radmindata = self.setdatadefaults(radmindata,len(zmindata)-1,None)
      self.rcmindata = self.setdatadefaults(rcmindata,len(zmindata)-1,None)
      self.zcmindata = self.setdatadefaults(zcmindata,len(zmindata)-1,None)
      self.checkarcs(self.zmindata,self.rminofzdata,self.radmindata,
                     self.zcmindata,self.rcmindata)
      zminmin = self.zmindata[0]
      zmaxmin = self.zmindata[-1]
    else:
      assert type(self.rminofz) in [FunctionType,StringType],\
             'The rminofz is not properly specified'
      self.lrminofz = true
      zminmin = zmin
      zmaxmin = zmax
      if isinstance(self.rminofz,StringType):
        # --- Check if the rofzfunc is in main. Complain if it is not.
        import __main__
        self.rminofz = __main__.__dict__[self.rminofz]
      self.rminofzdata = None
      self.zmindata = None
      self.radmindata = None
      self.rcmindata = None
      self.zcmindata = None

    if operator.isSequenceType(zmaxdata):
      self.lrmaxofz = false
      self.zmaxdata = zmaxdata
      self.rmaxofzdata = self.setdatadefaults(rmaxofzdata,len(zmaxdata),
                                              largepos)
      self.radmaxdata = self.setdatadefaults(radmaxdata,len(zmaxdata)-1,None)
      self.rcmaxdata = self.setdatadefaults(rcmaxdata,len(zmaxdata)-1,None)
      self.zcmaxdata = self.setdatadefaults(zcmaxdata,len(zmaxdata)-1,None)
      self.checkarcs(self.zmaxdata,self.rmaxofzdata,self.radmaxdata,
                     self.zcmaxdata,self.rcmaxdata)
      zminmax = self.zmaxdata[0]
      zmaxmax = self.zmaxdata[-1]
    else:
      assert type(self.rmaxofz) in [FunctionType,StringType],\
             'The rmaxofz is not properly specified'
      self.lrmaxofz = true
      zminmax = zmin
      zmaxmax = zmax
      if isinstance(self.rmaxofz,StringType):
        # --- Check if the rofzfunc is in main. Complain if it is not.
        import __main__
        self.rmaxofz = __main__.__dict__[self.rmaxofz]
      self.rmaxofzdata = None
      self.zmaxdata = None
      self.radmaxdata = None
      self.rcmaxdata = None
      self.zcmaxdata = None

    # --- If zmin or zmax were not specified, get the extremum from any tablized
    # --- data.
    if zmin is None:
      if zminmin is not None and zminmax is not None: zmin = min(zminmin,zminmax)
      elif zminmin is not None: zmin = zminmin
      elif zminmax is not None: zmin = zminmax
      else: raise 'zmin must be specified'
    if zmax is None:
      if zmaxmin is not None and zmaxmax is not None: zmax = min(zmaxmin,zmaxmax)
      elif zmaxmin is not None: zmax = zmaxmin
      elif zmaxmax is not None: zmax = zmaxmax
      else: raise 'zmax must be specified'
    self.zmin = zmin
    self.zmax = zmax

    if self.lrmaxofz: rmax = largepos
    else:             rmax = max(self.rmaxofzdata)
    self.createextent([-rmax,-rmax,self.zmin],
                      [+rmax,+rmax,self.zmax])

  def getkwlist(self):
    self.griddz = _griddzkludge[0]
    # --- If data arrays are specified, then put the data in the right place
    if self.lrminofz:
      rminofz.rminofz = self.rminofz
    else:
      f3d.npnts_srmin = len(self.zmindata)
      f3d.z_srmin = self.zmindata
      f3d.r_srmin = self.rminofzdata
      f3d.rad_srmin = self.radmindata
      f3d.zc_srmin = self.zcmindata
      f3d.rc_srmin = self.rcmindata

    if self.lrmaxofz:
      rmaxofz.rmaxofz = self.rmaxofz
    else:
      f3d.npnts_srmax = len(self.zmaxdata)
      f3d.z_srmax = self.zmaxdata
      f3d.r_srmax = self.rmaxofzdata
      f3d.rad_srmax = self.radmaxdata
      f3d.zc_srmax = self.zcmaxdata
      f3d.rc_srmax = self.rcmaxdata

    return Assembly.getkwlist(self)

  def draw(self,color='fg',filled=None,fullplane=1,nzpoints=None,**kw):
    """
Plots the r versus z
 - color='fg': color of outline, set to None to not plot the outline
 - filled=None: when set to an integer, fills the outline with the color
                specified from the current palette. Should be between 0 and 199.
 - fullplane=1: when true, plot the top and bottom, i.e. r vs z, and -r vs z.
 - narcpoints=40: number of points to draw along any circular arcs
    """
    narcpoints = kw.get('narcpoints',40)
    rminofz.rminofz = self.rminofz
    ri,zi = self.getplotdata(self.lrminofz,rminofz,
                             nzpoints,
                             self.rminofzdata,self.zmindata,self.radmindata,
                             self.rcmindata,self.zcmindata,narcpoints)
    rmaxofz.rmaxofz = self.rmaxofz
    ro,zo = self.getplotdata(self.lrmaxofz,rmaxofz,
                             nzpoints,
                             self.rmaxofzdata,self.zmaxdata,self.radmaxdata,
                             self.rcmaxdata,self.zcmaxdata,narcpoints)
    ro.reverse()
    zo.reverse()
    r,z = ri+ro,zi+zo
    r = r + [r[0]]
    z = z + [z[0]]
    if len(r) > 0:
      self.plotdata(r,z,color=color,filled=filled,fullplane=fullplane)

  def createdxobject(self,kwdict={},**kw):
    """
Creates internally the object to be used for visualization.
For options, see documentation of Opyndx.VisualRevolution.
    """
    kw.update(kwdict)
    if self.lrminofz:
      f3d.srfrv_z = self.zmin
      self.rminofz()
      rminzmin = f3d.srfrv_r
      f3d.srfrv_z = self.zmax
      self.rminofz()
      rminzmax = f3d.srfrv_r
    else:
      rminzmin = self.rminofzdata[0]
      rminzmax = self.rminofzdata[-1]
    if self.lrmaxofz:
      f3d.srfrv_z = self.zmin
      self.rmaxofz()
      rmaxzmin = f3d.srfrv_r
      f3d.srfrv_z = self.zmax
      self.rmaxofz()
      rmaxzmax = f3d.srfrv_r
    else:
      rmaxzmin = self.rmaxofzdata[0]
      rmaxzmax = self.rmaxofzdata[-1]
    rendzmin = 0.5*(rminzmin + rmaxzmin)
    rendzmax = 0.5*(rminzmax + rmaxzmax)

   #if not self.lrminofz or not self.lrmaxofz:

   #  --- This doesn't quite work and I didn't want to put the effort
   #  --- in to fix it.
   #  rr = concatenate((self.rmaxofzdata,array(self.rminofzdata)[::-1]))
   #  zz = concatenate((self.zmaxdata,array(self.zmindata)[::-1]))
   #  radmin = array(self.radmindata)
   #  radmin = where(radmin < largepos,-radmin,largepos)
   #  rad = concatenate((self.radmaxdata,[largepos],array(radmin)[::-1]))
   #  zc = concatenate((self.zcmaxdata,[0.],array(self.zcmindata)[::-1]))
   #  rc = concatenate((self.rcmaxdata,[0.],array(self.rcmindata)[::-1]))
   #  v = Opyndx.VisualRevolution(' ',self.zmin,self.zmax,
   #                   rendzmin=rendzmin,rendzmax=rendzmax,
   #                   xoff=self.xcent,yoff=self.ycent,zoff=self.zcent,
   #                   rofzdata=rr,zdata=zz,raddata=rad,zcdata=zc,rcdata=rc,
   #                   kwdict=kw)
   #else:
    if 1:
      vmin = Opyndx.VisualRevolution(self.rminofz,self.zmin,self.zmax,
                       rendzmin=rendzmin,rendzmax=rendzmax,
                       xoff=self.xcent,yoff=self.ycent,zoff=self.zcent,
                       rofzdata=self.rminofzdata,zdata=self.zmindata,
                       raddata=self.radmindata,zcdata=self.zcmindata,
                       rcdata=self.rcmindata,
                       normalsign=-1,
                       largepos=largepos,
                       kwdict=kw)
      vmax = Opyndx.VisualRevolution(self.rmaxofz,self.zmin,self.zmax,
                       rendzmin=rendzmin,rendzmax=rendzmax,
                       xoff=self.xcent,yoff=self.ycent,zoff=self.zcent,
                       rofzdata=self.rmaxofzdata,zdata=self.zmaxdata,
                       raddata=self.radmaxdata,zcdata=self.zcmaxdata,
                       rcdata=self.rcmaxdata,
                       largepos=largepos,
                       kwdict=kw)
      v = Opyndx.DXCollection(vmin,vmax)

    self.dxobject = v

#============================================================================
class ZSrfrvEllipticOut(ZSrfrvOut,EllipticAssembly):
  """
Outside of an elliptical surface of revolution
  - ellipticity: ratio of y radius to x radius
  - rofzfunc: name of python function describing surface
  - zmin,zmax: z-extent of the surface
  - rmax=largepos: max radius of the surface
  - voltage=0: conductor voltage
  - xcent=0.,ycent=0.,zcent=0.: center of conductor
  - condid=1: conductor id of conductor, must be integer, or can be 'next' in
              which case a unique ID is chosen
  - rofzdata=None: optional tablized data of radius of surface
  - zdata=None: optional tablized data of z locations of rofzdata
      raddata[i] is radius for segment from zdata[i] to zdata[i+1]
  - zcdata=None: z center of circle or curved segment
  - rcdata=None: r center of circle or curved segment
  - raddata=None: optional radius of curvature of segments
      The centers of the circles will be calculated automatically if
      not supplied, or if the centers are supplied, the radii will be
      calculated automatically.
      The length of the radii and centers lists is one less than the length
      of the list of r and z data.
    Note that if tablized data is given, the first argument is ignored.
  """
  def __init__(self,ellipticity,rofzfunc,zmin,zmax,rmax=largepos,
                    voltage=0.,xcent=0.,ycent=0.,zcent=0.,condid=1,
                    rofzdata=None,zdata=None,raddata=None,
                    zcdata=None,rcdata=None,**kw):
    ZSrfrvOut.__init__(self,rofzfunc,zmin,zmax,rmax,
                       voltage,xcent,ycent,zcent,condid,
                       rofzdata,zdata,raddata,
                       zcdata,rcdata)
    EllipticAssembly.__init__(self,ellipticity,voltage,xcent,ycent,zcent,
                              condid,self.kwlist,
                              self.generatorf,self.generatord,self.generatori,
                              kw=kw)

#============================================================================
class ZSrfrvEllipticIn(ZSrfrvIn,EllipticAssembly):
  """
Inside of an elliptical surface of revolution
  - ellipticity: ratio of y radius to x radius
  - rofzfunc: name of python function describing surface
  - zmin,zmax: z-extent of the surface
  - rmin=0: min radius of the surface
  - voltage=0: conductor voltage
  - xcent=0.,ycent=0.,zcent=0.: center of conductor
  - condid=1: conductor id of conductor, must be integer, or can be 'next' in
              which case a unique ID is chosen
  - rofzdata=None: optional tablized data of radius of surface
  - zdata=None: optional tablized data of z locations of rofzdata
  - raddata=None: optional radius of curvature of segments
      raddata[i] is radius for segment from zdata[i] to zdata[i+1]
  - zcdata=None: z center of circle or curved segment
  - rcdata=None: r center of circle or curved segment
      The centers of the circles will be calculated automatically if
      not supplied, or if the centers are supplied, the radii will be
      calculated automatically.
      The length of the radii and centers lists is one less than the length
      of the list of r and z data.
    Note that if tablized data is given, the first argument is ignored.
  """
  def __init__(self,ellipticity,rofzfunc,zmin,zmax,rmin=0,
                    voltage=0.,xcent=0.,ycent=0.,zcent=0.,condid=1,
                    rofzdata=None,zdata=None,raddata=None,
                    zcdata=None,rcdata=None,**kw):
    ZSrfrvIn.__init__(self,rofzfunc,zmin,zmax,rmin,
                      voltage,xcent,ycent,zcent,condid,
                      rofzdata,zdata,raddata,
                      zcdata,rcdata)
    EllipticAssembly.__init__(self,ellipticity,voltage,xcent,ycent,zcent,
                              condid,self.kwlist,
                              self.generatorf,self.generatord,self.generatori,
                              kw=kw)

#============================================================================
class ZSrfrvEllipticInOut(ZSrfrvInOut,EllipticAssembly):
  """
Between elliptical surfaces of revolution
  - ellipticity: ratio of y radius to x radius
  - rminofz,rmaxofz: names of python functions describing surfaces
  - zmin,zmax: z-extent of the surface
  - voltage=0: conductor voltage
  - xcent=0.,ycent=0.,zcent=0.: center of conductor
  - condid=1: conductor id of conductor, must be integer, or can be 'next' in
              which case a unique ID is chosen
  - rminofzdata,rmaxofzdata=None: optional tablized data of radii of surface
  - zmindata,zmaxdata=None: optional tablized data of z locations of r data
  - radmindata,radmaxdata=None: optional radius of curvature of segments
      radmindata[i] is radius for segment from zmindata[i] to zmindata[i+1]
  - zcmindata,zcmaxdata=None: z center of circle or curved segment
  - rcmindata,rcmaxdata=None: r center of circle or curved segment
      The centers of the circles will be calculated automatically if
      not supplied, or if the centers are supplied, the radii will be
      calculated automatically.
      The length of the radii and centers lists is one less than the length
      of the list of r and z data.
    Note that if tablized data is given, the first two arguments are ignored.
  """
  def __init__(self,ellipticity,rminofz,rmaxofz,zmin,zmax,
                    voltage=0.,xcent=0.,ycent=0.,zcent=0.,condid=1,
                    rminofzdata=None,zmindata=None,radmindata=None,
                    rcmindata=None,zcmindata=None,
                    rmaxofzdata=None,zmaxdata=None,radmaxdata=None,
                    rcmaxdata=None,zcmaxdata=None,**kw):
    ZSrfrvInOut.__init__(self,rminofz,rmaxofz,zmin,zmax,
                         voltage,xcent,ycent,zcent,condid,
                         rminofzdata,zmindata,radmindata,
                         rcmindata,zcmindata,
                         rmaxofzdata,zmaxdata,radmaxdata,
                         rcmaxdata,zcmaxdata)
    EllipticAssembly.__init__(self,ellipticity,voltage,xcent,ycent,zcent,
                              condid,self.kwlist,
                              self.generatorf,self.generatord,self.generatori,
                              kw=kw)

#============================================================================
class XSrfrvOut(ZSrfrvOut,XAssembly):
  """
Outside of an surface of revolution aligned along to X axis
  - rofxfunc: name of python function describing surface
  - xmin,xmax: x-extent of the surface
  - rmax=largepos: max radius of the surface
  - voltage=0: conductor voltage
  - xcent=0.,ycent=0.,zcent=0.: center of conductor
  - condid=1: conductor id of conductor, must be integer, or can be 'next' in
              which case a unique ID is chosen
  - rofxdata=None: optional tablized data of radius of surface
  - xdata=None: optional tablized data of x locations of rofxdata
      raddata[i] is radius for segment from xdata[i] to xdata[i+1]
  - xcdata=None: x center of circle or curved segment
  - rcdata=None: r center of circle or curved segment
  - raddata=None: optional radius of curvature of segments
      The centers of the circles will be calculated automatically if
      not supplied, or if the centers are supplied, the radii will be
      calculated automatically.
      The length of the radii and centers lists is one less than the length
      of the list of r and x data.
    Note that if tablized data is given, the first argument is ignored.
  """
  def __init__(self,rofxfunc,xmin,xmax,rmax=largepos,
                    voltage=0.,xcent=0.,ycent=0.,zcent=0.,condid=1,
                    rofxdata=None,xdata=None,raddata=None,
                    xcdata=None,rcdata=None,**kw):
    ZSrfrvOut.__init__(self,rofxfunc,xmin,xmax,rmax,
                       voltage,xcent,ycent,zcent,condid,
                       rofxdata,xdata,raddata,
                       xcdata,rcdata)
    XAssembly.__init__(self,voltage,xcent,ycent,zcent,
                            condid,self.kwlist,
                            self.generatorf,self.generatord,self.generatori,
                            kw=kw)

#============================================================================
class XSrfrvIn(ZSrfrvIn,XAssembly):
  """
Inside of a surface of revolution aligned along the X axis
  - rofxfunc: name of python function describing surface
  - xmin,xmax: x-extent of the surface
  - rmin=0: min radius of the surface
  - voltage=0: conductor voltage
  - xcent=0.,ycent=0.,zcent=0.: center of conductor
  - condid=1: conductor id of conductor, must be integer, or can be 'next' in
              which case a unique ID is chosen
  - rofxdata=None: optional tablized data of radius of surface
  - xdata=None: optional tablized data of x locations of rofxdata
  - raddata=None: optional radius of curvature of segments
      raddata[i] is radius for segment from xdata[i] to xdata[i+1]
  - xcdata=None: x center of circle or curved segment
  - rcdata=None: r center of circle or curved segment
      The centers of the circles will be calculated automatically if
      not supplied, or if the centers are supplied, the radii will be
      calculated automatically.
      The length of the radii and centers lists is one less than the length
      of the list of r and x data.
    Note that if tablized data is given, the first argument is ignored.
  """
  def __init__(self,rofxfunc,xmin,xmax,rmin=0,
                    voltage=0.,xcent=0.,ycent=0.,zcent=0.,condid=1,
                    rofxdata=None,xdata=None,raddata=None,
                    xcdata=None,rcdata=None,**kw):
    ZSrfrvIn.__init__(self,rofxfunc,xmin,xmax,rmin,
                      voltage,xcent,ycent,zcent,condid,
                      rofxdata,xdata,raddata,
                      xcdata,rcdata)
    XAssembly.__init__(self,voltage,xcent,ycent,zcent,condid,self.kwlist,
                            self.generatorf,self.generatord,self.generatori,
                            kw=kw)

#============================================================================
class XSrfrvInOut(ZSrfrvInOut,XAssembly):
  """
Between surfaces of revolution aligned along the X axis
  - rminofx,rmaxofx: names of python functions describing surfaces
  - xmin,xmax: x-extent of the surface
  - voltage=0: conductor voltage
  - xcent=0.,ycent=0.,zcent=0.: center of conductor
  - condid=1: conductor id of conductor, must be integer, or can be 'next' in
              which case a unique ID is chosen
  - rminofxdata,rmaxofxdata=None: optional tablized data of radii of surface
  - xmindata,xmaxdata=None: optional tablized data of x locations of r data
  - radmindata,radmaxdata=None: optional radius of curvature of segments
      radmindata[i] is radius for segment from xmindata[i] to xmindata[i+1]
  - xcmindata,xcmaxdata=None: x center of circle or curved segment
  - rcmindata,rcmaxdata=None: r center of circle or curved segment
      The centers of the circles will be calculated automatically if
      not supplied, or if the centers are supplied, the radii will be
      calculated automatically.
      The length of the radii and centers lists is one less than the length
      of the list of r and x data.
    Note that if tablized data is given, the first two arguments are ignored.
  """
  def __init__(self,rminofx,rmaxofx,xmin,xmax,
                    voltage=0.,xcent=0.,ycent=0.,zcent=0.,condid=1,
                    rminofxdata=None,xmindata=None,radmindata=None,
                    rcmindata=None,xcmindata=None,
                    rmaxofxdata=None,xmaxdata=None,radmaxdata=None,
                    rcmaxdata=None,xcmaxdata=None,**kw):
    ZSrfrvInOut.__init__(self,rminofx,rmaxofx,xmin,xmax,
                         voltage,xcent,ycent,zcent,condid,
                         rminofxdata,xmindata,radmindata,
                         rcmindata,xcmindata,
                         rmaxofxdata,xmaxdata,radmaxdata,
                         rcmaxdata,xcmaxdata)
    XAssembly.__init__(self,voltage,xcent,ycent,zcent,
                            condid,self.kwlist,
                            self.generatorf,self.generatord,self.generatori,
                            kw=kw)

#============================================================================
class YSrfrvOut(ZSrfrvOut,YAssembly):
  """
Outside of an surface of revolution aligned along to Y axis
  - rofyfunc: name of python function describing surface
  - ymin,ymax: y-extent of the surface
  - rmax=largepos: max radius of the surface
  - voltage=0: conductor voltage
  - xcent=0.,ycent=0.,zcent=0.: center of conductor
  - condid=1: conductor id of conductor, must be integer, or can be 'next' in
              which case a unique ID is chosen
  - rofydata=None: optional tablized data of radius of surface
  - ydata=None: optional tablized data of y locations of rofydata
      raddata[i] is radius for segment from ydata[i] to ydata[i+1]
  - ycdata=None: y center of circle or curved segment
  - rcdata=None: r center of circle or curved segment
  - raddata=None: optional radius of curvature of segments
      The centers of the circles will be calculated automatically if
      not supplied, or if the centers are supplied, the radii will be
      calculated automatically.
      The length of the radii and centers lists is one less than the length
      of the list of r and y data.
    Note that if tablized data is given, the first argument is ignored.
  """
  def __init__(self,rofyfunc,ymin,ymax,rmax=largepos,
                    voltage=0.,xcent=0.,ycent=0.,zcent=0.,condid=1,
                    rofydata=None,ydata=None,raddata=None,
                    ycdata=None,rcdata=None,**kw):
    ZSrfrvOut.__init__(self,rofyfunc,ymin,ymax,rmax,
                       voltage,xcent,ycent,zcent,condid,
                       rofydata,ydata,raddata,
                       ycdata,rcdata)
    YAssembly.__init__(self,voltage,xcent,ycent,zcent,
                            condid,self.kwlist,
                            self.generatorf,self.generatord,self.generatori,
                            kw=kw)

#============================================================================
class YSrfrvIn(ZSrfrvIn,YAssembly):
  """
Inside of a surface of revolution aligned along the Y axis
  - rofyfunc: name of python function describing surface
  - ymin,ymax: y-extent of the surface
  - rmin=0: min radius of the surface
  - voltage=0: conductor voltage
  - xcent=0.,ycent=0.,zcent=0.: center of conductor
  - condid=1: conductor id of conductor, must be integer, or can be 'next' in
              which case a unique ID is chosen
  - rofydata=None: optional tablized data of radius of surface
  - ydata=None: optional tablized data of y locations of rofydata
  - raddata=None: optional radius of curvature of segments
      raddata[i] is radius for segment from ydata[i] to ydata[i+1]
  - ycdata=None: y center of circle or curved segment
  - rcdata=None: r center of circle or curved segment
      The centers of the circles will be calculated automatically if
      not supplied, or if the centers are supplied, the radii will be
      calculated automatically.
      The length of the radii and centers lists is one less than the length
      of the list of r and y data.
    Note that if tablized data is given, the first argument is ignored.
  """
  def __init__(self,rofyfunc,ymin,ymax,rmin=0,
                    voltage=0.,xcent=0.,ycent=0.,zcent=0.,condid=1,
                    rofydata=None,ydata=None,raddata=None,
                    ycdata=None,rcdata=None,**kw):
    ZSrfrvIn.__init__(self,rofyfunc,ymin,ymax,rmin,
                      voltage,xcent,ycent,zcent,condid,
                      rofydata,ydata,raddata,
                      ycdata,rcdata)
    YAssembly.__init__(self,voltage,xcent,ycent,zcent,condid,self.kwlist,
                            self.generatorf,self.generatord,self.generatori,
                            kw=kw)

#============================================================================
class YSrfrvInOut(ZSrfrvInOut,YAssembly):
  """
Between surfaces of revolution aligned along the Y axis
  - rminofy,rmaxofy: names of python functions describing surfaces
  - ymin,ymax: y-extent of the surface
  - voltage=0: conductor voltage
  - xcent=0.,ycent=0.,zcent=0.: center of conductor
  - condid=1: conductor id of conductor, must be integer, or can be 'next' in
              which case a unique ID is chosen
  - rminofydata,rmaxofydata=None: optional tablized data of radii of surface
  - ymindata,ymaxdata=None: optional tablized data of y locations of r data
  - radmindata,radmaxdata=None: optional radius of curvature of segments
      radmindata[i] is radius for segment from ymindata[i] to ymindata[i+1]
  - ycmindata,ycmaxdata=None: y center of circle or curved segment
  - rcmindata,rcmaxdata=None: r center of circle or curved segment
      The centers of the circles will be calculated automatically if
      not supplied, or if the centers are supplied, the radii will be
      calculated automatically.
      The length of the radii and centers lists is one less than the length
      of the list of r and y data.
    Note that if tablized data is given, the first two arguments are ignored.
  """
  def __init__(self,rminofy,rmaxofy,ymin,ymax,
                    voltage=0.,xcent=0.,ycent=0.,zcent=0.,condid=1,
                    rminofydata=None,ymindata=None,radmindata=None,
                    rcmindata=None,ycmindata=None,
                    rmaxofydata=None,ymaxdata=None,radmaxdata=None,
                    rcmaxdata=None,ycmaxdata=None,**kw):
    ZSrfrvInOut.__init__(self,rminofy,rmaxofy,ymin,ymax,
                         voltage,xcent,ycent,zcent,condid,
                         rminofydata,ymindata,radmindata,
                         rcmindata,ycmindata,
                         rmayofxdata,ymaxdata,radmaxdata,
                         rcmaxdata,ycmaxdata)
    YAssembly.__init__(self,voltage,xcent,ycent,zcent,
                            condid,self.kwlist,
                            self.generatorf,self.generatord,self.generatori,
                            kw=kw)

#============================================================================
class ZAnnulus(ZSrfrvIn):
  """
Creates an Annulus as a surface of revolution.
  - rmin,rmax: Inner and outer radii
  - zmin,zmax: z-extent of the surface
  - voltage=0: conductor voltage
  - xcent=0.,ycent=0.,zcent=0.: center of conductor
  - condid=1: conductor id of conductor, must be integer, or can be 'next' in
              which case a unique ID is chosen
  """
  def __init__(self,rmin,rmax,length,
                    voltage=0.,xcent=0.,ycent=0.,zcent=0.,condid=1,**kw):

    self.rmin = rmin
    self.rmax = rmax
    self.length = length

    # --- Setup data for surface of revolution
    zmin = -length/2.
    zmax = +length/2.

    zdata = [zmin,zmax]
    rofzdata = [rmax,rmax]

    # --- ZSrfrvIn is a little faster than ZSrfrvInOut
    ZSrfrvIn.__init__(self,' ',zmin,zmax,rmin=rmin,
                      voltage=voltage,xcent=xcent,ycent=ycent,zcent=zcent,
                      condid=condid,
                      rofzdata=rofzdata,zdata=zdata,
                      kw=kw)

#============================================================================
class ZAnnulusElliptic(ZSrfrvEllipticIn,EllipticAssembly):
  """
Creates an Annulus as a surface of revolution.
  - ellipticity: ratio of y radius to x radius
  - rmin,rmax: Inner and outer radii
  - zmin,zmax: z-extent of the surface
  - voltage=0: conductor voltage
  - xcent=0.,ycent=0.,zcent=0.: center of conductor
  - condid=1: conductor id of conductor, must be integer, or can be 'next' in
              which case a unique ID is chosen
  """
  def __init__(self,ellipticity,rmin,rmax,length,
                    voltage=0.,xcent=0.,ycent=0.,zcent=0.,condid=1,**kw):

    self.rmin = rmin
    self.rmax = rmax
    self.length = length

    # --- Setup dat for surface of revolution
    zmin = -length/2.
    zmax = +length/2.

    zdata = [zmin,zmax]
    rofzdata = [rmax,rmax]

    # --- ZSrfrvEllipticIn is a little faster than ZSrfrvEllipticInOut
    ZSrfrvEllipticIn.__init__(self,ellipticity,' ',zmin,zmax,rmin=rmin,
                              voltage=voltage,
                              xcent=xcent,ycent=ycent,zcent=zcent,
                              condid=condid,
                              rofzdata=rofzdata,zdata=zdata,
                              kw=kw)

#============================================================================
class ZCone(ZSrfrvIn):
  """
Cone
  - r_zmin: radius at z min
  - r_zmax: radius at z max
  - length: length
  - voltage=0: cone voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cone
  - condid=1: conductor id of cone, must be integer, or can be 'next' in
              which case a unique ID is chosen
  """
  def __init__(self,r_zmin,r_zmax,length,
                    voltage=0.,xcent=0.,ycent=0.,zcent=0.,condid=1,**kw):

    self.r_zmin = r_zmin
    self.r_zmax = r_zmax
    self.length = length

    zmin = -length/2.
    zmax = +length/2.
    zdata = [zmin,zmax]
    rofzdata = [r_zmin,r_zmax]

    ZSrfrvIn.__init__(self,' ',zmin,zmax,
                      voltage=voltage,xcent=xcent,ycent=ycent,zcent=zcent,
                      condid=condid,
                      rofzdata=rofzdata,zdata=zdata,
                      kw=kw)

#============================================================================
class ZConeSlope(ZSrfrvIn):
  """
Cone
  - slope: ratio of radius at zmax minus radius at zmin over length
  - intercept: location where line defining cone crosses the axis, relative
               to zcent
  - length: length
  - voltage=0: cone voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cone
  - condid=1: conductor id of cone, must be integer, or can be 'next' in
              which case a unique ID is chosen
  """
  def __init__(self,slope,intercept,length,
                    voltage=0.,xcent=0.,ycent=0.,zcent=0.,condid=1,**kw):

    self.slope = slope
    self.intercept = intercept
    self.length = length

    r_zmin = self.slope*(-self.length/2. - self.intercept)
    r_zmax = self.slope*(+self.length/2. - self.intercept)
    zmin = -length/2.
    zmax = +length/2.
    zdata = [zmin,zmax]
    rofzdata = [r_zmin,r_zmax]

    ZSrfrvIn.__init__(self,' ',zmin,zmax,
                      voltage=voltage,xcent=xcent,ycent=ycent,zcent=zcent,
                      condid=condid,
                      rofzdata=rofzdata,zdata=zdata,
                      kw=kw)

#============================================================================
class ZConeOut(ZSrfrvOut):
  """
Cone outside
  - r_zmin: radius at z min
  - r_zmax: radius at z max
  - length: length
  - voltage=0: cone voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cone
  - condid=1: conductor id of cone, must be integer, or can be 'next' in
              which case a unique ID is chosen
  """
  def __init__(self,r_zmin,r_zmax,length,
                    voltage=0.,xcent=0.,ycent=0.,zcent=0.,condid=1,**kw):

    self.r_zmin = r_zmin
    self.r_zmax = r_zmax
    self.length = length

    zmin = -length/2.
    zmax = +length/2.
    zdata = [zmin,zmax]
    rofzdata = [r_zmin,r_zmax]

    ZSrfrvOut.__init__(self,' ',zmin,zmax,
                       voltage=voltage,xcent=xcent,ycent=ycent,zcent=zcent,
                       condid=condid,
                       rofzdata=rofzdata,zdata=zdata,
                       kw=kw)

#============================================================================
class ZConeOutSlope(ZSrfrvOut):
  """
Cone outside
  - slope: ratio of radius at zmax minus radius at zmin over length
  - intercept: location where line defining cone crosses the axis, relative
               to zcent
  - length: length
  - voltage=0: cone voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cone
  - condid=1: conductor id of cone, must be integer, or can be 'next' in
              which case a unique ID is chosen
  """
  def __init__(self,slope,intercept,length,
                    voltage=0.,xcent=0.,ycent=0.,zcent=0.,condid=1,**kw):

    self.slope = slope
    self.intercept = intercept
    self.length = length

    r_zmin = self.slope*(-self.length/2. - self.intercept)
    r_zmax = self.slope*(+self.length/2. - self.intercept)
    zmin = -length/2.
    zmax = +length/2.
    zdata = [zmin,zmax]
    rofzdata = [r_zmin,r_zmax]

    ZSrfrvOut.__init__(self,' ',zmin,zmax,
                       voltage=voltage,xcent=xcent,ycent=ycent,zcent=zcent,
                       condid=condid,
                       rofzdata=rofzdata,zdata=zdata,
                       kw=kw)

#============================================================================
#============================================================================
#============================================================================
def Quadrupole(ap=None,rl=None,rr=None,gl=None,gp=None,
               pa=None,pw=None,pr=None,vx=None,vy=None,
               glx=None,gly=None,axp=None,axm=None,ayp=None,aym=None,
               rxp=None,rxm=None,ryp=None,rym=None,
               vxp=None,vxm=None,vyp=None,vym=None,
               oxp=None,oxm=None,oyp=None,oym=None,
               pwl=None,pwr=None,pal=None,par=None,prl=None,prr=None,
               xcent=0.,ycent=0.,zcent=None,condid=None,splitrodids=false,
               elemid=None,elem='quad',**kw):
  """
Creates an interdigited quadrupole structure.
Either specify the quadrupole structure...
  - ap: pole tip aperture
  - rl: rod length
  - rr: rod radius
  - gl: gap length between rod end and end plate
  - gp: sign of gap location in x plane
  - pa: aperture of end plate
  - pw: width of end plate
  - pr: outer radius of end plate
  - vx: voltage of rod in x plane
  - vy: voltage of rod in y plane
  - glx: Change in gap length on x axis
  - gly: Change in gap length on y axis
  - axp: Change in aperture of rod on plus  x axis
  - axm: Change in aperture of rod on minus x axis
  - ayp: Change in aperture of rod on plus  y axis
  - aym: Change in aperture of rod on minus y axis
  - rxp: Change in radius of rod on plus  x axis
  - rxm: Change in radius of rod on minus x axis
  - ryp: Change in radius of rod on plus  y axis
  - rym: Change in radius of rod on minus y axis
  - vxp: Change in voltage of rod on plus  x axis
  - vxm: Change in voltage of rod on minus x axis
  - vyp: Change in voltage of rod on plus  y axis
  - vym: Change in voltage of rod on minus y axis
  - oxp: Perpendicular offset of rod on plus  x axis
  - oxm: Perpendicular offset of rod on minus x axis
  - oyp: Perpendicular offset of rod on plus  y axis
  - oym: Perpendicular offset of rod on minus y axis
  - pwl: Change on left  plate width
  - pwr: Change on right plate width
  - pal: Change on left  plate aperture
  - par: Change on right plate aperture
  - prl: Change on left  plate max radius
  - prr: Change on right plate max radius
  - xcent=0.,ycent=0.: transverse center of quadrupole
  - zcent: axial center of quadrupole, default taken from element
  - condid=1: conductor id of quadrupole, must be integer
  - splitrodids=false: when true, the condid's of the x and y rods are
                       different, y is the negative of x (which is condid)
Or give the quadrupole id to use...
  - elem='quad': element type to get data from
  - elemid: gets data from quad element. The above quantities can be specified
            as well to override the value from elemid.
  """
  if elemid is None:
    assert ap is not None,'ap must be specified'
    assert rl is not None,'rl must be specified'
    assert rr is not None,'rr must be specified'
    assert gl is not None,'gl must be specified'
    assert gp is not None,'gp must be specified'
    assert pa is not None,'pa must be specified'
    assert pw is not None,'pw must be specified'
    assert pr is not None,'pr must be specified'
    assert vx is not None,'vx must be specified'
    assert vy is not None,'vy must be specified'
    if zcent is None: zcent = 0.
    if condid is None: condid = 1
  else:
    if ap is None: ap = getattr(top,elem+'ap')[elemid]
    if rl is None: rl = getattr(top,elem+'rl')[elemid]
    if rr is None: rr = getattr(top,elem+'rr')[elemid]
    if gl is None: gl = getattr(top,elem+'gl')[elemid]
    if gp is None: gp = getattr(top,elem+'gp')[elemid]
    if pa is None: pa = getattr(top,elem+'pa')[elemid]
    if pw is None: pw = getattr(top,elem+'pw')[elemid]
    if elem == 'quad':
      if pr is None: pr = getattr(top,elem+'pr')[elemid]
      if vx is None: vx = top.quadvx[elemid]
      if vy is None: vy = top.quadvy[elemid]
      if xcent is None: xcent = top.qoffx[elemid]
      if ycent is None: ycent = top.qoffy[elemid]
    else:
      if vx is None: vx = 0.
      if vy is None: vy = 0.
      if xcent is None: xcent = getattr(top,elem+'ox')[elemid]
      if ycent is None: ycent = getattr(top,elem+'oy')[elemid]
    if zcent is None:
      zcent = 0.5*(getattr(top,elem+'zs')[elemid] +
                   getattr(top,elem+'ze')[elemid])
    if condid is None: condid = elemid

  dels = ['glx','gly','axp','axm','ayp','aym','rxp','rxm','ryp','rym',
          'vxp','vxm','vyp','vym','oxp','oxm','oyp','oym',
          'pwl','pwr','pal','par','prl','prr']
  if elemid is None or elem != 'quad':
    for d in dels:
      if locals()[d] is None: exec('%s = 0.'%d)
  else:
    edel = 'qdel'
    for d in dels:
      if locals()[d] is None: exec('%s = top.qdel%s[elemid]'%(d,d))

  if splitrodids:
    xidsign = +1
    yidsign = -1
  else:
    xidsign = +1
    yidsign = +1

  # --- Create x and y rods
  if ap > 0. and rr > 0. and rl > 0.:
    xrod1 = ZCylinder(rr+rxp,rl-glx,vx+vxp,xcent+ap+rr+axp,ycent+oxp,
                      zcent-gp*gl/2.,xidsign*condid,
                      kw=kw)
    xrod2 = ZCylinder(rr+rxm,rl-glx,vx+vxm,xcent-ap-rr-axm,ycent+oxm,
                      zcent-gp*gl/2.,xidsign*condid,
                      kw=kw)
    yrod1 = ZCylinder(rr+ryp,rl-gly,vy+vyp,xcent+oyp,ycent+ap+rr+ayp,
                      zcent+gp*gl/2.,yidsign*condid,
                      kw=kw)
    yrod2 = ZCylinder(rr+rym,rl-gly,vy+vym,xcent+oym,ycent-ap-rr-aym,
                      zcent+gp*gl/2.,yidsign*condid,
                      kw=kw)
    quad = xrod1 + xrod2 + yrod1 + yrod2
  else:
    quad = None

  # --- Add end plates
  if pw > 0. and (ap > 0. or pa > 0.):
    if pa == 0.: pa = ap
    if pr == 0.: pr = 2*w3d.xmmax
    if gp > 0.:
      v1 = vx+vxp
      v2 = vy+vyp
      gp = +1
    else:
      v1 = vy+vyp
      v2 = vx+vxp
      gp = -1
    if splitrodids:
      lidsign = +gp
      ridsign = -gp
    else:
      lidsign = +1
      ridsign = +1
    if pr < 1.4142*w3d.xmmax:
      plate1 = ZAnnulus(pa+pal,pr+prl,pw+pwl,v1,xcent,ycent,
                        zcent-0.5*(rl+gl)-pw/2.,lidsign*condid,
                        kw=kw)
      plate2 = ZAnnulus(pa+par,pr+prr,pw+pwr,v2,xcent,ycent,
                        zcent+0.5*(rl+gl)+pw/2.,ridsign*condid,
                        kw=kw)
    else:
      plate1 = ZCylinderOut(pa+pal,pw+pwl,v1,xcent,ycent,
                            zcent-0.5*(rl+gl)-pw/2.,lidsign*condid,
                            kw=kw)
      plate2 = ZCylinderOut(pa+par,pw+pwr,v2,xcent,ycent,
                            zcent+0.5*(rl+gl)+pw/2.,ridsign*condid,
                            kw=kw)
    quad = quad + plate1 + plate2

  return quad

#============================================================================
#============================================================================
#============================================================================
#============================================================================
#============================================================================
#============================================================================
try:
  enumerate
except:
  def enumerate(ll):
    tt = []
    for i in range(len(ll)):
      tt.append((i,ll[i]))
    return tt

class SRFRVLApart:
  """
Class for creating a surface of revolution conductor part using lines and arcs as primitives.
 - name:   name of conductor part
 - data:   list of lines and arcs (this must be a continuous curve)
            o the first element is the starting point
            o a line is described by the ending point
            o a circle is described by the ending point and the center
            o for both, the starting point is the previous ending point
            o each element is teminated by a letter 't' for top (or out)
              and 'b' for bottom (or in)
  """
  def __init__(self,name,data):
    self.name    = name
    self.data   = data

  def draw(self,ncirc=50,scx=1.,scy=1.,colort='blue',colorb='red',
                 color='none',signx=1.,width=1.):
    """
  Draw lines and arcs of SRFRVLA conductor part. Arcs are decomposed
  into ncirc segments.
    """
    if(color <> 'none'):
      colort=color
      colorb=color

    for j,d in enumerate(self.data):
      if(j==0):
        s=[d[0]*scx,d[1]*scy]
      else:
        if(d[-1]=='t'):
          colorn=colort
        else:
          colorn=colorb
        n=[d[0]*scx,d[1]*scy]
        if(len(d)==3):
          self.draw_line([s[0]],[signx*s[1]],[n[0]],[signx*n[1]],colorn,width)
        else:
          c=[d[2]*scx,d[3]*scy]
          self.draw_arc(s[0],signx*s[1],n[0],signx*n[1],c[0],signx*c[1],ncirc,colorn,width)
        s=n

  def shift(self,s,dir):
    """
  Shift SRFRVLA conductor part by s in the direction dir (0=Z, 1=R).
    """
    for d in self.data:
      d[dir]+=s
      if(len(d)>3):
        d[dir+2]+=s

  def draw_line(self,x0,y0,x1,y1,color='black',width=1.):
    """
  Draw a line.
    """
    pldj(x0,y0,x1,y1,color=color,width=width)

  def draw_arc(self,x0,y0,x1,y1,xc,yc,ncirc=50,color='black',width=1.):
    """
  Draw an arc as ncirc segments.
    """
    xy = self.get_xy_arc(x0,y0,x1,y1,xc,yc,ncirc)
    pldj(xy[0][0:ncirc-1],xy[1][0:ncirc-1],
         xy[0][1:ncirc],  xy[1][1:ncirc],
         color=color,width=width)

  def get_xy_arc(self,x0,y0,x1,y1,xc,yc,ncirc):
    """
  Retuns list of end points of ncirc segments along an arc of
  end points (x0,y0) and (x1,y1), and center (xc,yc).
    """
    xi=x0-xc
    yi=y0-yc
    xf=x1-xc
    yf=y1-yc
    r=sqrt(xi**2+yi**2)
    thetai=arctan2(yi,xi)
    thetaf=arctan2(yf,xf)
    if(abs(thetaf-thetai)>pi):
      if(thetaf-thetai<0.):
        thetaf=thetaf+2.*pi
      else:
        thetaf=thetaf-2.*pi
    dtheta=(thetaf-thetai)/(ncirc-1)
    x=zeros(ncirc+1,'d')
    y=zeros(ncirc+1,'d')
    x[0]=x0
    y[0]=y0
    for i in range(1,ncirc):
      theta=thetai+i*dtheta
      x[i]=xc+r*cos(theta)
      y[i]=yc+r*sin(theta)
    return [x,y]


class SRFRVLAcond:
  """
Class for creating Surface of revolution conductor using lines and arcs as primitives.
A conductor is a structure containing a list of parts. A part is a structure
containing a list of primitives.
 - name: name of conductor
 - parts: list of conductor parts
 - voltage: conductor part voltage
 - condid:  conductor part ID
 - install=1 (optional): flag for installation of conductors immediately after creation.
  """
  def __init__(self,name,parts,voltage,condid,install=1,l_verbose=1):
    """
  Register SRFRVLA conductor parts using ZSrfrvIn, ZSrfrvOut or ZSrfrvInOut.
    """
    self.name = name
    self.parts = parts
    self.voltage = voltage
    self.condid = condid
    for i,part in enumerate(self.parts):
    # Loop over parts and install them using ZSrfrvIn, ZSrfrvOut or ZSrfrvInOut.
    # The installed parts are stored in the list 'toinstall' and the corresponding
    # IDs and names in the lists 'condids' and 'cnames'.
      data = part.data
      it=0
      ib=0
      nsrmax=0
      nsrmin=0
      l_t=0
      l_b=0
      # Counts number of 'top' and 'bottom' primitives and allocate temporary arrays.
      for d in data[1:]:
        if(d[-1]=='t'):
            nsrmax += 1
        elif(d[-1]=='b'):
            nsrmin += 1
      if(nsrmax<>0):
        l_t = 1
        r_srmax = zeros(nsrmax+1,Float)
        z_srmax = zeros(nsrmax+1,Float)
        rc_srmax = zeros(nsrmax,Float)
        zc_srmax = zeros(nsrmax,Float)
      if(nsrmin<>0):
        l_b = 1
        r_srmin = zeros(nsrmin+1,Float)
        z_srmin = zeros(nsrmin+1,Float)
        rc_srmin = zeros(nsrmin,Float)
        zc_srmin = zeros(nsrmin,Float)

      # Fill arrays with datas from parts.
      do = data[0]
      for d in data[1:]:
        if(d[-1]=='t'):
          z_srmax[it]=do[0]
          r_srmax[it]=do[1]
          it=it+1
          z_srmax[it]=d[0]
          r_srmax[it]=d[1]
          if(len(d)>3):
            zc_srmax[it-1]=d[2]
            rc_srmax[it-1]=d[3]
          else:
            zc_srmax[it-1]=largepos
            rc_srmax[it-1]=largepos
        elif(d[-1]=='b'):
          z_srmin[ib]=do[0]
          r_srmin[ib]=do[1]
          ib=ib+1
          z_srmin[ib]=d[0]
          r_srmin[ib]=d[1]
          if(len(d)>3):
            zc_srmin[ib-1]=d[2]
            rc_srmin[ib-1]=d[3]
          else:
            zc_srmin[ib-1]=largepos
            rc_srmin[ib-1]=largepos
        do=d

      # Make sure arrays are in Z ascending order.
      if(z_srmin[0]>z_srmin[nsrmin]):
        args =argsort(z_srmin)
        argsc=argsort(z_srmin[1:])
        z_srmin=take(z_srmin,args)
        r_srmin=take(r_srmin,args)
        zc_srmin=take(zc_srmin,argsc)
        rc_srmin=take(rc_srmin,argsc)
      if(z_srmax[0]>z_srmax[nsrmax]):
        args =argsort(z_srmax)
        argsc=argsort(z_srmax[1:])
        z_srmax=take(z_srmax,args)
        r_srmax=take(r_srmax,args)
        zc_srmax=take(zc_srmax,argsc)
        rc_srmax=take(rc_srmax,argsc)

      # Register parts.
      if(l_t==1 and l_b==1):
        part.installed = ZSrfrvInOut('','',
                              min(z_srmax),
                              max(z_srmax),
                              voltage=voltage,
                              condid =condid,
                              rminofzdata=r_srmin,
                              rmaxofzdata=r_srmax,
                              zmindata=z_srmin,
                              zmaxdata=z_srmax,
                              rcmindata=rc_srmin,
                              rcmaxdata=rc_srmax,
                              zcmindata=zc_srmin,
                              zcmaxdata=zc_srmax,
                              )
        part.zmin = min(z_srmax)
        part.zmax = max(z_srmax)
        part.rmin = min(r_srmin)
        part.rmax = max(r_srmax)

      elif(l_b==1):
        part.installed = ZSrfrvOut('',
                            min(z_srmin),
                            max(z_srmin),
                            voltage=voltage,
                            condid =condid,
                            rofzdata=r_srmin,
                            zdata=z_srmin)
        part.zmin = min(z_srmin)
        part.zmax = max(z_srmin)
        part.rmin = min(r_srmin)
        part.rmax = max(r_srmin)

      elif(l_t==1):
        part.installed = ZSrfrvIn('',
                            min(z_srmax),
                            max(z_srmax),
                            voltage=voltage,
                            condid =condid,
                            rofzdata=r_srmax,
                            zdata=z_srmax)
        part.zmin = min(z_srmax)
        part.zmax = max(z_srmax)
        part.rmin = min(r_srmax)
        part.rmax = max(r_srmax)

      # store installed conductor in a list
      if i == 0:
        self.cond = self.parts[i].installed
      else:
        self.cond += self.parts[i].installed
    if(install):self.install(l_verbose=l_verbose)

  def install(self,grid=None,l_verbose=1,l_recursive=1):
    """
  Install SRFRVLA conductors.
    """
    if l_verbose:print 'installing',self.name,'( ID=',self.condid,')...'
    if w3d.solvergeom==w3d.RZgeom or w3d.solvergeom==w3d.XZgeom:
        if grid is None:grid=frz.basegrid
        print grid.gid
        for part in self.parts:
            installconductors(part.installed,xmin=part.rmin,xmax=part.rmax,
                        zmin=part.zmin,zmax=part.zmax,
                         nx=grid.nr,nzlocal=grid.nzlocal,nz=grid.nz,
                         xmmin=grid.xmin,xmmax=grid.xmax,
                         zmmin=grid.zmin,zmmax=grid.zmax,
                        gridrz=grid)
        if(l_recursive):
          try:
            self.install(grid.next,l_verbose=0)
          except:
            try:
              self.install(grid.down,l_verbose=0)
            except:
              pass
    else:
      for part in self.parts:
        installconductors(part.installed)

  def draw(self,ncirc=50,scx=1.,scy=1.,colort='blue',colorb='red',
                 color='none',signx=1.,width=1.):
    """
  Draws a list of conductor parts.
    """
    if(color <> 'none'):
      colort=color
      colorb=color
    for part in self.parts:
      part.draw(ncirc,scx,scy,colort,colorb,color,signx,width)

  def shift(self,s,dir):
    """
  Shifts conductor parts using a scalar shift of a list of shifts.
    - s: shift value (scalar or list same length as self.parts)
    - dir: direction of shift (0=Z, 1=R)
    """
    stype = 'list'
    try:
      l=len(s)
    except TypeError:
      stype = 'scalar'

    if(stype=='scalar'):
      for part in self.parts:
        part.shift(s,dir)
    else:
      assert len(s)==len(self.parts), 'Error in SRFRVLAcond.shift: shift must have same length as part or be a scalar.'
      for i,part in enumerate(self.parts):
        part.shift(s[i],dir)

class SRFRVLAsystem:
  """
Class for creating a SRFRVLAsystem of conductors (list of conductors).
A SRFRVLA contains two lists:
  - SRFRVLAconds which is the list of SRFRVLA conductors (SRFRVLAcond).
  - conds which is the list of conductors, as defined in generateconductors.py.
  """
  def __init__(self,SRFRVLAconds):
    self.conds = []
    self.SRFRVLAconds = SRFRVLAconds
    for SRFRVLAcond in SRFRVLAconds:
      self.conds += [SRFRVLAcond.cond]

  def save(self,filename):
    """
  Save conductor into external file:
    - filename: name of external file ('.wob' will be added')
    """
    f=open(filename+'.wob','w')
    f.write("Begin\n")

    for cond in self.SRFRVLAconds:
      f.write("  Conductor "+cond.name+"\n")
      for part in cond.parts:
        f.write("    Part "+part.name+"\n")
        d = part.data[0]
        f.write('     s   %13.6E %13.6E \n'%(d[0],d[1]))
        for d in part.data[1:]:
          if(len(d)==3): # line
            f.write('     l   %13.6E %13.6E                             %c \n'%(d[0],d[1],d[2]))
          else: # arc
            f.write('     a   %13.6E %13.6E %13.6E %13.6E %c \n'%(d[0],d[1],d[2],d[3],d[4]))
        f.write("    Endpart "+part.name+"\n")
      f.write("  End "+cond.name+"\n")

    f.write("End\n")
    f.close()

  def draw(self,ncirc=50,scx=1.,scy=1.,colort='blue',colorb='red',
                 color='none',signx=1.,width=1.):
    """
  Draws a list of conductors.
    """
    if(color <> 'none'):
      colort=color
      colorb=color
    for parts in self.SRFRVLAconds:
        parts.draw(ncirc,scx,scy,colort,colorb,color,signx,width)

  def getextent(self):
    return ConductorExtent([w3d.xmmin,w3d.ymmin,w3d.zmmin],[w3d.xmmax,w3d.ymmax,w3d.zmmax])

class SRFRVLAfromfile(SRFRVLAsystem):
  """
Class for reading SRFRVLA data from file.
  - filename: name of input file
  - voltages: list of voltages of conductor parts
  - condids: list of IDs of conductor parts
  - zshifts=None (optional): list of shifts in Z to apply to conductor parts
  - rshifts=None (optional): list of shifts in R to apply to conductor parts
  - install=1 (optional): flag for installation of conductors immediately after reading
  """
  def __init__(self,filename,voltages,condids,zshifts=None,rshifts=None,install=1,l_verbose=1):
    """
  Reads SRFRVLA conductors from external file. The series of conductors which is stored in a list.
    """
    f=open(filename+'.wob','r')
    cont=1
    ic = -1
    ip = -1
    # create list conductors
    conds=[]
    condnames = []
    while(cont==1):
      line=split(f.readline())
      if(line[0]=='Conductor'):
        ic += 1
        condnames += [line[1]]
        partnames = []
        parts = []
        # create list of data
        data = []
        ocont = 1
        while(ocont==1):
          line=split(f.readline())
          if(line[0]=='End'):
            ocont=0
          elif(line[0]=='Part'):
            ip += 1
            partnames += [line[1]]
          elif(line[0]=='Endpart'):
            parts += [SRFRVLApart(partnames[-1],data)]
            if zshifts is not None:parts[-1].shift(zshifts[ip],0)
            if rshifts is not None:parts[-1].shift(rshifts[ip],1)
            data = []
            partnames = []
          else:
            # append data to list
            if(line[0]=='S' or line[0]=='s'): #start
              data.append([float(line[1]),float(line[2])])
            if(line[0]=='L' or line[0]=='l'): #line
              data.append([float(line[1]),float(line[2]),line[3]])
            if(line[0]=='A' or line[0]=='a'): #arc
              data.append([float(line[1]),float(line[2]),\
                           float(line[3]),float(line[4]),line[5]])
        conds += [parts]
      elif(line[0]=='End'):
        cont=0
    f.close()
    self.conductors = conds
    self.condnames = condnames
    self.voltages = voltages
    self.condids = condids
    self.SRFRVLAconds = []
    for i,parts in enumerate(self.conductors):
      self.SRFRVLAconds += [SRFRVLAcond(self.condnames[i],parts,self.voltages[i],self.condids[i],install=0)]
    SRFRVLAsystem.__init__(self,self.SRFRVLAconds)
    if install:self.install(l_verbose=l_verbose)

  def install(self,l_verbose=1,grid=None,l_recursive=1):
    # Install conductors
    for conductor in self.SRFRVLAconds:
        conductor.install(l_verbose=l_verbose,grid=grid,l_recursive=l_recursive)

class SRFRVLA_circle(SRFRVLApart):
  """
Creates a circle.
  - name: name of conductor part
  - c: center
  - r: radius
  """
  def __init__(self,name,c,r):
    assert r>0., 'In SRFRVLA_circle, r must be >0.'
    x = c[0]
    y = c[1]
    data = []
    data.append([x-r,y])
    data.append([x+r,y,x,y,'t'])
    data.append([x,y-r,x,y,'b'])
    data.append([x-r,y,x,y,'b'])
    SRFRVLApart.__init__(self,name,data)

class SRFRVLA_rectangle(SRFRVLApart):
  """
Creates a rectangle.
  - name: name of conductor part
  - c: center
  - l: length
  - h: height
  """
  def __init__(self,name,c,l,h):
    assert l>0. and h>0., 'In SRFRVLA_rectangle, l, h must be >0.'
    x1 = c[0]-0.5*l
    x2 = c[0]+0.5*l
    y1 = c[1]-0.5*h
    y2 = c[1]+0.5*h
    data = []
    data.append([x1,y1])
    data.append([x1,y2,'t'])
    data.append([x2,y2,'t'])
    data.append([x2,y1,'t'])
    data.append([x1,y1,'b'])
    SRFRVLApart.__init__(self,name,data)

class SRFRVLA_rnd_rectangle(SRFRVLApart):
  """
Creates a rectangle with rounded edge.
  - name: name of conductor part
  - c: center
  - l: length
  - h: height
  - r: radius of rounded edges
  """
  def __init__(self,name,c,l,h,r):
    assert r>=0. and l>=0. and h>=0., 'In SRFRVLA_rnd_rectangle, r, l, h must be >=0.'
    l = max(l,2.*r)
    h = max(h,2.*r)
    x1 = c[0]-0.5*l
    x2 = c[0]+0.5*l
    y1 = c[1]-0.5*h
    y2 = c[1]+0.5*h
    data = []
    data.append([x1,y1+r])
    data.append([x1,y2-r,'t'])
    if(r>0.):data.append([x1+r,y2,x1+r,y2-r,'t'])
    data.append([x2-r,y2,'t'])
    if(r>0.):data.append([x2,y2-r,x2-r,y2-r,'t'])
    data.append([x2,y1+r,'t'])
    if(r>0.):data.append([x2-r,y1,x2-r,y1+r,'b'])
    data.append([x1+r,y1,'b'])
    if(r>0.):data.append([x1,y1+r,x1+r,y1+r,'b'])
    SRFRVLApart.__init__(self,name,data)

