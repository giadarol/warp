"""Utility routines for doing 3-d visualization using the VPython module.
Partly taken from VPython demos.
Modified by DPG

VisualMesh: can plot 3-D surfaces corresponding to meshed data.
"""
from warp import *
VPythonobjects_version = "$Id: VPythonobjects.py,v 1.10 2004/05/20 20:11:40 dave Exp $"

def VPythonobjectsdoc():
  import VPythonobjects
  print VPythonobjects.__doc__

##########################################################################

class VisualModel:
  def __init__(self,twoSided=1,normalsign=1,scene=None,
                    title='Visualization',labels=None,
                    vrange=None,viewer=None):
    if viewer is None: viewer = 'OpenDX'
    self.viewer = viewer
    self.triangles = []
    self.colors = []
    self.normals = []
    self.connections = None
    self.twoSided = twoSided  # add every face twice with opposite normals
    self.normalsign = normalsign
    self.defcolor = [0.5,0.7,1.0]
    if vrange is not None:
      self.vrange = vrange
      self.autoscale = 0
      self.uniform = 1
    else:
      self.vrange = [10.,10.,10.]
      self.autoscale = 1
      self.uniform = 1
    self.title = title
    self.labels = labels
    self.scene = scene
    self.dxobject = None

  def getdxobject(self):
    try:
      return self.dxobject
    except AttributeError:
      self.CreateDXObject()
      return dxobject

  def CreateDXObject(self):

    import pyOpenDX

    n = len(self.triangles)
    p = array(self.triangles).astype(Float32)
    ss = pyOpenDX.DXNewArray(pyOpenDX.TYPE_FLOAT,pyOpenDX.CATEGORY_REAL,1,3)
    pyOpenDX.DXAddArrayData(ss,0,n,p)

    normals = array(self.normals)*self.normalsign
    nn = pyOpenDX.DXNewArray(pyOpenDX.TYPE_FLOAT,pyOpenDX.CATEGORY_REAL,1,3)
    pyOpenDX.DXAddArrayData(nn,0,n,normals.astype(Float32))
    pyOpenDX.DXSetStringAttribute(nn,'dep','positions')

    if not self.colors: colors = array(n*[[0.5,0.7,1.0]])
    else:               colors = array(self.colors)
    co = pyOpenDX.DXNewArray(pyOpenDX.TYPE_FLOAT,pyOpenDX.CATEGORY_REAL,1,3)
    pyOpenDX.DXAddArrayData(co,0,n,colors.astype(Float32))
    pyOpenDX.DXSetStringAttribute(co,'dep','positions')

    if self.connections is None:
      self.connections = arange(n)
      self.connections.shape = (nint(n/3),3)
    cc = pyOpenDX.DXNewArray(pyOpenDX.TYPE_INT,pyOpenDX.CATEGORY_REAL,1,3)
    pyOpenDX.DXAddArrayData(cc,0,nint(n/3),self.connections.astype(Int))
    pyOpenDX.DXSetStringAttribute(cc,'ref','positions')
    pyOpenDX.DXSetStringAttribute(cc,'element type','triangles')
    pyOpenDX.DXSetStringAttribute(cc,'dep','connections')

    ff = pyOpenDX.DXNewField()
    pyOpenDX.DXSetComponentValue(ff,'positions',ss)
    pyOpenDX.DXSetComponentValue(ff,'colors',co)
    pyOpenDX.DXSetComponentValue(ff,'normals',nn)
    pyOpenDX.DXSetComponentValue(ff,'connections',cc)
    pyOpenDX.DXEndField(ff)

    self.dxobject = ff

  def Display(self,showgrid=0,labels=None):

    if self.viewer == 'VPython':
      import visual

      if self.scene is None:
        self.scene = visual.display(exit=0,width=500,height=500,
                                    uniform=self.uniform,title=self.title,
                                    autoscale=self.autoscale,range=self.vrange)

      self.frame = visual.frame(display=self.scene)
      if not self.colors: self.colors=None
      self.model = visual.faces(frame=self.frame,pos=self.triangles,
                                normal=self.normals,color=self.colors,
                                display=self.scene)
    else:
      import pyOpenDX
      if self.dxobject is None: self.CreateDXObject()
      ff = self.dxobject

      if showgrid:

        # --- Create a field containing two points on opposite side of the
        # --- mesh and use it to find the boinding box.
        n = 2
        p = zeros((n,3),'d')
        p[:,0] = [w3d.xmmin,w3d.xmmax]
        p[:,1] = [w3d.xmmin,w3d.xmmax]
        p[:,2] = [w3d.xmmin,w3d.xmmax]
        dxp = pyOpenDX.DXNewArray(pyOpenDX.TYPE_FLOAT,pyOpenDX.CATEGORY_REAL,1,3)
        pyOpenDX.DXAddArrayData(dxp,0,n,p.astype(Float32))
        # --- Create a DX array for the data
        dxd = pyOpenDX.DXNewArray(pyOpenDX.TYPE_DOUBLE,pyOpenDX.CATEGORY_REAL,0)
        pyOpenDX.DXAddArrayData(dxd,0,n,array([1.,1.]))
        pyOpenDX.DXSetStringAttribute(dxd,'dep','positions')
        # --- Create the field
        dxf = pyOpenDX.DXNewField()
        pyOpenDX.DXSetComponentValue(dxf,'positions',dxp)
        pyOpenDX.DXSetComponentValue(dxf,'data',dxd)
        pyOpenDX.DXEndField(dxf)

        minput = {'input':dxf}
        moutput = ['box']
        (box,) = pyOpenDX.DXCallModule('ShowBox',minput,moutput)

        minput = {'object':ff,'object1':box}
        moutput = ['group']
        (group,) = pyOpenDX.DXCallModule('Collect',minput,moutput)

      else:
        group = ff

      if labels is None: labels = self.labels

      pyOpenDX.DXImage(group,name=self.title,labels=labels)
      pyOpenDX.DXDelete(group)

  def FacetedTriangle(self, vv, nn = None, color=None):
    """Add a triangle to the model, apply faceted shading automatically"""
    if nn is None: nn = len(vv)*[None]
    if color is None: color = self.defcolor
    normal = self.Norm( self.Cross(vv[1]-vv[0], vv[2]-vv[0]) )
    for v,n in zip(vv,nn):
      self.triangles.append(v)
      self.colors.append(color)
      if n is None: self.normals.append(normal)
      else:         self.normals.append(n)
      #self.model.append( pos=v, color=color, normal=normal )
    if self.twoSided:
      for v,n in zip(vv,nn):
        #self.model.append( pos=v, color=color, normal=-normal )
        self.triangles.append(v)
        self.colors.append(color)
        if n is None: self.normals.append(-normal)
        else:         self.normals.append(-n)

  def FacetedPolygon(self, v, n=None, color=None):
    """Appends a planar polygon of any number of vertices to the model,
       applying faceted shading automatically."""
    if n is None: n = len(v)*[None]
    for t in range(len(v)-2):
      self.FacetedTriangle( vv=[v[0], v[t+1], v[t+2]],
                            nn=[n[0], n[t+1], n[t+2]], color=color)
   # --- This is an attempt to fill in areas that are concave without
   # --- having triangles sticking outside. It doesn't work yet.
   ## --- Get norm of first triangle
   #n1 = self.Norm( self.Cross(v[1]-v[0], v[2]-v[0]) )
   #starti = 0
   #startlen = len(v)
   #while 1:
   #  tilist = []
   #  badlist = []
   #  ilist = range(len(v))
   #  i = starti
   #  normsign = +1
   #  while len(ilist) > 2:
   #    n = len(ilist)
   #    i0 = ilist[(i  )%n]
   #    i1 = ilist[(i+1)%n]
   #    i2 = ilist[(i+2)%n]
   #    if [i0,i1,i2] in badlist:
   #      break
   #    ni = self.Norm( self.Cross(v[i1]-v[i0], v[i2]-v[i0]) )
   #    if self.Dot(n1,ni)*normsign >= 0.:
   #      tilist.append([i0,i1,i2])
   #      del ilist[(i+1)%n]
   #    else:
   #      badlist.append([i0,i1,i2])
   #      i += 1
   #  if len(ilist) < 3:
   #    break
   #  print starti,startlen,ilist
   #  starti += 1
   #  if starti == startlen:
   #    if normsign == -1: break
   #    starti = 0
   #    normsign = -1
   #for ti in tilist:
   #  self.FacetedTriangle( vv=[v[ti[0]], v[ti[1]], v[ti[2]]],
   #                        nn=[n[ti[0]], n[ti[1]], n[ti[2]]], color=color)
     

  def DoSmoothShading(self):
    rsq = sum(self.triangles**2,1)
    ii = argsort(rsq)
    # to be completed later

  def DoSmoothShading1(self):
    """Change a faceted model to smooth shaded, by averaging normals at
    coinciding vertices.
    
    This is a very slow and simple smooth shading
    implementation which has to figure out the connectivity of the
    model and does not attempt to detect sharp edges.

    It attempts to work even in two-sided mode where there are two
    opposite normals at each vertex.  It may fail somehow in pathological
    cases. """

    pos = self.model.pos
    normal = self.model.normal

    vertex_map = {}  # vertex position -> vertex normal
    vertex_map_backface = {}
    for i in range( len(pos) ):
      tp = tuple(pos[i])
      old_normal = vertex_map.get( tp, (0,0,0) )
      if dot(old_normal, normal[i]) >= 0:
        vertex_map[tp] = normal[i] + old_normal
      else:
        vertex_map_backface[tp] = normal[i] + vertex_map_backface.get(tp, (0,0,0))

    for i in range( len(pos) ):
      tp = tuple(pos[i])
      if dot(vertex_map[tp], normal[i]) >= 0:
        normal[i] = vertex_map[tp] and self.Norm( vertex_map[ tp ] )
      else:
        normal[i] = vertex_map_backface[tp] and self.Norm(vertex_map_backface[tp] )

  def DrawNormal(self, scale):
    pos = self.model.pos
    normal = self.model.normal
    for i in range(len(pos)):
      arrow(pos=pos[i], axis=normal[i]*scale)

  def Cross(self,v1,v2):
    return array([v1[1]*v2[2] - v1[2]*v2[1],
                  v1[2]*v2[0] - v1[0]*v2[2],
                  v1[0]*v2[1] - v1[1]*v2[0]])

  def Norm(self,v):
    magv = sqrt(sum(v**2))
    if magv != 0.: return v/magv
    else:          return v

  def Dot(self,v1,v2):
    return sum(array(v1)*array(v2))
########################################################################
class VisualMesh (VisualModel):
  """
xvalues, yvalues, zvalues: 2-D arrays containing the coordinates and data
twoSided=1: when true, surface is two sided
normalsign=1: when negative, show the under side (and twoSided == 0)
color=None: can be specified as an [r,g,b] list
scene=None: an already existing display scene. When None, create a new one.
title='Mesh': Display title - only used when new scene created.
  """
  def __init__(self, xvalues=None, yvalues=None, zvalues=None,
               xscaled=0,zscaled=1,
               twoSided=1,normalsign=1,color=None,color1=None,color2=None,
               scene=None,title=None,vrange=None,viewer=None):
    if not title: title = 'Mesh'
    VisualModel.__init__(self,twoSided=twoSided,normalsign=normalsign,
                              scene=scene,title=title,
                              vrange=vrange,viewer=viewer,display=1)

    assert zvalues is not None,"zvalues must be specified"

    s = shape(zvalues)
    if len(s) != 2:
      print 'First argument must be a 2-Dimensional array'
      return
    if xvalues is None:
      xvalues = arange(s[0])[:,NewAxis]*ones(s[1],'d')
    elif len(shape(xvalues))==1:
      xvalues = xvalues[:,NewAxis]*ones(s[1],'d')
    if yvalues is None:
      yvalues = arange(s[1])*ones(s[0],'d')[:,NewAxis]
    elif len(shape(yvalues))==1:
      yvalues = yvalues*ones(s[0],'d')[:,NewAxis]

    if xscaled:
      xrange = maxnd(xvalues) - minnd(xvalues)
      xvalues = xvalues/xrange
      yvalues = yvalues/xrange
    if zscaled:
      xrange = maxnd(xvalues) - minnd(xvalues)
      zrange = maxnd(zvalues) - minnd(zvalues)
      zvalues = zvalues/zrange*xrange/2.

    if color1 is not None:
      getcolor = 1
      zmin = minnd(zvalues)
      zmax = maxnd(zvalues)
      if color2 is None: color2 = zeros(3,'d')
    else:
      getcolor = 0

    points = zeros( xvalues.shape + (3,), Float )
    points[...,0] = xvalues
    points[...,1] = yvalues
    points[...,2] = zvalues

    for i in range(zvalues.shape[0]-1):
      for j in range(zvalues.shape[1]-1):
        if getcolor:
          color = color1 + color2*(points[i,j,2] - zmin)/(zmax - zmin)
        self.FacetedPolygon([points[i,j], points[i,j+1],
                             points[i+1,j+1], points[i+1,j]],
                            color=color)

    self.CreateDXObject()
    if display: self.Display()

########################################################################
class VisualRevolution(VisualModel):
  """
Visualize surface of revolution
 - srfrv: optional function defining surface
 - zzmin: min z extent of surface
 - zzmax: max z extent of surface
 - rendzmin: radius at zmin
 - rendzmax: radius at zmax
 - nz=20: if srfrv is given, number of z points radius sampled at
 - nth=20: number of theta angles points sampled at
 - phimin=0.: miminmum phi angle
 - phimax=2*pi: maximum phi angle
 - fillinends=1: if theta range < 2pi, fill in ends if true.
 - xoff=0: x offset
 - yoff=0: y offset
 - zoff=0: z offset
 - rofzdata=None: optional tablized radius defining surface
 - zdata=None: optional tablized z poins defining surface
 - raddata=None: optional tablized circle radius defining surface
 - zcdata=None: optional tablized circle z center defining surface
 - rcdata=None: optional tablized circle r center defining surface
 - ntpts=5: number of points sampled for circles
 - twoSided=1: if true, include both sides of the surface
 - normalsign=1: 1 when data is clockwise, -1 with counterclockwise
 - color=None: RGB color for surface, of form [r,g,b]
 - scene=None: include object in exising scene (only for VPython)
 - title=None: window title
 - vrange=None: set view range (only for VPython)
 - viewer=None: select viewer, either 'OpenDX' or 'VPython'
 - display=1: if 1, immeidately display object
  """
  def __init__(self,srfrv=None,zzmin=None,zzmax=None,
                    rendzmin=None,rendzmax=None,
                    nz=20,nth=20,phimin=None,phimax=None,fillinends=0,
                    xoff=0,yoff=0,zoff=0,
                    rofzdata=None,zdata=None,raddata=None,
                    zcdata=None,rcdata=None,ntpts=5,
                    twoSided=0,normalsign=1,color=None,
                    scene=None,title=None,vrange=None,
                    viewer=None,display=1,kwdict={}):
    kwdict.update(kwdict.setdefault('kwdict',{}))
    if 'kwdict' in kwdict: del kwdict['kwdict']
    for arg in kwdict.keys(): exec(arg+" = kwdict['"+arg+"']")
    if not title: title = 'Surface of revolution'
    VisualModel.__init__(self,twoSided=twoSided,normalsign=1,
                              scene=scene,title=title,
                              vrange=vrange,viewer=viewer)

    # --- If phimin and phimax are not given, then the ends never need
    # --- filling in
    if phimin is None and phimax is None: fillinends = 0
    if phimin is None: phimin = 0.
    if phimax is None: phimax = 2*pi

    if rofzdata is None and zdata is None:
      # --- Include extra points for the z end of the surface
      zz = arange(-1,nz+2)*(zzmax - zzmin)/nz + zzmin
      rr = ones(nz+3,'d')
      for i in range(nz+1):
        warp.f3d.srfrv_z = zz[i+1]
        srfrv()
        rr[i+1] = warp.f3d.srfrv_r
      zzleft = zz[:-1]
      zzrght = zz[1:]
      rrleft = rr[:-1]
      rrrght = rr[1:]
      ttleft = len(zzleft)*[None]
      ttrght = len(zzrght)*[None]
    else:
      zzleft = [0.]
      zzrght = [0.]
      rrleft = [0.]
      rrrght = [0.]
      ttleft = [0.]
      ttrght = [0.]
      for i in range(len(zdata)-1):
        z,r,rad,zc,rc = zdata[i],rofzdata[i],raddata[i],zcdata[i],rcdata[i]
        zp1,rp1 = zdata[i+1],rofzdata[i+1]
        if rad == largepos:
          zzleft.append(z)
          zzrght.append(zp1)
          rrleft.append(r)
          rrrght.append(rp1)
          ttleft.append(arctan2((zp1 - z),(r - rp1)))
          ttrght.append(ttleft[i])
        else:
          t1 = arctan2(r-rc,z-zc)
          t2 = arctan2(rp1-rc,zp1-zc)
          if t1 > t2 and rad < 0.: t2 += 2*pi
          if t1 < t2 and rad > 0.: t1 += 2*pi
          tt = t1 + (t2 - t1)*arange(ntpts+1)/ntpts
          zz = zc + abs(rad)*cos(tt)
          rr = rc + abs(rad)*sin(tt)
          zzleft += list(zz[:-1])
          zzrght += list(zz[1:])
          rrleft += list(rr[:-1])
          rrrght += list(rr[1:])
          if rad < 0.: addpi = pi
          else:        addpi = 0.
          ttleft += list(tt[:-1]+addpi)
          ttrght += list(tt[1:]+addpi)

      zzleft += [0.]
      zzrght += [0.]
      rrleft += [0.]
      rrrght += [0.]
      ttleft += [0.]
      ttrght += [0.]

    # --- Now set end points
    if rendzmin == largepos: rendzmin = 2.*max(max(rrleft),max(rrrght))
    if rendzmax == largepos: rendzmax = 2.*max(max(rrleft),max(rrrght))
    if rendzmin is not None:
      zzleft[0] = zzmin
      zzrght[0] = zzleft[1]
      rrleft[0] = rendzmin
      rrrght[0] = rrleft[1]
      ttleft[0] = pi *(normalsign > 0.)
      ttrght[0] = pi *(normalsign > 0.)
    if rendzmax is not None:
      zzleft[-1] = zzrght[-2]
      zzrght[-1] = zzmax
      rrleft[-1] = rrrght[-2]
      rrrght[-1] = rendzmax
      ttleft[-1] = 0. + pi*(normalsign < 0.)
      ttrght[-1] = 0. + pi*(normalsign < 0.)

    phi = phimin + (phimax - phimin)*arange(0,nth+1)/nth
    xx = cos(phi)
    yy = sin(phi)

    for i in range(len(zzleft)):
      for j in range(len(xx)-1):
        p1 = array([rrleft[i]*xx[j  ]+xoff, rrleft[i]*yy[j  ]+yoff,
                    zzleft[i] + zoff])
        p2 = array([rrleft[i]*xx[j+1]+xoff, rrleft[i]*yy[j+1]+yoff,
                    zzleft[i] + zoff])
        p3 = array([rrrght[i]*xx[j+1]+xoff, rrrght[i]*yy[j+1]+yoff,
                    zzrght[i] + zoff])
        p4 = array([rrrght[i]*xx[j  ]+xoff, rrrght[i]*yy[j  ]+yoff,
                    zzrght[i] + zoff])
        if ttleft[i] is not None:
          n1 = array([sin(ttleft[i])*cos(phi[j  ]),
                      sin(ttleft[i])*sin(phi[j  ]),cos(ttleft[i])])
          n2 = array([sin(ttleft[i])*cos(phi[j+1]),
                      sin(ttleft[i])*sin(phi[j+1]),cos(ttleft[i])])
          n3 = array([sin(ttrght[i])*cos(phi[j+1]),
                      sin(ttrght[i])*sin(phi[j+1]),cos(ttrght[i])])
          n4 = array([sin(ttrght[i])*cos(phi[j  ]),
                      sin(ttrght[i])*sin(phi[j  ]),cos(ttrght[i])])
        else:
          n1,n2,n3,n4 = None,None,None,None
        self.FacetedPolygon([p1,p2,p3,p4],[n1,n2,n3,n4],color=color)

    if fillinends:
      rr = array([rrleft[0]] + list(rrrght))
      zz = array([zzleft[0]] + list(zzrght)) + zoff
      n = len(rr)
      if phimin < phimax: nsign = -normalsign
      else:               nsign = +normalsign
      xx,yy = rr*cos(phimin)+xoff,rr*sin(phimin)+yoff
      pp = map(array,zip(xx,yy,zz))
      nx,ny = -nsign*sin(phimin),nsign*cos(phimin)
      nn = map(array,zip(n*[nx],n*[ny],n*[0.]))
      self.FacetedPolygon(pp,nn,color=color)
      xx,yy = rr*cos(phimax)+xoff,rr*sin(phimax)+yoff
      pp = map(array,zip(xx,yy,zz))
      nx,ny = -nsign*sin(phimax),nsign*cos(phimax)
      nn = map(array,zip(n*[nx],n*[ny],n*[0.]))
      self.FacetedPolygon(pp,nn,color=color)

    self.CreateDXObject()
    if display: self.Display()
