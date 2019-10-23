"""
Secondaries: class for generating secondaries
"""
from scipy.constants import m_e as e_mass
from scipy.constants import e as qe

from ..warp import *
from ..utils.appendablearray import AppendableArray

try:
    import pos
    from pos import pos as posC
    try:
        pos_version = pos.__version__
        pos.init_physical_constants() # needed for pos_version >= '17p3'
    except AttributeError:
        pos_version = '15p3'
except ImportError:
    pos = None

try:
    from txphysics import txigenelec, txstopping, txrand
    l_txphysics = 1
except:
    print 'WARNING: module txphysics is not accessible.'
    l_txphysics = 0

try:
    import desorb
    l_desorb = 1
except:
    print 'WARNING: module desorb is not accessible.'
    l_desorb = 0
import time

def coldotprod(a, b):
    return np.sum(a*b, axis=0)

def secondariesdoc():
    from ..particles import Secondaries
    print Secondaries.__doc__

class Secondaries:
    """Class for generating secondaries
     - isinc:      list of incident species
     - conductors: list of list of conductors
     - issec:      list of list of secondary species
     - type:       list of list of types of interaction
                   e- -> Cu              = 1
                   e- -> Stainless Steel = 2
                   H  -> Au              = 3
                   He -> Au              = 4
                   K  -> SS              = 5
     - min_age: this sets a minimum age for a macroparticle in order to
                qualify for emitting secondaries (must be old enough to procreate :-) ).
                This is in units of time steps. The default is None (no minimum age required).
     - set_params_user: optional function that sets SEY parameters provided by the user.
                        The default routine is always called, so this routine only
                        needs to set parameters which differ.
     - vmode: 1 (default) -> uses instantaneous velocity to compute angle to normal of conductor
              2           -> uses difference between current and old positions to compute angle to normal of conductor
     - l_verbose: sets verbosity (default=0).
     - lcallscrapercontrollers=false: After new particles are added, the
               particleboundaries3d routine is called. Normally, this would
               scrape particles on conductors, but since this can be slow and
               should not be necessary (since new particles are explicitly
               created outside of conductors), this call is skipped by telling
               particleboundaries3d not to call the controllers. If there are
               other conductors that the new particles could be lost in, then
               this option needs to be set to true for that scraping to happen.
     - luseoldpositionifout=False:
           If the new particle position is outside of the grid, use the old position of the lost particle
           (its position from the previous time step).
    """
    def __init__(self,isinc=None,conductors=None,issec=None,set_params_user=None,material=None,
                      xoldpid=None,yoldpid=None,zoldpid=None,min_age=None,vmode=1,l_verbose=0,
                      l_set_params_user_only=0,lcallscrapercontrollers=0,l_trackssnparents=0,
                      l_usenew=1,
                      luseoldpositionifout=False,
                      maxsec=None, pyecloud_secemi_object=None, pyecloud_fact_clean=None,
                      pyecloud_fact_split=None, pyecloud_nel_mp_ref=None):

        if pyecloud_secemi_object is not None:
            assert(l_set_params_user_only == 0)
            assert(l_usenew == 1)
            assert(luseoldpositionifout == False)
            assert(min_age is None) # Could be removed (but please test!)
            assert(set_params_user is None)
            assert(material is None)
            assert(maxsec is None) # Probably just ignored
           
            assert(pyecloud_fact_clean is not None)
            assert(pyecloud_fact_split is not None)
            assert(pyecloud_nel_mp_ref is not None)

            self.flag_pyecloud = True
            self.pyeclsecemi = pyecloud_secemi_object
            self.pyecloud_fact_clean = pyecloud_fact_clean
            self.pyecloud_fact_split = pyecloud_fact_split
            self.pyecloud_nel_mp_ref = pyecloud_nel_mp_ref
            
        else:
            self.flag_pyecloud = False


        self.totalcount = 0
        self.totallost = 0
        top.lresetlostpart=true
        self.inter={}
        self.outparts=[]
        # self.isinc = isinc
        # self.conductors = conductors
        # self.issec = issec
        # self.type = type
        self.vmode=vmode
        self.l_verbose=l_verbose
        self.l_record_timing=1
        if self.l_record_timing:self.timings=[]
        self.lcallscrapercontrollers = lcallscrapercontrollers
        self.l_trackssnparents=l_trackssnparents
        # self.condids={}
        # self.emitted={}
        if not self.flag_pyecloud:
            self.set_params_user=set_params_user
        self.l_set_params_user_only=l_set_params_user_only
        self.mat_number=1
        if not self.flag_pyecloud:
            if maxsec is None:
                if pos is None:
                    maxsec = 10
                else:
                    maxsec = posC.maxsec
            else:
                if pos is not None:
                    assert maxsec == posC.maxsec, Exception('maxsec must be the same as posC.maxsec')
        self.maxsec = maxsec
        if not self.flag_pyecloud:
            self.call_set_params_user(maxsec,self.mat_number)
        self.min_age=min_age
        if self.min_age is not None:
            w3d.l_inj_rec_inittime=true
            if top.tbirthpid==0:
                top.tbirthpid=nextpid()
                setuppgroup(top.pgroup)
        self.l_usenew=l_usenew
        self.luseoldpositionifout = luseoldpositionifout
        self.install()
        if xoldpid is None:
            self.xoldpid=top.npid-3
        else:
            self.xoldpid=xoldpid
        if yoldpid is None:
            self.yoldpid=top.npid-2
        else:
            self.yoldpid=yoldpid
        if zoldpid is None:
            self.zoldpid=top.npid-1
        else:
            self.zoldpid=zoldpid
        if self.l_trackssnparents:
            if top.ssnpid==0:
                top.ssnpid=nextpid()
                setuppgroup(top.pgroup)
            if top.ssnparentpid==0:
                top.ssnparentpid=nextpid()
                setuppgroup(top.pgroup)
        # --- set variables for secondary electrons routines
        self.secelec_ns = zeros(1,'l')
        self.secelec_un = zeros(maxsec,'d')
        self.secelec_ut = zeros(maxsec,'d')
        self.secelec_uz = zeros(maxsec,'d')
        self.secelec_ityps  = zeros(maxsec,'l')
        self.secelec_ekstot = zeros(maxsec,'d')
        self.secelec_dele   = zeros(1,'d')
        self.secelec_delr   = zeros(1,'d')
        self.secelec_delts  = zeros(1,'d')
        # --- set arrays for emitted particles
        self.npmax={}
        self.nps={}
        self.x={}
        self.y={}
        self.z={}
        self.ux={}
        self.uy={}
        self.uz={}
        self.pid={}
        # --- set history list
        self.htime=AppendableArray(typecode='d')
        self.ek0av=AppendableArray(typecode='d')     # average collision kinetic energy [eV]
        self.ek0max=AppendableArray(typecode='d')    # maximum collision kinetic energy [eV]
        self.costhav=AppendableArray(typecode='d')   # average collision angle
        self.pxsum=AppendableArray(typecode='d')     # x momentum deposited per step
        self.pysum=AppendableArray(typecode='d')     # y momentum deposited per step
        self.pzsum=AppendableArray(typecode='d')     # z momentum deposited per step
        self.power_dep=AppendableArray(typecode='d') # instantaneous power deposition [W]
        self.power_emit=AppendableArray(typecode='d') # instantaneous power emission [W]
        self.power_diff=AppendableArray(typecode='d') # instantaneous power deposition [W]
        self.piditype = 0

        if not self.flag_pyecloud:
            if pos is not None:
                if posC.nsteps_g != 0:
                    self.piditype = nextpid()

        self.lrecursivegenerate = 0

        if isinc is None:return
        for iis,js in enumerate(isinc):
            for ics,cond in enumerate(conductors[iis]):
                if material is None: m = None
                else:                m = material[iis][ics]
                self.add(js,cond,issec[iis][ics],m)

    def __getstate__(self):
        dict = self.__dict__.copy()

        # --- Functions cannot be pickled, so set_params_user has
        # --- to be specially handled. If it is defined, then only save
        # --- the name of the function so that it can be restored later.
        # --- This is dealt with in call_set_params_user so nothing needs
        # --- to happen in setstate. It is done there since this instance
        # --- may be restore before the function.
        if self.set_params_user is not None:
            try:
                dict['set_params_user'] = self.set_params_user.__name__
            except AttrbiuteError:
                print ("Warning: Secondaries set_params_user function '%s' is not "+
                       "a proper function and will not be saved")%self.set_params_user
                dict['set_params_user'] = None

        # --- These arrays are swig arrays and cannot be pickled.
        # --- They are only temporary arrays and don't need to be
        # --- restore so they can just be deleted. They are recreated
        # --- in the setstate below.
        try:
            del dict['emitted_e']
            del dict['emitted_bn']
            del dict['emitted_bt']
            del dict['emitted_bz']
        except KeyError:
            pass

        return dict

    def __setstate__(self,dict):
        'Restore the dictionary'
        self.__dict__.update(dict)

    def add(self,incident_species=None,conductor=None,emitted_species=None,material=None,interaction_type=None,
                 scale_factor=None,scale_factor_velocity=1.,forced_yield=None,init_position_offset=0.,
                 specular_fraction=None,conductor_temperature=None):
        if material is None: material = conductor.material
        if interaction_type is None and incident_species.type is Electron:
            if material=='Cu' and incident_species.type is Electron: interaction_type=1
            if material=='SS' and incident_species.type is Electron: interaction_type=2
            # if material=='Au' and incident_species.type is Hydrogen: interaction_type=3
            # if material=='Au' and incident_species.type is Helium:   interaction_type=4
            # if material=='SS' and incident_species.type is Potassium:interaction_type=5
            if interaction_type is None:raise Exception('Error in Secondaries: invalid material or incident species')
        isinc=incident_species
        if not isinstance(emitted_species,list):emitted_species=[emitted_species]
        issec=[]
        for e in emitted_species:
            issec.append(e.jslist[0])
        if isinc not in self.inter:
            self.inter[isinc]={}
            for key in ['emitted','absorbed','condids','issec','conductors','type','isneut','incident_species', \
                        'emitted_species','material','scale_factor','scale_factor_velocity','forced_yield', \
                        'init_position_offset','specular_fraction','conductor_temperature']:
                self.inter[isinc][key]=[]
            self.inter[isinc]['incident_species']=incident_species
        self.inter[isinc]['condids']               += [conductor.condid]
        self.inter[isinc]['conductors']            += [conductor]
        self.inter[isinc]['emitted']               += [[]]
        self.inter[isinc]['absorbed']               += [[]]
        for i in range(len(emitted_species)):
            self.inter[isinc]['emitted'][-1]         += [0.]
            self.inter[isinc]['absorbed'][-1]        += [0.]
        self.inter[isinc]['issec']                 += [issec]
        self.inter[isinc]['emitted_species']       += [emitted_species]
        self.inter[isinc]['material']              += [material]
        self.inter[isinc]['type']                  += [interaction_type]
        self.inter[isinc]['scale_factor']          += [scale_factor]
        self.inter[isinc]['scale_factor_velocity'] += [scale_factor_velocity]
        self.inter[isinc]['forced_yield']          += [forced_yield]
        self.inter[isinc]['init_position_offset']  += [init_position_offset]
        self.inter[isinc]['specular_fraction']     += [specular_fraction]
        self.inter[isinc]['conductor_temperature'] += [conductor_temperature]
        for e in emitted_species:
            js=e.jslist[0]
            if js not in self.x:
                self.nps[js]=0
                self.npmax[js]=4096
                self.allocate_temps(js)

    def allocate_temps(self,js):
        self.x[js]=fzeros(self.npmax[js],'d')
        self.y[js]=fzeros(self.npmax[js],'d')
        self.z[js]=fzeros(self.npmax[js],'d')
        self.ux[js]=fzeros(self.npmax[js],'d')
        self.uy[js]=fzeros(self.npmax[js],'d')
        self.uz[js]=fzeros(self.npmax[js],'d')
        if top.npid>0 :
            self.pid[js]=fzeros([self.npmax[js],top.npid],'d')

    def install(self):
        # --- The secondaries should be created immediately after particles
        # --- are scraped.
        if self.l_usenew:
            if not isinstalledafterscraper(self.generatenew):
                installafterscraper(self.generatenew)
        else:
            raise ValueError('Not supported anymore!')

    def addpart(self,nn,x,y,z,ux,uy,uz,js,weight=None,itype=None,ssnparent=None):
        if top.npid > self.pid[js].shape[1]:
            if self.nps[js] > 0:
                raise Exception('top.npid changed so existing particle pid data is corrupted')
            self.allocate_temps(js)
        if self.nps[js]+nn>self.npmax[js]:self.flushpart(js)
        if self.nps[js]+nn>self.npmax[js]:
            self.npmax[js] = nint(nn*1.2)
            self.allocate_temps(js)
        il=self.nps[js]
        iu=il+nn
        xx=fones(nn,'d')
        self.x[js][il:iu]=x*xx
        self.y[js][il:iu]=y*xx
        self.z[js][il:iu]=z*xx
        self.ux[js][il:iu]=ux
        self.uy[js][il:iu]=uy
        self.uz[js][il:iu]=uz
        if weight is not None:self.pid[js][il:iu,top.wpid-1]=weight
        if itype is not None:self.pid[js][il:iu,self.piditype-1]=itype.astype(float64)
        if ssnparent is not None:self.pid[js][il:iu,top.ssnparentpid-1]=ssnparent
        self.nps[js]+=nn

    def addparticles(self,nn,x,y,z,ux,uy,uz,js,weight=None,itype=None,ssnparent=None):
        if weight is not None or itype is not None or ssnparent is not None:
            pid = zeros((nn,top.npid))
            if weight is not None:pid[:,top.wpid-1]=weight
            if itype is not None:pid[:,self.piditype-1]=itype.astype(float64)
            if ssnparent is not None:pid[:,top.ssnparentpid-1]=ssnparent
        else:
            pid=0.
        gi=1./sqrt(1.+(ux**2+uy**2+uz**2)/clight**2)
        addparticles(x=x,
                     y=y,
                     z=z,
                     vx=ux,
                     vy=uy,
                     vz=uz,
                     gi=gi,
                     pid=pid,
                     js=js,
                     lmomentum=true,
                     lallindomain=true)

    def flushpart(self,js):
        if self.nps[js]>0:
            nn=self.nps[js]
            self.totalcount += nn
            if top.npid > 0:
                # --- Note that this includes the weights (wpid)
                pid = self.pid[js][:nn,:]
            else:
                pid = 0.
            ux=self.ux[js][:nn]
            uy=self.uy[js][:nn]
            uz=self.uz[js][:nn]
            gi=1./sqrt(1.+(ux**2+uy**2+uz**2)/clight**2)
            # --- get velocity in boosted frame if using a boosted frame of reference
            if top.boost_gamma>1.:
                uzboost = clight*sqrt(top.boost_gamma**2-1.)
                setu_in_uzboosted_frame3d(nn,ux,uy,uz,gi,
                                          uzboost,
                                          top.boost_gamma)
            addparticles(x=self.x[js][:nn],
                         y=self.y[js][:nn],
                         z=self.z[js][:nn],
                         vx=ux,
                         vy=uy,
                         vz=uz,
                         gi=gi,
                         pid=pid,
                         js=js,
                         lmomentum=true,
                         lallindomain=true)
            self.nps[js]=0

    def printall(self,l_cgm=0):
        title='        *** Particle/wall interactions ***\n'
        textblock=''
        swidth=0
        cwidth=0
        ewidth={}
        for js in self.inter:
            swidth=max(swidth,len(self.inter[js]['incident_species'].name))
            for ics,cond in enumerate(self.inter[js]['conductors']):
                cwidth=max(cwidth,len(cond.name))
                for ie,emitted_species in enumerate(self.inter[js]['emitted_species'][ics]):
                    if ie not in ewidth:
                        ewidth[ie]=len(emitted_species.name)
                    else:
                        ewidth[ie]=max(ewidth[ie],len(emitted_species.name))
        fs='%%-%gs'%swidth
        fc='%%-%gs'%cwidth
        fe={}
        for ie in ewidth:
            fe[ie]='%%-%gs'%ewidth[ie]
        for js in self.inter:
            sname=fs%self.inter[js]['incident_species'].name
            textblock+='\n'
            for ics,cond in enumerate(self.inter[js]['conductors']):
                cname=fc%cond.name
                for ie,emitted_species in enumerate(self.inter[js]['emitted_species'][ics]):
                    if ie==0:
                        ename=fe[ie]%emitted_species.name
                    else:
                        ename+=' + '+fe[ie]%emitted_species.name
                textblock += sname+' + '+cname+' => '+ename+'\n'
        if l_cgm:
            plt(title,0.22,0.905,justify="LT",height=14)
            plt(textblock,0.14,0.88,justify="LT",height=9,font='courier')
            fma()
        else:
            print title
            print textblock

    def generatenew(self,local=0,l_accumulate_hist=1):
        # --- theta and phi are angles from normal to surface with regard to z and x axis respectively
        # --- psi is angle to rotate warp local frame to Posinst local frame around normal
        # --- eta is angle between incident velocity vector and normal to surface (called theta in Posinst)

        if self.lrecursivegenerate: return

        if self.l_verbose>1:print 'start secondaries generation'

        if self.l_record_timing:t1 = time.clock()

        # --- reset 'emitted' list to zero
        for js in self.inter:
            for i in range(len(self.inter[js]['emitted'])):
                for j in range(len(self.inter[js]['emitted'][i])):
                    self.inter[js]['emitted'][i][j] = 0.
                    self.inter[js]['absorbed'][i][j] = 0.

        # --- set computing box mins and maxs
        if w3d.l4symtry:
            xmin=-w3d.xmmax
        else:
            xmin=w3d.xmmin
        xmax=w3d.xmmax
        if w3d.l2symtry or w3d.l4symtry:
            ymin=-w3d.ymmax
        else:
            ymin=w3d.ymmin
        ymax=w3d.ymmax
        zmin=w3d.zmmin+top.zgrid
        zmax=w3d.zmmax+top.zgrid

        # --- initializes history quantities
        weighttot=0.
        ek0sum=0.
        costhsum=0.
        ek0max=0.
        ek0emitsum=0.
        pxsum = 0.
        pysum = 0.
        pzsum = 0.
        if self.l_record_timing:t2 = time.clock()
        tinit=tgen=tprepadd=tadd=0.
        # compute number of secondaries and create them
        for ints in self.inter:
            incident_species=self.inter[ints]['incident_species']
            for js in incident_species.jslist:
                if self.l_verbose:print 'js',js
                if top.npslost[js]==0:continue
                # if top.npslost[js]==0 or top.it%top.pgroup.ndts[js] != 0:continue
                stride=top.pgroup.ndts[js]
                i1 = top.inslost[js] - 1
                i2 = top.inslost[js] + top.npslost[js] - 1
                for ics,cond in enumerate(self.inter[incident_species]['conductors']):
                    icond = cond.condid
                    if self.l_verbose:print 'ics',ics
                    iit = compress(top.pidlost[i1+top.it%stride:i2:stride,-1]==icond,arange(top.it%stride,top.npslost[js],stride))
                    n = len(iit)
                    if self.l_verbose:print 'nlost=',n
                    if n==0:continue
                    xplost = take(top.xplost[i1:i2],iit)
                    yplost = take(top.yplost[i1:i2],iit)
                    zplost = take(top.zplost[i1:i2],iit)
                    # --- exclude particles out of computational box
                    # if w3d.solvergeom==w3d.RZgeom:
                    #   condition = (sqrt(xplost**2+yplost**2)<xmax) & \
                    #               (zplost>zmin) & (zplost<zmax)
                    # else:
                    #   condition = (xplost>xmin) & (xplost<xmax) & \
                    #               (yplost>ymin) & (yplost<ymax) & \
                    #               (zplost>zmin) & (zplost<zmax)
                    # --- exclude particles recently created
                    if self.min_age is not None:
                        inittime = take(top.pidlost[i1:i2,top.tbirthpid-1],iit,0)
                        condition =  ((top.time-inittime)>self.min_age*top.dt)
                        # condition = condition & ((top.time-inittime)>self.min_age*top.dt)
                        iit2 = compress(condition,arange(n))
                    else:
                        iit2=arange(n)
                    n = len(iit2)
                    if self.l_verbose:print 'nlost=',n
                    if n==0:continue
                    self.totallost += n
                    xplost = take(xplost,iit2)
                    yplost = take(yplost,iit2)
                    zplost = take(zplost,iit2)
                    iit    = take(iit,iit2)
                    if self.l_record_timing:tstart=wtime()
                    uxplost = take(top.uxplost[i1:i2],iit).copy()
                    uyplost = take(top.uyplost[i1:i2],iit).copy()
                    uzplost = take(top.uzplost[i1:i2],iit).copy()
                    gaminvlost = take(top.gaminvlost[i1:i2],iit).copy()
                    # --- get velocity in lab frame if using a boosted frame of reference
                    if 0:#top.boost_gamma>1.:
                        uzboost = clight*sqrt(top.boost_gamma**2-1.)
                        setu_in_uzboosted_frame3d(n,uxplost,uyplost,uzplost,gaminvlost,
                                                  -uzboost,
                                                  top.boost_gamma)
                    if self.l_trackssnparents: ssnplost = take(top.pidlost[i1:i2,top.ssnpid-1],iit)
                    if self.vmode == 2 or self.luseoldpositionifout:
                        xplostold = take(top.pidlost[i1:i2,self.xoldpid],iit,0)
                        yplostold = take(top.pidlost[i1:i2,self.yoldpid],iit,0)
                        zplostold = take(top.pidlost[i1:i2,self.zoldpid],iit,0)
                    if self.vmode==1:
                        vxplost = uxplost*gaminvlost
                        vyplost = uyplost*gaminvlost
                        vzplost = uzplost*gaminvlost
                    elif self.vmode==2:
                        vxplost = (xplost - xplostold)/top.dt
                        vyplost = (yplost - yplostold)/top.dt
                        vzplost = (zplost - zplostold)/top.dt
                    else:
                        raise Exception('Error in Secondaries, one should have lmode=1 or 2, but have lmode=%g'%self.lmode)
                    # set energy of incident particle in eV
                    e0 = where(gaminvlost==1., \
                               0.5*top.pgroup.sm[js]*(uxplost**2+uyplost**2+uzplost**2)/jperev,
                               (1./gaminvlost-1.)*top.pgroup.sm[js]*clight**2/jperev)
                    if self.l_verbose:
                        print 'xplost',xplost
                        print 'yplost',yplost
                        print 'zplost',zplost
                        print 'e0',e0,gaminvlost,uxplost,uyplost,uzplost
                    v = array([vxplost,vyplost,vzplost])
                    # u = array([uxplost,uyplost,uzplost])
                    theta = take(top.pidlost[i1:i2,-3],iit,0)
                    phi   = take(top.pidlost[i1:i2,-2],iit,0)
                    if top.wpid>0:
                        weight = take(top.pidlost[i1:i2,top.wpid-1],iit,0)
                    else:
                        weight=1.
                    costheta = cos(theta)
                    sintheta = sin(theta)
                    cosphi   = cos(phi)
                    sinphi   = sin(phi)
                    # --- theta is relative to the z axis, phi to the x axis in the x-y plane
                    n_unit0 = array([sintheta*cosphi,sintheta*sinphi,costheta])
                    coseta = -sum(v*n_unit0,axis=0)/sqrt(sum(v*v,axis=0))
                    if top.wpid==0:
                        ek0sum+=sum(e0)*top.pgroup.sw[js]
                        costhsum+=sum(abs(coseta))*top.pgroup.sw[js]
                        weighttot+=n*top.pgroup.sw[js]
                        pxsum += sum(uxplost)*top.pgroup.sm[js]*top.pgroup.sw[js]
                        pysum += sum(uyplost)*top.pgroup.sm[js]*top.pgroup.sw[js]
                        pzsum += sum(uzplost)*top.pgroup.sm[js]*top.pgroup.sw[js]
                    else:
                        ek0sum+=sum(weight*e0)*top.pgroup.sw[js]
                        costhsum+=sum(weight*abs(coseta))*top.pgroup.sw[js]
                        weighttot+=sum(weight)*top.pgroup.sw[js]
                        pxsum += sum(weight*uxplost)*top.pgroup.sm[js]*top.pgroup.sw[js]
                        pysum += sum(weight*uyplost)*top.pgroup.sm[js]*top.pgroup.sw[js]
                        pzsum += sum(weight*uzplost)*top.pgroup.sm[js]*top.pgroup.sw[js]
                    ek0max=max(max(e0),ek0max)
                    
                    ## Collection of data
                    if 1:#cond.lcollectlpdata:
                        if js not in cond.lostparticles_angles:
                            cond.lostparticles_angles[js]=zeros(181,'d')
                        if js not in cond.lostparticles_energies:
                            cond.lostparticles_energies[js]=zeros(1001,'d')
                        e0min = min(e0)
                        e0max = max(e0)
                        # e0min=0.
                        # e0max=1.e6
                        if js not in cond.lostparticles_minenergy:
                            cond.lostparticles_minenergy[js]=e0min
                        if js not in cond.lostparticles_maxenergy:
                            cond.lostparticles_maxenergy[js]=e0max
                        l_rescale_energy_array=0
                        if e0min<cond.lostparticles_minenergy[js]:
                            e0minnew=0.9*e0min
                            l_rescale_energy_array=1
                        else:
                            e0minnew=cond.lostparticles_minenergy[js]
                        if e0max>cond.lostparticles_maxenergy[js]:
                            e0maxnew=1.1*e0max
                            l_rescale_energy_array=1
                        else:
                            e0maxnew=cond.lostparticles_maxenergy[js]
                        if l_rescale_energy_array:
                            newlostparticles_energies=zeros(1001,'d')
                            tmpcount=zeros(1001,'d')
                            e0minold = cond.lostparticles_minenergy[js]
                            e0maxold = cond.lostparticles_maxenergy[js]
                            e0old = e0minold+arange(1001)*(e0maxold-e0minold)/1000
                            deposgrid1d(1,1001,e0old,cond.lostparticles_energies[js],1000,newlostparticles_energies,tmpcount,e0minnew,e0maxnew)
                            try:
                                cond.itrescale.append(top.it)
                            except:
                                cond.itrescale=[top.it]
                            cond.lostparticles_minenergy[js]=e0minnew
                            cond.lostparticles_maxenergy[js]=e0maxnew
                            cond.lostparticles_energies[js]=newlostparticles_energies
                        if top.wpid >0:
                            eweights = weight*top.pgroup.sw[js]
                        else:
                            eweights = ones(n)*top.pgroup.sw[js]
                        setgrid1dw(shape(coseta)[0],arccos(coseta),eweights,180,cond.lostparticles_angles[js],0.,pi)
                        # --- Expand the range of energies by one approximately two cells.
                        # --- This is needed in case there is only one particle, where
                        # --- minenergy == maxenergy.
                        setgrid1dw(shape(e0)[0],e0,eweights,1000,cond.lostparticles_energies[js],0.999*cond.lostparticles_minenergy[js],1.001*cond.lostparticles_maxenergy[js])
                        if js==1:
                            try:
                                cond.sumlostw.append(sum(cond.lostparticles_energies[1]))
                            except:
                                cond.sumlostw = AppendableArray(typecode='d')
                                cond.sumlostw.append(sum(cond.lostparticles_energies[1]))
                            try:
                                cond.minlostw.append(sum(cond.lostparticles_minenergy[1]))
                            except:
                                cond.minlostw = AppendableArray(typecode='d')
                                cond.minlostw.append(sum(cond.lostparticles_minenergy[1]))
                            try:
                                cond.maxlostw.append(sum(cond.lostparticles_maxenergy[1]))
                            except:
                                cond.maxlostw = AppendableArray(typecode='d')
                                cond.maxlostw.append(sum(cond.lostparticles_maxenergy[1]))
                        # else:
                        #   setgrid1d(shape(coseta)[0],arccos(coseta),180,cond.lostparticles_angles[js],0.,pi)
                        #   # --- Expand the range of energies by one approximately two cells.
                        #   # --- This is needed in case there is only one particle, where
                        #   # --- minenergy == maxenergy.
                        #   setgrid1d(shape(e0)[0],e0,1000,cond.lostparticles_energies[js],0.999*cond.lostparticles_minenergy[js],1.001*cond.lostparticles_maxenergy[js])
                    for ie,emitted_species in enumerate(self.inter[incident_species]['emitted_species'][ics]):
                        js_new=emitted_species.jslist[0]
                        forced_yield = self.inter[incident_species]['forced_yield'][ics]
                        tstart = wtime()


                        # print 'v',v[0][i],v[1][i],v[2][i],i,iit[i],js
                        # print 'x',[xplost[i],yplost[i],zplost[i]]
                        # print 'xold',[xplostold[i],yplostold[i],zplostold[i]]
                        # print 'u',[uxplost[i],uyplost[i],uzplost[i]]
                        # print 'e0',e0[i]
                        # coseta[i] = -sum(v[0][i]*n_unit0[0][i]+v[1][i]*n_unit0[1][i]+v[2][i]*n_unit0[2][i])/sqrt(sum(v[0][i]*v[0][i]+v[1][i]*v[1][i]+v[2][i]*v[2][i]))
                        # print 'coseta',coseta[i]
                        l_warning=0
                        l_infinity=0
                        if any(coseta < 0.):
                            l_warning=1
                            swarn = 'WARNING issued by Secondaries.generate: coseta<0.'
                            inegative = (coseta < 0.)
                            coseta[inegative]=-coseta[inegative]
                            n_unit0[0][inegative]=-n_unit0[0][inegative]
                            n_unit0[1][inegative]=-n_unit0[1][inegative]
                            n_unit0[2][inegative]=-n_unit0[2][inegative]
                            costheta[inegative] = cos(pi+theta[inegative])
                            sintheta[inegative] = sin(pi+theta[inegative])
                            # print 'coseta 1,2 :',coseta[i],-sum(u[i]*n_unit0[i])/sqrt(sum(u[i]*u[i]))
                            # print 'n 1, 2',n_unit0[0][i],n_unit0[1][i],n_unit0[2][i],sintheta[i]*cosphi[i],sintheta[i]*sinphi[i],costheta[i]
                        if any(xplost==largepos) or any(yplost==largepos) or any(zplost==largepos):
                            l_warning=1
                            l_infinity=1
                            swarn = 'WARNING issued by Secondaries.generate: particle at infinity'
                        if l_warning and self.l_verbose:
                            print swarn
                            # print 'phi, theta',phi[i],theta[i]
                            # print 'n',n_unit0[0][i],n_unit0[1][i],n_unit0[2][i]
                            # print 'u',u[0][i],u[1][i],u[2][i]
                        if l_infinity:
                            continue
                        if self.l_record_timing:
                            tinit+=wtime()-tstart
                            tstart=wtime()
                        if self.l_verbose:
                            for i in range(n):
                                print 'e0, coseta',e0[i],coseta[i]

                        init_position_offset = self.inter[incident_species]['init_position_offset'][ics]
                        if forced_yield is not None:
                            ####################################################################
                            # emission with imposed yield (note that initial velocity is tiny) #
                            ####################################################################
                            for i in range(n):
                                if isinstance(forced_yield,float):
                                    ns=int(forced_yield)
                                    if ranf()<(forced_yield-ns):ns+=1
                                else:
                                    ns = forced_yield
                                if init_position_offset>0.:
                                    xnew = xplost[i]+n_unit0[0][i]*init_position_offset
                                    ynew = yplost[i]+n_unit0[1][i]*init_position_offset
                                    znew = zplost[i]+n_unit0[2][i]*init_position_offset
                                else:
                                    xnew = xplost[i]
                                    ynew = yplost[i]
                                    znew = zplost[i]
                                tmp = ones(ns,'d')
                                if self.l_trackssnparents:
                                    ssnparent = tmp*ssnplost[i]
                                else:
                                    ssnparent = None
                                self.inter[incident_species]['emitted'][ics][ie] += ns*top.pgroup.sq[js_new]*top.pgroup.sw[js_new]
                                self.inter[incident_species]['absorbed'][ics][ie] += top.pgroup.sq[js]*top.pgroup.sw[js]
                                if top.wpid==0:
                                    self.addpart(ns,xnew,ynew,znew,1.e-10*tmp,1.e-10*tmp,1.e-10*tmp,js_new,itype=None,ssnparent=ssnparent)
                                else:
                                    self.addpart(ns,xnew,ynew,znew,1.e-10*tmp,1.e-10*tmp,1.e-10*tmp,js_new,ones(ns)*weight[i],None,ssnparent=ssnparent)

                        elif emitted_species.type is Electron:
                            if incident_species.type is Electron:
                                ##########################################
                                # electron-induced emission of electrons #
                                ##########################################
                                if self.piditype>0:
                                    itypes=AppendableArray(typecode='i',autobump=100)
                                else:
                                    itypes=None
                                
                                itype=self.inter[incident_species]['type'][ics]
                                scale_factor=self.inter[incident_species]['scale_factor'][ics]
                                if scale_factor is None:scale_factor=1.
                                if not self.flag_pyecloud:
                                    self.prepare_secondaries(itype,posC.maxsec)
                                if top.wpid==0:weight=ones(n,'d')

                                if self.flag_pyecloud:
                                    if top.wpid == 0:
                                        raise ValueError('Not compatible with pyecloud')

                                    # DEBUG ##########################
                                    from collections import OrderedDict
                                    debug_dict = OrderedDict([
                                        ('uxplost', uxplost),
                                        ('uyplost', uyplost),
                                        ('uzplost', uzplost),
                                        ('gaminvlost', gaminvlost),
                                        ('e0', e0),
                                        ('Electron.mass', Electron.mass)
                                        ])
                                    for kk in debug_dict.keys():
                                        print('%s: %s'%(kk, repr(debug_dict[kk])))
                                    ####################################

                                    (xnew, ynew, znew,
                                     uxsec, uysec, uzsec,
                                       weightsec) = self.pyecloud_secondary_emission(
                                        sintheta[:n], costheta[:n], sinphi[:n], cosphi[:n],
                                        weight[:n],
                                        xplost[:n], yplost[:n], zplost[:n],
                                        vxplost[:n], vyplost[:n], vzplost[:n])
                                    ns = len(xnew)

                                else:
                                    xnew = zeros(n*posC.maxsec,'d')
                                    ynew = zeros(n*posC.maxsec,'d')
                                    znew = zeros(n*posC.maxsec,'d')
                                    uxsec = zeros(n*posC.maxsec,'d')
                                    uysec = zeros(n*posC.maxsec,'d')
                                    uzsec = zeros(n*posC.maxsec,'d')
                                    ns = array([n*posC.maxsec])
                                    pos.secelecarray(n,e0,coseta,weight,ns,xnew,ynew,znew,uxsec,uysec,uzsec,
                                                costheta,sintheta,sinphi,cosphi,n_unit0,xplost,yplost,zplost,vxplost,vyplost,vzplost,
                                                top.pgroup.sm[js_new],Electron.mass,scale_factor,init_position_offset)
                                    ns = ns[0]
                                    if ns>0:
                                        xnew=xnew[:ns]
                                        ynew=ynew[:ns]
                                        znew=znew[:ns]
                                        uxsec=uxsec[:ns]
                                        uysec=uysec[:ns]
                                        uzsec=uzsec[:ns]
                                        weightsec = weight[0] + 0 * xnew
                            else: # incidents are atoms or ions
                                ##########################################
                                # atom/ion-induced emission of electrons #
                                ##########################################
                                xnew=AppendableArray(typecode='d',autobump=100)
                                ynew=AppendableArray(typecode='d',autobump=100)
                                znew=AppendableArray(typecode='d',autobump=100)
                                uxsec=AppendableArray(typecode='d',autobump=100)
                                uysec=AppendableArray(typecode='d',autobump=100)
                                uzsec=AppendableArray(typecode='d',autobump=100)
                                if self.piditype>0:
                                    itypes=AppendableArray(typecode='i',autobump=100)
                                else:
                                    itypes=None
                                ns = 0
                                if incident_species.type.__class__ in [Atom,Molecule]:
                                    if self.inter[incident_species]['material'][ics]=='SS':target_num=10025
                                    scale_factor = self.inter[incident_species]['scale_factor'][ics]
                                    if scale_factor is None or scale_factor==1.:
                                        nbatches = 1
                                        l_scale_factor = false
                                    else:
                                        nbatches = int(scale_factor)+1
                                        l_scale_factor = true
                                    for i in range(n):
                                        untx=AppendableArray(typecode='d',autobump=10)
                                        uttx=AppendableArray(typecode='d',autobump=10)
                                        uztx=AppendableArray(typecode='d',autobump=10)
                                        for ibatch in range(nbatches):
                                            # --- version 0.2.1
                                            # ns=txphysics.ion_ind_elecs(e0[i]/(1.e6*incident_species.type.A),
                                            #                            max(0.04,coseta[i]),
                                            #                            float(incident_species.type.Z),
                                            #                            incident_species.type.A,
                                            #                            target_num,
                                            #                            self.emitted_e,
                                            #                            self.emitted_bn,
                                            #                            self.emitted_bt,
                                            #                            self.emitted_bz)
                                            # --- version 1.9.0
                                            ion_ind_e0 = zeros(1,'d')
                                            ion_ind_ct = zeros(1,'d')
                                            ion_ind_e0[0] = e0[i]/(1.e6)
                                            ion_ind_ct[0] = max(0.04,coseta[i])
                                            Te=0.
                                            Ne=0.
                                            self.emitted_e,self.emitted_bn,self.emitted_bt,self.emitted_bz = \
                                            txigenelec.ion_ind_elecs(ion_ind_e0,
                                                                    ion_ind_ct,
                                                                    float(incident_species.type.Z),
                                                                    incident_species.type.A,
                                                                    Te,Ne,
                                                                    target_num)
                                            nsemit = self.emitted_e.size
                                            if l_scale_factor:
                                                emitfrac=scale_factor-ibatch
                                            if l_scale_factor and emitfrac<1.:
                                                for iemit in range(nsemit):
                                                    if self.emitted_e[iemit]>=0. and ranf()<emitfrac:
                                                        gammac = clight/sqrt(1.-self.emitted_bn[iemit]**2 \
                                                                               +self.emitted_bt[iemit]**2 \
                                                                               +self.emitted_bz[iemit]**2)
                                                        untx.append(gammac*self.emitted_bn[iemit])
                                                        uttx.append(gammac*self.emitted_bt[iemit])
                                                        uztx.append(gammac*self.emitted_bz[iemit])
                                                        ns+=1
                                            else:
                                                for iemit in range(nsemit):
                                                    if self.emitted_e[iemit]>=0.:
                                                        gammac = clight/sqrt(1.-self.emitted_bn[iemit]**2 \
                                                                               +self.emitted_bt[iemit]**2 \
                                                                               +self.emitted_bz[iemit]**2)
                                                        untx.append(gammac*self.emitted_bn[iemit])
                                                        uttx.append(gammac*self.emitted_bt[iemit])
                                                        uztx.append(gammac*self.emitted_bz[iemit])
                                                        ns+=1

                                        un=untx.data()
                                        ut=uttx.data()
                                        uz=uztx.data()
                                        if self.l_verbose:
                                            print 'nb secondaries = ',ns,' from conductor ',icond, e0[i], coseta[i],i1,i2,iit[i],top.npslost
                                        x,y,z,ux,uy,uz =  self.getxv(i,costheta,sintheta,sinphi,cosphi,n_unit0,xplost,yplost,zplost,vxplost,vyplost,vzplost,un,ut,uz,init_position_offset)
                                        xnew.append(x)
                                        ynew.append(y)
                                        znew.append(z)
                                        uxsec.append(ux)
                                        uysec.append(uy)
                                        uzsec.append(uz)
                                xnew=xnew.data()
                                ynew=ynew.data()
                                znew=znew.data()
                                uxsec=uxsec.data()
                                uysec=uysec.data()
                                uzsec=uzsec.data()

                            if self.l_record_timing:
                                tgen+=wtime()-tstart
                                tstart=wtime()
                            if ns>0:
                                # Record number of emitted
                                self.inter[incident_species]['emitted'][ics][ie] += ns*top.pgroup.sq[js_new]*top.pgroup.sw[js_new]
                                self.inter[incident_species]['absorbed'][ics][ie] += top.pgroup.sq[js]*top.pgroup.sw[js]
                                
                                # pid[:,self.xoldpid]=xnew-vx*top.dt
                                # pid[:,self.yoldpid]=ynew-vy*top.dt
                                # pid[:,self.zoldpid]=znew-vz*top.dt
                                # --- apply perdiodic BC
                                # znew = where(znew<zmin,zmax-zmin+znew,znew)
                                # znew = where(znew>zmax,zmin-zmax+znew,znew)
                                if 0:
                                    if w3d.solvergeom==w3d.RZgeom:
                                        condition = (sqrt(xnew**2+ynew**2)>xmax) or \
                                                    (znew<zmin) or (znew>zmax)
                                    else:
                                        condition = (xnew<xmin) or (xnew>xmax) or \
                                                    (ynew<ymin) or (ynew>ymax) or \
                                                    (znew<zmin) or (znew>zmax)
                                    if condition:
                                        print 'WARNING from secondaries: new particle outside boundaries, skip creation',
                                        print '\nLost particle position: ',xplost[i],yplost[i],zplost[i],
                                        print '\nNew particle position: ',xnew,ynew,znew
                                        print 'XYZ min/max: ', xmin,xmax,ymin,ymax,zmin,zmax
                                        # self.outparts+=[[xnew,ynew,znew,xplost[i],yplost[i],zplost[i], \
                                        # xplostold[i],yplostold[i],zplostold[i],n_unit0[0][i],n_unit0[1][i],n_unit0[2][i],icond]]
                                        self.outparts+=[[xnew,ynew,znew,xplost[i],yplost[i],zplost[i], \
                                        n_unit0[0][i],n_unit0[1][i],n_unit0[2][i],icond]]
                                #else:
                                if self.l_trackssnparents:
                                    ssnparent = ones(ns)*ssnplost[i]
                                else:
                                    ssnparent = None
                                if top.wpid==0:
                                    e0emit = 0.5*top.pgroup.sm[js_new]*top.pgroup.sw[js_new]*(uxsec*uxsec+uysec*uysec+uzsec*uzsec)
                                    pxsum -= sum(top.pgroup.sm[js_new]*top.pgroup.sw[js_new]*uxsec)
                                    pysum -= sum(top.pgroup.sm[js_new]*top.pgroup.sw[js_new]*uysec)
                                    pzsum -= sum(top.pgroup.sm[js_new]*top.pgroup.sw[js_new]*uzsec)
                                    self.addparticles(ns,xnew,ynew,znew,uxsec,uysec,uzsec,js_new,itype=itypes,ssnparent=ssnparent)
                                else:
                                    e0emit = 0.5*top.pgroup.sm[js_new]*top.pgroup.sw[js_new]*weight[i]*(uxsec*uxsec+uysec*uysec+uzsec*uzsec)
                                    pxsum -= sum(top.pgroup.sm[js_new]*top.pgroup.sw[js_new]*weight[i]*uxsec)
                                    pysum -= sum(top.pgroup.sm[js_new]*top.pgroup.sw[js_new]*weight[i]*uysec)
                                    pzsum -= sum(top.pgroup.sm[js_new]*top.pgroup.sw[js_new]*weight[i]*uzsec)
                                    self.addparticles(ns,xnew,ynew,znew,uxsec,uysec,uzsec,js_new,weightsec,ssnparent=ssnparent)
                                ek0emitsum += sum(e0emit)

                            if self.l_record_timing:
                                tadd+=wtime()-tstart
                                tstart=wtime()

                        elif hasattr(emitted_species,'charge_state') and emitted_species.charge_state==0:
                            ########################
                            # emission of neutrals #
                            ########################
                            specular_fraction = self.inter[incident_species]['specular_fraction'][ics]
                            for i in range(n):
                                if incident_species.sq == 0. or specular_fraction is not None:
                                    ##################################################
                                    # emission of neutral from neutral or ion impact #
                                    ##################################################
                                    # --- This uses a model where a fraction is specularly reflected and the rest
                                    # --- are reemitted at the temperature of the conductor. An ion that impacts is
                                    # --- always neutralized before being emitted. One neutral is always emitted
                                    # --- for each incident ion or neutral.
                                    ns = 1
                                    if random.random() < specular_fraction:
                                        # --- Specular reflection
                                        vgt = cosphi[i]*costheta[i]*vxplost[i] + sinphi[i]*costheta[i]*vyplost[i] - sintheta[i]*vzplost[i]
                                        vgz =            -sinphi[i]*vxplost[i] +             cosphi[i]*vyplost[i]
                                        vgn = cosphi[i]*sintheta[i]*vxplost[i] + sinphi[i]*sintheta[i]*vyplost[i] + costheta[i]*vzplost[i]
                                        vt0 = vgt
                                        vz0 = vgz
                                        vn0 = -vgn
                                    else:
                                        # --- Thermal emission
                                        conductor_temperature = self.inter[incident_species]['conductor_temperature'][ics]
                                        vth = sqrt(2.*conductor_temperature*jperev/emitted_species.mass)
                                        vmag = vth*sqrt(2.*log(1./max(1.e-14, random.random())))
                                        thnew = pi/2.*random.random()
                                        phnew = 2.*pi*random.random()
                                        vgt = vmag*sin(thnew)*sin(phnew)
                                        vgn = vmag*cos(thnew)
                                        vgz = vmag*sin(thnew)*cos(phnew)
                                        if costheta[i] < 1.-1.e-10:
                                            z_unit0 = array([-sinphi[i],cosphi[i],0.])
                                            z = -array([vzplost[i]*n_unit0[1][i]-vyplost[i]*n_unit0[2][i],
                                                        vxplost[i]*n_unit0[2][i]-vzplost[i]*n_unit0[0][i],
                                                        vyplost[i]*n_unit0[0][i]-vxplost[i]*n_unit0[1][i]])
                                            z_unit  = z/sqrt(sum(z*z))
                                            cospsi  = sum(z_unit*z_unit0)
                                            sinpsi  = sqrt(max(0.,1.-cospsi*cospsi))
                                            vt0 = -(cospsi*vgt - sinpsi*vgz)
                                            vz0 = -(sinpsi*vgt + cospsi*vgz)
                                            vn0 = vgn
                                        else:
                                            vt0 = vgt
                                            vz0 = vgz
                                            vn0 = vgn

                                    vxnew = cosphi[i]*costheta[i]*vt0 - sinphi[i]*vz0 + cosphi[i]*sintheta[i]*vn0
                                    vynew = sinphi[i]*costheta[i]*vt0 + cosphi[i]*vz0 + sinphi[i]*sintheta[i]*vn0
                                    vznew =          -sintheta[i]*vt0                 +           costheta[i]*vn0

                                elif incident_species.sq != 0.:
                                    ############################################
                                    # emission of neutrals from neutral impact #
                                    ############################################
                                    # --- This is used for high energy incident ion.
                                    my_yield=1.+1.82e-4*exp(0.09*180./pi*arccos(coseta[i]))
                                    scale_factor = self.inter[incident_species]['scale_factor'][ics]
                                    if not (scale_factor is None or scale_factor==1.):
                                        my_yield*=scale_factor
                                    ns = int(my_yield)
                                    # --- The ns+1 is only a temporary fix to avoid an array out of
                                    # --- bounds errors. Once the desorb routine is fixed, this
                                    # --- code should be updated. Note that as the code is now,
                                    # --- in some cases, a desorbed particle will be thrown out.
                                    vx = desorb.floatArray(ns+1)
                                    vy = desorb.floatArray(ns+1)
                                    vz = desorb.floatArray(ns+1)
                                    vxnew = zeros(ns,'d')
                                    vynew = zeros(ns,'d')
                                    vznew = zeros(ns,'d')

                                    # --- compute the desorbed neutrals
                                    # --- note that the gamma0 (1.) and rel_weight (top.pgroup.sw[js_new]/top.pgroup.sw[js]) are actually not used
                                    desorb.desorb(my_yield,v[0][i],v[1][i],v[2][i],theta[i],phi[i],1.,0.4,top.pgroup.sm[js],top.pgroup.sw[js_new]/top.pgroup.sw[js],vx,vy,vz)
                                    scale_factor_velocity = self.inter[incident_species]['scale_factor_velocity'][ics]
                                    for ivnew in range(ns):
                                        vxnew[ivnew]=vx[ivnew]*scale_factor_velocity
                                        vynew[ivnew]=vy[ivnew]*scale_factor_velocity
                                        vznew[ivnew]=vz[ivnew]*scale_factor_velocity


                                self.inter[incident_species]['emitted'][ics][ie] += ns*top.pgroup.sq[js_new]*top.pgroup.sw[js_new]
                                self.inter[incident_species]['absorbed'][ics][ie] += top.pgroup.sq[js]*top.pgroup.sw[js]

                                init_position_offset = self.inter[incident_species]['init_position_offset'][ics]
                                if init_position_offset>0.:
                                    xnew = xplost[i]+n_unit0[0][i]*init_position_offset
                                    ynew = yplost[i]+n_unit0[1][i]*init_position_offset
                                    znew = zplost[i]+n_unit0[2][i]*init_position_offset
                                else:
                                    xnew = xplost[i]
                                    ynew = yplost[i]
                                    znew = zplost[i]
                                # pid[:,self.xoldpid]=xnew-vxnew*top.dt
                                # pid[:,self.yoldpid]=ynew-vynew*top.dt
                                # pid[:,self.zoldpid]=znew-vznew*top.dt
                                if w3d.solvergeom==w3d.RZgeom:
                                    condition = (sqrt(xnew**2+ynew**2)>xmax) or \
                                                (znew<zmin) or (znew>zmax)
                                else:
                                    condition = (xnew<xmin) or (xnew>xmax) or \
                                                (ynew<ymin) or (ynew>ymax) or \
                                                (znew<zmin) or (znew>zmax)
                                if condition and self.luseoldpositionifout:
                                    xnew = xplostold[i]
                                    ynew = yplostold[i]
                                    znew = zplostold[i]
                                    condition = False
                                if condition:
                                    print 'WARNING from secondaries: new neutral particle outside boundaries, skip creation',
                                    print '\nLost particle position: ',xplost[i],yplost[i],zplost[i],
                                    print '\nNew particle position: ',xnew,ynew,znew
                                    if self.vmode==1:
                                        self.outparts+=[[xnew,ynew,znew,xplost[i],yplost[i],zplost[i], \
                                      vxplost[i],vyplost[i],vzplost[i],n_unit0[0][i],n_unit0[1][i],n_unit0[2][i],icond]]
                                    if self.vmode==2:
                                        self.outparts+=[[xnew,ynew,znew,xplost[i],yplost[i],zplost[i], \
                                      xplostold[i],yplostold[i],zplostold[i],n_unit0[0][i],n_unit0[1][i],n_unit0[2][i],icond]]
                                else:
                                    if self.l_trackssnparents:
                                        ssnparent = ones(ns)*ssnplost[i]
                                    else:
                                        ssnparent = None
                                    if top.wpid==0:
                                        e0emit = 0.5*top.pgroup.sm[js_new]*top.pgroup.sw[js_new]*(vxnew*vxnew+vynew*vynew+vznew*vznew)
                                        pxsum -= sum(top.pgroup.sm[js_new]*top.pgroup.sw[js_new]*vxnew)
                                        pysum -= sum(top.pgroup.sm[js_new]*top.pgroup.sw[js_new]*vynew)
                                        pzsum -= sum(top.pgroup.sm[js_new]*top.pgroup.sw[js_new]*vznew)
                                        self.addpart(ns,xnew,ynew,znew,vxnew,vynew,vznew,js_new,itype=None,ssnparent=ssnparent)
                                    else:
                                        e0emit = 0.5*top.pgroup.sm[js_new]*top.pgroup.sw[js_new]*weight*(vxnew*vxnew+vynew*vynew+vznew*vznew)
                                        pxsum -= sum(top.pgroup.sm[js_new]*top.pgroup.sw[js_new]*weight[i]*vxnew)
                                        pysum -= sum(top.pgroup.sm[js_new]*top.pgroup.sw[js_new]*weight[i]*vynew)
                                        pzsum -= sum(top.pgroup.sm[js_new]*top.pgroup.sw[js_new]*weight[i]*vznew)
                                        self.addpart(ns,xnew,ynew,znew,vxnew,vynew,vznew,js_new,ones(ns)*weight[i],None,ssnparent=ssnparent)
                                    ek0emitsum += sum(e0emit)

                            # if self.l_record_timing:
                            #   tadd+=wtime()-tstart

        if self.l_record_timing:t3 = time.clock()
        # --- make sure that all particles are added
        for js in self.x:
            self.flushpart(js)
        # --- Check for particle out of bounds and exchange particles among
        # --- processors if needed. A call to particleboundaries3d is made
        # --- so that particles are scraped on any user defined conductors
        # --- as well as the grid boundaries.
        # --- Set the flag so that this generate routine is not called again
        # --- recursively (since generate is normally called at the end of
        # --- the scraping).
        # --- This is all needed in case lcallscrapercontrollers is true.
        # --- Also set the flag so that the lost particles are not reset.
        if not local:
            self.lrecursivegenerate = 1
            top.lresetlostpart = false
            particleboundaries3d(top.pgroup,-1,self.lcallscrapercontrollers)
            top.lresetlostpart = true
            self.lrecursivegenerate = 0

        # print "tinit,tgen,tadd:",tinit*1.e-6,tgen*1.e-6,tprepadd*1.e-6,tadd*1.e-6
        # --- append total emitted charge in conductors emitparticles_data arrays
        for js in self.inter:
            for ics,c in enumerate(self.inter[js]['conductors']):
                for ie in range(len(self.inter[js]['emitted'][ics])):
                    if local:
                        totemit = sum(self.inter[js]['emitted'][ics][ie])
                        totabsorb = sum(self.inter[js]['absorbed'][ics][ie])
                    else:
                        totemit = globalsum(self.inter[js]['emitted'][ics][ie])
                        totabsorb = globalsum(self.inter[js]['absorbed'][ics][ie])
                    if abs(totemit)>0.:
                        c.emitparticles_data.append(array([top.time,
                                                           totemit,
                                                           top.dt,
                                                           self.inter[js]['emitted_species'][ics][ie].jslist[0],
                                                           totabsorb]))

        # --- append history arrays
        if l_accumulate_hist:
            weighttot = globalsum(weighttot)
            ek0sum = globalsum(ek0sum)
            costhsum = globalsum(costhsum)
            ek0max = globalmax(ek0max)
            pxsum = globalsum(pxsum)
            pysum = globalsum(pysum)
            pzsum = globalsum(pzsum)
            if me==0:
                if weighttot != 0.:
                    self.htime.append(top.time)
                    self.ek0av.append(ek0sum/weighttot)      #cummulative collision kinetic energy [eV] this step
                    self.ek0max.append(ek0max)      #maximum collision kinetic energy [eV]
                    self.costhav.append(costhsum/weighttot)
                    self.pxsum.append(pxsum)
                    self.pysum.append(pysum)
                    self.pzsum.append(pzsum)
                    self.power_dep.append(ek0sum*jperev/top.dt)
                    self.power_emit.append(ek0emitsum/top.dt)
                    self.power_diff.append((ek0sum*jperev-ek0emitsum)/top.dt)
        # w3d.lcallscraper=0
        # particleboundaries3d(top.pgroup,-1,false)
        # w3d.lcallscraper=1
        # top.npslost=0
        if self.l_record_timing:t4 = time.clock()
        if self.l_record_timing:self.timings.append([t4-t1,t2-t1,t3-t2,t4-t3,tinit,tgen,tprepadd,tadd])
        # print 'time Secondaries = ',time.clock()-t1,'s',t2-t1,t3-t2,t4-t3
        if self.l_verbose>1:print 'secondaries generation finished'

    def call_set_params_user(self,maxsec,mat_num=None):
        
        if self.flag_pyecloud:
            raise ValueError('This method shound not be used with pyecloud')

        # --- Always call the default version of the routine. The user's routine
        # --- only needs to set the parameters that are different.
        if not self.l_set_params_user_only:self.set_params(maxsec,mat_num)
        # --- Now call the user's routine if there is one.
        if self.set_params_user is not None:
            if isinstance(self.set_params_user,basestring):
                    # --- This is needed primarily since a user defined function cannot be
                    # --- directly picklable and so only the name is saved.
                import __main__
                try:
                    self.set_params_user = __main__.__dict__[self.set_params_user]
                except KeyError:
                    # --- Maybe this should raise an error?
                    print "Warning: Secondaries set_params_user function '%s' is not defined."%self.set_params_user
            if callable(self.set_params_user):
                self.set_params_user(maxsec,mat_num)
            else:
                # --- Maybe this should raise an error?
                print "Warning: Secondaries set_params_user function '%s' is not callable."%self.set_params_user

    def getxv(self,i,costheta,sintheta,sinphi,cosphi,n_unit0,xplost,yplost,zplost,vxplost,vyplost,vzplost,un,ut,uz,init_position_offset):
        if costheta[i]<1.-1.e-10:
            z_unit0 = array([-sinphi[i],cosphi[i],0.])
            # z       = -array([uzplost[i]*n_unit0[1][i]-uyplost[i]*n_unit0[2][i],
            #                   uxplost[i]*n_unit0[2][i]-uzplost[i]*n_unit0[0][i],
            #                   uyplost[i]*n_unit0[0][i]-uxplost[i]*n_unit0[1][i]])
            z       = -array([vzplost[i]*n_unit0[1][i]-vyplost[i]*n_unit0[2][i],
                              vxplost[i]*n_unit0[2][i]-vzplost[i]*n_unit0[0][i],
                              vyplost[i]*n_unit0[0][i]-vxplost[i]*n_unit0[1][i]])
            z_unit  = z/sqrt(sum(z*z))
            cospsi  = sum(z_unit*z_unit0)
            sinpsi  = sqrt(max(0.,1.-cospsi*cospsi))
            # bt0 = (cospsi*bt - sinpsi*bz)
            # bz0 = (sinpsi*bt + cospsi*bz)
            ut0 = -(cospsi*ut - sinpsi*uz)
            uz0 = -(sinpsi*ut + cospsi*uz)
            un0 = un
        else:
            ut0 = ut
            uz0 = uz
            un0 = un
        uxsec = cosphi[i]*costheta[i]*ut0 - sinphi[i]*uz0 + cosphi[i]*sintheta[i]*un0
        uysec = sinphi[i]*costheta[i]*ut0 + cosphi[i]*uz0 + sinphi[i]*sintheta[i]*un0
        uzsec =          -sintheta[i]*ut0                 +           costheta[i]*un0
        tmp=ones(shape(un)[0])
        del un,ut,uz,un0,ut0,uz0
        xnew = xplost[i]*tmp
        ynew = yplost[i]*tmp
        znew = zplost[i]*tmp
        if init_position_offset>0.:
            xnew +=n_unit0[0][i]*init_position_offset
            ynew +=n_unit0[1][i]*init_position_offset
            znew +=n_unit0[2][i]*init_position_offset
        return xnew,ynew,znew,uxsec,uysec,uzsec

    def pyecloud_secondary_emission(self, sintheta, costheta, sinphi, cosphi,
                                        weightplost, 
                                        xplost, yplost, zplost,
                                        uxplost, uyplost, uzplost):

        #print('type(sintheta)', type(sintheta))
        #print('repr(sintheta)', repr(sintheta))
        # Build normal and tangent unit vectors 
        it1_impact = np.array([
            cosphi * costheta,
            sinphi * costheta,
            -sintheta])
        
        it2_impact = np.array([
            -sinphi,
            cosphi,
            0. * sintheta])
        
        in_impact = np.array([
            cosphi*sintheta,
            sinphi*sintheta, 
            costheta])

        # Define local reference system (_impact)
        ix_impact = it1_impact
        iy_impact = in_impact
        iz_impact = it2_impact

        # Transform velocities
        u_vect = np.array([
            uxplost,
            uyplost,
            uzplost,
            ])
        vx_impact = coldotprod(u_vect, ix_impact)
        vy_impact = coldotprod(u_vect, iy_impact)
        vz_impact = coldotprod(u_vect, iz_impact)

        # There is no need to transform the positions
        # (they are just copied from primary to secondary)
        x_impact = xplost
        y_impact = yplost
        z_impact = zplost

        # We defined our reference system so the the iy is the normal
        Norm_x = np.zeros_like(xplost)
        Norm_y = np.ones_like(xplost)

        # Some auxiliary quantity
        v_impact_n = vy_impact
        v_impact_mod = np.sqrt(np.sum(u_vect*u_vect, axis=0))
        costheta_impact = np.abs(v_impact_n / v_impact_mod)
        E_impact_eV = 0.5 * e_mass / qe * v_impact_mod * v_impact_mod
        i_found = np.arange(0, len(x_impact)) # I use this to trace the emitting MP
        flag_seg = True
        nel_mp_th = self.pyecloud_fact_split * self.pyecloud_nel_mp_ref
        nel_impact = weightplost
        
        (nel_emit_tot_events, event_type, event_info,
           nel_replace, x_replace, y_replace, z_replace,
           vx_replace, vy_replace, vz_replace, i_seg_replace,
           nel_new_MPs, x_new_MPs, y_new_MPs, z_new_MPs,
           vx_new_MPs, vy_new_MPs, vz_new_MPs, i_seg_new_MPs,
           ) = self.pyeclsecemi.impacts_on_surface(
           e_mass, nel_impact, x_impact, y_impact, z_impact,
           vx_impact, vy_impact, vz_impact, Norm_x, Norm_y, i_found,
           v_impact_n, E_impact_eV, costheta_impact, nel_mp_th, flag_seg)

        x_impact_new = np.concatenate([x_replace, x_new_MPs])
        y_impact_new = np.concatenate([y_replace, y_new_MPs])
        z_impact_new = np.concatenate([z_replace, z_new_MPs])
        vx_impact_new = np.concatenate([vx_replace, vx_new_MPs])
        vy_impact_new = np.concatenate([vy_replace, vy_new_MPs])
        vz_impact_new = np.concatenate([vz_replace, vz_new_MPs])
        weight_impact_new = np.concatenate([nel_replace, nel_new_MPs])

        if len(i_seg_new_MPs)>0:
            ix_impact_new = np.concatenate((ix_impact,
                    np.take(ix_impact, i_seg_new_MPs, axis=1)), axis=1)
            iy_impact_new = np.concatenate((iy_impact,
                    np.take(iy_impact, i_seg_new_MPs, axis=1)), axis=1)
            iz_impact_new = np.concatenate((iz_impact,
                    np.take(iz_impact, i_seg_new_MPs, axis=1)), axis=1)
        else:
            ix_impact_new = ix_impact
            iy_impact_new = iy_impact
            iz_impact_new = iz_impact



        # Remove small macroparticles
        mask_keep = weight_impact_new > (self.pyecloud_fact_clean * self.pyecloud_nel_mp_ref)
        
        xnew = x_impact_new[mask_keep]
        ynew = y_impact_new[mask_keep]
        znew = z_impact_new[mask_keep]
        vx_impact_new = vx_impact_new[mask_keep]
        vy_impact_new = vy_impact_new[mask_keep]
        vz_impact_new = vz_impact_new[mask_keep]
        weightnew = weight_impact_new[mask_keep]
        
        ix_impact_new = ix_impact_new[:, mask_keep] 
        iy_impact_new = iy_impact_new[:, mask_keep] 
        iz_impact_new = iz_impact_new[:, mask_keep] 


        # Transform velocities to warp frame
        u_new = vx_impact_new * ix_impact_new + vy_impact_new * iy_impact_new + vz_impact_new * iz_impact_new
        uxnew = u_new[0, :]
        uynew = u_new[1, :]
        uznew = u_new[2, :]

        # DEBUG ############################
        from collections import OrderedDict
        debug_dict = OrderedDict([
            ('sintheta', sintheta),
            ('costheta', costheta),
            ('sinphi', sinphi),
            ('cosphi', cosphi),
            ('weightplost', weightplost),
            ('xplost', xplost),
            ('yplost', yplost),
            ('zplost', zplost),
            ('uxplost', uxplost),
            ('uyplost', uyplost),
            ('uzplost', uzplost),
            ])
        for kk in debug_dict.keys():
            print('%s: %s'%(kk, repr(debug_dict[kk])))
        ####################################


        return xnew, ynew, znew, uxnew, uynew, uznew, weightnew


    def set_params(self,maxsec,mat_num=None):
        """This routine sets the parameters for a given material.
        Written by Peter Stoltz.  All values taken from Furman's paper

        Here is a list of ions & material numbers:
        e- -> Cu - mat_num=1 (default)
        e- -> Stainless Steel - mat_num=2
        H  -> Au - mat_num=3
        He -> Au - mat_num=4
        K  -> SS - mat_num=5
        """
        if self.flag_pyecloud:
            raise ValueError('This method shound not be used with pyecloud')
        
        
        if pos is None:
            return
        if mat_num is None:
            raise Exception("Error in set_params: mat_num has not been setup.")

        # --- Set the values
        # --- maxsec must be ten or greater
        # posC.th  = 0.
        # posC.pn = 0.
        # posC.pt = 0.
        # posC.pz = 0.
        # posC.ns = 0
        posC.matsurf =3
        # posC.ielswitch =1
        # posC.iredswitch =1
        # posC.itrueswitch =1
        # posC.irelk =0
        # posC.pmax =0.77939
        # posC.range = 0.
        # posC.freepath = 0.
        posC.iprob = 4
        # posC.ndelerm =0
        # posC.ndeltspm =0
        # posC.np0lt0 =0
        # posC.np1gt1 =0
        # posC.totsec1 =0.
        # posC.totsec2 =0.

        # --- Here we set material-dependent parameters...Cu (mat_num=1) is default
        posC.enpar[0] = 1.5
        posC.enpar[1] = 1.75
        posC.enpar[2] = 1.
        posC.enpar[3] = 3.75
        posC.enpar[4] = 8.5
        posC.enpar[5] = 11.5
        posC.enpar[6] = 2.5
        posC.enpar[7] = 3.0
        posC.enpar[8] = 2.5
        posC.enpar[9] = 3.0

        posC.pnpar[0] = 2.5
        posC.pnpar[1] = 3.3
        posC.pnpar[2] = 2.5
        posC.pnpar[3] = 2.5
        posC.pnpar[4] = 2.8
        posC.pnpar[5] = 1.3
        posC.pnpar[6] = 1.5
        posC.pnpar[7] = 1.5
        posC.pnpar[8] = 1.5
        posC.pnpar[9] = 1.5

        posC.dtspk = 1.8848
        posC.dtotpk = 2.1
        posC.pangsec =1.
        posC.pr =0.5
        posC.sige =2.
        posC.Ecr =0.0409
        posC.E0tspk =276.812
        posC.E0epk = 0.
        posC.E0w =60.8614
        posC.rpar1 =0.26
        posC.rpar2 =2.
        posC.tpar1 =0.66
        posC.tpar2 =0.8
        posC.tpar3 =0.7
        posC.tpar4 =1.
        posC.tpar5 =0.
        posC.tpar6 =0.
        posC.epar1 =0.26
        posC.epar2 =2.
        posC.P1rinf =0.2
        posC.P1einf =0.02
        posC.P1epk =0.49623
        posC.powts =1.54033
        posC.powe =1.
        posC.qr=0.104045


        # --- Now check for Stainless Steel
        if (mat_num == 2):

            posC.enpar[0] = 3.9
            posC.enpar[1] = 6.2
            posC.enpar[2] = 13.
            posC.enpar[3] = 8.8
            posC.enpar[4] = 6.25
            posC.enpar[5] = 2.25
            posC.enpar[6] = 9.2
            posC.enpar[7] = 5.3
            posC.enpar[8] = 17.8
            posC.enpar[9] = 10.

            posC.pnpar[0] = 1.6
            posC.pnpar[1] = 2.
            posC.pnpar[2] = 1.8
            posC.pnpar[3] = 4.7
            posC.pnpar[4] = 1.8
            posC.pnpar[5] = 2.4
            posC.pnpar[6] = 1.8
            posC.pnpar[7] = 1.8
            posC.pnpar[8] = 2.3
            posC.pnpar[9] = 1.8

            posC.pangsec =1.
            posC.pr =0.4
            posC.sige =1.9
            posC.dtspk = 1.22
            posC.dtotpk = 2.1
            posC.Ecr =40.0
            posC.E0tspk =310.
            posC.E0epk = 0.
            posC.E0w =100.
            posC.rpar1 =0.26
            posC.rpar2 =2.
            posC.tpar1 =0.66
            posC.tpar2 =0.8
            posC.tpar3 =0.7
            posC.tpar4 =1.
            posC.tpar5 =0.
            posC.tpar6 =0.
            posC.epar1 =0.26
            posC.epar2 =2.
            posC.P1rinf =0.74
            posC.P1einf =0.07
            posC.P1epk =0.5
            posC.powts =1.813
            posC.powe =0.9
            posC.qr=1.

        # --- Now check for Ion->Au
        # --- Because it is Ion, we set elastic and rediffused to zero
        #
        # --- First we take H ion (data from Eder, Rev Sci Inst 68 1(1997))
        if (mat_num == 3):
            posC.enpar[0] = 3.9
            posC.enpar[1] = 6.2
            posC.enpar[2] = 13.
            posC.enpar[3] = 8.8
            posC.enpar[4] = 6.25
            posC.enpar[5] = 2.25
            posC.enpar[6] = 9.2
            posC.enpar[7] = 5.3
            posC.enpar[8] = 17.8
            posC.enpar[9] = 10.

            posC.pnpar[0] = 1.6
            posC.pnpar[1] = 2.
            posC.pnpar[2] = 1.8
            posC.pnpar[3] = 4.7
            posC.pnpar[4] = 1.8
            posC.pnpar[5] = 2.4
            posC.pnpar[6] = 1.8
            posC.pnpar[7] = 1.8
            posC.pnpar[8] = 2.3
            posC.pnpar[9] = 1.8

            posC.pangsec =1.

            posC.dtspk = 2.382
            posC.dtotpk = 2.382
            posC.E0tspk = 104624.0
            posC.powts =1.44
            posC.tpar1 =0.66
            posC.tpar2 =0.8
            posC.tpar3 =0.7
            posC.tpar4 =1.
            posC.tpar5 =0.
            posC.tpar6 =0.

            posC.sige =0.0
            posC.pr =0.0
            posC.E0epk = 0.
            posC.E0w =0.
            posC.epar1 =0.0
            posC.epar2 =0.
            posC.P1einf =0.0
            posC.P1epk =0.0
            posC.powe =0.0

            posC.Ecr =0.0
            posC.rpar1 =0.0
            posC.rpar2 =0.
            posC.P1rinf =0.0
            posC.qr=0.

        # --- Now we take He ion (data from Eder, Rev Sci Inst 68 1(1997))
        if (mat_num == 4):
            posC.enpar[0] = 3.9
            posC.enpar[1] = 6.2
            posC.enpar[2] = 13.
            posC.enpar[3] = 8.8
            posC.enpar[4] = 6.25
            posC.enpar[5] = 2.25
            posC.enpar[6] = 9.2
            posC.enpar[7] = 5.3
            posC.enpar[8] = 17.8
            posC.enpar[9] = 10.

            posC.pnpar[0] = 1.6
            posC.pnpar[1] = 2.
            posC.pnpar[2] = 1.8
            posC.pnpar[3] = 4.7
            posC.pnpar[4] = 1.8
            posC.pnpar[5] = 2.4
            posC.pnpar[6] = 1.8
            posC.pnpar[7] = 1.8
            posC.pnpar[8] = 2.3
            posC.pnpar[9] = 1.8

            posC.pangsec =1.

            posC.dtspk = 5.68
            posC.dtotpk = 5.68
            posC.E0tspk = 1410000.0
            posC.powts =1.044
            posC.tpar1 =0.66
            posC.tpar2 =0.8
            posC.tpar3 =0.7
            posC.tpar4 =1.
            posC.tpar5 =0.
            posC.tpar6 =0.

            posC.sige =0.0
            posC.pr =0.0
            posC.E0epk = 0.
            posC.E0w =0.
            posC.epar1 =0.0
            posC.epar2 =0.
            posC.P1einf =0.0
            posC.P1epk =0.0
            posC.powe =0.0

            posC.Ecr =0.0
            posC.rpar1 =0.0
            posC.rpar2 =0.
            posC.P1rinf =0.0
            posC.qr=0.

        # --- This is K -> SS (see Phys. Rev. ST Accel. Beams 6, 054701 [2003])
        if (mat_num == 5):
            posC.enpar[0] = 3.9
            posC.enpar[1] = 6.2
            posC.enpar[2] = 13.
            posC.enpar[3] = 8.8
            posC.enpar[4] = 6.25
            posC.enpar[5] = 2.25
            posC.enpar[6] = 9.2
            posC.enpar[7] = 5.3
            posC.enpar[8] = 17.8
            posC.enpar[9] = 10.
            if (maxsec > 10):
                posC.enpar[10:] = 5.

            posC.pnpar[0] = 1.6
            posC.pnpar[1] = 2.
            posC.pnpar[2] = 1.8
            posC.pnpar[3] = 4.7
            posC.pnpar[4] = 1.8
            posC.pnpar[5] = 2.4
            posC.pnpar[6] = 1.8
            posC.pnpar[7] = 1.8
            posC.pnpar[8] = 2.3
            posC.pnpar[9] = 1.8
            if (maxsec > 10):
                posC.pnpar[10:] = 2.

            posC.pangsec =1.

            posC.dtspk = 55.
            posC.dtotpk = 2.382
            posC.E0tspk = 5.9e7
            posC.powts =1.25
            posC.tpar1 =0.66
            posC.tpar2 =0.8
            posC.tpar3 =0.7
            posC.tpar4 =1.
            posC.tpar5 =0.
            posC.tpar6 =0.

            posC.sige =0.0
            posC.pr =0.0
            posC.E0epk = 0.
            posC.E0w =0.
            posC.epar1 =0.0
            posC.epar2 =0.
            posC.P1einf =0.0
            posC.P1epk =0.0
            posC.powe =0.0

            posC.Ecr =0.0
            posC.rpar1 =0.0
            posC.rpar2 =0.
            posC.P1rinf =0.0
            posC.qr=0.

            posC.tpar1 = -1.
            posC.tpar2 = -1.
            posC.tpar3 = 0.
            posC.tpar4 = 0.


    def prepare_secondaries(self,itype,maxsec):

        if self.flag_pyecloud:
            raise ValueError('This method shound not be used with pyecloud')
        
        if(maxsec != posC.maxsec):
            posC.maxsec = maxsec
            posC.gchange("bincoeff")
            if pos_version >= '17p3':
                pos.init_pascal_triangle()
            else:
                pos.init_pascal_triangle(posC.nbc,posC.maxsec)

        if  itype != self.mat_number:
            self.mat_number=itype
            self.call_set_params_user(maxsec,self.mat_number)

    def getP1elast(self,E0,costheta,material,maxsec=10):
        
        if self.flag_pyecloud:
            raise ValueError('This method shound not be used with pyecloud')
        
        self.set_params(maxsec,material)
        return P1elast(E0,costheta,posC.P1einf,posC.P1epk,
                       posC.E0epk,posC.E0w,posC.powe,posC.epar1,posC.epar2)

    def getP1rediff(self,E0,costheta,material,maxsec=10):
        
        if self.flag_pyecloud:
            raise ValueError('This method shound not be used with pyecloud')
        
        self.set_params(maxsec,material)
        return P1rediff(E0,costheta,posC.Ecr,posC.rpar1,posC.rpar2,posC.qr,posC.P1rinf)

    def getdeltats(self,E0,costheta,material,maxsec=10):
              
        if self.flag_pyecloud:
            raise ValueError('This method shound not be used with pyecloud')
        
        self.set_params(maxsec,material)
        return deltats(E0,costheta,posC.dtspk,posC.E0tspk,
                       posC.powts,posC.tpar1,posC.tpar2,posC.tpar3,
                       posC.tpar4,posC.tpar5,posC.tpar6)

    def sey2(self,energy):
                
        if self.flag_pyecloud:
            raise ValueError('This method shound not be used with pyecloud')
        
        maxsec=10

        if(maxsec != posC.maxsec):
            posC.maxsec = maxsec
            posC.gchange("bincoeff")
            if pos_version >= '17p3':
                pos.init_pascal_triangle()
            else:
                pos.init_pascal_triangle(posC.nbc,posC.maxsec)

        # posC.enpar = zeros(maxsec,float64)
        # posC.pnpar = zeros(maxsec,float64)
        n=shape(energy)[0]
        s1=zeros(n,float64)
        s2=zeros(n,float64)
        s3=zeros(n,float64)
        for i in range(n):
            s1[i]=getdeltats(E0=energy[i],costheta=1.,material=2,maxsec=10)
            s2[i]=getP1rediff(E0=energy[i],costheta=1.,material=2,maxsec=10)
            s3[i]=getP1elast(E0=energy[i],costheta=1.,material=2,maxsec=10)
        return s1+s2+s3

    def getek0av(self):
        return self.htime[...],self.ek0av[...]

    def getek0max(self):
        return self.htime[...],self.ek0max[...]

    def getcosthav(self):
        return self.htime[...],self.costhav[...]

    def plek0av(self,color=black,width=1,type='solid',marker='o',marks=0,msize=1):
        if me != 0:return
        htime,ek0av=self.getek0av()
        plg(ek0av,htime,color=color,width=width,type=type,marker=marker,marks=marks,msize=msize)
        ptitles('','time (s)','ek0av (eV)')

    def plek0max(self,color=black,width=1,type='solid',marker='o',marks=0,msize=1):
        if me != 0:return
        htime,ek0max=self.getek0max()
        plg(ek0max,htime,color=color,width=width,type=type,marker=marker,marks=marks,msize=msize)
        ptitles('','time (s)','ek0max (eV)')

    def plcosthav(self,color=black,width=1,type='solid',marker='o',marks=0,msize=1):
        if me != 0:return
        htime,costhav=self.getcosthav()
        plg(costhav,htime,color=color,width=width,type=type,marker=marker,marks=marks,msize=msize)
        ptitles('','time (s)','costhav')

class PhotoElectrons:
    """Class for generating photo-electrons
    - posinst_file: name of Posinst input file
    - xfloor      : photo-electrons generated by Posinst that have x<xfloor will be forced to x=xfloor
    - xceiling    : photo-electrons generated by Posinst that have x>xceiling will be forced to x=xceiling
    - yfloor      : photo-electrons generated by Posinst that have y<yfloor will be forced to y=yfloor
    - yceiling    : photo-electrons generated by Posinst that have y>xceiling will be forced to y=yceiling
    - nz          : number of longitudinal slices (default=100)
    - l_xmirror   : turns mirroring of emitted photo-electrons with regard to x-axis on/off
    - l_verbose   : sets verbosity (default=0).
    """
    def __init__(self,posinst_file=None,xfloor=None,xceiling=None,yfloor=None,yceiling=None,
                 nz=100,l_xmirror=0,l_switchyz=0,l_verbose=0,lcallscrapercontrollers=0):
        self.totalcount = 0
        self.xfloor=xfloor
        self.xceiling=xceiling
        self.yfloor=yfloor
        self.yceiling=yceiling
        self.nz=nz
        self.l_xmirror=l_xmirror
        self.l_switchyz=l_switchyz
        self.l_verbose=l_verbose
        self.inter={}
        self.npmax={}
        self.nps={}
        self.x={}
        self.y={}
        self.z={}
        self.vx={}
        self.vy={}
        self.vz={}
        self.gi={}
        self.pid={}
        self.Lambda=0.
        self.lcallscrapercontrollers=lcallscrapercontrollers
        if posinst_file is not None:init_posinst_for_warp(posinst_file)
        self.install()

    def add(self,incident_species=None,emitted_species=None):
        isinc=incident_species
        issec=[]
        if isinc not in self.inter:
            self.inter[isinc]={}
            for key in ['incident_species','emitted_species']:
                self.inter[isinc][key]=[]
            self.inter[isinc]['incident_species']=incident_species
        self.inter[isinc]['emitted_species'] = emitted_species
        js=emitted_species.jslist[0]
        if js not in self.x:
            self.nps[js]=0
            self.npmax[js]=4096
            self.allocate_temps(js)

    def allocate_temps(self,js):
        self.x[js]=fzeros(self.npmax[js],'d')
        self.y[js]=fzeros(self.npmax[js],'d')
        self.z[js]=fzeros(self.npmax[js],'d')
        self.vx[js]=fzeros(self.npmax[js],'d')
        self.vy[js]=fzeros(self.npmax[js],'d')
        self.vz[js]=fzeros(self.npmax[js],'d')
        self.gi[js]=fzeros(self.npmax[js],'d')
        if top.wpid>0:
            self.pid[js]=fzeros([self.npmax[js],top.npid],'d')

    def install(self):
        if not isinstalleduserinjection(self.generate):
            installuserinjection(self.generate)

    def addpart(self,nn,x,y,z,vx,vy,vz,gi,js,weight=None):
        if self.nps[js]+nn>self.npmax[js]:self.flushpart(js)
        if self.nps[js]+nn>self.npmax[js]:
            self.npmax[js] = nint(nn*1.2)
            self.allocate_temps(js)
        il=self.nps[js]
        iu=il+nn
        self.x[js][il:iu]=x
        self.y[js][il:iu]=y
        self.z[js][il:iu]=z
        self.vx[js][il:iu]=vx
        self.vy[js][il:iu]=vy
        self.vz[js][il:iu]=vz
        self.gi[js][il:iu]=gi
        if weight is not None:self.pid[js][il:iu,top.wpid-1]=weight
        self.nps[js]+=nn

    def flushpart(self,js):
        if self.nps[js]>0:
            nn=self.nps[js]
            self.totalcount += nn
            if top.wpid==0:
                addparticles(x=self.x[js][:nn],
                             y=self.y[js][:nn],
                             z=self.z[js][:nn],
                             vx=self.vx[js][:nn],
                             vy=self.vy[js][:nn],
                             vz=self.vz[js][:nn],
                             gi=self.gi[js][:nn],
                             js=js,
                             lmomentum=true,
                             lallindomain=true)
            else:
                addparticles(x=self.x[js][:nn],
                             y=self.y[js][:nn],
                             z=self.z[js][:nn],
                             vx=self.vx[js][:nn],
                             vy=self.vy[js][:nn],
                             vz=self.vz[js][:nn],
                             gi=self.gi[js][:nn],
                             pid=self.pid[js][:nn,:],
                             js=js,
                             lmomentum=true,
                             lallindomain=true)
            self.nps[js]=0

    def generate(self):
        for ints in self.inter:
            incident_species=self.inter[ints]['incident_species']
            emitted_species=self.inter[incident_species]['emitted_species']
            if not isinstance(self.Lambda,ndarray):
                self.nz=0
                if self.l_switchyz:
                    self.ymin=top.ypminlocal
                    self.ymax=top.ypmaxlocal
                    self.dy=(top.ypmaxlocal-top.ypminlocal)
                else:
                    self.zmin=top.zpminlocal
                    self.zmax=top.zpmaxlocal
                    self.dz=(top.zpmaxlocal-top.zpminlocal)
            else:
                if incident_species is not None:
                    if self.l_switchyz:
                        self.ymin=min(incident_species.gety())
                        self.ymax=max(incident_species.gety())
                        self.dy=(ymax-ymin)/self.nz
                        self.Lambda = sum(sum(incident_species.get_density(nx=2,
                                                                           nz=2,
                                                                           ny=self.nz,
                                                                           ymin=self.ymin,
                                                                           ymax=self.ymax,
                                                                           l_minmax_grid=false,
                                                                           l_dividebyvolume=false,
                                                                           charge=1),2),0)
                    else:
                        self.zmin=min(incident_species.getz())
                        self.zmax=max(incident_species.getz())
                        self.dz=(self.zmax-self.zmin)/self.nz
                        self.Lambda = sum(sum(incident_species.get_density(nx=2,
                                                                           ny=2,
                                                                           nz=self.nz,
                                                                           zmin=self.zmin,
                                                                           zmax=self.zmax,
                                                                           l_minmax_grid=false,
                                                                           l_dividebyvolume=false,
                                                                           charge=1),0),0)
            weightemit=top.pgroup.sw[emitted_species.jslist[0]]*abs(top.pgroup.sq[emitted_species.jslist[0]])
            # for i in range(self.nz+1):
            for i in range(max(1,self.nz)):
                if self.nz<1:
                    rhel = self.Lambda*posC.queffp*posC.photpbppm*clight*top.dt/weightemit
                else:
                    rhel = self.Lambda[i]*posC.queffp*posC.photpbppm*clight*top.dt*self.dz/weightemit
                # rhel*=posC.slength
                # --- rhel is the number of photoelectrons created at each timestep
                # --- queffp is the quantum efficiency (photoelectrons produced per
                # --- photon) Miguel says queffp is between 0.1 and 1.0. Real
                # --- value not known.

                n=int(rhel)
                if ranf()<rhel-n:n+=1  # randomly add one electrons based on rhel fractional part
                if self.l_verbose:print ' *** i,rhel,nemit= ',i,rhel,n
                if n==0:continue
                posC.nphel[0]=n   # tells Posinst to emit n photoelectrons
                gen_photoelectrons(1) # number of beam slice in POSINST =1. Use only 1.
                if self.l_verbose:print 'nlast',posC.nlast,"nphel=",posC.nphel[0]

                if self.l_xmirror:
                    # put photons on both sides of the vacuum chamber
                    xran = ranf(posC.x[:posC.nlast])
                    xran = where(xran>0.5,1.,-1.)
                    posC.x[:posC.nlast] = posC.x[:posC.nlast]*xran
                    posC.vgx[:posC.nlast] = posC.vgx[:posC.nlast]*xran

                if self.l_verbose:print "min and max of photoelectrons=",min((posC.z[:posC.nlast]/posC.slength)*self.dz+i*self.dz),\
                                                                         max((posC.z[:posC.nlast]/posC.slength)*self.dz+i*self.dz)
                ns = posC.nlast
                js_new=emitted_species.jslist[0]
                x = posC.x[:posC.nlast]
                y = posC.y[:posC.nlast]
                z = posC.z[:posC.nlast]
                ux = posC.vgx[:posC.nlast]
                uy = posC.vgy[:posC.nlast]
                uz = posC.vgz[:posC.nlast]
                if self.xfloor is not None:
                    x=where(x>self.xfloor,x,self.xfloor)
                if self.xceiling is not None:
                    x=where(x<self.xceiling,x,self.xceiling)
                if self.yfloor is not None:
                    y=where(y>self.yfloor,y,self.yfloor)
                if self.yceiling is not None:
                    y=where(y<self.yceiling,y,self.yceiling)
                dt=ranf(x)*top.dt
                usq = (ux**2 + uy**2 + uz**2)/clight**2
                gaminv = 1./sqrt(1. + usq)
                if self.l_switchyz:
                    x = x+dt*ux*gaminv
                    y = (z/posC.slength)*self.dy+i*self.dy+self.ymin
                    z = y+dt*uy*gaminv
                else:
                    x = x+dt*ux*gaminv
                    y = y+dt*uy*gaminv
                    z = (z/posC.slength)*self.dz+i*self.dz+self.zmin
                xc = logical_and(x>=top.xpminlocal,x<top.xpmaxlocal)
                yc = logical_and(y>=top.ypminlocal,y<top.ypmaxlocal)
                ii = compress(logical_and(xc,yc),arange(ns))
                x = take(x,ii)
                y = take(y,ii)
                z = take(z,ii)
                ux = take(ux,ii)
                uy = take(uy,ii)
                uz = take(uz,ii)
                gaminv = take(gaminv,ii)
                np = shape(x)[0]
                if top.wpid==0:
                    weights = None
                else:
                    weights = ones(np,'d')
                if self.l_switchyz:
                    self.addpart(np,x,y,z,ux,uz,uy,gaminv,js_new,weights)
                else:
                    self.addpart(np,x,y,z,ux,uy,uz,gaminv,js_new,weights)
                posC.nlast=0

        # --- make sure that all particles are added
        for js in self.x:
            self.flushpart(js)
        # --- Check for particle out of bounds and exchange particles among
        # --- processors if needed. A call to particleboundaries3d is made
        # --- so that particles are scraped on any user defined conductors
        # --- as well as the grid boundaries.
        # --- Set the flag so that this generate routine is not called again
        # --- recursively (since generate is normally called at the end of
        # --- the scraping).
        # --- This is all needed in case lcallscrapercontrollers is true.
        # --- Also set the flag so that the lost particles are not reset.
        self.lrecursivegenerate = 1
        top.lresetlostpart = false
        particleboundaries3d(top.pgroup,-1,self.lcallscrapercontrollers)
        top.lresetlostpart = true
        self.lrecursivegenerate = 0
