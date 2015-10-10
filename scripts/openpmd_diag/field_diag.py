"""
This file defines the class FieldDiagnostic
"""
import os
import h5py
import numpy as np
from generic_diag import OpenPMDDiagnostic
from parallel import gatherarray
from data_dict import *

class FieldDiagnostic(OpenPMDDiagnostic) :
    """
    Class that defines the field diagnostics to be done.

    Usage
    -----
    After initialization, the diagnostic is called by using the
    `write` method.
    """

    def __init__(self, period, em, top, w3d, comm_world=None, 
                 fieldtypes=["rho", "E", "B", "J"], write_dir=None, 
                 lparallel_output=False) :
        """
        Initialize the field diagnostic.

        Parameters
        ----------
        period : int
            The period of the diagnostics, in number of timesteps.
            (i.e. the diagnostics are written whenever the number
            of iterations is divisible by `period`)
            
        em : an EM3D object (as defined in em3dsolver)
            Contains the fields data and the different methods to extract it

        top : the object representing the `top` package in Warp
            Contains information on the time.

        w3d : the object representing the `w3d` package in Warp
            Contains the dimensions of the grid.

        comm_world : a communicator object
            Either an mpi4py or a pyMPI object, or None (single-proc)
            
        fieldtypes : a list of strings, optional
            The strings are either "rho", "E", "B" or "J"
            and indicate which field should be written.
            Default : all fields are written
            
        write_dir : string, optional
            The POSIX path to the directory where the results are
            to be written. If none is provided, this will be the path
            of the current working directory.

        lparallel_output : boolean, optional
            Switch to set output mode (parallel or gathering)
            If "True" : Parallel output
        """
        # General setup
        OpenPMDDiagnostic.__init__(self, period, top, w3d, comm_world,
                                   lparallel_output, write_dir)
        
        # Register the arguments
        self.em = em
        self.fieldtypes = fieldtypes

    def setup_openpmd_meshes_group( self, dset ) :
        """
        Set the attributes that are specific to the mesh path
        
        Parameter
        ---------
        dset : an h5py.Group object that contains all the mesh quantities
        """
        # Field Solver
        dset.attrs["fieldSolver"] = field_solver_dict[ self.em.stencil ]
        # Field and particle boundary
        # - 2D and Circ
        if self.em.l_2dxz:
            dset.attrs["fieldBoundary"] = np.array([
                field_boundary_dict[ self.w3d.boundxy ],
                field_boundary_dict[ self.w3d.boundxy ],
                field_boundary_dict[ self.w3d.bound0 ],
                field_boundary_dict[ self.w3d.boundnz ] ])
            dset.attrs["particleBoundary"] = np.array([
                particle_boundary_dict[ self.top.pboundxy ],
                particle_boundary_dict[ self.top.pboundxy ],
                particle_boundary_dict[ self.top.pbound0 ],
                particle_boundary_dict[ self.top.pboundnz ] ])
        # - 3D
        else:
            dset.attrs["fieldBoundary"] = np.array([
                field_boundary_dict[ self.w3d.boundxy ],
                field_boundary_dict[ self.w3d.boundxy ],
                field_boundary_dict[ self.w3d.boundxy ],
                field_boundary_dict[ self.w3d.boundxy ],
                field_boundary_dict[ self.w3d.bound0 ],
                field_boundary_dict[ self.w3d.boundnz ] ])
            dset.attrs["particleBoundary"] = np.array([
                particle_boundary_dict[ self.top.pboundxy ],
                particle_boundary_dict[ self.top.pboundxy ],
                particle_boundary_dict[ self.top.pboundxy ],
                particle_boundary_dict[ self.top.pboundxy ],
                particle_boundary_dict[ self.top.pbound0 ],
                particle_boundary_dict[ self.top.pboundnz ] ])

        # Current Smoothing
        if np.all( self.em.npass_smooth == 0 ) :
            dset.attrs["currentSmoothing"] = np.string_("none")
        else :
            dset.attrs["currentSmoothing"] = np.string_("Binomial")
            dset.attrs["currentSmoothingParameters"] = str(self.em.npass_smooth)
        # Charge correction
        dset.attrs["chargeCorrection"] = np.string_("none")
        
    def setup_openpmd_mesh_record( self, dset, quantity ) :
        """
        Sets the attributes that are specific to a mesh record
        
        Parameter
        ---------
        dset : an h5py.Dataset or h5py.Group object

        quantity : string
           The name of the record (e.g. "rho", "J", "E" or "B")
        """
        # Generic record attributes
        self.setup_openpmd_record( dset, quantity )
        
        # Geometry parameters
        # - thetaMode
        if (self.em.l_2drz==True) :
            dset.attrs['geometry']  = np.string_("thetaMode")
            dset.attrs['geometryParameters'] = \
              np.string_("m=%d;imag=+" %(self.em.circ_m + 1))
            dset.attrs['gridSpacing'] = np.array([ self.em.dx, self.em.dz ])
            dset.attrs['axisLabels'] = np.array([ 'r', 'z' ])
            dset.attrs["gridGlobalOffset"] = np.array([
                self.w3d.xmmin, self.top.zgrid + self.w3d.zmmin])
        # - 2D Cartesian
        elif (self.em.l_2dxz==True) :
            dset.attrs['geometry'] = np.string_("cartesian")
            dset.attrs['gridSpacing'] = np.array([ self.em.dx, self.em.dz ])
            dset.attrs['axisLabels'] = np.array([ 'x', 'z' ])
            dset.attrs["gridGlobalOffset"] = np.array([
                self.w3d.xmmin, self.top.zgrid + self.w3d.zmmin])
        # - 3D Cartesian
        else :
            dset.attrs['geometry'] = np.string_("cartesian")
            dset.attrs['gridSpacing'] = np.array(
                [ self.em.dx, self.em.dy, self.em.dz ])
            dset.attrs['axisLabels'] = np.array([ 'x', 'y', 'z' ])
            dset.attrs["gridGlobalOffset"] = np.array([ self.w3d.xmmin,
                        self.w3d.ymmin, self.top.zgrid + self.w3d.zmmin])
            
        # Generic attributes
        dset.attrs["dataOrder"] = np.string_("C")
        dset.attrs["gridUnitSI"] = 1.
        dset.attrs["fieldSmoothing"] = np.string_("none")

    def setup_openpmd_mesh_component( self, dset, quantity ) :
        """
        Set up the attributes of a mesh component
    
        Parameter
        ---------
        dset : an h5py.Dataset or h5py.Group object
        
        quantity : string
            The field that is being written
        """
        # Generic setup of the component
        self.setup_openpmd_component( dset )
        
        # Field positions
        if (self.em.l_2dxz==True) :
            positions = np.array([0., 0.])
        else:
            positions = np.array([0.,0.,0.])
        # Along x
        positions[0] = x_offset_dict[ quantity ]
        # Along y (3D Cartesian only)
        if (self.em.l_2dxz==False) :
            positions[1] = y_offset_dict[quantity]
        # Along z
        positions[-1] = z_offset_dict[ quantity ]

        dset.attrs['position'] = positions
        
    def write_hdf5( self, iteration ) :
        """
        Write an HDF5 file that complies with the OpenPMD standard

        Parameter
        ---------
        iteration : int
             The current iteration number of the simulation.
        """
        # Create the file

        filename = "data%08d.h5" %iteration
        fullpath = os.path.join( self.write_dir, "diags/hdf5", filename )

        # In gathering mode, only the first proc creates the file.
        if self.lparallel_output == False and self.rank == 0 :
            # Create the filename and open hdf5 file
            f = h5py.File( fullpath, mode="a" )
            self.setup_openpmd_file( f )
            this_rank_writes = True
        # In parallel mode (lparallel_output=True), all proc create the file
        elif self.lparallel_output == True :
            # Create the filename and open hdf5 file
            f = h5py.File( fullpath, mode="a", driver='mpio',
                           comm=self.comm_world)
            self.setup_openpmd_file( f )
            this_rank_writes = True
        else:
            f = None
            this_rank_writes = False

        # Setup the fields group
        if this_rank_writes :
            field_path = "/data/%d/fields/" %iteration
            field_grp = f.require_group(field_path)
            self.setup_openpmd_meshes_group( field_grp )
        else:
            field_grp = None

        # Determine the components to be written (Cartesian or cylindrical)
        if (self.em.l_2drz == True) :
            coords = ['r', 't', 'z']
        else :
            coords = ['x', 'y', 'z']
            
        # Loop over the different quantities that should be written
        for fieldtype in self.fieldtypes :
            # Scalar field
            if fieldtype == "rho" :
                self.write_dataset( field_grp, "rho", "rho", this_rank_writes )
                if this_rank_writes :
                    self.setup_openpmd_mesh_record( field_grp["rho"], "rho" )
            # Vector field
            elif fieldtype in ["E", "B", "J"] :
                for coord in coords :
                    quantity = "%s%s" %(fieldtype, coord)
                    path = "%s/%s" %(fieldtype, coord)
                    self.write_dataset( field_grp, path, quantity,
                                        this_rank_writes )
                if this_rank_writes :
                    self.setup_openpmd_mesh_record(
                        field_grp[fieldtype], fieldtype )
            else :
                raise ValueError("Invalid string in fieldtypes: %s" %fieldtype)
        
        # Close the file
        if this_rank_writes :      
            f.close()

    def write_dataset( self, field_grp, path, quantity, this_rank_writes ) :
        """
        Write a given dataset
    
        Parameters
        ----------
        field_grp : an h5py.Group object
            The group that corresponds to the path indicated in meshesPath
        
        path : string
            The relative path where to write the dataset, in field_grp

        quantity : string
            Describes which field is being written.
            (Either rho, Er, Et, Ez, Br, Bz, Bt, Jr, Jt or Jz)

        this_rank_writes : bool
            Wether this proc participates in creating the dataset
            Parallel mode (lparallel_output=True) : all proc write
            Gathering mode : only the first proc writes
        """
        # Circ case
        if (self.em.l_2drz==True):
            self.write_circ_dataset( field_grp, path, quantity,
                                     this_rank_writes )
        # 2D Cartesian case
        elif (self.em.l_2dxz==True):
            self.write_cart2d_dataset( field_grp, path, quantity,
                                     this_rank_writes )
        # 3D Cartesian case
        else:
            self.write_cart3d_dataset( field_grp, path, quantity,
                                     this_rank_writes )

        
    def write_circ_dataset( self, field_grp, path, quantity,
                            this_rank_writes ) :
        """
        Write a dataset in Circ coordinates
        
        See the docstring of write_dataset for the parameters
        """
        # Create the dataset and setup its attributes
        if this_rank_writes :
            # Shape of the data : first write the real part mode 0
            # and then the imaginary part of the mode 1
            datashape = (3, self.em.nx+1, self.em.nz+1)
            dset = field_grp.require_dataset( path, datashape, dtype='f' )
            self.setup_openpmd_mesh_component( dset, quantity )
            
        # Fill the dataset with these quantities
        # Gathering mode
        if self.lparallel_output == False :
            F, F_circ, _ = self.get_circ_dataset( quantity, lgather=True )
            if self.rank == 0:
	            # Mode m=0
    	        dset[0,:,:] = F
    	        if F_circ is not None:
        	        # Mode m=1 (real and imaginary part)
            	    dset[1,:,:] = F_circ[:,:,0].real
            	    dset[2,:,:] = F_circ[:,:,0].imag
        # Parallel mode
        else:
            F, F_circ, bounds = self.get_circ_dataset( quantity, False )
            # Mode m=0
            dset[ 0, bounds[0,0]:bounds[1,0],
                     bounds[0,1]:bounds[1,1] ] = F
            if F_circ is not None:
                # Mode m=1 (real and imaginary part)
                dset[ 1, bounds[0,0]:bounds[1,0],
                         bounds[0,1]:bounds[1,1] ] = F_circ[:,:,0].real
                dset[ 2, bounds[0,0]:bounds[1,0],
                         bounds[0,1]:bounds[1,1] ] = F_circ[:,:,0].imag


    def write_cart2d_dataset( self, field_grp, path, quantity,
                            this_rank_writes ) :
        """
        Write a dataset in Cartesian coordinates
        
        See the docstring of write_dataset for the parameters
        """
        # Create the dataset and setup its attributes
        if this_rank_writes :
            datashape = (self.em.nx+1, self.em.nz+1)
            dset = field_grp.require_dataset( path, datashape, dtype='f' )
            self.setup_openpmd_mesh_component( dset, quantity )
            
        # Fill the dataset with these quantities
        # Gathering mode
        if self.lparallel_output == False :
            F, _ = self.get_cart_dataset( quantity, True )
            if self.rank == 0:
    	        dset[:,:] = F
        # Parallel mode
        else:
            F, bounds = self.get_cart_dataset( quantity, False )
            dset[ bounds[0,0]:bounds[1,0],
                    bounds[0,1]:bounds[1,1] ] = F

    def write_cart3d_dataset( self, field_grp, path, quantity,
                            this_rank_writes ) :
        """
        Write a dataset in Cartesian coordinates
        
        See the docstring of write_dataset for the parameters
        """
        # Create the dataset and setup its attributes
        if this_rank_writes :
            datashape = (self.em.nx+1, self.em.ny+1, self.em.nz+1)
            dset = field_grp.require_dataset( path, datashape, dtype='f' )
            self.setup_openpmd_mesh_component( dset, quantity )
            
        # Fill the dataset with these quantities
        # Gathering mode
        if self.lparallel_output == False :
            F, _ = self.get_cart_dataset( quantity, True )
            if self.rank == 0:
    	        dset[:,:,:] = F
        # Parallel mode
        else:
            F, bounds = self.get_cart_dataset( quantity, False )
            dset[ bounds[0,0]:bounds[1,0],
                  bounds[0,1]:bounds[1,1],
                  bounds[0,2]:bounds[1,2] ] = F
                     
    def get_circ_dataset( self, quantity, lgather) :
        """
        Get a given quantity in Circ coordinates

        Parameters
        ----------
        quantity : string
            Describes which field is being written.
            (Either rho, Er, Et, Ez, Br, Bz, Bt, Jr, Jt or Jz)

        lgather : boolean
            Defines if data is gathered on me (process) = 0
            If "False": No gathering is done
        """
        F_circ = None
        em = self.em
        
        # Treat the currents in a special way
        if quantity in ['Jr', 'Jt', 'Jz'] :
            # Get the array index that corresponds to that component
            i = circ_dict_Jindex[ quantity ]
            # Extract mode 0
            F = em.getarray( em.fields.J[:,:,:,i] )
            # Extract higher modes
            if em.circ_m > 0 :
                F_circ = em.getarray_circ( em.fields.J_circ[:,:,i,:] )
                
        # Treat the fields E, B, rho in a more systematic way
        elif quantity in ['Er', 'Et', 'Ez', 'Br', 'Bt', 'Bz', 'rho' ] :
            # Get the field name in Warp
            field_name = circ_dict_quantity[ quantity ]
            # Extract mode 0
            field_array = getattr( em.fields, field_name )
            F = em.getarray( field_array )
            # Extract higher modes
            if em.circ_m > 0 :
                field_array = getattr( em.fields, field_name + '_circ')
                F_circ = em.getarray_circ( field_array )

        # Gather array if lgather = True 
        # (Mutli-proc operations using gather)
        # Only done in non-parallel case
        if lgather == True:
            F = em.gatherarray( F )
            if em.circ_m > 0 :
                F_circ = em.gatherarray( F_circ )
        
        # Get global positions (indices) of local domain
        # Only needed for parallel output
        if lgather == False :
            nx, nz = np.shape(F)
            bounds = np.zeros([2,2], dtype = np.int)
            bounds[0,0] = int(round((em.block.xmin - em.xmmin) / em.dx))
            bounds[1,0] = bounds[0,0] + nx
            bounds[0,1] = int(round((em.block.zmin - em.zmmin) / em.dz))
            bounds[1,1] = bounds[0,1] + nz
        else :
            bounds = None

        return( F, F_circ, bounds )

        
    def get_cart_dataset( self, quantity, lgather) :
        """
        Get a given quantity in Cartesian coordinates

        Parameters
        ----------
        quantity : string
            Describes which field is being written.
            (Either rho, Er, Et, Ez, Br, Bz, Bt, Jr, Jt or Jz)

        lgather : boolean
            Defines if data is gathered on me (process) = 0
            If "False": No gathering is done
        """
        em = self.em
        
        # Treat the currents in a special way
        if quantity in ['Jx', 'Jy', 'Jz'] :
            # Get the array index that corresponds to that component
            i = cart_dict_Jindex[ quantity ]
            F = em.getarray( em.fields.J[:,:,:,i] )

        # Treat the fields E, B, rho in a more systematic way
        elif quantity in ['Ex', 'Ey', 'Ez', 'Bx', 'By', 'Bz', 'rho' ] :
            # Get the field name in Warp
            field_name = cart_dict_quantity[ quantity ]
            field_array = getattr( em.fields, field_name )
            F = em.getarray( field_array )

        # Gather array if lgather = True 
        # (Mutli-proc operations using gather)
        # Only done in non-parallel case
        if lgather == True:
            F = em.gatherarray( F )
        
        # Get global positions (indices) of local domain
        # Only needed for parallel output
        if lgather == False :
            if F.ndim == 2:
                nx, nz = np.shape(F)
                bounds = np.zeros([2,2], dtype = np.int)
                bounds[0,0] = int(round((em.block.xmin - em.xmmin) / em.dx))
                bounds[1,0] = bounds[0,0] + nx
                bounds[0,1] = int(round((em.block.zmin - em.zmmin) / em.dz))
                bounds[1,1] = bounds[0,1] + nz
            elif F.ndim == 3:
                nx, ny, nz = np.shape(F)
                bounds = np.zeros([2,3], dtype = np.int)
                bounds[0,0] = int(round((em.block.xmin - em.xmmin) / em.dx))
                bounds[1,0] = bounds[0,0] + nx
                bounds[0,1] = int(round((em.block.ymin - em.ymmin) / em.dy))
                bounds[1,1] = bounds[0,1] + ny
                bounds[0,2] = int(round((em.block.zmin - em.zmmin) / em.dz))
                bounds[1,2] = bounds[0,2] + nz
        else :
            bounds = None

        return( F, bounds )