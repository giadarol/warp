"""
This file defines the class ProbeParticleDiagnostic

This diagnostic records the particles that go through a given plane
in time and saves them in a single openPMD file (while doing caching, in
order not to access the disk at every single timestep)

Note: this is a hack the openPMD standard in that it stores particles from
different time in the same openPMD file. (The particles have an addiditional
attribute `t` which is stored in the openPMD file.)
"""
import os
import numpy as np
import time
from scipy.constants import c
from particle_diag import ParticleDiagnostic
from parallel import gatherarray

class ProbeParticleDiagnostic(ParticleDiagnostic):
    """
    Class that writes the particles that go cross a given plane, in
    the direction given by `plane_normal_vector`
    (The particles that cross the plane in the other direction are not saved.)

    Usage
    -----
    After initialization, the diagnostic is called by using
    the 'write' method.
    """
    def __init__(self, plane_position, plane_normal_vector,
                 period, top, w3d, comm_world=None,
                 particle_data=["position", "momentum", "weighting", "t"],
                 select=None, write_dir=None, lparallel_output=False,
                 species={"electrons": None}):
        """
        Initialize diagnostics that retrieve the particles crossing a given
        plane.

        Parameters
        ----------
        plane_position: a list of 3 floats (in meters)
            The position (in x, y, z) of one of the points of the plane

        plane_normal_vector: a list of 3 floats
            The coordinates (in x, y, z) of one of the vectors of the plane

        period: int
            Number of iterations for which the data is accumulated in memory,
            before finally writing it to the disk.

        See the documentation of ParticleDiagnostic for the other parameters
        """
        # Do not leave write_dir as None, as this may conflict with
        # the default directory ('./diags')
        if write_dir is None:
            write_dir = 'probe_diags'

        # Initialize Particle diagnostic normal attributes
        ParticleDiagnostic.__init__(self, period, top, w3d, comm_world,
            species=species, particle_data=particle_data, select=select,
            write_dir=write_dir, lparallel_output=lparallel_output)

        # Initialize proper helper objects
        self.particle_storer = ParticleStorer( top.dt, self.write_dir,
            self.species_dict, self.lparallel_output, self.rank )
        self.particle_catcher = ParticleCatcher( top, plane_position,
                                                plane_normal_vector )
        self.particle_catcher.allocate_previous_instant()

        # Initialize a corresponding empty file
        if self.lparallel_output == False and self.rank == 0:
            self.create_file_empty_particles(
                self.particle_storer.filename, 0, 0, self.top.dt)

    def write( self ):
        """
        Redefines the method write of the parent class ParticleDiagnostic

        Should be registered with installafterstep in Warp
        """
        # At each timestep, store new particles in memory buffers
        self.store_new_particles()

        # Every self.period, write the buffered slices to disk
        if self.top.it % self.period == 0:
            self.flush_to_disk()

    def store_new_particles( self ):
        """
        Store new of the particles in the memory buffers of the
        particle storer

        The stored particles are also selected in accordance with
        the selection rules provided as argument to ProbeParticleDiagnostic
        """
        # Loop through the particle species and register the
        # particle arrays in the particle storer object (buffering)
        for species_name, species in self.species_dict.iteritems():

            slice_array = self.particle_catcher.extract_slice(
                        species, self.select )
            self.particle_storer.register_slice( slice_array, species_name )

    def flush_to_disk(self):
        """
        Writes the buffered slices of particles to the disk. Erase the
        buffered slices of the ParticleStorer object

        Notice: In parallel version, data are gathered to proc 0
        before being saved to disk
        """
        # Compact the successive slices that have been buffered
        # over time into a single array
        for species_name in self.species_dict:

            particle_array = self.particle_storer.compact_slices(species_name)

            if self.comm_world is not None:
                # In MPI mode: gather an array containing the number
                # of particles on each process
                n_rank = self.comm_world.allgather(np.shape(particle_array)[1])

                # Note that gatherarray routine in parallel.py only works
                # with 1D array. Here we flatten the 2D particle arrays
                # before gathering.
                g_curr = gatherarray(particle_array.flatten(), root=0,
                    comm=self.comm_world)

                if self.rank == 0:
                    # Get the number of quantities
                    nquant = np.shape(
                        self.particle_catcher.particle_to_index.keys())[0]

                    # Prepare an empty array for reshaping purposes. The
                    # final shape of the array is (8, total_num_particles)
                    p_array = np.empty((nquant, 0))

                    # Index needed in reshaping process
                    n_ind = 0

                    # Loop over all the processors, if the processor
                    # contains particles, we reshape the gathered_array
                    # and reconstruct by concatenation
                    for i in xrange(self.top.nprocs):

                        if n_rank[i] != 0:
                            p_array = np.concatenate((p_array, np.reshape(
                                g_curr[n_ind:n_ind+nquant*n_rank[i]],
                                (nquant,n_rank[i]))),axis=1)

                            # Update the index
                            n_ind += nquant*n_rank[i]

            else:
                p_array = particle_array

            # Write this array to disk (if this self.particle_storer has new slices)
            if self.rank == 0 and p_array.size:
                self.write_slices(p_array, species_name, self.particle_storer,
                    self.particle_catcher.particle_to_index)

        # Erase the buffers
        self.particle_storer.buffered_slices[species_name] = []

    def write_probe_dataset(self, species_grp, path, data, quantity):
        """
        Writes each quantity of the buffered dataset to the disk, the
        final step of the writing
        """
        dset = species_grp[path]
        index = dset.shape[0]

        # Resize the h5py dataset
        dset.resize(index+len(data), axis=0)

        # Write the data to the dataset at correct indices
        dset[index:] = data

    def write_slices( self, particle_array, species_name, particle_storer, p2i ):
        """
        Write the slices of the different species to an openPMD file

        Parameters
        ----------
        particle_array: array of reals
            Array of shape (8, num_part)

        species_name: String
            A String that acts as the key for the buffered_slices dictionary

        particle_storer: a ParticleStorer object

        p2i: dict
            Dictionary of correspondance between the particle quantities
            and the integer index in the particle_array
        """
        # Open the file without parallel I/O in this implementation
        f = self.open_file( particle_storer.filename, parallel_open=False )
        particle_path = "/data/%d/particles/%s" %(particle_storer.iteration,
                                                    species_name)
        species_grp = f[particle_path]

        # Loop over the different quantities that should be written
        for particle_var in self.particle_data:

            if particle_var == "position":
                for coord in ["x","y","z"]:
                    quantity= coord
                    path = "%s/%s" %(particle_var, quantity)
                    data = particle_array[ p2i[ quantity ] ]
                    self.write_probe_dataset(
                            species_grp, path, data, quantity)

            elif particle_var == "momentum":
                for coord in ["x","y","z"]:
                    quantity= "u%s" %coord
                    path = "%s/%s" %(particle_var,coord)
                    data = particle_array[ p2i[ quantity ] ]
                    self.write_probe_dataset(
                            species_grp, path, data, quantity)

            elif particle_var == "t":
               quantity= "t"
               path = "t"
               data = particle_array[ p2i[ quantity ] ]
               self.write_probe_dataset(species_grp, path, data, quantity)

            elif particle_var == "weighting":
               quantity= "w"
               path = "weighting"
               data = particle_array[ p2i[ quantity ] ]
               self.write_probe_dataset(species_grp, path, data, quantity)

        # Close the file
        f.close()

class ParticleStorer:
    """
    Class that stores data relative to the particles that are crossing the plane
    """
    def __init__(self, dt, write_dir, species_dict, lparallel_output, rank):
        """
        Initialize a ParticleStorer object

        Parameters
        ----------
        write_dir: string
            Absolute path to the directory where the data for
            this snapshot is to be written

        species_dict: dict
            Contains all the species name of the species object
            (inherited from Warp)
        """
        # Deduce the name of the filename where this snapshot writes
        if lparallel_output == False and rank == 0:
            self.filename = os.path.join( write_dir, 'hdf5/data%08d.h5' %0)
        self.iteration = 0
        self.dt = dt

        # Prepare buffered slices
        self.buffered_slices = {}
        for species_name in species_dict:
            self.buffered_slices[species_name] = []

    def register_slice(self, slice_array, species):
        """
        Store the slice of particles represented by slice_array

        Parameters
        ----------
        slice_array: array of reals
            An array of packed fields that corresponds to one slice,
            as given by the ParticleCatcher object

        species: String, key of the species_dict
            Act as the key for the buffered_slices dictionary
        """
        # Store the values
        self.buffered_slices[species].append(slice_array)

    def compact_slices(self, species):
        """
        Compact the successive slices that have been buffered
        over time into a single array.

        Parameters
        ----------
        species: String, key of the species_dict
            Act as the key for the buffered_slices dictionary

        Returns
        -------
        paticle_array: an array of reals of shape (9, numPart)
        regardless of the dimension

        Returns None if the slices are empty
        """
        if self.buffered_slices[species] != []:
            particle_array = np.concatenate(
                self.buffered_slices[species], axis=1)
        else:
            particle_array = np.empty((8,0))

        return particle_array

class ParticleCatcher:
    """
    Class that extracts, interpolates and gathers particles
    """
    def __init__(self, top, plane_position, plane_normal_vector ):
        """
        Initialize the ParticleCatcher object

        Parameters
        ----------
        plane_position: a list of 3 floats (in meters)
            The position (in x, y, z) of one of the points of the plane

        plane_normal_vector: a list of 3 floats
            The coordinates (in x, y, z) of one of the vectors of the plane

        top: WARP object
        """
        # Some attributes neccessary for particle selections
        self.plane_position = plane_position
        self.plane_normal_vector = plane_normal_vector
        self.top = top

        # Create a dictionary that contains the correspondance
        # between the particles quantity and array index
        self.particle_to_index = {'x':0, 'y':1, 'z':2, 'ux':3,
                'uy':4, 'uz':5, 'w':6, 't':7}

    def get_particle_slice( self, species ):
        """
        Select the particles for the current slice, and extract their
        positions and momenta at the current and previous timestep

        Parameters
        ----------
        species: a Species object of Warp
            Contains the particle data from which one slice will be extracted

        Returns
        -------
        num_part: int
            Number of selected particles
        """
        # Quantities at current time step
        current_x = self.get_quantity( species, "x" )
        current_y = self.get_quantity( species, "y" )
        current_z = self.get_quantity( species, "z" )
        current_ux = self.get_quantity( species, "ux" )
        current_uy = self.get_quantity( species, "uy" )
        current_uz = self.get_quantity( species, "uz" )
        weights = self.get_quantity( species, "w" )

        # Quantities at previous time step
        previous_x = self.get_quantity( species, "x", l_prev=True )
        previous_y = self.get_quantity( species, "y", l_prev=True )
        previous_z = self.get_quantity( species, "z", l_prev=True )
        previous_ux = self.get_quantity( species, "ux", l_prev=True )
        previous_uy = self.get_quantity( species, "uy", l_prev=True )
        previous_uz = self.get_quantity( species, "uz", l_prev=True )

        # A particle array for mapping purposes
        particle_indices = np.arange( len(current_z) )

        # For this snapshot:
        # - check if the particles where before the plane at the previous timestep
        # - check if the particle are beyond the plane at the current timestep
        r = self.plane_position
        n = self.plane_normal_vector
        previous_position_relative_to_plane = \
              n[0]*(previous_x - r[0]) \
            + n[1]*(previous_y - r[1]) \
            + n[2]*(previous_z - r[2])
        current_position_relative_to_plane = \
              n[0]*(current_x - r[0]) \
            + n[1]*(current_y - r[1]) \
            + n[2]*(current_z - r[2])
        selected_indices = np.compress(
            (previous_position_relative_to_plane <= 0 ) &
            (current_position_relative_to_plane > 0 )  , particle_indices)

        num_part = np.shape(selected_indices)[0]

        ## Select the particle quantities that satisfy the
        ## aforementioned condition
        self.mass = species.mass
        self.w_captured = np.take(weights, selected_indices)

        current_x = np.take(current_x, selected_indices)
        current_y = np.take(current_y, selected_indices)
        current_z = np.take(current_z, selected_indices)
        current_ux = np.take(current_ux, selected_indices)
        current_uy = np.take(current_uy, selected_indices)
        current_uz = np.take(current_uz, selected_indices)
        current_position_relative_to_plane = np.take(
            current_position_relative_to_plane, selected_indices )

        previous_x = np.take(previous_x, selected_indices)
        previous_y = np.take(previous_y, selected_indices)
        previous_z = np.take(previous_z, selected_indices)
        previous_ux = np.take(previous_ux, selected_indices)
        previous_uy = np.take(previous_uy, selected_indices)
        previous_uz = np.take(previous_uz, selected_indices)
        previous_position_relative_to_plane = np.take(
            previous_position_relative_to_plane, selected_indices )

        # Interpolate particle quantity to the time when they cross the plane
        norm_factor = 1 / ( np.abs(previous_position_relative_to_plane) \
                + current_position_relative_to_plane )
        interp_current = np.abs(previous_position_relative_to_plane) * norm_factor
        interp_previous = current_position_relative_to_plane * norm_factor

        self.t_captured = interp_current * self.top.time + \
                            interp_previous * (self.top.time - self.top.dt)
        self.x_captured = interp_current * current_x + \
                            interp_previous * previous_x
        self.y_captured = interp_current * current_y + \
                            interp_previous * previous_y
        self.z_captured = interp_current * current_z + \
                            interp_previous * previous_z
        self.ux_captured = interp_current * current_ux + \
                            interp_previous * previous_ux
        self.uy_captured = interp_current * current_uy + \
                            interp_previous * previous_uy
        self.uz_captured = interp_current * current_uz + \
                            interp_previous * previous_uz

        return( num_part )

    def gather_array(self, quantity):
        """
        Gather the quantity arrays and normalize the momenta
        Parameters
        ----------
        quantity: String
            Quantity of the particles that is wished to be gathered

        Returns
        -------
        ar: array of reals
            An array of gathered particle's quantity
        """
        ar = np.zeros(np.shape(self.x_captured)[0])

        if quantity == "x":
            ar = np.array(self.x_captured)
        elif quantity == "y":
            ar = np.array(self.y_captured)
        elif quantity == "z":
            ar = np.array(self.z_captured)
        elif quantity == "ux":
            ar = np.array(self.ux_captured)
        elif quantity == "uy":
            ar = np.array(self.uy_captured)
        elif quantity == "uz":
            ar = np.array(self.uz_captured)
        elif quantity == "w":
            ar = np.array(self.w_captured)
        elif quantity == "t":
            ar = np.array(self.t_captured)
        return ar

    def extract_slice(self, species, select ):
        """
        Extract a slice of the particles

        If select is present, extract only particles that satisfy the criteria

        Parameters
        ----------
        species: a Species object of Warp
            Contains the particle data from which one slice will be extracted

        select: dict
            A set of rules defined by the users in selecting the particles
            Ex: {"uz": [50, 100]} for particles which have normalized
            values between 50 and 100

        Returns
        -------
        slice_array: An array of reals of shape (8, num_part)
            An array that packs together the different particle quantities
            (x, y, z, ux, uy, uz, weight, t)
        """
        # Declare an attribute for convenience
        p2i = self.particle_to_index

        # Get the particles
        num_part = self.get_particle_slice( species )
        slice_array = np.empty((np.shape(p2i.keys())[0], num_part,))

        # Get the particle quantities
        for quantity in self.particle_to_index.keys():
            # Here typical values for 'quantity' are e.g. 'z', 'ux', 'gamma'
            # you should just gather array locally
            slice_array[ p2i[quantity], ... ] = self.gather_array(quantity)

        # Choose the particles based on the select criteria defined by the
        # users.
        if (select is not None) and slice_array.size:
            select_array = self.apply_selection(select, slice_array)
            row, column =  np.where(select_array==True)
            temp_slice_array = slice_array[row,column]

            # Temp_slice_array is a 1D numpy array, we reshape it so that it
            # has the same size as slice_array
            slice_array = np.reshape(
                temp_slice_array,(np.shape(p2i.keys())[0],-1))

        # Multiplying momenta by the species mass to make them unitless
        for quantity in ["ux", "uy", "uz"]:
            slice_array[p2i[quantity]] *= species.mass

        return slice_array

    def get_quantity(self, species, quantity, l_prev=False):
        """
        Get a given particle quantity

        Parameters
        ----------
        species: a Species object of Warp
            Contains the particle data from which the quantity is extracted

        quantity: string
            Describes which quantity is queried
            Either "x", "y", "z", "ux", "uy", "uz", "w"

        l_prev: boolean
            If True, then return the quantities of the previous timestep;
            else return quantities of the current timestep
        """
        # Extract the chosen quantities

        # At current timestep
        if not(l_prev):
            if quantity == "x":
                quantity_array = species.getx( gather=False )
            elif quantity == "y":
                quantity_array = species.gety( gather=False )
            elif quantity == "z":
                quantity_array = species.getz( gather=False )
            elif quantity == "ux":
                quantity_array = species.getux( gather=False )
            elif quantity == "uy":
                quantity_array = species.getuy( gather=False )
            elif quantity == "uz":
                quantity_array = species.getuz( gather=False )
            elif quantity == "w":
                quantity_array = species.getweights( gather=False )

        # Or at previous timestep
        else:
            if quantity == "x":
                quantity_array = species.getpid( id=self.top.xoldpid-1,
                    gather=0, bcast=0)
            elif quantity == "y":
                quantity_array = species.getpid( id=self.top.yoldpid-1,
                    gather=0, bcast=0)
            elif quantity == "z":
                quantity_array = species.getpid( id=self.top.zoldpid-1,
                    gather=0, bcast=0)
            elif quantity == "ux":
                quantity_array = species.getpid( id=self.top.uxoldpid-1,
                    gather=0, bcast=0)
            elif quantity == "uy":
                quantity_array = species.getpid( id=self.top.uyoldpid-1,
                    gather=0, bcast=0)
            elif quantity == "uz":
                quantity_array = species.getpid( id=self.top.uzoldpid-1,
                    gather=0, bcast=0)

        return( quantity_array )

    def allocate_previous_instant(self):
        """
        Allocate the top.'quantity'oldpid arrays. This is used to store
        the previous values of the quantities.
        """
        if not self.top.xoldpid:
            self.top.xoldpid = self.top.nextpid()
        if not self.top.yoldpid:
            self.top.yoldpid = self.top.nextpid()
        if not self.top.zoldpid:
            self.top.zoldpid = self.top.nextpid()
        if not self.top.uxoldpid:
            self.top.uxoldpid = self.top.nextpid()
        if not self.top.uyoldpid:
            self.top.uyoldpid = self.top.nextpid()
        if not self.top.uzoldpid:
            self.top.uzoldpid = self.top.nextpid()

    def apply_selection(self, select, slice_array):
        """
        Apply the rules of self.select to determine which
        particles should be written

        Parameters
        ----------
        select: a dictionary that defines all selection rules based
        on the quantities

        Returns
        -------
        A 1d array of the same shape as that particle array
        containing True for the particles that satify all
        the rules of self.select
        """
        p2i = self.particle_to_index

        # Initialize an array filled with True
        select_array = np.ones( np.shape(slice_array), dtype='bool' )

        # Apply the rules successively
        # Go through the quantities on which a rule applies
        for quantity in select.keys():
            # Lower bound
            if select[quantity][0] is not None:
                select_array = np.logical_and(
                    slice_array[p2i[quantity]] >\
                     select[quantity][0], select_array )
            # Upper bound
            if select[quantity][1] is not None:
                select_array = np.logical_and(
                    slice_array[p2i[quantity]] <\
                    select[quantity][1], select_array )

        return select_array
