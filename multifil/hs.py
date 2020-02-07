#!/usr/bin/env python
# encoding: utf-8
"""
hs.py - A half-sarcomere model with multiple thick and thin filaments

Created by Dave Williams on 2009-12-31.
"""

import multiprocessing as mp
import sys
import time

import numpy as np

from . import af
from . import mf
from . import ti

# Parameter Standards
HS_MF_K = 0
HS_AF_K = 1
HS_TI_A = 0
HS_TI_B = 1


class hs:
    """The half-sarcomere and ways to manage it"""

    def __init__(self, lattice_spacing=None, z_line=None, poisson=None,
                 pCa=None, timestep_len=1, time_dependence=None, starts=None, **kwargs):
        """ Create the data structure that is the half-sarcomere model

        Parameters:
            lattice_spacing: the surface-to-surface distance (14.0)
            z_line: the length of the half-sarcomere (1250)
            poisson: poisson ratio obeyed when z-line changes. Significant
                values are:
                    * 0.5 - constant volume
                    * 0.0 - constant lattice spacing, default value
                    * any negative value - auxetic lattice spacing
            pCa: pCa, controlling tropomyosin movement and thus binding site
                availability, positive by convention (4.0)
            timestep_len: how many ms per timestep (1)
            time_dependence: a dictionary to override the initial lattice
                spacing, sarcomere length, and pCa at each
                timestep. Each key may contain a list of the values, to be
                iterated over as timesteps proceed. The first entry in these
                lists will override passed initial values. The valid keys
                time_dependence can control are:
                    * "lattice_spacing"
                    * "z_line"
                    * "pCa"
            starts: starting polymer/orientation for thin/thick filaments in
                form ((rand(0,25), ...), (rand(0,3), ...))
        Returns:
            None

        This is the organizational basis that the rest of the model, and
        classes representing the other structures therein will use when
        running. It contains the following properties:

        ## Half-sarcomere properties: these are properties that can be
        interpreted as belonging to the overall model, not to any thick or
        thin filament.

        lattice_spacing:
            the face to face lattice spacing for the whole model
        m_line:
            x axis location of the m line
        h_line:
            x axis location of the h line
        hiding_line:
            x axis location below which actin sites are hidden by actin
            overlap (crossing through the m-line from adjacent half sarc)

        ## Thick Filament Properties: each is a tuple of thick filaments
        (filament_0, filament_1, filament_2, filament_3) where each
        filament_x is giving the actual properties of that particular
        filament.

        thick_location:
            each tuple location is a list of x axis locations
        thick_crowns:
            each tuple location is a tuple of links to crown instances
        thick_link:
            each tuple location is a list consisting of three (one for each
            myosin head in the crown) of either None or a link to a thin_site
        thick_adjacent:
            each tuple location is a tuple of links to adjacent thin filaments
        thick_face:
            each tuple location is a tuple of length six, each location of
            which contains a tuple of links to myosin heads that are facing
            each surrounding thin filament
        thick_bare_zone:
            a single value, the length of each filament before the first crown
        thick_crown_spacing:
            a single value, the distance between two crowns on a single filament
        thick_k:
            a single value, the spring constant of the thick filament between
            any given pair of crowns

        ## Thin Filament Properties: arranged in the same manner as the
        thick filament properties, but for the eight thin filaments

        thin_location:
            each tuple location is a list of x axis locations
        thin_link:
            each tuple location is a list consisting of entries (one for each
            thin_site on the thin_filament) of either a None or a link to a
            thick_crown
        thin_adjacent:
            each tuple location is a tuple of links to adjacent thick filaments
        thin_face:
            each tuple location is a tuple of length three, each location of
            which contains a tuple of links to thin filament sites that are
            facing each surrounding thick filament
        thin_site_spacing:
            the axial distance from one thin filament binding site to another
        thin_k:
            a single value, the spring constant of the thin filament between
            any given pair of thin binding sites

        """
        # Versioning, to be updated when backwards incompatible changes to the
        # data structure are made, not on release of new features
        self.version = 1.4  # Includes support for tropomyosin AND titin
        # Begin handling kwargs
        titin_params = None
        if 'ti_params' in kwargs.keys():
            titin_params = kwargs['ti_params']
        # Parse initial LS and Z-line
        if time_dependence is not None:
            if 'lattice_spacing' in time_dependence:
                lattice_spacing = time_dependence['lattice_spacing'][0]
            if 'z_line' in time_dependence:
                z_line = time_dependence['z_line'][0]
            if 'pCa' in time_dependence:
                pCa = time_dependence['pCa'][0]
        self.time_dependence = time_dependence
        # The next few lines use detection of None rather than a sensible
        # default value as a passed None is an explicit selection of default
        if lattice_spacing is None:
            lattice_spacing = 14.0
        if z_line is None:
            z_line = 1250
        if pCa is None:
            pCa = 4.0
        if poisson is None:
            poisson = 0.0
        # Record initial values for use with poisson driven ls
        self._initial_z_line = z_line
        self._initial_lattice_spacing = lattice_spacing
        self.poisson_ratio = poisson
        # Store these values for posterity
        self.lattice_spacing = lattice_spacing
        self.z_line = z_line
        self.pCa = pCa
        # Create the thin filaments, unlinked but oriented on creation.
        thin_orientations = ([4, 0, 2], [3, 5, 1], [4, 0, 2], [3, 5, 1],
                             [3, 5, 1], [4, 0, 2], [3, 5, 1], [4, 0, 2])
        np.random.seed(None)
        if starts is None:
            thin_starts = []
            for i in range(len(thin_orientations)):
                thin_starts.append(np.random.randint(25))
        else:
            thin_starts = starts[0]
        self._thin_starts = thin_starts
        thin_ids = range(len(thin_orientations))
        new_thin = lambda thin_id: af.ThinFilament(self, thin_id, thin_orientations[thin_id],
                                                   thin_starts[thin_id])
        self.thin = tuple([new_thin(thin_id) for thin_id in thin_ids])
        # Determine the hiding line
        self.hiding_line = None
        self.update_hiding_line()

        '''Begin connecting things'''
        # Create the thick filaments, remembering they are arranged thus:
        # ----------------------------
        # |   Actin around myosin    |
        # |--------------------------|
        # |      a1      a3          |
        # |  a0      a2      a0      |
        # |      M0      M1          |
        # |  a4      a6      a4      |
        # |      a5      a7      a5  |
        # |          M2      M3      |
        # |      a1      a3      a1  |
        # |          a2      a0      |
        # ----------------------------
        # and that when choosing which actin face to link to which thick
        # filament face, use these face orders:
        # ----------------------------------------------------
        # | Myosin face order  |       Actin face order      |
        # |--------------------|-----------------------------|
        # |         a1         |                             |
        # |     a0      a2     |  m0      m1         m0      |
        # |         mf         |      af      OR             |
        # |     a5      a3     |                     af      |
        # |         a4         |      m2         m2      m1  |
        # ----------------------------------------------------
        if starts is None:
            thick_starts = []
            for i in range(4):
                thick_starts.append(np.random.randint(1, 4))
        else:
            thick_starts = starts[1]
        self._thick_starts = thick_starts
        self.thick = (
            mf.ThickFilament(self, 0, (
                self.thin[0].thin_faces[1], self.thin[1].thin_faces[2],
                self.thin[2].thin_faces[2], self.thin[6].thin_faces[0],
                self.thin[5].thin_faces[0], self.thin[4].thin_faces[1]),
                             thick_starts[0]),
            mf.ThickFilament(self, 1, (
                self.thin[2].thin_faces[1], self.thin[3].thin_faces[2],
                self.thin[0].thin_faces[2], self.thin[4].thin_faces[0],
                self.thin[7].thin_faces[0], self.thin[6].thin_faces[1]),
                             thick_starts[1]),
            mf.ThickFilament(self, 2, (
                self.thin[5].thin_faces[1], self.thin[6].thin_faces[2],
                self.thin[7].thin_faces[2], self.thin[3].thin_faces[0],
                self.thin[2].thin_faces[0], self.thin[1].thin_faces[1]),
                             thick_starts[2]),
            mf.ThickFilament(self, 3, (
                self.thin[7].thin_faces[1], self.thin[4].thin_faces[2],
                self.thin[5].thin_faces[2], self.thin[1].thin_faces[0],
                self.thin[0].thin_faces[0], self.thin[3].thin_faces[1]),
                             thick_starts[3])
        )
        # Now the thin filaments need to be linked to thick filaments, use
        # the face orders from above and the following arrangement:
        # ----------------------------
        # |   Myosin around actin    |
        # |--------------------------|
        # |      m3      m2      m3  |
        # |          A1      A3      |
        # |      A0      A2          |
        # |  m1      m0      m1      |
        # |      A4      A6          |
        # |          A5      A7      |
        # |      m3      m2      m3  |
        # ----------------------------
        # The following may be hard to read, but it has been checked and
        # may be moderately trusted. CDW-20100406
        self.thin[0].set_thick_faces((self.thick[3].thick_faces[4],
                                      self.thick[0].thick_faces[0], self.thick[1].thick_faces[2]))
        self.thin[1].set_thick_faces((self.thick[3].thick_faces[3],
                                      self.thick[2].thick_faces[5], self.thick[0].thick_faces[1]))
        self.thin[2].set_thick_faces((self.thick[2].thick_faces[4],
                                      self.thick[1].thick_faces[0], self.thick[0].thick_faces[2]))
        self.thin[3].set_thick_faces((self.thick[2].thick_faces[3],
                                      self.thick[3].thick_faces[5], self.thick[1].thick_faces[1]))
        self.thin[4].set_thick_faces((self.thick[1].thick_faces[3],
                                      self.thick[0].thick_faces[5], self.thick[3].thick_faces[1]))
        self.thin[5].set_thick_faces((self.thick[0].thick_faces[4],
                                      self.thick[2].thick_faces[0], self.thick[3].thick_faces[2]))
        self.thin[6].set_thick_faces((self.thick[0].thick_faces[3],
                                      self.thick[1].thick_faces[5], self.thick[2].thick_faces[1]))
        self.thin[7].set_thick_faces((self.thick[1].thick_faces[4],
                                      self.thick[3].thick_faces[0], self.thick[2].thick_faces[2]))
        # Create the titin filaments and link them from thick
        # faces to thin faces
        # |--------------------------------------------------|
        # |            Actin & titin around myosin           |
        # |--------------------------------------------------|
        # |           a1               a3                    |
        # |                                                  |
        # |  a0       t1      a2       t4       a0           |
        # |       t0     t2        t3      t5                |
        # |           M0               M1                    |
        # |       t6     t8        t9      t11               |
        # |  a4       t7      a6       t10      a4           |
        # |                                                  |
        # |           a5     t13       a7       t16      a5  |
        # |               t12    t14        t15    t17       |
        # |                   M2                M3           |
        # |               t18    t20        t21    t23   a1  |
        # |           a1      t19      a3       t22          |
        # |                                                  |
        # |                   a2                a0           |
        # |--------------------------------------------------|
        # ## CHECK_JDP ## Link Thick filament to titin
        # ## Checked - AMA 13JUN19
        """ titin connection loop format:
        link_list = []
        num = 0
        an_list = [0, 1, 2, 2, 3, 0]
        for half in range(0, 2):
            mf_list = [0, 1, 2]
            af_list = [1, 2, 2]

            for quart in range(0, 2):
                for eighth in range(0, 2):
                    for triple in range(0, 3):
                        m_n = eighth + half * 2
                        m_f = mf_list[triple]
                        a_n = an_list[triple + eighth * 3]
                        a_f = af_list[triple]
                        link_list.append((num, m_n, m_f, a_n, a_f))
                        num += 1
                mf_list = [5, 4, 3]
                af_list = [1, 0, 0]
                for i in range(0, len(an_list)):
                    an_list[i] = an_list[i] - half * 8 + 4
            an_list = [5, 6, 7, 7, 4, 5]

        for item in link_list:
            print("ti.Titin(self, " + str(item[0]) + ", ti_thick(" +
                  str(item[1]) + ", " + str(item[2]) + "), ti_thin(" +
                  str(item[3]) + ", " + str(item[4]) + "), a=a, b=b),")
        # """

        a = None
        b = None
        if titin_params is not None:
            a = titin_params[HS_TI_A]
            b = titin_params[HS_TI_B]
        ti_thick = lambda thick_i, j: self.thick[thick_i].thick_faces[j]
        ti_thin = lambda thin_i, j: self.thin[thin_i].thin_faces[j]
        self.titin = (
            ti.Titin(self, 0, ti_thick(0, 0), ti_thin(0, 1), a=a, b=b),
            ti.Titin(self, 1, ti_thick(0, 1), ti_thin(1, 2), a=a, b=b),
            ti.Titin(self, 2, ti_thick(0, 2), ti_thin(2, 2), a=a, b=b),
            ti.Titin(self, 3, ti_thick(1, 0), ti_thin(2, 1), a=a, b=b),
            ti.Titin(self, 4, ti_thick(1, 1), ti_thin(3, 2), a=a, b=b),
            ti.Titin(self, 5, ti_thick(1, 2), ti_thin(0, 2), a=a, b=b),

            ti.Titin(self, 6, ti_thick(0, 5), ti_thin(4, 1), a=a, b=b),
            ti.Titin(self, 7, ti_thick(0, 4), ti_thin(5, 0), a=a, b=b),
            ti.Titin(self, 8, ti_thick(0, 3), ti_thin(6, 0), a=a, b=b),
            ti.Titin(self, 9, ti_thick(1, 5), ti_thin(6, 1), a=a, b=b),
            ti.Titin(self, 10, ti_thick(1, 4), ti_thin(7, 0), a=a, b=b),
            ti.Titin(self, 11, ti_thick(1, 3), ti_thin(4, 0), a=a, b=b),

            ti.Titin(self, 12, ti_thick(2, 0), ti_thin(5, 1), a=a, b=b),
            ti.Titin(self, 13, ti_thick(2, 1), ti_thin(6, 2), a=a, b=b),
            ti.Titin(self, 14, ti_thick(2, 2), ti_thin(7, 2), a=a, b=b),
            ti.Titin(self, 15, ti_thick(3, 0), ti_thin(7, 1), a=a, b=b),
            ti.Titin(self, 16, ti_thick(3, 1), ti_thin(4, 2), a=a, b=b),
            ti.Titin(self, 17, ti_thick(3, 2), ti_thin(5, 2), a=a, b=b),

            ti.Titin(self, 18, ti_thick(2, 5), ti_thin(1, 1), a=a, b=b),
            ti.Titin(self, 19, ti_thick(2, 4), ti_thin(2, 0), a=a, b=b),
            ti.Titin(self, 20, ti_thick(2, 3), ti_thin(3, 0), a=a, b=b),
            ti.Titin(self, 21, ti_thick(3, 5), ti_thin(3, 1), a=a, b=b),
            ti.Titin(self, 22, ti_thick(3, 4), ti_thin(0, 0), a=a, b=b),
            ti.Titin(self, 23, ti_thick(3, 3), ti_thin(1, 0), a=a, b=b),
        )
        '''Initialize the last few variables'''
        # Set the timestep for all our new cross-bridges
        self.timestep_len = timestep_len
        # Track how long we've been running
        self.current_timestep = 0

        # ## Tropomyosin species concentrations
        self.c_tn = None
        self.c_ca = None,
        self.c_tnca = None
        # TODO determine if we need to self.update_concentrations()

        # ## variables previously initialized in methods (hiding line included above)
        self.last_transitions = None
        self._volume = None
        self.update_volume()

    def to_dict(self):
        """Create a JSON compatible representation of the thick filament

        Example usage: json.dumps(sarc.to_dict(), indent=1)

        Current output includes:
            version: version of the sarcomere model
            timestep_len: the length of the timestep in ms
            current_timestep: time to get a watch
            lattice_spacing: the thick to thin distance
            z_line: the z_line location
            pCa: the calcium level
            hiding_line: where binding sites become unavailable due to overlap
            time_dependence: dictionary of how "lattice_spacing", "z_line", 
                and "pCa" can change
            last_transitions: keeps track of the last state change by thick
                filament and by crown
            thick: the structures for the thick filaments
            thin: the structures for the thin filaments
        """
        sd = self.__dict__.copy()  # sarc dict
        sd['current_timestep'] = self.current_timestep
        # set act_perm as mean since prop access returns values at every point
        sd['actin_permissiveness'] = np.mean(self.actin_permissiveness)
        sd['thick'] = [t.to_dict() for t in sd['thick']]
        sd['thin'] = [t.to_dict() for t in sd['thin']]
        return sd

    def from_dict(self, sd):
        """ Load values from a sarcomere dict. Values read in correspond to
        the current output documented in to_dict.
        """
        # Warn of possible version mismatches
        read, current = sd['version'], self.version
        if read != current:
            import warnings
            warnings.warn("Versioning mismatch, reading %0.1f into %0.1f."
                          % (read, current))
        # Get filaments in right orientations
        self.__init__(
            lattice_spacing=sd['_initial_lattice_spacing'],
            z_line=sd['_initial_z_line'],
            poisson=sd['poisson_ratio'],
            pCa=sd['pCa'],
            timestep_len=sd['timestep_len'],
            time_dependence=sd['time_dependence'],
            starts=(sd['_thin_starts'], sd['_thick_starts'])
        )
        # Local keys
        self.current_timestep = sd['current_timestep']
        self._z_line = sd['_z_line']
        self._lattice_spacing = sd['_lattice_spacing']
        self.hiding_line = sd['hiding_line']
        if 'last_transitions' in sd.keys():
            self.last_transitions = sd['last_transitions']
        # Sub-structure keys
        for data, thick in zip(sd['thick'], self.thick):
            thick.from_dict(data)
        for data, thin in zip(sd['thin'], self.thin):
            thin.from_dict(data)

    def run(self, time_steps=100, callback=None, bar=True, every=1, print_every=10):
        """Run the model for the specified number of timesteps

        Parameters:
            time_steps: number of time steps to run the model for (100)
            callback: function to be executed after each time step to
                collect data. The callback function takes the sarcomere
                in its current state as its only argument. (Defaults to
                the axial force at the M-line if not specified.)
            bar: progress bar control,False means don't display, True
                means give us the basic progress reports, if a function
                is passed, it will be called as f(completed_steps,
                total_steps, sec_left, sec_passed, process_name).
                (Defaults to True)
            every: how many timesteps to update after
            print_every: how many timesteps to print update after
        Returns:
            output: the results of the callback after each timestep
        """
        # Callback defaults to the axial force at the M-line
        if callback is None:
            callback = self.tm_report

        # ## logic to handle bar is type(True || False || Function)
        use_bar = False
        update_bar = self.print_bar
        if isinstance(bar, bool):
            use_bar = bar
        elif isinstance(bar, type(lambda x: x)):
            use_bar = True
            update_bar = bar
        # Create a place to store callback information and note the time
        output = {}
        tic = time.time()
        # Run through each timestep
        for i in range(time_steps):
            try:
                self.timestep(print_state=i % print_every == 0)
                output = self._append(callback(), output)
                # Update us on how it went
                toc = int((time.time() - tic) / (i + 1) * (time_steps - i - 1))
                proc_name = mp.current_process().name

                if use_bar and i % every == 0:
                    update_bar(i, time_steps, toc, time.time() - tic, proc_name)
            except KeyboardInterrupt:
                return output, 130
            except Exception as e:
                import traceback
                print("/n")
                print(e)
                traceback.print_exc()
                return output, 1
        return output, 0

    @staticmethod
    def print_bar(i, time_steps, toc, tic, proc_name):
        if tic < -1:
            print('Causality has failed')
        sys.stdout.write("\n" + proc_name +
                         " finished timestep %i of %i, %ih%im%is left"
                         % (i + 1, time_steps, toc / 60 / 60, toc / 60 % 60, toc % 60))
        sys.stdout.flush()

    @staticmethod
    def _append(data_dict, dictionary):
        for key, value in data_dict.items():
            if key in dictionary.keys():
                dictionary[key].append(value)
            else:
                dictionary.update({key: [value]})
        return dictionary

    def timestep(self, current=None, print_state=True):
        """Move the model one step forward in time, allowing the
        myosin heads a chance to bind and then balancing forces
        """
        # Record our passage through time
        if current is not None:
            self.current_timestep = current
        else:
            self.current_timestep += 1
        # Update bound states
        if print_state:
            print(" F", end=" ")
        self.last_transitions = [thick.transition() for thick in self.thick]
        if print_state:
            print(" T", end=" ")
        self.update_concentrations()
        tm_activation = [thin.transition() for thin in self.thin]  # TODO: work into storage
        active = 0
        total = 0
        for i in range(0, len(tm_activation)):
            active += tm_activation[i][0]
            total += tm_activation[i][1]
        if print_state:
            print(active / total)
        # Settle forces
        self.settle()

    @property
    def current_timestep(self):
        """Return the current timestep"""
        return self._current_timestep

    @current_timestep.setter
    def current_timestep(self, new_timestep):
        """Set the current timestep"""
        # Update boundary conditions
        self.update_hiding_line()
        td = self.time_dependence
        i = new_timestep
        if td is not None:
            if 'lattice_spacing' in td:
                self.lattice_spacing = td['lattice_spacing'][i]
            if 'z_line' in td:
                self.z_line = td['z_line'][i]
            if 'pCa' in td:
                self.pCa = td['pCa'][i]
        self._current_timestep = i
        return

    @property
    def actin_permissiveness(self):
        """How active & open to binding, 0 to 1, are binding sites?"""
        return [thin.permissiveness for thin in self.thin]

    @property
    def z_line(self):
        """Axial location of the z-line, length of the half sarcomere"""
        return self._z_line

    @z_line.setter
    def z_line(self, new_z_line):
        """Set a new z-line, updating the lattice spacing at the same time"""
        self._z_line = new_z_line
        self.update_ls_from_poisson_ratio()
        self.update_volume()

    @property
    def lattice_spacing(self):
        """Return the current lattice spacing"""
        return self._lattice_spacing

    @lattice_spacing.setter
    def lattice_spacing(self, new_lattice_spacing):
        """Assign a new lattice spacing"""
        self._lattice_spacing = new_lattice_spacing

    @property
    def _tn_count(self):
        return self._tn_total - self._tnca_count

    @property
    def _tnca_count(self):
        bound = 0
        for thin in self.thin:
            for tm in thin.tm:
                for tm_site in tm.sites:
                    if tm_site.state != 0:
                        bound += 1
        return bound

    @property
    def _tn_total(self):
        total_tn = 0
        for thin in self.thin:
            for tm in thin.tm:
                total_tn += len(tm.sites)
        return total_tn

    def update_concentrations(self):
        self.c_tn = self._concentration(self._tn_count)
        self.c_ca = 10.0 ** (-self.pCa)
        self.c_tnca = self._concentration(self._tnca_count)

    @property
    def concentrations(self):
        return {"free_tm": self.c_tn,
                "free_ca": self.c_ca,
                "bound_tm": self.c_tnca}

    def tm_report(self):
        report = {}
        # calculate average rates across all binding sites and average cooperativity state
        rates = [0, 0, 0, 0, 0, 0]
        count = 0
        coop = 0
        for thin in self.thin:
            for tm in thin.tm:
                for tm_site in tm.sites:
                    if tm_site.subject_to_cooperativity:
                        coop += 1
                    new_rates = tm_site.get_rates
                    for i in range(len(rates)):
                        rates[i] += new_rates[i]
                    count += 1
        for i in range(len(rates)):
            rates[i] /= count

        rate_keys = ["r_12", "r_21", "r_23", "r_32", "r_31", "r_13"]
        for i in range(len(rates)):
            report.update({rate_keys[i]: rates[i]})

        xb_fracs = self.get_frac_in_states()

        report.update({"coop": coop / count,
                       "axial_force": self.axial_force(),
                       "ca": self.c_ca,
                       "xb_fraction_free": xb_fracs[0],
                       "xb_fraction_loose": xb_fracs[1],
                       "xb_fraction_tight": xb_fracs[2]})
        report.update(self.concentrations)
        return report

    @property
    def volume(self):
        """return the current fluid volume of the half-sarcomere AMA-11JAN2020"""
        return self._volume

    # @volume.setter
    def update_volume(self):
        """re-calculate the fluid volume of the half sarcomere - ASSUMING CONSTANT LATTICE SPACING
        returns volume in L(AMA-3FEB2020)
        (nm3 AMA-11JAN2020"""
        ls = self._lattice_spacing
        length = self._z_line

        # calculate area of 4 hexagons, with edge length 9/2 + 16/2 + ls
        edge = 9 / 2 + 16 / 2 + ls
        area = 4 * 3 / 2 * np.sqrt(3) * edge * edge

        thin_volume = np.pi * 22659.75  # radius 4.5 * radius 4.5 * length 1119
        thick_volume = np.pi * 58624  # radius 8 * radius 8 * length 916
        filament_volume = 10 * thin_volume + 4 * thick_volume
        whole_volume = area * length
        fluid_volume = whole_volume - filament_volume
        nm3_p_L = 1e-24
        fluid_volume *= nm3_p_L
        self._volume = fluid_volume

    def _concentration(self, count):
        return count / self.volume

    @staticmethod
    def ls_to_d10(face_dist):
        """Convert face-to-face lattice spacing to d10 spacing.

        Governing equations:
            ls = ftf, the face to face distance
            filcenter_dist = face_dist + .5 * dia_actin + .5 * dia_myosin
            d10 = 1.5 * filcenter_dist
        Values:
            dia_actin: 9nm [1]_
            dia_myosin: 16nm [2]_
            example d10: 37nm for cardiac muscle at 2.2um [3]_
        References:
            .. [1] Egelman 1985, The structure of F-actin.
                   J Muscle Res Cell Motil, Pg 130, values from 9 to 10 nm
            .. [2] Woodhead et al. 2005, Atomic model of a myosin filament in
                   the relaxed state. Nature, Pg 1195, in tarantula filament
            .. [3] Millman 1998, The filament lattice of striated muscle.
                   Physiol Rev,  Pg 375
        Note: Arguably this should be moved to a support class as it really
        isn't something the half-sarcomere knows about or does. I'm leaving it
        here as a convenience for now.

        Parameters:
            face_dist: face to face lattice spacing in nm
        Returns:
            d10: d10 spacing in nm
        """
        filcenter_dist = face_dist + 0.5 * 9 + 0.5 * 16
        d10 = 1.5 * filcenter_dist
        return d10

    @staticmethod
    def d10_to_ls(d10):
        """Convert d10 spacing to face-to-face lattice spacing

        Governing equations: See ls_to_d10
        Values: See ls_to_d10

        Parameters:
            d10: d10 spacing in nm
        Returns:
            face_dist: face to face lattice spacing in nm
        """
        filcenter_dist = d10 * 2 / 3
        face_dist = filcenter_dist - 0.5 * 9 - 0.5 * 16
        return face_dist

    def axial_force(self):
        """Sum of each thick filament's axial force on the M-line """
        return sum([thick.effective_axial_force() for thick in self.thick])

    def radial_tension(self):
        """The sum of the thick filaments' radial tensions"""
        return sum([t.radial_tension() for t in self.thick])

    def radial_force(self):
        """The sum of the thick filaments' radial forces, as a (y,z) vector"""
        return np.sum([t.radial_force_of_filament() for t in self.thick], 0)

    def _single_settle(self, factor=0.95):
        """Settle down now, just a little bit"""
        thick = [thick.settle(factor) for thick in self.thick]
        thin = [thin.settle(factor) for thin in self.thin]
        return np.max((np.max(np.abs(thick)), np.max(np.abs(thin))))

    def settle(self):
        """Jiggle those locations around until the residual forces are low

        We choose the convergence limit so that 95% of thermal forcing events
        result in a deformation that produces more axial force than the
        convergence value, 0.12pN.
        """
        converge_limit = 0.12  # see doc string
        converge = self._single_settle()
        while converge > converge_limit:
            converge = self._single_settle()

    def _get_residual(self):
        """Get the residual force at every point in the half-sarcomere"""
        thick_f = np.hstack([t.axial_force() for t in self.thick])
        thin_f = np.hstack([t.axial_force() for t in self.thin])
        mash = np.hstack([thick_f, thin_f])
        return mash

    def get_frac_in_states(self):
        """Calculate the fraction of cross-bridges in each state"""
        nested = [t.get_states() for t in self.thick]
        xb_states = [xb for fil in nested for face in fil for xb in face]
        num_in_state = [xb_states.count(state) for state in range(3)]
        frac_in_state = [n / float(len(xb_states)) for n in num_in_state]
        return frac_in_state

    def update_ls_from_poisson_ratio(self):
        """Update the lattice spacing consistent with the poisson ratio,
        initial lattice spacing, current z-line, and initial z-line

        Governing equations
        ===================
        Poisson ratio := ν
            ν = dε_r/dε_z = Δr/r_0 / Δz/z_0
        From Mathematica derivation
        γ := center to center distance between filaments
            γ(ν, γ_0, z_0, Δz) = γ_0 (z_0/(z_0+Δz))^ν
        And since we want the face-to-face distance, aka ls, we convert with:
            γ = ls + 0.5 (dia_actin + dia_myosin)
        and
            γ_0 = ls_0 + 0.5 (dia_actin + dia_myosin)
        and the simplifying
            β = 0.5 (dia_actin + dia_myosin)
        to get
            ls = (ls_0 + β) (z_0/(z_0 + Δz))^ν - β
        which is what we implement below.
        Note: this is a novel derivation and so there is no current
            citation to be invoked.

        Values: See ls_to_d10

        Parameters:
            self: half-sarcomere, automatically passed
        Returns:
            None
        """
        beta = 0.5 * (9 + 16)
        ls_0 = self._initial_lattice_spacing
        z_0 = self._initial_z_line
        nu = self.poisson_ratio
        dz = self.z_line - z_0
        ls = (ls_0 + beta) * (z_0 / (z_0 + dz)) ** nu - beta
        self.lattice_spacing = ls
        return

    def update_hiding_line(self):
        """Update the line determining which actin sites are unavailable"""
        farthest_actin = min([min(thin.axial) for thin in self.thin])
        self.hiding_line = -farthest_actin

    def resolve_address(self, address):
        """Give back a link to the object specified in the address
        Addresses are formatted as the object type (string) followed by a list
        of the indices that the object occupies in each level of organization.
        Valid string values are:
            thin_fil
            thin_face
            bs
            thick_fil
            crown
            thick_face
            xb
        and an example valid address would be ('bs', 1, 14) for the binding
        site at index 14 on the thin filament at index 1.
        """
        if address[0] == 'thin_fil':
            return self.thin[address[1]]
        elif address[0] in ['thin_face', 'bs', 'tm', 'tm_site']:
            return self.thin[address[1]].resolve_address(address)
        elif address[0] == 'thick_fil':
            return self.thick[address[1]]
        elif address[0] in ['crown', 'thick_face', 'xb']:
            return self.thick[address[1]].resolve_address(address)
        import warnings
        warnings.warn("Unresolvable address: %s" % str(address))


sarc = hs()
