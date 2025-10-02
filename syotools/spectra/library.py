#!/usr/bin/env python
"""
Created on Tue Oct 18 11:19:05 2016
@author: gkanarek
"""
import os
import astropy.units as u
import synphot as syn
import specutils as specu
from syotools.spectra.spec_defaults import default_spectra, syn_spectra_library 
from pathlib import Path

class _spec_library(object):
    """
    This is a container object for spectra, which can handle spectra stored
    in the data/ folder, as well as accept new spectra from user upload.
    
    synphot is the package which will handle all the gruntwork of spectrum
    processing; we also use specutils for somewhat more forgiving file IO.
    
    Each spectrum has an id (given by the user for non-default spectra), which
    is used to access the spectrum and its optional description.

    the property _available_spectra is a dictionary with keys for each spectrum 
    
    PLEASE NOTE:
        Accessing a spectrum with object dot notation (using its spec id) will
    return the spectrum object itself; HOWEVER, accessing via indexing will
    return the spectrum's description string (if any) instead.
    """
    
    _available_spectra = {} 
    _descriptions = {}
    
    def __init__(self):
        """
        Load the default spectra.
        """
        self._available_spectra.update(default_spectra['specs'])
        self._descriptions.update(default_spectra['descs'])
    
    #Implement dot and dict index access to the items stored in the container:
    def __getattr__(self, key):
        return self._available_spectra.get(key)

    def __setattr__(self, key, value):
        """
        Add or change an item in the _available_spectra dict. We do some type
        checking & parsing to keep conversions under the hood; everything
        should eventually be stored as an SourceSpectrum or FileSpectrum.
        """
        if isinstance(value, syn.spectrum.SourceSpectrum):
            self._available_spectra[key] = value
        elif isinstance(value, specu.Spectrum1D):
            self.add_spec_from_spectrum1d(key, value)
        elif all(isinstance(value, (list,tuple)), len(value) == 2):
            self.add_spec_from_arrays(key, *value)
        else:
            name = lambda cls: '{}.{}'.format(cls.__module__, cls.__name__)
            allowed = ', '.join([name(s) for s in (syn.spectrum.SourceSpectrum,
                                                 specu.Spectrum1D, u.Quantity, 
                                                 list, tuple)])
            raise TypeError('Only the following types are supported: '+allowed)
        self._descriptions[key] = ''

    def __delattr__(self, key):
        del self._available_spectra[key]
        if key in self._descriptions:
            del self._descriptions[key]

    def __getitem__(self, key):
        return self._descriptions.get(key, '')
    
    def __setitem__(self, key, value):
        if key not in self._available_spectra:
            raise IndexError("Spectrum with id '{}' does not exist".format(key))
        if not isinstance(key, str):
            raise TypeError("Only string descriptions allowed.")
        self._descriptions[key] = value

    def __delitem__(self, key):
        if key in self._descriptions:
            self._descriptions[key] = ''

    #More dict-style access
    def update(self, new_dict):
        for key, value in new_dict.items():
            try:
                self.key = value
            except TypeError as err:
                msg = 'For spectrum id {}, o'.format(key) + err[1:]
                raise TypeError(msg)
    
    def keys(self):
        for k in sorted(self._available_spectra):
            yield k
    
    def values(self):
        for k in sorted(self._available_spectra):
            yield self._available_spectra[k]

    def items(self):
        for k in sorted(self._available_spectra):
            yield (k, self._available_specta[k])
    
    def descriptions(self):
        for k in sorted(self._available_spectra):
            yield self._descriptions[k]
    
    def get(self, key, default=None):
        return self._available_spectra.get(key, default)

    #I/O methods
    def load_spec_from_file(self, filepath, specid, waveunits=None, 
                            fluxunits=None, multispec=False):
        """
        Load a spectrum from a FITS or ascii file.
        
        FITS loading is accomplished by specutils.io.readfits, while ascii
        (and final spectrum format) is handled by synphot.
        
        Arguments:
            filepath  - the file path to be loaded. pathlib.Path.resolve() is
                        is used to convert to an absolute path. (string)
            specname  - the id key for the spectrum in the library, either a 
                        string or (if a multispec file) a list of strings. all
                        strings must be valid for dot access -- i.e., 
                        alphanumeric characters & underscores only.
            waveunits - units for the wavelength array. if None, the units are
                        read from the header, if any. Must be included for an 
                        ascii file, ignored for fits.
                        (astropy.units.Unit, None)
            fluxunits - units for the flux array. if None, the units are read
                        from the header, if any. Must be included for an ascii
                        file, ignored for fits. (astropy.units.Unit, 
                        None)
            multispec - is the file a multispec file? if so, specid MUST be
                        a list of the same length as the number of spectra.
                        (bool, False)
        """
        
        if multispec:
            if not isinstance(specid, (list, tuple)):
                raise TypeError('specid must be a list or tuple when multispec=True')
        
        use_pathlib = True
        if use_pathlib:
            path = Path(filepath).resolve() # using pathlib
            extensions = path.suffixes
            abspath = str(path)
        else:
            path = os.path.expanduser(filepath)
            tmp, extensions = os.path.splitext(path)
            abspath = os.path.abspath(path)
        
        #We load fits files using specutils, because it's more forgiving.
        
        if extensions[-1] == '.fits' or (len(extensions) >= 2 and
                                     extensions[-2:] == ['.fits','.gz']):

            #JT 
            new_spec = syn.spectrum.SourceSpectrum.from_file(abspath, ext=1, keep_neg=True)

            #handle a multispec file
            if multispec:
                if len(new_spec) != len(specid):
                    raise ValueError('Length of specid must equal number of spectra')
                for sid, spec in zip(specid, new_spec):
                    self.add_spec_from_spectrum1d(sid, spec)
            else:
                self.add_spec_from_spectrum1d(specid, new_spec)
        else:
            if waveunits is None or fluxunits is None:
                raise ValueError('Wavelength and flux units must be specified for ascii files')
            sp = syn.spectrum.SourceSpectrum.from_file(abspath, wave_unit=waveunits, flux_unit=fluxunits, keep_neg=True)
            self.add_spec_from_spectrum1d(specid, sp)
    
    def add_spec_from_arrays(self, specid, wavelength, flux):
        """
        Accept Quantity arrays for wavelength and flux, store them in a synphot
        SourceSpectrum, and add them to the library.
        """
        
        if not (isinstance(wavelength, u.Quantity) and isinstance(flux, u.Quantity)):
            raise TypeError('wavelength and flux must be Quantity arrays')
        sp = syn.spectrum.SourceSpectrum(syn.models.Empirical1D, points=wavelength, lookup_table=flux, keep_neg=True)
        self._available_spectra[specid] = sp
    
    def add_spec_from_spectrum1d(self, specid, spectrum):
        """
        Accept specutils.Spectrum1D object, store it in a synphot 
        SourceSpectrum, and add it to the library.
        """
        self.add_spec_from_arrays(specid, spectrum.wavelength, spectrum.flux)

    def save_spec_to_file(self, fname, specid, **options):
        """
        Save a library spectrum to a FITS file, using synphot's to_fits.
        Any options to write_fits_spec may be passed via **options; see
        https://synphot.readthedocs.io/en/latest/api/synphot.spectrum.SourceSpectrum.html#synphot.spectrum.SourceSpectrum.to_fits
        """
        if use_pathlib:
            path = Path(fname)
            abspath = str(path.resolve())
        else:
            abspath = os.path.abspath(os.path.expanduser(fname))
        
        outspec = self._available_spectra[specid]
        outspec.to_fits(abspath, **options)

#Initialize the default library:
SpectralLibrary = _spec_library()
