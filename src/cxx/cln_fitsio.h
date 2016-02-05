/*
 * cln_fitsio.h - This file is part of MRS3D
 * 
 * Copyright 2011 CEA
 * Contributor : François Lanusse (francois.lanusse@gmail.com)
 * 
 * Created on 16/05/11
 * 
 * This software is a computer program whose purpose is to apply mutli-
 * resolution signal processing algorithms on spherical 3D data.
 * 
 * This software is governed by the CeCILL  license under French law and
 * abiding by the rules of distribution of free software.  You can  use, 
 * modify and/ or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info". 
 * 
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability. 
 *
 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or 
 * data to be ensured and,  more generally, to use and operate it in the 
 * same conditions as regards security. 
 * 
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 *
 */

/*
 * This file is inspired from the Healpix_cxx software package distributed
 * under GNU General Public License.
 * 
 * For more information about HEALPix, see http://healpix.jpl.nasa.gov
 * 
 */

#ifndef CLN_FITSIO_H
#define CLN_FITSIO_H

#include <string>
#include <datatypes.h>
#include <fitshandle.h>
class fitshandle;

template<typename T> class Cln;

/*! \defgroup cln_fitsio_group FITS-based I/O of power spectra */
/*! \{ */

/*! Returns the maximum \a l and \a n moments found in the FITS HDU
 pointed to be \a inp in \a lmax and \a nmax. */
bool get_clnsize(fitshandle &inp, int &lmax, int &nmax);
/*! Returns the maximum \a l and \a n moments found in the HDU
 \a hdunum of file \a filename in \a lmax and \a nmax. */
bool get_clnsize(const std::string &filename, int &lmax, int &nmax, int hdunum=2);

/*! Reads the c_ln of the FITS binary table pointed to by \a inp into
 \a cln. \a cln is reallocated with the parameters \a lmax and \a nmax.
 Values not present in the FITS table are set to zero; values outside
 the requested (l,n) range are ignored. */
template<typename T> void read_Cln_from_fits(fitshandle &inp, Cln<T> &cln, int lmax, int nmax);
/*! Opens the FITS file \a filename, jumps to the HDU \a hdunum, then reads
 the c_ln from the FITS binary table there into \a cln. \a cln is
 reallocated with the parameters \a lmax and \a nmax.
 Values not present in the FITS table are set to zero; values outside
 the requested (\a l,\a n) range are ignored. */
template<typename T> void read_Cln_from_fits(const std::string &filename, Cln<T> &cln, int lmax, int nmax, int hdunum=2);

/*! Inserts a new binary table into \a out, which contains two columns
 of type PLANCK_INT32, \a datatype and \a datatype, respectively.
 The data in \a cln is written into this table; values outside
 the requested (\a lmax , \a nmax) range are omitted. */
template<typename T> void write_Cln_to_fits(fitshandle &out, const Cln<T> &cln, int lmax, int nmax, PDT datatype);

template<typename T> inline void write_Cln_to_fits(const std::string &outfile, const Cln<T> &cln, int lmax, int nmax, PDT datatype)
{
	fitshandle out;
	out.create(outfile);
	write_Cln_to_fits (out, cln, lmax, nmax, datatype);
}

/*! \} */

#endif
