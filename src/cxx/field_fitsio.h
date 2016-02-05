/*
 * field_fitsio.h - This file is part of MRS3D
 * 
 * Copyright 2011 CEA
 * Contributor : François Lanusse (francois.lanusse@gmail.com)
 * 
 * Created on 19/06/11
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

#ifndef FIELD_FITSIO_H
#define FIELD_FITSIO_H

#include <string>
#include <datatypes.h>
#include <fitshandle.h>
class fitshandle;

template<typename T> class Field;

/*! \defgroup field_fitsio_group FITS-based I/O of 3D density field */
/*! \{ */

/*! Returns the number \a Nside of pixels along one edge of the cube and \a size the physical size of the cube in Mpc/h */
template<typename T> bool get_fieldsize(fitshandle &inp, int &Nside, T &size);
/*! Returns the number of pixels \a Nside and \a size of one edge of the cube found in the HDU
 \a hdunum of file \a filename in \a Nside and \a size. */
template<typename T> bool get_fieldsize(const std::string &filename, int &Nside, T &size, int hdunum=2);

/*! Reads the Field of the FITS binary table pointed to by \a inp into
 \a field. \a field is reallocated with the parameters \a Nside and \a size.
 Values not present in the FITS table are set to zero. */
template<typename T> void read_Field_from_fits
(fitshandle &inp, Field<T> &field);
/*! Opens the FITS file \a filename, jumps to the HDU \a hdunum, then reads
 the Field from the FITS binary table there into \a field. \a field is
 reallocated with the parameters \a Nside and \a size.
 Values not present in the FITS table are set to zero. */
template<typename T> void read_Field_from_fits
(const std::string &filename, Field<T> &field, int hdunum=2);

/*! Inserts a new binary table into \a out, which contains three columns
 of type PLANCK_INT32, \a datatype and \a datatype, respectively.
 The data in \a field is written into this table. */
template<typename T> void write_Field_to_fits
(fitshandle &out, const Field<T> &field, PDT datatype);

template<typename T> inline void write_Field_to_fits
(const std::string &outfile, const Field<T> &field, PDT datatype)
{
	fitshandle out;
	out.create(outfile);
	write_Field_to_fits (out, field,datatype);
}

/*! \} */

#endif
