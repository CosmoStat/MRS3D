/*
 * almn_cln_tools.h - This file is part of MRS3D
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

#ifndef ALMN_POWSPEC_TOOLS_H
#define ALMN_POWSPEC_TOOLS_H

#include "xcomplex.h"

template<typename T> class Almn;
template<typename T> class Cln;
class planck_rng;

/*! \defgroup almn_ps_group Conversions between a_lmn and power spectra */
/*! \{ */

/*! Creates a Gaussian realisation of the input power spectrum \a cln,
 using the random number generator \a rng, and returns the result
 in \a almn. */
template<typename T1, typename T2> void create_almn (const Cln<T1> &powspec, Almn<xcomplex<T2> > &almn, planck_rng &rng);

/*! Returns the power spectrum of \a almn in \a cln. */
template<typename T1, typename T2> void extract_cln
(const Almn<xcomplex<T1> > &almn, Cln<T2> &cln);
/*! Returns the cross power spectrum of \a almn1 and \a almn2 in \a cln. */
template<typename T1, typename T2> void extract_crosscln
(const Almn<xcomplex<T1> > &alm1, const Almn<xcomplex<T1> > &alm2, Cln<T2> &cln);

/*! \} */

#endif
