/*
 * sbtools.h - This file is part of FastDSBT
 * 
 * Copyright 2011 CEA
 * Contributor : François Lanusse (francois.lanusse@gmail.com)
 * 
 * Created on 04/05/11
 * 
 * This software is a computer program whose purpose is to compute the
 * Spherical Bessel Transform.
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

#ifndef ST_TOOLS_H
#define ST_TOOLS_H

#include <stdint.h>

extern "C" {
	void __sbtools_MOD_brootfinder(double* qln, int64_t *nnmax, int64_t *nlmax);
	void __sbtools_MOD_bjl(int64_t* l, double* x, double* jl);
}

/*! \defgroup dsbt_sbtools_group Spherical Bessel utilities used by FastDSBT */
/*! \{ */

/*! Computes the first \a nmax roots of the bessel functions of order up to \a lmax.
 *  The results are stored in a pre-allocated array \a qln of size nmax*(lmax +1). */
void BesselRoots(int64_t nmax, int64_t lmax, double* qln);

/*! Returns the spherical Bessel function of order \a l evaluated on \a x */ 
double sphericalBesselJ(int64_t l, double x);

/*! Save the zeros of bessel functions from the array \a qln to a fits file named
 *  \a tableFileName. \a nnmax and \a nlmax specify the size of the \a qln array */
void write_qln(char* tableFileName, int64_t nnmax, int64_t nlmax, double *qln);

/*! Load the zeros of bessel functions from the fits file \a tableFileName into 
 * the array \a qln. \a nnmax and \a nlmax specify the size of the \a qln array */
int load_qln(const char* tableFileName, int64_t nnmax, int64_t nlmax, double *qln);

/*! \} */
#endif
