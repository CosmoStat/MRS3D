/*
 * almn_healpix_tools.h - This file is part of MRS3D
 *
 * Copyright 2011 CEA
 * Contributor : François Lanusse (francois.lanusse@gmail.com)
 *
 * Created on 22/04/11
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

#ifndef ALMN_HEALPIX_TOOLS_H
#define ALMN_HEALPIX_TOOLS_H

#include "xcomplex.h"
#include "arr.h"
#include "fastDSBT.h"
#include "sbtools.h"
#include <alm.h>
#include <alm_healpix_tools.h>
#include <healpix_base.h>
#include <algorithm>
#include <math.h>
#include "almn.h"

/*! \defgroup almn_healpix_group Conversions between a_lmn and Healpix maps */
/*! \{ */

/*!
 * Computes a Healpix \a map from the specified \a almn at a radial distance \a rho (in Mpc/h).
 * \a dsbt is a previously initialized FastDSBT object to be reused by almn2rmap.
 *
 */
template<typename T> void almn2rmap(FastDSBT &dsbt, const Almn<xcomplex<T> > &almn,  Healpix_Map<T> &map, T rho) {

    double Rmax = almn.Rmax();
    double constFactor = sqrt(2.0 * M_PI)/(Rmax*Rmax*Rmax);
    double Kmax ;
    Alm<xcomplex<T> > alm(almn.Lmax(),almn.Mmax());
    alm.SetToZero();
    double coeff = sqrt(2.0*M_PI)/(Rmax*Rmax*Rmax);

    Kmax = dsbt.get_Kln(0, almn.Nmax()-1, Rmax);

    // Inverse Spherical Bessel Transfom
#pragma omp parallel
    for (int n=0; n < almn.Nmax(); ++n) {
#pragma omp for
        for (int l=0; l <= almn.Lmax(); ++l) {
            double normfactor = sphericalBesselJ(l+1, dsbt.get_Kln(l,n,1.0));
            normfactor = coeff/(normfactor*normfactor);
            if (dsbt.get_Kln(l, n, Rmax) <= Kmax) {
                double factor = sphericalBesselJ(l,dsbt.get_Kln(l,n,Rmax)*rho)*normfactor;
                for (int m=0 ; m <= std::min(almn.Mmax(),l); ++m) {
                    alm(l,m) += ((T) factor) * almn(l,m,n);
                }
            }
        }
    }
    alm2map(alm,map);
}

/*!
 * Computes a portion of the density field for a given set of \a almn coefficients by using
 * the Healpix pixelisation scheme. Starting at radial distance \a rmin, a healpix is reconstructed
 * and the first of the twelve base sectors of the map (therefore a 2D image of size \a nside x \a nside)
 * is place in \a cube.
 * The process is repeated \a nr times, each time increasing the value of the radial distance by \a rstep.
 * At the end, \a cube is a 3D array containing the same area of the density field projected on a sphere
 * for different radial values.
 *
 * \a dsbt is a previously initialized FastDSBT object to be reused by almn2rmap.
 * \a nside specifies the resolution of the healpix maps
 */
template<typename T> void almn2cube(FastDSBT &dsbt, const Almn<xcomplex<T> > &almn, T* cube, T rmin, T rstep,int64 nr, int64 nside) {

    Healpix_Map<T> map(nside,RING,SET_NSIDE);
    T Rmax = almn.Rmax();

    /*Setup a matrix with the pixel numbers corresponding to (x,y) coordinates*/
    int64** indexCoords = new int64 * [nside];
    for (int64 i = 0 ; i < nside; ++i) {
        indexCoords[i] = new int64 [nside];
    }

    indexCoords[0][0] = 0;
    for (int64 x = 0; x<nside; x++)
    {
        for (int64 y=0 ; y<nside; y++)
        {

            if (x % 2 == 0)
            {
                if (y % 2 == 0)
                {
                    indexCoords[x][y] = 4 * indexCoords[x/2][y/2];
                }
                else
                {
                    indexCoords[x][y] = 4 * indexCoords[x/2][y/2]+2;
                }
            }
            else
            {
                if (y % 2 == 0)
                {
                    indexCoords[x][y] = 4 * indexCoords[x/2][y/2]+1;
                }
                else
                {
                    indexCoords[x][y] = 4 * indexCoords[x/2][y/2]+3;
                }
            }
        }
    }


    for (int64 r=0; r < nr ; r++) {
        almn2rmap(dsbt, almn, map, rmin + r*rstep);

        map.swap_scheme();

        //Extract only the portion of the map we are going to use
        for (int64 x =0; x < nside; ++x) {
            for (int64 y = 0; y < nside; ++y) {
                cube[r*nside*nside + y*nside + x] = map[indexCoords[x][y]];
            }
        }
        std::cout << " r : " << r*rstep + rmin << std::endl;
        map.SetNside(nside,RING);
    }

}
/*! \} */
#endif
