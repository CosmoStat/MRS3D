/*
 * almn_field_tools.h - This file is part of MRS3D
 *
 * Copyright 2011 CEA
 * Contributor : François Lanusse (francois.lanusse@gmail.com)
 *
 * Created on 20/06/11
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

#ifndef ALMN_FIELD_TOOLS_H
#define ALMN_FIELD_TOOLS_H

#include "field.h"
#include "almn.h"
#include "string.h"
#include "almn_healpix_tools.h"
#include <xcomplex.h>
#include <healpix_map.h>
#include <healpix_map_fitsio.h>
#include "fastDSBT.h"
#include <mpi.h>

using namespace MPI;

/*! \defgroup almn_field_groop Conversion between cartesian fields and a_lmn */
/*! /{ */

/*! Computes a_lmn coefficients \a almn from a given cartesian density field \a field.
 *  \a tableFile refers to a fits file containing q_ln coefficients.
 *  \a Healpix_Nside specifies the resolution of the Healpix maps to be used by the algorithm.
 *  \a useMPI, if set to \b true enables the use of MPI.
 */
template<typename T> void field2almn(std::string tableFile , Field<T> &field, Almn<xcomplex<T> > &almn, int64 Healpix_Nside, bool useMPI=false) {

    int64 index;
    int64 lmax = almn.Lmax();
    int64 mmax = almn.Mmax();
    int64 nmax = almn.Nmax();
    T Rmax = almn.Rmax();

    arr<double> weight(2*Healpix_Nside);
    weight.fill(1.0);
    Alm<xcomplex<T> > alm(lmax,mmax);
    Healpix_Map<T> map(Healpix_Nside, RING,SET_NSIDE);

    int NbProc = 1;
    int ProcId = 0;

    if (useMPI) {
        std::cout << "Using MPI !" << std::endl;
        NbProc = COMM_WORLD.Get_size();
        ProcId = COMM_WORLD.Get_rank();
    }
    int nmax_local = nmax/NbProc;
    int nOffset = ProcId*nmax_local;

    //Initialize FastDSBT
    FastDSBT dsbt(nmax,lmax,tableFile.c_str());
    double Kmax;
    double normfactor = 1.0/(Rmax*Rmax*Rmax);

    //Temporary array to compute the Discrete Spherical Bessel Transform
    xcomplex<T> *out = (xcomplex<T> *) malloc(almn.Size()/NbProc*sizeof(xcomplex<T>));
    double      *transformMatrix = (double *) malloc(nmax*nmax_local*(lmax+1)*sizeof(double));
    int l_ext=0;

    // First step : reverse DSBT for all the almns (if actually needed)
    normfactor = 1.0/(Rmax*Rmax*Rmax);


    for (int64 n=0; n < nmax_local; ++n) {

        T rho = dsbt.get_Rln(l_ext,n + nOffset,Rmax);
        //We build an healpix map from the field for the given distance
#pragma omp parallel for
        for (int64 npix=0; npix < map.Npix(); ++npix) {
            pointing angle = map.pix2ang(npix);
            map[npix] = field.getValueFromSphCoords(rho, angle.theta, angle.phi);
        }

        //Now compute alms
        map2alm(map, alm, weight, false);

        for (int64 l = 0; l <= lmax;  ++l) {
            for (int64 m = 0; m <= min(mmax,l); ++m) {
                index = l*(l+1)/2 +m;
                out[index*nmax_local + n] = alm(l,m);
            }
        }
    }

    // Last step : DSBT to compute only the right a(l_ext)mn for which we have the correct grid
    double *matrix = dsbt.getTransformMatrix(l_ext,l_ext);

    Kmax = dsbt.Kmax(l_ext,Rmax);
    normfactor = 1.0/(Kmax*Kmax*Kmax);

    for (int i=0; i < NbProc ; ++i) {
        for (int64 l=0; l <= lmax; ++l) {

#pragma omp parallel for schedule(dynamic,1) private(index)
            for (int64 m = 0; m <= min(mmax,l); ++m) {
                index = l*(l+1)/2 +m;
                for (int64_t p=0; p < nmax_local; ++p) {
                    if (i == 0) {
                        almn(l,m,p+nOffset) = 0.0;
                    }
                    for (int64_t q=0; q < nmax_local; ++q) {
                        almn(l,m,p+nOffset)  +=  out[index*nmax_local + q] * matrix[(p + nOffset)*nmax + q + ((ProcId + NbProc - i )%NbProc )*nmax_local];
                    }
                    if (i == NbProc - 1) {
                        almn(l,m,p+nOffset) = almn(l,m,p+nOffset)*normfactor;
                    }
                }
            }
        }

        if (useMPI) {
            COMM_WORLD.Sendrecv_replace(out , sizeof(xcomplex<T>)*almn.Size() /NbProc, BYTE, (ProcId + 1 )% NbProc , i, (ProcId - 1 + NbProc) % NbProc, i);
        }

    }

    //Now we convert all the coefficients to the correct l

    //Before hand, we compute all the necessary matrices.
    for (int64 l=0; l<(lmax+1)/NbProc; ++l) {
        double *matrix = dsbt.getConversionMatrix(l_ext,l+ ProcId*(lmax+1)/NbProc);
	
	if(useMPI){
        for (int64 i=0 ; i < NbProc ; ++i) {
            COMM_WORLD.Sendrecv(&(matrix[((i + ProcId) % NbProc)*nmax*nmax_local]), nmax*nmax_local, DOUBLE, (i+ProcId) % NbProc, i, &(transformMatrix[(((ProcId + NbProc-i)%NbProc)*(lmax+1)/NbProc + l)*nmax*nmax_local]), nmax*nmax_local, DOUBLE, (ProcId + NbProc-i)%NbProc, i);
        }
	}else{
	    memcpy(&(transformMatrix[l*nmax*nmax]),matrix,nmax*nmax);
	}
    }

    for (int i=0; i<NbProc; ++i) {

#pragma omp parallel for schedule(dynamic,1) private(index)
        for (int64 l = 0; l <= lmax;  ++l) {

            for (int64 m = 0; m <= min(mmax,l); ++m) {
                index = l*(l+1)/2 +m;
                for (int64_t p=0; p < nmax_local; ++p) {
                    if (i == 0) {
                        out[index*nmax_local + p] = 0.0;
                    }

                    for (int64_t q=0; q < nmax_local; ++q) {
                        out[index*nmax_local + p] += almn(l,m,q+nOffset) * transformMatrix[l*nmax*nmax_local + p*nmax + q  + ((ProcId - i + NbProc) % NbProc)*nmax_local];
                    }
                }
            }

        }
        if (useMPI) {
            COMM_WORLD.Sendrecv_replace(&(almn(0,0,nOffset)), sizeof(xcomplex<T>)*almn.Size()/NbProc, BYTE, (ProcId + 1 )% NbProc , i, (ProcId - 1 + NbProc) % NbProc, i);
        }
    }

    for (int64 l=0; l <= lmax ; ++l) {
        for (int64 m=0; m <= min(mmax,l); ++m) {
            index = l*(l+1)/2 +m;
            for (int64 n =0; n < nmax_local; ++n) {
                almn(l,m,n+nOffset) = out[index*nmax_local + n];
            }
        }
    }

    if (useMPI) {
        //Recover results from all the processes
        for (int i = 1; i < NbProc; ++i) {
            COMM_WORLD.Sendrecv(&(almn(0,0,nOffset)), almn.Size()/NbProc*sizeof(xcomplex<T>), BYTE, (ProcId+i)%NbProc, i, &(almn(0,0, nmax_local*((ProcId - i + NbProc) % NbProc))) ,almn.Size()/NbProc*sizeof(xcomplex<T>), BYTE, (ProcId - i + NbProc) % NbProc, i);
        }
    }

    free(out);
    free(transformMatrix);

}

/*! Reconstructs a cartesian density field \a field from given a_lmn coefficients \a almn.
 *  \a tableFile refers to a fits file containing q_ln coefficients.
 *  \a Healpix_Nside specifies the resolution of the Healpix maps to be used by the algorithm.
 *  \a useMPI, if set to \b true enables the use of MPI.
 */
template<typename T> void almn2field(std::string tableFile ,  Almn<xcomplex<T> > &almn, Field<T> &field, int64 Healpix_Nside, bool useMPI=false) {
    Healpix_Map<T> map1(Healpix_Nside,RING,SET_NSIDE);
    Healpix_Map<T> map2(Healpix_Nside,RING,SET_NSIDE);
    FastDSBT dsbt(almn.Nmax(),almn.Lmax(),tableFile.c_str());
    int Nside = field.Nside;
    T r1 = 0;
    T r2 = 0;

    for (int64 n=0; n < 2*Nside; ++n) {
        almn2rmap(dsbt, almn, map2, ((T) n)*field.size/((T) 2*Nside));
        r2 =  ((T) n)*field.size/((T) 2*Nside);

#pragma omp parallel for
        for (int x=0;  x < Nside; ++x) {
            for (int y=0; y < Nside; ++y) {
                for (int z=0; z < Nside; ++z) {

                    T r, theta, phi;
                    field.getSphericalCoordinates(x,y,z,r,theta,phi);

                    if (r > r1 && r <= r2) {
                        int64 npix   = map2.ang2pix(pointing(theta,phi));
                        T u ;
                        if (n == 0) {
                            u = 0;
                        } else {
                            u = (r - r1)/(r2 - r1);
                        }

                        field(x,y,z) = map1[npix]*( 1.0 - u) + map2[npix]*u;
                    }
                }
            }
        }

        map1.Import_nograde(map2);
        r1 = r2;
        std::cout << "n > : " << n << std::endl;
    }
}

/*! \} */


#endif
