/*
 * waveletDecomposition.h - This file is part of MRS3D
 * 
 * Copyright 2011 CEA
 * Contributor : François Lanusse (francois.lanusse@gmail.com)
 * 
 * Created on 14/04/11
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


#ifndef WAVELET_DECOMPOSITION_H
#define WAVELET_DECOMPOSITION_H

#include "almn.h"
#include "almn_fitsio.h"
#include <xcomplex.h>
#include <healpix_map.h>
#include <cstdlib>
#include "almn_healpix_tools.h"
#include "fastDSBT.h"
#include <planck_rng.h>
#include "waveletTreatment.h"
#include <math.h>
#include <string.h>
#include <mpi.h>
#include <omp.h>
#include <algorithm>

using namespace MPI;
using namespace std;

/*! Class for the processing and storage of the wavelet decomposition. */
template<typename T> class waveletDecomposition {

public:
    int64 nside,lmax,mmax,nmax,Jmax;
    double Rmax;
    T kc;
    std::string tableFile_;

    /*! Computes the wavelet decomposition from the input \a almn. 
     *  \a nside_ specifies the resolution of HealpixMaps used to process the Discrete Spherical Fourier-Bessel Transform
     * when applying treatments on the wavelet coefficients.
     *  \a Jmax_ is the number of wavelet scales in the decomposition
     *  \a kc_   is the cutoff frequency of the lowPass filter used in the decomposition 
     * 	\a tableFile is the name of the fits file of tabulated Bessel function zeros
     *  \a useMPI_ indicates whether the MPI implementation should be used */
    waveletDecomposition(int64 nside_, int64 Jmax_, T kc_, Almn<xcomplex<T> > &almn, std::string tableFile="", bool useMPI_=false)
            : nside(nside_),lmax(almn.Lmax()),mmax(almn.Mmax()),nmax(almn.Nmax()),Rmax(almn.Rmax()),nmax_local(nmax),
            Jmax(Jmax_),kc(kc_),wlmns(Jmax_+1),
            useMPI(useMPI_),
            dsbt(nmax,lmax,tableFile.c_str()),
            LowPass(Jmax_,lmax+1,1),
            Normalization(Jmax_+1)
    {
        tableFile_ = tableFile;

        //	if (useMPI) {
        NbProc = COMM_WORLD.Get_size();
        ProcId = COMM_WORLD.Get_rank();
        //	}
        nmax_local = nmax/NbProc;
        nOffset = ProcId*nmax_local;

        // Setup the sizes of the wlmns
        for (int64 j=0; j <= Jmax; ++j) {
            wlmns[j].Set(lmax,mmax,nmax_local,Rmax);
        }
        
     //   kc = 2.0*dsbt.get_Kln(0,nmax-1,Rmax);

        // Computes lowPass filters
        computeLowPass();
	
	//Copy almn to wlmns(0,0)
	#pragma omp parallel for
        for (int64 n=0; n < nmax_local; ++n) {
            for (int64 l=0 ; l <= lmax ; l++) {
                for (int64 m=0; m <= l ; ++m) {
                    wlmns[0](l,m,n) = almn(l,m,n+nOffset);
                }
            }
        }
	
	//Process the wavelet decomposition
	computeTransform();
	
    }

    /*! Performs the inverse wavelet decomposition and returns in \a almn 
     *  the spherical Fourier-Bessel coefficients of the reconstructed data. */
    void reconstructionAlmn(Almn<xcomplex<T> > &almn) {
      almn.Set(lmax,mmax,nmax,Rmax);

#pragma omp parallel
        {
            //Copy rlmn to almn
#pragma omp for
            for (int64 n=0; n < nmax_local; ++n) {
                for (int64 l=0 ; l <= lmax ; ++l) {
                    for (int64 m=0; m <= l ; ++m) {
                        almn(l,m,n+nOffset) = wlmns[Jmax](l,m,n);
                    }
                }
            }

            for (int64 j=Jmax-1; j >= 0; --j) {
#pragma omp for
                for (int64 n=0; n < nmax_local; ++n) {
                    for (int64 l=0 ; l <= lmax ; ++l) {
                        T lowpass = LowPass(j,l,n);
                        for (int64 m=0; m <= l ; ++m) {
                            almn(l,m,n+nOffset) = (almn(l,m,n+nOffset) * lowpass + (((T) 1.0) - lowpass)*wlmns[j](l,m,n)*sqrt(Normalization[j]))/(lowpass*lowpass + (((T) 1.0) - lowpass)*(((T) 1.0) - lowpass));
                        }
                    }
                }
            }
        }

        if (useMPI) {
            //Recover results from all the processes
            for (int i = 1; i < NbProc; ++i) {
                COMM_WORLD.Sendrecv(&(almn(0,0,nOffset)), wlmns[0].Size()*sizeof(xcomplex<T>), BYTE, (ProcId+i)%NbProc, i, &(almn(0,0, nmax_local*((ProcId - i + NbProc) % NbProc))) ,wlmns[0].Size()*sizeof(xcomplex<T>), BYTE, (ProcId - i + NbProc) % NbProc, i);
            }
        }
    }

    /*! Returns in \a almn a copy of the Spherical Fourier-Bessel coefficients of the wavelet scale \a j when
     *  \a j varies from \b 0 to \b Jmax-1. For \b j = \b Jmax, this function returns the coefficient of the smoothed 
     *  away approximation. */  
    void getWlmn(Almn<xcomplex<T> > &almn, int64 j) {
        almn.Set(lmax,mmax,nmax,Rmax);
      
#pragma omp parallel for
        for (int64 n=0; n < nmax_local; ++n) {
            for (int64 l=0 ; l <= lmax ; ++l) {
                for (int64 m=0; m <= l ; ++m) {
                    almn(l,m,n+nOffset) = wlmns[j](l,m,n);
                }
            }
        }

        if (useMPI) {
            //Recover results from all the processes
            for (int i = 1; i < NbProc; ++i) {
                COMM_WORLD.Sendrecv(&(almn(0,0,nOffset)), wlmns[0].Size()*sizeof(xcomplex<T>), BYTE, (ProcId+i)%NbProc, i, &(almn(0,0, nmax_local*((ProcId - i + NbProc) % NbProc))) ,wlmns[0].Size()*sizeof(xcomplex<T>), BYTE, (ProcId - i + NbProc) % NbProc, i);
            }
        }
    }

    /*! Apply a treatment specified through the object \a treatment to the wavelet scale \a j .
     * \a Nside specifies the resolution of Healpix maps used during the Discrete Spherical Fourier-Bessel Transform.
     * waveletTreatment is a generic class which can be specialised to implement a specific operation. See waveletTreatment.h */
    void applyTreatment(int64 Nside, int64 j, waveletTreatment<T> *treatment) {
        int64 index;
        arr<double> weight(2*Nside);
        weight.fill(1.0);
        Alm<xcomplex<T> > alm(lmax,mmax);

        //Initialize FastDSBT
        double Kmax;
        double normfactor = 1.0/(Rmax*Rmax*Rmax);

        //Temporary array to compute the Discrete Spherical Bessel Transform
        xcomplex<T> *out = (xcomplex<T> *) malloc(wlmns[0].Size()*sizeof(xcomplex<T>));
        double      *transformMatrix = (double *) malloc(nmax*nmax_local*(lmax+1)*sizeof(double));

        //Map Buffer for block processing
        arr<Healpix_Map<T> *> mapBuffer(treatment->numberOfMaps());
        for (int i =0; i < treatment->numberOfMaps(); ++i) {
            mapBuffer[i] = new Healpix_Map<T>(Nside,RING,SET_NSIDE);
        }

        // First step : reverse DSBT for all the almns (if actually needed)
        normfactor = 1.0/(Rmax*Rmax*Rmax);
        int64 l_ext = 0;

        if (treatment->type() != WRITE_ONLY) {

            //Before hand, we compute all the necessary matrices : saves time ^^ but not space :(
            for (int64 l=0; l<(lmax+1)/NbProc; ++l) {
                double *matrix = dsbt.getTransformMatrix(l+ ProcId*(lmax+1)/NbProc,l_ext);

                for (int64 i=0 ; i < NbProc ; ++i) {
                    COMM_WORLD.Sendrecv(&(matrix[((i + ProcId) % NbProc)*nmax*nmax_local]), nmax*nmax_local, DOUBLE, (i+ProcId) % NbProc, i, &(transformMatrix[(((ProcId + NbProc-i)%NbProc)*(lmax+1)/NbProc + l)*nmax*nmax_local]), nmax*nmax_local, DOUBLE, (ProcId + NbProc-i)%NbProc, i);
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
                                out[index*nmax_local + p] += (wlmns[j])(l,m,q) * transformMatrix[l*nmax*nmax_local + p*nmax + q  + ((ProcId - i + NbProc) % NbProc)*nmax_local];
                            }
                            if (i == NbProc - 1) {
                                out[index*nmax_local + p] = out[index*nmax_local + p]*normfactor;
                            }
                        }
                    }

                }
                if (useMPI) {
                    COMM_WORLD.Sendrecv_replace(&(wlmns[j](0,0,0)), sizeof(xcomplex<T>)*wlmns[j].Size(), BYTE, (ProcId + 1 )% NbProc , i, (ProcId - 1 + NbProc) % NbProc, i);
                }
            }
        }

        for (int64 i=0; i < nmax_local / treatment->numberOfMaps(); ++i) {

            int64 blockOffset = i* treatment->numberOfMaps();
            arr<T> radii(treatment->numberOfMaps());

            //Now, out contains the alm(r), let's rebuild an actual map, a block at a time
            for (int64 n=0; n < treatment->numberOfMaps(); ++n) {
                if (treatment->type() != WRITE_ONLY) {
                    //First let's put back the alm in an Alm object to use Healpix
                    for (int64 l = 0; l <= lmax;  ++l) {
                        for (int64 m = 0; m <= min(mmax,l); ++m) {
                            index = l*(l+1)/2 +m;
                            alm(l,m) = out[index*nmax_local + n + blockOffset];
                        }
                    } 
                    alm2map(alm, * mapBuffer[n]);
                } else {
                    mapBuffer[n]->fill(0.0);
                }
                radii[n] = dsbt.get_Rln(l_ext,n + nOffset + blockOffset,Rmax);
            }
	    
            //Here if we feel like it we can actually apply treatments to the map, why not...
           treatment->operation(mapBuffer,radii,j);

            if (treatment->type() != READ_ONLY) {
                for (int64 n=0; n < treatment->numberOfMaps() ; ++n) {

                    //Now inverse operation
                    map2alm( * mapBuffer[n], alm, weight, false);

                    for (int64 l = 0; l <= lmax;  ++l) {
                        for (int64 m = 0; m <= min(mmax,l); ++m) {
                            index = l*(l+1)/2 +m;
                            out[index*nmax_local + n + blockOffset] = alm(l,m);
                        }
                    }
                }
            }
        }
	

        // Last step : DSBT to compute only the right a(l_ext)mn for which we have the correct grid
        if (treatment->type() != READ_ONLY) {
            double *matrix = dsbt.getTransformMatrix(l_ext,l_ext);

            Kmax = dsbt.Kmax(l_ext,Rmax);
            normfactor = 1.0/(Kmax*Kmax*Kmax);

            for (int i=0; i < NbProc ; ++i) {

#pragma omp parallel for schedule(dynamic,1) private(index)
                for (int64 l = 0;  l <= lmax; ++l) {


                    for (int64 m = 0; m <= min(mmax,l); ++m) {
                        index = l*(l+1)/2 +m;
                        for (int64_t p=0; p < nmax_local; ++p) {
                            if (i == 0) {
                                (wlmns[j])(l,m,p) = 0.0;
                            }
                            for (int64_t q=0; q < nmax_local; ++q) {
                                (wlmns[j])(l,m,p)  +=  out[index*nmax_local + q] * matrix[(p + nOffset)*nmax + q + ((ProcId + NbProc - i )%NbProc )*nmax_local];
                            }
                            if (i == NbProc - 1) {
                                (wlmns[j])(l,m,p) = (wlmns[j])(l,m,p)*normfactor;
                            }
                        }
                    }
                }

                if (useMPI) {
                    COMM_WORLD.Sendrecv_replace(out , sizeof(xcomplex<T>)*wlmns[j].Size() , BYTE, (ProcId + 1 )% NbProc , i, (ProcId - 1 + NbProc) % NbProc, i);
                }

            }


	    // Now we convert all the coefficients to the correct l

           // Before hand, we compute all the necessary matrices
            for (int64 l=0; l<(lmax+1)/NbProc; ++l) {
                double *matrix = dsbt.getConversionMatrix(l_ext,l+ ProcId*(lmax+1)/NbProc);

                for (int64 i=0 ; i < NbProc ; ++i) {
                    COMM_WORLD.Sendrecv(&(matrix[((i + ProcId) % NbProc)*nmax*nmax_local]), nmax*nmax_local, DOUBLE, (i+ProcId) % NbProc, i, &(transformMatrix[(((ProcId + NbProc-i)%NbProc)*(lmax+1)/NbProc + l)*nmax*nmax_local]), nmax*nmax_local, DOUBLE, (ProcId + NbProc-i)%NbProc, i);
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
                                out[index*nmax_local + p] += (wlmns[j])(l,m,q) * transformMatrix[l*nmax*nmax_local + p*nmax + q  + ((ProcId - i + NbProc) % NbProc)*nmax_local];
                            }
                        }
                    }

                }
                if (useMPI) {
                    COMM_WORLD.Sendrecv_replace(&(wlmns[j](0,0,0)), sizeof(xcomplex<T>)*wlmns[j].Size(), BYTE, (ProcId + 1 )% NbProc , i, (ProcId - 1 + NbProc) % NbProc, i);
                }
            }

	    //  Ok, now we put back the coefficients in wlmns
            for (int64 l = 0; l <= lmax;  ++l) {

                for (int64 m = 0; m <= min(mmax,l); ++m) {
                    index = l*(l+1)/2 +m;
                    for (int64_t p=0; p < nmax_local; ++p) {
                        wlmns[j](l,m,p) = out[index*nmax_local + p] ;
			 xcomplex<T> toto = (wlmns[j])(l,m,p) ;
                    }
                }
            }
      }


        free(out);
        free(transformMatrix);

        for (int64 i=0; i < treatment->numberOfMaps(); ++i) {
            delete mapBuffer[i];
        }
    }

private:

    bool useMPI;
    int NbProc;
    int ProcId;
    int64 nmax_local;
    int64 nOffset;

    arr<Almn<xcomplex<T> > > wlmns;

    FastDSBT dsbt;

    void computeTransform() {
#pragma omp parallel
        {
            for (int64 j=0; j < Jmax ; ++j) {
#pragma omp for
                for (int64 l=0 ; l <= lmax ; ++l) {
                    for (int64 n=0; n < nmax_local; ++n) {
                        T lowpass = LowPass(j,l,n);
                        for (int64 m=0; m <= l ; ++m) {
                            (wlmns[j+1])(l,m,n) = (wlmns[j])(l,m,n) * lowpass;
			    if(Normalization[j] > 0.0){
			      (wlmns[j])(l,m,n) = (wlmns[j])(l,m,n) * (((T)1.0) - lowpass)/sqrt(Normalization[j]);
			    }
                        }
                    }
                }
            }

        }
    }
    
    /*! Computes the lowPass filter used in the wavelet decomposition.
     *  To implement another wavelet transform, this function should be modified to
     *  compute a different filter. */
    void computeLowPass() {
        double x, phi_top, phi_bot;

        LowPass.alloc(Jmax+1,lmax+1,nmax_local);
        Normalization.alloc(Jmax+1);
        Normalization.fill(0.0);

        for (int j=0; j <= Jmax;  j++) {
            double cutoff = kc /((double) pow(2.0,j+1.0));

            for (int l=0; l <= lmax; ++l) {
  
                for (int n=0; n < nmax_local ; ++n) {

                    x =  2.0*((double) dsbt.get_Kln(l,n+nOffset,Rmax))/cutoff;

                    phi_top =  1.0/8.0 * ( pow(fabs(x - 2.0),3.0) - 4.0*pow( fabs(x - 1.0),3.0) + 6.0 *pow( fabs(x),3.0)
                                           - 4.0 * pow( fabs(x + 1.0),3.0) + pow( fabs(x + 2.0),3.0));

                    x = ((double) dsbt.get_Kln(l,n+nOffset,Rmax))/((double) cutoff);
                    phi_bot =  1.0/8.0 * ( pow( fabs(x - 2.0),3.0) - 4.0*pow( fabs(x - 1.0),3.0) + 6.0 *pow( fabs(x),3.0)
                                           - 4.0 * pow( fabs(x + 1.0),3.0) + pow( fabs(x + 2.0),3.0));
                    if (x > 1.00 ) {
                        LowPass(j,l,n) = 0.0;
                    } else {
                        LowPass(j,l,n) = phi_top/phi_bot;
                        if (LowPass(j,l,n) <=0.0) {
                            LowPass(j,l,n) =0.0;
                        }
                    }

                    if (j==0) {
                        Normalization[0] += pow(fabs(1.0 - LowPass(0,l,n)),2);
                    } else {
                        Normalization[j] += pow(fabs(LowPass(j-1,l,n) - LowPass(j,l,n)), 2);
                    }

                }

            }
            Normalization[j] = Normalization[j]/(12.0*nside*nside/(4.0*M_PI));
	    MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE, &(Normalization[j]),1,MPI_DOUBLE,MPI_SUM);
        }

    }

    arr3<T>			LowPass;
    arr<T>			Normalization;

};

#endif
