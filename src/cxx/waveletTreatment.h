/*
 * waveletTreatment.h - This file is part of MRS3D
 * 
 * Copyright 2011 CEA
 * Contributor : François Lanusse (francois.lanusse@gmail.com)
 * 
 * Created on 17/05/11
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

#ifndef WAVELET_TREATMENT_H
#define WAVELET_TREATMENT_H

#include <healpix_map.h>
#include <healpix_map_fitsio.h>
#include <mpi.h>
#include <lsconstants.h>
#include <cmath>
#include "waveletDecomposition.h"
#include <algorithm>

using namespace std;

/*! Type of wavelet treatment.
 *  READ_WRITE : The coefficients are read and modified
 *  READ_ONLY  : The coefficients are only read, only an inverse Discrete Spherical Fourier-Bessel Transform is required.
 *  WRITE_ONLY : All the coefficients are modified independently of their previous value (e.g. set to zero), only a direct Discrete Spherical Fourier-Bessel Transform is required. */
enum WAVELET_TREATMENT_TYPE {
		READ_WRITE,
		READ_ONLY,
		WRITE_ONLY
};

/*! Generic class that represents a treatment on wavelet coefficients in \b waveletDecomposition */
template<typename T> class waveletTreatment{

	public:
	
	/*! Apply an operation on the healpix maps in the \a mapBuffer array. These maps correspond
	 * to "shells" of the spherical 3D grid used in the Discrete Spherical Fourier-Bessel Transform, their
	 * radii are specified in the same order in the array \a radii. The size of the arrays \a mapBuffer and 
	 * \a radii is the value returned by the function \b numberOfMaps (see below).
	 * \a j is the scale of the wavelet coefficients bein processed. */
	virtual void operation(arr<Healpix_Map<T>* > &mapBuffer,arr<T> &radii ,int64 j){
		return;
	}
	
	/*! Returns the type of wavelet treatment as defined in \b WAVELET_TREATMENT_TYPE.
	 * This can be used to avoid unnecessary computations. For instance, if \b operation() only involves
	 * reading the values of the wavelet coefficients, they will not be modified and there is no need to update their almns. */
	virtual int type(){
		return READ_WRITE;
	}
	
	/*! The value returned by this method indicates to the method \b applyTreatment in \b waveletDecomposition,
	 * how many "shells" of the spherical 3D grid should be given at each call of \b operation().
	 * This is usefull for instance if the operation on the wavelet coefficients requires to now the values of the 
	 * neighbouring coefficients along the radial direction. */
	virtual int numberOfMaps(){
		return 1;
	}
	
	/*! Specifies the overlap in the \a mapBuffer given to \b operation between two calls.
	 * 0 means that \a mapBuffer will contain different maps. */
	virtual int overlap(){
		return 0;
	}
};

/*! Implementation of the abstract waveletTreatment class.
 *  \b backProjection will set all the wavelet coefficients to 0 except one.
 *  This can be used to illustrate the shape of a wavelet */
template<typename T> class backProjection : public waveletTreatment<T>{
	
public:
	/*! Constructor used to specify which coefficient should be non zero.
	 *  \a r, \a theta, \a phi give the spherical coordinates of the wavelet coefficient
	 *  and \a j specifies on which wavelet scale it is. */
	backProjection(T r, T theta, T phi, int64 j) : 
	r_(r),j_(j),angle_(theta,phi),lastR(0.0)
	{}
	
	
	void operation(arr<Healpix_Map<T> *> &mapBuffer,arr<T> &radii,int64 j){
	  //It is not necessary to set all the coefficients ot zero here, this is done by default in waveletDecomposition when
	  // the treatment applied is of type WRITE_ONLY.
		if(j == j_){
				if (r_ > lastR  && r_ <= radii[0]){
					(* mapBuffer[0])[mapBuffer[0]->ang2pix(angle_)] = 10.0;
				}
				lastR = radii[0];
		}

	}
	
	int type(){
		return WRITE_ONLY;
	}
	
	T r_;
	int64 j_;
	pointing angle_;
	
private:
	T lastR;
};

/*! Implementation of the abstract waveletTreatment class.
 *  \b HardThresholding will apply a threshold on every wavelet coefficients */
template<typename T> class HardThresholding : public waveletTreatment<T>{

public:
	/*! \a k is the threshold parameter and \a sigma is the noise variance */
	HardThresholding(T k, T sigma):
	k_(k),sigma_(sigma)
	{
	}
	
	void operation(arr<Healpix_Map<T> *> &mapBuffer,arr<T> &radii,int64 j){
		for (int i=0; i < numberOfMaps(); ++i) {
			#pragma omp parallel for
		for (int64 n=0; n < mapBuffer[i]->Npix(); ++n) {
			if (fabs((* mapBuffer[i])[n]) < k_*sigma_) {
				(* mapBuffer[i])[n]=0.0;
			}
		}
		}
		
	}
	
	int type(){
		return READ_WRITE;
	}
	
	virtual int numberOfMaps(){
		return 1;
	}
	
	virtual int overlap(){
		return 0;
	}
	
private:
	T k_;
	T sigma_;
};

/*! Implementation of the abstract waveletTreatment class.
 *  exportMaps can be used for debugging, it will export every "shell" of the spherical 3D grid
 *  in a healpix fits file */
template<typename T> class exportMaps : public waveletTreatment<T>{

public:
	/*! \a filenameBase is the generic name for the exported maps */
	exportMaps(std::string filenameBase){
		filenameBase_ = filenameBase;
	}
	
	void operation(arr<Healpix_Map<T> *> &mapBuffer,arr<T> &radii,int64 j){
		if (MPI::COMM_WORLD.Get_rank() == 0) {
		std::stringstream filename;
		filename << filenameBase_ <<  "_" << (int) radii[0] << "_" << j <<".fits" ;
			write_Healpix_map_to_fits(filename.str(), * (mapBuffer[0]), planckType<float>());
		}
	}
	
	int type(){
		return READ_ONLY;
	}
	
	int numberOfMaps(){
		return 1;
	}
	
	int overlap(){
		return 0;
	}
	
private:
	std::string filenameBase_;
};


#endif