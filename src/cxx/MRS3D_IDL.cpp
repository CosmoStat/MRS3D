/*
 * MRS3D.cpp - This file is part of MRS3D
 * 
 * Copyright 2011 CEA
 * Contributor : François Lanusse (francois.lanusse@gmail.com)
 * 
 * Created on 11/04/11
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

#include "almn.h"
#include "almn_fitsio.h"
#include "cln.h"
#include "cln_fitsio.h"
#include "almn_cln_tools.h"
#include "waveletDecomposition.h"
#include "almn_healpix_tools.h"
#include "field.h"
#include "field_fitsio.h"
#include "almn_field_tools.h"
#include "selectionFunction.h"
#include <cstdlib>
#include <cstdio>
#include "idl_export.h"

/*! \defgroup dsbt_idl_group IDL interface to FastDSBT */
/*! \{ */

extern "C" {
	  
	/*! Save almn coefficients to a fits file.
	 * Arguments in argv:
	 * 0 : IDL_STRING *\a almnFile  -> name of the file 
	 * 1 : IDL_LONG64 *\a lmax      -> lmax of almn array
	 * 2 : IDL_LONG64 *\a mmax      -> mmax of almn array
	 * 3 : IDL_LONG64 *\a nmax      -> nmax of almn array
	 * 4 : float      *\a rmax      -> value of Rmax used to compute almns
	 * 5 : float      *\a almn_real -> pointer to the real component of almn 
	 * 6 : float      *\a almn_imag -> pointer to the imaginary component of almn */
	void MRS3D_saveAlmn_IDL(int argc, void *argv[]){
		char buffer[256];
		sprintf(buffer, "Error in MRS3D_loadAlmn_IDL : Invalid number of arguments %d instead of 7\n",argc);
		if(argc != 7){
			IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, buffer);
		}
		
		IDL_STRING *almnFile = (IDL_STRING *) argv[0];
		IDL_LONG64 *lmax = (IDL_LONG64 *) argv[1];
		IDL_LONG64 *mmax = (IDL_LONG64 *) argv[2];
		IDL_LONG64 *nmax = (IDL_LONG64 *) argv[3];
		float     *rmax  = (float *) argv[4];
		float *almn_real = (float *) argv[5];
		float *almn_imag = (float *) argv[6];
		
		Almn<xcomplex<float> > almn(*lmax,*mmax,*nmax,*rmax);
		
		for (int64 n=0; n < *nmax; ++n) {
			for (int64 l=0; l <= *lmax; ++l) {
				for (int64 m=0; m <= min((IDL_LONG64)l , *mmax); ++m) {
					almn(l,m,n).re = almn_real[n*(*lmax+1)*(*mmax+1) + m*(*lmax+1) + l];
					almn(l,m,n).im = almn_imag[n*(*lmax+1)*(*mmax+1) + m*(*lmax+1) + l];
				}
			}
		}
		write_Almn_to_fits(std::string(almnFile->s), almn, *lmax, *mmax, *nmax, planckType<float>());
	}
	
	/*! load almn coefficients from a fits file.
	 * Arguments in argv:
	 * 0 : IDL_STRING *\a almnFile  -> name of the file 
	 * 1 : IDL_LONG64 *\a lmax      -> lmax of almn array
	 * 2 : IDL_LONG64 *\a mmax      -> mmax of almn array
	 * 3 : IDL_LONG64 *\a nmax      -> nmax of almn array
	 * 4 : float      *\a almn_real -> pointer to the real component of almn
	 * 5 : float      *\a almn_imag -> pointer to the imaginary component of almn */
	void MRS3D_loadAlmn_IDL(int argc, void *argv[]){
		char buffer[256];
		sprintf(buffer, "Error in MRS3D_loadAlmn_IDL : Invalid number of arguments %d instead of 6\n",argc);
		if(argc != 6){
			IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, buffer);
		}
		
		IDL_STRING *almnFile = (IDL_STRING *) argv[0];
		IDL_LONG64 *lmax	 = (IDL_LONG64 *) argv[1];
		IDL_LONG64 *mmax	 = (IDL_LONG64 *) argv[2];
		IDL_LONG64 *nmax	 = (IDL_LONG64 *) argv[3];
		float *almn_real    = (float *) argv[4];
		float *almn_imag    = (float *) argv[5];
		
		Almn<xcomplex<float> > almn;
		
		read_Almn_from_fits(std::string(almnFile->s), almn, *lmax, *mmax, *nmax);
				
		for (int64 n=0; n < *nmax; ++n) {
			for (int64 l=0; l <= *lmax; ++l) {
				for (int64 m=0; m <= *mmax; ++m) {
					almn_real[n*(*lmax+1)*(*mmax+1) + m*(*lmax+1) + l] = almn(l,m,n).re;
					almn_imag[n*(*lmax+1)*(*mmax+1) + m*(*lmax+1) + l] = almn(l,m,n).im;
				}
			}
		}
	}
	
	/*! load cln coefficients from a fits file.
	 * Arguments in argv:
	 * 0 : IDL_STRING *\a clnFile   -> name of the file 
	 * 1 : IDL_LONG64 *\a lmax      -> lmax of cln array
	 * 2 : IDL_LONG64 *\a nmax      -> nmax of cln array
	 * 3 : float      *\a data      -> pointer to the array in which to store the cln coefficients */
	void MRS3D_loadCln_IDL(int argc, void *argv[]){
		char buffer[256];
		sprintf(buffer, "Error in MRS3D_loadCln_IDL : Invalid number of arguments %d instead of 4\n",argc);
		if(argc != 4){
			IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, buffer);
		}
		
		IDL_STRING *clnFile  = (IDL_STRING *) argv[0];
		IDL_LONG64 *lmax	 = (IDL_LONG64 *) argv[1];
		IDL_LONG64 *nmax	 = (IDL_LONG64 *) argv[2];
		float      *data	 = (float *) argv[3];
		
		Cln<float > cln(*lmax,*nmax);
		
		read_Cln_from_fits(std::string(clnFile->s), cln, *lmax, *nmax);
		
		for (int64 n=0; n < *nmax; ++n) {
			for (int64 l=0; l <= *lmax; ++l) {
				data[n*(*lmax+1) + l] = cln(l,n);
			}
		}
	}
	
	/*! Save cln coefficients to a fits file.
	 * Arguments in argv:
	 * 0 : IDL_STRING *\a clnFile   -> name of the file 
	 * 1 : IDL_LONG64 *\a lmax      -> lmax of cln array
	 * 2 : IDL_LONG64 *\a nmax      -> nmax of cln array
	 * 3 : float      *\a rmax	-> rmax used to dicretize the power spectrum
	 * 4 : float      *\a data      -> pointer to the array in which the cln coefficients are stored */
	void MRS3D_saveCln_IDL(int argc, void *argv[]){
		char buffer[256];
		sprintf(buffer, "Error in MRS3D_saveAlmn_IDL : Invalid number of arguments %d instead of 5\n",argc);
		if(argc != 5){
			IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, buffer);
		}
		
		IDL_STRING *clnFile  = (IDL_STRING *) argv[0];
		IDL_LONG64 *lmax	 = (IDL_LONG64 *) argv[1];
		IDL_LONG64 *nmax	 = (IDL_LONG64 *) argv[2];
		float *rmax		 = (float *) argv[3] ;
		float *data	 	 = (float *) argv[4];
		
		Cln<float > cln(*lmax,*nmax,*rmax);
		
		for (int64 n=0; n < *nmax; ++n) {
			for (int64 l=0; l <= *lmax; ++l) {
				cln(l,n) = data[n*(*lmax+1) + l];
			}
		}
		
		write_Cln_to_fits(std::string(clnFile->s), cln, *lmax, *nmax, planckType<float>());
	}
	
	/*! Loads a almn fits file and returns the power spectrum.
	 * Arguments in argv:
	 * 0 : IDL_STRING *\a AlmnFile  -> name of the file 
	 * 1 : IDL_LONG64 *\a lmax      -> lmax of almn to use
	 * 2 : IDL_LONG64 *\a mmax      -> nmax of almn to use
	 * 3 : IDL_LONG64 *\a nmax      -> nmax of almn to use
	 * 4 : float      *\a data      -> pointer to the array in which to store the cln coefficients */
	void MRS3D_AlmnFileToCln_IDL(int argc, void *argv[]){
		char buffer[256];
		sprintf(buffer, "Error in MRS3D_AlmnFileToCln_IDL : Invalid number of arguments %d instead of 5\n",argc);
		if(argc != 5){
			IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, buffer);
		}
		
		IDL_STRING *AlmnFile  = (IDL_STRING *) argv[0];
		IDL_LONG64 *lmax	 = (IDL_LONG64 *) argv[1];
		IDL_LONG64 *mmax	 = (IDL_LONG64 *) argv[2];
		IDL_LONG64 *nmax	 = (IDL_LONG64 *) argv[3];
		float *data		 = (float *) argv[4];
		
		Almn<xcomplex<float> > almn(*lmax,*mmax,*nmax);
		
		read_Almn_from_fits(std::string(AlmnFile->s), almn, *lmax, *mmax, *nmax);
		
		
		Cln<float > cln(*lmax,*nmax);
				
		extract_cln(almn, cln);
				
		for (int64 n=0; n < *nmax; ++n) {
			for (int64 l=0; l <= *lmax; ++l) {
				data[n*(*lmax+1) + l] = cln(l,n);
			}
		}
	}
	
	/*! Convert an almn fits file to a cln fits file.
	 * Arguments in argv:
	 * 0 : IDL_STRING *\a AlmnFile  -> name of the file 
	 * 1 : IDL_LONG64 *\a lmax      -> lmax of almn to use
	 * 2 : IDL_LONG64 *\a mmax      -> nmax of almn to use
	 * 3 : IDL_LONG64 *\a nmax      -> nmax of almn to use
	 * 4 : float      *\a data      -> pointer to the array in which to store the cln coefficients */
	void MRS3D_AlmnFileToClnFile_IDL(int argc, void *argv[]){
		char buffer[256];
		sprintf(buffer, "Error in MRS3D_AlmnFileToClnFile_IDL : Invalid number of arguments %d instead of 5\n",argc);
		if(argc != 5){
			IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, buffer);
		}
		
		IDL_STRING *AlmnFile  = (IDL_STRING *) argv[0];
		IDL_LONG64 *lmax	 = (IDL_LONG64 *) argv[1];
		IDL_LONG64 *mmax	 = (IDL_LONG64 *) argv[2];
		IDL_LONG64 *nmax	 = (IDL_LONG64 *) argv[3];
		IDL_STRING *ClnFile  = (IDL_STRING *) argv[4];
		
		Almn<xcomplex<float> > almn(*lmax,*mmax,*nmax);
		
		read_Almn_from_fits(std::string(AlmnFile->s), almn, *lmax, *mmax, *nmax);
		
		Cln<float > cln(*lmax,*nmax);
		
		extract_cln(almn, cln);
	
		write_Cln_to_fits(std::string(ClnFile->s), cln, *lmax, *nmax, planckType<float>());
	}
	
	
	/*! Extract the same portion of Healpix maps of the density field reconstructed from almns
	 *  and store them in a 3D array (each layer is a map).
	 *  Arguments in argv:
	 *  0 : IDL_STRING *\a tableFile  -> name of the qln table file
	 *  1 : IDL_STRING *\a AlmnFile   -> name of the almn fits file
	 *  2 : IDL_LONG64 *\a lmax       -> lmax of almn to use
	 *  3 : IDL_LONG64 *\a mmax       -> mmax of almn to use
	 *  4 : IDL_LONG64 *\a nmax       -> nmax of almn to use
	 *  5 : IDL_LONG64 *\a nside      -> nside parameter of the Healpix maps
	 *  6 : IDL_LONG64 *\a nr         -> number of healpix maps to compute
	 *  7 : float      *\a rmin       -> radius of the first healpix map
	 *  8 : float      *\a rstep      -> step between to consecutive radii
	 *  9 : float      *\a cube       -> pointer to the 3D array in which to store the healpix map (size: nside x nside x nr) */
	void MRS3D_ExtractCubeFromAlmn_IDL(int argc, void *argv[]){
		char buffer[256];
		sprintf(buffer, "Error in MRS3D_ExtractCubeFromAlmn_IDL : Invalid number of arguments %d instead of 10\n",argc);
		if(argc != 10){
			IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, buffer);
		}
		
		IDL_STRING *tableFile  = (IDL_STRING *) argv[0];
		IDL_STRING *AlmnFile  = (IDL_STRING *) argv[1];
		IDL_LONG64 *lmax	 = (IDL_LONG64 *) argv[2];
		IDL_LONG64 *mmax	 = (IDL_LONG64 *) argv[3];
		IDL_LONG64 *nmax	 = (IDL_LONG64 *) argv[4];
		IDL_LONG64 *nside    = (IDL_LONG64 *) argv[5];
		IDL_LONG64 *nr		 = (IDL_LONG64 *) argv[6];
		float      *rmin     = (float *)	  argv[7];
		float      *rstep	 = (float *)      argv[8];
		float      *cube     = (float *)      argv[9];
		
		Almn<xcomplex<float> > almn(*lmax,*mmax,*nmax);
		
		read_Almn_from_fits(std::string(AlmnFile->s), almn, *lmax, *mmax, *nmax);
		
		FastDSBT dsbt(*nmax, *lmax,tableFile->s);
		almn2cube(dsbt, almn, cube, *rmin,*rstep, *nr, *nside);
	}

	/*! Load cartesian field from a fits file.
	 * Arguments in argv:
	 * 0 : IDL_STRING *\a fieldFile -> name of the file 
	 * 1 : IDL_LONG64 *\a Nside     -> number of pixels along one side of the cube
	 * 2 : float      *\a field      -> pointer to the array in which to store the field */
	void MRS3D_loadField_IDL(int argc, void *argv[]){
		char buffer[256];
		sprintf(buffer, "Error in MRS3D_loadField_IDL : Invalid number of arguments %d instead of 3\n",argc);
		if(argc != 3){
			IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, buffer);
		}
		
		IDL_STRING *fieldFile = (IDL_STRING *) argv[0];
		IDL_LONG64 *Nside	  = (IDL_LONG64 *) argv[1];
		float      *field_     = (float *) argv[2];
		
		int actualNside;
		float size;
		
		get_fieldsize(std::string(fieldFile->s), actualNside, size);
		
		if (actualNside != *Nside) {
			sprintf(buffer, "Error in MRS3D_loadField_IDL : Invalid size of field. \n");
			IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, buffer);
		}
		
		Field<float> field(actualNside,size);
		
		read_Field_from_fits(std::string(fieldFile->s), field);
				
		for (int x=0; x < actualNside; ++x) {
			for (int y=0; y < actualNside; ++y) {
				for (int z=0; z < actualNside; ++z) {
					field_[z*(actualNside*actualNside) + y*(actualNside) + x] = field(x,y,z);
				}
			}
		}
	}
	
	/*! Save a cartesian field to a fits file.
	 * Arguments in argv:
	 * 0 : IDL_STRING *\a fieldFile -> name of the file 
	 * 1 : IDL_LONG64 *\a Nside     -> number of pixels along one side of the cube
	 * 2 : float      *\a size	-> physical size of the field (in Mpc/h)
	 * 2 : float      *\a field     -> pointer to the array in which the field is stored */
	void MRS3D_saveField_IDL(int argc, void *argv[]){
		char buffer[256];
		sprintf(buffer, "Error in MRS3D_saveField_IDL : Invalid number of arguments %d instead of 4\n",argc);
		if(argc != 4){
			IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, buffer);
		}
		
		IDL_STRING *fieldFile  = (IDL_STRING *) argv[0];
		IDL_LONG64 *Nside	   = (IDL_LONG64 *) argv[1];
		float      *size	   = (float *) argv[2]; 
		float      *field_     = (float *) argv[3];

		
		Field<float> field(*Nside,*size);
		
		for (int64 x=0; x < *Nside; ++x) {
			for (int64 y=0; y < *Nside; ++y) {
				for (int64 z=0; z < *Nside; ++z) {
					field(x,y,z) = field_[z*((*Nside)*(*Nside)) + y*(*Nside) + x];
				}
			}
		}
		write_Field_to_fits(std::string(fieldFile->s), field,  planckType<float>());
	}
	
	/*! Computes almn coefficients from cartesian field
	 *  Arguments in argv:
	 *  0 : IDL_STRING *\a tableFile  -> name of the qln table file
	 *  1 : IDL_STRING *\a fieldFile  -> name of the field fits file
	 *  2 : IDL_STRING *\a AlmnFile   -> name of the almn fits file
	 *  3 : IDL_LONG64 *\a Nside      -> nside parameter of the Healpix maps
	 *  4 : IDL_LONG64 *\a lmax       -> lmax of almn to use
	 *  5 : IDL_LONG64 *\a mmax       -> mmax of almn to use
	 *  6 : IDL_LONG64 *\a nmax       -> nmax of almn to use
	 *  7 : float      *\a Rmax       -> Rmax value to use for the computation of almns */
	void MRS3D_FieldFileToAlmnFile_IDL(int argc, void *argv[]){
		char buffer[256];
		sprintf(buffer, "Error in MRS3D_FieldFileToAlmnFile_IDL : Invalid number of arguments %d instead of 8.\n",argc);
		if(argc != 8){
			IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, buffer);
		}
		
		IDL_STRING *tableFile  = (IDL_STRING *) argv[0];
		IDL_STRING *fieldFile  = (IDL_STRING *) argv[1];
		IDL_STRING *almnFile   = (IDL_STRING *) argv[2];
		IDL_LONG64 *Nside	   = (IDL_LONG64 *) argv[3];
		IDL_LONG64 *lmax	   = (IDL_LONG64 *) argv[4];
		IDL_LONG64 *mmax	   = (IDL_LONG64 *) argv[5];
		IDL_LONG64 *nmax	   = (IDL_LONG64 *) argv[6];
		float      *Rmax	   = (float *)      argv[7];
		
		int FieldNside;
		float size;
		
		get_fieldsize(std::string(fieldFile->s), FieldNside, size);
		
		Field<float> field(FieldNside,size);
		Almn<xcomplex<float> > almn(*lmax,*mmax,*nmax,*Rmax);
		
		read_Field_from_fits(std::string(fieldFile->s), field);
			
		field2almn(std::string(tableFile->s), field,almn,*Nside);
		
		write_Almn_to_fits(std::string(almnFile->s), almn, *lmax, *mmax, *nmax, planckType<float>());

	}
	
	/*! Reconstructs a cartesian density field from almns.
	 *  Arguments in argv:
	 *  0 : IDL_STRING *\a tableFile  -> name of the qln table file
	 *  1 : IDL_STRING *\a fieldFile  -> name of the field fits file
	 *  2 : IDL_STRING *\a AlmnFile   -> name of the almn fits file
	 *  3 : IDL_LONG64 *\a fieldNside -> number of pixels along one side of the cube
	 *  4 : float      *\a fieldSize  -> physical size of the field to reconstruct
	 *  5 : IDL_LONG64 *\a Nside      -> nside parameter of the Healpix maps */
	void MRS3D_AlmnFileToFieldFile_IDL(int argc, void *argv[]){
		char buffer[256];
		sprintf(buffer, "Error in MRS3D_AlmnFileToFieldFile_IDL : Invalid number of arguments %d instead of 6.\n",argc);
		if(argc != 6){
			IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, buffer);
		}
		
		IDL_STRING *tableFile   = (IDL_STRING *) argv[0];
		IDL_STRING *fieldFile   = (IDL_STRING *) argv[1];
		IDL_STRING *almnFile    = (IDL_STRING *) argv[2];
		IDL_LONG64 *fieldNside  = (IDL_LONG64 *) argv[3];
		float      *fieldSize   = (float *) argv[4];
		IDL_LONG64 *Nside	    = (IDL_LONG64 *) argv[5];
		
		int lmax,mmax,nmax;
		get_almnsize(std::string(almnFile->s),lmax, mmax, nmax);
		
		Almn<xcomplex<float> > almn(lmax,mmax,nmax);
		Field<float> field(*fieldNside,*fieldSize);
		
		read_Almn_from_fits(std::string(almnFile->s), almn, lmax, mmax, nmax);
		
		almn2field(std::string(tableFile->s), almn, field, *Nside);

		write_Field_to_fits(std::string(fieldFile->s), field,  planckType<float>());
	}
	
}
/*! \} */