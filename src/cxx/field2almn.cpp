/*
 * Field2Almn.cpp - This file is part of MRS3D
 * 
 * Copyright 2011 CEA
 * Contributor : François Lanusse (francois.lanusse@gmail.com)
 * 
 * Created on 28/06/11
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
#include "almn_field_tools.h"
#include "field_fitsio.h"
#include <mpi.h>

int main(int argc, char* argv[]){
  
	std::cout << std::endl;
	std::cout << "+------------+" << std::endl;
	std::cout << "| field2almn |" << std::endl;
	std::cout << "+------------+" << std::endl;
	std::cout << std::endl;
	std::cout << "Computes the Discrete Spherical Fourier-Bessel Transform from a 3D cartesian field .fits file." << std::endl;
	std::cout << std::endl;
	std::cout << "Usage : field2almn qln_Table_File input_field output_almn Healpix_Nside Rmax Lmax Mmax Nmax" << std::endl;
	std::cout << std::endl;
	
	if(argc != 9){
		cout << "Wrong number of arguments" << endl;
		return 0;
	}	
	
	MPI::Init(argc,argv);
	char *qlnTableFile = argv[1];
	char *fieldFile   = argv[2];
	char *almnFile     = argv[3];
	
	int64  nside = atoi(argv[4]);
	double Rmax  = atof(argv[5]);
	
	int64 lmax	 =  atoi(argv[6]);
	int64 mmax	 =  atoi(argv[7]);
	int64 nmax	 =  atoi(argv[8]);
	
	int NbProc = MPI::COMM_WORLD.Get_size();
	int ProcId = MPI::COMM_WORLD.Get_rank();
	
	int nmax_local = nmax/NbProc;
	int nOffset = ProcId*nmax_local;
	
	
	int FieldNside;
	double size;
	
	get_fieldsize(std::string(fieldFile), FieldNside, size);
	std::cout << "Parameters of the input field: " << std::endl;
	std::cout << "    Nside = " << FieldNside << std::endl;
	std::cout << "    size  = " << size << std::endl;
	std::cout << std::endl;
	
	Field<double> field(FieldNside,size);
	
	Almn<xcomplex<double> > almn(lmax,mmax,nmax,Rmax);
	
	std::cout << "Loading field  ";
	if (MPI::COMM_WORLD.Get_rank() == 0) {
		read_Field_from_fits(std::string(fieldFile), field);
	}
	
	MPI::COMM_WORLD.Bcast(&(field(0,0,0)),field.Nside*field.Nside*field.Nside,MPI::DOUBLE,0);
	
	MPI::COMM_WORLD.Barrier();
	std::cout << " Done" << std::endl;
	std::cout << std::endl;
	
	std::cout << "Computing Discrete Spherical Fourier-Bessel Transform... ";
	field2almn(std::string(qlnTableFile), field,almn,nside,true);
	std::cout << " Done" << std::endl;
	std::cout << std::endl;
	
	if (ProcId == 0) {
		std::cout << "Saving almn... ";
		write_Almn_to_fits(std::string(almnFile), almn, lmax, mmax, nmax, planckType<double>());
		std::cout << " Done" << std::endl;
	}
	
	MPI::COMM_WORLD.Barrier();
	MPI::Finalize();
}	
	
