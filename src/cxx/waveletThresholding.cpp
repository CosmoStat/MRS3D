/*
 * waveletThresholding.cpp - This file is part of MRS3D
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

#include "almn.h"
#include "waveletDecomposition.h"
#include "almn_fitsio.h"
#include "cln.h"
#include "cln_fitsio.h"
#include "almn_cln_tools.h"
#include <healpix_map.h>
#include <healpix_map_fitsio.h>
#include "waveletTreatment.h"
#include <mpi.h>

int main(int argc, char* argv[]){
  
	std::cout << std::endl;
	std::cout << "+---------------------+" << std::endl;
	std::cout << "| waveletThresholding |" << std::endl;
	std::cout << "+---------------------+" << std::endl;
	std::cout << std::endl;
	std::cout << "Applies wavelet hard thresholding to the input almn .fits file." << std::endl;
	std::cout << std::endl;
	std::cout << "Usage : waveletThresholding qln_Table_File input_almn output_almn Nside Jmax Kc k sigmaN" << std::endl;
	std::cout << std::endl;
  
	if(argc != 9){
		cout << "Wrong number of arguments" << endl;
		return 0;	
	}	
	
	MPI::Init(argc,argv);
	
	char *qlnTableFile = argv[1];
	char *almnFileIn   = argv[2];
	char *almnFileOut  = argv[3];
	int64 Nside        = atoi(argv[4]);
	int64 Jmax	       = atoi(argv[5]);
	double Kc		   = atof(argv[6]);
	double k		   = atof(argv[7]);
	double sigmaN      = atof(argv[8]);
	
	int lmax,mmax,nmax;
	get_almnsize(std::string(almnFileIn),lmax, mmax, nmax);
	
	std::cout << "Parameters of the input almns: " << std::endl;
	std::cout << "    Lmax = " << lmax << std::endl;
	std::cout << "    Mmax = " << mmax << std::endl;
	std::cout << "    Nmax = " << nmax << std::endl;
	std::cout << std::endl;
	
	int NbProc = MPI::COMM_WORLD.Get_size();
	int ProcId = MPI::COMM_WORLD.Get_rank();
	int nmax_local = nmax/NbProc;
	int nOffset = ProcId*nmax_local;

	Almn<xcomplex<double> > almn(lmax,mmax,nmax);
	
	std::cout << "Loading almns  ";
	if (MPI::COMM_WORLD.Get_rank() == 0) {
		read_Almn_from_fits(std::string(almnFileIn), almn, lmax, mmax, nmax);
		
		for (int i = 1; i < MPI::COMM_WORLD.Get_size(); ++i) {
			MPI::COMM_WORLD.Send(&(almn(0,0,i*nmax/MPI::COMM_WORLD.Get_size())), almn.Size()*sizeof(xcomplex<double>)/MPI::COMM_WORLD.Get_size(), MPI::BYTE, i, i);
			std::cout << ".";
		}
	}else {
		MPI::COMM_WORLD.Recv(&(almn(0,0,MPI::COMM_WORLD.Get_rank()*nmax/MPI::COMM_WORLD.Get_size())), almn.Size()*sizeof(xcomplex<double>)/MPI::COMM_WORLD.Get_size() , MPI::BYTE, 0,MPI::COMM_WORLD.Get_rank());
	}
	
	MPI::COMM_WORLD.Barrier();
	std::cout << " Done" << std::endl;
	std::cout << std::endl;
	
	std::cout << "Computing Wavelet Transform... ";
	waveletDecomposition<double> wd(Nside, Jmax, Kc, almn ,std::string(qlnTableFile),true);
	std::cout << " Done" << std::endl;
	std::cout << std::endl;
	
	
	almn.SetToZero();
	
	waveletTreatment<double> *operation;
	
	operation = new HardThresholding<double>(k,sigmaN);
	
	std::cout << "Applying Hard Thresholding to wavelet scales... ";
	for (int64 j=0; j < Jmax; ++j) {
		wd.applyTreatment(Nside, j, operation);
		wd.getWlmn(almn,j);
		std::stringstream filename;
		filename << std::string(almnFileOut) << "_" << j << ".fits";
		
		if (MPI::COMM_WORLD.Get_rank() == 0) {
			write_Almn_to_fits(filename.str(), almn, lmax, mmax, nmax,planckType<double>());
		}
	}
	std::cout << " Done" << std::endl;
	std::cout << std::endl;
	
	std::cout << "Computing inverse Wavelet Transform... ";
	wd.reconstructionAlmn(almn);
	std::cout << " Done" << std::endl;
	std::cout << std::endl;
	
	std::stringstream filename;
	filename << std::string(almnFileOut);
	
	if (MPI::COMM_WORLD.Get_rank() == 0) {
	std::cout << "Saving almn... ";
	  write_Almn_to_fits( std::string(almnFileOut), almn, lmax, mmax, nmax,planckType<float>());
	std::cout << " Done" << std::endl;
	  
	}
	
	MPI::COMM_WORLD.Barrier();
	
	MPI::Finalize();
	
}