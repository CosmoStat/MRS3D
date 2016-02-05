/*
 * backProjection.cpp - This file is part of MRS3D
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

int main(int argc, char* argv[]) {

    std::cout << std::endl;
    std::cout << "+----------------+" << std::endl;
    std::cout << "| backprojection |" << std::endl;
    std::cout << "+----------------+" << std::endl;
    std::cout << std::endl;
    std::cout << "Reconstruct a density cube from only one non-zero wavelet coefficient in the wavelet decomposition" << std::endl;
    std::cout << "For illustration purposes." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage : backprojection qlnTableFile backProjectionField_fits NSide lmax mmax nmax Rmax Jmax Kc j r theta phi" << std::endl;
    std::cout << std::endl;

    if (argc != 14) {
        cout << "Expecting following arguments: qlnTableFile backProjection_fits NSide lmax mmax nmax Rmax Jmax Kc j r theta phi" << endl;
        return 0;
    }
    MPI::Init(argc,argv);

    char *qlnTableFile = argv[1];
    char *almnFileOut  = argv[2];
    int64 Nside        = atoi(argv[3]);
    int lmax           = atoi(argv[4]);
    int mmax           = atoi(argv[5]);
    int nmax	   = atoi(argv[6]);
    double Rmax        = atof(argv[7]);
    int64  Jmax	   = atoi(argv[8]);
    double Kc  	   = atof(argv[9]);
    int64  j	   = atoi(argv[10]);
    double r	   = atof(argv[11]);
    double theta	   = atof(argv[12]);
    double phi         = atof(argv[13]);



    int NbProc = MPI::COMM_WORLD.Get_size();
    int ProcId = MPI::COMM_WORLD.Get_rank();
    int nmax_local = nmax/NbProc;
    int nOffset = ProcId*nmax_local;

    Almn<xcomplex<double> > almn(lmax,mmax,nmax,Rmax);
    almn.SetToZero();

    waveletDecomposition<double> wd(Nside, Jmax , Kc, almn ,std::string(qlnTableFile),true);

    waveletTreatment<double> *operation;
    operation = new backProjection<double>(r,theta,phi,j);

    almn.SetToZero();
    wd.applyTreatment(Nside, j, operation);

    wd.getWlmn(almn,j);
    std::stringstream filename1;
    filename1 << std::string(almnFileOut) << "_wavelet_" << j <<".fits";
    if (ProcId == 0) {
        write_Almn_to_fits(filename1.str(), almn, lmax, mmax, nmax, planckType<double>());
    }
    MPI::COMM_WORLD.Barrier();

    wd.reconstructionAlmn(almn);

    std::stringstream filename;
    filename << std::string(almnFileOut);

    if (MPI::COMM_WORLD.Get_rank() == 0) {
        write_Almn_to_fits( std::string(almnFileOut), almn, lmax, mmax, nmax,planckType<float>());
    }

    MPI::COMM_WORLD.Barrier();

    MPI::Finalize();

}