/*
 * Almn2Field.cpp - This file is part of MRS3D
 *
 * Copyright 2011 CEA
 * Contributor : François Lanusse (francois.lanusse@gmail.com)
 *
 * Created on 01/08/11
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
#include "field.h"
#include "field_fitsio.h"
#include "almn_field_tools.h"
#include "selectionFunction.h"

int main(int argc, char *argv[]) {


    std::cout << std::endl;
    std::cout << "+------------+" << std::endl;
    std::cout << "| almn2field |" << std::endl;
    std::cout << "+------------+" << std::endl;
    std::cout << std::endl;
    std::cout << "Reconstructs a cartesian 3D field from an almn .fits file." << std::endl;
    std::cout << "Warning: For data visualization only !"<< std::endl;
    std::cout << std::endl;
    std::cout << "Usage : almn2field qln_Table_File input_almn output_field field_Nside field_size Healpix_Nside" << std::endl;
    std::cout << std::endl;

    if (argc != 7) {
        cout << "Wrong number of arguments" << endl;
        return 0;
    }

    char *tableFile   =  argv[1];
    char *fieldFile   =  argv[2];
    char *almnFile    =  argv[3];
    int  fieldNside  =   atoi(argv[4]);
    float fieldSize   = atof(argv[5]);
    int   Nside	    = atoi(argv[6]);

    int lmax,mmax,nmax;
    get_almnsize(std::string(almnFile),lmax, mmax, nmax);

    std::cout << "Parameters of the input almns: " << std::endl;
    std::cout << "    Lmax = " << lmax << std::endl;
    std::cout << "    Mmax = " << mmax << std::endl;
    std::cout << "    Nmax = " << nmax << std::endl;
    std::cout << std::endl;

    Almn<xcomplex<double> > almn;
    Field<double> field(fieldNside,fieldSize);

    std::cout << "Loading almns  ";
    read_Almn_from_fits(std::string(almnFile), almn, lmax, mmax, nmax);
    std::cout << " Done" << std::endl;
    std::cout << std::endl;

    std::cout << "Computing reconstruction... ";
    almn2field(std::string(tableFile), almn, field, Nside);
    std::cout << " Done" << std::endl;
    std::cout << std::endl;

    std::cout << "Saving field... ";
    write_Field_to_fits(std::string(fieldFile), field,  planckType<double>());
    std::cout << " Done" << std::endl;

    return 0;
}
