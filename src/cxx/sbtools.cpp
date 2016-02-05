/*
 * sbtools.cpp - This file is part of FastDSBT
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

#include "sbtools.h"

#include <fitshandle.h>
#include <iostream>
#include <fstream>

void BesselRoots(int64_t nmax, int64_t lmax, double* qln){
	__sbtools_MOD_brootfinder(qln,&nmax,&lmax);
}

double sphericalBesselJ(int64_t l, double x){
	double result;
	
	__sbtools_MOD_bjl(&l, &x, &result);
	
	return result;
}

void write_qln(char* tableFileName, int64_t nnmax, int64_t nlmax, double *qln){
	fitshandle handle;
	handle.create(std::string(tableFileName));
	
	//creation des colonnes 
	std::vector<fitscolumn> cols;
	cols.push_back(fitscolumn(std::string("Spherical Bessel Zeros (n,l) (FORTRAN ARRAY INDEXATION)"),"unknown",1,planckType<double>()));
	
	handle.insert_bintab(cols);
	
	//Parameters
	handle.set_key("NLMAX",(int64) nlmax, "");
	handle.set_key("NNMAX",(int64) nnmax, "");
	
	handle.write_column_raw_void(1, qln, planckType<double>(), (nlmax+1)*nnmax, 0);
	
	handle.close();
}

int load_qln(const char* tableFileName, int64_t nnmax, int64_t nlmax, double *qln){
	
	//Check if the fitsFile exists
	std::ifstream fitsFile;
	fitsFile.open(tableFileName);
	
	if(fitsFile.is_open()){
		fitsFile.close();
		fitshandle handle;
		
		handle.open(std::string(tableFileName));
		
		handle.goto_hdu(2);
		
		if(handle.key_present("NLMAX") && handle.key_present("NNMAX")){
			
			int nlmax_table, nnmax_table; 
			handle.get_key("NLMAX",nlmax_table);
			handle.get_key("NNMAX",nnmax_table);
			
			if(nnmax_table >= nnmax && nlmax_table >= nlmax){
				
				for (int64_t n = 0; n < nnmax; ++n) {
					handle.read_column_raw_void(1, &(qln[(nlmax+1)*n]) ,planckType<double>(), (nlmax+1), (nlmax_table+1)*n);
				}
				
				handle.close();
				return 0;
			}
		}
		
		handle.close();
	}
	
	return -1;
}
