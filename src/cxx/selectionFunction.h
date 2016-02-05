/*
 * selectionFunction.h - This file is part of MRS3D
 * 
 * Copyright 2011 CEA
 * Contributor : François Lanusse (francois.lanusse@gmail.com)
 * 
 * Created on 16/05/11
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


#ifndef SELECTION_FUNCTION_H
#define SELECTION_FUNCTION_H

#include <math.h>

/*! Abstract class for a selection function */
template<typename T> class selectionFunction{
	
public:
	
	virtual T evaluate(T rho){
		return 1;
	}
};

/*! Implementation of a normalized gaussian selection function */
template<typename T> class gaussianSelection : public selectionFunction<T>{
	
	
public:
	gaussianSelection(T rho0) : 
	rho0_(rho0)
	{}
	
	T evaluate(T rho){
		return 4.0/(sqrt(M_PI)*rho0_*rho0_*rho0_) * exp(-pow(rho/rho0_, 2));
	}
	
private:
	T rho0_;
};

#endif