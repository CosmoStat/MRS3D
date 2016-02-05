/*
 * fastDSBT.h - This file is part of FastDSBT
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


#ifndef FAST_DSBT_H
#define FAST_DSBT_H

#include "sbtools.h"
#include <math.h>
#include <cstdlib>
#include <cstdio>
#include <iostream>

/*! Class for processing the Discrete Spherical Bessel Transform.
 *  Stores in memory the values of the transform matrices to allow
 *  fast evalution of the transform at a given order.*/
class FastDSBT {

public:

    /*! Initializes the transform for given \b maximum values of \a Nmax and \a Lmax.
     * The values of the roots of the bessel functions are loaded from the fits file \a qlnTableFile */
    FastDSBT(int64_t Nmax, int64_t Lmax,const char* qlnTableFile) {
        N=Nmax;
        L=Lmax;

        // If the initialization has already been done we free the arrays and start fresh
        //		if (initialized) {
        //			free(qln);
        //			free(transformMatrix);
        //		}

        // Allocate the qln table
        qln = (double*) malloc((L+1)*N*sizeof(double));

        //try to load the qln table from specified file
        if (load_qln(qlnTableFile,N, L,qln)) {
            //An error occured, either the file does not exist or it does not contain a sufficient number of zeros
            //We have to compute them ourselves (it may take quite some time)
            BesselRoots(N,L, qln);
        }

        //Allocate array to store the matrix array in case it is used more than once
        transformMatrix  = (double*) malloc(N*N*sizeof(double));
        conversionMatrix = (double*) malloc(N*N*sizeof(double));
        matrixComputed = false;
        conversionMatrixComputed = false;
        initialized = true;
    }
    
    /*! Frees the space allocated for the transform matrices and qln value */
    ~FastDSBT() {
        if (initialized) {
            free(qln);
            free(transformMatrix);
        }
        initialized = false;
        matrixComputed = false;
    }
    
    /*! Computes the inverse Discrete Spherical Bessel Transform of order \a l for
     *  the specified \a Rmax value. \a inputArray contains the values of the Spherical 
     *  Bessel Coefficients sampled on a kln grid (which can be obtained using \a computeKGrid)
     *  and \a outputArray will contain the result of the inverse transform i.e. a function
     *  on the corresponding rln grid (which can be obtained using \a computeKGrid).*/
    template<typename T1,typename T2>   void computeInverseTransform(int64_t l, T1 Rmax, T2 *inputArray, T2 *outputArray) {
        double normFactor = 1.0/(Rmax*Rmax*Rmax);
        computeTransform(l,l, normFactor, inputArray, outputArray);
    }

    /*! Same as for \b computeInverseTransform(int64_t l, T1 Rmax, T2 *inputArray, T2 *outputArray)
     * except the result of the inverse transform of order \a l is provided on a radial grid of
     * different order \a lgrid.*/
    template<typename T1,typename T2>   void computeInverseTransform(int64_t lj, int64_t lgrid, T1 Rmax, T2 *inputArray, T2 *outputArray) {
        double normFactor = 1.0/(Rmax*Rmax*Rmax);
        computeTransform(lj, lgrid, normFactor, inputArray, outputArray);
    }

    /*! Computes the direct Discrete Spherical Bessel Transform of order \a l for
     *  the specified \a Rmax value. \a inputArray contains an input function sampled on the rln grid
     *  of order \a l (which can be obtained using \a computeRGrid) and \a outputArray will contain the
     *  result of the direct transform i.e. the Spherical Bessel coefficients sampled on the 
     *  corresponding kln grid (which can be obtained using \a computeKGrid).*/
    template<typename T1,typename T2>   void computeDirectTransform(int64_t l, T1 Rmax, T2 *inputArray, T2 *outputArray) {
        double normFactor = 1.0/pow(Kmax(l,(double) Rmax), 3.0);
        computeTransform(l,l, normFactor, inputArray, outputArray);
    }

    /*! Same as for \b computeDirectTransform(int64_t l, T1 Rmax, T2 *inputArray, T2 *outputArray)
     * except the result of the direct transform of order \a l is provided on a kln grid of
     * different order \a lgrid.*/
    template<typename T1,typename T2>   void computeDirectTransform(int64_t lj, int64_t lgrid, T1 Rmax, T2 *inputArray, T2 *outputArray) {
        double normFactor = 1.0/pow(Kmax(lj,(double) Rmax), 3.0);
        computeTransform(lj,lgrid, normFactor, inputArray, outputArray);
    }

    /*! Returns the maximum value of k accessible at order \a l for the specific \a Rmax
     *  from the qln coefficients loaded in memory at the initialization of the object.*/ 
    template<typename T>   T Kmax(int64_t l, T Rmax) {
        if (! initialized) {
            return 0;
        } else {
            return qln[(L+1)*(N-1) + l]/Rmax;
        }
    }

    /*! Returns the value of the rln grid for \a l and \n using the specified \a Rmax value*/ 
    template<typename T>   T get_Rln(int64_t l, int64_t n, T Rmax) {
        if (! (initialized && n < N && l <= L)) {
            return 0;
        } else {
            return Rmax * qln[(L+1)*n + l]/qln[(L+1)*(N-1) + l];
        }
    }

     /*! Returns the value of the kln grid for \a l and \n using the specified \a Rmax value*/ 
    template<typename T>   T get_Kln(int64_t l, int64_t n, T Rmax) {
        if (! (initialized && n < N && l <= L)) {
            return 0;
        } else {
            return qln[(L+1)*n + l]/Rmax;
        }
    }

    /*! Computes all the values of the radial grid of order \a l for a given value of \a Rmax.
     *  The results are stored in \a RGrid which must point to an allocated array of size Nmax.*/
    template<typename T>   void computeRGrid(int64_t l, T Rmax, T *RGrid) {
        if (! (initialized  && l <= L)) {
            return;
        }

        double Kmax = qln[(L+1)*(N-1) + l]/( (double) Rmax);
        for (int64_t n=0 ; n < N; ++n) {
            RGrid[n] = qln[(L+1)*n + l]/Kmax;
        }
    }

     /*! Computes all the values of the kln grid of order \a l for a given value of \a Rmax.
     *  The results are stored in \a KGrid which must point to an allocated array of size Nmax.*/
    template<typename T>   void computeKGrid(int64_t l, T Rmax, T *KGrid) {
        if (! (initialized  && l <= L)) {
            return;
        }

        for (int64_t n=0 ; n < N; ++n) {
            KGrid[n] = qln[(L+1)*n + l]/Rmax;
        }
    }

    /*! Convert the coefficients computed for a given order \a l1 to the coefficients of order \a l2.
     *  Coefficients (which can be either a sampled function or its spherical bessel coefficients) must
     *  be provided in \a inputArray on the grid of order \a l1 and the converted coefficients will be 
     *  sampled on the gird of order \a l2 and stored in \a outputArray which must point to an 
     *  allocated array of size Nmax.*/
    template<typename T>   void convertTransform(int64_t l1, int64_t l2, T *inputArray, T *outputArray) {

        if (! initialized) {
            std::cout << "FastDSBT is not initialized" << std::endl;
            return;
        }

        if (!(conversionMatrixComputed && l1 == conversionMatrixL1 && l2 == conversionMatrixL2)) {
            computeConversionMatrix(l1,l2);
        }

#pragma omp parallel for
        for (int64_t p=0; p < N; ++p) {
            outputArray[p] = 0.0;

            for (int64_t q=0; q < N; ++q) {
                outputArray[p] += inputArray[q] * conversionMatrix[p*N+q];
            }
        }

    }
    
    /*! Returns a pointer to the 2D array containing the transform matrix of order \a lj and
     * which provides results on the grid of order \a lgrid.
     * This matrix is of size Nmax x Nmax. */ 
    double* getTransformMatrix(int64_t lj, int64_t lgrid) {
        if (!(matrixComputed && lj == matrixLJ && lgrid == matrixLGrid)) {
            computeMatrix(lj,lgrid);
        }
        return transformMatrix;
    }

    /*! Returns a pointer to the 2D array containing the conversion matrix between orders \a l1 and
     * l2.
     * This matrix is of size Nmax x Nmax. */ 
    double* getConversionMatrix(int64_t l1, int64_t l2) {
        if (!(conversionMatrixComputed && l1 == conversionMatrixL1 && l2 == conversionMatrixL2)) {
            computeConversionMatrix(l1,l2);
        }

        return conversionMatrix;

    }

private:

    /*! Generic function used to compute the direct of inverse transform with the same transform matrix.
     *  The results must be rescaled by \a normFactor to obtain the correct transform, which is done in 
     *  \b computeInverseTransform and \b computeDirectTransform.
     *  \a lj specifies the order of the transform
     * 	\a lgrid specifies the order of the grid on which are provided the results
     * 	\a normFactor specifies the normalization factor to use which depends on whether the direct or 
     * inverse transform is computed
     *  \a inputArray points to an input array of size Nmax.
     *  \a outputArray points to an allocated array of size Nmax in which are stored the results.*/
    template<typename T>   void computeTransform(int64_t lj, int64_t lgrid, double normFactor, T *inputArray, T *outputArray) {

        if (! initialized) {
            std::cout << "FastDSBT is not initialized" << std::endl;
            return;
        }

        if (!(matrixComputed && lj == matrixLJ && lgrid == matrixLGrid)) {
            computeMatrix(lj,lgrid);
        }

#pragma omp parallel for
        for (int64_t p=0; p < N; ++p) {
            outputArray[p] = 0.0;

            for (int64_t q=0; q < N; ++q) {
                outputArray[p] += inputArray[q] * transformMatrix[p*N+q];
            }

            outputArray[p] = outputArray[p]*normFactor;
        }
    }

    /*! Computes the conversion matrix between orders \a l1 and \a l2. */
    void computeConversionMatrix(int64_t l1, int64_t l2) {

        if (! initialized) {
            std::cout << "FastDSBT is not initialized" << std::endl;
            return;
        }

        if (!(matrixComputed && l2 == matrixLJ && l2 == matrixLGrid)) {
            computeMatrix(l2,l2);
        }

        double normFactor = sqrt(2.0 * M_PI);
        double tempvalue;

#pragma omp parallel for private(tempvalue)
        for (int64_t q=0; q < N; ++q) {
            for (int64_t n=0; n < N; ++n) {
                tempvalue = sphericalBesselJ(l1,qln[(L+1)*n + l2]*qln[(L+1)*q + l1]/qln[(L+1)*(N-1) + l2]);
                for (int64_t p=0; p < N ; ++p) {
                    if (n==0) {
                        conversionMatrix[p*N +q] =0;
                    }
                    conversionMatrix[p*N+q] += tempvalue * transformMatrix[p*N+n];


                    if (n == N -1) {
                        conversionMatrix[p*N + q] = conversionMatrix[p*N + q]  * normFactor/(pow(qln[(L+1)*(N-1) + l2], 3.0) * pow(sphericalBesselJ(l1+1, qln[(L+1)*q + l1]),2.0));
                    }
                }
            }
        }

        conversionMatrixL1 = l1;
        conversionMatrixL2 = l2;
        conversionMatrixComputed = true;
    }

    /*! Computes the transform matrix of order \a lj which provides results on the grid of order \a lgrid. */
    void computeMatrix(int64_t lj, int64_t lgrid) {

        if (! initialized)
            return;

        double normFactor = sqrt(2.0 * M_PI);
        double *weights   = (double *) malloc(N*sizeof(double));
        double tempValue;


#pragma omp parallel private(tempValue)
        {

            if (lj == lgrid) {
#pragma omp for
                for (int64_t p=0; p < N ; ++p) {
                    for (int64_t q=0; q <= p ; ++q) {
                        tempValue = sphericalBesselJ(lj, qln[(L+1)*p + lj]*qln[(L+1)*q + lj]/qln[(L+1)*(N-1) + lj]);

                        transformMatrix[p*N + q] = tempValue;
                        if ( p != q)
                            transformMatrix[q*N + p] = tempValue;
                    }

                    tempValue  = sphericalBesselJ(lj+1, qln[(L+1)*p + lj]);
                    weights[p] = normFactor/(tempValue*tempValue);
                }
            } else {

#pragma omp for
                for (int64_t p=0; p < N ; ++p) {
                    for (int64_t q=0; q < N ; ++q) {
                        tempValue = sphericalBesselJ(lj, qln[(L+1)*p + lgrid]*qln[(L+1)*q + lj]/qln[(L+1)*(N-1) + lgrid]);
                        transformMatrix[p*N + q] = tempValue;
                    }

                    tempValue  = sphericalBesselJ(lj+1, qln[(L+1)*p + lj]);
                    weights[p] = normFactor/(tempValue*tempValue);
                }

            }


#pragma omp for
            for (int64_t p=0; p < N ; ++p) {
                for (int64_t q=0; q < N ; ++q) {
                    transformMatrix[p*N + q] = transformMatrix[p*N + q]*weights[q];
                }
            }

        }

        matrixLJ = lj;
        matrixLGrid = lgrid;
        matrixComputed = true;

        free(weights);
    }

    bool matrixComputed;
    bool conversionMatrixComputed;
    bool initialized;
    int64_t N,L;
    int64_t matrixLJ;
    int64_t matrixLGrid;
    int64_t conversionMatrixL1;
    int64_t conversionMatrixL2;
    double *qln;
    double *transformMatrix;
    double *conversionMatrix;
};

#endif
