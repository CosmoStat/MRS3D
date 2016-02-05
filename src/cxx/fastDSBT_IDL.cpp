/*
 * fastDSBT_IDL.cpp - This file is part of FastDSBT
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

#include "fastDSBT.h"
#include "sbtools.h"
#include "idl_export.h"
#include <iostream>

/*! \defgroup dsbt_idl_group IDL interface to FastDSBT */
/*! \{ */

extern "C" {

    static FastDSBT *dsbt = NULL;

    /*! Computes the direct transform.
     * Arguments in argv:
     *  0 : IDL_LONG64 *\a l -> order of the transform
     *  1 : double *\a Rmax  -> value of Rmax
     *  2 : double *\a inputArray  -> pointer to the array of input coefficients (size Nmax)
     *  3 : double *\a outputArray -> pointer to the array which will contain output coefficients (size Nmax)
     *
     * or alternatively:
     *  0 : IDL_LONG64 *\a lj    -> order of the transform
     *  1 : IDL_LONG64 *\a lgrid -> order of the grid of the output coefficients
     *  2 : double *\a Rmax  -> value of Rmax
     *  3 : double *\a inputArray  -> pointer to the array of input coefficients (size Nmax)
     *  4 : double *\a outputArray -> pointer to the array which will contain output coefficients (size Nmax)
     */
    void Compute_DiscreteSBT_IDL(int argc, void * argv[]) {
        char buffer[256];
        sprintf(buffer, "Error in DiscreteSBT_IDL : Invalid number of arguments %d instead of 4\n",argc);

        if (argc == 4) {
            IDL_LONG64 *l = (IDL_LONG64 *)		argv[0];
            double *Rmax  = (double *)		argv[1];
            double *inputArray  = (double *)	argv[2];
            double *outputArray = (double *)	argv[3];

            if (dsbt != NULL) {
                dsbt->computeDirectTransform((int64_t) *l,*Rmax, inputArray,outputArray);
            }
        } else if (argc == 5) {
            IDL_LONG64 *lj = (IDL_LONG64 *)		argv[0];
            IDL_LONG64 *lgrid = (IDL_LONG64 *)	argv[1];
            double *Rmax  = (double *)		argv[2];
            double *inputArray  = (double *)	argv[3];
            double *outputArray = (double *)	argv[4];

            if (dsbt != NULL) {
                dsbt->computeDirectTransform((int64_t) *lj, (int64_t) *lgrid,*Rmax, inputArray,outputArray);
            }
        } else {
            IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, buffer);
        }
    }

    /*! Convert coefficients from order l1 to l2
     * Arguments in argv:
     *  0 : IDL_LONG64 *\a l1    -> order of the input transform
     *  1 : IDL_LONG64 *\a l2    -> order of the output transform
     *  2 : double *\a inputArray  -> pointer to the array of input coefficients (size Nmax)
     *  3 : double *\a outputArray -> pointer to the array which will contain output coefficients (size Nmax)
     */
    void Convert_Transform_IDL(int argc, void * argv[]) {
        char buffer[256];
        sprintf(buffer, "Error in ConvertTransform_IDL : Invalid number of arguments %d instead of 4\n",argc);
        if (argc != 4) {
            IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, buffer);
        }

        IDL_LONG64 *l1 = (IDL_LONG64 *)		argv[0];
        IDL_LONG64 *l2 = (IDL_LONG64 *)		argv[1];
        double *inputArray  = (double *)	argv[2];
        double *outputArray = (double *)	argv[3];

        if (dsbt != NULL) {
            dsbt->convertTransform((int64_t) *l1,(int64_t) *l2,inputArray,outputArray);
        }
    }

    /*! Computes the inverse transform.
     * Arguments in argv:
     *  0 : IDL_LONG64 *\a l -> order of the transform
     *  1 : double *\a Rmax  -> value of Rmax
     *  2 : double *\a inputArray  -> pointer to the array of input coefficients (size Nmax)
     *  3 : double *\a outputArray -> pointer to the array which will contain output coefficients (size Nmax)
     *
     * or alternatively:
     *  0 : IDL_LONG64 *\a lj    -> order of the transform
     *  1 : IDL_LONG64 *\a lgrid -> order of the grid of the output coefficients
     *  2 : double *\a Rmax  -> value of Rmax
     *  3 : double *\a inputArray  -> pointer to the array of input coefficients (size Nmax)
     *  4 : double *\a outputArray -> pointer to the array which will contain output coefficients (size Nmax)
     */
    void Compute_DiscreteISBT_IDL(int argc, void * argv[]) {
        char buffer[256];
        sprintf(buffer, "Error in DiscreteISBT_IDL : Invalid number of arguments %d instead of 4 or 5\n",argc);

        if (argc == 4) {
            IDL_LONG64 *l = (IDL_LONG64 *)		argv[0];
            double *Rmax  = (double *)		argv[1];
            double *inputArray  = (double *)	argv[2];
            double *outputArray = (double *)	argv[3];

            if (dsbt != NULL) {
                dsbt->computeInverseTransform((int64_t) *l,*Rmax, inputArray,outputArray);
            }
        } else if (argc == 5) {
            IDL_LONG64 *lj = (IDL_LONG64 *)		argv[0];
            IDL_LONG64 *lgrid = (IDL_LONG64 *)	argv[1];
            double *Rmax  = (double *)		argv[2];
            double *inputArray  = (double *)	argv[3];
            double *outputArray = (double *)	argv[4];

            if (dsbt != NULL) {
                dsbt->computeInverseTransform((int64_t) *lj, (int64_t) *lgrid,*Rmax, inputArray,outputArray);
            }
        } else {
            IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, buffer);
        }
    }

    /*! Compute the radial grid of order l for a given value of Rmax.
     * Arguments in argv:
     *  0 : IDL_LONG64 *\a l    -> order of grid
     *  1 : double *\a Rmax     -> value of Rmax
     *  2 : double *\a outputArray -> pointer to the array which will contain output coefficients (size Nmax)
     */
    void Compute_RGrid_IDL(int argc, void * argv[]) {
        char buffer[256];
        sprintf(buffer, "Error in RGrid_IDL : Invalid number of arguments %d instead of 3\n",argc);
        if (argc != 3) {
            IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, buffer);
        }

        IDL_LONG64 *l = (IDL_LONG64 *)		argv[0];
        double *Rmax  = (double *)		argv[1];
        double *outputArray  = (double *)	argv[2];

        if (dsbt != NULL) {
            dsbt->computeRGrid((int64_t) *l,*Rmax, outputArray);
        }
    }

    /*! Compute the kln grid of order l for a given value of Rmax.
     * Arguments in argv:
     *  0 : IDL_LONG64 *\a l    -> order of grid
     *  1 : double *\a Rmax     -> value of Rmax
     *  2 : double *\a outputArray -> pointer to the array which will contain output coefficients (size Nmax)
     */
    void Compute_KGrid_IDL(int argc, void * argv[]) {
        char buffer[256];
        sprintf(buffer, "Error in KGrid_IDL : Invalid number of arguments %d instead of 3\n",argc);
        if (argc != 3) {
            IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, buffer);
        }

        IDL_LONG64 *l = (IDL_LONG64 *)		argv[0];
        double *Rmax  = (double *)		argv[1];
        double *outputArray  = (double *)	argv[2];

        if (dsbt != NULL) {
            dsbt->computeKGrid((int64_t) *l,*Rmax, outputArray);
        }
    }

    /*! Computes the n first zeros of bessel function of order l and stores them in the file tableFile.
     * Arguments in argv:
     *  0 : IDL_LONG64 *\a n   -> number of zeros
     *  1 : IDL_LONG64 *\a l   -> order of bessel function
     *  2 : IDL_STRING *\a tableFile   -> Name of the output fits file storing the qln coefficients
     */
    void Compute_BesselZeros_IDL(int argc, void * argv[]) {
        char buffer[256];
        sprintf(buffer, "Error in BesselZeros_IDL : Invalid number of arguments %d instead of 3\n",argc);
        if (argc != 3) {
            IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, buffer);
        }

        IDL_LONG64 *n = (IDL_LONG64 *)		argv[0];
        IDL_LONG64 *l = (IDL_LONG64 *)		argv[1];
        IDL_STRING *tableFile = (IDL_STRING *)	argv[2];

        double *qln_ = (double*) malloc(((*l)+1)*(*n)*sizeof(double));

        BesselRoots(*n,*l,qln_);
        write_qln(tableFile->s, *n, *l, qln_);
        free(qln_);

    }

    /*! Initializes the transform with maximum values for L and N and reads the qln values from tablefile
     * Arguments in argv:
     *  0 : IDL_LONG64 *\a N   -> maximum value of N
     *  1 : IDL_LONG64 *\a L   -> maximum value of L
     *  2 : IDL_STRING *\a tableFile   -> Name of the fits file containing qln coefficients
     */
    void Init_DiscreteSBT_IDL(int argc, void * argv[]) {
        char buffer[256];
        sprintf(buffer, "Error in DiscreteSBT_IDL : Invalid number of arguments %d instead of 3\n",argc);
        if (argc != 3) {
            IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, buffer);
        }

        IDL_LONG64 *N = (IDL_LONG64 *)		argv[0];
        IDL_LONG64 *L = (IDL_LONG64 *)		argv[1];
        IDL_STRING *tableFile = (IDL_STRING *)	argv[2];

        if (dsbt != NULL) {
            delete dsbt;
        }
        dsbt = new FastDSBT(*N, *L, tableFile->s);
    }


    /*! Frees the C++ object */
    void Free_DiscreteSBT_IDL(int argc, void * argv[]) {
        if (dsbt != NULL) {
            delete dsbt;
        }
    }
}

/*! \} */