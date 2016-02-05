/*
 * almn_cln_tools.cpp - This file is part of MRS3D
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

/*
 * This file is inspired from the Healpix_cxx software package distributed
 * under GNU General Public License.
 *
 * For more information about HEALPix, see http://healpix.jpl.nasa.gov
 *
 */


#include "almn_cln_tools.h"
#include "almn.h"
#include "planck_rng.h"
#include "cln.h"
#include "xcomplex.h"
#include "openmp_support.h"
#include "wigner.h"
#include "lsconstants.h"

using namespace std;
template<typename T1, typename T2> void create_almn(const Cln<T1> &cln, Almn<xcomplex<T2> > &almn, planck_rng &rng)
{
    int lmax = almn.Lmax();
    int mmax = almn.Mmax();
    int nmax = almn.Nmax();
    const double hsqrt2 = 1/sqrt(2.);

    for (int l=0; l<=lmax; ++l)
    {
        for (int n=0; n < nmax; ++n) {

            double rms_tt = sqrt(cln(l,n));
            double zeta1_r = rng.rand_gauss();
            almn(l,0,n) = T2(zeta1_r * rms_tt);
            for (int m=1; m<=min(l,mmax); ++m)
            {
                zeta1_r = rng.rand_gauss()*hsqrt2;
                double zeta1_i = rng.rand_gauss()*hsqrt2;
                almn(l,m,n).Set (T2(zeta1_r*rms_tt), T2(zeta1_i*rms_tt));
            }
        }
    }
}

template void create_almn (const Cln<float> &cln, Almn<xcomplex<float> > &almn, planck_rng &rng);
template void create_almn (const Cln<float> &cln, Almn<xcomplex<double> > &almn, planck_rng &rng);
template void create_almn (const Cln<double> &cln, Almn<xcomplex<float> > &almn, planck_rng &rng);
template void create_almn (const Cln<double> &cln, Almn<xcomplex<double> > &almn, planck_rng &rng);


template<typename T1, typename T2> void extract_cln
(const Almn<xcomplex<T1> > &almn1,
 const Almn<xcomplex<T1> > &almn2,Cln<T2> &cln)
{
    planck_assert (almn1.conformable(almn2), "a_lmn are not conformable");
    cln.Set(almn1.Lmax(),almn1.Nmax(),almn1.Rmax());
    for (int l=0; l<=almn1.Lmax(); ++l)
    {
        for (int n=0; n < almn1.Nmax(); ++n) {
            cln(l,n) = almn1(l,0,n).re*almn2(l,0,n).re;
            int limit = min(l,almn1.Mmax());
            for (int m=1; m<=limit; ++m)
                cln(l,n) += 2.0 * (almn1(l,m,n).re*almn2(l,m,n).re + almn1(l,m,n).im*almn2(l,m,n).im);
            cln(l,n) /= (T2) (2*l+1);
        }
    }
}

template void extract_cln
(const Almn<xcomplex<float> > &almn1,
 const Almn<xcomplex<float> > &almn2, Cln<float> &cln);
template void extract_cln
(const Almn<xcomplex<double> > &almn1,
 const Almn<xcomplex<double> > &almn2, Cln<float> &cln);
template void extract_cln
(const Almn<xcomplex<float> > &almn1,
 const Almn<xcomplex<float> > &almn2, Cln<double> &cln);
template void extract_cln
(const Almn<xcomplex<double> > &almn1,
 const Almn<xcomplex<double> > &almn2, Cln<double> &cln);


template<typename T1, typename T2> void extract_cln
(const Almn<xcomplex<T1> > &almn, Cln<T2> &cln)
{
    extract_cln (almn,almn,cln);
}

template void extract_cln
(const Almn<xcomplex<float> > &almn, Cln<float> &cln);
template void extract_cln
(const Almn<xcomplex<double> > &almn, Cln <float> &cln);
template void extract_cln
(const Almn<xcomplex<float> > &almn, Cln<double> &cln);
template void extract_cln
(const Almn<xcomplex<double> > &almn, Cln <double> &cln);

