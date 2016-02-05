/*
 * almn.h - This file is part of MRS3D
 *
 * Copyright 2011 CEA
 * Contributor : François Lanusse (francois.lanusse@gmail.com)
 *
 * Created on 13/04/11
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

#ifndef ALMN_H
#define ALMN_H


#include <alm.h>
#include <arr.h>

/*! Class for storing Spherical Fourier-Bessel coefficients */
template<typename T> class Almn {

private:
    int lmax,mmax,nmax,tval;
    double rmax;
    arr2<T> almn;

public:
    /*! Returns the number of alm coefficients with maximum
        quantum numbers \a l and \a m. */
    static long Num_Alms (int l, int m)
    {
        planck_assert(m<=l,"mmax must not be larger than lmax");
        return ((m+1)*(m+2))/2 + (m+1)*(l-m);
    }

    /*! Returns \a true, if both objects have the same \b lmax, \b mmax, \b nmax  and \b rmax,
     else  \a false. */
    bool conformable (const Almn &other) const
    {
        return ((lmax==other.lmax) && (mmax==other.mmax) && (nmax==other.nmax) && (rmax == other.rmax));
    }



    /*! Constructs an Almn object with given \a lmax, \a mmax and \a nmax, for a given value of \a Rmax. */
    Almn (int lmax_ = 0, int mmax_=0, int nmax_=0, double Rmax_=1)
            : lmax(lmax_), mmax(mmax_), nmax(nmax_),rmax(Rmax_),tval(2*lmax+1),
            almn(nmax,Num_Alms(lmax,mmax))
    {}


    /*! Deletes the old coefficients and allocates storage according to
      \a lmax , \a mmax and \a nmax. Updates the value of Rmax if needed. */
    void Set (int lmax_, int mmax_, int nmax_,  double Rmax_=-1)
    {
        lmax=lmax_;
        mmax=mmax_;
        nmax=nmax_;
        tval=2*lmax+1;
        if (Rmax_ > 0) {
            rmax = Rmax_;
        }
        almn.alloc(nmax,Num_Alms(lmax,mmax));
    }

    /*! Sets all coefficients to zero. */
    void SetToZero() {
        almn.fill(0);
    }

    /*! Multiplies all coefficients by \a factor. */
    template<typename T2> void Scale (const T2 &factor) {
        for (int i =0; i <almn.size() ; ++i) {
            *(&almn(0,0)+i) *= factor;
        }
    }

    /*! \a a(l,m,n) *= \a factor[l] for all \a l,m,n. */
    template<typename T2> void ScaleL (const arr<T2> &factor)
    {
        planck_assert(factor.size()>lmax, "almn.ScaleL: factor array too short");
        for (int n=0; n<nmax; ++n)
            for (int m=0; m<=mmax; ++m)
                for (int l=m; l<=lmax; ++l)
                    operator()(l,m,n)*=factor[l];
    }

    /*! \a a(l,m,n) *= \a factor[n] for all \a l,m,n. */
    template<typename T2> void ScaleN (const arr<T2> &factor)
    {
        planck_assert(factor.size()>nmax, "almn.ScaleN: factor array too short");
        for (int n=0; n<nmax; ++n) {
            for (int i= 0; i < almn.size2(); ++i)
                *(&almn(n,0)+i)*=factor[n];
        }
    }

    /*! Returns a reference to the specified coefficient. */
    T &operator() (int l, int m,int n)
    {
        return almn(n,((m*(tval-m))>>1) + l);
    }

    /*! Returns a constant reference to the specified coefficient. */
    const T &operator() (int l, int m, int n) const
    {
        return almn(n, ((m*(tval-m))>>1) + l);
    }

    /*! Returns the maximum \a l */
    int Lmax() const {
        return lmax;
    }

    /*! Returns the maximum \a m */
    int Mmax() const {
        return mmax;
    }

    /*! Returns the maximum \a n */
    int Nmax() const {
        return nmax;
    }

    /*! Returns the value of \a Rmax */
    double Rmax() const {
        return rmax;
    }

    /*! Returns the size of the full \a n */
    int Size() const {
        return almn.size();
    }
};


#endif // ALMN_H
