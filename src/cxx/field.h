/*
 * field.h - This file is part of MRS3D
 *
 * Copyright 2011 CEA
 * Contributor : François Lanusse (francois.lanusse@gmail.com)
 *
 * Created on 17/06/11
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

#ifndef FIELD_H
#define FIELD_H

#include <arr.h>
#include <paramfile.h>
#include <trafos.h>

#include "selectionFunction.h"

using namespace std;

/*! Class for storing a 3D density field on a cartesian grid. Provides cartesian to spherical coordinates tools.*/
template<typename T> class Field {

public:
    int Nside;    // Number of pixel along on side if the cube
    T size;         // in Mpc/h

    arr3<T> field;

    /*! Creates and allocate a new Field of size \a Nside (number of pixels along one side of the cube)
      and physical size \a size in Mpc/h */
    Field(int Nside_, T size_) :
            Nside(Nside_), size(size_), field(Nside,Nside,Nside)
    {
        //We do not allow totally empty cubes
        if (Nside <= 0) {
            Nside = 1;
        }
        if (size <= 0.0) {
            size = 1.0;
        }

        pixelSize = size/Nside;
        halfSize  = size/2.0;
        Rmax = halfSize;
    }

    /*! Sets all coefficients to zero. */
    void SetToZero() {
        field.fill(0);
    }

    /*! Multiplies whole field by \a factor. */
    template<typename T2> void Scale (const T2 &factor) {
        for (int i =0; i <field.size() ; ++i) {
            *(&field(0,0,0)+i) *= factor;
        }
    }

    /*! Returns a reference to the specified coefficient. */
    T &operator() (int x, int y,int z)
    {
        return field(x,y,z);
    }

    /*! Returns a constant reference to the specified coefficient. */
    const T &operator() (int x, int y,int z) const
    {
        return field(x,y,z);
    }

    /*! Returns the value of the field at the cartesian coordinates (x,y,z) by performing trilinear interpolation */
    const T getValueFromCartCoords(T x, T y, T z) {

        x = (x + halfSize)/pixelSize ;
        y = (y + halfSize)/pixelSize ;
        z = (z + halfSize)/pixelSize ;

        int X = max((int) floor(x),0);
        int Y = max((int) floor(y),0);
        int Z = max((int) floor(z),0);

        if (X >= Nside - 1 || Y >= Nside - 1|| Z >= Nside - 1) {
            return 0;
        }

        T t = (x - floor(x));
        T u = (y - floor(y));
        T v = (z - floor(z));

        return	  field(X,Y,Z)       * (1.0 - t) * (1.0 - u) * (1.0 - v)
                 + field(X,Y,Z+1)     * (1.0 - t) * (1.0 - u) *     v
                 + field(X,Y+1,Z)     * (1.0 - t) *     u     * (1.0 - v)
                 + field(X,Y+1,Z+1)   * (1.0 - t) *     u     *     v
                 + field(X+1,Y,Z)     *     t     * (1.0 - u) * (1.0 - v)
                 + field(X+1,Y,Z+1)   *     t     * (1.0 - u) *     v
                 + field(X+1,Y+1,Z)   *     t     *     u     * (1.0 - v)
                 + field(X+1,Y+1,Z+1) *     t     *     u     *     v    ;

    }

    /*! Returns the value of the field at the spherical coordinates (\a rho,\a theta,\a phi) by performing trilinear interpolation */
    const T getValueFromSphCoords(T rho, T theta, T phi) {
        if (rho > Rmax) {
            return 0.0;
        }

        T x = rho*sin(theta)*cos(phi);
        T y = rho*sin(theta)*sin(phi);
        T z = rho*cos(theta);

        return getValueFromCartCoords(x, y, z);
    }

    /*! Applies a radial selection function to the field. */
    void applySelectionFunction(selectionFunction<T> &function) {
        for (int x=0; x < Nside; ++x) {
            for (int y =0; y < Nside; ++y) {
                for (int z =0; z < Nside; ++z) {
                    field(x,y,z) = field(x,y,z) * function.evaluate(getR(x,y,z));
                }
            }
        }
    }

    /*! Deletes the old field and allocates storage according to
     \a Nside_*/
    void Set (int Nside_, T size_)
    {
        Nside = Nside_;
        size = size_;
        field.alloc(Nside_,Nside_,Nside_);
    }

    /*! Returns the Radial value at point specified in cartesian coordinates (\a x, \a y, \a z). */
    T getR(int x, int y, int z) {
        float centeredX,centeredY,centeredZ;
        centeredX = x - ((T) Nside)/2.0;
        centeredY = y - ((T) Nside)/2.0;
        centeredZ = z - ((T) Nside)/2.0;
        return sqrt(centeredX*centeredX + centeredY*centeredY + centeredZ*centeredZ)*pixelSize;
    }

    /*! Converts the cartesian coordinates (\a x, \a y, \a z) to spherical coordinates (\a rho, \a theta, \a phi)*/
    T getSphericalCoordinates(int x, int y, int z, T &rho, T &theta, T &phi) {

        rho   = getR(x,y,z);
        theta = acos((((T) z)*pixelSize - halfSize)/rho);
        phi   = atan2(((T) y)*pixelSize - halfSize, (((T) x)*pixelSize - halfSize));
    }

private:

    T pixelSize;
    T halfSize;
    T Rmax;

};

#endif // FIELD_H