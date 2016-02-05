/*
 * almn_fitsio.cpp - This file is part of MRS3D
 *
 * Copyright 2011 CEA
 * Contributor : François Lanusse (francois.lanusse@gmail.com)
 *
 * Created on 22/04/11
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

#include "almn_fitsio.h"

#include "almn.h"
#include <fitshandle.h>
#include <share_utils.h>
#include <xcomplex.h>
#include <safe_cast.h>

using namespace std;

bool get_almnsize(fitshandle &inp, int &lmax, int &mmax, int &nmax)
{
    if (inp.key_present("MAX-LPOL") && inp.key_present("MAX-MPOL") && inp.key_present("MAX-NPOL"))
    {
        inp.get_key ("MAX-LPOL",lmax);
        inp.get_key ("MAX-MPOL",mmax);
        inp.get_key ("MAX-NPOL",nmax);
        return(true);
    }
    return(false);
}

bool get_almnsize(const string &filename, int &lmax, int &mmax, int &nmax, int hdunum)
{
    fitshandle inp;
    inp.open (filename);
    inp.goto_hdu(hdunum);
    return(get_almnsize (inp, lmax, mmax, nmax));
}

template<typename T> void read_Almn_from_fits
(fitshandle &inp, Almn<xcomplex<T> >&almns, int lmax, int mmax, int nmax)
{
    T rmax;
    int n_almns = safe_cast<int>(inp.nelems(1));
    arr<int> indexl;
    arr<int> indexm;
    arr<int> indexn;

    arr<double> re, im;

    if (inp.key_present("RMAX")) {
        inp.get_key("RMAX",rmax);
    } else {
        std::cout << "Warning: The almn fits file does not contain an Rmax value, set to 1 Mpc/h by default." << std::endl;
        rmax = 1;
    }

    almns.Set(lmax, mmax,nmax,rmax);
    almns.SetToZero();

    // Read the columns headers
    arr<int> colNum;
    int arrFull = 0;
    colNum.alloc(5);
    colNum.fill(-1);
    for (int i=1; i <= inp.ncols() ; i++) {
        if (inp.colname(i).rfind("l index") != -1) {
            if (colNum[0] == -1) {
                arrFull++;
            }
            colNum[0]=i;
        } else if (inp.colname(i).rfind("m index") != -1) {
            if (colNum[1] == -1) {
                arrFull++;
            }
            colNum[1]=i;
        } else if (inp.colname(i).rfind("n index") != -1) {
            if (colNum[2] == -1) {
                arrFull++;
            }
            colNum[2]=i;
        } else if (inp.colname(i).rfind("real") != -1) {
            if (colNum[3] == -1) {
                arrFull++;
            }
            colNum[3]=i;
        } else if (inp.colname(i).rfind("imag") != -1) {
            if (colNum[4] == -1) {
                arrFull++;
            }
            colNum[4]=i;
        }
    }

    planck_assert(arrFull == 5,"The five required columns were not found in the almn fits file.");

    chunkMaker cm(n_almns,inp.efficientChunkSize(1));
    uint64 offset,ppix;
    while (cm.getNext(offset,ppix))
    {
        indexl.alloc(ppix);
        indexm.alloc(ppix);
        indexn.alloc(ppix);

        re.alloc(ppix);
        im.alloc(ppix);
        inp.read_column(colNum[0],indexl,offset);
        inp.read_column(colNum[1],indexm,offset);
        inp.read_column(colNum[2],indexn,offset);
        inp.read_column(colNum[3],re,offset);
        inp.read_column(colNum[4],im,offset);

        for (tsize i=0; i<ppix; ++i)
        {
            planck_assert(indexm[i]>=0,"negative m encountered");
            planck_assert(indexn[i]>0,"n <= 0 encountered");
            planck_assert(indexl[i]>=indexm[i], "wrong l,m combination");
            if ((indexl[i]<=lmax) && (indexm[i]<=mmax) && (indexn[i] <= nmax))
                almns(indexl[i],indexm[i],indexn[i]-1).Set ((T) re[i],(T) im[i]);
        }
    }
}

template void read_Almn_from_fits (fitshandle &inp, Almn<xcomplex<double> > &almns, int lmax, int mmax, int nmax);
template void read_Almn_from_fits (fitshandle &inp, Almn<xcomplex<float> > &almns, int lmax, int mmax, int nmax);


template<typename T> void read_Almn_from_fits
(const string &filename, Almn<xcomplex<T> >&almns, int lmax, int mmax, int nmax,
 int hdunum)
{
    fitshandle inp;
    inp.open (filename);
    inp.goto_hdu(hdunum);
    read_Almn_from_fits(inp,almns,lmax,mmax,nmax);
}

template void read_Almn_from_fits (const string &filename, Almn<xcomplex<double> > &almns, int lmax, int mmax, int nmax, int hdunum);
template void read_Almn_from_fits (const string &filename, Almn<xcomplex<float> > &almns, int lmax, int mmax, int nmax, int hdunum);


template<typename T> void write_Almn_to_fits
(fitshandle &out, const Almn<xcomplex<T> > &almns, int lmax, int mmax, int nmax,
 PDT datatype)
{
    vector<fitscolumn> cols;
    cols.push_back (fitscolumn("l index","unknown",1,PLANCK_INT32));
    cols.push_back (fitscolumn("m index","unknown",1,PLANCK_INT32));
    cols.push_back (fitscolumn("n index","unknown",1,PLANCK_INT32));
    cols.push_back (fitscolumn("real","unknown",1,datatype));
    cols.push_back (fitscolumn("imag","unknown",1,datatype));

    out.insert_bintab(cols);

    arr<int64> indexl;
    arr<int64> indexm;
    arr<int64> indexn;
    arr<double> re, im;

    int64 l=0, m=0, n=1;

    int64 lm=almns.Lmax(), mm=almns.Mmax(),nm = almns.Nmax();

    int n_almns = nmax*(((mmax+1)*(mmax+2))/2 + (mmax+1)*(lmax-mmax));
    chunkMaker cm(n_almns,out.efficientChunkSize(1));

    uint64 offset,ppix;
    while (cm.getNext(offset,ppix))
    {
        indexl.alloc(ppix);
        indexm.alloc(ppix);
        indexn.alloc(ppix);

        re.alloc(ppix);
        im.alloc(ppix);

        for (tsize i=0; i<ppix; ++i)
        {
            indexl[i] = l;
            indexm[i] = m;
            indexn[i] = n;
            if ((l<=lm) && (m<=mm) && (n <= nm))
            {
                re[i] = almns(l,m,n-1).re;
                im[i] = almns(l,m,n-1).im;
            }
            else
            {
                re[i] = 0;
                im[i] = 0;
            }
            ++m;

            if ((m>l) || (m>mmax)) {
                ++l;
                m=0;
            }
            if ((l>lmax)) {
                ++n;
                l=0;
            }
        }
        out.write_column(1,indexl,offset);
        out.write_column(2,indexm,offset);
        out.write_column(3,indexn,offset);
        out.write_column(4,re,offset);
        out.write_column(5,im,offset);
    }
    out.set_key("MAX-LPOL",lmax,"Maximum L multipole order");
    out.set_key("MAX-MPOL",mmax,"Maximum M multipole order");
    out.set_key("MAX-NPOL",nmax,"Maximum N multipole order");
    out.set_key("RMAX",almns.Rmax(),"Boundary Rmax used for quantization (Mpc/h)");
}

template void write_Almn_to_fits(fitshandle &out, const Almn<xcomplex<double> > &almns, int lmax, int mmax, int nmax, PDT datatype);
template void write_Almn_to_fits(fitshandle &out, const Almn<xcomplex<float> > &almns, int lmax, int mmax, int nmax, PDT datatype);
