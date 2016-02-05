/*
 * cln_fitsio.cpp - This file is part of MRS3D
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

#include "cln_fitsio.h"

#include <string>
#include "cln.h"
#include <fitshandle.h>
#include <safe_cast.h>
#include <share_utils.h>

using namespace std;

bool get_clnsize(fitshandle &inp, int &lmax, int &nmax)
{
    if (inp.key_present("MAX-LPOL") && inp.key_present("MAX-NPOL"))
    {
        inp.get_key ("MAX-LPOL",lmax);
        inp.get_key ("MAX-NPOL",nmax);
        return(true);
    }
    return(false);
}

bool get_clnsize(const string &filename, int &lmax, int &nmax, int hdunum)
{
    fitshandle inp;
    inp.open (filename);
    inp.goto_hdu(hdunum);
    return(get_clnsize (inp, lmax, nmax));
}

template<typename T> void read_Cln_from_fits
(fitshandle &inp, Cln<T>&cln, int lmax, int nmax)
{
    T rmax;
    int n_clns = safe_cast<int>(inp.nelems(1));
    arr<int> indexl;
    arr<int> indexn;

    arr<double> val;

    if (inp.key_present("RMAX")) {
        inp.get_key("RMAX",rmax);
    } else {
        std::cout << "Warning: The cln fits file does not contain an Rmax value, set to 1 Mpc/h by default." << std::endl;
        rmax = 1;
    }

    cln.Set(lmax,nmax,rmax);
    cln.SetToZero();

    // Read the columns headers
    arr<int> colNum;
    int arrFull = 0;
    colNum.alloc(3);
    colNum.fill(-1);
    for (int i=1; i <= inp.ncols() ; i++) {
        if (inp.colname(i).rfind("l index") != -1) {
            if (colNum[0] == -1) {
                arrFull++;
            }
            colNum[0]=i;
        } else if (inp.colname(i).rfind("n index") != -1) {
            if (colNum[1] == -1) {
                arrFull++;
            }
            colNum[1]=i;
        } else if (inp.colname(i).rfind("cln value") != -1) {
            if (colNum[2] == -1) {
                arrFull++;
            }
            colNum[2]=i;
        }
    }

    planck_assert(arrFull == 3,"The three required columns were not found in the cln fits file.");

    chunkMaker cm(n_clns,inp.efficientChunkSize(1));
    uint64 offset,ppix;
    while (cm.getNext(offset,ppix))
    {
        indexl.alloc(ppix);
        indexn.alloc(ppix);

        val.alloc(ppix);
        inp.read_column(colNum[0],indexl,offset);
        inp.read_column(colNum[1],indexn,offset);
        inp.read_column(colNum[2],val,offset);

        for (tsize i=0; i<ppix; ++i)
        {
            planck_assert(indexn[i]>0,"n <= 0 encountered");
            planck_assert(indexl[i]>=0,"l < 0 encountered");
            if ((indexl[i]<=lmax) && (indexn[i] <= nmax))
                cln(indexl[i],indexn[i]-1) = val[i];
        }
    }
}

template void read_Cln_from_fits (fitshandle &inp, Cln<double> &cln, int lmax,int nmax);
template void read_Cln_from_fits (fitshandle &inp, Cln<float> &cln, int lmax, int nmax);


template<typename T> void read_Cln_from_fits(const string &filename, Cln<T>&cln, int lmax, int nmax, int hdunum)
{
    fitshandle inp;
    inp.open (filename);
    inp.goto_hdu(hdunum);
    read_Cln_from_fits(inp,cln,lmax,nmax);
}

template void read_Cln_from_fits (const string &filename, Cln<double> &cln, int lmax, int nmax, int hdunum);
template void read_Cln_from_fits (const string &filename, Cln<float> &cln, int lmax, int nmax, int hdunum);


template<typename T> void write_Cln_to_fits(fitshandle &out, const Cln<T> &cln, int lmax, int nmax, PDT datatype)
{
    vector<fitscolumn> cols;
    cols.push_back (fitscolumn("l index","unknown",1,PLANCK_INT32));
    cols.push_back (fitscolumn("n index","unknown",1,PLANCK_INT32));
    cols.push_back (fitscolumn("cln value","unknown",1,datatype));

    out.insert_bintab(cols);

    arr<int64> indexl;
    arr<int64> indexn;
    arr<double> val;

    int64 l=0, n=1;

    int64 lm=cln.Lmax(),nm = cln.Nmax();

    int n_clns = nmax*(lmax +1);
    chunkMaker cm(n_clns,out.efficientChunkSize(1));

    uint64 offset,ppix;
    while (cm.getNext(offset,ppix))
    {
        indexl.alloc(ppix);
        indexn.alloc(ppix);

        val.alloc(ppix);

        for (tsize i=0; i<ppix; ++i)
        {
            indexl[i] = l;
            indexn[i] = n;
            if ((l<=lm) && (n <= nm))
            {
                val[i] = cln(l,n-1);
            }
            else
            {
                val[i] = 0;
            }
            ++l;

            if ((l>lmax)) {
                ++n;
                l=0;
            }
        }
        out.write_column(1,indexl,offset);
        out.write_column(2,indexn,offset);
        out.write_column(3,val,offset);
    }
    out.set_key("MAX-LPOL",lmax,"Maximum L multipole order");
    out.set_key("MAX-NPOL",nmax,"Maximum L multipole order");
    out.set_key("RMAX",cln.Rmax(),"Boundary Rmax used for quantization (Mpc/h)");
}

template void write_Cln_to_fits(fitshandle &out, const Cln<double> &cln, int lmax, int nmax, PDT datatype);
template void write_Cln_to_fits(fitshandle &out, const Cln<float> &cln, int lmax, int nmax, PDT datatype);
