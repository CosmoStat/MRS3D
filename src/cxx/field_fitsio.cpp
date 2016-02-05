/*
 * field_fitsio.h - This file is part of MRS3D
 * 
 * Copyright 2011 CEA
 * Contributor : François Lanusse (francois.lanusse@gmail.com)
 * 
 * Created on 19/06/11
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

#include "field_fitsio.h"

#include <string>
#include "field.h"
#include <fitshandle.h>
#include <safe_cast.h>
#include <share_utils.h>

using namespace std;

template<typename T> bool get_fieldsize(fitshandle &inp, int &Nside, T &size)
{
	if (inp.key_present("NSIDE") && inp.key_present("SIZE"))
    {
		inp.get_key ("NSIDE",Nside);
		inp.get_key ("SIZE",size);
		return(true);
    }
	return(false);
}

template bool get_fieldsize(fitshandle &inp, int &Nside, float &size);
template bool get_fieldsize(fitshandle &inp, int &Nside, double &size);


template<typename T> bool get_fieldsize(const string &filename, int &Nside, T &size, int hdunum)
{
	fitshandle inp;
	inp.open (filename);
	inp.goto_hdu(hdunum);
	return(get_fieldsize (inp, Nside, size));
}

template bool get_fieldsize(const string &filename, int &Nside, float &size, int hdunum);
template bool get_fieldsize(const string &filename, int &Nside, double &size, int hdunum);

template<typename T> void read_Field_from_fits(fitshandle &inp, Field<T>&field)
{
	int64 n_field = safe_cast<int64>(inp.nelems(1));
	arr<int> indexx;
	arr<int> indexy;
	arr<int> indexz;
	
	arr<double> val;
	
	int Nside = 1;
	T size = 1;
	get_fieldsize(inp, Nside, size);
	
	field.Set(Nside,size);
	field.SetToZero();
	
	chunkMaker cm(n_field,inp.efficientChunkSize(1));
	uint64 offset,ppix;
	while(cm.getNext(offset,ppix))
    {
		indexx.alloc(ppix);
		indexy.alloc(ppix);
		indexz.alloc(ppix);
		
		val.alloc(ppix);
		inp.read_column(1,indexx,offset);
		inp.read_column(2,indexy,offset);
		inp.read_column(3,indexz,offset);
		inp.read_column(4,val,offset);
		
		for (tsize i=0; i<ppix; ++i)
		{
			planck_assert(indexx[i]>=0,"x < 0 encountered");
			planck_assert(indexy[i]>=0,"y < 0 encountered");
			planck_assert(indexz[i]>=0,"z < 0 encountered");
			if ((indexx[i]<Nside) && (indexy[i] < Nside) && (indexz[i] < Nside))
				field(indexx[i],indexy[i],indexz[i]) = val[i];
		}
    }
}

template void read_Field_from_fits (fitshandle &inp, Field<double> &field);
template void read_Field_from_fits (fitshandle &inp, Field<float> &field);


template<typename T> void read_Field_from_fits
(const string &filename, Field<T>&field, int hdunum)
{
	fitshandle inp;
	inp.open (filename);
	inp.goto_hdu(hdunum);
	read_Field_from_fits(inp,field);
}

template void read_Field_from_fits (const string &filename, Field<double> &field, int hdunum);
template void read_Field_from_fits (const string &filename, Field<float> &field,  int hdunum);


template<typename T> void write_Field_to_fits(fitshandle &out, const Field<T> &field, PDT datatype)
{
	vector<fitscolumn> cols;
	cols.push_back (fitscolumn("x index","unknown",1,planckType<int64>()));
	cols.push_back (fitscolumn("y index","unknown",1,planckType<int64>()));
	cols.push_back (fitscolumn("z index","unknown",1,planckType<int64>()));
	cols.push_back (fitscolumn("value","unknown",1,datatype));
	
	out.insert_bintab(cols);
	
	arr<int64> indexx;
	arr<int64> indexy;
	arr<int64> indexz;
	arr<double> val;
	
	int64 Nside=1,x=0,y=0,z=0;
	T size =1;
	
	Nside = field.Nside;
	size  = field.size;
	
	int64 n_clns = Nside*Nside*Nside;
	chunkMaker cm(n_clns,out.efficientChunkSize(1));
	
	uint64 offset,ppix;
	while(cm.getNext(offset,ppix))
    {
		indexx.alloc(ppix);
		indexy.alloc(ppix);
		indexz.alloc(ppix);
		
		val.alloc(ppix);
		
		for (tsize i=0; i<ppix; ++i)
		{
			indexx[i] = x;
			indexy[i] = y;
			indexz[i] = z;
			val[i] = field(x,y,z);
			++z;
			if (z >= Nside) { ++y; z=0; }
			if (y >= Nside) { ++x; y=0; }
		}
		out.write_column(1,indexx,offset);
		out.write_column(2,indexy,offset);
		out.write_column(3,indexz,offset);
		out.write_column(4,val,offset);
    }
	out.set_key("NSIDE",Nside,"Number of pixels along one edge");
	out.set_key("SIZE",size,"size of one edge in Mpc/h");
}

template void write_Field_to_fits(fitshandle &out, const Field<double> &field, PDT datatype);
template void write_Field_to_fits(fitshandle &out, const Field<float>  &field, PDT datatype);
