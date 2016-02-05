;--------------------------------------------------------------------------
; mrs3d_save_Almn.pro - This file is part of MRS3D
; 
; Copyright 2011 CEA
; Contributor : François Lanusse (francois.lanusse@gmail.com)
; 
; Created on 16/05/11
; 
; This software is a computer program whose purpose is to apply mutli-
; resolution signal processing algorithms on spherical 3D data.
; 
; This software is governed by the CeCILL  license under French law and
; abiding by the rules of distribution of free software.  You can  use, 
; modify and/ or redistribute the software under the terms of the CeCILL
; license as circulated by CEA, CNRS and INRIA at the following URL
; "http://www.cecill.info". 
; 
; As a counterpart to the access to the source code and  rights to copy,
; modify and redistribute granted by the license, users are provided only
; with a limited warranty  and the software's author,  the holder of the
; economic rights,  and the successive licensors  have only  limited
; liability. 
;
; In this respect, the user's attention is drawn to the risks associated
; with loading,  using,  modifying and/or developing or reproducing the
; software by the user in light of its specific status of free software,
; that may mean  that it is complicated to manipulate,  and  that  also
; therefore means  that it is reserved for developers  and  experienced
; professionals having in-depth computer knowledge. Users are therefore
; encouraged to load and test the software's suitability as regards their
; requirements in conditions enabling the security of their systems and/or 
; data to be ensured and,  more generally, to use and operate it in the 
; same conditions as regards security. 
; 
; The fact that you are presently reading this means that you have had
; knowledge of the CeCILL license and that you accept its terms.
;
;------------------------------------------------------------------------
pro mrs3d_save_Almn,filename,rmax,almn,HELP=help
;+
; NAME:
;    mrs3d_save_almn
;
; PURPOSE:
;    Save Spherical Fourier-Bessel coefficients in a FITS file.
;
; CATEGORY:
;    IDL Interface to external facility
;
; CALLING SEQUENCE:
;    mrs3d_save_almn, almn_out, rmax, almn_in
;
; INPUTS:
;   almn_out : FITS file in which to store the almn coefficients
;   rmax     : Value of rmax used for the computation of the almn
;   almn_in  : Input complex array of size [lmax+1,mmax+1,nmax]
;
; OUTPUTS:
;
;
; KEYWORD PARAMETERS:
;   HELP     = If set, this documentation header is shown
;
; COMMON BLOCKS:
;
;
; SIDE EFFECTS:
;    Creates a new FITS file to store the computed almn.
;
;
; RESTRICTIONS:
;
;
; PROCEDURE:
;   MRS3D_saveAlmn_IDL
;
; EXAMPLE:
;    almn = mrs3d_load_almn( 'almn.fits' )
;    almn[*,0,0] = 0
;    mrs3d_save_almn, 'almn_modified.fits', 250.0, almn
;
; will load a set of almn coefficients, modify them and store them
;back in a FITS file. NOTE : the value of Rmax (here 250.0) is normally
;displayed by mrs3d_load_almn.
;-

if keyword_set(help) then begin
    doc_library,'mrs3d_save_almn'
    return
endif

almn_real = re(almn)
almn_imag = im(almn)

almn_size = size(almn_real)
nlmax = almn_size[1]-1
nmmax = almn_size[2]-1
nnmax = almn_size[3]

S = CALL_EXTERNAL('libmrs3d_IDL.dylib','MRS3D_saveAlmn_IDL', string(filename), long64(nlmax),long64(nmmax),long64(nnmax),float(rmax),float(almn_real),float(almn_imag), /UNLOAD)

end
