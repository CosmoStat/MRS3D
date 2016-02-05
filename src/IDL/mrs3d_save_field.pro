;--------------------------------------------------------------------------
; mrs3d_save_field.pro - This file is part of MRS3D
; 
; Copyright 2011 CEA
; Contributor : François Lanusse (francois.lanusse@gmail.com)
; 
; Created on 21/06/11
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
pro mrs3d_save_Field,filename,size_,field, HELP=help
;+
; NAME:
;    mrs3d_save_field
;
; PURPOSE:
;    Save density field in a FITS file.
;
; CATEGORY:
;    IDL Interface to external facility
;
; CALLING SEQUENCE:
;    mrs3d_save_field, field_out, field_size, field_in
;
; INPUTS:
;   field_out  : FITS file in which to store the density field
;   field_size : Physical size of one edge of the cube (Mpc/h)
;   cln_in     : Input array of size [nside, nside,nside]
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
;    Creates a new FITS file to store the density field.
;
;
; RESTRICTIONS:
;
;
; PROCEDURE:
;   MRS3D_saveField_IDL
;
; EXAMPLE:
;    field  = mrs3d_load_field('field.fits')
;    field  = field + randomn(1,256,256,256)   ;in the case of 256^3 field
;    mrs3d_save_field, 'field_noisy.fits', 250.0, field
;
; will load a density field, add a gaussian noise and store it back
; in a different file.
;-

if keyword_set(help) then begin
    doc_library,'mrs3d_save_field'
    return
endif

field_size = size(field)
Nside = field_size[1]
print,Nside

S = CALL_EXTERNAL('libmrs3d_IDL.dylib','MRS3D_saveField_IDL', string(filename), long64(Nside),float(size_),float(field), /UNLOAD)

end
