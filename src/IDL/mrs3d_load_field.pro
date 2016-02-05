;--------------------------------------------------------------------------
; mrs3d_load_Field.pro - This file is part of MRS3D
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
function mrs3d_load_Field,filename, HELP=help
;+
; NAME:
;    mrs3d_load_field
;
; PURPOSE:
;    Loads density field stored in a FITS file.
;
; CATEGORY:
;    IDL Interface to external facility
;
; CALLING SEQUENCE:
;    field = mrs3d_load_field( field_in )
;
; INPUTS:
;   field_in  : Input density field stored in a FITS file
;
; OUTPUTS:
;   field     : array of size [nside,nside,nside] with the loaded field.
;
; KEYWORD PARAMETERS:
;    HELP     = If set, this documentation header is shown
;
; COMMON BLOCKS:
;
;
; SIDE EFFECTS:
;
;
; RESTRICTIONS:
;
;
; PROCEDURE:
;   MRS3D_loadField_IDL
;
; EXAMPLE:
;    mrs3d_field2almn, 'qln.fits', 'field.fits', 'almn.fits', 512, 511, 511, 512, 250.0
;    mrs3d_almn2field, 'qln.fits', 'almn.fits', 'field2.fits', 512, 250.0, 512
;    field  = mrs3d_load_field('field.fits')
;    field2 = mrs3d_load_field('field2.fits')
;    ivolume, field
;    ivolume, field2
;
; will plot the original 3D density field and the reconstructed density field after 
;direct and inverse Spherical Fourier-Bessel Transform.
; 
;-

if keyword_set(help) then begin
    doc_library,'mrs3d_load_field'
    return,-1
endif

; open the file
FXBOPEN, unit, filename, 'XTENSION'

; read and parse the header
header = FXBHEADER(unit)
Nside  = FXPAR(header, 'NSIDE')
size_  = FXPAR(header, 'SIZE')

print,"Field size : ",size_, "(Mpc/h)^3"

field=fltarr(Nside,Nside,Nside)

S = CALL_EXTERNAL('libmrs3d_IDL.dylib','MRS3D_loadField_IDL', string(filename), long64(Nside),float(field), /UNLOAD)

return,field

end
