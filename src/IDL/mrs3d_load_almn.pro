;--------------------------------------------------------------------------
; mrs3d_load_Almn.pro - This file is part of MRS3D
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
function mrs3d_load_Almn, filename, HELP=help
;+
; NAME:
;    mrs3d_load_Almn
;
; PURPOSE:
;    Loads Spherical Fourier-Bessel coefficients stored in a FITS file.
;
; CATEGORY:
;    IDL Interface to external facility
;
; CALLING SEQUENCE:
;    almn = mrs3d_load_almn( almn_in)
;
; INPUTS:
;   almn_in  : Input almn stored in a FITS file
;
; OUTPUTS:
;   almn     : array of size [lmax+1,mmax+1,nmax] with the loaded almn.
;
; KEYWORD PARAMETERS:
;
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
;   MRS3D_loadAlmn_IDL
;
; EXAMPLE:
;    mrs3d_field2almn, 'qln.fits', 'field.fits', 'almn.fits', 512, 511, 511, 512, 250.0
;    almn = mrs3d_load_almn( 'almn.fits' )
;    isurface, almn[*,0,*]
;
; will plot the almn coefficients for m = 0 as a surface.
; 
;-

if keyword_set(help) then begin
    doc_library,'mrs3d_load_almn'
    return,-1
endif

; open the file
FXBOPEN, unit, filename, 'XTENSION'

; read and parse the header
header = FXBHEADER(unit)
nlmax  = FXPAR(header, 'MAX-LPOL')
nmmax  = FXPAR(header, 'MAX-MPOL')
nnmax  = FXPAR(header, 'MAX-NPOL')
rmax   = FXPAR(header, 'RMAX')

print,"Almn size : "
print,"lmax : ",nlmax
print,"mmax : ",nmmax
print,"nmax : ",nnmax
print,"Rmax : ",rmax 

almn_real = fltARR(nlmax+1,nmmax+1,nnmax)
almn_imag = fltARR(nlmax+1,nmmax+1,nnmax)

S = CALL_EXTERNAL('libmrs3d_IDL.dylib','MRS3D_loadAlmn_IDL', string(filename), nlmax,nmmax,nnmax,float(almn_real),float(almn_imag), /UNLOAD)

almn = complex(almn_real,almn_imag)
return, almn
end
