;--------------------------------------------------------------------------
; mrs3d_almn2cln.pro - This file is part of MRS3D
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
function mrs3d_almn2cln, filename, HELP=help
;+
; NAME:
;    mrs3d_almn2cln
;
; PURPOSE:
;    Computes and returns power spectrum from almn FITS file
;
; CATEGORY:
;    IDL Interface to external facility
;
; CALLING SEQUENCE:
;    cln = mrs3d_almn2cln( almn_in );
;
; INPUTS:
;   almn_in : input almn stored in a FITS file
;
; OUTPUTS:
;    cln : output discretized power spectrum
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
;    MRS3D_AlmnFileToCln_IDL
;
; EXAMPLE:
;    cln = mrs3d_almn2cln('almn.fits')
;    isurface,cln
;
; will plot the 2D power spectrum computed from almn.fits.
; However, the n axis should be normalized using bessel zeros.
;-

if keyword_set(help) then begin
    doc_library,'mrs3d_almn2cln'
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

print,"Cln size : "
print,"lmax : ",nlmax
print,"nmax : ",nnmax
print,"Rmax : ",rmax 

cln = fltArr(nlmax+1,nnmax)

S = CALL_EXTERNAL('libmrs3d_IDL.dylib','MRS3D_AlmnFileToCln_IDL', string(filename), long64(nlmax), long64(nmmax), long64(nnmax),float(cln), /UNLOAD)
return,cln
end