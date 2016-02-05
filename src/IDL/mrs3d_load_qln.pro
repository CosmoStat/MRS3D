;--------------------------------------------------------------------------
; mrs3d_load_qln.pro - This file is part of MRS3D
; 
; Copyright 2011 CEA
; Contributor : François Lanusse (francois.lanusse@gmail.com)
; 
; Created on 21/11/11
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
function mrs3d_load_qln,filename, HELP=help
;+
; NAME:
;    mrs3d_load_qln
;
; PURPOSE:
;    Loads the zeros of Bessel fucntions stored in a FITS table
;
; CATEGORY:
;    IDL Interface to external facility
;
; CALLING SEQUENCE:
;    qln = mrs3d_load_qln( qlnTable )
;
; INPUTS:
;   qlnTable  : Input qln Table stored in a FITS file
;
; OUTPUTS:
;   qln     : array of size [lmax+1,nmax] with the loaded qlns.
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
;  
;
; EXAMPLE:
;    qln = mrs3d_load_qln('qlnTable.fits')
;    plot,qln[12,*]
;
; will load the zeros in qln and then plot the zeros of the Bessel function
; of order 12.
; 
;-

if keyword_set(help) then begin
    doc_library,'mrs3d_load_qln'
    return,-1
endif

; open the file
FXBOPEN,unit,filename,'XTENSION'

; read and parse the header
header = FXBHEADER(unit)
nlmax  = FXPAR(header,'NLMAX')
nnmax  = FXPAR(header,'NNMAX')

; reshape array into a qln(l,n) matrix
qln_array = fltarr(nlmax+1,nnmax)

; read the binary table
FXBREAD,unit, qln_array,1

return,qln_array

END
