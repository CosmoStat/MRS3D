;--------------------------------------------------------------------------
; mrs3d_save_cln.pro - This file is part of MRS3D
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
pro mrs3d_save_cln, filename, rmax, cln, HELP=help
;+
; NAME:
;    mrs3d_save_cln
;
; PURPOSE:
;    Save discretized power spectrum in a FITS file.
;
; CATEGORY:
;    IDL Interface to external facility
;
; CALLING SEQUENCE:
;    mrs3d_save_cln, cln_out, rmax, cln_in
;
; INPUTS:
;   cln_out : FITS file in which to store the cln coefficients
;   rmax    : Value of rmax used for the computation of the cln
;   cln_in  : Input array of size [lmax+1,nmax]
;
; OUTPUTS:
;
;
; KEYWORD PARAMETERS:
;    HELP     = If set, this documentation header is shown
;
; COMMON BLOCKS:
;
;
; SIDE EFFECTS:
;    Creates a new FITS file to store the cln.
;
;
; RESTRICTIONS:
;
;
; PROCEDURE:
;   MRS3D_saveCln_IDL
;
; EXAMPLE:
;    cln = mrs3d_almn2cln('almn.fits')
;    mrs3d_save_cln, 'cln.fits', 250.0, cln
;
; will compute the power spectrum from almns and then store it in a FITS file.
; NOTE : the value of Rmax (here 250.0) is normally displayed by mrs3d_load_almn.
;-

if keyword_set(help) then begin
    doc_library,'mrs3d_save_cln'
    return
endif

S = SIZE(cln)

;Check if two dimensional.  
IF S[0] NE 2 THEN $  
   ;Print error message.  
   MESSAGE, 'Cln a is not two dimensional.'  
  
;Get number of columns and rows.  
nlmax = S[1] & nnmax = S[2]  

S = CALL_EXTERNAL('libmrs3d_IDL.dylib','MRS3D_saveCln_IDL', string(filename), long64(nlmax-1), long64(nnmax),float(rmax),float(cln), /UNLOAD)

end
