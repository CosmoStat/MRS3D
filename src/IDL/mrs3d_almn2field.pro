;--------------------------------------------------------------------------
; mrs3d_almn2field.pro - This file is part of MRS3D
; 
; Copyright 2011 CEA
; Contributor : François Lanusse (francois.lanusse@gmail.com)
; 
; Created on 24/06/11
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
pro mrs3d_almn2field,tableFilename,almnFilename,fieldFilename,Field_Nside,Field_Size,Healpix_nside,HELP=help
;+
; NAME:
;    mrs3d_almn2field
;
; PURPOSE:
;    Reconstructs a cartesian density field from Spherical Fourier-Bessel
;coefficients.
;
; CATEGORY:
;    IDL Interface to external facility
;
; CALLING SEQUENCE:
;    mrs3d_almn2field, qlnTable, almn_in, field_out, field_nside, field_size, healpix_nside 
;
; INPUTS:
;   qlnTable      : FITS file with tabulated values of Bessel zeros  
;   almn_in       : Input almn stored in a FITS file
;   field_out     : Output field FITS file
;   field_nside   : Number of voxels along one side of the cube
;   field_size    : Physical size of one edge of  the field (Mpc/h)
;   healpix_nside : Resolution of the Healpix maps used for the reconstruction
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
;   Creates a new FITS file to store the reconstructed field.
;
; RESTRICTIONS:
;    The name specified in field_out must be available. The reconstructed field
;is not perfect due to losses in the spherical to cartesian conversion, the quality
;of the reconstruction improves with healpix_nside.
;
; PROCEDURE:
;   MRS3D_AlmnFileToFieldFile_IDL
;
; EXAMPLE:
;    mrs3d_almn2field, 'qln.fits', 'almn.fits', 'field.fits', 512, 250.0, 512
;    field = mrs3d_load_field('field.fits')
;    ivolume, field
;
; will plot the 3D density cube reconstructed from almn.fits.
; 
;-

if keyword_set(help) then begin
    doc_library,'mrs3d_almn2field'
    return
endif

  S = CALL_EXTERNAL('libmrs3d_IDL.dylib','MRS3D_AlmnFileToFieldFile_IDL', string(tablefilename), string(fieldfilename), string(almnfilename), long64(Field_nside),float(Field_Size),long64(Healpix_nside), /UNLOAD)
end
