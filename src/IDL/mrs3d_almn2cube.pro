;--------------------------------------------------------------------------
; mrs3d_almn2cube.pro - This file is part of MRS3D
; 
; Copyright 2011 CEA
; Contributor : François Lanusse (francois.lanusse@gmail.com)
; 
; Created on 20/05/11
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
function mrs3d_almn2cube,qlnTable,filename,nside,nr,rmin,rstep,HELP=help
;+
; NAME:
;    mrs3d_almn2cube
;
; PURPOSE:
;    Reconstructs the density field from almns on one of the 12 base areas
;    of the Healpix sphere (which gives 2D pictures) for increasing radial 
;    values.
;
; CATEGORY:
;    IDL Interface to external facility
;
; CALLING SEQUENCE:
;    cube = mrs3d_almn2cube( qlnTable, almn_in, nside, nr, rmin, rstep )
;
; INPUTS:
;   qlnTable : FITS file with tabulated values of Bessel zeros  
;   almn_in  : Input almn stored in a FITS file
;   nside    : Resolution of the Healpix maps, and size of the 2D images
;in the output data
;   nr       : Number of layers of reconstructed maps to compute
;   rmin     : Starting value for the radius of the reconstructed maps
;   rstep    : Distance between two consecutive layers
;
; OUTPUTS:
;    cube : output data of size nside*nside*nr. Each cube[*,*,i] is a 2D
;image of size nside*nside.
;
; KEYWORD PARAMETERS:
;     HELP     = If set, this documentation header is shown
;
; COMMON BLOCKS:
;
;
; SIDE EFFECTS:
;
;
; RESTRICTIONS:
;   Only compute the first base area of the Healpix Sphere.
;
; PROCEDURE:
;   MRS3D_ExtractCubeFromAlmn_IDL
;
; EXAMPLE:
;    cube = mrs3d_almn2cube( 'qln.fits', 'almn.fits', 256, 128, 100, 1 )
;    tvscl, cube[*,*,0]
;    ivolume,cube
;
; will plot the density extracted at radius 100 Mpc/h from almn.fits on
; 1/12 of the sphere. Then it will plot a 3D density in a cube (scales 
; respected)
; 
;-

if keyword_set(help) then begin
    doc_library,'mrs3d_almn2cube'
    return, -1
endif

; open the file
FXBOPEN, unit, filename, 'XTENSION'

; read and parse the header
header = FXBHEADER(unit)
nlmax  = FXPAR(header, 'MAX-LPOL')
nmmax  = FXPAR(header, 'MAX-MPOL')
nnmax  = FXPAR(header, 'MAX-NPOL')

cube = fltArr(nside,nside,nr)

S = CALL_EXTERNAL('libmrs3d_IDL.dylib','MRS3D_ExtractCubeFromAlmn_IDL', string(qlnTable), string(filename), long64(nlmax), long64(nmmax), long64(nnmax),long64(nside),long64(nr), float(rmin), float(rstep),float(cube), /UNLOAD)
return,cube
end
