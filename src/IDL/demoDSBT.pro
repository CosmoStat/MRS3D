;----------------------
; How to use fastDSBT
;----------------------

; Parameters for the demonstration
Rmax=10.0d*!pi
Nmax = 100
Lmax = 10
l =1


; Create a table of zeros of Bessel functions for l <= 10 n < 100. The table if saved in 'qlnTable_100_10.fits'
; dsbt_besselZeros,nmax,lmax,qlnTableFileName
dsbt_besselZeros,Nmax,Lmax,'qlnTable_100_10.fits'

; Initialize the Discrete Spherical Bessel Transform Library using the generated qln table
; dsbt_init,nmax,lmax,qlnTableFileName
dsbt_init,Nmax,Lmax,'qlnTable_100_10.fits'


; Compute the r-grid on which to evaluate our test function
; dsbt_rgrid(nmax, l, Rmax)
r = dsbt_rgrid(Nmax,l,Rmax)

; Now we evaluate our test function on this grid
signal = sin(r)

; Compute the direct DSBT of signal
; dsbt_direct(l, Rmax, f)
signal_t = dsbt_direct(l,Rmax,signal)

; To plot the transformed signal we can retrieve the corresponding k-grid
; dsbt_kgrid(nmax, l, Rmax)
k = dsbt_kgrid(Nmax,l,Rmax)

window,0 & plot,k,signal_t

; Now let us compute the inverse transform
; dsbt_inverse(l, Rmax, f)
signal_2 = dsbt_inverse(l,Rmax,signal_t)

; And we plot the original and reconstructed signals
window,1 & plot,r,signal & oplot,r,signal_2
window,2 & plot,r,abs(signal - signal_2)

; Now that we are done with the DSBT we free the environment
;dsbt_free
dsbt_free

END
