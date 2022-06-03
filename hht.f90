! Hilbert-Huang transform for HEALPix maps

program hht

use healpix_modules

implicit none

!===============================================================
real(DP), allocatable          :: map_in(:,:), map_out(:,:,:)
integer                        :: nside, npix, nmaps, ord, n, imfmax, itermax, i
real(DP)                       :: radius, phi
character(len=80)              :: fin, arg, header(43)
character(len=80), allocatable :: fout(:)
!===============================================================
type it1; real(DP), allocatable :: RES(:), IMF(:); end type
type it2; real(DP), allocatable :: MI(:), MH(:), Emax(:), Emin(:), Emean(:); end type
!===============================================================

! Show help information
if (nArguments() == 1) then
	call getArgument(1, fin)
	if (fin == "-help" .or. fin == "-h") then; call info(); stop; end if
	call fatal_error("Cannot parse argument.")
end if

! Show error message if there are not enough arguments
if (nArguments() <= 5) then; call info(); call fatal_error('Not enough arguments.'); end if

! Reading input parameters
call getArgument(1, fin)                      ! Input file name
call getArgument(2, arg); read(arg,*) radius  ! Radius in arcminutes
call getArgument(3, arg); read(arg,*) imfmax  ! Number of IMFs
call getArgument(4, arg); read(arg,*) itermax ! Number of iterations in the EMD
call getArgument(5, arg); read(arg,*) phi     ! Tension parameter

! Show error message if there are not enough arguments
if (nArguments() < 5+2*imfmax+2) call fatal_error('Not enough arguments.')
if (nArguments() > 5+2*imfmax+2) call fatal_error('Too many arguments.')

! Reading output file names
allocate(fout(2*imfmax+2)); do i = 1, 2*imfmax+2; call getArgument(i+5, fout(i)); end do

! Radius from arcminutes to radians
radius = (radius/60)*(pi/180) 

! Parameters of the FITS file containing the input map
npix = getsize_fits(fin, nmaps=nmaps, nside=nside, ordering=ord); n = nside2npix(nside) - 1

! Allocating arrays
allocate(map_in(0:n,nmaps), map_out(0:n,nmaps,2*imfmax), source=0.0)

! Reading input map
call input_map(fin, map_in, npix, nmaps)
write (*,*) "Map read successfully."

! All subroutines are designed for maps with nested ordering
if (ord == 1) call convert_ring2nest(nside, map_in)

! Empirical mode decomposition
call emd(map_in(:,1), nside, radius, imfmax, itermax, phi, map_out(:,1,:))
write (*,*) "Empirical Mode Decomposition completed successfully."

! Go back to ring ordering if necessary
if (ord == 1) then; do i = 1, 2*imfmax; call convert_nest2ring(nside, map_out(:,:,i)); end do; end if

! Generating output files
call write_minimal_header(header, 'map', nside=nside, order=ord)
do i = 1, 2*imfmax; call output_map(map_out(:,:,i), header, fout(i)); end do

! Generating file containing the sum of IMFs
call output_map(map_out(:,:,2*imfmax)+sum(map_out(:,:,1:imfmax),dim=3), header, fout(2*imfmax+1))

! Generating residual file
call output_map(map_in-map_out(:,:,2*imfmax)-sum(map_out(:,:,1:imfmax),dim=3), header, fout(2*imfmax+2))

write (*,*) "Output files generated successfully."

deallocate(map_in, map_out, fout)

contains

! Empirical Mode Decomposition process
subroutine emd(map_in, nside, radius, imfmax, itermax, phi, map_out)
	integer, intent(in)   :: nside, imfmax, itermax
	real(DP), intent(in)  :: map_in(0:12*nside**2-1), radius, phi
	real(DP), intent(out) :: map_out(0:12*nside**2-1,imfmax)
	integer               :: i, j, n, l, p, q, nlist, nmax, nmin, imax(12*nside**2/8), imin(12*nside**2/8), list(8)
	real(DP)              :: eps, neigh(8), z1(0:12*nside**2-1), z2(0:12*nside**2-1)
	type(it1)             :: EMD1(imfmax)
	type(it2)             :: EMD2(imfmax,itermax)
	
	n   = nside2npix(nside) - 1
	l   = nside2npix(nside) * sin(radius/2.0)**2
	eps = 0.01 * sum(map_in) / (n+1)
	
	! IMF number
	do i = 1, imfmax
		write (*,'("- Computing Intrinsic Mode Function", X, I1, X, "out of", X, I1, ".")') i, imfmax
		
		allocate(EMD1(i)%IMF(0:n), EMD1(i)%RES(0:n), source=0.0)
		
		! Iteration number
		do j = 1, itermax
			write (*,'("-- Iteration in progress:", X, I1, X, "out of", X, I1, ".")') j, itermax
			
			allocate(EMD2(i,j)%MI(0:n), EMD2(i,j)%MH(0:n), source=0.0)
			allocate(EMD2(i,j)%Emax(0:n), EMD2(i,j)%Emin(0:n), EMD2(i,j)%Emean(0:n), source=0.0)
			
			! To find the first IMF, start with the map itself as input signal
			if (i == 1 .and. j == 1) EMD2(i,j)%MI = map_in
			! To find the next IMF, now use the residue as input signal
			if (i /= 1 .and. j == 1) EMD2(i,j)%MI = EMD1(i-1)%RES
			! Repeat the process with the resulting signal from the last step as input signal
			if (j /= 1) EMD2(i,j)%MI = EMD2(i,j-1)%MH
			
			!EMD2(i,j)%Emin = HPX_DBADVAL; EMD2(i,j)%Emax = HPX_DBADVAL
			
			! Local extrema counter
			nmin = 0; nmax = 0
			
			! Finding the positions and values of the local extrema of the input map
			do p = 0, n
				! Pixel indices of a neighborhood around "p"
				call neighbours_nest(nside, p, list, nlist)
				
				! Pixel values of the neighborhood around "p"
				forall (q = 1:nlist) neigh(q) = EMD2(i,j)%MI(list(q))
				
				! Indices of the data for the upper envelope
				if (EMD2(i,j)%MI(p) >= maxval(neigh(1:nlist))) then; nmax = nmax + 1; imax(nmax) = p; end if
				
				! Indices of the data for the lower envelope
				if (EMD2(i,j)%MI(p) <= minval(neigh(1:nlist))) then; nmin = nmin + 1; imin(nmin) = p; end if
				
			end do
			
			! Interpolation by convolution with the inverse Laplacian of a Gaussian beam (not finished)
			!call interpolate(nside, 4*nside, nmax, EMD2(i,j)%MI(imax(1:nmax)), 120.0, imax(1:nmax), EMD2(i,j)%Emax)
			!call interpolate(nside, 3*nside, nmin, EMD2(i,j)%MI(imin(1:nmin)), 120.0, imin(1:nmin), EMD2(i,j)%Emin)
			
			! Interpolation subroutine to compute the smooth envelopes
			!$OMP PARALLEL DO
			do p = 0, n
				call rst_interp(nside, EMD2(i,j)%MI(imax(1:nmax)), nmax, imax(1:nmax), p, radius*i, phi, EMD2(i,j)%Emax(p))
				call rst_interp(nside, EMD2(i,j)%MI(imin(1:nmin)), nmin, imin(1:nmin), p, radius*i, phi, EMD2(i,j)%Emin(p))
			end do
			!$OMP END PARALLEL DO
			
			! Computing the mean of the upper envelope and the lower envelope
			EMD2(i,j)%Emean = (EMD2(i,j)%Emax + EMD2(i,j)%Emin) / 2.0
			
			! Subtracting the envelope mean signal from the input signal
			EMD2(i,j)%MH = EMD2(i,j)%MI - EMD2(i,j)%Emean
		end do
		
		! The IMF is defined as the result from the last iteration over "j"
		EMD1(i)%IMF = EMD2(i,itermax)%MH
		
		! The residue is defined as
		EMD1(i)%RES = EMD2(i,1)%MI - EMD1(i)%IMF
		
	end do
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Minima and maxima, after defining the arrays minima(0:12*nside**2-1) and maxima(0:12*nside**2-1)
	!minima = HPX_DBADVAL; do p = 1, nmin; minima(imin(p)) = map_in(imin(p)); end do
	!maxima = HPX_DBADVAL; do p = 1, nmax; maxima(imax(p)) = map_in(imax(p)); end do
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	do i = 1, imfmax; map_out(:,i) = EMD1(i)%IMF; map_out(:,i+4) = EMD1(i)%RES; end do
	
end subroutine emd

! Local interpolation by regularized spline with tension (Mitášová & Mitáš, 1993)
subroutine rst_interp(nside, lut, next, iext, pix, radius, phi, z)
	integer, intent(in)   :: pix, nside, next, iext(next)
	real(DP), intent(in)  :: lut(next), radius, phi
	real(DP), intent(out) :: z
	real(DP)              :: vpix(3)
	real(DP), allocatable :: A(:,:), B(:), v(:,:)
	integer               :: le(next), nle, i, j, l, nlist, info
	integer, allocatable  :: list(:), lutext(:)
	
	l = nside2npix(nside) * sin(radius/2.0)**2
	
	! If the value at pixel "pix" is already an extremum, then return it and stop this subroutine
	!if (any(iext == pix)) then; z = lut(findloc(iext, value=pix, dim=1)); return; end if
	
	! Pixel indices of a neighborhood around "pix"
	allocate(list(0:3*l/2)); call pix2vec_nest(nside, pix, vpix); call query_disc(nside, vpix, radius, list, nlist, nest=1)
	
	! Local extrema counter
	nle = 0
	
	! If "list(i)" is a local extremum, then store it into the array "le"
	do i = 0, nlist-1; if (any(iext == list(i))) then; nle = nle + 1; le(nle) = list(i); end if; end do
	
	! Indices of "lut" corresponding to the local extrema
	allocate(lutext(nle)); forall (i = 1:nle) lutext(i) = findloc(iext, value=le(i), dim=1)
	
	! Computing the interpolation coefficients stored in "X" by solving "A*X=B"
	allocate(A(nle,nle), B(nle), v(3,nle))
	
	! Vector coordinates of the extrema in "le"
	do i = 1, nle; call pix2vec_nest(nside, iext(lutext(i)), v(:,i)); end do
	
	! Computing "A" and "B"
	A(:,1) = 1.0; B = lut(lutext); forall (i=1:nle, j=2:nle) A(i,j) = R(v(:,i),v(:,j-1),phi) - R(v(:,i),v(:,nle),phi)
	
	! Interpolation coefficients "X=A^(-1)*B"
	call lsolve(nle, A, B)
	
	! Interpolated value at pixel "pix"
	z = B(1); do i = 1, nle-1; z = z + B(i+1)*R(vpix,v(:,i),phi); end do
	z = z - sum(B(2:nle))*R(vpix,v(:,nle),phi)
	
	deallocate(list, lutext, A, B, v)
	
end subroutine rst_interp

! Generate the beam window function in multipole space of the inverse Laplacian of a Gaussian beam parametrized by its FWHM
subroutine il_gaussbeam(fwhm, lmax, wlm)
	integer, intent(in)   :: lmax
	real(DP), intent(in)  :: fwhm ! In arcminutes
	real(DP), intent(out) :: wlm(0:lmax)
	integer               :: l
	real(DP)              :: gaussian_wlm(0:lmax,1)
	
	call gaussbeam(fwhm, lmax, gaussian_wlm)
	
	wlm(0) = 0.0; forall (l=1:lmax) wlm(l) = gaussian_wlm(l,1) / (l*(l+1))
end subroutine il_gaussbeam

! Sample beam function on a specified grid of "x=cos(theta)"
subroutine sample_beam(lmax, bl, n, f, nside)
	integer, intent(in)  :: nside, lmax, n
	real(DP), intent(in) :: bl(0:lmax)
	real(DP)             :: f(n), wlm(0:lmax,1), x(n), P(n,0:2), S(n)
	integer              :: i, l
	
	! Initialize sampling grid
	forall (i=1:n) x(i) = dble(2*i-n-1) / (n-1)
	
	! Initialize beam convolved with delta function and pixel window
	wlm = 1.0; call pixel_window(wlm, nside, "pixel_window_n0064.fits")
	forall (l=0:lmax) wlm(l,:) = (2*l+1) / (4.0*pi) * bl(l) * wlm(l,:)
	
	! Start Legendre recursion
	P(:,0) = 1.0; P(:,1) = x
	S = wlm(0,1)*P(:,0) + wlm(1,1)*P(:,1)
	
	! Sum up "Y_{l0}" harmonic series
	do l = 2, lmax
		P(:,mod(l,3)) = ((2*l-1)*x*P(:,mod(l-1,3)) - (l-1)*P(:,mod(l-2,3))) / l
		S = S + wlm(l,1)*P(:,mod(l,3))
	end do
	
	! Return normalized beam
	f = S
	
end subroutine sample_beam

! Linear interpolation on interval "[a,b]" using a uniformly spaced LUT
pure function linterp1d(n, lut, a, b, x)
	integer, intent(in)  :: n
	real(DP), intent(in) :: lut(n), a, b, x
	integer              :: i
	real(DP)             :: w, linterp1d
	
	! Index in lookup table
	w = (x-a)/(b-a)*(n-1)
	i = min(max(floor(w)+1, 1), n-1)
	w = w - dble(i-1)
	
	! Linear combination of LUT values
	linterp1d = (1.0-w) * lut(i) + w * lut(i+1)
	
end function

! Interpolate scattered values using specified beam and locations
subroutine interpolate(nside, lmax, n, lut, fwhm, pixels, map_out)
	integer, intent(in)            :: nside, lmax, n, pixels(n)
	real(DP), intent(in)           :: lut(n), fwhm
	real(DP), intent(out)          :: map_out(0:12*nside**2-1)
	integer                        :: i, j, l, samples, stat
	real(DP)                       :: avg
	integer, allocatable           :: p(:)
	complex(DPC), allocatable      :: alm(:,:,:)
	real(DP), allocatable          :: bl(:), kernel(:), v(:,:), A(:,:), B(:)
	real, parameter                :: rad2arcmin = 60*180/pi
	
	! Allocate storage
	allocate(bl(0:lmax), alm(1,0:lmax,0:lmax), v(3,n), A(n,n), B(n))
	
	! Process interpolation kernel
	call il_gaussbeam(fwhm, lmax, bl)
	
	samples = 16*nside; allocate(kernel(samples))
	call sample_beam(lmax, bl, samples, kernel, nside)
	
	! Sample positions
	do i = 1, n; call pix2vec_nest(nside, pixels(i), v(:,i)); end do
	
	! Build kernel overlap matrix
	B = lut; avg = sum(B)/n; B = (B - avg) * (12*nside**2) / (4.0*pi)
	forall (i=1:n, j=1:n) A(i,j) = linterp1d(samples, kernel, -1.0, 1.0, sum(v(:,i)*v(:,j)))
	
	! Solve "A*X=B"
	call lsolve(n, A, B)
	
	! Inject kernel amplitudes into an emty map
	map_out = 0.0; map_out(pixels) = B
	
	! Convolve injected pixels with intepolation kernel
	call map2alm(nside, lmax, lmax, map_out, alm)
	forall (l=0:lmax) alm(:,l,0:l) = bl(l) * alm(:,l,0:l)
	alm(:,0,0) = sqrt(4.0*pi) * avg
	call alm2map(nside, lmax, lmax, alm, map_out)
	
	deallocate(alm, bl, kernel, v, A, B)

end subroutine interpolate

! Radial basis function for interpolation
pure function R(v1, v2, phi) 
	real(DP), intent(in) :: v1(3), v2(3), phi
	real(DP)             :: v3(3), r12, R

	if (all(v1 == v2)) then; R = 0.0; return; end if
	
	! Vectorial product 
	v3(1) = v1(2) * v2(3) - v1(3) * v2(2)
	v3(2) = v1(3) * v2(1) - v1(1) * v2(3)
	v3(3) = v1(1) * v2(2) - v1(2) * v2(1)
	
	! Angular distance
	r12 = atan2(norm2(v3), dot_product(v1,v2))
	
	!call angdist(v1, v2, r12)
	
	R = - log((phi*r12/2.0)**2) - expint(1,(phi*r12/2.0)**2) - EULER

end function R

! Exponential integral function "E_n(x)" from Numerical Recipes
pure function expint(n, x)
	integer, intent(in)  :: n
	real(DP), intent(in) :: x
	real(DP)             :: expint
	integer, parameter   :: MAXIT = 100
	real(DP), parameter  :: EPS = epsilon(x), BIG = huge(x)*EPS
	integer              :: i, nm1
	real(DP)             :: a, b, c, d, del, fact, h

	!if (n < 0 .or. x < 0.0 .or. (x <= 0.0 .and. n <= 1)) call fatal_error('EXPINT: Cannot parse arguments.')

	if (n == 0) then; expint = exp(-x)/x; return; end if

	nm1 = n - 1

	if (x == 0.0) then
		expint = 1.0/nm1
	else if (x > 1.0) then
		b = x + n
		c = BIG
		d = 1.0/b
		h = d
		
		do i = 1, MAXIT
			a = -i*(nm1 + i)
			b = b + 2.0
			d = 1.0/(a*d + b)
			c = b + a/c
			del = c*d
			h = h*del
			if (abs(del-1.0) <= EPS) exit
		end do
		
		expint = h*exp(-x)
	else
		if (nm1 /= 0) then
			expint = 1.0/nm1
		else
			expint = -log(x)-EULER
		end if
		
		fact = 1.0
		
		do i = 1, MAXIT
			fact = -fact*x/i
			
			if (i /= nm1) then
				del = -fact/(i-nm1)
			else
				del = fact*(-log(x)-EULER+sum(1.0/arth(1,1,nm1)))
			end if
			
			expint = expint + del
			
			if (abs(del) < abs(expint)*EPS) exit
		end do
	end if

end function expint

! Array of given length containing an arithmetic progression
pure function arth(first, increment, length)
	integer, intent(in) :: first, increment, length
	integer             :: arth(length)
	integer             :: k

	if (n > 0) arth(1) = first
	
	do k = 2, length; arth(k) = arth(k-1) + increment; end do

end function arth

! Solve the system of "n" linear equations in "n" unknowns in the form "A*X=B"
subroutine lsolve(n, A, B)
	integer, intent(in)     :: n
	real(DP), intent(in)    :: A(n,n)
	real(DP), intent(inout) :: B(n)
	integer                 :: pivot(n), stat
	
	! Status indicator ("stat /= 0" indicates an error)
	stat = 0
	
	! LAPACK subroutine to compute the solution to a real system of linear equations
	call dgesv(n, 1, A, n, pivot, B, n, stat)
	
	! Stop the program if necessary
	if(stat /= 0) call fatal_error('Singular matrix in lsolve()')
end subroutine lsolve

subroutine info()
	write (*,*) "Usage: hht IFN RAD IMF NIT PHI OF1 OF2 OF3 ... SUM RES"
	write (*,*) "IFN = Input file name"
	write (*,*) "RAD = Radius in arcminutes"
	write (*,*) "IMF = Number of Intrinsic Mode Functions"
	write (*,*) "NIT = Number of iterations in the Empirical Mode Decomposition"
	write (*,*) "PHI = Tension parameter"
	write (*,*) "OFN = Nth output file name"
	write (*,*) "SUM = Name of the output file containing the sum of the Intrinsic Mode Functions"
	write (*,*) "RES = Residual file name"
end subroutine info

end program
