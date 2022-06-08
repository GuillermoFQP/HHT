! Hilbert-Huang transform for HEALPix maps

program hht

use healpix_modules

implicit none

!===============================================================
real(DP), allocatable          :: map_in(:,:), map_out(:,:,:)
integer                        :: nside, npix, nmaps, ord, n, imf, nit, ch, i, stiff
character(len=80)              :: fin, arg, header(43)
character(len=80), allocatable :: fout(:)
!===============================================================
type it1; real(DP), allocatable :: RES(:), IMF(:); end type
type it2; real(DP), allocatable :: MI(:), Emax(:), Emin(:); end type
!===============================================================

! Show help information
select case (nArguments())
	case (1)
		call getArgument(1, arg); if (arg == "-h") then; call usage(); stop; end if
		call fatal_error("Cannot parse argument.")
	case (3)
		call getArgument(1, arg); if (arg /= "-lext") call fatal_error("Cannot parse argument.")
		call getArgument(2, fin)
		call getArgument(3, arg); read(arg,*) stiff ! Stiffness
		
		! Output file names
		ch = 4; allocate(fout(ch))
		do i = 1, ch; write (fout(i),'("extrema",I0,".fits")') i; end do
	case (4)
		! Input parameters
		call getArgument(1, fin)                    ! Input file name
		call getArgument(2, arg); read(arg,*) imf   ! Number of IMFs
		call getArgument(3, arg); read(arg,*) nit   ! Number of iterations in the EMD
		call getArgument(4, arg); read(arg,*) stiff ! Stiffness
		
		! Output file names
		ch = 2*imf; allocate(fout(ch))
		do i = 1, imf; write (fout(i),'("imf",I0,".fits")') i; end do
		do i = imf+1, ch; write (fout(i),'("res",I0,".fits")') i; end do
	case default
		call fatal_error("Invalid number of arguments.")
end select

! Parameters of the FITS file containing the input map
npix = getsize_fits(fin, nmaps=nmaps, nside=nside, ordering=ord); n = nside2npix(nside) - 1

! Allocating arrays
allocate(map_in(0:n,nmaps), map_out(0:n,nmaps,ch), source=0.0)

! Reading input map
call input_map(fin, map_in, npix, nmaps)
write (*,*) "Map read successfully."

! All subroutines are designed for maps with nested ordering
if (ord == 1) call convert_ring2nest(nside, map_in)

! Find extrema or perforn empirical mode decomposition
if (nArguments() == 3) then
	call extrema(nside, map_in(:,1), stiff, map_out(:,1,:))
	write (*,*) "Extrema and envelopes obtained successfully."
else
	call emd(nside, map_in(:,1), imf, nit, stiff, map_out(:,1,:))
	write (*,*) "Empirical Mode Decomposition completed successfully."
end if

! Go back to ring ordering if necessary
if (ord == 1) then
	call convert_nest2ring(nside, map_in)
	do i = 1, ch; call convert_nest2ring(nside, map_out(:,:,i)); end do
end if

! Generating output files
call write_minimal_header(header, 'map', nside=nside, order=ord)
do i = 1, ch; call output_map(map_out(:,:,i), header, fout(i)); end do

! Generating file containing the EMD and the residual
if (nArguments() == 4) then
	call output_map(map_out(:,:,2*imf)+sum(map_out(:,:,1:imf),dim=3), header, "emd.fits")
	call output_map(map_in-map_out(:,:,2*imf)-sum(map_out(:,:,1:imf),dim=3), header, "res.fits")
end if
write (*,*) "Output files generated successfully."

deallocate(map_in, map_out, fout)

contains

! Empirical Mode Decomposition process
subroutine emd(nside, map_in, imf, nit, stiff, map_out)
	integer, intent(in)   :: nside, imf, nit, stiff
	real(DP), intent(in)  :: map_in(0:12*nside**2-1)
	real(DP), intent(out) :: map_out(0:12*nside**2-1,2*imf)
	integer               :: i, j, n, nlist, nmax, nmin, imax(12*nside**2/8), imin(12*nside**2/8)
	real(DP)              :: eps
	type(it1)             :: EMD1(imf)
	type(it2)             :: EMD2(imf,nit)
	
	n   = nside2npix(nside) - 1
	eps = 0.01 * sum(map_in) / (n+1)
	
	! IMF number
	do i = 1, imf
		write (*,'("- Computing Intrinsic Mode Function", X, I2, X, "out of", X, I2, ".")') i, imf
		
		allocate(EMD1(i)%IMF(0:n), source=0.0)
		
		! Iteration number
		do j = 1, nit
			write (*,'("-- Iteration in progress:", X, I2, X, "out of", X, I2, ".")') j, nit
			
			allocate(EMD2(i,j)%MI(0:n), EMD2(i,j)%Emax(0:n), EMD2(i,j)%Emin(0:n), source=0.0)
			
			! To find the first IMF, start with the map itself as input signal
			if (i == 1 .and. j == 1) EMD2(i,j)%MI = map_in
			! To find the next IMF, now use the residue as input signal
			if (i /= 1 .and. j == 1) EMD2(i,j)%MI = EMD2(i-1,1)%MI - EMD1(i-1)%IMF
			! Repeat the process with the resulting signal from the last step as input signal
			if (j /= 1) EMD2(i,j)%MI = EMD2(i,j-1)%MI - (EMD2(i,j-1)%Emax + EMD2(i,j-1)%Emin) / 2.0
			
			! Finding the positions and values of the local extrema of the input map
			call local_extrema(nside, EMD2(i,j)%MI, nmax, nmin, imax, imin)
			
			! Compute the smooth envelopes (interpolation by spherical spline with stiffness)
			call ss_interp(nside, 100, stiff, nmax, EMD2(i,j)%MI(imax(1:nmax)), imax(1:nmax), EMD2(i,j)%Emax)
			call ss_interp(nside, 100, stiff, nmin, EMD2(i,j)%MI(imin(1:nmin)), imin(1:nmin), EMD2(i,j)%Emin)
		end do
		
		! The IMF is defined as the result from the last iteration over "j"
		EMD1(i)%IMF = EMD2(i,nit)%MI - (EMD2(i,nit)%Emax + EMD2(i,nit)%Emin) / 2.0
		
	end do
	
	do i = 1, imf; map_out(:,i) = EMD1(i)%IMF; map_out(:,i+imf) = EMD2(i,1)%MI - EMD1(i)%IMF; end do
	
end subroutine emd

! Empirical Mode Decomposition process
subroutine extrema(nside, map_in, stiff, map_out)
	integer, intent(in)   :: nside, stiff
	real(DP), intent(in)  :: map_in(0:12*nside**2-1)
	real(DP), intent(out) :: map_out(0:12*nside**2-1,4)
	integer               :: nmax, nmin, imax(12*nside**2/8), imin(12*nside**2/8)
	
	! Local extrema counter
	nmin = 0; nmax = 0
	
	! Finding the positions and values of the local extrema of the input map
	call local_extrema(nside, map_in, nmax, nmin, imax, imin)
	
	! Minima and maxima
	map_out(:,1) = HPX_DBADVAL; forall (i=1:nmax) map_out(imax(i),2) = map_in(imax(i))
	map_out(:,2) = HPX_DBADVAL; forall (i=1:nmin) map_out(imin(i),1) = map_in(imin(i))
	
	! Compute the smooth envelopes (interpolation by spherical spline with stiffness)
	call ss_interp(nside, 100, stiff, nmax, map_in(imax(1:nmax)), imax(1:nmax), map_out(:,3))
	call ss_interp(nside, 100, stiff, nmin, map_in(imin(1:nmin)), imin(1:nmin), map_out(:,4))
	
end subroutine extrema

! Positions and values of the local extrema of an input map
subroutine local_extrema(nside, map_in, nmax, nmin, imax, imin)
	integer, intent(in)  :: nside
	real(DP), intent(in) :: map_in(0:12*nside**2-1)
	integer, intent(out) :: nmax, nmin, imax(12*nside**2/8), imin(12*nside**2/8)
	integer              :: n, p, q, nlist, list(8)
	real(DP)             :: neigh(8)
	
	n = nside2npix(nside) - 1
	
	! Local extrema counter
	nmin = 0; nmax = 0
	
	do p = 0, n
		! Pixel indices of a neighborhood around "p"
		call neighbours_nest(nside, p, list, nlist)
		
		! Pixel values of the neighborhood around "p"
		forall (q = 1:nlist) neigh(q) = map_in(list(q))
		
		! Find local extrema
		if (map_in(p) >= maxval(neigh(1:nlist))) then; nmax = nmax + 1; imax(nmax) = p; end if
		if (map_in(p) <= minval(neigh(1:nlist))) then; nmin = nmin + 1; imin(nmin) = p; end if
	end do
	
end subroutine local_extrema

! Local interpolation by regularized spline with tension (Mitášová & Mitáš, 1993)
subroutine rst_interp(nside, lut, next, iext, pix, radius, phi, z)
	integer, intent(in)   :: pix, nside, next, iext(next)
	real(DP), intent(in)  :: lut(next), radius, phi
	real(DP), intent(out) :: z
	real(DP)              :: vpix(3)
	real(DP), allocatable :: A(:,:), B(:), v(:,:)
	integer               :: le(next), nle, i, j, l, nlist
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

! Global interpolation by spherical spline (Perrin et al, 1988)
subroutine ss_interp(nside, lmax, stiff, next, lut, iext, map_out)
	integer, intent(in)   :: nside, lmax, stiff, next, iext(next)
	real(DP), intent(in)  :: lut(next)
	real(DP), intent(out) :: map_out(0:12*nside**2-1)
	real(DP)              :: A(next,next), B(next), v(3,next)
	integer               :: i, j, p, n
	real                  :: start, finish
	
	n = nside2npix(nside) - 1
	
	! Vector coordinates of the extrema in "le"
	do i = 1, next; call pix2vec_nest(nside, iext(i), v(:,i)); end do
	
	! Computing "A" and "B"
	A(:,1) = 1.0; B = lut; forall (i=1:next, j=2:next) A(i,j) = G(v(:,i),v(:,j-1),lmax,stiff) - G(v(:,i),v(:,next),lmax,stiff)
	
	! Interpolation coefficients "X=A^(-1)*B"
	call lsolve(next, A, B)
	
	! Interpolated value at pixel "p"
	call cpu_time(start)
	!$OMP PARALLEL DO
	do p = 0, n; call ss_interp_val(nside, p, lmax, stiff, next, A, B, v, map_out); end do
	!$OMP END PARALLEL DO
	call cpu_time(finish)
	
	! Show execution time
	write (*,*) "Time to compute interpolated values (in seconds):", finish - start
	
end subroutine ss_interp

! Interpolated value at pixel "p"
subroutine ss_interp_val(nside, p, lmax, stiff, next, A, B, v, map_out)
	integer, intent(in)   :: nside, p, lmax, stiff, next
	real(DP), intent(in)  :: A(next,next), B(next), v(3,next)
	real(DP), intent(out) :: map_out(0:12*nside**2-1)
	real(DP)              :: vp(3)
	integer               :: i
	
	call pix2vec_nest(nside, p, vp)
	
	map_out(p) = B(1); do i = 1, next-1; map_out(p) = map_out(p) + B(i+1) * G(vp,v(:,i),lmax,stiff); end do
	map_out(p) = map_out(p) - sum(B(2:next)) * G(vp,v(:,next),lmax,stiff)
	
end subroutine ss_interp_val

! Basis function for interpolation by spherical spline
pure function G(v1, v2, lmax, stiff)
	integer, intent(in)   :: lmax, stiff
	real(DP), intent(in)  :: v1(3), v2(3)
	real(DP)              :: v3(3), P(0:lmax), S(1:lmax), theta, G
	integer               :: i
	
	! Vectorial product 
	v3(1) = v1(2) * v2(3) - v1(3) * v2(2)
	v3(2) = v1(3) * v2(1) - v1(1) * v2(3)
	v3(3) = v1(1) * v2(2) - v1(2) * v2(1)
	
	! Angular distance
	theta = atan2(norm2(v3), dot_product(v1,v2))
	
	! Array containing Legrendre polynomials up to order "lmax" evaluated at "cos(theta)"
	P = fleg(cos(theta), lmax+1)
	
	! Summands
	forall (i=1:lmax) S(i) = P(i) * dble(2*i+1) / (i*(i+1))**stiff
	
	! Sum
	G = (1.0/4*pi) * sum(S)
	
end function G

! Radial basis function for interpolation by regularized spline
pure function R(v1, v2, phi) 
	real(DP), intent(in) :: v1(3), v2(3), phi
	real(DP)             :: v3(3), theta, R

	if (all(v1 == v2)) then; R = 0.0; return; end if
	
	! Vectorial product 
	v3(1) = v1(2) * v2(3) - v1(3) * v2(2)
	v3(2) = v1(3) * v2(1) - v1(1) * v2(3)
	v3(3) = v1(1) * v2(2) - v1(2) * v2(1)
	
	! Angular distance
	theta = atan2(norm2(v3), dot_product(v1,v2))
	
	R = - log((phi*theta/2.0)**2) - expint(1,(phi*theta/2.0)**2) - EULER

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

! Fitting routine for an expansion with "nl" Legendre polynomials evaluated at "x" from Numerical Recipes
pure function fleg(x, nl)
	real(DP), intent(in) :: x
	integer, intent(in)  :: nl
	real(DP)             :: fleg(nl), d, f1, f2, twox
	integer              :: j
	
	fleg(1) = 1.0; fleg(2) = x
	
	if (nl > 2) then
		twox = 2.0 * x; f2 = x; d = 1.0
		
		do j = 3, nl
			f1 = d; f2 = f2 + twox; d = d + 1.0
			fleg(j) = (f2 * fleg(j-1) - f1 * fleg(j-2)) / d
		end do
	end if
	
end function fleg

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

subroutine usage()
	write (*,'(/,A)') "Usage:"
	write (*,*) "For Hilbert-Huang transform: hht IFN IMF NIT STF"
	write (*,*) "For finding local extrema: hht -lext IFN STF"
	write (*,*) "IFN = Input file name"
	write (*,*) "IMF = Number of Intrinsic Mode Functions"
	write (*,*) "NIT = Number of iterations in the Empirical Mode Decomposition"
	write (*,*) "STF = Stiffness"
	
end subroutine usage

end program
