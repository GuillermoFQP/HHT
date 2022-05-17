! Hilbert-Huang transform for HEALPix maps

program hht

use healpix_modules

implicit none

!===============================================================
real(DP), allocatable :: map_in(:,:), map_out_1(:,:), map_out_2(:,:), map_out_3(:,:), map_out_4(:,:)
integer               :: nside, npix, nmaps, ord, n, imfmax, itermax
real(DP)              :: radius, phi
character(len=80)     :: fin, fout1, fout2, fout3, fout4, arg, header(43)
!===============================================================
type it1; real(DP), allocatable :: RES(:), IMF(:); end type
type it2; real(DP), allocatable :: MI(:), MH(:), Emax(:), Emin(:), Emean(:); end type
!===============================================================

! Input parameters
if (nArguments() /= 9) then
	write (*,*) "Usage: hht IFN OF1 OF2 OF3 OF4 RAD IMF NIT PHI"
	write (*,*) "IFN = Input file name"
	write (*,*) "OF1 = Output file name 1"
	write (*,*) "OF1 = Output file name 2"
	write (*,*) "OF3 = Output file name 1"
	write (*,*) "OF4 = Output file name 2"
	write (*,*) "RAD = Radius in arcminutes"
	write (*,*) "IMF = Number of Intrinsic Mode Functions"
	write (*,*) "NIT = Number of iterations in the Empirical Mode Decomposition"
	write (*,*) "PHI = Tension parameter"
	call fatal_error('Not enough arguments')
end if

call getArgument(1, fin)   ! Input file name
call getArgument(2, fout1) ! Output 1 file name
call getArgument(3, fout2) ! Output 2 file name
call getArgument(4, fout3) ! Output 3 file name
call getArgument(5, fout4) ! Output 4 file name
call getArgument(6, arg); read(arg,*) radius  ! Radius in arcminutes
call getArgument(7, arg); read(arg,*) imfmax  ! Number of IMFs
call getArgument(8, arg); read(arg,*) itermax ! Number of iterations in the EMD
call getArgument(9, arg); read(arg,*) phi     ! Tension parameter

radius = (radius/60)*(pi/180) ! Arcminutes to radians

npix = getsize_fits(fin, nmaps=nmaps, nside=nside, ordering=ord)
n    = nside2npix(nside) - 1

! Allocating arrays
allocate(map_in(0:n,nmaps), map_out_1(0:n,nmaps), map_out_2(0:n,nmaps), map_out_3(0:n,nmaps), map_out_4(0:n,nmaps), source=0.0)

! Reading input map
call input_map(fin, map_in, npix, nmaps)
write (*,*) "Map read successfully."

! All subroutines are designed for maps with nested ordering
if (ord == 1) call convert_ring2nest(nside, map_in)

! Empirical mode decomposition
call emd(map_in(:,1), nside, radius, imfmax, itermax, phi, map_out_1(:,1), map_out_2(:,1), map_out_3(:,1), map_out_4(:,1))
write (*,*) "Empirical Mode Decomposition completed successfully."

! Go back to ring ordering if necessary
if (ord == 1) then
	call convert_nest2ring(nside, map_out_1)
	call convert_nest2ring(nside, map_out_2)
	call convert_nest2ring(nside, map_out_3)
	call convert_nest2ring(nside, map_out_4)
end if

! Generating output files
call write_minimal_header(header, 'map', nside=nside, order=ord)
call output_map(map_out_1, header, fout1)
call output_map(map_out_2, header, fout2)
call output_map(map_out_3, header, fout3)
call output_map(map_out_4, header, fout4)
write (*,*) "Output files generated successfully."

deallocate(map_in, map_out_1, map_out_2, map_out_3, map_out_4)

contains

! Empirical Mode Decomposition process
subroutine emd(map_in, nside, radius, imfmax, itermax, phi, map_out_1, map_out_2, map_out_3, map_out_4)
	integer, intent(in)   :: nside, imfmax, itermax
	real(DP), intent(in)  :: map_in(0:12*nside**2-1), radius, phi
	real(DP), intent(out) :: map_out_1(0:12*nside**2-1), map_out_2(0:12*nside**2-1), map_out_3(0:12*nside**2-1), map_out_4(0:12*nside**2-1)
	integer               :: i, j, n, l, p, q, nlist, nmax, nmin, imax(12*nside**2), imin(12*nside**2), list(8)
	real(DP)              :: eps, vec(3), neigh(8)
	type(it1)             :: EMD1(imfmax)
	type(it2)             :: EMD2(imfmax,itermax)
	
	n   = nside2npix(nside) - 1
	l   = nside2npix(nside) * sin(radius/2.0)**2
	eps = 0.01 * sum(map_in) / (n+1)
	
	! IMF number
	do i = 1, imfmax
		allocate(EMD1(i)%IMF(0:n), EMD1(i)%RES(0:n))
		
		! Iteration number
		do j = 1, itermax
			allocate(EMD2(i,j)%MI(0:n), EMD2(i,j)%MH(0:n), EMD2(i,j)%Emax(0:n), EMD2(i,j)%Emin(0:n), EMD2(i,j)%Emean(0:n))
			
			! To find the first IMF, start with the map itself as input signal
			if (i == 1 .and. j == 1) EMD2(i,j)%MI = map_in
			! To find the next IMF, now use the residue as input signal
			if (i /= 1 .and. j == 1) EMD2(i,j)%MI = EMD1(i-1)%RES
			! Repeat the process with the resulting signal from the last step as input signal
			if (j /= 1) EMD2(i,j)%MI = EMD2(i,j-1)%MH
			
			! EMD2(i,j)%Emin = HPX_DBADVAL
			! EMD2(i,j)%Emax = HPX_DBADVAL
			
			! Local extrema counter
			nmin = 0; nmax = 0
			
			! Finding the positions and values of the local extrema of the input map
			do p = 0, n
				! Pixel indices of a neighborhood around "p"
				call neighbours_nest(nside, p, list, nlist)
				
				! Pixel values of the neighborhood around "p"
				forall (q = 1:nlist) neigh(q) = EMD2(i,j)%MI(list(q))
				
				! Indices of the data for the upper envelope
				if (EMD2(i,j)%MI(p) >= maxval(neigh(1:nlist))) then
					nmax = nmax + 1
					imax(nmax) = p
				end if
				
				! Indices of the data for the lower envelope
				if (EMD2(i,j)%MI(p) <= minval(neigh(1:nlist))) then
					nmin = nmin + 1
					imin(nmin) = p
				end if
				
			end do
			
			! Interpolation subroutine to compute the smooth envelopes
			!$OMP PARALLEL DO
			do p = 0, n
				call rst_interpolation(EMD2(i,j)%MI, nside, imax, nmax, p, radius, phi, EMD2(i,j)%Emax(p))
				call rst_interpolation(EMD2(i,j)%MI, nside, imin, nmin, p, radius, phi, EMD2(i,j)%Emin(p))
			end do
			!$OMP END PARALLEL DO
			
			write (*,*) minval(EMD2(i,j)%Emax)
			
			! Computing the mean of the upper envelope and the lower envelope
			EMD2(i,j)%Emean = (EMD2(i,j)%Emax + EMD2(i,j)%Emin) / 2.0
			
			! Subtracting the envelope mean signal from the input signal
			EMD2(i,j)%MH = EMD2(i,j)%MI - EMD2(i,j)%Emean
		end do
		
		! The IMF is defined as the result from the last iteration over "j"
		EMD1(i)%IMF = EMD2(i,itermax)%MH
		
		! The residue is defined as
		EMD1(i)%RES = EMD2(i,1)%MI - EMD1(i)%IMF
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Minima and lower envelope
		map_out_1 = HPX_DBADVAL
		do p = 1, nmin; map_out_1(imin(p)) = map_in(imin(p)); end do
		map_out_2 = EMD2(1,itermax)%Emin
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Maxima and upper envelope
		map_out_3 = HPX_DBADVAL
		do p = 1, nmax; map_out_3(imax(p)) = map_in(imax(p)); end do
		map_out_4 = EMD2(1,itermax)%Emax
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	end do
	
end subroutine emd

! Local interpolation by regularized spline with tension (Mitášová & Mitáš, 1993)
subroutine rst_interpolation(map_in, nside, extlabels, extnum, pix, radius, phi, z)
	integer, intent(in)   :: pix, nside, extnum, extlabels(12*nside**2)
	real(DP), intent(in)  :: map_in(0:12*nside**2-1), radius, phi
	real(DP), intent(out) :: z
	real(DP)              :: vec(3)
	real(DP), allocatable :: A(:,:), B(:), extvec(:,:)
	integer               :: extlist(12*nside**2), extcount, i, j, l, p, nlist, info
	integer, allocatable  :: list(:)
	
	l   = nside2npix(nside) * sin(radius/2.0)**2
	
	allocate(list(0:3*l/2))
	
	! If the value at pixel "pix" is already an extremum, then return it and stop this subroutine
	!if (any(extlabels == pix)) then; z = map_in(pix); return; end if
	
	! Pixel indices of a neighborhood around "pix"
	call pix2vec_nest(nside, pix, vec); call query_disc(nside, vec, radius, list, nlist, nest=1)
	
	! Local extrema counter
	extcount = 0
	
	! If "list(p)" is a local extremum, then store "p" into the array "extlist"
	do p = 0, nlist-1
		if (any(extlabels == list(p))) then
			extcount          = extcount + 1
			extlist(extcount) = list(p)
		end if
	end do
	
	! Write the number of extrema around the input pixel
	!write (*,*) extcount, "extrema around pixel", pix 
	!write (*,*) extlist(1:extcount)
	
	! Computing the interpolation coefficients stored in "X" by solving "A*X=B"
	allocate(A(extcount,extcount), B(extcount), extvec(extcount,3))
	
	! Vector coordinates of the extrema in "extlist"
	do i = 1, extcount; call pix2vec_nest(nside, extlist(i), extvec(i,:)); end do
	
	! Computing "A" and "B"
	do i = 1, extcount
		! Interpolation function evaluated at the local extrema
		B(i) = map_in(extlist(i))
		! Coefficients of the system of linear equations
		A(i,1) = 1.0
		do j = 2, extcount; A(i,j) = R(extvec(i,:),extvec(j-1,:),phi) - R(extvec(i,:),extvec(extcount,:),phi); end do
	end do
	
	! Interpolation coefficients "X=A^(-1)*B"
	call lsolve(extcount, A, B)
	
	! First component of the interpolated value at pixel "pix"
	z = B(1)
	
	! Second component of the interpolated value at pixel "pix"
	do i = 1, extcount-1; z = z + B(i+1)*R(vec,extvec(i,:),phi); end do
	
	z = z - sum(B(2:extcount))*R(vec,extvec(extcount,:),phi)
	
	! Write the difference between the true value and the interpolated value at an extremum
	!if (any(extlabels == pix)) write (*,*) z - map_in(pix)
	
	deallocate(list, A, B, extvec)
	
end subroutine rst_interpolation

! Radial basis function for interpolation
function R(v1, v2, phi) 
	real(DP), intent(in) :: v1(3), v2(3), phi
	real(DP)             :: r12, R

	if (all(v1 == v2)) then; R = 0.0; return; end if
	
	call angdist(v1, v2, r12)
	
	R = - log((phi*r12/2.0)**2) - expint(1,(phi*r12/2.0)**2) - EULER

end function R

! Exponential integral function "E_n(x)" from Numerical Recipes
function expint(n, x)
	integer, intent(in)  :: n
	real(DP), intent(in) :: x
	real(DP)             :: expint
	integer, parameter   :: MAXIT = 100
	real(DP), parameter  :: EPS = epsilon(x), BIG = huge(x)*EPS
	integer              :: i, nm1
	real(DP)             :: a, b, c, d, del, fact, h

	if (n < 0 .or. x < 0.0 .or. (x <= 0.0 .and. n <= 1)) call fatal_error('EXPINT: Cannot parse arguments.')

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
		
		if (i > MAXIT) call fatal_error('EXPINT: Continued fraction failed.')
		
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
		
		if (i > MAXIT) call fatal_error('EXPINT: Series failed.')
	end if

end function expint

! Array of given length containing an arithmetic progression
function arth(first, increment, length)
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

end program

