! Empirical Mode Decomposition for HEALPix maps 

program hht

use healpix_modules

implicit none

!===============================================================
real(DP), allocatable          :: map_in(:,:), map_out(:,:,:)
integer                        :: nside, npix, nmaps, ord, n, imf, nit, ch, i, stiff
character(len=80)              :: fin, arg, header(43)
character(len=80), allocatable :: fout(:)
integer, parameter             :: nLeg = 50   ! Number of Legendre polynomials for the spherical spline interpolation.
integer, parameter             :: nlemin = 30 ! Number of minimum local extrema for local interpolation
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
		do i = imf+1, ch; write (fout(i),'("res",I0,".fits")') i - imf; end do
	case default
		call fatal_error("Invalid number of arguments.")
end select

! Parameters of the FITS file containing the input map
npix = getsize_fits(fin, nmaps=nmaps, nside=nside, ordering=ord); n = nside2npix(nside) - 1

write (*,'(X, "NSIDE = ", I0)') nside

! Allocating arrays
allocate(map_in(0:n,nmaps), map_out(0:n,nmaps,ch), source=0.0)

! Reading input map
call input_map(fin, map_in, npix, nmaps)
write (*,'(/,X,A)') "Map read successfully."

! All subroutines are designed for maps with nested ordering
if (ord == 1) call convert_ring2nest(nside, map_in)

! Find extrema or perform empirical mode decomposition
write (*,'(X,"Using ",I0," Legendre polynomials for the interpolation with stiffness parameter ",I0,".")') nLeg, stiff
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
	integer               :: p, i, j, n, nlist, nmax, nmin, imax(12*nside**2/8), imin(12*nside**2/8)
	real(DP), allocatable :: inp(:), Emax(:), Emin(:)
	
	n   = nside2npix(nside) - 1
	
	! IMF number
	do i = 1, imf
		write (*,'(/,X, "- Computing Intrinsic Mode Function ", I0, " out of ", I0, ".")') i, imf
		
		! Initialize residue
		if (i == 1) map_out(:,imf+i) = map_in
		if (i /= 1) map_out(:,imf+i) = map_out(:,imf+i-1)
		
		! Iteration number
		do j = 1, nit
			write (*,'(/,X, "-- Iteration in progress: " I0, " out of ", I0, ".")') j, nit
			
			allocate(inp(0:n), Emax(0:n), Emin(0:n), source=0.0)
			
			! Start sifting process
			if (j == 1) inp = map_out(:,imf+i)
			if (j /= 1) inp = map_out(:,i)
			
			! Finding the positions and values of the local extrema of the input map
			call local_extrema(nside, inp, nmax, nmin, imax, imin)
			!if (nmax < ishft(n+1,-4) .or. nmin < ishft(n+1,-4)) then
            !    write (*,'(X, "Npix = ", I0, X, "&", X, "Nmax = ", I0, X, "&", X, "Nmin = ", I0)') n + 1, nmax, nmin
            !    call fatal_error("Too few extrema on the map.")
            !end if

            ! Compute the smooth envelopes (2 ways)
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! Interpolation by spherical spline with stiffness
			!call ss_interp(nside, nLeg, stiff, nmax, inp(imax(1:nmax)), imax(1:nmax), Emax)
			!call ss_interp(nside, nLeg, stiff, nmin, inp(imin(1:nmin)), imin(1:nmin), Emin)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Local interpolation by spherical spline with stiffness
            !$OMP PARALLEL DO
            do p = 0, n; call local_interp(nside, p, nLeg, stiff, inp(imax(1:nmax)), nmax, imax(1:nmax), Emax(p)); end do
            !$OMP END PARALLEL DO
            !$OMP PARALLEL DO
            do p = 0, n; call local_interp(nside, p, nLeg, stiff, inp(imin(1:nmin)), nmin, imin(1:nmin), Emin(p)); end do
            !$OMP END PARALLEL DO
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			write (*,*) "    STOP CRITERIA:", maxval(abs((Emax + Emin) / 2.0))
			
			! Update IMF
			map_out(:,i) = inp - (Emax + Emin) / 2.0
			
			deallocate(inp, Emax, Emin)
		
		end do
		
		! Update residue
		map_out(:,imf+i) = map_out(:,imf+i) - map_out(:,i)
		
	end do
	
end subroutine emd

! Empirical Mode Decomposition process
subroutine extrema(nside, map_in, stiff, map_out)
	integer, intent(in)   :: nside, stiff
	real(DP), intent(in)  :: map_in(0:12*nside**2-1)
	real(DP), intent(out) :: map_out(0:12*nside**2-1,4)
	integer               :: nmax, nmin, imax(12*nside**2/8), imin(12*nside**2/8)
	
	! Finding the positions and values of the local extrema of the input map
	call local_extrema(nside, map_in, nmax, nmin, imax, imin)
	
	! Minima and maxima
	write (*,'(/,X,A)') "- Finding local maxima."
	map_out(:,1) = HPX_DBADVAL; forall (i=1:nmax) map_out(imax(i),1) = map_in(imax(i))
	write (*,*) "- Finding local minima."
	map_out(:,2) = HPX_DBADVAL; forall (i=1:nmin) map_out(imin(i),2) = map_in(imin(i))
	
	! Compute the smooth envelopes (interpolation by spherical spline with stiffness)
	write (*,*) "- Computing upper envelope."
	call ss_interp(nside, nLeg, stiff, nmax, map_in(imax(1:nmax)), imax(1:nmax), map_out(:,3))
	write (*,*) "- Computing lower envelope."
	call ss_interp(nside, nLeg, stiff, nmin, map_in(imin(1:nmin)), imin(1:nmin), map_out(:,4))
	
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

! Local interpolation by spherical spline (Perrin et al, 1988)
subroutine local_interp(nside, pix, lmax, stiff, lut, next, iext, z)
	integer, intent(in)   :: nside, pix, lmax, stiff, next, iext(next)
	real(DP), intent(in)  :: lut(next)
	real(DP), intent(out) :: z
	real(DP)              :: vpix(3), rpix, radius
	real(DP), allocatable :: A(:,:), B(:), v(:,:)
	integer               :: nle, i, j, l, nlist
	integer, allocatable  :: list(:), lutext(:), le(:)
	
    ! Initial number of local extrema, radius of one single pixel and initial radius
    nle = 0; rpix = sqrt(4.0*pi/nside2npix(nside)) * sqrt(2.0)/2; radius = rpix
    
    ! Find the adequate number of local extrema around pixel "pix"
    do while (nle == 0)
        ! List size
        l = nside2npix(nside) * sin(radius/2.0)**2
        
        ! Pixel indices of a neighborhood around pixel "pix"
        allocate(list(0:3*l/2), le(1:3*l/4)); call pix2vec_nest(nside, pix, vpix); call query_disc(nside, vpix, radius, list, nlist, nest=1)
        
        ! If "list(i)" is a local extremum, then store it into the array "le"
        do i = 0, nlist-1; if (any(iext == list(i))) then; nle = nle + 1; le(nle) = list(i); end if; end do
        
        ! Increase radius and reset parameters if there are too few local extrema
        if (nle < nlemin) then; radius = radius + rpix; nle = 0; deallocate(list, le); end if
        
    end do
    
    write (*,'(X,"--- Computing kernel matrix of dimension ",I0,".")') nle
	
	! Indices of "lut" corresponding to the local extrema
	allocate(lutext(nle)); forall (i = 1:nle) lutext(i) = findloc(iext, value=le(i), dim=1)
	
	! Computing the interpolation coefficients stored in "X" by solving "A*X=B"
	allocate(A(nle,nle), B(nle), v(3,nle))
	
	! Vector coordinates of the extrema in "le"
	do i = 1, nle; call pix2vec_nest(nside, iext(lutext(i)), v(:,i)); end do
	
	! Computing "A" and "B"
	A(:,1) = 1.0; B = lut(lutext); forall (i=1:nle, j=2:nle) A(i,j) = G(v(:,i),v(:,j-1),lmax,stiff) - G(v(:,i),v(:,nle),lmax,stiff)
	
	! Interpolation coefficients "X=A^(-1)*B"
	call lsolve(nle, A, B)
	
	! Interpolated value at pixel "pix"
	z = B(1); do i = 1, nle-1; z = z + B(i+1)*G(vpix,v(:,i),lmax,stiff); end do
	z = z - sum(B(2:nle))*G(vpix,v(:,nle),lmax,stiff)
	
	deallocate(list, le, lutext, A, B, v)
	
end subroutine local_interp

! Global interpolation by spherical spline (Perrin et al, 1988)
subroutine ss_interp(nside, lmax, stiff, next, lut, iext, map_out)
	integer, intent(in)   :: nside, lmax, stiff, next, iext(next)
	real(DP), intent(in)  :: lut(next)
	real(DP), intent(out) :: map_out(0:12*nside**2-1)
	real(DP)              :: A(next,next), B(next), v(3,next)
	integer               :: i, j, p, n
	
	n = nside2npix(nside) - 1
	
	! Vector coordinates of the extrema in "le"
	do i = 1, next; call pix2vec_nest(nside, iext(i), v(:,i)); end do
	
	! Computing "A" and "B"
	write (*,'(/,X,"--- Computing kernel matrix of dimension ",I0,".")') next
	A(:,1) = 1.0; B = lut
	!$OMP PARALLEL DO
	do i = 1, next; do j = 2, next; A(i,j) = G(v(:,i),v(:,j-1),lmax,stiff) - G(v(:,i),v(:,next),lmax,stiff); end do; end do
	!$OMP END PARALLEL DO
	
	! Interpolation coefficients "X=A^(-1)*B"
	write (*,'(X,A)') "--- Solving system of equations."; call lsolve(next, A, B)
	
	! Interpolated value at pixel "p"
	write (*,'(X,A)') "--- Interpolation started."
	!$OMP PARALLEL DO
	do p = 0, n; call ss_interp_val(nside, p, lmax, stiff, next, B, v, map_out(p)); end do
	!$OMP END PARALLEL DO
	write (*,'(X,A,/)') "--- Interpolation finished."
	
end subroutine ss_interp

! Interpolated value at pixel "p"
subroutine ss_interp_val(nside, p, lmax, stiff, next, B, v, z)
	integer, intent(in)   :: nside, p, lmax, stiff, next
	real(DP), intent(in)  :: B(next), v(3,next)
	real(DP), intent(out) :: z
	real(DP)              :: vp(3)
	integer               :: i
	
	call pix2vec_nest(nside, p, vp)
	
	z = B(1); do i = 1, next-1; z = z + B(i+1) * G(vp,v(:,i),lmax,stiff); end do
	z = z - sum(B(2:next)) * G(vp,v(:,next),lmax,stiff)
	
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
	if (stat /= 0) call fatal_error('Singular matrix in LSOLVE()')
	
end subroutine lsolve

pure function dilog(x)
	real(DP), intent(in) :: x
	integer              :: i
	real(DP)             :: dilog
	
	! Initialize series
	i = 1; dilog = 0.0
	
	do i = 1, 100
	!do while (x**i / dble(i)**2 >= epsilon(x))
		dilog = dilog + x**i / dble(i)**2
		!i = i + 1
	end do
	
end function dilog

subroutine usage()
	write (*,'(/,X,A,/)') "Usage:"
	write (*,*) "For Hilbert-Huang transform: hht IFN IMF NIT STF"
	write (*,'(X,A,/)') "For finding local extrema: hht -lext IFN STF"
	write (*,*) "IFN = Input file name"
	write (*,*) "IMF = Number of Intrinsic Mode Functions"
	write (*,*) "NIT = Number of iterations in the Empirical Mode Decomposition"
	write (*,'(X,A,/)') "STF = Stiffness"
	
end subroutine usage

end program
