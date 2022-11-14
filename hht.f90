! Empirical Mode Decomposition for HEALPix maps 

program hht

use healpix_modules

implicit none

!===============================================================
real(DP), allocatable          :: map_in(:,:), map_out(:,:,:)
integer                        :: nside, npix, nmaps, ord, n, imf, nit, ch, i, stiff, tens, nLeg, fwhm
character(len=80)              :: fin, arg, header(43)
character(len=80), allocatable :: fout(:)
integer, parameter             :: nlemin = 100 ! Number of minimum local extrema for local interpolation
integer, parameter             :: gridres = 64 ! Coarse grid resolution for local interpolation
!===============================================================

! Show help information
select case (nArguments())
	case (1)
		call getArgument(1, arg); if (arg == "-h") then; call usage(); stop; end if
		call fatal_error("Cannot parse argument.")
	case (5)
		call getArgument(1, arg); if (arg /= "-lext") call fatal_error("Cannot parse argument.")
		call getArgument(2, fin)
		call getArgument(3, arg); read(arg,*) stiff ! Stiffness
        call getArgument(4, arg); read(arg,*) tens  ! Tension parameter
        call getArgument(5, arg); read(arg,*) fwhm  ! FWHM for Gaussian smoothing
		
		! Output file names
		ch = 4; allocate(fout(ch))
		do i = 1, ch; write (fout(i),'("extrema",I0,".fits")') i; end do
	case (6)
		! Input parameters
		call getArgument(1, fin)                    ! Input file name
		call getArgument(2, arg); read(arg,*) imf   ! Number of IMFs
		call getArgument(3, arg); read(arg,*) nit   ! Number of iterations in the EMD
		call getArgument(4, arg); read(arg,*) stiff ! Stiffness
        call getArgument(5, arg); read(arg,*) tens  ! Tension parameter
        call getArgument(6, arg); read(arg,*) fwhm  ! FWHM for Gaussian smoothing
		
		! Output file names
		ch = 2*imf; allocate(fout(ch))
		do i = 1, imf; write (fout(i),'("imf",I0,".fits")') i; end do
		do i = imf+1, ch; write (fout(i),'("res",I0,".fits")') i - imf; end do
	case default
		call fatal_error("Invalid number of arguments.")
end select

! Parameters of the FITS file containing the input map
npix = getsize_fits(fin, nmaps=nmaps, nside=nside, ordering=ord); n = nside2npix(nside) - 1; nLeg = 2 * nside

write (*,'(/,X, "Input map Nside = ", I0)') nside

! Allocating arrays
allocate(map_in(0:n,nmaps), map_out(0:n,nmaps,ch), source=0.0)

! Reading input map
call input_map(fin, map_in, npix, nmaps)
write (*,'(/,X,A)') "Map read successfully."

! Gaussian smoothing the input map
if (fwhm /= 0) call smoothing(nside, ord, map_in(:,1), fwhm)

! All the following subroutines are designed for maps in NESTED ordering
if (ord == 1) call convert_ring2nest(nside, map_in)

! Find extrema or perform empirical mode decomposition
write (*,'(/,X,"Using ",I0," Legendre polynomials for the interpolation with stiffness parameter ",I0," and tension parameter ",I0,".")') nLeg, stiff, tens
if (nArguments() == 5) then
	call extrema(nside, map_in(:,1), stiff, tens, map_out(:,1,:))
	write (*,*) "Extrema and envelopes obtained successfully."
else
	call emd(nside, map_in(:,1), imf, nit, stiff, tens, map_out(:,1,:))
	write (*,*) "Empirical Mode Decomposition completed successfully."
end if

! Go back to RING ordering if necessary
if (ord == 1) then
	call convert_nest2ring(nside, map_in)
	do i = 1, ch; call convert_nest2ring(nside, map_out(:,:,i)); end do
end if

! Generating output files
call write_minimal_header(header, 'map', nside=nside, order=ord)
do i = 1, ch; call output_map(map_out(:,:,i), header, fout(i)); end do

! Generating file containing the EMD and the residual
if (nArguments() == 6) then
	call output_map(map_out(:,:,2*imf)+sum(map_out(:,:,1:imf),dim=3), header, "emd.fits")
	call output_map(map_in-map_out(:,:,2*imf)-sum(map_out(:,:,1:imf),dim=3), header, "res.fits")
end if
write (*,*) "Output files generated successfully."

deallocate(map_in, map_out, fout)

contains

! Gaussian smoothing process
subroutine smoothing(nside, ord, map_in, fwhm)
	integer, intent(in)       :: nside, ord, fwhm
	real(DP), intent(inout)   :: map_in(0:12*nside**2-1)
	complex(DPC), allocatable :: alm(:,:,:)
	integer                   :: lmax
	
	lmax = 3*nside - 1
	
	write (*,'(/,X,"Smoothing input map with Gaussian beam.")')
	
	allocate(alm(1,0:lmax,0:lmax))
	
	! The subroutines used here are designed for maps in RING ordering
	if (ord == 2) call convert_nest2ring(nside, map_in)
	
	call map2alm(nside, lmax, lmax, map_in, alm)
	call alter_alm(nside, lmax, lmax, dble(fwhm), alm)
	call alm2map(nside, lmax, lmax, alm, map_in)
	
	! Go back to NESTED ordering if necessary
	if (ord == 2) call convert_ring2nest(nside, map_in)
	
end subroutine

! Empirical Mode Decomposition process
subroutine emd(nside, map_in, imf, nit, stiff, tens, map_out)
	integer, intent(in)   :: nside, imf, nit, stiff, tens
	real(DP), intent(in)  :: map_in(0:12*nside**2-1)
	real(DP), intent(out) :: map_out(0:12*nside**2-1,2*imf)
	integer               :: p, i, j, n, nlist, nmax, nmin, imax(12*nside**2/9), imin(12*nside**2/9)
	real(DP), allocatable :: inp(:), Emax(:), Emin(:)
	
	n = nside2npix(nside) - 1
	
	! IMF number
	do i = 1, imf
		write (*,'(/, X, "- Computing Intrinsic Mode Function ", I0, "/", I0, ".")') i, imf
		
		! Initialize residue
		if (i == 1) map_out(:,imf+i) = map_in
		if (i /= 1) map_out(:,imf+i) = map_out(:,imf+i-1)
		
		! Iteration number
		do j = 1, nit
			write (*,'(/, X, "-- Iteration in progress: " I0, "/", I0, ".")') j, nit
			
			allocate(inp(0:n), Emax(0:n), Emin(0:n), source=0.0)
			
			! Start sifting process
			if (j == 1) inp = map_out(:,imf+i)
			if (j /= 1) inp = map_out(:,i)
			
			! Finding the positions and values of the local extrema of the input map
			call local_extrema(nside, inp, nmax, nmin, imax, imin)

			! Compute the smooth envelopes (3 ways)
            
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
			! Interpolation by spherical spline with stiffness
			call ss_interp(nside, nLeg, stiff, tens, nmax, inp(imax(1:nmax)), imax(1:nmax), Emax)
			call ss_interp(nside, nLeg, stiff, tens, nmin, inp(imin(1:nmin)), imin(1:nmin), Emin)
			
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
			! Local interpolation by spherical spline with stiffness
			!!$OMP PARALLEL DO
			!do p = 0, n; call loc_interp(nside, p, nLeg, stiff, tens, inp(imax(1:nmax)), nmax, imax(1:nmax), Emax(p)); end do
			!!$OMP END PARALLEL DO
			!!$OMP PARALLEL DO
			!do p = 0, n; call loc_interp(nside, p, nLeg, stiff, tens, inp(imin(1:nmin)), nmin, imin(1:nmin), Emin(p)); end do
			!!$OMP END PARALLEL DO
			
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
			! Faster local interpolation by spherical spline with stiffness
			!call loc_interp_mg(nside, gridres, nLeg, stiff, tens, inp(imax(1:nmax)), nmax, imax(1:nmax), Emax)
			!call loc_interp_mg(nside, gridres, nLeg, stiff, tens, inp(imin(1:nmin)), nmin, imin(1:nmin), Emin)
			
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
			write (*,'(X, X, X, "Iteration finished with stop criteria = ", E10.4)') maxval(abs((Emax + Emin) / 2.0))
			
			! Update IMF
			map_out(:,i) = inp - (Emax + Emin) / 2.0
			
			deallocate(inp, Emax, Emin)
		
		end do
		
		! Update residue
		map_out(:,imf+i) = map_out(:,imf+i) - map_out(:,i)
		
	end do
	
end subroutine emd

! Empirical Mode Decomposition process
subroutine extrema(nside, map_in, stiff, tens, map_out)
	integer, intent(in)   :: nside, stiff, tens
	real(DP), intent(in)  :: map_in(0:12*nside**2-1)
	real(DP), intent(out) :: map_out(0:12*nside**2-1,4)
	integer               :: nmax, nmin, imax(12*nside**2/9), imin(12*nside**2/9)
	
	! Finding the positions and values of the local extrema of the input map
	call local_extrema(nside, map_in, nmax, nmin, imax, imin)
	
	! Minima and maxima
	write (*,'(/,X,A)') "- Finding local maxima."
	map_out(:,1) = HPX_DBADVAL; forall (i=1:nmax) map_out(imax(i),1) = map_in(imax(i))
	write (*,*) "- Finding local minima."
	map_out(:,2) = HPX_DBADVAL; forall (i=1:nmin) map_out(imin(i),2) = map_in(imin(i))
	
	! Compute the smooth envelopes (interpolation by spherical spline with stiffness)
	write (*,*) "- Computing upper envelope."
	call ss_interp(nside, nLeg, stiff, tens, nmax, map_in(imax(1:nmax)), imax(1:nmax), map_out(:,3))
	write (*,*) "- Computing lower envelope."
	call ss_interp(nside, nLeg, stiff, tens, nmin, map_in(imin(1:nmin)), imin(1:nmin), map_out(:,4))
	
end subroutine extrema

! Positions and values of the local extrema of an input map
subroutine local_extrema(nside, map_in, nmax, nmin, imax, imin)
	integer, intent(in)  :: nside
	real(DP), intent(in) :: map_in(0:12*nside**2-1)
	integer, intent(out) :: nmax, nmin, imax(12*nside**2/9), imin(12*nside**2/9)
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
subroutine loc_interp(nside, pix, lmax, stiff, tens, lut, next, iext, z)
	integer, intent(in)   :: nside, pix, lmax, stiff, tens, next, iext(next)
	real(DP), intent(in)  :: lut(next)
	real(DP), intent(out) :: z
	real(DP)              :: vpix(3), rpix, radius, const
	real(DP), allocatable :: A(:,:), B(:), v(:,:)
	integer               :: nle, i, j, l, nlist
	integer, allocatable  :: list(:), lutext(:), le(:)
	
	! Initial number of local extrema, radius of one single pixel and initial radius
	nle = 0; rpix = sqrt(4.0*pi/nside2npix(nside)) * sqrt(2.0)/2; radius = rpix
	
	! Cartesian vector coordinates of pixel "pix" on the coarse map
	call pix2vec_nest(nside, pix, vpix)
	
	! Find the adequate number of local extrema around pixel "pix"
	do while (nle == 0)
		! List size
		l = nside2npix(nside) * sin(radius/2.0)**2
		
		! Pixel indices of a neighborhood around pixel "pix"
		allocate(list(0:3*l/2), le(1:3*l/4)); call query_disc(nside, vpix, radius, list, nlist, nest=1)
		
		! If "list(i)" is a local extremum, then store it into the array "le"
		do i = 0, nlist-1; if (any(iext == list(i))) then; nle = nle + 1; le(nle) = list(i); end if; end do
		
		! Increase radius and reset parameters if there are too few local extrema
		if (nle < nlemin) then; radius = radius + rpix; nle = 0; deallocate(list, le); end if
		
	end do
	
	!write (*,'(X,"--- Computing kernel matrix of dimension ",I0,".")') nle
	
	! Indices of "lut" corresponding to the local extrema
	allocate(lutext(nle)); forall (i = 1:nle) lutext(i) = findloc(iext, value=le(i), dim=1)
	
	! Computing the interpolation coefficients stored in "X" by solving "A*X=B"
	allocate(A(nle,nle), B(nle), v(3,nle))
	
	! Constant parameter in interpolation function
	if (tens == 0) const = sum(lut(lutext)) / nle; if (tens /= 0) const = 0.0
	
	! Vector coordinates of the extrema in "le"
	do i = 1, nle; call pix2vec_nest(nside, iext(lutext(i)), v(:,i)); end do
	
	! Computing "A" and "B"
	forall (i=1:nle, j=1:nle) A(i,j) = G(dot_product(v(:,i),v(:,j)),lmax,stiff,tens); B = lut(lutext) - const
	
	! Interpolation coefficients "X=A^(-1)*B"
	call lsolve(nle, A, B)
	
	! Interpolated value at pixel "pix"
	z = const; do j = 1, nle; z = z + B(j) * G(dot_product(vpix,v(:,j)),lmax,stiff,tens); end do
	
	deallocate(list, le, lutext, A, B, v)
	
end subroutine loc_interp

! Multigrid local interpolation by spherical spline (Perrin et al, 1988)
subroutine loc_interp_mg(nside, nsidec, lmax, stiff, tens, lut, next, iext, map_out)
	integer, intent(in)   :: nside, nsidec, lmax, stiff, tens, next, iext(next)
	real(DP), intent(in)  :: lut(next)
	real(DP), intent(out) :: map_out(0:12*nside**2-1)
	real(DP)              :: ang(2), outvals((nside/nsidec)**2)
	integer               :: n, nc, pc, nle, i
	integer, allocatable  :: list(:), lutext(:), le(:), coarse(:,:), count(:)
	
	n = nside2npix(nside) - 1; nc = nside2npix(nsidec) - 1
	
	write (*,'(/,X, "Parameters for multigrid local interpolation:")')
	write (*,'(X, "Coarser Nside = ", I0)') gridres
	write (*,'(X, "Fine pixels per coarse pixel = ", I0)') (nside/gridres)**2
	
	! Coarse map
	allocate(coarse(0:nc,(nside/nsidec)**2), count(0:nc)); count = 0
	
	! Map pixels indices from original Nside to pixels in lower Nside
	do i = 0, n
		! Position of pixel "i" on the original map and corresponding pixel "pc" on the coarse map
		call pix2ang_nest(nside, i, ang(1), ang(2)); call ang2pix_nest(nsidec, ang(1), ang(2), pc)
		
		! Count and store pixels on the original map within pixel "pc" on the coarse map
		count(pc) = count(pc) + 1; coarse(pc,count(pc)) = i
	end do
	
	!write (*,'(X,X,X,X, "COUNT min = ", I0, " max = ", I0)') maxval(count), minval(count)
	
	! Calculate interpolated values for pixels of the original map within pixel "i" in the coarse map
	!$OMP PARALLEL DO PRIVATE(outvals)
	do i = 0, nc; call loc_interp_mg_val(nside, nsidec, coarse(i,:), lmax, stiff, tens, lut, next, iext, outvals, i); map_out(coarse(i,:)) = outvals; end do
	!$OMP END PARALLEL DO
	
	deallocate(coarse, count)
	
end subroutine loc_interp_mg

! Pixel value for multigrid local interpolation by spherical spline
subroutine loc_interp_mg_val(nside, nsidec, ind, lmax, stiff, tens, lut, next, iext, outvals, pix)
	integer, intent(in)   :: nside, nsidec, ind((nside/nsidec)**2), lmax, stiff, tens, next, iext(next), pix
	real(DP), intent(in)  :: lut(next)
	real(DP), intent(out) :: outvals((nside/nsidec)**2)
	real(DP)              :: ang(2), vpix(3), rpix, radius, const
	real(DP), allocatable :: A(:,:), B(:), v(:,:)
	integer               :: nle, i, j, k, l, nlist
	integer, allocatable  :: list(:), lutext(:), le(:)
	
	! Initial number of local extrema, radius of one single pixel and initial radius
	nle = 0; rpix = sqrt(4.0*pi / nside2npix(nsidec)) * sqrt(2.0)/2; radius = rpix
	
	! Cartesian vector coordinates of pixel "pix" on the coarse map
	call pix2vec_nest(nsidec, pix, vpix)
	
	! Find the adequate number of local extrema around pixel "pix" in the coarse map
	do while (nle == 0)
		! List size
		l = nside2npix(nside) * sin(radius/2.0)**2
		
		! Pixel indices of a neighborhood around pixel "pix"
		allocate(list(0:3*l/2), le(1:3*l/4)); call query_disc(nside, vpix, radius, list, nlist, nest=1)
		
		! If "list(i)" is a local extremum, then store it into the array "le"
		do i = 0, nlist-1; if (any(iext == list(i))) then; nle = nle + 1; le(nle) = list(i); end if; end do
		
		! Increase radius and reset parameters if there are too few local extrema
		if (nle < nlemin) then; radius = radius * 2; nle = 0; deallocate(list, le); end if
		
	end do
	
	write (*,'(X,"--- Computing kernel matrix of dimension ",I0,".")') nle
	
	allocate(lutext(nle), A(nle,nle), B(nle), v(3,nle))
	
	! Indices of "lut" corresponding to the local extrema
	forall (i = 1:nle) lutext(i) = findloc(iext, value=le(i), dim=1)
	
	! Constant parameter in interpolation function
    	if (tens == 0) const = sum(lut(lutext)) / nle; if (tens /= 0) const = 0.0
	
    	! Vector coordinates of the extrema in "le"
    	do i = 1, nle; call pix2vec_nest(nside, iext(lutext(i)), v(:,i)); end do
	
    	! Computing the interpolation coefficients stored in "X" by solving "A*X=B"
    	forall (i=1:nle, j=1:nle) A(i,j) = G(dot_product(v(:,i),v(:,j)),lmax,stiff,tens); B = lut(lutext) - const
	
    	! Interpolation coefficients "X=A^(-1)*B"
    	call lsolve(nle, A, B)
	
    	! Interpolated value for all pixels on the original map within pixel "pix" on the coarse map
    	do i = 1, (nside/nsidec)**2
		! Position of the pixel on the original map
		call pix2vec_nest(nside, ind(i), vpix)
        	
        	! Calculation of the interpolated value
        	outvals(i) = const; do j = 1, nle; outvals(i) = outvals(i) + B(j) * G(dot_product(vpix,v(:,j)),lmax,stiff,tens); end do
    	end do

    	deallocate(list, le, lutext, A, B, v)
    	
end subroutine loc_interp_mg_val

! Global interpolation by spherical spline (Perrin et al, 1988)
subroutine ss_interp(nside, lmax, stiff, tens, next, lut, iext, map_out)
	integer, intent(in)   :: nside, lmax, stiff, tens, next, iext(next)
	real(DP), intent(in)  :: lut(next)
	real(DP), intent(out) :: map_out(0:12*nside**2-1)
	real(DP)              :: A(next,next), B(next), v(3,next), vi(3), const
	integer               :: i, j, p, n
	
	n = nside2npix(nside) - 1
    
    	! Constant parameter in interpolation function
    	if (stiff == 0 .and. tens /= 0) then; const = 0.0; else; const = sum(lut) / next; end if
	
	! Vector coordinates of the extrema in "le"
	do i = 1, next; call pix2vec_nest(nside, iext(i), v(:,i)); end do
    	
    	write (*,'(/,X,"--- Computing kernel matrix of dimension ",I0,".")') next
	
	! Computing "A" and "B"
    	!$OMP PARALLEL DO
    	do j = 1, next; do i = 1, next; A(i,j) = G(dot_product(v(:,i),v(:,j)),lmax,stiff,tens); end do; end do
    	!$OMP END PARALLEL DO
    	B = lut - const
    
    	write (*,'(X,A)') "--- Solving system of equations."
    	
	! Interpolation coefficients "X=A^(-1)*B"
	call lsolve(next, A, B)
	
	write (*,'(X,A)') "--- Interpolation started."
    	
	!$OMP PARALLEL DO PRIVATE(vi)
	do i = 0, n
        	! Position of the pixel on the map
        	call pix2vec_nest(nside, i, vi)
        
        	! Calculation of the interpolated value
        	map_out(i) = const; do j = 1, next; map_out(i) = map_out(i) + B(j) * G(dot_product(vi,v(:,j)),lmax,stiff,tens); end do
    	end do
	!$OMP END PARALLEL DO
    
	write (*,'(X,A,/)') "--- Interpolation finished."
	
end subroutine ss_interp

! Basis function for interpolation by spherical spline
pure function G(x, lmax, stiff, tens)
	integer, intent(in)   :: lmax, stiff, tens
	real(DP), intent(in)  :: x
	real(DP)              :: v3(3), P(0:2), S, theta, G
	integer               :: i
	
    	G = - 1.0 / (4*pi); P(0) = 1.0; P(1) = x
    	
    	if (stiff == 0 .and. tens /= 0) then
        	S = P(0) * dble(2*0+1) / ((0*(0+1))+tens) + P(1) * dble(2*1+1) / ((1*(1+1))+tens)
        	do i = 2, lmax
            		P(mod(i,3)) = ((2*i-1) * x * P(mod(i-1,3)) - (i-1) * P(mod(i-2,3))) / i
            		S = S + P(mod(i,3)) * dble(2*i+1) / ((i*(i+1))+tens)
        	end do
        	G = G * S
    	else
        	S = P(1) * dble(2*1+1) / (((1*(1+1))+tens) * (1*(1+1))**stiff)
        	do i = 2, lmax
            		P(mod(i,3)) = ((2*i-1) * x * P(mod(i-1,3)) - (i-1) * P(mod(i-2,3))) / i
            		S = S + P(mod(i,3)) * dble(2*i+1) / (((i*(i+1))+tens) * (i*(i+1))**stiff)
        	end do
        	G = (-1)**stiff * G * S
    	end if
	
end function G

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

subroutine usage()
    	write (*,'(/,X,A,/)') "Usage:"
    	write (*,*) "For Hilbert-Huang transform: hht IFN IMF NIT STF TNS GSP"
    	write (*,'(X,A,/)') "For finding local extrema: hht -lext IFN STF TNS SMP"
    	write (*,*) "IFN = Input file name"
    	write (*,*) "IMF = Number of Intrinsic Mode Functions"
    	write (*,*) "NIT = Number of iterations in the Empirical Mode Decomposition"
    	write (*,*) "TNS = Tension parameter"
    	write (*,*) "GSP = FWHM parameter for Gaussian smoothing (in arcminutes)"
    	write (*,'(X,A,/)') "STF = Stiffness"
	
end subroutine usage

end program
