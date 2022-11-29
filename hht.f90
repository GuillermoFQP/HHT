! Empirical Mode Decomposition for HEALPix maps 

program hht

use healpix_modules

implicit none

!======================================================================================
real(DP), allocatable          :: map_in(:,:), map_out(:,:,:)
integer                        :: nside, npix, nmaps, ord, n, imf, nit, ch, i, stff, tens, fwhm, lmax
character(len=80)              :: fin, arg, header(43)
character(len=80), allocatable :: fout(:)
!======================================================================================
! Parameters for the calculations of spherical harmonics
integer, parameter  :: LOG2LG = 100, RSMAX = 20, RSMIN = -20
real(DP), parameter :: FL_LARGE = 2.0**LOG2LG, FL_SMALL = 2.0**(-LOG2LG), ALN2_INV = 1.4426950408889634073599246810 ! 1/log(2)
real(DP)            :: rescale_tab(RSMIN:RSMAX)
!======================================================================================

select case (nArguments())
	case (1)
		! Display usage information
		call getArgument(1, arg)
		if (arg == "-h") then
			call usage()
			stop
		end if
		call fatal_error("Cannot parse argument.")
	case (5)
		! Generate maps of local extrema and their corresponding interpolations
		call getArgument(1, arg)
		if (arg /= "-lext") call fatal_error("Cannot parse argument.")
		call getArgument(2, fin)
		call getArgument(3, arg) ! Stiffness
		read(arg,*) stff
        call getArgument(4, arg) ! Tension parameter
		read(arg,*) tens
        call getArgument(5, arg) ! FWHM for Gaussian smoothing
		read(arg,*) fwhm
		
		! Output file names
		ch = 4
		allocate(fout(ch))
		write (fout(1),*) "loc_max.fits"
		write (fout(2),*) "loc_min.fits"
		write (fout(3),*) "int_max.fits"
		write (fout(4),*) "int_min.fits"
	case (6)
		! Hilbert-Huang transform
		call getArgument(1, fin) ! Input file name
		call getArgument(2, arg) ! Number of IMFs
		read(arg,*) imf
		call getArgument(3, arg) ! Number of iterations in the EMD
		read(arg,*) nit
		call getArgument(4, arg) ! Stiffness
		read(arg,*) stff
        call getArgument(5, arg) ! Tension parameter
		read(arg,*) tens
        call getArgument(6, arg) ! FWHM for Gaussian smoothing
		read(arg,*) fwhm
		
		! Output file names
		ch = imf
		allocate(fout(2*ch))
		
		do i = 1, ch
			write (fout(i),'("imf",I0,".fits")') i
			write (fout(ch+i),'("res",I0,".fits")') i
		end do
	case default
		call fatal_error("Invalid number of arguments.")
end select

! Parameters of the FITS file containing the input map
npix = getsize_fits(fin, nmaps=nmaps, nside=nside, ordering=ord)
n = nside2npix(nside) - 1
lmax = min(3*nside - 1, 2*2048)

write (*,'(/,X, "Input map Nside = ", I0)') nside

! Allocating arrays
allocate(map_in(0:n,nmaps), map_out(0:n,nmaps,ch), source=0.0)

! Reading input map
call input_map(fin, map_in, npix, nmaps)
write (*,'(/,X,A)') "Map read successfully."

! Gaussian smoothing the input map
if (fwhm /= 0) call smoothing(nside, ord, lmax, map_in(:,1), fwhm)

! All the following subroutines are designed for maps in NESTED ordering
if (ord == 1) call convert_ring2nest(nside, map_in)

! Find extrema or perform empirical mode decomposition
write (*,'(/,X,"Interpolation with stiffness parameter ",I0," and tension parameter ",I0,".")') stff, tens
if (nArguments() == 5) then
	call extrema(nside, map_in(:,1), stff, tens, lmax, fwhm, map_out(:,1,:))
	write (*,'(/, X, A)') "Extrema and envelopes obtained successfully."
else
	call emd(nside, map_in(:,1), imf, nit, stff, tens, lmax, fwhm, map_out(:,1,:))
	write (*,'(/, X, A)') "Empirical Mode Decomposition completed successfully."
end if

! Go back to RING ordering if necessary
if (ord == 1) then
	call convert_nest2ring(nside, map_in)
	do i = 1, ch
		call convert_nest2ring(nside, map_out(:,:,i))
	end do
end if

! Generating output files
call write_minimal_header(header, 'map', nside=nside, order=ord)
do i = 1, ch
	call output_map(map_out(:,:,i), header, fout(i))
end do

! Generating files containing the residues and the EMD
if (nArguments() == 6) then
	do i = 1, ch
		call output_map(map_in - sum(map_out(:,:,1:i),dim=3), header, fout(ch+i))
	end do
end if
write (*,*) "Output files generated successfully."

deallocate(map_in, map_out, fout)

contains

! Gaussian smoothing process
subroutine smoothing(nside, ord, lmax, map_in, fwhm)
	integer, intent(in)       :: nside, ord, fwhm
	real(DP), intent(inout)   :: map_in(0:12*nside**2-1)
	complex(DPC), allocatable :: alm(:,:,:)
	integer                   :: lmax
	
	write (*,'(/,X,"Smoothing input map with Gaussian beam.")')
	
	allocate(alm(1,0:lmax,0:lmax))
	
	! The subroutines used here are designed for maps in RING ordering
	if (ord == 2) call convert_nest2ring(nside, map_in)
	
	call map2alm(nside, lmax, lmax, map_in, alm)
	call alter_alm(nside, lmax, lmax, dble(fwhm), alm)
	call alm2map(nside, lmax, lmax, alm, map_in)
	
	deallocate(alm)
	
	! Go back to NESTED ordering if necessary
	if (ord == 2) call convert_ring2nest(nside, map_in)
	
end subroutine smoothing

! Positions and values of the local extrema of an input map
subroutine local_extrema(nside, map_in, nmax, nmin, imax, imin)
	integer, intent(in)   :: nside
	real(DP), intent(in)  :: map_in(0:12*nside**2-1)
	integer, intent(out)  :: nmax, nmin, imax(12*nside**2/9), imin(12*nside**2/9)
	integer               :: n, i, j, l, nlist
	integer, allocatable  :: list(:)
	real(DP)              :: rpix, radius, vi(3)
	logical               :: use_disc
	
	n = nside2npix(nside) - 1
	
	nmin = 0 ! Local minima counter
	nmax = 0 ! Local maxima counter
	
	! Use "query_disc" subroutine to find local extrema, otherwise use "neighbours_nest"
	use_disc = .false.
	
	if (use_disc) then
		!======================================================================================
		! Finding local extrema using disk with given radius
		!======================================================================================
		rpix = sqrt(4.0*pi/(n+1)) / sqrt(2.0) ! Approximate radius of one pixel
		radius = 5.0*rpix                       ! Disk radius
		l = (2.0*radius)**2 / (4.0*pi/(n+1))    ! Approximate number of pixels inside the disk
		
		allocate(list(0:l))
		
		do i = 0, n
			! Pixel indices of a neighborhood around pixel "i"
			call pix2vec_nest(nside, i, vi)
			call query_disc(nside, vi, radius, list, nlist, nest=1)
			
			! Find local extrema
			if (map_in(i) >= maxval(map_in(list(1:nlist)))) then
				nmax = nmax + 1
				imax(nmax) = i
			end if
			if (map_in(i) <= minval(map_in(list(1:nlist)))) then
				nmin = nmin + 1
				imin(nmin) = i
			end if
		end do
		!======================================================================================
	else
		!======================================================================================
		! Finding local extrema using nearest neighbours
		!======================================================================================
		allocate(list(8))
		
		do i = 0, n
			! Pixel indices of a neighborhood around pixel "i"
			call neighbours_nest(nside, i, list, nlist)
			
			! Find local extrema
			if (map_in(i) >= maxval(map_in(list(1:nlist)))) then
				nmax = nmax + 1
				imax(nmax) = i
			end if
			if (map_in(i) <= minval(map_in(list(1:nlist)))) then
				nmin = nmin + 1
				imin(nmin) = i
			end if
		end do
		!======================================================================================
	end if
	
	deallocate(list)
	
end subroutine local_extrema

! Solve the real system of "n" symmetric linear equations in "n" unknowns in the form "A*X=B"
subroutine lsolve(n, A, B)
	integer, intent(in)     :: n
	real(DP), intent(in)    :: A(n,n)
	real(DP), intent(inout) :: B(n)
	real(DP), allocatable   :: work(:)
	real(DP)                :: get_lwork(1)
	integer                 :: pivot(n), stat, lwork
	
	stat = 0 ! Status indicator ("stat/=0" indicates an error)
	
	! Returns the optimal size of the WORK array as the first entry of the GET_LWORK array
	call dsysv('U', n, 1, A, n, pivot, B, n, get_lwork, -1, stat)
	
	lwork = get_lwork(1) ! Parameter to calculate the optimal size of the WORK array
	
	allocate(work(lwork))
	
	! LAPACK subroutine to get "X=A^(-1)*B" for a symmetric matrix "A"
	call dsysv('U', n, 1, A, n, pivot, B, n, work, lwork, stat)
	
	! Stop the program if necessary
	if (stat /= 0) call fatal_error('Singular matrix in subroutine LSOLVE()')
	
	deallocate(work)
	
end subroutine lsolve

! First factor for interpolation in harmonic space
subroutine interp_alms(c, ang, lmax, mmax, alm)
	integer, intent(in)       :: lmax, mmax
	real(DP), intent(in)      :: c(:), ang(size(c),2)
	complex(DPC), intent(out) :: alm(1,0:lmax,0:mmax)
	real(DP)                  :: mfac(0:mmax), recfac(0:1,0:lmax), lam_lm(0:lmax)
	integer                   :: l, m
	
	alm = (0.0, 0.0)
	
	! Recursion factors used in "lambda_mm" calculation for all "m" in "0<=m<=m_max"
	call gen_mfactor(mmax, mfac)
	
	!$OMP PARALLEL PRIVATE(m, i, recfac, lam_lm, l) SHARED(mmax, lmax, c, ang, mfac, alm)
	!$OMP DO SCHEDULE(DYNAMIC)
	do m = 0, mmax
		! Generate recursion factors useful for "lambda_lm" for a given "m"
		call gen_recfactor(lmax, m, recfac)
		
		do i = 1, size(c)
			! Compute "lam_lm(theta_i)" for all "l>=m" for a given "m"
			call do_lambda_lm(lmax, m, abs(cos(ang(i,1))), sin(ang(i,1)), mfac(m), recfac, lam_lm)
			if (cos(ang(i,1)) < 0.0) forall (l=m:lmax) lam_lm(l) = (-1.0)**(l+m) * lam_lm(l)
			
			! Compute numerators for all "l>=m" for a given "m"
			alm(1,m:lmax,m) = alm(1,m:lmax,m) + cmplx(c(i) * lam_lm(m:lmax), kind=DPC) * cdexp((0.0, -1.0) * m * ang(i,2))
		end do
	end do
	!$OMP END DO
	!$OMP END PARALLEL
	
end subroutine interp_alms

! Global interpolation by spherical spline (Perrin et al, 1988)
subroutine ss_interp(nside, lmax, stff, tens, next, LUT, iext, map_out)
	integer, intent(in)       :: nside, lmax, stff, tens, next, iext(next)
	real(DP), intent(in)      :: LUT(next)
	real(DP), intent(out)     :: map_out(0:12*nside**2-1)
	real(DP)                  :: vec(next,3), ang(next,2), vi(3), fwhm
	integer                   :: i, j, n, lmin, omega_pix
	real(DP), allocatable     :: A(:,:), B(:), bl1(:,:), bl2(:,:), wl(:,:)
	complex(DPC), allocatable :: alm(:,:,:)
	logical                   :: harmonic_space
	integer(I8B)              :: l
	
	n = nside2npix(nside) - 1
	fwhm = 2.0 * acos(1.0 - 2.0/next) * 180.0 * 60.0 / pi
	
	allocate(A(next,next), B(next), source=0.0)                      ! System of linear equations
	allocate(bl1(0:lmax,1), bl2(0:lmax,1), wl(0:lmax,1), source=1.0) ! Beams
	
	! Positions of the local extrema
	do i = 1, next
		call pix2vec_nest(nside, iext(i), vec(i,:))
		call pix2ang_nest(nside, iext(i), ang(i,1), ang(i,2))
	end do
	
	! Initialize RESCALE_TAB useful for calculation of spherical harmonics
	call init_rescale_tab()
	
	! Inverse of differential operator in harmonic space
	if (stff /= 0 .or. tens == 0) then
		lmin = 1
		bl1(0,1) = 1.0
	else
		lmin = 0
		bl1(0,1) = -1.0 / tens
	end if
	
	forall (l=1:lmax) bl1(l,1) = dble((-1)**(stff+1)) / dble((l*(l+1) + tens) * (l*(l+1))**stff)
	
	! Pixel window function
	call pixel_window(wl, nside)
	
	! Gaussian beam
	call gaussbeam(dble(fwhm), lmax, bl2)
	
	write (*,'(/, X, "--- Computing interpolation matrix of dimension ", I0, ".")') next
	
	! Computing "A"
	!$OMP PARALLEL PRIVATE(i, j) SHARED(A, vec, lmin, lmax)
	!$OMP DO SCHEDULE(DYNAMIC)
	do j = 1, next
		do i = 1, j
			!A(i,j) = Aij(dot_product(vec(i,:),vec(j,:)), wl*bl1, lmin, lmax)
			A(i,j) = G(dot_product(vec(i,:),vec(j,:)), wl*bl1, lmin, lmax)
			if (i /= j) A(j,i) = A(i,j)
		end do
	end do
	!$OMP END DO
	!$OMP END PARALLEL
	
	! Computing "B"
	B = LUT - lmin * sum(LUT)/next

	write (*,'(X, A)') "--- Solving system of equations."

	! Calculating and storing the interpolation coefficients "X=A^(-1)*B" in array "B"
	call lsolve(next, A, B)
	
	deallocate(A)
	
	harmonic_space = .true.
	
	if (harmonic_space) then
		!======================================================================================
		! Interpolation in spherical harmonic space
		!======================================================================================
		write (*,'(X, A)') "--- Calculating spherical harmonic coefficients."
		
		allocate(alm(1,0:lmax,0:lmax), source=(0.0, 0.0))
		
		! Calculating the sum of "C_i*Y^(*)_lm(theta_i,phi_i)" over all pixels "i" with values in LUT
		call interp_alms(B, ang, lmax, lmax, alm)
		
		! Convolution in harmonic space and generation of output map
		forall (l=lmin:lmax) alm(1,l,0:l) = alm(1,l,0:l) * wl(l,1) * bl1(l,1) ! * bl2(l,1)
		
		! "a_00" takes a different value when stff>0 or tens=0
		if (lmin == 1) alm(1,0,0) = sqrt(4.0*pi) * sum(LUT) / next
		
		write (*,'(X, A)') "--- Generating map."
		
		call alm2map(nside, lmax, lmax, alm, map_out) ! Generates map in RING ordering
		call convert_ring2nest(nside, map_out)        ! Goes back to NESTED ordering

		write (*,'(X, "--- Max. and min. interpolation error: ", E10.4, X, E10.4)') maxval(abs(LUT-map_out(iext))), minval(abs(LUT-map_out(iext)))
		
		deallocate(B, wl, bl1, bl2, alm)
		!======================================================================================
	else
		!======================================================================================
		! Interpolation in real space
		!======================================================================================
		write (*,'(X, A)') "--- Interpolation in real space started."
		
		map_out = lmin * sum(LUT)/next
		
		!$OMP PARALLEL PRIVATE (i, vi) SHARED(n, nside, map_out, lmin, lmax)
		!$OMP DO SCHEDULE(DYNAMIC)
		do i = 0, n
			! Position of the pixel on the map
			call pix2vec_nest(nside, i, vi)
			
			! Calculation of the interpolated value
			do j = 1, next
				map_out(i) = map_out(i) + B(j) * Aij(dot_product(vi,vec(j,:)), wl*bl1*bl2, lmin, lmax)
				!map_out(i) = map_out(i) + B(j) * G(dot_product(vi,vec(j,:)), wl*bl1*bl2, lmin, lmax)
			end do
		end do
		!$OMP END DO
		!$OMP END PARALLEL
		
		deallocate(B, wl, bl1, bl2)
		
		write (*,'(X, "--- Max. and min. interpolation error: ", E10.4, X, E10.4)') maxval(abs(LUT-map_out(iext))), minval(abs(LUT-map_out(iext)))
		!======================================================================================
	end if
	
end subroutine ss_interp

! Computes the upper and lower envelopes of an input map
subroutine extrema(nside, map_in, stff, tens, lmax, fwhm, map_out)
	integer, intent(in)   :: nside, stff, tens, lmax, fwhm
	real(DP), intent(in)  :: map_in(0:12*nside**2-1)
	real(DP), intent(out) :: map_out(0:12*nside**2-1,4)
	integer               :: nmax, nmin, imax(12*nside**2/9), imin(12*nside**2/9)
	
	! Finding the positions and values of the local extrema of the input map
	write (*,'(/, X, A)') "- Finding local extrema."
	call local_extrema(nside, map_in, nmax, nmin, imax, imin)
	
	map_out(:,1:2) = HPX_DBADVAL
	
	map_out(imax(1:nmax),1) = map_in(imax(1:nmax))
	map_out(imin(1:nmin),2) = map_in(imin(1:nmin))
	
	write (*,'(/, X, A)') "- Computing upper envelope."
	call ss_interp(nside, lmax, stff, tens, nmax, map_in(imax(1:nmax)), imax(1:nmax), map_out(:,3))
	write (*,'(/, X, A)') "- Computing lower envelope."
	call ss_interp(nside, lmax, stff, tens, nmin, map_in(imin(1:nmin)), imin(1:nmin), map_out(:,4))
	
end subroutine extrema

! Empirical Mode Decomposition process
subroutine emd(nside, map_in, nimf, nitr, stff, tens, lmax, fwhm, imf)
	integer, intent(in)   :: nside, nimf, nitr, stff, tens, lmax, fwhm
	real(DP), intent(in)  :: map_in(0:12*nside**2-1)
	real(DP), intent(out) :: imf(0:12*nside**2-1,nimf)
	integer               :: i, j, k, n, nmax, nmin, imax(12*nside**2/9), imin(12*nside**2/9)
	real(DP), allocatable :: inp(:), Emax(:), Emin(:)
	
	n = nside2npix(nside) - 1
	
	allocate(inp(0:n), Emax(0:n), Emin(0:n), source=0.0)
	
	! IMF number
	do i = 1, nimf
		write (*,'(/, X, "- Computing Intrinsic Mode Function ", I0, "/", I0, ".")') i, nimf
		
		! Iteration number
		do j = 1, nitr
			write (*,'(/, X, "-- Iteration in progress: " I0, "/", I0, ".")') j, nitr
			
			! Start sifting process
			if (j == 1) then
				inp = map_in
				if (i /= 1) then
					do k = 1, i - 1
						inp = inp - imf(:,k)
					end do
				end if
			else
				inp = imf(:,i)
			end if
			
			! Finding the positions and values of the local extrema of the input map
			call local_extrema(nside, inp, nmax, nmin, imax, imin)

			! Computing smooth upper and lower envelopes
			call ss_interp(nside, lmax, stff, tens, nmax, inp(imax(1:nmax)), imax(1:nmax), Emax)
			call ss_interp(nside, lmax, stff, tens, nmin, inp(imin(1:nmin)), imin(1:nmin), Emin)
			
			write (*,'(/, X, "-- Stop criteria = ", E10.4)') sqrt(sum(abs((Emax + Emin) / 2.0)**2) / dble(n+1))
			
			! Update IMF
			imf(:,i) = inp - (Emax + Emin) / 2.0
		end do
	end do
	
	deallocate(inp, Emax, Emin)
	
end subroutine emd

! Calculation of the elements of the interpolation matrix (fast)
pure function G(cth, bl, lmin, lmax)
	integer, intent(in)  :: lmin, lmax
	real(DP), intent(in) :: cth, bl(0:lmax,1)
	real(DP)             :: theta_ij, P(0:2), G
	integer(I8B)         :: l
	
	G = 0.0
	
	if (cth < 1.0 - epsilon(cth)) then
		P(0) = 1.0
		P(1) = cth
		
		! Initial terms of summation
		if (lmin == 0) G = G + P(0) * bl(0,1) / (4.0*pi)
		if (lmin <= 1) G = G + P(1) * bl(1,1) * 3.0 / (4.0*pi)
		
		! Summation
		do l = 2, lmax
			P(mod(l,3)) = (dble(2*l-1) * cth * P(mod(l-1,3)) - dble(l-1) * P(mod(l-2,3))) / l
			if (l >= lmin) G = G + P(mod(l,3)) * bl(l,1) * dble(2*l+1) / (4.0*pi)
		end do
	else
		! Initial terms of summation
		if (lmin == 0) G = G + bl(0,1) / (4.0*pi)
		if (lmin <= 1) G = G + bl(1,1) * 3.0 / (4.0*pi)
		
		! Summation
		do l = 2, lmax
			if (l >= lmin) G = G + bl(l,1) * dble(2*l+1) / (4.0*pi)
		end do
	end if
	
end function G

! Calculation of the elements of the interpolation matrix (slow)
function Aij(cth, bl, lmin, lmax)
	integer, intent(in)  :: lmin, lmax
	real(DP), intent(in) :: cth, bl(0:lmax,1)
	real(DP)             :: sth, mfac(0:0), recfac(0:1,0:lmax), lam_lm(0:lmax), Aij
	integer(I8B)         :: l
	
	Aij = 0.0
	sth = sqrt(1.0 - cth**2)
	
	! Recursion factor used in "lambda_00" calculation
	call gen_mfactor(0, mfac)
	
	! Generate recursion factors useful for "lambda_l0"
	call gen_recfactor(lmax, 0, recfac)
	
	if (cth < 1.0 - epsilon(cth)) then
		! Compute "lambda_l0(theta_ij)=sqrt((2*l+1)/4*pi)*Pl(cos(theta_ij))" for all "l" for "m=0"
		call do_lambda_lm(lmax, 0, abs(cth), sth, mfac(0), recfac, lam_lm)
		if (cth < 0.0) forall (l=0:lmax) lam_lm(l) = dble((-1)**l) * lam_lm(l)
	else
		forall (l=0:lmax) lam_lm(l) = sqrt(dble(2*l+1) / (4.0*pi))
	end if
	
	! Compute summation
	do l = lmin, lmax
		Aij = Aij + bl(l,1) * sqrt(dble(2*l+1) / (4.0*pi)) * lam_lm(l)
	end do
	
end function Aij

! Generates factor used in "lambda_mm" calculation for all "m" in "0<=m<=m_max"
subroutine gen_mfactor(m_max, m_fact)
	integer, intent(in)   :: m_max
	real(DP), intent(out) :: m_fact(0:m_max)
	integer               :: m

	! fact(m) = fact(m-1) * sqrt((2m+1)/(2m))
	m_fact(0) = 1.0
	do m = 1, m_max
		m_fact(m) = m_fact(m-1) * sqrt(dble(2*m+1) / dble(2*m))
	end do

	! Log_2 ( fact(m) / sqrt(4 Pi) )
	do m=0,m_max
		m_fact(m) = log(SQ4PI_INV * m_fact(m)) * ALN2_INV
	enddo

end subroutine gen_mfactor

! Generates recursion factors used to computes the "Y_lm" of degree "m" for all "l" in "m<=l<=l_max"
subroutine gen_recfactor( l_max, m, recfac)
	integer, intent(in)   :: l_max, m
	real(DP), intent(out) :: recfac(0:1, 0:l_max)
	real(DP)              :: fm2, fl2
	integer               :: l

	recfac(0:1,0:m) = 0.0
	fm2 = dble(m)**2
	do l = m, l_max
		fl2 = dble(l+1)**2
		recfac(0,l) = sqrt((4.0*fl2-1.0) / (fl2-fm2))
	end do
	! Put outside the loop because of problem on some compilers
	recfac(1,m:l_max) = 1.0 / recfac(0,m:l_max)

end subroutine gen_recfactor

! Initialize RESCALE_TAB array
subroutine init_rescale_tab()
	integer           :: s, smax
	real(DP)          :: logOVFLOW
	
	logOVFLOW=log(FL_LARGE)
	smax = int(log(MAX_DP) / logOVFLOW)
	rescale_tab(RSMIN:RSMAX) = 0.0
	do s = -smax, smax
		rescale_tab(s) = FL_LARGE**s
	end do
	rescale_tab(0) = 1.0
	
end subroutine init_rescale_tab

! Computes scalar "lambda_lm(theta)" for all "l" in "[m,lmax]" for a given "m", and given theta
subroutine do_lambda_lm(lmax, m, cth, sth, mfac, recfac, lam_lm)
	integer, intent(in)   :: lmax,  m
	real(DP), intent(in)  :: cth, sth, mfac, recfac(0:1,0:lmax)
	real(DP), intent(out) :: lam_lm(0:lmax)
	real(DP)              :: log2val, dlog2lg, ovflow, unflow, corfac, lam_mm, lam_0, lam_1, lam_2
	integer               :: scalel, l, l_min
	
	! Define constants
	ovflow = rescale_tab(1)
	unflow = rescale_tab(-1)
	l_min = l_min_ylm(m, sth)
	dlog2lg = real(LOG2LG, kind=DP)

	! Computes "lambda_mm"
	log2val = mfac + m*log(sth) * ALN2_INV     ! "log_2(lambda_mm)"
	scalel = int(log2val / dlog2lg)
	corfac = rescale_tab(max(scalel,RSMIN))
	lam_mm = 2.0**(log2val - scalel * dlog2lg) ! Rescaled "lambda_mm"
	if (IAND(m,1)>0) lam_mm = -lam_mm          ! Negative for odd "m"

	lam_lm(0:lmax) = 0.0
	
	! "l=m"
	lam_lm(m) = lam_mm * corfac

	! "l>m"
	lam_0 = 0.0
	lam_1 = 1.0
	lam_2 = cth * lam_1 * recfac(0,m)
	do l = m+1, lmax
		! Do recursion
		if (l >= l_min) then
			lam_lm(l) = lam_2 * corfac * lam_mm
		end if
		lam_0 = lam_1 * recfac(1,l-1)
		lam_1 = lam_2
		lam_2 = (cth * lam_1 - lam_0) * recfac(0,l)

		! Do dynamic rescaling
		if (abs(lam_2) > ovflow) then
			lam_1 = lam_1 * unflow
			lam_2 = lam_2 * unflow
			scalel = scalel + 1
			corfac = rescale_tab(max(scalel,RSMIN))
		else if (abs(lam_2) < unflow .and. abs(lam_2) /= 0.0) then
			lam_1 = lam_1 * ovflow
			lam_2 = lam_2 * ovflow
			scalel = scalel - 1
			corfac = rescale_tab(max(scalel,RSMIN))
		end if

	end do
end subroutine do_lambda_lm

! Display usage information
subroutine usage()
	write (*,'(/, X, A, /)') "Usage:"
	write (*,*) "For Hilbert-Huang transform: hht IFN IMF NIT STF TNS GSP"
	write (*,'(X, A, /)') "For finding local extrema: hht -lext IFN STF TNS SMP"
	write (*,*) "IFN = Input file name"
	write (*,*) "IMF = Number of Intrinsic Mode Functions"
	write (*,*) "NIT = Number of iterations in the Empirical Mode Decomposition"
	write (*,*) "TNS = Tension parameter"
	write (*,*) "GSP = FWHM parameter for Gaussian smoothing (in arcminutes)"
	write (*,'(X, A, /)') "STF = Stiffness"
	
end subroutine usage

end program
