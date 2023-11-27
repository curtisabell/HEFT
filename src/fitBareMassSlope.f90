! TODO:
! 1. generalise enabledMasses to allow for N masses,
!    and read in which is enabled from slopeFitMasses.data
! 2. probably a config file option for test mode so dont have to recompile
program fitBareMassSlope
    use kinds
    use numFort
    use heft
    use hamiltonian
    use bobyqa_module
    implicit none

    real(DP) :: testVar, tempVar

    ! Pion mass variation
    real(DP) :: m_pi2
    ! real(DP) :: m_pi2_max = m_pi0**2
    real(DP) :: m_pi2_max = 0.4_DP
    real(DP) :: m_pi2_increment
    integer  :: m_pi2_counter, m_pi2_counter_nImg
    integer  :: m_pi2_counter_end

    logical :: useStored = .true.
    ! real(DP), dimension(:), allocatable :: m_bare ! bare mass
    character (len=8) :: dummyString
    character (len=9) :: wrt_fmt, wrt_fmtInt

    real(DP) :: L0 = 2.9_DP, Lmax = 10.0_DP
    real(DP) :: L, L_increment ! fm
    integer  :: n_L_points = 100

    ! integer :: n_bare = 1! number of bare states
    ! integer :: n_ch = 2! number of channels
    integer :: n_k ! number of momentum states
    integer :: n_min_forced = 128 ! forced min matrix size
    integer :: n_mesh ! n_ch * n_k
    ! integer :: n_init_k ! 0 for odd parity, 1 for even
    integer :: N_H ! full Hamiltonian size
    integer :: ch_i, ch_j, bare_i, bare_j
    integer :: ik, ich, i, j, ii, jk, jj, i_bare, j_bare ! Loop indicies
    real(DP) :: Lambda_max
    integer :: jImg

    ! Slope fitting vars
    integer :: iBare_fit
    integer :: nMasses, iMasses
    integer :: dof
    integer :: maxFitIters
    real(DP) :: bare_chi2
    real(DP), parameter :: bareFitStepSize = 1.0_DP
    real(DP), dimension(:,:), allocatable :: alpha_bare
    real(DP), dimension(:), allocatable :: alpha_bare_lin
    real(DP), dimension(:), allocatable :: slope_low, slope_high
    real(DP), dimension(:), allocatable :: bareFitError
    real(DP), dimension(:,:), allocatable :: m_pi_lat
    real(DP), dimension(:,:), allocatable :: m_Bar_lat
    real(DP), dimension(:,:), allocatable :: dm_Bar_lat
    real(DP), dimension(:), allocatable :: m_Bare_fit
    logical :: doBareFitPrinting = .false.

    ! real(DP), parameter :: u_kmax = 0.01_DP
    integer  :: n_max ! max number of possible momenta1
    integer  :: Neigs = 24 ! Num. of eigenvalues for luscher
    integer  :: L_counter ! Loop variables
    integer  :: mat_size ! n_max*channels
    integer  :: C_3zeros
    real(DP) :: k_max_sqrd ! GeV**2

    real(DP), dimension(:,:), allocatable :: H, H_orig ! Hamiltonian
    real(DP), dimension(:), allocatable   :: omega ! H0 eigenvalues
    integer, dimension(:), allocatable    :: ch_num
    real(DP), dimension(:), allocatable   :: E_int ! eigenvalues
    real(DP), dimension(:), allocatable   :: E_temp
    real(DP), dimension(:), allocatable   :: k_allowed
    integer,  dimension(:), allocatable   :: C_3n!, C_3packed
    integer, dimension(:,:), allocatable  :: bare_index
    integer,  dimension(:), allocatable   :: index_arr

    real(DP) :: V_ij_temp
    real(DP), dimension(1,1) :: V_ij_temp11
    real(DP) :: g_i_temp

    real(DP) :: t_init, t_mid, t_end

    logical :: readBareSlopesFromFile = .false.
    logical :: testMode = .false.
    logical :: saveResult = .true.
    logical :: debugMode = .false.
    integer :: bareSlopeFitOrder
    integer :: nActiveBareFits
    logical, dimension(:), allocatable :: activeBareFits

    real(DP), dimension(:,:), allocatable :: Efit

    ! bobyqa stuff (bq) (turn off limitFit)
    integer :: bq_npt
    real(DP) :: bq_rhobeg, bq_rhoend
    integer :: bq_iprint
    integer :: bq_maxfun
    logical :: useBOBYQA = .true.

    ! logical for disabling some points
    logical, dimension(5) :: enabledMasses = .true.


    call initialiseHEFT()
    call initialiseHEFTFinite()
    call printCurrentParameters(iParamChoice)
    allocate( bare_index(n_bare,2) )

    if (n_bare.eq.0) then
        if (IamRoot) write(*,*) 'No bare state present, stopping...'
        stop
    end if

    if (IamRoot) then
        allocate(activeBareFits(n_bare))
        activeBareFits = .true.
        bareSlopeFitOrder = bareSlopeOrder

        nActiveBareFits = count(activeBareFits)
        allocate(alpha_bare(nActiveBareFits,bareSlopeFitOrder) &
            & , bareFitError(bareSlopeFitOrder*nActiveBareFits) &
            & , alpha_bare_lin(bareSlopeFitOrder*nActiveBareFits) &
            & , slope_low(bareSlopeFitOrder*nActiveBareFits) &
            & , slope_high(bareSlopeFitOrder*nActiveBareFits))

        open(101, file='slopeFitMasses.data', action='read')
        read(101,*) dummyStr, nMasses
        read(101,*) dummyStr, L_m_pi
        read(101,*)

        allocate( m_pi_lat(nMasses,n_bare) &
            & , m_Bar_lat(nMasses,n_bare), dm_Bar_lat(nMasses,n_bare))

        do i = 1, n_bare
           read(101,*)
           do iMasses = 1, nMasses
              read(101,*) testVar, m_pi_lat(iMasses,i) &
                  & , m_Bar_lat(iMasses,i) &
                  & , dm_Bar_lat(iMasses,i)
           end do
        end do
        close(101)

        dof = nMasses - bareSlopeFitOrder*n_bare
        allocate(Efit(nMasses,n_bare))

        ! --------------------Setup Hamiltonian variables-------------------
        L = L_m_pi
        Lambda_max = maxval( (/maxval( Lambda(:,:), dim=2 ), Lambda_v(:)/) )
        k_max_sqrd = inverse_u_k(u_kmax, Lambda_max)**2
        n_max = k_max_sqrd * (L/2.0_DP/pi/hbar_c)**2
        n_max = maxval( (/ n_max, n_min_forced /) )

        allocate( C_3n(0:n_max) )
        call count_3_integers(C_3n, n_max)

        ! Get rid of all the zero states
        C_3zeros = count(C_3n(:) .eq. 0)
        ! n_k = n_max - C_3zeros + 1 - n_init_k

        ! n_k = n_max - C_3zeros
        ! allocate( C_3packed(0:n_k) )

        ! Because the index ik in hamiltonian.f90 always
        ! begins at 1, for S-wave when we have k=0, need
        ! to start C_3packed at 1 to access k=0.
        ! For P-wave, need to start at 0 to ignore k=0.
        ! Hacky and needs to be done better
        if (ch_pWave(ch_onshell) .eq. 'S') then
            n_k = n_max - C_3zeros + 1
            allocate( C_3packed(n_k) )
        else
            n_k = n_max - C_3zeros
            allocate( C_3packed(0:n_k) )
        end if

        C_3packed(:) = pack(C_3n(:), C_3n(:) .ne. 0)
        n_mesh = n_k*n_ch
        N_H = n_bare + n_mesh

        allocate( H(N_H,N_H), omega(N_H), E_temp(N_H) )
        allocate( index_arr(N_H), k_allowed(n_k) )
        allocate( E_int(N_H) )

        ! Excludes momentum when C_3(n) = 0, e.g. n=7
        k_allowed(:) = 2.0_DP*pi/L * sqrt( real(pack( (/ (i,i=n_init_k,n_max) /), &
            & C_3n(n_init_k:n_max) .ne. 0), DP) ) * hbar_c

        write(*,'(a,f0.2,a)') 'Bare mass slope fitting,   L = ', L_m_pi, ' fm'

        ! Set initial guess
        do ii = 1, nActiveBareFits
           alpha_bare(ii,:) = slp_bare_default
        end do
        alpha_bare_lin(:) = reshape(transpose(alpha_bare), [size(alpha_bare_lin)])

        do ii = 1, nActiveBareFits
           write(*,'(a,i1,a3,2f10.3)') 'Initial mass slope guess: alpha_' &
               & , ii, ' = ', alpha_bare(ii,:)
        end do
        ! write(*,'(a,2f0.3)') 'Initial mass slope guess: alpha = ', alpha_bare(1,:)

        if (testMode) then
            write(*,*)
            write(*,*) 'Test Mode'
            call fitBareSlope(1, alpha_bare_lin(:), bare_chi2)

        else
            if (useBOBYQA) then
                ! Use bobyqa for fitting with bounds
                bq_npt = 2*size(alpha_bare_lin(:)) + 1
                bq_rhobeg = 1.0d-4
                bq_rhoend = 1.0d-7
                bq_iprint = 2
                bq_maxfun = 100000
                if (bareSlopeFitOrder.eq.1) then
                    slope_low(:) = 0.1_DP
                    slope_high(:) = 1.6_DP
                else if (bareSlopeFitOrder.eq.2) then
                    if (n_bare.eq.1) then
                        slope_low(:)  = [0.1_DP, -1.0_DP]
                        slope_high(:) = [1.6_DP, 1.0_DP]
                    else if (n_bare.eq.2) then
                        slope_low(:) = [0.7_DP, -0.5_DP, 0.6_DP, -0.5_DP]
                        slope_high(:) = [1.2_DP, 0.0_DP, 1.2_DP, 0.0_DP]
                    end if
                else
                    slope_low(:) = -1.6_DP
                    slope_high(:) = 1.6_DP
                end if
                bare_chi2 = HUGE(1.0d0)

                call bobyqa(size(alpha_bare_lin), bq_npt, alpha_bare_lin &
                    & , slope_low, slope_high, bq_rhobeg, bq_rhoend, bq_iprint &
                    & , bq_maxfun, fitBareSlope_bq)
            else
                ! Use minfun for fitting (unbounded)
                bareFitError = 0.0001_DP
                doBareFitPrinting = .true.
                maxFitIters = 1000
                call Minimize( fitBareSlope &
                    & , size(alpha_bare_lin), alpha_bare_lin(:) &
                    & , bareFitError(:), bareFitStepSize, bare_chi2 &
                    & , doBareFitPrinting )
            end if

            write(*,*)
            do ii = 1, nActiveBareFits
               write(*,'(a,i1,a3,2f10.3)') 'Fitted mass slope: alpha_' &
                   & , ii, ' = ', alpha_bare(ii,:)
            end do
            write(*,'(a,f0.3)') 'chi2 = ', bare_chi2
            write(*,'(a,f0.3)') 'chi2/dof = ', bare_chi2 / real(dof,DP)
        end if


        deallocate(m_pi_lat, m_Bar_lat, dm_Bar_lat)
        deallocate(Efit)

        if (.not.testMode .and. saveResult) then
            ! Write slopes for this fit to file
            open(101, file='data/bare_mass_'//to_string(bareSlopeFitOrder)//'slope_' &
                & //trim(fileNameSuffix)//'.fit', action='write')
            write(101,'(a)') ' Bare Mass Slopes,  ' // fileNameSuffix
            do ii = 1, n_bare
               write(101,'(1x,a,10f10.3)') bare_labels(ii), alpha_bare(ii,:)
            end do
            write(101,'(a,i4)')   ' dof     ', dof
            write(101,'(a,f0.3)') ' chi2         ', bare_chi2
            write(101,'(a,f0.3)') ' chi2/dof     ', bare_chi2 / real(dof,DP)
            close(101)
        end if
        deallocate(alpha_bare, bareFitError)
    end if
    call finaliseHEFT()



contains


    subroutine count_3_integers(C_3, n_max)
        implicit none
        integer, intent(in) :: n_max
        integer, dimension(0:n_max), intent(inout) :: C_3

        integer :: nx, ny, nz, i, sqrt_i

        C_3(:) = 0
        ! C_3(0) = 4
        do i = 0, n_max
           sqrt_i = int( sqrt(real(i,DP)) )
           do nx = -sqrt_i, sqrt_i
              do ny = -sqrt_i, sqrt_i
                 do nz = -sqrt_i, sqrt_i
                    if (nx**2 + ny**2 + nz**2 .eq. i) then
                        C_3(i) = C_3(i) + 1
                    end if
                 end do
              end do
           end do
        end do

    end subroutine count_3_integers


    function fliplr(array)
        implicit none
        real(DP), dimension(:), intent(in) :: array
        real(DP), dimension(size(array)) :: fliplr

        fliplr = array(size(array):1:-1)
    end function fliplr


    subroutine fitBareSlope_bq(n, x, f)
        implicit none
        integer, intent(in) :: n
        real(DP), dimension(:), intent(in) :: x
        real(DP), intent(out) :: f ! chi^2

        call fitBareSlope(n, x, f)

        if (f.lt.bare_chi2) then
            bare_chi2 = f
        end if
    end subroutine fitBareSlope_bq


    subroutine fitBareSlope(n, x, f)
        implicit none
        integer, intent(in) :: n ! n = 1
        real(DP), dimension(n), intent(in) :: x ! alpha_Bare
        real(DP), intent(out) :: f ! chi2

        integer :: iM, iB, xIndexStart, bareCounter, iOrder

        ! integer, dimension() ::

        ! Unpack the x array into alpha_bare
        do iB = 1, nActiveBareFits
           ! how the fuck does this indexing work??
           ! xIndexStart = nActiveBareFits*(iB-1)+1
           ! alpha_bare(iB,:) = x(xIndexStart:xIndexStart+bareSlopeFitOrder)

           xIndexStart = bareSlopeFitOrder*(iB-1) + 1
           if (testMode) write(*,*) 'xIndexStart', xIndexStart
           alpha_bare(iB,:) = x(xIndexStart:xIndexStart+bareSlopeFitOrder)
        end do

        ! Only want to change the slopes if they're active for fitting
        bareCounter = 0
        do iB = 1, n_bare
           if (activeBareFits(iB)) then
               bareCounter = bareCounter+1
               slp_bare(iB,:) = alpha_bare(bareCounter,:)
           end if
        end do

        do iM = 1, nMasses
           m_pi = m_pi_lat(iM,1)
           m_pi2 = m_pi**2

           do jj = 1,n_ch
              m_mes(jj) = sqrt( m_mes_phys(jj)**2 + slp_mes(jj)*(m_pi2 - m_pi0**2) )
              m_bar(jj) = m_bar_phys(jj) + slp_bar(jj)*(m_pi2 - m_pi0**2)
           end do

           m_bare(:) = m_bare_phys(:)
           do iOrder = 1, bareSlopeFitOrder
              m_bare(:) = m_bare(:) + slp_bare(:,iOrder) &
                  &  * (m_pi**(2*iOrder) - m_pi0**(2*iOrder))
           end do

           call generateHamiltonian(H, omega, k_allowed, L)

           if (iM.eq.1 .and. testMode) then
               write(*,*) '---------------Initial Hamiltonian---------------'
               do i = 1,7
                  write(*,'(10f7.3)') H(i,1:7)
               end do
               write(*,*) '-------------------------------------------------'
               write(*,*)
           end if

           call syevd( H(:,:),  E_int(:), jobz='V' )

           ! Info about the states dominated by the bare state
           do ii = 1, n_bare
              bare_index(ii,1) = maxloc(H(ii,:)**2, 1)
              bare_index(ii,2) = maxloc(H(ii,:)**2 &
                  & , 1, (index_arr(:) .ne. bare_index(ii,1)))


              Efit(iM,ii) = E_int(bare_index(ii,1))
              ! write(*,*) iM, ii, Efit(iM,ii)
              ! write(*,*) bare_index(ii,1)
           end do


           ! ! cheeky swaparoo
           ! tempVar = Efit(4,1)
           ! Efit(4,1) = Efit(4,2)
           ! Efit(4,2) = tempVar

           ! force which eigenvalue

           ! Efit(iM) = E_int(1)

           H(:,:) = 0.0_DP
           omega(:) = 0.0_DP
        end do

        f = 0.0_DP
        do ii = 1, n_bare
           ! f = f + sum( (Efit(:,ii) - m_Bar_lat(:,ii))**2 / dm_Bar_lat(:,ii)**2 )

           ! enabledMasses = [1, 1, 1, 1, 1]
           ! enabledMasses = [0, 1, 1, 1, 1]
           ! enabledMasses = [0, 1, 0, 1, 1]
           ! enabledMasses = [0, 1, 1, 0, 1]
           ! enabledMassyes = [0, 0, 1, 1, 1]
           ! enabledMasses = [0, 0, 1, 0, 1]
           enabledMasses = .true.

           f = f + sum( (Efit(:,ii) - m_Bar_lat(:,ii))**2 / dm_Bar_lat(:,ii)**2 &
               & , mask=enabledMasses )

        end do

        if (testMode .or. debugMode) then
            write(*,*) 'x:'
            write(*,*) alpha_bare_lin
            write(*,*) 'x unpacked:'
            do iB = 1, n_bare
               write(*,*) iB, alpha_bare(iB, :)
            end do
            if (debugMode) write(*,*) 'chi2', f
        end if

        if (testMode) then
            do ii = 1, n_bare
               write(*,*) 'Bare', ii
               write(*,*) ' iM     m_pi        E_data      dE_data      E_test'
               do iM = 1, nMasses
                  write(*,'(i4,4f12.5)') iM, m_pi_lat(iM,1), m_Bar_lat(iM,ii) &
                      & , dm_Bar_lat(iM,ii), Efit(iM,ii)
               end do
            end do
            write(*,'(a,f0.2)') 'chi2 = ', f
        end if

    end subroutine fitBareSlope


end program fitBareMassSlope
