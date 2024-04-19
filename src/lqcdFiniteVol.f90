! TODO:
! 1. k_allowed should have an index for channel, currently
!    if I have both S and D wave contributions there's gonna
!    be a k=0 state in D wave which should be forbidden
! 2. Correlation function stuff needs to be generalised better,
!    really should have some higher rank object for alpha instead
!    of alpha1, beta2, etc. Should maybe even be in its own program
program lqcdFiniteVol
    use kinds
    use numFort
    use heft
    use hamiltonian
    implicit none

    real(DP) :: testVar

    ! Pion mass variation
    real(DP) :: m_pi2[*]
    ! real(DP) :: m_pi2_max = m_pi0**2
    real(DP) :: m_pi2_max = 0.4_DP
    real(DP) :: m_pi2_increment
    integer  :: i_m_pi2, m_pi2_counter_nImg
    integer  :: m_pi2_counter_end
    real(DP), dimension(:), allocatable :: m_pi2_array

    logical :: useStored = .true.
    ! real(DP), dimension(:), allocatable :: m_bare ! bare mass
    character (len=8) :: dummyString
    character (len=9) :: wrt_fmt, wrt_fmtInt

    ! Box size variation
    real(DP) :: L[*]
    real(DP), dimension(:), allocatable :: L_array

    ! real(DP) :: L0 = 2.9_DP, Lmax = 10.0_DP
    ! real(DP) :: L, L_increment ! fm
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
    integer :: nSpacings, iSpacing
    real(DP), parameter :: bareFitStepSize = 10.0_DP
    real(DP), dimension(1) :: alpha_bare
    real(DP), dimension(:), allocatable :: bareFitError
    real(DP), dimension(:), allocatable :: a_fm
    real(DP), dimension(:), allocatable :: m_pi_lat
    real(DP), dimension(:), allocatable :: dm_pi_lat
    real(DP), dimension(:), allocatable :: m_Bar_lat
    real(DP), dimension(:), allocatable :: dm_Bar_lat
    real(DP), dimension(:), allocatable :: m_Bare_fit
    real(DP), dimension(:), allocatable :: bare_chi2
    logical :: doBareFitPrinting = .false.

    integer  :: n_max ! max number of possible momenta
    integer  :: Neigs = 24 ! Num. of eigenvalues to output
    integer  :: C_3zeros
    real(DP) :: k_max_sqrd ! GeV**2

    real(DP), dimension(:,:), allocatable :: H[:], H_orig ! Hamiltonian
    real(DP), dimension(:), allocatable   :: omega[:] ! H0 eigenvalues
    integer, dimension(:), allocatable    :: ch_num
    real(DP), dimension(:), allocatable   :: E_int[:] ! eigenvalues
    real(DP), dimension(:), allocatable   :: E_temp
    real(DP), dimension(:,:), allocatable :: k_allowed
    integer,  dimension(:), allocatable   :: C_3n!, C_3packed
    integer, dimension(:,:), allocatable  :: bare_index[:]
    integer,  dimension(:), allocatable   :: index_arr

    real(DP) :: V_ij_temp
    real(DP), dimension(1,1) :: V_ij_temp11
    real(DP) :: g_i_temp

    real(DP) :: t_init, t_mid, t_end

    logical  :: doUnphysicalPion = .true.
    real(DP) :: m_pi2_unphys = 0.0001_DP
    real(DP) :: m_pi2_min
    integer  :: n_inc_unphys, init_mpi2_counter
    integer  :: mpi2_counter_offset

    logical :: slopeFileExists
    logical :: verboseFileRead = .false.
    logical :: readBareSlopesFromFile = .false.
    logical :: readMassSlopesManually = .false.

    integer :: nKappa, nSites
    real(DP) :: m2_low, m2_high
    real(DP), dimension(:), allocatable :: lqcd_a, lqcd_m_pi, lqcd_L

    ! correlation stuff
    logical :: calculateCorrelations = .true.
    logical :: useHEFTEigenvectors = .true.
    logical :: removeSecondEigenvalue = .false.
    logical :: contaminationRedLine = .true.

    real(DP) :: secondEvecThres = 0.01_DP ! |<B_i|E_j>|^2
    integer, parameter :: nt = 64 ! number of time steps
    real(DP) :: t_min = 0.0_DP, t_max, dt
    real(DP) :: alpha1, alpha2, beta1, beta2, scaling_G ! 2 bare states
    real(DP), dimension(:), allocatable :: t[:]
    real(DP), dimension(:), allocatable :: Gt[:] ! correlation
    real(DP), dimension(:,:), allocatable :: Gt_multi[:] ! correlation (multiple)
    real(DP), dimension(:), allocatable :: Ct[:] ! contamination
    real(DP), dimension(:,:), allocatable :: Ct_multi[:] ! contamination (multiple)
    real(DP), dimension(:,:), allocatable :: alpha_interp, beta_interp

    integer :: file_correlation = 151
    integer :: file_contamination = 152

    ! --------------------------Initialise HEFT-------------------------
    call initialiseHEFT()
    call initialiseHEFTFinite()
    call printCurrentParameters(iParamChoice)
    allocate( bare_index(n_bare,3)[*] )

    ! Check the number of points is a multiple of the number of images
    if (mod(m_pi_points,nImg).ne.0) then
        if (IamRoot) then
            write(*,*) 'nPoints = ', m_pi_points
            write(*,*) 'nImages = ', nImg
            write(*,*) 'Make sure the number of points is a multiple'
            write(*,*) " of the number of images you're running it with"
        end if
        stop
    end if
    ! ------------------------------------------------------------------


    ! ------------------Get slope(s) for the bare state-----------------
    if (IamRoot) then
        write(*,'(a,f5.2,a)') ' L_m_pi: ', L_m_pi, ' fm'
        if (readMassSlopesManually) then
            write(*,*) 'Enter mass slopes manually:'
            write(*,*) 'Bare slope order = ', bareSlopeOrder
            do i = 1, n_bare
               write(*,'(a,i1,a)') 'Enter slope(s) for bare state ', i, ': '
               read(*,*) slp_bare(i,:)
               write(*,*) 'You entered: '
               write(*,*) slp_bare(i,:)
            end do
            write(*,*)

        else if (readBareSlopesFromFile) then
            ! For doing gridSearch
            slp_bare(:,:) = 0.0_DP
            if (verboseFileRead) then
                write(*,*) 'Reading mass slopes from bareMassSlopes.in'
                write(*,*) "If you're not doing a grid search, "
                write(*,*) "you probably want to disable readBareSlopesFromFile"
            end if

            open(101, file='bareMassSlopes.in', action='read')
            do i = 1, n_bare
               read(101,*) slp_bare(i,:)
               if (verboseFileRead) write(*,*) 'Slope ', i, '=', slp_bare(1,:)
            end do
            close(101)

        else if (doMassSlopeFitting .and. .not.(readBareSlopesFromFile)) then
            ! If the bare slope fitting program has been run for this fit
            inquire(file = 'data/bare_mass_'//trim(int2str(bareSlopeOrder))//'slope_' &
                & //trim(fileNameSuffix)//'.fit', exist=slopeFileExists)

            if (slopeFileExists) then
                open(101, file='data/bare_mass_'//trim(int2str(bareSlopeOrder))//'slope_' &
                    & //trim(fileNameSuffix)//'.fit', action='read')
                read(101,*)
                do i = 1, n_bare
                   read(101,*) dummyStr, slp_bare(i,:)
                   write(*,'(a,2f10.3)') ' Setting alpha_' // trim(dummyStr) &
                       & // ' = ', slp_bare(i,:)
                end do
                close(101)

            else
                write(*,*) 'Slope data file for ' // trim(fileNameSuffix) &
                    & // ' does not exist. Please run fitBare.x'
                write(*,*) 'Setting slope to default value of', slp_bare_default
                slp_bare(:,:) = 0.0_DP
                do i = 1, n_bare
                   slp_bare(i,:) = slp_bare_default(:)
                end do
            end if

        else
            write(*,'(a,5f8.3)') 'Setting slope(s) to default value of', slp_bare_default
            do i = 1, n_bare
               slp_bare(i,:) = slp_bare_default(:)
            end do
        end if
        write(*,*)
    end if
    if (doCA) call co_broadcast(slp_bare, 1)
    ! ------------------------------------------------------------------


    ! ------------------------Initialise File I/o-----------------------
    if (IamRoot) then
        open(101, file='data/H_eigenvalues_lqcd_' &
            & //trim(fileNameSuffix)//'.out', action='write')
        open(102, file='data/H0_eigenvalues_lqcd_' &
            & //trim(fileNameSuffix)//'.out', action='write')
        open(105, file='data/H_eigenvectors_lqcd_' &
            & //trim(fileNameSuffix)//'.out', action='write')
        open(103, file='data/bare_state_lqcd_' &
            & //trim(fileNameSuffix)//'.out', action='write')
        open(104, file='data/finiteParams_lqcd_' &
            & //trim(fileNameSuffix)//'.out', action='write')
        open(142, file='data/initialHamiltonian_lqcd_' &
            & //trim(fileNameSuffix)//'.out', action='write')
        write(wrt_fmt,'(a1,i2,a6)') '(', Neigs+1, 'f14.6)'

        if (calculateCorrelations) then
            open(file_correlation, file='data/correlation_lqcd_' &
                & //trim(fileNameSuffix)//'.out', action='write')
            open(file_contamination, file='data/contamination_lqcd_' &
                & //trim(fileNameSuffix)//'.out', action='write')
        end if


        ! quick file to tell the plot script if lqcd or mpi was run
        open(143, file='data/latest_finiteProgram_' &
            & //trim(fileNameSuffix)//'.out', action='write')
        write(143,*) 'lqcd'
        close(143)

        ! Write the channel numbers to the file
        write(wrt_fmtInt,'(a1,i2,a6)') '(', Neigs+1, 'i14.1)'
        allocate(ch_num(Neigs))
        do ii = 1,Neigs
           if (ii.le.n_bare) then
               ch_num(ii) = 0
           else
               ch_num(ii) = mod(ii-n_bare-1,n_ch)+1
           end if
        end do
        write(102,wrt_fmtInt) 0, ch_num(:Neigs)
        deallocate(ch_num)
    end if
    ! ------------------------------------------------------------------

    ! Find maximum Lambda to use for Hamiltonian cutoff at each L
    Lambda_max = maxval( (/maxval( Lambda(:,:), dim=2 ), Lambda_v(:)/) )
    k_max_sqrd = inverse_u_k(u_kmax,Lambda_max)**2

    ! array of quark masses
    m_pi2_max = m_pi_max**2
    m_pi2_min = m_pi_min**2
    allocate(m_pi2_array(m_pi_points))
    m_pi2_array = linspace(m_pi2_min, m_pi2_max, m_pi_points)

    ! array of lattice sizes
    allocate(L_array(m_pi_points))

    ! Use lQCD eigenvectors for correlation functions
    !  - This code is just for the odd-parity nucleons
    !    TODO: generalise/input this from a data file
    if (.not. useHEFTEigenvectors) then
        ! alpha beta interpolants
        ! --structure and flow of the nucleon eigenstates in lattice qcd
        ! alpha 1: 0.7 - 0.8
        ! beta 1: 0.25 - 0.55
        ! alpha 2: 0.6 - 0.8
        ! beta 2: -0.6 - -0.45
        if (n_bare.gt.0) then
            allocate( alpha_interp(n_bare,m_pi_points), beta_interp(n_bare,m_pi_points) )
            alpha_interp(1,:) = real([(i,i=0,m_pi_points-1)],DP) &
                & * (0.8_DP-0.7_DP) / real(m_pi_points-1,DP) + 0.7_DP
            beta_interp(1,:) = real([(i,i=0,m_pi_points-1)],DP) &
                & * (0.55_DP-0.25_DP) / real(m_pi_points-1,DP) + 0.25

            if (n_bare.eq.2) then
                alpha_interp(2,:) = real([(i,i=0,m_pi_points-1)],DP) &
                    & * (0.8_DP-0.6_DP) / real(m_pi_points-1,DP) + 0.6_DP
                beta_interp(2,:) = real([(i,i=0,m_pi_points-1)],DP) &
                    & * (-0.45_DP + 0.6_DP) / real(m_pi_points-1,DP) - 0.6_DP
            end if
        end if
    end if

    ! Need to correctly vary the lattice extent L as a function of m_pi**2
    !   TODO: need to generalise this a bit for various lattice sizes
    if (IamRoot) then
        open(167, file=fileName_lQCD, action='read')
        read(167,*) dummyStr, nKappa
        read(167,*)
        read(167,*) dummyStr, nSites
        read(167,*)
        allocate(lqcd_a(nKappa), lqcd_m_pi(nKappa), lqcd_L(nKappa))

        do i = 1, nKappa
           read(167,*) lqcd_a(i), lqcd_m_pi(i)
        end do
        close(167)
        lqcd_L(:) = lqcd_a * real(nSites,DP) ! convert spacings to lattice extents

        ! for points between lQCD points, vary L linearly
        do i = 1, nKappa-1
           m2_low = lqcd_m_pi(i)**2
           m2_high = lqcd_m_pi(i+1)**2

           where((m_pi2_array(:).ge.m2_low) .and. (m_pi2_array(:).lt.m2_high))
               L_array(:) = (lqcd_L(i+1) - lqcd_L(i)) / (m2_high - m2_low) &
                   & * (m_pi2_array(:) - m2_low) + lqcd_L(i)
           end where
        end do

        ! for points outside the bounds of the lQCD points, keep lattice size constant
        where(m_pi2_array(:) .lt. lqcd_m_pi(1)**2) L_array(:) = lqcd_L(1)
        where(m_pi2_array(:) .ge. lqcd_m_pi(nKappa)**2) L_array(:) = lqcd_L(nKappa)

    end if

    call co_broadcast(L_array, 1)

    allocate(t(nt)[*], Ct(nt)[*], Gt(nt)[*])

    sync all
    if (IamRoot) call cpu_time(t_init)

    ! --------------------Main loop over pion masses--------------------
    massloop: do m_pi2_counter_nImg = 1, m_pi_points, nImg
       i_m_pi2 = m_pi2_counter_nImg + (iImg - 1)

       m_pi2 = m_pi2_array(i_m_pi2)
       m_pi = sqrt(m_pi2)
       L = L_array(i_m_pi2)

       ! Alter the hadron masses as a function of quark mass
       do jj = 1,n_ch
          m_mes(jj) = sqrt( m_mes_phys(jj)**2 + slp_mes(jj)*(m_pi2 - m_pi0**2) )
          m_bar(jj) = m_bar_phys(jj) + slp_bar(jj)*(m_pi2 - m_pi0**2)
       end do

       ! Bare mass extrapolation up to Order(m_{\pi}^{2*bareSlopeOrder})
       if (n_bare.gt.0) then
           m_bare(:) = m_bare_phys(:)
           do jj = 1, bareSlopeOrder
              m_bare(:) = m_bare(:) + slp_bare(:,jj) &
                  &  * (m_pi**(2*jj) - m_pi0**(2*jj))
           end do
       end if

       ! Get the size of the Hamiltonian
       ! To make general plot scripts a bit easier, force it to have
       !    some minimum size
       n_max = k_max_sqrd * (L/2.0_DP/pi/hbar_c)**2
       n_max = maxval( [n_max, n_min_forced] )

       allocate( C_3n(0:n_max) )
       call count_3_integers(C_3n, n_max)

       ! Get rid of all the zero states
       C_3zeros = count(C_3n(:) .eq. 0)

       ! n_k = n_max - C_3zeros
       ! allocate( C_3packed(0:n_k) )

       ! Because the index ik in hamiltonian.f90 always
       ! begins at 1, for S-wave when we have k=0, need
       ! to start C_3packed at 1 to access k=0.
       ! For P-wave, need to start at 0 to ignore k=0.
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

       allocate( H(N_H,N_H)[*], omega(N_H)[*], E_temp(N_H) )
       allocate( index_arr(N_H), k_allowed(n_k,n_ch) )
       allocate( E_int(N_H)[*] )

       ! Excludes momentum when C_3(n) = 0, e.g. n=7
       do i = 1, n_ch
          k_allowed(:,i) = 2.0_DP*pi/L * sqrt(real(pack([(i,i=n_init_k,n_max-(1-n_init_k))] &
              & , C_3n(n_init_k:(n_max-(1-n_init_k))) .ne. 0), DP) ) * hbar_c
       end do
       index_arr(:) = [ (i,i=1,N_H) ]


       ! Generate the hamiltonian, see hamiltonian.f90
       call generateHamiltonian(H, omega, k_allowed, L)

       if (i_m_pi2 .eq. 1) then
           write(*,*) '---------------Initial Hamiltonian---------------'
           do i = 1,7
              write(*,'(10f7.3)') H(i,1:7)
           end do
           write(*,*) '-------------------------------------------------'
           write(*,*)
           write(*,*) '        -------------------------------------------'
           write(*,*) '       |        n       |    m_pi^2   |     L      |'
           write(*,*) '       |----------------|-------------|------------|'

           ! -------------------Write the initial H to a file------------------
           do i = 1, N_H
              write(142,'('//trim(int2str(N_H))//'f12.6)') H(i, :)
           end do
       end if

       ! Solves for eigenvalues of H
       call syevd( H(:,:),  E_int(:), jobz='V' )

       ! Info about the 3 states with largest bare state(s) contributions
       do ii = 1, n_bare
          bare_index(ii,1) = maxloc(H(ii,:)**2, 1)
          bare_index(ii,2) = maxloc(H(ii,:)**2 &
              & , 1, (index_arr(:) .ne. bare_index(ii,1)))
          bare_index(ii,3) = maxloc(H(ii,:)**2 &
              & , 1, (index_arr(:) .ne. bare_index(ii,1)) &
              & .and. (index_arr(:) .ne. bare_index(ii,2)))
       end do

       ! calculate correlation function
       if (calculateCorrelations) then
           t_max = L ! isotropic
           dt = (t_max - t_min) / real(nt-1,DP) + t_min
           t(:) = [(real([(i,i=0,nt-1)], DP) * dt)]

           Gt(:) = 0.0_DP
           Ct(:) = 0.0_DP
           if (n_bare.eq.1) then
               do i = 1, N_H
                  Gt(:) = Gt(:) + abs(H(1,i))**2 &
                      & * exp(-E_int(i) * t(:) / GeVfm)
                  if (i.ne.bare_index(1,1)) then
                      Ct(:) = Ct(:) + abs(H(1,i))**2 &
                          & * exp(-E_int(i) * t(:) / GeVfm)
                  end if
               end do
               Ct(:) = Ct(:) / Gt(:)

           else if (n_bare.eq.2) then
               if (useHEFTEigenvectors) then
                   alpha1 = H(1, bare_index(1,1))
                   beta1 = H(2, bare_index(1,1))

                   alpha2 = H(1, bare_index(2,1))
                   beta2 = H(2, bare_index(2,1))
               else
                   alpha1 = alpha_interp(1,i_m_pi2)
                   beta1 = beta_interp(1,i_m_pi2)

                   alpha2 = alpha_interp(2,i_m_pi2)
                   beta2 = beta_interp(2,i_m_pi2)
               end if

               ! scaling_G = 1.0_DP / (alpha1**2 + beta1**2)
               ! alpha1 = alpha1 * sqrt(scaling_G)
               ! beta1 = beta1 * sqrt(scaling_G)

               ! scaling_G = 1.0_DP / (alpha2**2 + beta2**2)
               ! alpha2 = alpha2 * sqrt(scaling_G)
               ! beta2 = beta2 * sqrt(scaling_G)

               write(*,*)
               write(*,'(5f6.2,f8.3)') m_pi2, alpha1, H(1,bare_index(1,1)), beta1, H(2,bare_index(1,1)), alpha1*H(1,bare_index(1,1)) &
                   & + beta1*H(2,bare_index(1,1))
               write(*,'(5f6.2,f8.3)') m_pi2, alpha1, H(1,bare_index(2,1)), beta1, H(2,bare_index(2,1)), alpha1*H(1,bare_index(2,1)) &
                   & + beta1*H(2,bare_index(2,1))

               ! TODO: dear god this needs a cleanup and rewrite
               do i = 1, N_H
                  if (contaminationRedLine) then
                      ! first contamination function (red)
                      Gt(:) = Gt(:) + (alpha1*H(1,i) + beta1*H(2,i))**2 &
                          & * exp(-E_int(i)*t(:) / GeVfm)

                      Ct(:) = Gt(:)
                      do bare_i = 1, n_bare
                         ! remove contribution from bare dominated state
                         ! for the second bare state, check it hasn't already been removed
                         ! i.e. one eigenstate has largest contribution from both
                         !      bare states
                         if (.not.((bare_i.eq.2) .and. (bare_index(1,1).eq.bare_index(2,1)))) then
                             Ct(:) = Ct(:) - (alpha1*H(1,bare_index(bare_i,1)) &
                                 & + beta1*H(2,bare_index(bare_i,1)))**2 &
                                 & * exp(-E_int(bare_index(bare_i,1))*t(:) / GeVfm)
                         end if

                         ! check if the contribution from the eigenvalue with
                         ! 2nd largest bare_i component is greater than secondEvecThres%
                         ! also check that the state with 2nd highest hasnt already been removed
                         ! if so, remove it from the correlation function
                         if (removeSecondEigenvalue .and. &
                             & (H(bare_i,bare_index(bare_i,2))**2).ge.secondEvecThres &
                             & .and. (bare_index(bare_i,2).ne.bare_index(mod(bare_i,2)+1,1))) then

                             ! also make sure that the state with second largest hasn't already
                             ! been removed when you look at n_bare=2
                             if (.not.((bare_i.eq.2) &
                                 & .and. (bare_index(2,2).eq.bare_index(1,2)))) then
                                 Ct(:) = Ct(:) - (alpha1*H(1,bare_index(bare_i,2)) &
                                     & + beta1*H(2,bare_index(bare_i,2)))**2 &
                                     & * exp(-E_int(bare_index(bare_i,2))*t(:) / GeVfm)
                             end if
                         end if
                      end do

                  else
                      ! second contamination function (blue)
                      Gt(:) = Gt(:) + (alpha2*H(1,i) + beta2*H(2,i))**2 &
                          & * exp(-E_int(i)*t(:) / GeVfm)
                      ! don't add the states dominated by a bare state
                      ! if ((i.ne.bare_index(1,1)) .and. (i.ne.bare_index(2,1)) ) then
                      !     Ct(:) = Ct(:) + (alpha2*H(1,i) + beta2*H(2,i))**2 &
                      !         & * exp(-E_int(i)*t(:) / GeVfm)
                      ! end if

                      Ct(:) = Gt(:)
                      do bare_i = 1, n_bare
                         ! remove contribution from bare dominated state
                         if (.not.((bare_i.eq.2) .and. (bare_index(1,1).eq.bare_index(2,1)))) then
                             Ct(:) = Ct(:) - (alpha2*H(1,bare_index(bare_i,1)) &
                                 & + beta2*H(2,bare_index(bare_i,1)))**2 &
                                 & * exp(-E_int(bare_index(bare_i,1))*t(:) / GeVfm)
                         end if

                         if (bare_i.eq.1) then
                             write(*,*) alpha2*H(1,bare_index(bare_i,1)) &
                                 & + beta2*H(2,bare_index(bare_i,1))
                         end if

                         ! check if the contribution from the eigenvalue with
                         ! 2nd largest bare_i component is greater than secondEvecThres%
                         ! also check that the state with 2nd highest hasnt already been removed
                         ! if so, remove it from the correlation function
                         if (removeSecondEigenvalue .and. &
                             & (H(bare_i,bare_index(bare_i,2))**2).ge.secondEvecThres &
                             & .and. (bare_index(bare_i,2).ne.bare_index(mod(bare_i,2)+1,1))) then
                             ! write(*,*) m_pi2, bare_i, (H(bare_i,bare_index(bare_i,1))**2)
                             Ct(:) = Ct(:) - (alpha2*H(1,bare_index(bare_i,2)) &
                                 & + beta2*H(2,bare_index(bare_i,2)))**2 &
                                 & * exp(-E_int(bare_index(bare_i,2))*t(:) / GeVfm)
                         end if
                      end do
                  end if
               end do
               Ct(:) = Ct(:) / Gt(:)
           end if
       end if

       sync all

       ! ------------------------Write data to files-----------------------
       if (IamRoot) then
           do jImg = 1, nImg
              E_temp = E_int(:)[jImg]
              ! H eigenvalues
              write(101,wrt_fmt) m_pi2[jImg], E_temp(:Neigs)
              E_temp = omega(:)[jImg]
              ! H0 eigenvalues
              write(102,wrt_fmt) m_pi2[jImg], E_temp(:Neigs)

              ! Eigenvectors
              do ii = 1, 11
                 write(105,'(12f14.9)') m_pi2[jImg], abs(H(ii,:8)[jImg])**2
              end do

              do ii = 1, n_bare
                 write(103,'(f8.4,3f12.6,3i4)') m_pi2[jImg] &
                     & , E_int(bare_index(ii,1)[jImg])[jImg] &
                     & , E_int(bare_index(ii,2)[jImg])[jImg] &
                     & , E_int(bare_index(ii,3)[jImg])[jImg] &
                     & , bare_index(ii,1)[jImg] &
                     & , bare_index(ii,2)[jImg] &
                     & , bare_index(ii,3)[jImg]
              end do

              ! write correlation stuff
              if (calculateCorrelations) then
                  write(file_correlation,'(101f12.6)') m_pi2, t(:)[jImg]
                  write(file_correlation,'(101f12.6)') L, Gt(:)[jImg]
                  write(file_contamination,'(101f12.6)') m_pi2, t(:)[jImg]
                  write(file_contamination,'(101f12.6)') L, Ct(:)[jImg]
              end if
           end do
       end if



       ! Write initial/zeroth/final point to terminal
       if ((i_m_pi2.eq.1) &
           & .or. (i_m_pi2.eq.m_pi_points)) then
           write(*,'(8x,a1,2x,i4,a2,i4,4x,a1,f9.3,4x,a1,f8.2,4x,a1)') &
               & '|', i_m_pi2, '/', m_pi_points, '|', m_pi2 &
               & , '|', L, '|'
       end if

       ! Write progress to terminal
       if (iImg.eq.nImg .and. i_m_pi2.gt.1) then
           write(*,'(8x,a1,2x,i4,a2,i4,4x,a1,f9.3,4x,a1,f8.2,4x,a1,a,\)') &
               & '|', i_m_pi2, '/', m_pi_points, '|' &
               & , m_pi2, '|', L, '|', char(13)
       end if
       if (mod(i_m_pi2,m_pi_points/10) == 0) then
           write(*,'(8x,a1,2x,i4,a2,i4,4x,a1,f9.3,4x,a1,f8.2,4x,a1)') &
               & '|', i_m_pi2, '/', m_pi_points, '|', m_pi2 &
               & , '|', L, '|'
       end if

       deallocate(C_3n, H, omega, E_temp, index_arr &
           & , k_allowed, E_int, C_3packed)
       sync all
    end do massloop


    if (IamRoot) then
        call cpu_time(t_end)
        write(*,*) '        ------------------------------------------- '
        write(*,*)
        write(*,*) 'Time taken: ' // timePrint(t_end-t_init)
        write(*,*)
        write(*,*)

        write(104,'(4a14)') 'n_ch', 'n_bare', 'L', 'Lambda_max'
        write(104,'(2i14,2f14.7)') n_ch, n_bare, L, Lambda_max
    end if

    deallocate(m_pi2_array, L_array, bare_index)
    if (IamRoot) then
        deallocate(lqcd_a, lqcd_m_pi, lqcd_L)
    end if
    call finaliseHEFT()


    close(101)
    close(102)
    close(103)
    close(104)
    close(105)
    !---------------------------------------------------------------------!
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

end program lqcdFiniteVol
