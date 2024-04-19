program mpiFiniteVol
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

    ! real(DP), parameter :: u_kmax = 0.01_DP
    integer  :: n_max ! max number of possible momenta1
    integer  :: Neigs = 24 ! Num. of eigenvalues for luscher
    integer  :: L_counter ! Loop variables
    integer  :: mat_size ! n_max*channels
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

    logical  :: doUnphysicalPion = .false.
    real(DP) :: m_pi2_unphys = 0.0001_DP
    real(DP) :: m_pi2_min
    integer  :: n_inc_unphys, init_mpi2_counter
    integer  :: mpi2_counter_offset

    logical :: slopeFileExists
    logical :: verboseFileRead = .false.
    logical :: readBareSlopesFromFile = .false.
    logical :: readMassSlopesManually = .false.

    ! correlation function stuff
    real(DP), dimension(:), allocatable   :: t ! Time (fm)
    real(DP), dimension(:), allocatable   :: Gt[:], Ct[:] ! Correl., Contamin. Fns

    call initialiseHEFT()
    call initialiseHEFTFinite()
    call printCurrentParameters(iParamChoice)
    allocate( bare_index(n_bare,3)[*] )

    if (IamRoot) then
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
               read(101,*) slp_bare(i,1)
               if (verboseFileRead) write(*,*) 'Slope ', i, '=', slp_bare(1,1)
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
    end if

    if (doCA) call co_broadcast(slp_bare, 1)

    ! ------------------------Initialise File I/o-----------------------
    if (IamRoot) then
        open(101, file='data/H_eigenvalues_m_pi_' &
            & //trim(fileNameSuffix)//'.out', action='write')
        open(102, file='data/H0_eigenvalues_m_pi_' &
            & //trim(fileNameSuffix)//'.out', action='write')
        open(105, file='data/H_eigenvectors_m_pi_' &
            & //trim(fileNameSuffix)//'.out', action='write')
        open(106, file='data/bare_state_eigenvectors_m_pi_' &
            & //trim(fileNameSuffix)//'.out', action='write')
        open(103, file='data/bare_state_m_pi_' &
            & //trim(fileNameSuffix)//'.out', action='write')
        open(104, file='data/finiteParams_m_pi_' &
            & //trim(fileNameSuffix)//'.out', action='write')
        open(142, file='data/initialHamiltonian_' &
            & //trim(fileNameSuffix)//'.out', action='write')
        write(wrt_fmt,'(a1,i2,a6)') '(', Neigs+1, 'f14.6)'

        ! quick file to tell the plot script if lqcd or mpi was run
        open(143, file='data/latest_finiteProgram_' &
            & //trim(fileNameSuffix)//'.out', action='write')
        write(143,*) 'mpi'
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

    Lambda_max = maxval( [maxval( Lambda(:,:), dim=2 ), Lambda_v(:)] )
    L = L_m_pi

    m_pi2_max = m_pi_max**2
    m_pi2_min = m_pi_min**2
    m_pi2_increment = (m_pi2_max - m_pi2_min)/real(m_pi_points-1,DP)

    ! -----------------Smaller than physical pion stuff-----------------
    if (doUnphysicalPion) then
        init_mpi2_counter = -int( (m_pi0**2 - m_pi2_unphys) / m_pi2_increment )
        mpi2_counter_offset = -init_mpi2_counter + 1
        m_pi2_min = m_pi2_unphys
    else
        init_mpi2_counter = 1
        m_pi2_min = m_pi2_min
        ! m_pi2_min = m_pi0**2
        mpi2_counter_offset = 1
    end if
    m_pi2_counter_end = init_mpi2_counter + m_pi_points - 1

    ! ---------------------Initial Hamiltonian Setup--------------------
    k_max_sqrd = inverse_u_k(u_kmax, Lambda_max)**2
    n_max = k_max_sqrd * (L/2.0_DP/pi/hbar_c)**2 + 1
    n_max = maxval( (/ n_max, n_min_forced /) )

    allocate( C_3n(0:n_max) )
    call count_3_integers(C_3n, n_max)

    ! Get rid of all the zero states
    C_3zeros = count(C_3n(:) .eq. 0)

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

    allocate( H(N_H,N_H)[*], omega(N_H)[*], E_temp(N_H) )
    allocate( index_arr(N_H), k_allowed(n_k,n_ch) )
    allocate( E_int(N_H)[*] )

    ! Excludes momentum when C_3(n) = 0, e.g. n=7
    do i = 1, n_ch
       k_allowed(:,i) = 2.0_DP*pi/L * sqrt(real(pack([(i,i=n_init_k,n_max-(1-n_init_k))] &
           & , C_3n(n_init_k:(n_max-(1-n_init_k))) .ne. 0), DP) ) * hbar_c
    end do
    index_arr(:) = [ (i,i=1,N_H) ]

    if (IamRoot) then
        write(*,'(a10,f10.3,a10,f6.2,a10,i8)') 'm_pi0: ', m_pi0, 'L: ', L, 'n_H: ', n_H
        write(*,*)
    end if


    ! ------------------------------------------------------------------
    ! -------------------------Do pion mass loop------------------------
    ! ------------------------------------------------------------------
    m_pi2 = m_pi2_min
    sync all
    if (IamRoot) call cpu_time(t_init)

    massloop : do m_pi2_counter_nImg = init_mpi2_counter, m_pi2_counter_end, nImg

       m_pi2_counter = m_pi2_counter_nImg + iImg - 1

       ! need -1 (unphys) vs -2 (phys)
       ! please fix at some point
       if (doUnphysicalPion) then
           m_pi2 = m_pi2_min + (m_pi2_counter-1+mpi2_counter_offset) * m_pi2_increment
       else
           m_pi2 = m_pi2_min + (m_pi2_counter-2+mpi2_counter_offset) * m_pi2_increment
       end if
       m_pi = sqrt(m_pi2)

       do jj = 1,n_ch
          m_mes(jj) = sqrt( m_mes_phys(jj)**2 + slp_mes(jj)*(m_pi2 - m_pi0**2) )
          m_bar(jj) = m_bar_phys(jj) + slp_bar(jj)*(m_pi2 - m_pi0**2)
       end do

       if (n_bare.gt.0) then
           m_bare(:) = m_bare_phys(:)
           do jj = 1, bareSlopeOrder
              m_bare(:) = m_bare(:) + slp_bare(:,jj) &
                  &  * (m_pi**(2*jj) - m_pi0**(2*jj))
           end do
       end if

       ! Generate the hamiltonian, see hamiltonian.f90
       call generateHamiltonian(H, omega, k_allowed, L)

       if (m_pi2_counter .eq. init_mpi2_counter) then
           write(*,*) '---------------Initial Hamiltonian---------------'
           do i = 1,7
              write(*,'(10f7.3)') H(i,1:7)
           end do
           write(*,*) '-------------------------------------------------'
           write(*,*)
           write(*,*) '        ------------------------------'
           write(*,'(8x,a1,2x,a10,4x,a1,1x,a9,3x,a1)') &
               & '|', '   n   ', '|', 'm_pi^2', '|'
           write(*,*) '       |----------------|-------------|'

           ! -------------------Write the initial H to a file------------------
           do i = 1, N_H
              write(142,'('//trim(int2str(N_H))//'f12.6)') H(i, :)
           end do
       end if

       ! Solves for eigenvalues of H
       call syevd( H(:,:),  E_int(:), jobz='V' )

       ! Info about the states dominated by the bare state
       do ii = 1, n_bare
          bare_index(ii,1) = maxloc(H(ii,:)**2, 1)
          bare_index(ii,2) = maxloc(H(ii,:)**2 &
              & , 1, (index_arr(:) .ne. bare_index(ii,1)))
          bare_index(ii,3) = maxloc(H(ii,:)**2 &
              & , 1, (index_arr(:) .ne. bare_index(ii,1)) &
              & .and. (index_arr(:) .ne. bare_index(ii,2)))
       end do

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

              do ii = 1, n_bare
                 write(103,'(f8.4,3f12.6,3i4)') m_pi2[jimg] &
                     & , E_int(bare_index(ii,1)[jImg])[jImg] &
                     & , E_int(bare_index(ii,2)[jImg])[jImg] &
                     & , E_int(bare_index(ii,3)[jImg])[jImg] &
                     & , bare_index(ii,1)[jImg] &
                     & , bare_index(ii,2)[jImg] &
                     & , bare_index(ii,3)[jImg]
              end do

              ! Write eigenvectors
              do ii = 1, 11
                 write(105,'(12f14.9)') m_pi2[jImg], abs(H(ii,:8)[jImg])**2
              end do

              ! write eigenvector of the bare dominated value
              if (n_bare.eq.1) then
                  write(106, '(f8.4,16f14.9)') m_pi2[jImg] &
                      & , abs(H(:16,bare_index(1,1)[jImg])[jImg]**2)
              end if
           end do
       end if

       ! Write initial/zeroth/final point to terminal
       if (m_pi2_counter.eq.init_mpi2_counter &
           & .or. (m_pi2_counter.eq.m_pi2_counter_end)&
           & .or. (m_pi2_counter.eq.0)) then
           write(*,'(8x,a1,2x,i4,a2,i4,4x,a1,f9.2,4x,a1)') &
               & '|', m_pi2_counter, '/', m_pi2_counter_end, '|', m_pi2 &
               & , '|'
       end if

       ! Write progress to terminal
       if (iImg.eq.nImg .and. m_pi2_counter.gt.1) then
           write(*,'(8x,a1,2x,i4,a2,i4,4x,a1,f9.2,4x,a1,a,\)') &
               & '|', m_pi2_counter, '/', m_pi2_counter_end, '|' &
               & , m_pi2, '|', char(13)
       end if
       if (mod(m_pi2_counter,m_pi_points/10) == 0) then
           write(*,'(8x,a1,2x,i4,a2,i4,4x,a1,f9.2,4x,a1)') &
               & '|', m_pi2_counter, '/', m_pi2_counter_end, '|', m_pi2 &
               & , '|'
       end if
       sync all

    end do massloop

    if (IamRoot) then
        call cpu_time(t_end)
        write(*,*) '        ------------------------------ '
        write(*,*)
        write(*,*) 'Time taken: ' // timePrint(t_end-t_init)
        write(*,*)
        write(*,*)

        write(104,'(4a14)') 'n_ch', 'n_bare', 'L', 'Lambda_max'
        write(104,'(2i14,2f14.7)') n_ch, n_bare, L, Lambda_max
    end if

    deallocate( H, omega, k_allowed, E_int &
        & , C_3n, C_3packed, index_arr, bare_index )
    deallocate(slp_bare)
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


    function fliplr(array)
        implicit none
        real(DP), dimension(:), intent(in) :: array
        real(DP), dimension(size(array)) :: fliplr

        fliplr = array(size(array):1:-1)
    end function fliplr


    subroutine fitBareSlope(n, x, f)
        ! Subroutine for fitting the mass slope of the bare mass
        !    as a function of the pion mass (from LQCD)
        implicit none
        integer, intent(in) :: n ! n = 1a
        real(DP), dimension(n), intent(in) :: x ! alpha_Bare
        real(DP), intent(out) :: f ! chi2

        real(DP) :: alpha ! slope of the bare mass
        real(DP) :: yint

        alpha = x(1)
        ! Fix the curve such that it passes through the physical point correctly
        yint = m_bare(iBare_fit) - alpha*m_pi0**2
        m_Bare_fit(:) = alpha*m_pi_lat**2 + yint

        ! 2D chi2
        ! https://www.astro.rug.nl/software/kapteyn/kmpfittutorial.html#fitting-data-when-both-variables-have-uncertainties

        f = sum( (m_Bar_lat(:) - m_Bare_fit(:))**2 &
            & / (dm_pi_lat(:)**2 + dm_Bar_lat(:)**2*alpha**2) )

        ! Only fit to the last few
        f = sum( (m_Bar_lat(2:) - m_Bare_fit(2:))**2 &
            & / (dm_pi_lat(2:)**2 + dm_Bar_lat(2:)**2*alpha**2) )

    end subroutine fitBareSlope


end program mpiFiniteVol
