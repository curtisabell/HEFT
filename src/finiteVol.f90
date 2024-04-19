program finiteVol
    use kinds
    use NumFort
    use heft
    use hamiltonian
    implicit none

    real(DP) :: testVar

    ! Pion mass variation
    real(DP) :: m_pi2
    real(DP) :: m_pi2_increment
    integer  :: m_pi2_counter

    ! Box size variation
    real(DP) :: L[*]
    real(DP) :: L_increment
    integer  :: L_counter, L_counter_nImg

    logical :: useStored = .true.
    logical :: doPrinting = .true.
    character(len=8) :: dummyString
    character(len=9) :: wrt_fmt, wrt_fmtInt
    character(len=1) :: n_digits ! digits of n_max for H_fmt
    character(:), allocatable :: H_fmt

    integer :: n_k ! number of momentum states
    integer :: n_mesh ! n_ch * n_k
    integer :: N_H ! full Hamiltonian size
    integer :: ch_i, ch_j, bare_i, bare_j
    integer :: ik, ich, ii, jk, jj, i_bare, j_bare ! Loop indicies
    real(DP) :: Lambda_max

    integer  :: n_max ! max number of possible momenta
    integer  :: n_max_forced = 1024 ! forced max matrix size
    integer  :: n_min_forced = 128 ! forced min matrix size
    integer  :: Neigs = 40 ! Num. of eigenvalues for luscher

    integer  :: i, j, jImg ! Loop variables
    integer  :: mat_size ! n_max*channels
    integer  :: C_3zeros
    real(DP) :: k_max_sqrd ! GeV**2
    real(DP) :: m_pi_diff_tol = 1.0d-6
    logical  :: slopeFileExists
    integer, dimension(:,:), allocatable :: bare_index[:]

    real(DP), dimension(:,:), allocatable :: H[:], H_orig ! Hamiltonian
    real(DP), dimension(:), allocatable   :: omega[:] ! H0 eigenvalues
    integer , dimension(:), allocatable   :: ch_num
    real(DP), dimension(:), allocatable   :: E_int[:] ! eigenvalues
    real(DP), dimension(:), allocatable   :: E_temp
    real(DP), dimension(:,:), allocatable :: k_allowed
    integer,  dimension(:), allocatable   :: C_3n !, C_3packed
    integer,  dimension(:), allocatable   :: index_arr

    real(DP) :: V_ij_temp
    real(DP), dimension(1,1) :: V_ij_temp11
    real(DP) :: g_i_temp

    ! File labels
    integer :: file_eigenvalues  = 101
    integer :: file_basisEnergy  = 102
    integer :: file_bareStates   = 103
    integer :: file_finiteParams = 104

    call initialiseHEFT('fin')
    call initialiseHEFTFinite()
    call printCurrentParameters(iParamChoice)

    allocate( bare_index(n_bare,3)[*] )
    E_gl = m_mes_phys(ch_onshell) + m_bar_phys(ch_onshell)

    ! m_pi = 0.388_DP
    ! if (m_pi.ne.m_pi0) then
    if (abs(m_pi - m_pi0) .gt. m_pi_diff_tol) then
        if (IamRoot) write(*,*) 'Unphysical pion mass, m_pi = ', m_pi
        do jj = 1,n_ch
           m_mes(jj) = sqrt( m_mes_phys(jj)**2 + slp_mes(jj)*(m_pi**2 - m_pi0**2) )
           m_bar(jj) = m_bar_phys(jj) + slp_bar(jj)*(m_pi**2 - m_pi0**2)
        end do

        ! Read in bare slope from file/use default
        if (IamRoot) then
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
                do i = 1, n_bare
                   slp_bare(i,:) = [slp_bare_default, 0.0_DP]
                end do
            end if
        end if
        if (doCA) call co_broadcast(slp_bare, 1)

        if (n_bare.gt.0) then
            m_bare(:) = m_bare_phys(:)
            do jj = 1, bareSlopeOrder
               m_bare(:) = m_bare(:) + slp_bare(:,jj) &
                   &  * (m_pi**(2*jj) - m_pi0**(2*jj))
            end do
        end if
    end if


    ! ------------------------Initialise File I/O-----------------------
    if (IamRoot) then
        open(file_eigenvalues, file='data/H_eigenvalues_L_' &
            & //trim(fileNameSuffix)//'.out', action='write')
        open(file_basisEnergy, file='data/H0_eigenvalues_L_' &
            & //trim(fileNameSuffix)//'.out', action='write')
        open(file_bareStates, file='data/bare_state_L_' &
            & //trim(fileNameSuffix)//'.out', action='write')
        open(file_finiteParams, file='data/finiteParams_L_' &
            & //trim(fileNameSuffix)//'.out', action='write')
        open(108,file='data/bare_state_eigenvectors2_L_' &
            & //trim(fileNameSuffix)//'.out', action='write')
        open(109,file='data/H_eigenvectors2_L_' &
            & //trim(fileNameSuffix)//'.out', action='write')

        ! Format for writing eigenvalues to file
        write(wrt_fmt,'(a1,i2,a6)') '(', Neigs+1, 'f11.6)'

        ! Write the channel numbers to the file
        write(wrt_fmtInt,'(a1,i2,a6)') '(', Neigs+1, 'i11.6)'
        allocate(ch_num(Neigs))
        do ii = 1,Neigs
           if (ii.le.n_bare) then
               ch_num(ii) = 0
           else
               ch_num(ii) = mod(ii-n_bare-1,n_ch)+1
           end if
        end do
        write(file_basisEnergy,wrt_fmtInt) 135, ch_num(:Neigs)
        deallocate(ch_num)
    end if

    ! Find maximum Lambda to use for Hamiltonian cutoff at each L
    Lambda_max = maxval( (/maxval( Lambda(:,:), dim=2 ), Lambda_v(:)/) )
    ! L = L_min
    L_increment = (L_max - L_min)/real(L_points-1,DP)
    k_max_sqrd = inverse_u_k(u_kmax,Lambda_max)**2

    ! ------------------------------------------------------------------
    ! ----------------------Do lattice volume loop----------------------
    ! ------------------------------------------------------------------
    do L_counter_nImg = 1, L_points, nImg
       L_counter = L_counter_nImg + iImg - 1
       L = L_min + (L_counter-1) * L_increment

       ! n_max = int(k_max_sqrd * (L/2.0_DP/pi/hbar_c)**2) + n_bare

       ! Add 1 since this allows us to subtract 1 in the case
       !    where we have an S-wave channel with k=0, giving
       !    the same number of momenta in k_allowed for each channel
       n_max = int(k_max_sqrd * (L/2.0_DP/pi/hbar_c)**2) + 1
       n_max = maxval( (/ n_max, n_min_forced /) )

       ! Need to reduce the size of the matrix Hamiltonian by the number of times
       !    C_3n == 0. This is to avoid the contributions of "false" basis states, such as n=7
       !    where there are no ways to sum three square integers to produce 7 etc.
       allocate( C_3n(0:n_max) )
       call count_3_integers(C_3n, n_max)
       C_3zeros = count(C_3n(:) .eq. 0)

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
       n_mesh = (1 - n_init_k + n_k)*n_ch
       N_H = n_bare + n_mesh

       ! -----------------------Generate Hamiltonian-----------------------

       allocate( H(N_H,N_H)[*] )
       allocate( E_int(N_H)[*],  omega(N_H)[*], E_temp(N_H) )
       allocate( index_arr(N_H), k_allowed(n_k,n_ch) )

       ! Excludes momentum when C_3(n) = 0, e.g. n=7
       do i = 1, n_ch
          k_allowed(:,i) = 2.0_DP*pi/L * sqrt(real(pack([(i,i=n_init_k,n_max-(1-n_init_k))] &
              & , C_3n(n_init_k:(n_max-(1-n_init_k))) .ne. 0), DP) ) * hbar_c
       end do
       index_arr(:) = [ (i,i=1,N_H) ]

       ! See hamiltonian.f90
       call generateHamiltonian(H, omega, k_allowed, L)
       ! Print the top left corner of the first H
       if (L_counter .eq. 1 .and. IamRoot) then
           write(*,*) '---------------Initial Hamiltonian---------------'
           do i = 1,7
              write(*,'(10f7.3)') H(i,1:7)
           end do
           write(*,*) '-------------------------------------------------'
           write(*,*)
           write(*,'(21x,a8,a7)') 'L (fm)', 'n_H'
       end if

       ! -------------------Solve the eigenvalue equation------------------
       call syevd( H(:,:),  E_int(:), jobz='V' )

       if (L_counter.eq.1 .and. IamRoot) then
           write(*,*) 'L: ', L
           ! write(*,'(a7,f8.4,3f18.10,3x)') 'u_min: ', u_kmax, E_int(1), E_int(2), E_int(3)
       end if

       ! Calculate the positions of states with largest bare componants
       do ii = 1, n_bare
          bare_index(ii,1) = maxloc(H(ii,:)**2, 1)
          bare_index(ii,2) = maxloc(H(ii,:)**2 &
              & , 1, (index_arr(:) .ne. bare_index(ii,1)))
          bare_index(ii,3) = maxloc(H(ii,:)**2 &
              & , 1, (index_arr(:) .ne. bare_index(ii,1)) &
              & .and. (index_arr(:) .ne. bare_index(ii,2)))
       end do

       ! Sync images so that the correct info is written to file
       sync all

       ! Using the root image, write all info to files
       if (IamRoot) then
           do jImg = 1, nImg
              ! Write energy eigenvalues to file
              ! For some reason need this E_temp to avoid a weird
              ! crash because of "undefined variables???" ifort err maybe
              E_temp = E_int(:)[jImg]
              write(file_eigenvalues, wrt_fmt) L[jImg], E_temp(:Neigs)
              E_temp = omega(:)[jImg]
              write(file_basisEnergy, wrt_fmt) L[jImg], E_temp(:Neigs)

              ! Write location of states with largest bare componants
              do ii = 1, n_bare
                 write(file_bareStates,'(f8.4,3f9.6,3i4)') L[jImg] &
                     & , E_int(bare_index(ii,1)[jImg])[jImg] &
                     & , E_int(bare_index(ii,2)[jImg])[jImg] &
                     & , E_int(bare_index(ii,3)[jImg])[jImg] &
                     & , bare_index(ii,1)[jImg] &
                     & , bare_index(ii,2)[jImg] &
                     & , bare_index(ii,3)[jImg]
              end do

              ! write eigenvector of the bare dominated value
              if (n_bare.eq.1) then
                  write(108, '(f8.4,16f14.9)') L[jImg] &
                      & , abs(H(:16,bare_index(1,1)[jImg])[jImg]**2)
              end if

              ! write Eigenvectors
              do ii = 1, 16
                 write(109,'(f8.4,11f14.9)') L[jImg], abs(H(ii,:8)[jImg])**2
              end do

           end do

           ! Write first point to terminal
           if (L_counter.eq.1) then
               write(*,'(10x,i4,a2,i4,f8.2,i8)') L_counter, '/', L_points, L, N_H
           end if
       end if


       ! Write progress to terminal
       sync all
       if (iImg.eq.nImg .and. L_counter.gt.1) then
           write(*,'(10x,i4,a2,i4,f8.2,i8,a,\)') L_counter, '/', L_points, L, N_H, char(13)
           if ((mod(L_counter,L_points/maxval((/4,nImg/))) == 0) &
               & .or. (L_counter.eq.1)) then
               write(*,'(10x,i4,a2,i4,f8.2,i8)') L_counter, '/', L_points, L, N_H
           end if
       end if

       deallocate( H, omega, k_allowed, E_int &
           & , C_3n, C_3packed, index_arr, E_temp)
    end do

    if (IamRoot) then
        write(*,*)
        write(file_finiteParams,*) Lambda_max
        close(file_eigenvalues)
        close(file_basisEnergy)
        close(file_bareStates)
        close(file_finiteParams)
        close(108)
        close(109)
    end if

    call finaliseHEFT()
    deallocate(bare_index)

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

end program finiteVol
