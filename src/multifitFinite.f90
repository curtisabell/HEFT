program mpiFiniteVol
    use kinds
    use numFort
    use heft
    use hamiltonian
    use, intrinsic :: iso_fortran_env, only : iostat_end
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
    character (len=9) :: wrt_fmtInt
    character (len=13) :: wrt_fmt

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
    integer :: initial_N_H[*]
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
    real(DP), dimension(:), allocatable   :: H_row[:] ! temp array
    real(DP), dimension(:), allocatable   :: E_temp
    real(DP), dimension(:), allocatable   :: k_allowed
    integer,  dimension(:), allocatable   :: C_3n!, C_3packed
    integer, dimension(:,:), allocatable  :: bare_index[:]
    integer,  dimension(:), allocatable   :: index_arr

    real(DP) :: V_ij_temp
    real(DP), dimension(1,1) :: V_ij_temp11
    real(DP) :: g_i_temp

    real(DP) :: t_init, t_mid, t_end
    real(DP) :: m_pi_diff_tol = 1.0d-6

    logical :: slopeFileExists

    integer :: io_error
    logical :: io_stop = .false.
    logical :: paramsFileExists
    integer :: nParamFits ! how many fits to loop over
    integer :: iFit[*], iFit_nImg
    real(DP), dimension(:), allocatable :: thisParam
    real(DP), dimension(:,:), allocatable :: paramArray
    character(len=128) :: multifit_params_file, multifit_params_file_full

    call initialiseHEFT()
    call initialiseHEFTFinite()
    call printCurrentParameters(iParamChoice)
    allocate( bare_index(n_bare,3)[*] )
    m_pi = m_pi0
    m_pi2 = m_pi**2

    ! -----------------------Unphysical pion mass-----------------------
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

    ! ----------------------Read in fits from file----------------------
    if (IamRoot) then
        call get_command_argument(1, multifit_params_file)
        ! check if there is an argument
        if (multifit_params_file.eq.'') then
            write(*,*) 'Enter a param file as an argument'
            write(*,*) 'Stopping...'
            io_stop = .true.

        else
            ! if there is an argument, check it is valid
            multifit_params_file_full = trim(multifit_params_file) // '.params'
            inquire(file=trim(multifit_params_file_full), exist=paramsFileExists)
            if (paramsFileExists) then
                ! invalid if called one of these
                if (trim(multifit_params_file).eq.'allFits') then
                    write(*,*) 'allFits is not a valid input'
                    write(*,*) 'Stopping...'
                    io_stop = .true.
                else if (trim(multifit_params_file).eq.'singleFit') then
                    write(*,*) 'singleFit is not a valid input'
                    write(*,*) 'Stopping...'
                    io_stop = .true.
                end if
            else
                ! Invalid if the file doesn't even exist
                write(*,*) 'That .params file does not exist'
                write(*,*) 'Stopping...'
                io_stop = .true.
            end if
        end if
    end if

    ! Stop the program if no params file was given as an argument
    sync all
    if (doCA) call co_broadcast(io_stop, 1)
    if (io_stop) stop

    ! Read params from the file now we know its a valid input
    if (IamRoot) then
        write(*,*) 'Params file:  ' // trim(multifit_params_file_full)
        open(192, file=trim(multifit_params_file_full), action='read')
        read(192,*)
        read(192,*) dummyString, nParamFits
    end if

    allocate(paramArray(nParamFits,nParamsTotal))

    if (IamRoot) then
        do ii = 1, nParamFits
           read(192,*) paramArray(ii,:)
        end do
        close(192)
    end if

    sync all
    if (doCA) call co_broadcast(paramArray, 1)

    ! ------------------------Initialise File I/o-----------------------
    if (IamRoot) then
        open(101, file='data/H_eigenvalues_multifit_' &
            & //trim(multifit_params_file)//'.out', action='write')
        open(102, file='data/H0_eigenvalues_multifit_' &
            & //trim(multifit_params_file)//'.out', action='write')
        open(105, file='data/H_eigenvectors_multifit_' &
            & //trim(multifit_params_file)//'.out', action='write')
        open(103, file='data/bare_state_multifit_' &
            & //trim(multifit_params_file)//'.out', action='write')
        write(wrt_fmt,'(a5,i2,a6)') '(i14,', Neigs, 'f14.6)'

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

    L = L_multifit
    allocate(E_temp(n_min_forced), thisParam(nParamsTotal))

    if (IamRoot) then
        write(*,'(a,i0)') 'nFits: ', nParamFits
        write(*,'(a,f8.3)') 'm_pi: ', m_pi
        write(*,'(a,f8.3)') 'L (fm): ', L
        write(*,*)
        write(*,*) '        ------------------------------'
    end if

    ! ------------------------------------------------------------------
    ! -------------------------Do loop over fits------------------------
    ! ------------------------------------------------------------------
    sync all
    fitLoop : do iFit_nImg = 1, nParamFits, nImg
       iFit = iFit_nImg + iImg - 1
       thisParam(:) = paramArray(iFit,:)
       call unpackFitParams(thisParam)

       ! do m_pi2_counter_nImg = init_mpi2_counter, m_pi2_counter_end, nImg
       ! m_pi2_counter = m_pi2_counter_nImg + iImg - 1

       Lambda_max = maxval( [maxval( Lambda(:,:), dim=2 ), Lambda_v(:)] )
       k_max_sqrd = inverse_u_k(u_kmax, Lambda_max)**2
       n_max = k_max_sqrd * (L/2.0_DP/pi/hbar_c)**2
       n_max = maxval( (/ n_max, n_min_forced /) )

       allocate( C_3n(0:n_max) )
       call count_3_integers(C_3n, n_max)

       ! Get rid of all the zero states
       C_3zeros = count(C_3n(:) .eq. 0)
       ! n_k = n_max - C_3zeros + 1 - n_init_k
       n_k = n_max - C_3zeros
       allocate( C_3packed(0:n_k) )
       C_3packed(:) = pack(C_3n(:), C_3n(:) .ne. 0)
       n_mesh = n_k*n_ch
       N_H = n_bare + n_mesh

       ! this is used for writing the eigenvectors to file
       ! needs to be the same every time so set it to the first one
       ! if (iFit==1) then
           ! initial_N_H = N_H
           ! call co_broadcast(initial_N_H, iImg)
       ! end if

       allocate( H(N_H,N_H)[*], omega(N_H)[*] )
       allocate( index_arr(N_H), k_allowed(n_k) )
       allocate( E_int(N_H)[*], H_row(N_H)[*] )

       ! Excludes momentum when C_3(n) = 0, e.g. n=7
       k_allowed(:) = 2.0_DP*pi/L * sqrt( real(pack( (/ (i,i=n_init_k,n_max) /), &
           & C_3n(n_init_k:n_max) .ne. 0), DP) ) * hbar_c

       index_arr(:) = [(i,i=1,N_H)]

       sync all
       if (IamRoot) call cpu_time(t_init)

       ! Generate the hamiltonian, see hamiltonian.f90
       call generateHamiltonian(H, omega, k_allowed, L)

       ! Solves for eigenvalues/eigenvectors of H
       call syevd( H(:,:),  E_int(:), jobz='V' )

       ! Get location of states dominated by the bare state(s)
       do ii = 1, n_bare
          H_row(:) = H(ii,:)**2
          bare_index(ii,1) = maxloc(H_row, 1)
          H_row(bare_index(ii,1)) = 0.0_DP
          bare_index(ii,2) = maxloc(H_row, 1)
          H_row(bare_index(ii,2)) = 0.0_DP
          bare_index(ii,3) = maxloc(H_row, 1)
          ! write(*,*) H_row(:4)
          ! write(*,*) bare_index(ii,:)
          ! stop
       end do

       ! ! Get location of states dominated by the bare state(s)
       ! do ii = 1, n_bare
       !    bare_index(ii,1) = maxloc(H(ii,:)**2, 1)
       !    bare_index(ii,2) = maxloc(H(ii,:)**2 &
       !        & , 1, (index_arr(:) .ne. bare_index(ii,1)))
       !    bare_index(ii,3) = maxloc(H(ii,:)**2 &
       !        & , 1, (index_arr(:) .ne. bare_index(ii,1)) &
       !        & .and. (index_arr(:) .ne. bare_index(ii,2)))
       ! end do

       sync all
       ! ------------------------Write data to files-----------------------
       if (IamRoot) then
           do jImg = 1, nImg
              ! H eigenvalues
              E_temp(:) = E_int(:n_min_forced)[jImg]
              write(101,wrt_fmt) iFit[jImg], E_temp(:Neigs)

              ! H0 eigenvalues
              E_temp = omega(:n_min_forced)[jImg]
              write(102,wrt_fmt) iFit[jImg], E_temp(:Neigs)


              ! if (Lambda(1,1).gt.7.9_DP) then
              !     write(*,*) N_H
              !     do ii = 1, 11
              !        write(*,*) ii, H(1,ii)**2
              !     end do
              ! end if

              ! print all the Eigenvectors
              do ii = 1, Neigs
                 write(105,'(i8,11f14.9)') iFit[jImg], H(ii,:11)[jImg]**2
              end do

              ! ! Eigenvectors, and sum the final one
              ! do ii = 1, 4
              !    write(105,'(i8,11f14.9)') iFit[jImg], H(ii,:11)[jImg]**2
              ! end do
              ! write(105,'(i8,11f14.9)') iFit[jImg], sum(H(5:,:11)[jImg]**2, dim=1)

              do ii = 1, n_bare
                 write(103,'(i8,3f12.6,3i4)') iFit[jimg] &
                     & , E_int(bare_index(ii,1)[jImg])[jImg] &
                     & , E_int(bare_index(ii,2)[jImg])[jImg] &
                     & , E_int(bare_index(ii,3)[jImg])[jImg] &
                     & , bare_index(ii,1)[jImg] &
                     & , bare_index(ii,2)[jImg] &
                     & , bare_index(ii,3)[jImg]
              end do

              ! do ii = 1, n_bare
              !    write(103,'(i8,3f12.6,3i4)') iFit[jimg] &
              !        & , E_temp(bare_index(ii,1)[jImg]) &
              !        & , E_temp(bare_index(ii,2)[jImg]) &
              !        & , E_temp(bare_index(ii,3)[jImg]) &
              !        & , bare_index(ii,1)[jImg] &
              !        & , bare_index(ii,2)[jImg] &
              !        & , bare_index(ii,3)[jImg]
              ! end do
           end do
       end if

       ! --------------------Write progress to terminal--------------------
       if ((iFit.eq.1) .or. (iFit.eq.nParamFits) &
           & .or. (mod(iFit,nParamFits/10) == 0)) then
           write(*,'(8x,a1,2x,i4,a2,i4,4x,a1,f9.2,4x,a1,i9,a3)') &
               & '|', iFit, '/', nParamFits, '|', m_pi2 &
               & , '|', N_H, '  |'
       end if

       if (iImg.eq.nImg .and. iFit.gt.1) then
           write(*,'(8x,a1,2x,i4,a2,i4,4x,a1,f9.2,4x,a1,i9,a3,a,\)') &
               & '|', iFit, '/', nParamFits, '|' &
               & , m_pi2, '|', N_H, '  |', char(13)
       end if

       ! if (mod(iFit,nParamFits/10) == 0) then
       !     write(*,'(8x,a1,2x,i4,a2,i4,4x,a1,f9.2,4x,a1)') &
       !         & '|', iFit, '/', nParamFits, '|', m_pi2 &
       !         & , '|'
       ! end if

       deallocate( H, omega, k_allowed, E_int, H_row &
           & , C_3n, C_3packed, index_arr )
       sync all
    end do fitloop


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

    deallocate(slp_bare, bare_index)
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


end program mpiFiniteVol
