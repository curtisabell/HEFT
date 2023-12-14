module heft
    use kinds
    use NumFort
    implicit none

    real(DP), parameter :: f_pi = 0.0924_DP ! GeV
    real(DP), parameter :: m_pi0 = 0.1385_DP ! GeV
    real(DP), parameter :: hbar_c = 0.1973_DP ! GeVfm
    real(DP), parameter :: GeVfm  = 0.1973_DP ! hbar * c
    real(DP), parameter :: GeV2_mb = 0.3894_DP ! hbar * c (mb = millibarns)
    real(DP), parameter :: degrees = 180.0_DP / pi ! / rad
    real(DP), parameter :: M_octet = 0.920_DP
    real(DP), parameter :: u_kmax = 0.01_DP ! Limits k_max

    real(DP) :: beta = 0.67_DP

    ! Coarray variables
    integer :: iImg, nImg
    logical :: IamRoot, doCA
    integer, parameter :: root_image = 1

    ! SU6 couplings
    real(DP), parameter :: D_SU6 = 0.76_DP
    real(DP), parameter :: F_SU6 = 0.50_DP
    real(DP) :: C_SU6 = -2.0_DP * D_SU6
    ! real(DP) :: chi_Delta = 3.0_DP/32.0_DP/pi/f_pi**2 * 2.0_DP/9.0_DP * C_SU6**2

    integer  :: n_bare
    integer  :: n_ch
    integer  :: n_f
    integer  :: n_init_k
    real(DP) :: E_gl ! Global energy
    real(DP) :: E_init, E_fin
    integer  :: nPoints_inf
    real(DP) :: m_pi
    logical  :: useDataPoints ! If energy array should be read or generated from E_init, E_fin
    integer,  dimension(:), allocatable   :: C_3packed

    real(DP) :: L_min, L_max
    real(DP) :: L_multifit
    real(DP) :: L_m_pi
    real(DP) :: fit_chi2
    integer  :: L_points

    real(DP) :: m_pi_max
    real(DP) :: m_pi_min
    integer  :: m_pi_points

    ! Ends of each section of the parameter array
    integer :: bare_end, bare_end_file
    integer :: g_end, g_end_file
    integer :: Lam_end, Lam_end_file
    integer :: v_end, v_end_file
    integer :: Lam_v_end, Lam_v_end_file

    ! Fitting variables
    integer :: phaseIndex
    logical :: doFitting
    logical :: doEtaFitting ! inelasticity
    logical :: doCSFitting ! cross section
    logical :: printFittingOutput
    integer :: nParamsTotal
    integer :: nParams ! Only active parameters
    real(DP), dimension(:), allocatable :: fitParams, activeFitParams, fitParamsOrig
    real(DP), dimension(:), allocatable :: fitParamErrors, activeFitParamErrors
    logical,  dimension(:), allocatable :: isActiveParam
    logical,  dimension(:), allocatable :: isVaryingParam
    integer,  dimension(:), allocatable :: fitParamIndex, varyingParamIndex
    real(DP), dimension(:), allocatable :: fitParamLow, fitParamHigh
    real(DP), dimension(:), allocatable :: activeFitParamLow, activeFitParamHigh
    real(DP), dimension(:), allocatable :: fitParamGuess, activeFitParamGuess
    real(DP), dimension(:), allocatable :: varyingParamMin, varyingParamMax

    real(DP), dimension(:), allocatable   :: m_bare ! bare mass
    real(DP), dimension(:), allocatable   :: m_bare_phys
    real(DP), dimension(:), allocatable   :: m_mes, m_bar
    real(DP), dimension(:), allocatable   :: m_mes_phys, m_bar_phys
    real(DP), dimension(:,:), allocatable :: gBare ! bare state -> channel couplings
    real(DP), dimension(:,:), allocatable :: vCh ! channel -> channel couplings
    real(DP), dimension(:,:), allocatable :: Lambda ! bare->channel cutoffs
    real(DP), dimension(:),   allocatable :: Lambda_v ! channel->channel cutoffs
    real(DP), dimension(:,:), allocatable :: id ! identity matrix

    ! Mass slopes
    real(DP), dimension(:,:), allocatable :: slp_bare ! bare mass
    real(DP), dimension(:), allocatable :: slp_mes  ! meson mass
    real(DP), dimension(:), allocatable :: slp_bar  ! baryon mass
    integer :: bareSlopeOrder
    real(DP), dimension(:), allocatable :: slp_bare_default
    ! Mass slope fitting variables
    logical :: doMassSlopeFitting
    character(len=32), dimension(:), allocatable :: massSlopeFileNames ! dim=nBare

    logical :: isEnergyDepPotential
    character(len=1), dimension(:), allocatable :: ch_AM ! angular momentum in each ch
    character(len=9), dimension(:), allocatable :: ch_labels
    character(len=9), dimension(:), allocatable :: bare_labels
    character(len=1), dimension(:), allocatable :: ch_pWave
    integer :: ch_onshell
    integer :: onshell_AM ! onshell angular momentum
    character(len=1) :: onshell_pWave

    character(len=32) :: fitParamsFileName
    character(len=32) :: fileNameSuffix
    character(len=32) :: dummyStr
    character(len=32) :: projName
    character(len=4) :: nbmc
    logical :: useCustomMasses
    character(len=32) :: customMassFile

    integer :: unit_paramOutput = 201
    character(len=128) :: fileName_paramOutput
    character(len=128) :: fileName_particles = '../src/particles.in'
    character(len=128) :: fileName_lQCD = 'slopeFitMasses.data'

    ! Pole search variables
    complex(DP) :: E_pole
    complex(DP), dimension(:), allocatable :: k_pole
    complex(DP), dimension(:), allocatable :: E_poles, E_poles_init
    integer :: nPoles
    integer, dimension(:), allocatable :: poles_bare_loc

    real(DP), dimension(:), allocatable :: pole_rots
    complex(DP) :: k_complex_thetas

    ! Stuff for choosing potentials
    character(len=1) :: f_choice
    character(len=1) :: g_choice
    character(len=1) :: u_choice
    character(len=3) :: pot_choices

    ! The fit number
    integer :: iParamChoice
    logical :: readParamsFromFile

    type particle
        character(len=16) :: name
        logical :: isMeson
        real(DP) :: mass
        real(DP) :: spin
        integer :: sparity
        real(DP) :: slope
    end type particle

    integer :: nParticles
    type(particle), dimension(:), allocatable :: particles
    type(particle), dimension(:,:), allocatable :: ch

    logical :: verbose = .false.

    interface u_k
        module procedure u_k_cmplx, u_k_real
    end interface u_k

    interface f_i
        module procedure f_i_cmplx, f_i_real
    end interface f_i

    interface V_ij
        module procedure V_ij_cmplx, V_ij_real
    end interface V_ij

    interface g_i
        module procedure g_i_cmplx, g_i_real
    end interface g_i

    interface pWaveConvert
        module procedure pWave_to_AM, AM_to_pWave
    end interface pWaveConvert

contains


    subroutine initialiseHEFT(volume)
        implicit none
        character(len=3), optional :: volume ! must be inf or fin
        ! Loop variables
        integer :: i, kk, iPart

        ! IO variables
        integer :: file_particles
        integer :: file_system
        character(len=128) :: fileName
        character(len=128) :: string
        character(len=64), dimension(:), allocatable :: chNameA, chNameB

        ! Default parameter values
        real(DP) :: def_m_bare   = 2.0_DP
        real(DP) :: def_gBare    = 0.1_DP
        real(DP) :: def_Lambda   = 0.8_DP
        real(DP) :: def_vCh      = -0.1_DP
        real(DP) :: def_Lambda_v = 0.8_DP

        ! Variables for reading parameters from allFits.params
        integer :: iFileParams, iParams
        character(len=32) :: file_all_fits = 'allFits.params'
        character(len=16) :: fmt_readFileParams
        real(DP), dimension(:), allocatable :: fileParams

        ! -----------------------Set coarray variables----------------------
        iImg = this_image()
        nImg = num_images()
        IamRoot = (iImg .eq. root_image)
        doCA = (nImg .gt. 1)

        ! -------------------Read scattering system config------------------
        if (IamRoot) then
            file_system = 142
            fileName = 'HEFT.config'
            open(file_system, file=fileName, action='read')
            do i = 1, 2
               read(file_system,*)
            end do
            read(file_system,*) string, n_ch
            read(file_system,*) string, n_bare
            read(file_system,*)

            write(nbmc, '(i1,a1,i1,a1)') n_bare, 'b', n_ch, 'c'

            ! Write details to terminal
            write(*,*) '--------------------HEFT '// nbmc &
                & //'--------------------'
            write(*,*) 'nImages: ' // trim(int2str(nImg))
        end if

        if (doCA) then
            call co_broadcast(n_ch, 1)
            call co_broadcast(n_bare, 1)
        end if

        ! Allocate parameter data
        allocate(m_mes(n_ch), m_bar(n_ch), m_mes_phys(n_ch), m_bar_phys(n_ch) &
            &, gBare(n_ch,n_bare), vCh(n_ch,n_ch), Lambda(n_ch,n_bare) &
            &, Lambda_v(n_ch), m_bare(n_bare), m_bare_phys(n_bare) &
            &, slp_bar(n_ch), slp_mes(n_ch), id(n_ch,n_ch) &
            &, ch_labels(n_ch), bare_labels(n_bare) &
            &, ch_AM(n_ch), ch_pWave(n_ch), chNameA(n_ch), chNameB(n_ch))

        ! Get the scattering channels for this system
        if (IamRoot) then
            read(file_system,*)
            do i = 1, n_ch
               read(file_system,*) chNameA(i), chNameB(i), ch_pWave(i)
            end do

            ! Read which channel number is onshell, and which partial wave
            read(file_system,*)
            read(file_system,*) string, ch_onshell
            onshell_pWave = ch_pWave(ch_onshell)
            ! read(file_system,*) string, onshell_pWave

            ! Logical for using default/custom masses
            read(file_system,*) string, useCustomMasses
            read(file_system,*) customMassFile

            ! If doing infinite volume calcs, always want physical masses
            if (present(volume)) then
                if (volume.eq.'inf') useCustomMasses = .false.
            end if

            ! Read info about the bare state(s)
            if (n_bare .gt. 0) then
                read(file_system,*)
                read(file_system,*)
                do i = 1, n_bare
                   read(file_system,*) bare_labels(i)
                end do
            end if
            read(file_system,*)
            read(file_system,*)

            if (n_bare .eq. 0) then
                write(projName, *) nbmc
            else
                write(projName, *) trim(adjustl(bare_labels(1))), '_', nbmc
            end if

            ! Read which potentials to use
            read(file_system,*) g_choice, f_choice, u_choice
            pot_choices = g_choice//f_choice//u_choice
            close(file_system)
        end if

        if (doCA) then
            call co_broadcast(ch_onshell, 1)
            call co_broadcast(onshell_pWave, 1)
            call co_broadcast(bare_labels, 1)
            call co_broadcast(g_choice, 1)
            call co_broadcast(f_choice, 1)
            call co_broadcast(u_choice, 1)
            call co_broadcast(pot_choices, 1)
        end if

        ! Set n_f
        ! Define the dimension of the potential functions f
        if (f_choice.eq.'W') then
            ! Weinberg-Tamozawa potential
            n_f = 2
        else
            n_f = 1
        end if

        ! -------------------------Read in particles------------------------
        if (IamRoot) then
            file_particles = 151
            if (useCustomMasses) then
                fileName = trim(adjustl(customMassFile))
                if (trim(fileName).eq.'particles.in') then
                    fileName = fileName_particles
                else
                    write(*,*) 'Using custom particle file: ' // trim(fileName)
                end if
            else
                fileName = fileName_particles
                ! fileName = '../src/particles.in'
            end if
            open(file_particles, file=fileName, action='read')
            read(file_particles, *) string, nParticles
        end if

        if (doCA) call co_broadcast(nParticles, 1)

        allocate(particles(nParticles), ch(2,nParticles))
        if (IamRoot) then
            ! Skip blank space
            do i = 1,6
               read(file_particles,*)
            end do

            ! Read particle data
            do i = 1, nParticles
               read(file_particles,*)
               read(file_particles,*) particles(i)%name
               read(file_particles,*) particles(i)%mass
               read(file_particles,*) particles(i)%spin
               read(file_particles,*) particles(i)%sparity
               read(file_particles,*) particles(i)%slope
            end do
            close(file_particles)


            do i = 1, n_ch
               do iPart = 1, nParticles
                  if (chNameA(i).eq.particles(iPart)%name) then
                      ch(1,i) = particles(iPart)
                  end if

                  if (chNameB(i).eq.particles(iPart)%name) then
                      ch(2,i) = particles(iPart)
                  end if
               end do
            end do
        end if

        if (doCA) then
            call co_broadcast(particles, 1)
            call co_broadcast(ch, 1)
        end if

        ! Convert partial wave to angular momentum
        onshell_AM = pWaveConvert(onshell_pWave)
        if (onshell_AM .eq. -1) then
            if (IamRoot) then
                write(*,*)
                write(*,*) 'Enter S, P, D, F, G or H as the onshell partial wave in HEFT.config.'
                write(*,*) 'Exiting...'
            end if
            stop
        end if

        ! Set initial value of allowed k
        if (onshell_pWave.eq.'S' .or. onshell_pWave.eq.'s') then
            n_init_k = 0
        else
            n_init_k = 1
        end if

        ! Print scattering channels and bare states
        if (IamRoot) then
            if (n_bare.gt.0) then
                write(*,*) 'Bare states: '
                do i = 1, n_bare
                   write(*,*) '     ', trim(bare_labels(i))
                end do
            end if
            write(*,*) 'Scattering Channels: '
            do i = 1, n_ch
               write(*,*) '     ', trim(ch(1,i)%name) // '-' // trim(ch(2,i)%name) &
                   &  // ' (' // ch_pWave(i) // '-wave)'
            end do
            write(*,*) '-------------------------------------------------'
            write(*,*)
        end if

        ! Distribute channel info to some global variables
        do i = 1, n_ch
           ch_labels(i)  = trim(ch(1,i)%name) // '-' // trim(ch(2,i)%name)
           m_mes_phys(i) = ch(1,i)%mass
           m_bar_phys(i) = ch(2,i)%mass
           slp_mes(i)    = ch(1,i)%slope
           slp_bar(i)    = ch(2,i)%slope
        end do

        m_bar(:) = m_bar_phys(:)
        m_mes(:) = m_mes_phys(:)

        ! Set some misc variables
        id(:,:) = real(identityMat(n_ch),DP)
        E_gl = m_mes(ch_onshell) + m_bar(ch_onshell)

        ! ------------------Read in parameters from allFits-----------------
        if (IamRoot) then
            readParamsFromFile = .true.
            file_all_fits = 'allFits.params'
            iParamChoice = 1
            if (readParamsFromFile) then
                open(108, file='allFits.params', action='read')
                read(108,*) dummyStr, iFileParams
                read(108,*) dummyStr, iParamChoice
                read(108,*) dummyStr, bare_end_file, g_end_file &
                    & , Lam_end_file, v_end_file, Lam_v_end_file
                read(108,*)

                allocate(fileParams(iFileParams+1))
                do
                   read(108,*) iParams, fileParams
                   if (iParams.eq.iParamChoice) then
                       if (n_bare.gt.0) then
                           m_bare_phys = fileParams(1:bare_end_file)
                           m_bare = m_bare_phys
                           gBare  = vectorToMat(fileParams(bare_end_file+1 : g_end_file) &
                               & , n_ch, n_bare)
                           Lambda = vectorToMat(fileParams(g_end_file+1 : Lam_end_file) &
                               & ,  n_ch, n_bare)
                       end if
                       vCh      = vectorToSymMat(fileParams(Lam_end_file+1 : v_end_file))
                       Lambda_v = fileParams(v_end_file+1 : Lam_v_end_file)

                       fit_chi2 = fileParams(Lam_v_end_file+1)
                       exit

                   else if (iParams.eq.0) then
                       write(*,'(x,a,i0,a)') 'Parameter set ', iParamChoice, ' does not exist.'
                       write(*,*) 'Please choose a correct value for iFileParams.'
                       stop
                   end if
                end do
                deallocate(fileParams)

                fileNameSuffix = 'fit'//int2str(iParamChoice)

                ! Write current parameter set to a file so it can be more easily parsed
                !    by the python plotting scripts
                fileName_paramOutput = 'data/params_' // trim(fileNameSuffix) &
                    & // '.out'
                open(unit_paramOutput, file=trim(fileName_paramOutput), action='write')
                if (n_bare.ge.1) then
                    do i = 1, n_bare
                       write(unit_paramOutput, '(a10, f12.6)') &
                           & 'm_b'//trim(int2str(i)), m_bare(i)
                    end do
                    do i = 1, n_bare
                       write(unit_paramOutput, '(a10, 10f12.6)') &
                           & 'g_b'//trim(int2str(i)), gBare(:,i)
                    end do
                    do i = 1, n_bare
                       write(unit_paramOutput, '(a10, 10f12.6)') &
                           & 'Lam_b'//trim(int2str(i)), Lambda(:,i)
                    end do
                end if
                do i = 1, n_ch
                   write(unit_paramOutput, '(a10, 10f12.6)') &
                       & 'v_c'//trim(int2str(i)), vCh(:,i)
                end do
                write(unit_paramOutput, '(a10, 10f12.6)') 'Lam_v', Lambda_v(:)
                close(unit_paramOutput)

            else
                ! Otherwise set to default values
                m_bare   = def_m_bare
                gBare    = def_gBare
                Lambda   = def_Lambda
                vCh      = def_vCh
                Lambda_v = def_Lambda_v
                fit_chi2 = HUGE(1.0_DP)

                fileNameSuffix = 'fitDefault'
            end if
        end if

        if (doCA) then
            call co_broadcast(ch_pWave,    1)
            call co_broadcast(ch_labels,   1)
            call co_broadcast(m_bare,      1)
            call co_broadcast(m_bare_phys, 1)
            call co_broadcast(gBare,       1)
            call co_broadcast(Lambda,      1)
            call co_broadcast(vCh,         1)
            call co_broadcast(Lambda_v,    1)
        end if

        ! ----------Info for packing parameters in/out of an array----------
        nParamsTotal = n_bare + 2*n_ch*n_bare + (n_ch*(n_ch+1))/2 + n_ch
        bare_end     = n_bare
        g_end        = bare_end + n_ch*n_bare
        Lam_end      = g_end + n_ch*n_bare
        v_end        = Lam_end + (n_ch*(n_ch+1))/2 ! Size of upper triangle
        Lam_v_end    = v_end + n_ch

        if (verbose .and. IamRoot) then
            write(*,*) 'Parameter ends: '
            write(*,*) bare_end, g_end, Lam_end, v_end, Lam_v_end
        end if

        ! Set pion mass to physical point
        m_pi = m_pi0

    end subroutine initialiseHEFT


    subroutine finaliseHEFT()
        implicit none
        deallocate(m_mes, m_bar, m_mes_phys, m_bar_phys &
            & , gBare, vCh, Lambda &
            & , Lambda_v, m_bare, m_bare_phys &
            & , slp_bar, slp_mes, id &
            & , ch_labels, bare_labels, ch, ch_AM, ch_pWave)
    end subroutine finaliseHEFT


    subroutine initialiseHEFTFinite()
        implicit none
        integer :: i
        integer :: file_finite = 128
        character(len=16) :: string

        ! -------------Read finite volume HEFT configs from file------------
        if (IamRoot) then
            open(file_finite, file='HEFTFinite.config', action='read')
            read(file_finite,*)
            read(file_finite,*)
            read(file_finite,*)

            read(file_finite,*) string, m_pi

            read(file_finite,*)
            read(file_finite,*)

            read(file_finite,*) string, L_min
            read(file_finite,*) string, L_max
            read(file_finite,*) string, L_points

            read(file_finite,*)
            read(file_finite,*)

            read(file_finite,*) string, m_pi_min
            read(file_finite,*) string, m_pi_max
            read(file_finite,*) string, m_pi_points
            read(file_finite,*) string, L_m_pi
            read(file_finite,*) string, doMassSlopeFitting
            read(file_finite,*) string, bareSlopeOrder
        end if

        if (doCA) call co_broadcast(bareSlopeOrder, 1)
        allocate(slp_bare_default(bareSlopeOrder))

        if (IamRoot) then
            read(file_finite,*) string, slp_bare_default(:)

            read(file_finite,*)
            read(file_finite,*)

            read(file_finite,*) fileName_lQCD

            read(file_finite,*)
            read(file_finite,*)

            read(file_finite,*) string, L_multifit
            close(file_finite)
        end if


        if (doCA) then
            call co_broadcast(m_pi,               1)
            call co_broadcast(L_min,              1)
            call co_broadcast(L_max,              1)
            call co_broadcast(L_points,           1)
            call co_broadcast(m_pi_min,           1)
            call co_broadcast(m_pi_max,           1)
            call co_broadcast(m_pi_points,        1)
            call co_broadcast(L_m_pi,             1)
            call co_broadcast(doMassSlopeFitting, 1)
            call co_broadcast(slp_bare_default,   1)
            call co_broadcast(L_multifit,         1)
        end if

        if (IamRoot) then
            allocate( massSlopeFileNames(n_bare) )
            do i = 1, n_bare
               massSlopeFileNames(i) = 'mass_' // trim(bare_labels(i)) &
                   & // '.in'
            end do
        end if

        allocate(slp_bare(n_bare,bareSlopeOrder))

    end subroutine initialiseHEFTFinite


    subroutine initialiseHEFTInfinite()
        implicit none

        integer :: i
        integer :: file_infinite = 128
        character(len=16) :: string

        if (IamRoot) then
            open(file_infinite, file='HEFTInfinite.config', action='read')

            read(file_infinite,*)
            read(file_infinite,*)
            read(file_infinite,*)

            read(file_infinite,*) string, E_init
            read(file_infinite,*) string, E_fin
            read(file_infinite,*) string, nPoints_inf
            read(file_infinite,*) string, useDataPoints

            close(file_infinite)
        end if

        if (doCA) then
            call co_broadcast(E_init, 1)
            call co_broadcast(E_fin, 1)
            call co_broadcast(nPoints_inf, 1)
            call co_broadcast(useDataPoints, 1)
        end if

        nPoles = 4
        allocate(E_poles(nPoles), E_poles_init(nPoles))

    end subroutine initialiseHEFTInfinite


    subroutine initialiseHEFTPoles()
      implicit none

    end subroutine initialiseHEFTPoles


    subroutine initialiseHEFTFitting()
        implicit none
        integer :: file_fitConfig
        integer :: kk
        real(DP) :: tempError
        real(DP), dimension(2) :: tempBounds
        character(len=4) :: tempStr
        integer :: ib

        if (IamRoot) then
            if (verbose) write(*,*) '------------Initialising HEFTFitting-------------'
            allocate( fitParams(nParamsTotal), isActiveParam(nParamsTotal) &
                & , fitParamLow(nParamsTotal), fitParamHigh(nParamsTotal) &
                & , fitParamsOrig(nParamsTotal), fitParamErrors(nParamsTotal) )

            if (verbose) write(*,'(a,\)') ' Opening fitting config'
            file_fitConfig = 301
            open(file_fitConfig, file='HEFTFitting.config', action='read')
            if (verbose) write(*,'(a,i0)') ' as ID: ', file_fitConfig
            read(file_fitConfig,*)

            doEtaFitting = .true.

            ! ----------------Read bounds and errors for each var---------------
            if (verbose) write(*,*) '     Bounds                 Error'
            if (n_bare.gt.0) then
                ! Bare mass(es))
                do ib = 1, n_bare
                   read(file_fitConfig,*) tempStr, tempBounds(:), tempError
                   if (verbose) write(*,*) tempStr, tempBounds, tempError
                   fitParamLow(ib) = tempBounds(1)
                   fitParamHigh(ib) = tempBounds(2)
                   fitParamErrors(ib) = tempError
                end do
                ! g
                read(file_fitConfig,*) tempStr, tempBounds(:), tempError
                if (verbose) write(*,*) tempStr, tempBounds, tempError
                fitParamLow(bare_end+1:g_end) = tempBounds(1)
                fitParamHigh(bare_end+1:g_end) = tempBounds(2)
                fitParamErrors(bare_end+1:g_end) = tempError
                ! Lam
                read(file_fitConfig,*) tempStr, tempBounds(:), tempError
                if (verbose) write(*,*) tempStr, tempBounds, tempError
                fitParamLow(g_end+1:Lam_end) = tempBounds(1)
                fitParamHigh(g_end+1:Lam_end) = tempBounds(2)
                fitParamErrors(g_end+1:Lam_end) = tempError
            end if

            ! v
            read(file_fitConfig,*) tempStr, tempBounds(:), tempError
            if (verbose) write(*,*) tempStr, tempBounds, tempError
            fitParamLow(Lam_end+1:v_end) = tempBounds(1)
            fitParamHigh(Lam_end+1:v_end) = tempBounds(2)
            fitParamErrors(Lam_end+1:v_end) = tempError
            ! Lam
            read(file_fitConfig,*) tempStr, tempBounds(:), tempError
            if (verbose) write(*,*) tempStr, tempBounds, tempError
            fitParamLow(v_end+1:Lam_v_end) = tempBounds(1)
            fitParamHigh(v_end+1:Lam_v_end) = tempBounds(2)
            fitParamErrors(v_end+1:Lam_v_end) = tempError

            if (verbose) then
                write(*,*) 'Fit param low'
                write(*,*) fitParamLow
                write(*,*) 'Fit param high'
                write(*,*) fitParamHigh
                write(*,*) 'Fit param error'
                write(*,*) fitParamErrors
            end if

            ! ------------------Read in which params are active-----------------
            read(file_fitConfig,*)
            read(file_fitConfig,*)
            read(file_fitConfig,*)
            if (verbose) write(*,*)
            if (verbose) write(*,*) 'Reading which parameters are active'
            read(file_fitConfig,*) isActiveParam(:)

            ! -------------------Set the active fit variables-------------------
            nParams = count(isActiveParam)
            if (verbose) write(*,'(x,i0,a)') nParams, ' parameters active'
            allocate( activeFitParams(nParams), fitParamIndex(nParams) )
            allocate( activeFitParamErrors(nParams), activeFitParamGuess(nParams) )
            allocate( activeFitParamLow(nParams), activeFitParamHigh(nParams) )
            fitParamIndex = pack( (/ (kk,kk=1,nParamsTotal) /), isActiveParam )

            activeFitParamErrors = pack( fitParamErrors, isActiveParam )
            activeFitParamLow    = pack( fitParamLow,    isActiveParam )
            activeFitParamHigh   = pack( fitParamHigh,   isActiveParam )
            ! activeFitParamGuess  = pack( fitParamGuess,  isActiveParam )

            if (verbose) then
                write(*,*) 'Active param low'
                write(*,*) activeFitParamLow
                write(*,*) 'Active param high'
                write(*,*) activeFitParamHigh
                write(*,*) 'Active param errors'
                write(*,*) activeFitParamErrors
            end if

            if (verbose) write(*,*) '-------------------------------------------------'
            if (verbose) write(*,*)
        end if

    end subroutine initialiseHEFTFitting


    ! ------------------------------------------------------------------
    ! -------------------------Useful functions-------------------------
    ! ------------------------------------------------------------------


    subroutine triangularise(M)
        ! Copies the upper triangle of a matrix to the lower triangle
        implicit none
        real(DP), dimension(:,:), intent(inout) :: M
        ! Local variables
        integer :: n
        integer :: ii, jj

        n = size(M,1)
        do jj= 1,n
           do ii = 1,n
              if (jj.gt.ii) then
                  M(jj,ii) = M(ii,jj)
              end if
           end do
        end do

    end subroutine triangularise


    function vectorToMat(vec, iM, jM) result(M)
        ! Converts a length N vector into an iM*jM matrix where N = iM*jM
        implicit none
        real(DP), dimension(:) :: vec
        integer :: iM, jM
        real(DP), dimension(iM, jM) :: M
        integer :: kk, iiM, jjM
        integer :: N

        N = size(vec)
        do kk = 1,N
           iiM = mod(kk-1,N) + 1
           jjM = (kk-1)/N + 1
           M(iiM,jjM) = vec(kk)
        end do

    end function vectorToMat


    function vectorToSymMat(vec) result(M)
        ! Converts a vector into a symmetric matrix
        implicit none
        real(DP), dimension(:), intent(in) :: vec
        real(DP), dimension(int(0.5*sqrt(8.0*size(vec)+1.0) - 0.5) &
            & , int(0.5*sqrt(8.0*size(vec)+1.0) - 0.5)) :: M
        integer :: size_v, size_M
        integer :: ii, jj, start, fin
        size_v = size(vec)
        size_M = int(0.5*sqrt(8.0*size_v+1.0) - 0.5)

        do jj = 1, size_M
           ii = size_M - jj + 1
           start = (ii**2 - ii) / 2 + 1
           fin   = (ii**2 + ii) / 2
           M(jj:, jj) = vec(size_v-fin+1 : size_v-start+1)
           M(jj, jj:) = M(jj:, jj)
        end do
    end function vectorToSymMat


    pure function int2str(int) result (str)
        ! Converts an integer into a string
        implicit none
        integer, intent(in) :: int
        character(len=10) :: str
        write(str,'(i0)') int
    end function int2str


    function timePrint(t_in) result(timeString)
        ! Gives a human readable time from an integer
        !    given by cpu_time
        implicit none
        real(DP), intent(in) :: t_in
        character(len=6) :: timeString
        real(DP) :: dt

        dt = t_in
        if (dt.le.120.0_DP) then
            write(timeString,'(f5.1,a1)') dt, 's'
        else
            dt = dt/60.0_DP

            if (dt.le.120.0_DP) then
                write(timeString,'(f5.1,a1)') dt, 'm'
            else
                dt = dt/60.0_DP

                if (dt.le.24.0_DP) then
                    write(timeString,'(f5.1,a1)') dt, 'h'
                else
                    dt = dt/24.0_DP
                    write(timeString,'(f5.1,a1)') dt, 'd'
                end if
            end if
        end if
    end function timePrint


    function symMatToVector(M) result(vec)
        ! Converts a symmetric matrix into an array, column-wise
        !  i.e. takes the upper triangle of the matrix and converts to a vector
        implicit none
        real(DP), dimension(:,:) :: M
        real(DP), dimension((size(M,1)*(size(M,1)+1))/2) :: vec
        ! Local variables
        integer :: dimM
        integer :: iM, jM

        dimM = size(M,1)
        do jM = 1,dimM
           do iM = 1,jM
              vec((jM*(jM-1))/2 + iM) = M(iM,jM)
           end do
        end do

    end function symMatToVector


    function vectorToSqrMat(vec) result(M)
        ! Converts a square vector to a square matrix
        !  e.g. a length 9 vector to a 3x3 matrix
        !  please dont use non square length vectors :(
        implicit none
        real(DP), dimension(:) :: vec
        real(DP), dimension(int(sqrt(real(size(vec)))) &
            & , int(sqrt(real(size(vec))))) :: M
        integer :: kk, iM, jM
        integer :: N

        N = size(vec)
        do kk = 1,N
           iM = mod(kk-1,N) + 1
           jM = (kk-1)/N + 1
           M(iM,jM) = vec(kk)
        end do

    end function vectorToSqrMat


    subroutine packFitParams()
        ! Packs the variables for fitting into a rank 1 array
        implicit none
        integer :: kk

        ! Bare state parameters
        if (n_bare.gt.0) then
            fitParams(1:bare_end) = m_bare
            fitParams(bare_end+1 : g_end) &
                & = reshape(gBare, (/n_ch*n_bare/))
            fitParams(g_end+1 : Lam_end) &
                & = reshape(Lambda, (/n_ch*n_bare/))
        end if

        fitParams(Lam_end+1 : v_end) &
            & = symMatToVector(vCh)
        fitParams(v_end+1 : Lam_v_end) = Lambda_v

        activeFitParams = pack( fitParams, isActiveParam )

    end subroutine packFitParams


    subroutine unpackFitParams(params)
        ! Unpacks the fit parameters back to the actual variables
        implicit none
        real(DP), dimension(:), intent(inout) :: params
        integer :: kk

        ! ! If one of the parameters depends on the value of another, set that here
        ! do kk = 1, size(params)
        !     params(kk) = params(correlationArray(kk)) * correlationValues(kk)
        ! end do

        if (n_bare.gt.0) then
            m_bare = params(1:bare_end)
            gBare  = vectorToMat(params(bare_end+1 : g_end), n_ch, n_bare)
            Lambda = vectorToMat(params(g_end+1 : Lam_end),  n_ch, n_bare)
        end if

        vCh = vectorToSymMat(params(Lam_end+1 : v_end))
        Lambda_v = params(v_end+1 : Lam_v_end)

    end subroutine unpackFitParams



    subroutine scrambleActiveParams(params, paramLow, paramHigh)
        implicit none
        ! Used to randomise the initial guess for parameter fitting
        real(DP), dimension(:), intent(out) :: params
        real(DP), dimension(:), intent(in) :: paramLow, paramHigh
        real(DP) :: randN
        integer :: iP

        do iP = 1, nParams
           call random_number(randN)
           params(iP) = paramLow(iP) &
               & + randN*(paramHigh(iP) - paramLow(iP))
        end do
        call unpackFitParams(params)

    end subroutine scrambleActiveParams


    function vectoriseUpperTriangle(M) result(vec)
        implicit none
        ! Converts an upper triangular matrix into an array, column-wise
        real(DP), dimension(:,:) :: M
        real(DP), dimension((size(M,1)*(size(M,1)+1))/2) :: vec
        ! Local variables
        integer :: dimM
        integer :: iM, jM

        dimM = size(M,1)
        do jM = 1,dimM
           do iM = 1,jM
              vec((jM*(jM-1))/2 + iM) = M(iM,jM)
           end do
        end do

    end function vectoriseUpperTriangle


    function pWave_to_AM(pW) result(AM)
        character(len=1), intent(in) :: pW
        integer :: AM

        select case(pW)
        case('s', 'S')
            AM = 0
        case('p', 'P')
            AM = 1
        case('d', 'D')
            AM = 2
        case('f', 'F')
            AM = 3
        case('g', 'G')
            AM = 4
        case('h', 'H')
            AM = 5
        case default
            AM = -1
        end select
    end function pWave_to_AM


    function AM_to_pWave(AM) result(pW)
        implicit none
        integer, intent(in) :: AM
        character(len=1) :: pW

        select case(AM)
        case(0)
            pW = 'S'
        case(1)
            pW = 'P'
        case(2)
            pW = 'D'
        case(3)
            pW = 'F'
        case(4)
            pW = 'G'
        case(5)
            pW = 'H'
        case default
            pW = 'Z'
        end select
    end function AM_to_pWave


    subroutine printCurrentParameters(iChoice)
        implicit none
        integer :: ii
        integer, optional :: iChoice
        if (IamRoot) then
            write(*,*) '--------------Current Parameter Set--------------'
            if (present(iChoice)) then
                write(*,'(a,i4)') ' Fit number: ', iChoice
            end if
            write(*,'(a4,20f12.5)') 'm',  m_bare_phys
            do ii = 1, n_bare
               write(*,'(a3,i1,20f12.5)') 'g',  ii, gBare(:,ii)
            end do
            do ii = 1, n_bare
               write(*,'(a3,i1,20f12.5)') 'L',  ii, Lambda(:,ii)
            end do
            do ii = 1, n_ch
               write(*,'(a3,i1,20f12.5)') 'v',  ii, vCh(ii,:)
            end do
            write(*,'(a4,20f12.5)') 'Lv', Lambda_v
            write(*,*) 'Potential Choices: ' // pot_choices
            write(*,*) '-------------------------------------------------'
            write(*,*)
        end if
    end subroutine printCurrentParameters



    ! ------------------------------------------------------------------
    ! -----------------------Potentials/Regulators----------------------
    ! ------------------------------------------------------------------


    ! -------------------------Actual potentials------------------------

    function g_i_A(k, ich, ib)
        ! Delta1232
        implicit none
        complex(DP), intent(in) :: k
        integer,  intent(in) :: ich
        integer,  intent(in) :: ib
        complex(DP) :: g_i_A
        ! Local variables
        complex(DP) :: omega_mes
        integer :: angMom

        angMom = pWaveConvert(ch_pWave(ich))
        omega_mes = sqrt(k**2 + m_mes(ich)**2)

        g_i_A = 1.0_DP/m_pi0 * k**angMom / sqrt(omega_mes) &
            & * u_k(k, Lambda(ich,ib), ch_pWave(ich))

    end function g_i_A

    function g_i_B(k, ich, ib)
        ! Standard S-wave (odd-parity nucleons, Lambda1405)
        implicit none
        complex(DP), intent(in) :: k
        integer,  intent(in) :: ich
        integer,  intent(in) :: ib
        complex(DP) :: g_i_B
        ! Local variables
        complex(DP) :: omega_mes, u_tilde
        integer :: angMom

        angMom = pWaveConvert(ch_pWave(ich))
        omega_mes = sqrt(k**2 + m_mes(ich)**2)

        g_i_B = sqrt(3.0_DP)/(2.0_DP*pi*f_pi) &
            & * k**angMom * sqrt(omega_mes) &
            & * u_k(k, Lambda(ich,ib), ch_pWave(ich))

    end function g_i_B

    function g_i_C(k, ich, ib)
        ! Roper with converted couplings
        implicit none
        complex(DP), intent(in) :: k
        integer,  intent(in) :: ich
        integer,  intent(in) :: ib
        complex(DP) :: g_i_C
        ! Local variables
        complex(DP) :: omega_mes
        integer :: angMom

        angMom = pWaveConvert(ch_pWave(ich)) ! convert P to 1, D to 2 etc.
        omega_mes = sqrt(k**2 + m_mes(ich)**2)

        g_i_C = &
            & sqrt(3.0_DP)/(2.0_DP*pi*f_pi) &
            & * k**angMom &
            & * u_k(k, Lambda(ich,ib), ch_pWave(ich)) &
            & / sqrt(omega_mes)

        ! if (ich.eq.1) then
        !     g_i_C = &
        !         & sqrt(3.0_DP)/(2.0_DP*pi*f_pi) &
        !         & * k &
        !         & * u_k(k, Lambda(ich,ib), ch_pWave(ich)) &
        !         & / sqrt(omega_mes)
        ! else if (ich.eq.2) then
        !     g_i_C = &
        !         & 1.0_DP/(sqrt(3.0_DP)*pi*f_pi) &
        !         & * k &
        !         & * u_k(k, Lambda(ich,ib), ch_pWave(ich)) &
        !         & / sqrt(omega_mes)
        ! else if (ich.eq.3) then
        !     g_i_C = &
        !         & 1.0_DP/(2.0_DP*pi) &
        !         & * u_k(k, Lambda(ich,ib), ch_pWave(ich)) &
        !         & / sqrt(omega_mes)
        ! end if

    end function g_i_C

    function g_i_D(k, ich, ib)
        ! Delta(1232) but with f_pi
        implicit none
        complex(DP), intent(in) :: k
        integer,  intent(in) :: ich
        integer,  intent(in) :: ib
        complex(DP) :: g_i_D
        ! Local variables
        complex(DP) :: omega_mes
        integer :: angMom

        angMom = pWaveConvert(ch_pWave(ich))
        omega_mes = sqrt(k**2 + m_mes(ich)**2)
        g_i_D = 1.0_DP/f_pi * k**angMom / sqrt(omega_mes) &
            & * u_k(k, Lambda(ich,ib), ch_pWave(ich))

    end function g_i_D

    function g_i_E(k, ich, ib)
        ! Roper with Zhan-Wei's couplings
        implicit none
        complex(DP), intent(in) :: k
        integer,  intent(in) :: ich
        integer,  intent(in) :: ib
        complex(DP) :: g_i_E
        ! Local variables
        complex(DP) :: omega_mes

        omega_mes = sqrt(k**2 + m_mes(ich)**2)

        if (ich.eq.1) then
            g_i_E = &
                & sqrt(3.0_DP)/(2.0_DP*pi*f_pi) &
                & * k &
                & * u_k(k, Lambda(ich,ib), ch_pWave(ich)) &
                & / sqrt(omega_mes)
        else if (ich.eq.2) then
            g_i_E = &
                & 1.0_DP/(sqrt(3.0_DP)*pi*f_pi) &
                & * k &
                & * u_k(k, Lambda(ich,ib), ch_pWave(ich)) &
                & / sqrt(omega_mes)
        else if (ich.eq.3) then
            g_i_E = &
                & 1.0_DP/(2.0_DP*pi) &
                & * u_k(k, Lambda(ich,ib), ch_pWave(ich)) &
                & / sqrt(omega_mes)
        end if

    end function g_i_E


    function f_i_A(k, ich)
        implicit none
        complex(DP), intent(in) :: k
        integer,  intent(in) :: ich
        ! Local variables
        complex(DP), dimension(1,1) :: f_i_A
        complex(DP) :: omega_mes
        integer :: angMom

        angMom = pWaveConvert(ch_pWave(ich)) ! convert P to 1, D to 2 etc.
        omega_mes = sqrt(k**2 + m_mes(ich)**2)

        f_i_A(1,1) = 1.0_DP/m_pi0 * k**angMom / omega_mes &
            & * u_k(k, Lambda_v(ich), ch_pWave(ich))

    end function f_i_A

    function f_i_B(k, ich)
        implicit none
        complex(DP), intent(in) :: k
        integer,  intent(in) :: ich
        ! Local variables
        complex(DP), dimension(1,1) :: f_i_B
        complex(DP) :: u_tilde
        integer :: angMom

        angMom = pWaveConvert(ch_pWave(ich))
        u_tilde = u_k(k, Lambda_v(ich), ch_pWave(ich)) &
            & * (sqrt(k**2 + m_pi0**2) + m_pi0) &
            & / sqrt(k**2 + m_pi0**2)
        f_i_B(1,1) = sqrt(3.0_DP)/(2.0_DP*pi*f_pi) &
            & * k**angMom * u_tilde

    end function f_i_B

    function f_i_C(k, ich)
        ! Similar to f_i_B, except with an additional low energy term for S-wave
        implicit none
        complex(DP), intent(in) :: k
        integer,  intent(in) :: ich
        ! Local variables
        complex(DP), dimension(1,1) :: f_i_C
        complex(DP) :: omega_mes
        complex(DP) :: u_tilde
        integer :: angMom

        angMom = pWaveConvert(ch_pWave(ich))
        u_tilde = u_k(k, Lambda_v(ich), ch_pWave(ich)) &
            & * (sqrt(k**2 + m_pi0**2) + m_pi0) &
            & / sqrt(k**2 + m_pi0**2)
        f_i_C(1,1) = sqrt(3.0_DP)/(2.0_DP*pi*f_pi) &
            & * k**angMom * u_tilde

        beta = 0.6_DP
        f_i_C(1, 1) = f_i_C(1, 1) / (k**2 + beta**2)**2

    end function f_i_C

    function f_i_D(k, ich)
        ! f_i_C, but with a different low energy enhancement term
        implicit none
        complex(DP), intent(in) :: k
        integer,  intent(in) :: ich
        ! Local variables
        complex(DP), dimension(1,1) :: f_i_D
        complex(DP) :: omega_mes
        integer :: angMom

        angMom = pWaveConvert(ch_pWave(ich))
        f_i_D(1,1) = sqrt(3.0_DP)/(2.0_DP*pi*f_pi) &
            & * u_k(k, Lambda_v(ich), ch_pWave(ich)) * k**angMom

        beta = 0.6_DP
        ! beta = 0.0_DP
        f_i_D(1, 1) = f_i_D(1, 1) / (k**2 + beta**2)**2

    end function f_i_D

    function f_i_E(k, ich)
        ! Roper with converted couplings
        implicit none
        complex(DP), intent(in) :: k
        integer,  intent(in) :: ich
        ! Local variables
        complex(DP), dimension(1,1) :: f_i_E
        complex(DP) :: omega_mes
        integer :: angMom

        angMom = pWaveConvert(ch_pWave(ich))
        omega_mes = sqrt(k**2 + m_mes(ich)**2)

        f_i_E(1,1) = &
            & sqrt(3.0_DP)/(2.0_DP*pi*f_pi) &
            & * k**angMom &
            & * u_k(k, Lambda_v(ich), ch_pWave(ich)) &
            & / omega_mes

    end function f_i_E

    function f_i_F(k, ich)
        implicit none
        complex(DP), intent(in) :: k
        integer,  intent(in) :: ich
        ! Local variables
        complex(DP), dimension(1,1) :: f_i_F
        complex(DP) :: omega_mes
        integer :: angMom

        angMom = pWaveConvert(ch_pWave(ich)) ! convert P to 1, D to 2 etc.
        omega_mes = sqrt(k**2 + m_mes(ich)**2)
        f_i_F(1,1) = 1.0_DP/f_pi * k**angMom / omega_mes &
            & * u_k(k, Lambda_v(ich), ch_pWave(ich))

    end function f_i_F


    function f_i_G(k, ich)
        ! Roper paper potential
        implicit none
        complex(DP), intent(in) :: k
        integer,  intent(in) :: ich
        ! Local variables
        complex(DP), dimension(1,1) :: f_i_G
        complex(DP) :: omega_mes
        integer :: angMom

        omega_mes = sqrt(k**2 + m_mes(ich)**2)

        ! f_i_G(1,1) = &
        !     & sqrt(3.0_DP)/(2.0_DP*pi*f_pi) &
        !     & / omega_mes &
        !     & * k**angMom &
        !     & * u_k(k, Lambda_v(ich), ch_pWave(ich))

        if (ich.eq.1) then
            f_i_G(1,1) = &
                & sqrt(3.0_DP)/(2.0_DP*pi*f_pi) &
                & / omega_mes &
                & * k &
                & * u_k(k, Lambda_v(ich), ch_pWave(ich))
        else if (ich.eq.2) then
            f_i_G(1,1) = &
                & 1.0_DP/(sqrt(3.0_DP)*pi*f_pi) &
                & / omega_mes &
                & * k &
                & * u_k(k, Lambda_v(ich), ch_pWave(ich))
        else if (ich.eq.3) then
            f_i_G(1,1) = &
                & 1.0_DP/(2.0_DP*pi) &
                & / omega_mes &
                & * u_k(k, Lambda_v(ich), ch_pWave(ich))
        end if

    end function f_i_G


    function f_i_W(k, ich)
        ! Weinberg-Tamozawa potential
        implicit none
        complex(DP), intent(in) :: k
        integer,  intent(in) :: ich
        ! Local variables
        complex(DP), dimension(2,1) :: f_i_W
        complex(DP) :: omega_mes

        omega_mes = sqrt(k**2 + m_mes(ich)**2)

        f_i_W(1,1) = omega_mes / sqrt(2.0_DP*omega_mes)
        f_i_W(2,1) = 1.0_DP / sqrt(2.0_DP*omega_mes)

        f_i_W(:,:) = f_i_W(:,:) / (sqrt(8.0_DP)*pi*f_pi) &
            & * u_k(k, Lambda_v(ich), ch_pWave(ich))

    end function f_i_W


    ! Dipole regulator
    function u_k_A(k, Lam, pWave)
        implicit none
        complex(DP) :: k
        real(DP) :: Lam
        complex(DP) :: u_k_A
        character(len=1) :: pWave
        integer :: angMom

        angMom = pWaveConvert(pWave)
        if (angMom.le.1) then
            u_k_A = 1.0_DP &
                & / (1.0_DP + (k/Lam)**2)**2
        else
            u_k_A = 1.0_DP &
                & / (1.0_DP + (k/Lam)**2)**3
        end if

    end function u_k_A

    function inverse_u_k_A(u, Lam)
        implicit none
        real(DP), intent(in) :: u, Lam
        ! Local variables
        real(DP) :: inverse_u_k_A

        inverse_u_k_A = Lam * sqrt(1.0_DP/sqrt(u) - 1.0_DP)

    end function inverse_u_k_A


    ! Guassian regulator
    function u_k_B(k, Lam, pWave)
        implicit none
        complex(DP) :: k
        real(DP) :: Lam
        character(len=1) :: pWave
        complex(DP) :: u_k_B

        if ((pWave .eq. 'd') .or. (pWave .eq. 'D')) then
            u_k_B = exp( -k**4 / Lam**4 )
        else
            u_k_B = exp( -abs(k/Lam)**2 )
        end if

    end function u_k_B

    function inverse_u_k_B(u, Lam)
        implicit none
        real(DP), intent(in) :: u, Lam
        ! Local variables
        real(DP) :: inverse_u_k_B

        inverse_u_k_B = Lam * sqrt(-log(u))

    end function inverse_u_k_B


    ! Step function reg
    function u_k_C(k, Lam, pWave)
        implicit none
        complex(DP) :: k
        real(DP) :: Lam
        character(len=1) :: pWave
        complex(DP) :: u_k_C

        if (abs(k) .le. Lam) then
            u_k_C = 1.0_DP
        else
            u_k_C = 0.0_DP
        end if

    end function u_k_C

    function inverse_u_k_C(u, Lam)
        implicit none
        real(DP), intent(in) :: u, Lam
        ! Local variables
        real(DP) :: inverse_u_k_C

        inverse_u_k_C = Lam

    end function inverse_u_k_C


    ! Sharpe function reg
    function u_k_D(k, Lam, pWave)
        implicit none
        complex(DP) :: k
        real(DP) :: Lam
        character(len=1) :: pWave
        complex(DP) :: u_k_D
        complex(DP) :: zReg
        integer :: nPow

        zReg = abs(k / Lam)
        nPow = 1

        if (abs(zReg) .le. 0.0_DP) then
            u_k_D = 1.0_DP
        else if (abs(zReg) .ge. 1.0_DP) then
            u_k_D = 0.0_DP
        else
            u_k_D = 1.0_DP - exp( -exp(-(1.0_DP/(1.0_DP-zReg)**nPow)) &
                & / zReg**nPow )
        end if

    end function u_k_D

    function inverse_u_k_D(u, Lam)
        implicit none
        real(DP), intent(in) :: u, Lam
        ! Local variables
        real(DP) :: inverse_u_k_D

        inverse_u_k_D = Lam

    end function inverse_u_k_D

    ! ------------------------Choosing functions------------------------
    ! Dipole regulator
    function u_k_cmplx(k, Lam, pWave)
        implicit none
        complex(DP) :: k
        real(DP) :: Lam
        character(len=1), optional :: pWave
        complex(DP) :: u_k_cmplx

        if (.not. present(pWave)) then
            pWave = 'Z'
        end if

        select case (u_choice)
        case('A')
            u_k_cmplx = u_k_A(k, Lam, pWave)
        case('B')
            u_k_cmplx = u_k_B(k, Lam, pWave)
        case('C')
            u_k_cmplx = u_k_C(k, Lam, pWave)
        case('D')
            u_k_cmplx = u_k_D(k, Lam, pWave)
        end select
    end function u_k_cmplx

    function u_k_real(k, Lam, pWave)
        implicit none
        real(DP), intent(in) :: k, Lam
        character(len=1), optional :: pWave
        ! Local variables
        real(DP) :: u_k_real

        if (.not. present(pWave)) then
            pWave = 'Z'
        end if

        u_k_real = real(u_k_cmplx(cmplx(k,0.0_DP,DP) &
            & , Lam, pWave),DP)
    end function u_k_real

    function inverse_u_k(u, Lam)
        implicit none
        real(DP), intent(in) :: u, Lam
        ! Local variables
        real(DP) :: inverse_u_k
        select case (u_choice)
        case('A')
            inverse_u_k = inverse_u_k_A(u, Lam)
        case('B')
            inverse_u_k = inverse_u_k_B(u, Lam)
        end select
    end function inverse_u_k


    function f_i_cmplx(k, ich)
        implicit none
        complex(DP), intent(in) :: k
        integer,  intent(in) :: ich
        ! Local variables
        complex(DP), dimension(n_f,1) :: f_i_cmplx
        complex(DP) :: omega_mes

        select case (f_choice)
        case('A')
            f_i_cmplx = f_i_A(k, ich)
        case('B')
            f_i_cmplx = f_i_B(k, ich)
        case('C')
            f_i_cmplx = f_i_C(k, ich)
        case('D')
            f_i_cmplx = f_i_D(k, ich)
        case('E')
            f_i_cmplx = f_i_E(k, ich)
        case('F')
            f_i_cmplx = f_i_F(k, ich)
        case('G')
            f_i_cmplx = f_i_G(k, ich)
        case('W')
            f_i_cmplx = f_i_W(k, ich)
        end select

    end function f_i_cmplx


    function f_i_real(k, ich)
        implicit none
        real(DP), intent(in) :: k
        integer,  intent(in) :: ich
        ! Local variables
        real(DP), dimension(n_f,1) :: f_i_real
        real(DP) :: omega_mes

        f_i_real = real(f_i_cmplx(cmplx(k,0.0_DP,DP), ich), DP)

    end function f_i_real


    function V_ij_cmplx(ki, i, kj, j)
        implicit none
        complex(DP), intent(in) :: ki, kj
        integer,  intent(in) :: i, j
        complex(DP) :: V_ij_cmplx
        complex(DP), dimension(1,1) :: V_ij_cmplx_matrix

        if (n_f.eq.1) then
            V_ij_cmplx_matrix = f_i(ki, i) * vCh(i,j) * f_i(kj, j)

        else if (n_f.eq.2) then
            ! Weinberg-Tamozawa potential (potential "W")
            ! (f^T*v)*f = f^T*(v*f)
            V_ij_cmplx_matrix = matmul( transpose(f_i(ki,i)), &
                & matmul(reshape([ 0.0_DP, vCh(i,j) &
                & , vCh(i,j), 0.0_DP ], [2,2]), f_i(kj,j)) )
        end if

        V_ij_cmplx = V_ij_cmplx_matrix(1,1)

    end function V_ij_cmplx


    function V_ij_real(ki, i, kj, j)
        implicit none
        real(DP), intent(in) :: ki, kj
        integer,  intent(in) :: i, j
        real(DP) :: V_ij_real

        V_ij_real = real(V_ij_cmplx(cmplx(ki,0.0_DP,DP), i &
            & , cmplx(kj,0.0_DP,DP), j), DP)

    end function V_ij_real


    function g_i_cmplx(k, ich, ib)
        implicit none
        complex(DP), intent(in) :: k
        integer,  intent(in) :: ich
        integer,  intent(in) :: ib
        complex(DP) :: g_i_cmplx

        select case (g_choice)
        case('A')
            g_i_cmplx = g_i_A(k, ich, ib)
        case('B')
            g_i_cmplx = g_i_B(k, ich, ib)
        case('C')
            g_i_cmplx = g_i_C(k, ich, ib)
        case('D')
            g_i_cmplx = g_i_D(k, ich, ib)
        case('E')
            g_i_cmplx = g_i_E(k, ich, ib)
        end select

    end function g_i_cmplx


    function g_i_real(k, ich, ib)
        implicit none
        real(DP), intent(in) :: k
        integer,  intent(in) :: ich
        integer,  intent(in) :: ib
        real(DP) :: g_i_real

        g_i_real = real(g_i_cmplx(cmplx(k,0.0_DP,DP), ich, ib), DP)

    end function g_i_real


end module heft
