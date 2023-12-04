module potentials
  use kinds
  use numFort
  implicit none

  real(DP), parameter :: f_pi = 0.0924_DP ! GeV
  real(DP), parameter :: m_pi0 = 0.1385_DP ! GeV
  real(DP), parameter :: u_kmax = 0.01_DP ! Limits k_max
  real(DP), parameter :: hbar_c = 0.1973_DP ! GeVfm
  real(DP), parameter :: GeVfm  = 0.1973_DP ! hbar * c
  real(DP), parameter :: GeV2_mb = 0.3894_DP ! hbar * c (mb = millibarns)
  real(DP), parameter :: degrees = 180.0_DP / pi ! / rad
  real(DP), parameter :: M_octet = 0.920_DP

  ! Coarray variables
  integer :: iImg, nImg
  logical :: IamRoot
  integer, parameter :: root_image = 1

  ! ! Couplings SU6
  ! real(DP), parameter :: D_SU6 = 0.76_DP
  ! real(DP) :: C_SU6 = -2.0_DP * D_SU6
  ! real(DP) :: chi_Delta = 3.0_DP/32.0_DP/pi/f_pi**2 * 2.0_DP/9.0_DP * C_SU6**2

  integer  :: n_bare
  integer  :: n_ch
  integer  :: n_init_k
  real(DP) :: E_gl ! Global energy
  real(DP) :: E_init, E_fin
  real(DP) :: m_pi
  logical  :: readData ! If energy array should be read or generated from E_init, E_fin
  real(DP) :: L_min, L_max
  real(DP) :: L_renorm
  real(DP) :: L_m_pi
  real(DP) :: fit_chi2
  integer  :: L_points

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
  real(DP), dimension(:), allocatable :: fitParams, activeFitParams
  real(DP), dimension(:), allocatable :: fitParamErrors, activeFitParamErrors
  logical,  dimension(:), allocatable :: isActiveParam
  logical,  dimension(:), allocatable :: isVaryingParam
  integer,  dimension(:), allocatable :: fitParamIndex, varyingParamIndex
  real(DP), dimension(:), allocatable :: fitParamLow, fitParamHigh
  real(DP), dimension(:), allocatable :: activeFitParamLow, activeFitParamHigh
  real(DP), dimension(:), allocatable :: fitParamGuess, activeFitParamGuess
  real(DP), dimension(:), allocatable :: varyingParamMin, varyingParamMax
  integer,  dimension(:), allocatable :: correlationArray
  real(DP), dimension(:), allocatable :: correlationValues
  logical, dimension(:), allocatable  :: printedParams
  character(len=64) :: printedParamsLabel

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
  real(DP), dimension(:), allocatable :: slp_bare ! bare mass
  real(DP), dimension(:), allocatable :: slp_mes  ! meson mass
  real(DP), dimension(:), allocatable :: slp_bar  ! baryon mass
  ! Mass slope fitting variables
  logical :: doMassSlopeFitting
  character(len=32), dimension(:), allocatable :: massSlopeFileNames ! dim=nBare

  logical :: isEnergyDepPotential
  character(len=9), dimension(:), allocatable :: ch_labels
  character(len=9), dimension(:), allocatable :: bare_labels
  integer :: ch_onshell

  character(len=32) :: fitParamsFileName
  character(len=32) :: dummyStr

  ! Pole search variables
  complex(DP) :: E_pole
  complex(DP), dimension(:), allocatable :: k_pole
  complex(DP), dimension(:), allocatable :: E_poles, E_poles_init
  integer :: nPoles
  integer, dimension(:), allocatable :: poles_bare_loc


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

contains




  subroutine initialiseParams()
    implicit none
    ! Physical masses (GeV)
    real(DP) :: m_N0 = 0.9385_DP
    real(DP) :: m_Delta0 = 1.232_DP
    ! real(DP) :: m_Delta0 = 1.193_DP
    real(DP) :: slp_N = 1.435_DP
    real(DP) :: slp_Delta = 1.079758_DP

    ! ! ----------------------------Good Stuff----------------------------
    ! Fixed Lambda 0.8 GeV, chi2 = 290
    real(DP), parameter :: fm_bare           =  1.3766735565540680_DP
    real(DP), parameter :: fg_piN            =  0.1363879430157954_DP
    real(DP), parameter :: fg_piDelta        =  0.1111802369967659_DP
    real(DP), parameter :: fLambda_piN       =  0.8000000000000000_DP
    real(DP), parameter :: fLambda_piDelta   =  0.8000000000000000_DP
    real(DP), parameter :: fv_piNpiN         = -0.0131918928769021_DP
    real(DP), parameter :: fv_piNpiDelta     = -0.0735282162476920_DP
    real(DP), parameter :: fv_piDeltapiDelta = -0.0201388085295084_DP
    real(DP), parameter :: fLambda_v_piN     =  0.8000000000000000_DP
    real(DP), parameter :: fLambda_v_piDelta =  0.8000000000000000_DP
    ! ! ------------------------------------------------------------------

    integer :: kk
    integer :: iFileParams, iParams, iParamChoice
    logical :: readParamsFromFile
    character(len=16) :: fmt_readFileParams
    real(DP), dimension(:), allocatable :: fileParams

    ! Delta deets
    n_ch = 2
    n_bare = 1
    n_init_k = 1
    E_init = 1.10364_DP
    E_fin = 1.64331_DP
    isEnergyDepPotential = .false.
    readData = .true.

    L_min = 3.0_DP
    L_max = 8.0_DP
    L_points = 100
    L_renorm = 2.98556_DP
    ! L_renorm = 4.0_DP
    ! L_m_pi = 6.54_DP
    ! L_m_pi = 3.27_DP
    L_m_pi = 2.9856_DP

    allocate(m_mes(n_ch), m_bar(n_ch), m_mes_phys(n_ch), m_bar_phys(n_ch) &
         & , gBare(n_ch,n_bare), vCh(n_ch,n_ch), Lambda(n_ch,n_bare) &
         & , Lambda_v(n_ch), m_bare(n_bare), m_bare_phys(n_bare) &
         & , slp_bare(n_bare), slp_bar(n_ch), slp_mes(n_ch), id(n_ch,n_ch) &
         & , ch_labels(n_ch), bare_labels(n_bare))

    ch_labels = (/ 'piN', 'piDelta' /)
    ch_onshell = 1
    bare_labels = (/ 'Delta' /)

    m_bar_phys = (/ m_N0, m_Delta0 /)
    m_mes_phys = (/ m_pi0, m_pi0 /)
    m_bare_phys = fm_bare

    slp_bar = (/ slp_N, slp_Delta /)
    slp_mes = 1.0_DP
    slp_bare = slp_Delta

    m_pi = m_pi0
    m_bar = m_bar_phys
    m_mes = m_mes_phys
    m_bare = m_bare_phys

    ! For g_piDelta, see Young et al. (2002) (Chiral Analysis...)
    gBare(:,1) = (/ fg_piN, fg_piDelta /)
    ! gBare(:,1) = (/ fg_piN, 25.0_DP/8.0_DP*fg_piN /)
    vCh(1,:)   = (/ fv_piNpiN, fv_piNpiDelta /)
    vCh(2,:)   = (/ fv_piNpiDelta, fv_piDeltapiDelta /)

    Lambda(:,1)  = (/ fLambda_piN, fLambda_piDelta /)
    Lambda_v(:) = (/ fLambda_v_piN, fLambda_v_piDelta /)

    id(:,:) = real(identityMat(n_ch),DP)
    E_gl = m_mes(ch_onshell) + m_bar(ch_onshell)

    ! phaseIndex = 1 ! which channel is on shell

    ! Info about parameter locations when packing into an array
    nParamsTotal = n_bare + 2*n_ch*n_bare + (n_ch*(n_ch+1))/2 + n_ch
    bare_end  = n_bare
    g_end     = bare_end + n_ch*n_bare
    Lam_end   = g_end + n_ch*n_bare
    v_end     = Lam_end + (n_ch*(n_ch+1))/2 ! Size of upper triangle
    Lam_v_end = v_end + n_ch

    ! -----------Read in the parameters from the allFits file-----------
    readParamsFromFile = .true.
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
                m_bare = fileParams(1:bare_end_file)
                m_bare_phys = m_bare
                gBare  = vectorToMat(fileParams(bare_end_file+1 : g_end_file) &
                     & , n_bare, n_ch)
                Lambda = vectorToMat(fileParams(g_end_file+1 : Lam_end_file) &
                     & ,  n_bare, n_ch)
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
    end if

    ! ----------------------Bare mass slope fitting---------------------
    doMassSlopeFitting = .true.
    allocate(massSlopeFileNames(n_bare))
    massSlopeFileNames(1) = 'DeltaMass.in'
    ! ------------------------------------------------------------------

    ! ----------------------Fitting Parameter Stuff---------------------
    ! Fitting initialisations
    doFitting = .true.
    doEtaFitting = .true.
    doCSFitting = .false.

    ! bare_mass + (g + Lambda) + vCh + Lambda_v
    nParamsTotal = n_bare + 2*n_ch*n_bare + (n_ch*(n_ch+1))/2 + n_ch
    allocate( fitParams(nParamsTotal), isActiveParam(nParamsTotal) &
         & , fitParamLow(nParamsTotal), fitParamHigh(nParamsTotal) &
         & , isVaryingParam(nParamsTotal), varyingParamMin(nParamsTotal) &
         & , varyingParamMax(nParamsTotal), correlationArray(nParamsTotal) &
         & , correlationValues(nParamsTotal), printedParams(nParamsTotal) )

    ! Only change this
    !                  m, g, g, L, L, v11, v12, v22, Lv, Lv
    ! isActiveParam = (/ 1, 1, 1, 1, 1, 1,   1,   1,   1,  1 /)
    ! isActiveParam = (/ 0, 1, 1, 1, 1, 1,   1,   1,   1,  1 /)
    ! isActiveParam = (/ 1, 1, 1, 0, 0, 1,   1,   1,   1,  1 /)
    ! isActiveParam = (/ 1, 1, 0, 1, 1, 1,   1,   1,   1,  1 /)
    isActiveParam = (/ 1, 1, 1, 0, 0, 1,   1,   1,   0,  0 /)
    ! isActiveParam = (/ 1, 1, 0, 0, 0, 1,   1,   1,   0,  0 /)
    ! isActiveParam = (/ 1, 1, 1, 0, 0, 0,   0,   0,   0,  0 /)
    ! isActiveParam = (/ 0, 1, 0, 0, 0, 0,   0,   0,   0,  0 /)

    !                   m, g, g, L, L, v11, v12, v22, Lv, Lv
    isVaryingParam = (/ 0, 0, 0, 1, 1, 0,   0,   0,   1,  1 /)
    ! isVaryingParam = (/ 0, 0, 0, 1, 0, 0,   0,   0,   0,  0 /)
    varyingParamMin = (/ 0.0_DP, 0.0_DP, 0.0_DP, fLambda_piN, fLambda_piDelta &
         & , 0.0_DP, 0.0_DP, 0.0_DP, fLambda_v_piN, fLambda_v_piDelta /)
    ! varyingParamMax = (/ 0.0_DP, 0.0_DP, 8.0_DP, 8.0_DP, 0.0_DP &
    !      & , 0.0_DP, 0.0_DP, 8.0_DP, 8.0_DP /)
    ! varyingParamMax = (/ 0.0_DP, 0.0_DP, 4.0_DP, 4.0_DP, 0.0_DP &
    !      & , 0.0_DP, 0.0_DP, 4.0_DP, 4.0_DP /)
    ! varyingParamMax = (/ 0.0_DP, 0.0_DP, 0.0_DP, 1.2_DP, 1.2_DP, 0.0_DP &
    !      & , 0.0_DP, 0.0_DP, 1.2_DP, 1.2_DP /)
    varyingParamMax = (/ 0.0_DP, 0.0_DP, 0.0_DP, 1.0_DP, 1.0_DP, 0.0_DP &
         & , 0.0_DP, 0.0_DP, 1.0_DP, 1.0_DP /)
    ! varyingParamMax = (/ 0.0_DP, 0.0_DP, 0.0_DP, 0.85_DP, 0.85_DP, 0.0_DP &
    !      & , 0.0_DP, 0.0_DP, 0.85_DP, 0.85_DP /)

    ! Gives the index of which parameter each one depends one
    !    e.g. for g_piDelta = C*g_piN, we have a 2 in i=3
    correlationArray(:)  = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 /)
    correlationValues(:) = 1.0_DP
    ! correlationArray(:)  = (/ 1, 2, 2, 4, 5, 6, 7, 8, 9, 10 /)
    ! correlationValues(:) = (/ 1.0_DP, 1.0_DP, 25.0_DP/8.0_DP, 1.0_DP, 1.0_DP &
    !      & , 1.0_DP, 1.0_DP, 1.0_DP, 1.0_DP, 1.0_DP /)

    !                  m, g, g, L, L, v11, v12, v22, Lv, Lv
    printedParams = (/ 1, 1, 1, 0, 0, 1,   0,   1,   0,  0 /)
    write(printedParamsLabel,'(5a9)') 'm', 'g1', 'g2', 'v11', 'Lambda'

    nParams = count(isActiveParam)
    allocate( activeFitParams(nParams), fitParamIndex(nParams) &
         & , fitParamErrors(nParams), activeFitParamErrors(nParams) &
         & , activeFitParamLow(nParams), activeFitParamHigh(nParams) &
         & , activeFitParamGuess(nParams), varyingParamIndex(nParams))
    fitParamIndex = pack( (/ (kk,kk=1,nParamsTotal) /), isActiveParam )
    varyingParamIndex = pack( (/ (kk,kk=1,nParamsTotal) /), isVaryingParam )


    !                  m, g, g, L, L, v11, v12, v22, Lv, Lv
    ! fitParamErrors = (/ 1e-3_DP, 1e-3_DP, 1e-3_DP &
    !      & , 1e-3_DP, 1e-3_DP, 1e-3_DP, 1e-3_DP, 1e-3_DP &
    !      & , 1e-3_DP, 1e-3_DP /)
    ! fitParamErrors = (/ 1e-4_DP, 1e-4_DP, 1e-4_DP &
    !      & , 1e-4_DP, 1e-4_DP, 1e-4_DP, 1e-4_DP, 1e-4_DP &
    !      & , 1e-4_DP, 1e-4_DP /)
    fitParamErrors = (/ 1e-5_DP, 1e-4_DP, 1e-4_DP &
         & , 1e-5_DP, 1e-5_DP, 1e-3_DP, 1e-3_DP, 1e-3_DP &
         & , 1e-5_DP, 1e-5_DP /)
    ! fitParamErrors = (/ 1e-6_DP, 1e-6_DP, 1e-6_DP &
    !      & , 1e-6_DP, 1e-6_DP, 1e-6_DP, 1e-6_DP, 1e-6_DP &
    !      & , 1e-6_DP, 1e-6_DP /)
    ! fitParamErrors = (/ 1.0_DP, 1e-3_DP, 1e-3_DP &
    !      & , 1.0_DP, 1.0_DP, 1e-3_DP, 1e-3_DP, 1e-3_DP &
    !      & , 1.0_DP, 1.0_DP /)
    fitParamLow = (/ 1.2_DP, 0.0001_DP, 0.0001_DP, 0.6_DP, 0.6_DP &
         & , -0.9_DP, -0.9_DP, -0.9_DP, 0.6_DP, 0.6_DP /)
    fitParamHigh = (/ 5.0_DP, 0.9_DP, 0.9_DP, 3.0_DP, 3.0_DP &
         & , -0.00001_DP, -0.00001_DP, -0.00001_DP, 3.0_DP, 3.0_DP /)
    fitParamGuess = (/ 1.5_DP, 0.1_DP, 0.1_DP, 0.8_DP, 0.8_DP &
         & , -0.1_DP, -0.1_DP, -0.1_DP, 0.8_DP, 0.8_DP /)


    activeFitParamErrors = pack( fitParamErrors, isActiveParam )
    activeFitParamLow    = pack( fitParamLow,    isActiveParam )
    activeFitParamHigh   = pack( fitParamHigh,   isActiveParam )
    activeFitParamGuess  = pack( fitParamGuess,  isActiveParam )
    printFittingOutput   = .false.

    fitParamsFileName = 'fitParams.params'


  end subroutine initialiseParams



  subroutine finaliseParams()
    implicit none
    deallocate(m_mes, m_bar, m_bar_phys, m_mes_phys, gBare &
         & , vCh, Lambda, Lambda_v, m_bare, m_bare_phys &
         & , slp_bare, slp_bar, slp_mes, id, ch_labels &
         & , bare_labels, fitParams, isActiveParam, activeFitParams &
         & , fitParamIndex, fitParamErrors, activeFitParamErrors &
         & , fitParamLow, fitParamHigh, activeFitParamLow &
         & , activeFitParamHigh, activeFitParamGuess, isVaryingParam &
         & , varyingParamMin, varyingParamMax, varyingParamIndex &
         & , correlationValues, correlationArray, printedParams )

  end subroutine finaliseParams



  function C_ij_E()
    implicit none
    real(DP), dimension(2,2) :: C_ij_E

    ! ------------------------------------------------------------------
    real(DP), parameter :: fv_piNpiN         = -0.0103014637890360_DP
    real(DP), parameter :: fv_piNpiDelta     = -0.0810565209400564_DP
    real(DP), parameter :: fv_piDeltapiDelta = -0.0014744087735186_DP
    ! ------------------------------------------------------------------

    C_ij_E(1,:)   = (/ fv_piNpiN, fv_piNpiDelta /)
    C_ij_E(2,:)   = (/ fv_piNpiDelta, fv_piDeltapiDelta /)

  end function C_ij_E


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
    real(DP), dimension(:) :: vec
    ! This is disgusting, please redo this
    real(DP), dimension(int(0.5*sqrt(8.0*size(vec)+1.0) - 0.5) &
         & , int(0.5*sqrt(8.0*size(vec)+1.0) - 0.5)) :: M
    integer :: size_M
    integer :: kk
    integer :: vecStart, vecEnd

    size_M = int(0.5_DP*sqrt(8.0_DP*size(vec)+1.0_DP) - 0.5_DP)
    do kk = 1,size_M
       vecStart = kk*(kk-1)/2 + 1
       vecEnd = kk*(kk+1)/2
       M(1:kk,kk) = vec(vecStart:vecEnd)
       M(kk,1:kk) = vec(vecStart:vecEnd)
    end do

  end function vectorToSymMat


  pure function int2str(int) result (str)
    implicit none
    integer, intent(in) :: int
    character(len=10) :: str
    write(str,'(i0)') int
  end function int2str




  ! Dipole regulator
  function u_k_cmplx(k, Lam)
    implicit none
    complex(DP) :: k
    real(DP) :: Lam
    complex(DP) :: u_k_cmplx

    u_k_cmplx = 1.0_DP &
         & / (1.0_DP + (k/Lam)**2)**2

  end function u_k_cmplx

  function u_k_real(k, Lam)
    implicit none
    real(DP), intent(in) :: k, Lam
    ! Local variables
    real(DP) :: u_k_real

    u_k_real = real(u_k_cmplx(cmplx(k,0.0_DP,DP), Lam),DP)

  end function u_k_real

  function inverse_u_k(u, Lam)
    implicit none
    real(DP), intent(in) :: u, Lam
    ! Local variables
    real(DP) :: inverse_u_k

    inverse_u_k = Lam * sqrt(1.0_DP/sqrt(u) - 1.0_DP)

  end function inverse_u_k

  function f_i_cmplx(k, ich)
    implicit none
    complex(DP), intent(in) :: k
    integer,  intent(in) :: ich
    ! Local variables
    complex(DP), dimension(1,1) :: f_i_cmplx
    complex(DP) :: omega_mes

    omega_mes = sqrt(k**2 + m_mes(ich)**2)
    f_i_cmplx(1,1) = 1.0_DP/m_pi * k / omega_mes * u_k(k,Lambda_v(ich))

  end function f_i_cmplx

  function f_i_real(k, ich)
    implicit none
    real(DP), intent(in) :: k
    integer,  intent(in) :: ich
    ! Local variables
    real(DP), dimension(1,1) :: f_i_real
    real(DP) :: omega_mes

    f_i_real = real(f_i_cmplx(cmplx(k,0.0_DP,DP), ich), DP)

  end function f_i_real


  function V_ij_cmplx(ki, i, kj, j)
    implicit none
    complex(DP), intent(in) :: ki, kj
    integer,  intent(in) :: i, j
    complex(DP) :: V_ij_cmplx
    complex(DP), dimension(1,1) :: V_ij_cmplx_matrix

    V_ij_cmplx_matrix = f_i(ki, i) * vCh(i,j) * f_i(kj, j)
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
    ! Local variables
    complex(DP) :: omega_mes

    omega_mes = sqrt(k**2 + m_mes(ich)**2)
    g_i_cmplx = 1.0_DP/m_pi * k/sqrt(omega_mes) * u_k(k,Lambda(ich,ib))

  end function g_i_cmplx

  function g_i_real(k, ich, ib)
    implicit none
    real(DP), intent(in) :: k
    integer,  intent(in) :: ich
    integer,  intent(in) :: ib
    real(DP) :: g_i_real

    g_i_real = real(g_i_cmplx(cmplx(k,0.0_DP,DP), ich, ib), DP)

  end function g_i_real

  ! function g_i_cmplx(k, ich)
  !   implicit none
  !   complex(DP), intent(in) :: k
  !   integer,  intent(in) :: ich
  !   complex(DP) :: g_i_cmplx
  !   ! Local variables
  !   complex(DP) :: omega_mes

  !   omega_mes = sqrt(k**2 + m_mes(ich)**2)
  !   g_i_cmplx = 1.0_DP/m_pi * k/sqrt(omega_mes) * u_k(k,Lambda(ich,1))

  ! end function g_i_cmplx

  ! function g_i_real(k, ich)
  !   implicit none
  !   real(DP), intent(in) :: k
  !   integer,  intent(in) :: ich
  !   real(DP) :: g_i_real

  !   ! omega_mes = sqrt(k**2 + m_mes(ich)**2)
  !   ! g_i_real = 1.0_DP/m_pi * k/sqrt(omega_mes) * u_k(k,Lambda(ich,1))
  !   g_i_real = real(g_i_cmplx(cmplx(k,0.0_DP,DP), ich), DP)

  ! end function g_i_real

end module potentials
