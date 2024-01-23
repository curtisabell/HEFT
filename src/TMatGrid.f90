program TMatGrid
    use kinds
    use numFort
    use heft
    use SMatrix
    implicit none

    integer :: i_re, i_im ! loop vars

    ! corners of the grid
    real(DP) :: E_real, E_real_min, E_real_max
    real(DP) :: E_imag, E_imag_min, E_imag_max
    integer :: nPoints = 128 ! points in each dim

    complex(DP) :: E
    real(DP), dimension(2) :: E_2d
    real(DP) :: invQuantity ! Tmat for 0b, A for nb

    character(len=:), allocatable :: fileName_grid

    call initialiseHEFT()
    call printCurrentParameters(iParamChoice)
    call initialiseHEFTInfinite()


    if (IamRoot) then
        allocate( lowBound(n_ch), midBound(n_ch) )
        lowBound(:) = m_mes(:) + m_bar(:)
        absErr = 1.0E-6_DP
        relErr = 0.0_DP


        ! TODO: generalise/config file
        ! E_real_min = (m_mes(ch_onshell) + m_bar(ch_onshell))*1.01_DP
        E_real_min = 1.4_DP
        E_imag_min = -0.01_DP
        E_real_max = 1.8_DP
        E_imag_max = -0.1_DP

        fileName_grid = 'data/TMatGrid_fit' &
            & // trim(int2str(iParamChoice)) // '.out'
        open(101, file=trim(fileName_grid), action='write')


        do i_im = 0, nPoints
           E_imag = E_imag_min + (E_imag_max-E_imag_min) &
               & / real(nPoints-1,DP) * real(i_im-1,DP)
           do i_re = 0, nPoints
              E_real = E_real_min + (E_real_max-E_real_min) &
                  & / real(nPoints-1,DP) * real(i_re-1,DP)

              if (i_im.eq.0) then
                  ! write the energy values on the first line
                  if (i_re.eq.0) then
                      write(101, '(f16.11)', advance='no') 0.0_DP
                  else
                      write(101, '(f16.11)', advance='no') E_real
                  end if

              else
                  if (i_re.eq.0) then
                      write(101, '(f16.11)', advance='no') E_imag
                  else
                      E = cmplx(E_real, E_imag, DP)
                      E_2d(:) = [E_real, E_imag]
                      call calcSMatrixPole(2, E_2d, invQuantity)
                      write(101, '(f16.11)', advance='no') invQuantity
                  end if
              end if

           end do
           write(101,*)

           write(*,'(f9.4,4x,i0,a,i0)') E_imag, i_im, '/', nPoints

        end do

        close(101)
        deallocate( lowBound, midBound )
    end if

    call finaliseHEFT()
    !-----------------------------------------------------------------------!

end program TMatGrid
