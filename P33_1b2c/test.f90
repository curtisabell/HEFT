program test
  use kinds
  use heft
  implicit none

  ! ----------------------Setup coarray variables---------------------
  nImg = num_images()
  iImg = this_image()
  IamRoot = (iImg .eq. root_image)
  doCA = (num_images() .gt. 1)

  ! if (IamRoot) call initialiseHEFT()
  call initialiseHEFT()

  sync all

  write(*,*) iImg, gBare


  if (iImg.eq.3) then
     write(*,*) iImg, g_choice, f_choice, u_choice
  end if


  call finaliseHEFT()

contains

end program test
