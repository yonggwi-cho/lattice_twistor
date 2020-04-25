! //  global variable===============================================
module params
  implicit none
  integer,parameter :: Nspn = 4   ! number of Dirac indies
  integer,parameter :: Ndir = 4   ! number of lorentz indies
  integer,parameter :: Ns = 8 ! number of s length
  integer,parameter :: Nvol = Ns**4  ! vox size
  integer,parameter :: Nint = Nvol*Nspn   ! internal degree of freedom
  real(8),parameter :: a = 1.0d0
end module params
! ==================================================================
