! //  global variable===============================================
module params
  implicit none
  integer,parameter :: Nspn = 2   ! number of Dirac indies
  integer,parameter :: Ndir = 4   ! number of lorentz indies
  integer,parameter :: Ns = 2   ! number of s length
  integer,parameter :: Nvol = Ns**4 ! numeber of lattice points
  integer,parameter :: Nint = Nvol*Nspn ! internal degree of freedom 
  real(8),parameter :: a = 1.0d0 ! lattice spacing of t direction
end module params
! ==================================================================