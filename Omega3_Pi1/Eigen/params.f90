! //  global variable===============================================
module params
  implicit none
  integer,parameter :: Nspn = 2   ! number of Dirac indies
  integer,parameter :: Ndir = 4   ! number of lorentz indies
  integer,parameter :: Nt = 4  ! number of t length
  integer,parameter :: Ns = 2  ! number of s length
  integer,parameter :: Nvol = Nt**3*Ns ! numeber of lattice points
  integer,parameter :: Nint = Nvol*Nspn ! internal degree of freedom 
  real(8),parameter :: a = 11.0d0 ! lattice spacing of t direction
  real(8),parameter :: b = 3.0d0 ! lattice spacing of s direction
  real(8),parameter :: c = 1.0d00  ! constant value of t2
end module params
! ==================================================================
