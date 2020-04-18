! //  global variable===============================================
module params
  implicit none
  integer,parameter :: Nspn = 2   ! number of Dirac indies
  integer,parameter :: Ndir = 4   ! number of lorentz indies
  integer,parameter :: Nt = 4  ! number of t spacetime's length
  integer,parameter :: Ns = 4  ! number of s spacetime's length
  integer,parameter :: Nvol = Nt**2*Ns**2 ! numeber of lattice points
  integer,parameter :: Nint = Nvol*Nspn ! internal degree of freedom 
  real(8),parameter :: a = 1.0d0 ! lattice spacing of t direction
  real(8),parameter :: b = 2.0d0 ! lattice spacing of s direction
  real(8),parameter :: c1 = 2.0d0  ! constant value of t2
  real(8),parameter :: c2 = 1.0d0  ! constant value of t4 
  real(8),parameter :: PI = 6.0d0*asin(0.5d0) ! PI 
end module params
! ==================================================================
! //  function======================================================
module subprog
  implicit none
contains
! // impose periodic boundary condtion-----------------------------
  subroutine periodic(i,m)
    use params
    integer,intent(in)::i
    integer,intent(out):: m
    if (i==-1) then
       m = Nt
    else if(i==Nt) then
       m = 1
    else
       m = i
    endif
  end subroutine periodic
!------------------------------------------------------------------
! // impose dirichlet boundary condtion-----------------------------
  subroutine dirichlet(i,m) 
    use params
    integer,intent(in)::i
    integer,intent(out):: m
    if(i==0)then
       m = -1000
    else if(i==Ns+1) then
       m = -10000
    else
       m = i
    endif
  end subroutine dirichlet
!------------------------------------------------------------------
! // klonecker delta -------------
  function delta1(i,j) result(n)
    use params
    integer,intent(in) :: i,j
    integer k,l
    real(8) n
    call periodic(i,k)
    call periodic(j,l)
    if(k==l)then
       n = 1.0d0
    else if(k.ne.l)then 
       n = 0.0d0
    else
       stop "something is wrong!"
    endif
  end function delta1
!--------------------------------
! // klonecker delta -------------
  function delta2(i,j) result(n)
    use params
    integer,intent(in) :: i,j
    integer k,l
    real(8) n
    call dirichlet(i,k)
    call dirichlet(j,l)
    if(k==l)then
       n = 1.0d0
    else if(k.ne.l)then 
       n = 0.0d0
    else
       stop "something is wrong!"
    endif
  end function delta2
!--------------------------------
! // pauli matrix-----------------------------------------
  function sgm(mu,alpha,beta) result(s)
    integer,intent(in)::mu,alpha,beta
    complex(kind(0d0)) A(2,2,4)
    complex(kind(0d0)) s
    A(1,1:2,1) = (/(0.0d0,0.0d0),(1.0d0,0.0d0)/)
    A(2,1:2,1) = (/(1.0d0,0.0d0),(0.0d0,0.0d0)/)
    A(1,1:2,2) = (/(0.0d0,0.0d0),(0.0d0,-1.0d0)/)
    A(2,1:2,2) = (/(0.0d0,1.0d0),(0.0d0,0.0d0)/)
    A(1,1:2,3) = (/(1.0d0,0.0d0),(0.0d0,0.0d0)/)
    A(2,1:2,3) = (/(0.0d0,0.0d0),(-1.0d0,0.0d0)/)
    A(1,1:2,4) = (/(1.0d0,0.0d0),(0.0d0,0.0d0)/)
    A(2,1:2,4) = (/(0.0d0,0.0d0),(1.0d0,0.0d0)/)
    s = A(alpha,beta,mu)
  end function sgm
! --------------------------------------------------------
! // coordinate ------------------------------------------
  function x(mu,n1,n3,m1,m2) result(y)
    use params
    integer,intent(in) :: mu,n1,n3,m1,m2
    real(8) y
    if(mu==1) then
       y = (b*c1*m1 - a*b*n1*m2 -  c2)   /dble(1.0d0+b**2*(m1**2+m2**2))
    else if(mu==2) then
       y = (a*n1 - a*b*n3*m1 + b*c2*m2)/ dble(1.0d0+b**2*(m1**2+m2**2))
    else if(mu==3) then
       y = (c1 - a*b*n3*m2 + b*c2*m1)   /dble(1.0d0+b**2*(m1**2+m2**2))
    else if(mu==4) then
       y = (a*b*n1*m1 + b*c1*m2 + a*n3)/ dble(1.0d0+b**2*(m1**2+m2**2))
    endif
!    write(*,*) y
  end function x
! --------------------------------------------------------
! // Delivative term --------------------------------------------------------------------------------------
  function del(mu,n1,n1p,n3,n3p,m1,m1p,m2,m2p) result(dc)
    use params
    integer,intent(in) :: mu,n1,n1p,n3,n3p,m1,m1p,m2,m2p
    real(8) dr
    complex(kind(0d0)) dc
    if(mu==1) then
       dr = -(delta1(n1+1,n1p)-delta1(n1-1,n1p))*delta1(n3,n3p)*delta2(m1,m1p)*delta2(m2,m2p)/dble(2.0d0*a*b*m2) &
            + (delta2(m1+1,m1p)-delta2(m1-1,m1p))*delta1(n3,n3p)*delta1(n1,n1p)*delta2(m2,m2p) &
            /(dble(-2.0d0*b*m1*x(1,n1,n3,m1,m2)+c1)*dble(2.0d0*b)) &
            + (delta2(m2+1,m2p)-delta2(m2-1,m2p))*delta1(n3,n3p)*delta1(n1,n1p)*delta2(m1,m1p) &
            /(dble(-2.0d0*b*m2*x(1,n1,n3,m1,m2)-a*n1)*dble(2.0d0*b))
    else if(mu==2) then
       dr = (delta1(n1+1,n1p)-delta1(n1-1,n1p))*delta1(n3,n3p)*delta2(m1,m1p)*delta2(m2,m2p)/dble(2.0d0*a) &
            - (delta1(n3+1,n3p)-delta1(n3-1,n3p))*delta1(n1,n1p)*delta2(m1,m1p)*delta2(m2,m2p)/dble(2.0d0*m1*a*b)&
            + (delta2(m1+1,m1p)-delta2(m1-1,m1p))*delta1(n3,n3p)*delta1(n1,n1p)*delta2(m2,m2p) &
            /(dble(-2.0d0*b*m1*x(2,n1,n3,m1,m2)-a*n3)*dble(2.0d0*b)) &
            + (delta2(m2+1,m2p)-delta2(m2-1,m2p))*delta1(n3,n3p)*delta1(n1,n1p)*delta2(m1,m1p) &
            /(dble(-2.0d0*b*m2*x(2,n1,n3,m1,m2)-c2)*dble(2.0d0*b))
    else if(mu==3) then
       dr = -(delta1(n3+1,n3p)-delta1(n3-1,n3p))*delta1(n1,n1p)*delta2(m1,m1p)*delta2(m2,m2p)/dble(2.0d0*a*b*m2) &
            + (delta2(m1+1,m1p)-delta2(m1-1,m1p))*delta1(n3,n3p)*delta1(n1,n1p)*delta2(m2,m2p) &
            /(dble(-2.0d0*b*m1*x(3,n1,n3,m1,m2)+c2)*dble(2.0d0*b)) &
            + (delta2(m2+1,m2p)-delta2(m2-1,m2p))*delta1(n3,n3p)*delta1(n1,n1p)*delta2(m1,m1p) &
            /(dble(-2.0d0*b*m2*x(3,n1,n3,m1,m2)-a*n3)*dble(2.0d0*b))
    else if(mu==4) then
       dr = (delta1(n1+1,n1p)-delta1(n1-1,n1p))*delta1(n3,n3p)*delta2(m1,m1p)*delta2(m2,m2p)/dble(2.0d0**a*b*m1) &
           + (delta1(n3+1,n3p)-delta1(n3-1,n3p))*delta1(n1,n1p)*delta2(m1,m1p)*delta2(m2,m2p)/dble(2.0d0*b) &
           + (delta2(m1+1,m1p)-delta2(m1-1,m1p))*delta1(n3,n3p)*delta1(n1,n1p)*delta2(m2,m2p) &
           /(dble(-2.0d0*b*m1*x(4,n1,n3,m1,m2)+a*n1)*dble(2.0d0*b)) &
           + (delta2(m2+1,m2p)-delta2(m2-1,m2p))*delta1(n3,n3p)*delta1(n1,n1p)*delta2(m1,m1p) &
           /(dble(-2.0d0*b*m2*x(4,n1,n3,m1,m2)+c1)*dble(2.0d0*b))
    endif
    dr = dble(1+b**2*(m1**2+m2**2))*dr/4.0d0
    dc = cmplx(dr,0,kind(0d0))
  end function del
! ---------------------------------------------------------------------------------------------------------
end module subprog
!=============================================================================================================
! // Fermions on the twister lattice space       
! // Solve Dirac operator                
! // impose periodic boundary condtion (s,t)  
! // we use weyl fermion                        
! // t2, t4 set to constant
! // lattice spacing is 1 (a,b)
!=============================================================================================================
program main
  use params
  use subprog
  implicit none
  ! Parameters
  Integer nb, nmax
  Parameter (nb=64, nmax=Nint)
  Integer lda, ldvr, lwork
  Parameter (lda=nmax, ldvr=nmax, lwork=(1+nb)*nmax)
  ! Local Scalars
  Integer i,j,ifail,info,lwkopt
  integer  mu,n1,n1p,n3,n3p,m1,m1p,m2,m2p,alpha,beta
  real(8) time1,time2
  ! Local Arrays
  Complex (kind(0D0)) dummy(1, 1),work(lwork)
  Double Precision rwork(2*nmax)
  ! set Dirac operator
  open(100,file="./coordinate.d")
  do n1 = 1, Nt
!     write(100,*)x(1,0,0,m1,1),x(2,0,0,m1,1),m1
     write(100,*) x(1,n1-1,0,1,1),x(2,n1-1,0,1,1),n1-1
  enddo
  close(100)
End Program main
!=============================================================================================================
