! //  function======================================================
module subprog
  implicit none
contains
! // bubble sort---------------------------------
  subroutine sort(Eigen)
    use params
    complex(kind(0d0)),intent(inout) :: Eigen(:)
    integer i,j
    complex(kind(0d0)) temp
    do i = 1, Nint
       do j = 1, Nint-i
          if(abs(Eigen(j))>abs(Eigen(j+1)))then
             temp = Eigen(j+1)
             Eigen(j+1) = Eigen(j)
             Eigen(j) = temp
          endif
       enddo
    enddo
  end subroutine sort
! // -------------------------------------------
! // impose periodic boundary condtion-----------------------------
  subroutine periodic(i,m)
    use params
    integer,intent(in)::i
    integer,intent(out):: m
    if (i==-1) then
       m = Nt-1
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
    complex(kind(0d0)),parameter::I = (0.0d0,1.0d0)
    A(1,1:2,3) = (/(0.0d0,0.0d0),(1.0d0,0.0d0)/)
    A(2,1:2,3) = (/(1.0d0,0.0d0),(0.0d0,0.0d0)/)
    A(1,1:2,2) = (/(0.0d0,0.0d0),(0.0d0,-1.0d0)/)
    A(2,1:2,2) = (/(0.0d0,1.0d0),(0.0d0,0.0d0)/)
    A(1,1:2,1) = (/(1.0d0,0.0d0),(0.0d0,0.0d0)/)
    A(2,1:2,1) = (/(0.0d0,0.0d0),(-1.0d0,0.0d0)/)
    A(1,1:2,4) = (/(1.0d0,0.0d0),(0.0d0,0.0d0)/)
    A(2,1:2,4) = (/(0.0d0,0.0d0),(1.0d0,0.0d0)/)
    if (mu.ne.4) then
       s = -I*A(alpha,beta,mu)
    else
       s = A(alpha,beta,mu)
    endif
  end function sgm
! --------------------------------------------------------
! // coordinate -------------------------------------------------------
  function x(mu,n1,n3,n4,m) result(y)
    use params
    integer,intent(in) :: mu,n1,n3,n4,m
    real(8) y
    if(mu==1) then
       y = -c*b*m + a*n4
    else if(mu==2) then
       y = -a*n1 + a*b*n3*m
    else if(mu==3) then
       y = -c - a*b*n4*m
    else if(mu==4) then
       y = -a*b*n1*m - a*n3
    endif
    y = y*dsqrt(2.0d0)/dble(1.0d0 + b**2*m**2)
  end function x
!--------------------------------------------------------------------------------------------------------------
! // Delivative term ------------------------------------------------------------------------------------------
  function del(mu,n1,n1p,n3,n3p,n4,n4p,m,mp) result(dc)
    use params
    integer,intent(in) :: mu,n1,n1p,n3,n3p,n4,n4p,m,mp
    real(8) dr
    complex(kind(0d0)) dc
    if(mu==1) then
       dr = - m*(delta2(m+1,mp)-delta2(m-1,mp))*delta2(n1,n1p)*delta2(n3,n3p)*delta2(n4,n4p) &
            - (a*n1 -x(2,n1,n3,n4,m)/dsqrt(2.0d0))*(delta2(n1+1,n1p)-delta2(n1-1,n1p))&
            *delta2(n3,n3p)*delta2(n4,n4p)*delta2(m,mp)/dble(a) &
            + (a*n3 -x(4,n1,n3,n4,m)/dsqrt(2.0d0))*(delta2(n3+1,n3p)-delta2(n3-1,n3p))&
            *delta2(n1,n1p)*delta2(n4,n4p)*delta2(m,mp)/dble(a) &
            + (a*n4 -x(3,n1,n3,n4,m)/dsqrt(2.0d0))*(delta2(n4+1,n4p)-delta2(n4-1,n4p))&
            *delta2(n1,n1p)*delta2(n3,n3p)*delta2(m,mp)/dble(a)
       dr = dr/x(1,n1,n3,n4,m)
    else if(mu==2) then
       dr = (delta2(n1+1,n1p)-delta2(n1-1,n1p))*delta2(n1,n1p)*delta2(n3,n3p)*delta2(n4,n4p)/dble(dsqrt(2.0d0)*a) &
            + (a*n3 -x(4,n1,n3,n4,m)/dsqrt(2.0d0))*(delta2(n3+1,n3p)-delta2(n3-1,n3p))&
            *delta2(n1,n1p)*delta2(n4,n4p)*delta2(m,mp)/dble(a*x(2,n1,n3,n4,m))
    else if(mu==3) then
       dr = - (delta2(m+1,mp)-delta2(m-1,mp))*delta2(n1,n1p)*delta2(n3,n3p)*delta2(n4,n4p)/dble(b) &
            - x(4,n1,n3,n4,m)*(delta2(n1+1,n1p)-delta2(n1-1,n1p))&
            *delta2(n3,n3p)*delta2(n4,n4p)*delta2(m,mp)/dble(a) &
            + x(2,n1,n3,n4,m)*(delta2(n3+1,n3p)-delta2(n3-1,n3p))&
            *delta2(n1,n1p)*delta2(n4,n4p)*delta2(m,mp)/dble(a) &
            - (c-dsqrt(2.0d0)*x(3,n1,n3,n4,m))*(delta2(n4+1,n4p)-delta2(n4-1,n4p))&
            *delta2(n1,n1p)*delta2(n3,n3p)*delta2(m,mp)/dble(a)
       dr = dr/dble(dsqrt(2.0d0)*x(1,n1,n3,n4,m))
    else if(mu==4) then
       dr = (a*n1-x(2,n1,n3,n4,m)/dsqrt(2.0d0))*(delta2(n1+1,n1p)-delta2(n1-1,n1p))&
            *delta2(n3,n3p)*delta2(n4,n4p)*delta2(m,mp)/dble(a*x(4,n1,n3,n4,m)) &
            +(delta2(n3+1,n3p)-delta2(n3-1,n3p))*delta2(n1,n1p)*delta2(n4,n4p)*delta2(m,mp)/dble(a*dsqrt(2.0d0)) 
    endif
    dc = cmplx(dr,0,kind(0d0))
  end function del
  ! ----------------------------------------------------------------------------------------------------------
  function jac(n1,n3,n4,m) result(g)
    use params
    integer,intent(in)::n1,n3,n4,m
    integer info,ipiv(Ndir),i
    real(8) g, mat(Ndir,Ndir)
    mat(:,:) = 0.0d0
    ipiv(:) = 0.0d0
    g = 1.0d0
    call matrix(n1,n3,n4,m,mat) ! set matrix
    call dgetrf(Ndir,Ndir,mat,Ndir,ipiv,info) ! LU decomposition by using LAPACK
    if(info /= 0)then
       write(*,*) "info =", info
       stop "LAPACK error!"
    endif
!//calculate the determinat
    do i = 1, Ndir
       g = g * mat(i,i)
       if(ipiv(i) /= i) then
          g = -g
       endif
    enddo
    g = g*(dsqrt(2.0d0)/(1.0d0+b**2*m**2))
    g = abs(g)
  end function jac
! ----------------------------------------------------------------------------------------------------------
  subroutine matrix(n1,n3,n4,m,J)
    use params
    integer,intent(in)::n1,n3,n4,m
    real(8),intent(inout):: J(:,:)
    J(1,1:4) =(/-b*m, -1.0d0, 0.0d0, -dsqrt(2.0d0)*b*m*x(0,n1,n3,n4,m)-a*n1/)
    J(2,1:4) =(/0.0d0, 0.0d0, 1.0d0, -dsqrt(2.0d0)*b*m*x(1,n1,n3,n4,m)-c/)
    J(3,1:4) =(/-1.0d0, b*m, 0.0d0,-dsqrt(2.0d0)*b*m*x(2,n1,n3,n4,m)+a*n3/)
    J(4,1:4) =(/0.0d0, 0.0d0, -b*m, -dsqrt(2.0d0)*b*m*x(3,n1,n3,n4,m)-a*n4/)
  end subroutine matrix
! ----------------------------------------------------------------------------------------------------------
end module subprog
