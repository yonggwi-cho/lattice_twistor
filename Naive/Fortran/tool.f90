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
       m = Ns-1
    else if(i==Ns) then
       m = 1
    else
       m = i
    endif
  end subroutine periodic
!------------------------------------------------------------------
! // impose dirichlet boundary condtion-----------------------------
  subroutine dirichlet(i,m) 
    use params
    integer,intent(in) :: i
    integer,intent(out) :: m
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
  function delta(i,j) result(n)
    use params
    integer,intent(in) :: i,j
    integer k,l
    real(8) n
    if(i == j)then
       n = 1.0d0
    else if(i /= j)then 
       n = 0.0d0
    else
       stop "something is wrong!"
    endif
  end function delta
!---------------------------------------------------------------------
! naive fermion------------------------------------------------------
  function naive(i,j,alpha,beta,mu) result(d)
    use params
    integer,intent(in) :: i,j,alpha,beta,mu
    integer k,l,m
    complex(kind(0d0)) d
    call periodic(i,k)
    call periodic(j,l)
    call periodic(j-1,m)
    d = (delta(k,l)-delta(k,m))*gamma(mu,alpha,beta)
  end function naive
!---------------------------------------------------------------
! // klonecker delta -------------------------------------------
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
  !-------------------------------------------------------------------------
  !--- euclidian gamma matrix ----------------------------------------
  function gamma(mu,a,b)  result (d)
    integer,intent(in) :: a,b,mu
    complex(kind(0d0)) G(4,4,4)
    complex(kind(0d0)) d
    !--mu=1-----------
    G(1,1:4,1) = (/(0.0d0,0.0d0),(0.0d0,0.0d0),(0.0d0,0.0d0),(0.0d0,-1.0d0)/)
    G(2,1:4,1) = (/(0.0d0,0.0d0),(0.0d0,0.0d0),(0.0d0,-1.0d0),(0.0d0,0.0d0)/)
    G(3,1:4,1) = (/(0.0d0,0.0d0),(0.0d0,1.0d0),(0.0d0,0.0d0),(0.0d0,0.0d0)/)
    G(4,1:4,1) = (/(0.0d0,1.0d0),(0.0d0,0.0d0),(0.0d0,0.0d0),(0.0d0,0.0d0)/)
    !--mu=2-----------
    G(1,1:4,2) = (/(0.0d0,0.0d0),(0.0d0,0.0d0),(0.0d0,0.0d0),(-1.0d0,0.0d0)/)
    G(2,1:4,2) = (/(0.0d0,0.0d0),(0.0d0,0.0d0),(1.0d0,0.0d0),(0.0d0,0.0d0)/)
    G(3,1:4,2) = (/(0.0d0,0.0d0),(1.0d0,0.0d0),(0.0d0,0.0d0),(0.0d0,0.0d0)/)
    G(4,1:4,2) = (/(-1.0d0,0.0d0),(0.0d0,0.0d0),(0.0d0,0.0d0),(0.0d0,0.0d0)/)
    !--mu=3-----------
    G(1,1:4,3) = (/(0.0d0,0.0d0),(0.0d0,0.0d0),(0.0d0,-1.0d0),(0.0d0,0.0d0)/)
    G(2,1:4,3) = (/(0.0d0,0.0d0),(0.0d0,0.0d0),(0.0d0,0.0d0),(0.0d0,1.0d0)/)
    G(3,1:4,3) = (/(0.0d0,1.0d0),(0.0d0,0.0d0),(0.0d0,0.0d0),(0.0d0,0.0d0)/)
    G(4,1:4,3) = (/(0.0d0,0.0d0),(0.0d0,-1.0d0),(0.0d0,0.0d0),(0.0d0,0.0d0)/)
    !--mu=4-----------
    G(1,1:4,4) = (/(0.0d0,0.0d0),(0.0d0,0.0d0),(1.0d0,0.0d0),(0.0d0,0.0d0)/)
    G(2,1:4,4) = (/(0.0d0,0.0d0),(0.0d0,0.0d0),(0.0d0,0.0d0),(1.0d0,0.0d0)/)
    G(3,1:4,4) = (/(1.0d0,0.0d0),(0.0d0,0.0d0),(0.0d0,0.0d0),(0.0d0,0.0d0)/)
    G(4,1:4,4) = (/(0.0d0,0.0d0),(1.0d0,0.0d0),(0.0d0,0.0d0),(0.0d0,0.0d0)/)
    !----------------
    d = G(a,b,mu)
  end function gamma
!---------------------------------------------------------------------
! // pauli matrix----------------------------------------------------------
  function sgm(mu,alpha,beta) result(s)
    integer,intent(in)::mu,alpha,beta
    complex(kind(0d0)) A(2,2,4)
    complex(kind(0d0)) s
    complex(kind(0d0)),parameter::I = (0.0d0,1.0d0)
    A(1,1:2,1) = (/(0.0d0,0.0d0),(1.0d0,0.0d0)/)
    A(2,1:2,1) = (/(1.0d0,0.0d0),(0.0d0,0.0d0)/)
    A(1,1:2,2) = (/(0.0d0,0.0d0),(0.0d0,-1.0d0)/)
    A(2,1:2,2) = (/(0.0d0,1.0d0),(0.0d0,0.0d0)/)
    A(1,1:2,3) = (/(1.0d0,0.0d0),(0.0d0,0.0d0)/)
    A(2,1:2,3) = (/(0.0d0,0.0d0),(-1.0d0,0.0d0)/)
    A(1,1:2,4) = (/(1.0d0,0.0d0),(0.0d0,0.0d0)/)
    A(2,1:2,4) = (/(0.0d0,0.0d0),(1.0d0,0.0d0)/)
    if (mu.ne.4) then
       s = -I*A(alpha,beta,mu)
    else
       s = A(alpha,beta,mu)
    endif
  end function sgm
! ---------------------------------------------------------------------
  function dis(n1,n2,n3,n4) result(r)
    use params
    integer,intent(in)::n1,n2,n3,n4
    real(8) r
    r = n1**2 + n2**2 + n3**2 + n4**2
  end function dis
! // Delivative term --------------------------------------------------------------------------------------------------
  function del(mu,n1,n1p,n2,n2p,n3,n3p,n4,n4p) result(dc)
    use params
    integer,intent(in) :: mu,n1,n1p,n2,n2p,n3,n3p,n4,n4p
    real(8) dr,dp
    complex(kind(0d0)) dc
    if(mu==1) then
       dp = n1*n1*(delta2(n1+1,n1p)-delta2(n1-1,n1p))*delta2(n2,n2p)*delta2(n3,n3p)*delta2(n4,n4p) &
            +  n1*n2*(delta2(n2+1,n2p)-delta2(n2-1,n2p))*delta2(n1,n1p)*delta2(n3,n3p)*delta2(n4,n4p) &
            +  n1*n3*(delta2(n3+1,n3p)-delta2(n3-1,n3p))*delta2(n1,n1p)*delta2(n2,n2p)*delta2(n4,n4p) &
            +  n1*n4*(delta2(n4+1,n4p)-delta2(n4-1,n4p))*delta2(n1,n1p)*delta2(n2,n2p)*delta2(n3,n3p)
       dr = dis(n1,n2,n3,n4)*(delta2(n1+1,n1p)-delta2(n1-1,n1p))*delta2(n2,n2p)*delta2(n3,n3p)*delta2(n4,n4p) - 2.0*dp
    else if(mu==2) then
       dp =  n2*n1*(delta2(n1+1,n1p)-delta2(n1-1,n1p))*delta2(n2,n2p)*delta2(n3,n3p)*delta2(n4,n4p) &
            +  n2*n2*(delta2(n2+1,n2p)-delta2(n2-1,n2p))*delta2(n1,n1p)*delta2(n3,n3p)*delta2(n4,n4p) & 
            +  n2*n3*(delta2(n3+1,n3p)-delta2(n3-1,n3p))*delta2(n1,n1p)*delta2(n2,n2p)*delta2(n4,n4p) & 
            +  n2*n4*(delta2(n4+1,n4p)-delta2(n4-1,n4p))*delta2(n1,n1p)*delta2(n2,n2p)*delta2(n3,n3p)
       dr = dis(n1,n2,n3,n4)*(delta2(n2+1,n2p)-delta2(n2-1,n2p))*delta2(n1,n1p)*delta2(n3,n3p)*delta2(n4,n4p) - 2.0*dp
    else if(mu==3) then
       dp=  n3*n1*(delta2(n1+1,n1p)-delta2(n1-1,n1p))*delta2(n2,n2p)*delta2(n3,n3p)*delta2(n4,n4p) &
            +  n3*n2*(delta2(n2+1,n2p)-delta2(n2-1,n2p))*delta2(n1,n1p)*delta2(n3,n3p)*delta2(n4,n4p) &
            +  n3*n3*(delta2(n3+1,n3p)-delta2(n3-1,n3p))*delta2(n1,n1p)*delta2(n2,n2p)*delta2(n4,n4p) &
            +  n3*n4*(delta2(n4+1,n4p)-delta2(n4-1,n4p))*delta2(n1,n1p)*delta2(n2,n2p)*delta2(n3,n3p)
       dr = dis(n1,n2,n3,n4)*(delta2(n3+1,n3p)-delta2(n3-1,n3p))*delta2(n1,n1p)*delta2(n2,n2p)*delta2(n4,n4p) - 2.0*dp
    else if(mu==4) then
       dp =  n4*n1*(delta2(n1+1,n1p)-delta2(n1-1,n1p))*delta2(n2,n2p)*delta2(n3,n3p)*delta2(n4,n4p) &
            +  n4*n2*(delta2(n2+1,n2p)-delta2(n2-1,n2p))*delta2(n1,n1p)*delta2(n3,n3p)*delta2(n4,n4p) &
            +  n4*n3*(delta2(n3+1,n3p)-delta2(n3-1,n3p))*delta2(n1,n1p)*delta2(n2,n2p)*delta2(n4,n4p) &
            +  n4*n4*(delta2(n4+1,n4p)-delta2(n4-1,n4p))*delta2(n1,n1p)*delta2(n2,n2p)*delta2(n3,n3p)
       dr = dis(n1,n2,n3,n4)*(delta2(n4+1,n4p)-delta2(n4-1,n4p))*delta2(n1,n1p)*delta2(n2,n2p)*delta2(n3,n3p) - 2.0*dp
    endif
    dr = a*dr
    dc = cmplx(dr,0,kind(0d0))    
  end function del
!-----------------------------------------------------------------------------------------------------------------------
  function jac(n1,n2,n3,n4) result(g)
    use params
    integer,intent(in)::n1,n2,n3,n4
    integer info,ipiv(Ndir),i
    real(8) g, mat(Ndir,Ndir)
!    mat(:,:) = 0.0d0
!    ipiv(:) = 0.0d0
    g = 1.0d0
!    call matrix(n1,n2,n3,n4,mat) ! set matrix
!    call dgetrf(Ndir,Ndir,mat,Ndir,ipiv,info) ! LU decomposition by using LAPACK
!    if(info /= 0)then
!       write(*,*) "info =", info
!       stop "LAPACK error!"
!    endif
!    write(*,*) mat(:,:)
!    stop
!//calculate the determinat
!    do i = 1, Ndir
!       g = g * mat(i,i)
!       if(ipiv(i) /= i) then
!          g = -g
!       endif
!    enddo
    g = 4*3*2*1
    g = g/(dis(n1,n2,n3,n4)**4)
    g = abs(g)
  end function jac
! ----------------------------------------------------------------------------------------------------------
  subroutine matrix(n1,n2,n3,n4,mat)
    use params
    integer,intent(in)::n1,n2,n3,n4
    real(8),intent(inout) :: mat(:,:)
    integer i,j
    do i = 1, Ndir
       do j = 1, Ndir
          mat(i,j)  = delta(i,j)-2.0*nn(i,j,n1,n2,n3,n4)/dis(n1,n2,n3,n4)
       enddo
    enddo
  end subroutine matrix
! ----------------------------------------------------------------------------------------------------------
  function nn(mu,nu,n1,n2,n3,n4) result(r)
    integer,intent(in):: mu,nu,n1,n2,n3,n4
    real(8) r
    if(mu == nu)then
       if(mu == 1) then
          r = n1*n1
       else if(mu == 2) then
          r = n2*n2
       else if(mu == 3) then
          r = n3*n3
       else if(mu == 4) then
          r = n4*n4
       end if
    else if((mu == 1 .and. nu ==2) .or. (mu ==2 .and. nu == 1)) then
       r = n1*n2
    else if((mu == 1 .and. nu ==3) .or. (mu ==3 .and. nu == 1)) then
       r = n1*n3
    else if((mu == 1 .and. nu ==4) .or. (mu ==4 .and. nu == 1)) then
       r = n1*n4
    else if((mu == 2 .and. nu ==3) .or. (mu ==3 .and. nu == 2)) then
       r = n2*n3
    else if((mu == 2 .and. nu ==4) .or. (mu ==4 .and. nu == 2)) then
       r = n2*n4
    else if((mu == 3 .and. nu ==4) .or. (mu ==4 .and. nu == 3)) then
       r = n3*n4
    else
       stop "something is wrong 1!"
    endif
  end function nn
! ----------------------------------------------------------------------------------------------------------
end module subprog
