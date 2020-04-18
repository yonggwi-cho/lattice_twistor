! //  global variable===============================================
module params
  implicit none
  integer,parameter :: Nspn = 2   ! number of Dirac indies
  integer,parameter :: Ndir = 4   ! number of lorentz indies
  integer,parameter :: Nt = 2  ! number of t length+1
  integer,parameter :: Ns = 3  ! number of s length
  integer,parameter :: Nvol = Nt**2*Ns**2 ! numeber of lattice points
  integer,parameter :: Nint = Nvol*Nspn ! internal degree of freedom 
  real(8),parameter :: a = 11.0d0 ! lattice spacing of t direction
  real(8),parameter :: b = 2.3d-01 ! lattice spacing of s direction
  real(8),parameter :: c1 = 1.0d00  ! constant value of t2
  real(8),parameter :: c2 = 1.0d00  ! constant value of t4 
end module params
! ==================================================================
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
! --------------------------------------------------------
! // coordinate -------------------------------------------------------
  function x(mu,n1,n3,m1,m2) result(y)
    use params
    integer,intent(in) :: mu,n1,n3,m1,m2
    real(8) y
    if(mu==1) then
       y = (b*c1*m1 - a*b*n1*m2 - c2)   /dble(1.0d0+b**2*(m1**2+m2**2))
    else if(mu==2) then
       y = (a*n1 - a*b*n3*m1 + b*c2*m2) /dble(1.0d0+b**2*(m1**2+m2**2))
    else if(mu==3) then
       y = (c1 - a*b*n3*m2 + b*c2*m1)   /dble(1.0d0+b**2*(m1**2+m2**2))
    else if(mu==4) then
       y = (a*b*n1*m1 + b*c1*m2 + a*n3) /dble(1.0d0+b**2*(m1**2+m2**2))
    endif
  end function x
! ----------------------------------------------------------------------
! // Delivative term ------------------------------------------------------------------------------------------
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
! ----------------------------------------------------------------------------------------------------------
end module subprog
!=========================================================================================================
! // Fermions on the twister lattice space       
! // Solve Dirac operator                
! // impose periodic boundary condtion on s  
! // impose dirichlet boundary condtion on t  
! // we use weyl fermion                        
! // t2, t4 set to constant
! // lattice spacing is (a,b)
!=========================================================================================================
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
  character(100) output
  ! Local Arrays
  Complex (kind(0D0)) dummy(1, 1),work(lwork)
  Double Precision rwork(2*nmax)
  complex(kind(0d0)),allocatable :: Dirac(:,:),Eigen(:)
  complex(kind(0d0)) d
  allocate(Dirac(Nint,Nint))
  ! zero clear
  Dirac(:,:) = 0.0d0
  write(output,12) Nt,Ns
12 format("./Output/eigen_Nt_",i2.2,"_Ns_",i2.2,".d")
  open(10,file=output,status="replace")
  ! set Dirac operator
  write(*,*)"==================================================================================================="
  write(10,*)"==================================================================================================="
  write(*,*)"eigensolver of the Dirac operator(Weyl fermion) on the twister-like lattice space"
  write(10,*)"eigensolver of the Dirac operator(Weyl fermion) on the twister-like lattice space"
  write(*,*)"use LAPACK!"
  write(*,*)"Nt :",Nt
  write(10,*)"Nt :",Nt
  write(*,*)"Ns :",Ns
  write(10,*)"Ns :",Ns
  write(*,'(a,f5.2)')"Lattice spacing a : ",a
  write(10,'(a,f5.2)')"Lattice spacing a : ",a
  write(*,'(a,f5.2)')"Lattice spacing b : ",b
  write(10,'(a,f5.2)')"Lattice spacing b  :",b
  write(*,'(a,f5.2)')"constant value of t2,t4 : ",c1
  write(10,'(a,f5.2)')"constant value of t2,t4 : ",c1
  write(*,*)"==================================================================================================="
  write(10,*)"==================================================================================================="
  call cpu_time(time1)
  write(*,*) "Start making a Dirac operator!"
  do alpha = 1, Nspn
  do beta = 1, Nspn
     write(*,150) alpha,beta
     do n1 = 1, Nt
     do n1p = 1, Nt
        do n3 = 1, Nt
        do n3p = 1, Nt
           do m1 = 1, Ns
           do m1p = 1, Ns
              do m2 = 1, Ns
              do m2p = 1, Ns
                 do mu = 1, Ndir
                    i =  alpha + Nspn*(n1-1) + Nspn*Nt*(n3-1) + Nspn*Nt**2*(m1-1) + Nspn*Nt**2*Ns*(m2-1) 
                    j =  beta  + Nspn*(n1p-1) + Nspn*Nt*(n3p-1) + Nspn*Nt**2*(m1p-1) + Nspn*Nt**2*Ns*(m2p-1)
                    d = sgm(mu,alpha,beta)*del(mu,n1-1,n1p-1,n3-1,n3p-1,m1,m1p,m2,m2p)
                    Dirac(i,j) = Dirac(i,j) + d
                    if((abs(Dirac(i,j))*0.0d0)/=0.0d0) then
                       write(*,140)mu,alpha,n1-1,n3-1,m1,m2,Dirac(i,j)
                       stop 
                    endif
                 enddo
              enddo
              enddo
           enddo
           enddo
        enddo
        enddo
     enddo
     enddo
  enddo
  enddo
  write(*,*) "Finished making a Dirac operator!"
  write(*,*) "Start solving a Dirac operator!"
  ! Compute the eigenvalues by using the lapack
  allocate(Eigen(Nint))
  Eigen(:) = 0.0d0
  Call zgeev('N','N',Nint,Dirac,Nint,Eigen,dummy,1,dummy,1,work,lwork,rwork,info)
  deallocate(Dirac)
  call cpu_time(time2)
  write(*,*) "Finished solving a Dirac operator!"
  lwkopt = work(1)
  ! sort eigenvlues
  call sort(Eigen)
  !  Print solution
  write(10,*) " arg(eingenvalue)      |     absolute value     |       eigenvalue"
  if(info == 0) then
     do j = 1, Nint
        write(10, *) atan(aimag(Eigen(j))/real(Eigen(j))),abs(Eigen(j)),Eigen(j)
     enddo
     ifail = 0
  else
     Write (*, *)
     Write (*, 110) 'Failure in ZGEEV.  INFO = ', info
  endif
  ! Print workspace information
  If (lwork<lwkopt) Then
     Write (*, *)
     Write (*, 120) 'Optimum workspace required = ', lwkopt, &
          'Workspace provided         = ', lwork
  EndIf
  close(10)
  deallocate(Eigen)
  write(*,130) time2-time1
100 Format ((3X,4(' (',F7.4,',',F7.4,')',:)))
110 Format (1X, A, I4)
120 Format (1X, A, I5, /1X, A, I5)
130 format("Solving time [sec]=",f15.4)
140 format("(mu,alpha,n1,n3,m1,m2,Dirac(i,j))=(",i1.1,","i1.1,",",i2.2,",",i2.2,",",i2.2,",",i2.2,",(",f3.3,",",f3.3"))")
150  format("(alpha,beta)=(",i1.1,",",i1.1,")")
End Program main
!=============================================================================================================
