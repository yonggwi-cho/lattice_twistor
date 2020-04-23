!=========================================================================================================
! // naive Fermions
! // Solve eigen values for a Dirac operator              
!=========================================================================================================
program main
  use params
  use subprog
  implicit none
  ! Parameters
  integer,parameter :: nb=64, nmax=Nint
  integer,parameter :: lda=nmax, ldvr=nmax, lwork=(1+nb)*nmax
  ! Local Scalars
  Integer i,j,ifail,info,lwkopt
  integer  mu,n1,n1p,n2,n2p,n3,n3p,n4,n4p,alpha,beta
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
  write(output,12) Ns
12 format("./Output/eigen_Ns_",i2.2,".d")
  open(10,file=output,status="replace")
  ! set Dirac operator
  write(*,*)"================================================================================================="
  write(10,*)"================================================================================================"
  write(*,*)" eigensolver of the Dirac operator (Naive fermion)"
  write(*,*)"use LAPACK!"
  write(*,*)"Ns :",Ns
  write(10,*)"Ns :",Ns
  write(*,'(a,f5.2)')"Lattice spacing a : ",a
  write(10,'(a,f5.2)')"Lattice spacing a : ",a
  write(*,*)"================================================================================================="
  write(10,*)"================================================================================================"
  call cpu_time(time1)
  write(*,*) "Start making a Dirac operator!"
  do alpha = 1, Nspn
     do beta = 1, Nspn
        write(*,150) alpha,beta
        do n1 = 1, Ns
           do n1p = 1, Ns
              do n2 = 1, Ns
                 do n2p = 1, Ns
                    do n3 = 1, Ns
                       do n3p = 1, Ns
                          do n4 = 1, Ns
                             do n4p = 1, Ns
                                do mu = 1, Ndir
                                   i =  alpha + Nspn*(n1-1) + Nspn*Ns*(n2-1) + Nspn*Ns**2*(n3-1) + Nspn*Ns**3*(n4-1)
                                   j =  beta  + Nspn*(n1p-1) + Nspn*Ns*(n2p-1) + Nspn*Ns**2*(n3p-1) + Nspn*Ns**3*(n4p-1)
                                   d = naive(i,j,alpha,beta,mu)
                                   Dirac(i,j) = Dirac(i,j) + d
                                   if((abs(Dirac(i,j))*0.0d0)/=0.0d0) then
                                      write(*,140)mu,alpha,n1,n2,n3,n4,Dirac(i,j)
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
  write(10,*) " eigenvalue(Real)       |eigenvalue(imag)       |     absolute value "
  if(info == 0) then
     do j = 1, Nint
        write(10, *) Real(Eigen(j)),aimag(Eigen(j)),abs(Eigen(j)),j
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
