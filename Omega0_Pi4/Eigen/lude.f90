program main
  implicit none
  integer,parameter :: Ndir = 3 
  integer info,ipiv(Ndir),i
  real(8) g, mat(Ndir,Ndir)
  mat(:,:) = 0.0d0
  ipiv(:) = 0.0d0
  g = 1.0d0
  mat(1:3,1) = (/4,1,1/)
  mat(1:3,2) = (/1,3,1/)
  mat(1:3,3) = (/2,1,5/)
  write(*,*) "A = "
  write(*,*) mat(:,:)
  call dgetrf(Ndir,Ndir,mat,Ndir,ipiv,info) ! LU decomposition by using LAPACK
  if(info /= 0)then
     write(*,*) "info =", info
     stop "LAPACK error!"
  endif
  write(*,*) "L,U = "
  write(*,*) mat(:,:)
  do i = 1, Ndir
     g = g * mat(i,i)
     write(*,*) ipiv(i),i
     if(ipiv(i) /= i) then
        g = -g
     endif
  enddo
  write(*,*) "det(A) = ", g
end program main
