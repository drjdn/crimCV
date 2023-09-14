!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Copyright (C) 2023 Jason D. Nielsen
!!!
!!! This program is free software! you can redistribute it and/or modify
!!! it under the terms of the GNU General Public License as published by
!!! the Free Software Foundation! either version 2 of the License, or
!!! (at your option) any later version.
!!!
!!! This program is distributed in the hope that it will be useful,
!!! but WITHOUT ANY WARRANTY! without even the implied warranty of
!!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!! GNU General Public License for more details.
!!!
!!! You should have received a copy of the GNU General Public License
!!! along with this program! if not, write to the Free Software
!!! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111, USA
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module mytype
  implicit none
  integer, parameter:: dp=kind(1.0D0)
  integer, parameter:: qp=2*dp
end module mytype

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Copyright (C) 2023 Jason D. Nielsen
!!!
!!! This program is free software! you can redistribute it and/or modify
!!! it under the terms of the GNU General Public License as published by
!!! the Free Software Foundation! either version 2 of the License, or
!!! (at your option) any later version.
!!!
!!! This program is distributed in the hope that it will be useful,
!!! but WITHOUT ANY WARRANTY! without even the implied warranty of
!!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!! GNU General Public License for more details.
!!!
!!! You should have received a copy of the GNU General Public License
!!! along with this program! if not, write to the Free Software
!!! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111, USA
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module dat_mod
  use mytype
  implicit none
  integer, save:: g_nn, g_ni, g_no, g_npp, g_npl, g_ng, g_npg
  integer, parameter:: nmaxit=10000
  real(dp), allocatable, save:: g_Dat(:,:), g_y(:), g_X(:,:), g_Z(:,:), g_expd(:), g_zero(:), &
       & g_nzero(:), g_llc(:), g_miss(:), g_offt(:)
  integer, allocatable, save:: g_indi(:)
  real(dp), allocatable, save, target:: g_pr_wt(:,:), g_llike_t(:,:)
  real(dp), pointer:: g_llikei(:) => null(), g_gwt(:) => null()
  logical:: isEM
end module dat_mod

module gamma_mod
  use mytype
  implicit none

  interface gammafn
     module procedure gamma_d, gamma_dv
  end interface

  interface lgammafn
     module procedure lgamma_d, lgamma_dv
  end interface

contains      

  function gamfn(x,kf) result(gl)
  !
  !       ==================================================
  !       Purpose: Compute gamma (x) or (x)]
  !       Input:   x  --- Argument (x) ( x > 0 )
  !                KF --- Function code
  !                       KF=1 a(x); KF=0 for ln[a(x)]
  !       Output:  GL --- a(x) or ln[a(x)] 
  !       ==================================================
  !
    implicit none
    real(dp), intent(in):: x
    integer, intent(in):: kf
    real(dp):: gl
    real(dp), dimension(10),parameter:: a=(/8.333333333333333e-02_dp, &
         & -2.777777777777778e-03_dp, 7.936507936507937e-04_dp, &
         & -5.952380952380952e-04_dp, 8.417508417508418e-04_dp, &
         & -1.917526917526918e-03_dp, 6.410256410256410e-03_dp, &
         & -2.955065359477124e-02_dp, 1.796443723688307e-01_dp, &
         & -1.39243221690590_dp/)
    real(dp):: x0, x2, xp, glo
    integer:: k, n
    x0=x
    if (x == 1.0 .OR. x == 2.0) then
       gl=0.0_dp
       if (kf == 1) gl=exp(gl)
       return
    else if (x <= 7.0) then
       n=int(7-x)
       x0=x+n
    end if
    x2=1.0_dp/(x0*x0)
    xp=6.283185307179586477_dp
    glo=A(10)
    do k=9,1,-1
       glo=glo*x2+a(k)
    end do
    gl=glo/x0+0.5_dp*log(xp)+(x0-.5_dp)*log(x0)-x0
    if (x <= 7.0) then
       do k=1,n
          gl=gl-log(x0-1.0_dp)
          x0=x0-1.0_dp
       end do
    end if
    if (kf == 1) gl=exp(gl)
    return
  end function gamfn

  function lgamma_d(x) result(gmln)
    implicit none
    real(dp), intent(in):: x
    real(dp):: gmln
    gmln=gamfn(x,0)       
  end function lgamma_d

  function gamma_d(x) result(gmln)
    implicit none
    real(dp), intent(in):: x
    real(dp):: gmln
    gmln=gamfn(x,1)
  end function gamma_d

  function lgamma_dv(x) result(gmln)
    implicit none
    real(dp), intent(in):: x(:)
    real(dp):: gmln(size(x))
    integer:: nx, i
    nx=size(x)
    do i=1,nx
       gmln(i)=gamfn(x(i),0)
    end do
  end function lgamma_dv

  function gamma_dv(x) result(gmln)
    implicit none
    real(dp), intent(in):: x(:)
    real(dp):: gmln(size(x))
    integer:: nx, i
    nx=size(x)
    do i=1,nx
       gmln(i)=gamfn(x(i),1)
    end do
  end function gamma_dv
  
end module gamma_mod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Copyright (C) 2023 Jason D. Nielsen
!!!
!!! This program is free software! you can redistribute it and/or modify
!!! it under the terms of the GNU General Public License as published by
!!! the Free Software Foundation! either version 2 of the License, or
!!! (at your option) any later version.
!!!
!!! This program is distributed in the hope that it will be useful,
!!! but WITHOUT ANY WARRANTY! without even the implied warranty of
!!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!! GNU General Public License for more details.
!!!
!!! You should have received a copy of the GNU General Public License
!!! along with this program! if not, write to the Free Software
!!! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111, USA
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module matrix
  use mytype
  implicit none

  interface operator (.mm.)
     module procedure mm, vm, mv, vv
  end interface

  interface operator (.tmm.)
     module procedure tmm
  end interface

  interface operator (.mmt.)
     module procedure mmt
  end interface

  interface operator (.mmd.)
     module procedure mmdiag_l, mmdiag_r
  end interface

  interface operator (.kp.)
     module procedure kron_mm, kron_mv, kron_vm
  end interface

  interface fsymsolve
     module procedure fsymsolve_v, fsymsolve_m
  end interface

  interface rsymsolve
     module procedure rsymsolve_v, rsymsolve_m
  end interface

contains      

!  subroutine printmat(A)
!    implicit none
!    real(dp), dimension(:,:):: A
!    integer:: i, j
!    do i=1,size(A,1)
!       write(*,fmt="(a2)",advance="no") "|"
!       do j=1,size(A,2)
!          write(*,fmt="(f9.4)",advance="no") A(i,j)
!       end do
!       write(*,fmt="(a2)") " |"
!    end do
!    write(*,*)
!  end subroutine printmat

  function mm(A,B) result(out)
    implicit none
    real(dp), dimension(:,:), intent(in):: A, B
    real(dp), dimension(size(A,1),size(B,2)):: out
    call dgemm_f95(A,B,out,'N','N')
  end function mm

  function tmm(A,B) result(out)
    implicit none
    real(dp), dimension(:,:), intent(in):: A, B
    real(dp), dimension(size(A,2),size(B,2)):: out
    call dgemm_f95(A,B,out,'T','N')
  end function tmm

  function mmt(A,B) result(out)
    implicit none
    real(dp), dimension(:,:), intent(in):: A, B
    real(dp), dimension(size(A,1),size(B,1)):: out
    call dgemm_f95(A,B,out,'N','T')
  end function mmt

  function mv(A,B) result(out)
    implicit none
    real(dp), dimension(:,:), intent(in):: A
    real(dp), dimension(:), intent(in):: B
    real(dp), dimension(size(A,1)):: out
    call dgemv_f95(A,B,out,'N')
  end function mv

  function vm(A,B) result(out)
    implicit none
    real(dp), dimension(:,:), intent(in)::  B
    real(dp), dimension(:), intent(in):: A
    real(dp), dimension(size(B,2)):: out
    call dgemv_f95(B,A,out,'T')
  end function vm

  function vv(A,B) result(out)
    implicit none
    real(dp), dimension(:), intent(in)::  B
    real(dp), dimension(:), intent(in):: A
    real(dp):: out
    out=dot_product(A,B)
  end function vv

  SUBROUTINE dgemm_f95(a,b,c,t1,t2)
    ! .. Use Statements ..
    ! .. Parameters ..
    REAL (dp), PARAMETER :: zero = 0.0_dp
    REAL (dp), PARAMETER :: one = 1.0_dp
    ! .. Array Arguments ..
    REAL (dp), INTENT (IN) :: a(:,:), b(:,:)
    CHARACTER, INTENT (IN) :: t1, t2
    REAL (dp), INTENT (OUT) :: c(:,:)
    ! .. Local Scalar ..
    INTEGER :: k, m, n
    ! .. External Procedures ..
    EXTERNAL dgemm  
    m = SIZE(c,1)
    n = SIZE(c,2)
    IF (t1 == 'T'.and. t2 == 'T') then
       k = SIZE(a,1)
       CALL dgemm(t1,t2,m,n,k,one,a,k,b,n,zero,c,m)
    ELSEIF (t1 == 'T'.and. t2 == 'N') then
       k = SIZE(a,1)
       CALL dgemm(t1,t2,m,n,k,one,a,k,b,k,zero,c,m)
    ELSEIF (t1 == 'N'.and. t2 == 'T') then
       k = SIZE(a,2)
       CALL dgemm(t1,t2,m,n,k,one,a,m,b,n,zero,c,m)
    ELSE
       k = SIZE(a,2)
       CALL dgemm(t1,t2,m,n,k,one,a,m,b,k,zero,c,m)
    END IF
  END SUBROUTINE dgemm_f95
  
  SUBROUTINE dgemv_f95(a,b,c,t1)
    ! .. Use Statements ..
    ! .. Parameters ..
    REAL (dp), PARAMETER :: one = 1.0_dp
    REAL (dp), PARAMETER :: zero = 0.0_dp
    ! .. Array Arguments ..
    CHARACTER, INTENT (IN) :: t1
    REAL (dp), INTENT (IN) :: a(:,:), b(:)
    REAL (dp), INTENT (OUT) :: c(:)
    ! .. Local Scalar ..
    INTEGER :: m, n
    ! .. External Procedures ..
    EXTERNAL dgemv  
    m = SIZE(a,1)
    n = SIZE(a,2)
    CALL dgemv(t1,m,n,one,a,m,b,1,zero,c,1)
  END SUBROUTINE dgemv_f95

  function mmdiag_l(A,diag)
    implicit none
    real(dp), dimension(:,:), intent(in):: A
    real(dp), dimension(size(A,2)), intent(in):: diag
    real(dp), dimension(size(A,1),size(A,2)):: mmdiag_l
    integer:: i
    do i=1,size(A,2)
       mmdiag_l(:,i)=A(:,i)*diag(i)
    end do
  end function mmdiag_l

  function mmdiag_r(diag,A)
    implicit none
    real(dp), dimension(:,:), intent(in):: A
    real(dp), dimension(size(A,1)), intent(in):: diag
    real(dp), dimension(size(A,1),size(A,2)):: mmdiag_r
    integer:: i
    do i=1,size(A,1)
       mmdiag_r(i,:)=A(i,:)*diag(i)
    end do
  end function mmdiag_r

  function ddiag(A)
    implicit none
    real(dp), dimension(:,:), intent(in):: A
    integer:: i
    real(dp), dimension(size(A,1)):: ddiag
    do i=1,size(A,1)
       ddiag(i)=A(i,i)
    end do
  end function ddiag

  subroutine fsymsolve_v(A,b,syslv,info,ldet)
    implicit none
    real(dp), dimension(:,:), intent(in):: A
    real(dp), dimension(size(A,2)), intent(in):: b
    real(dp), dimension(size(A,2)), intent(out):: syslv
    integer, intent(out):: info
    real(dp), optional, intent(out):: ldet
    real(dp), dimension(size(A,1),size(A,2)):: D
    real(dp), dimension(size(A,2),1):: c
    integer:: n, nrhs, ipiv(size(A,1)), lwrk
    real(dp), allocatable:: wrk(:)
    real(dp):: w(1)
    n=size(A,1)
    nrhs=1
!    if(n/=size(A,2)) then
!       stop 'Matrix must be square'
!    end if
    D=A
    c=reshape(b,(/n,1/))
    lwrk=-1
    call dsysv('U',n,nrhs,D,n,ipiv,c,n,w,lwrk,info)
    lwrk=w(1)
    allocate(wrk(lwrk))
    call dsysv('U',n,nrhs,D,n,ipiv,c,n,wrk,lwrk,info)
    if (present(ldet)) then
       ldet=sum(log(ddiag(D)))
    end if
    syslv=pack(c,.TRUE.)
    deallocate(wrk)
  end subroutine fsymsolve_v

  subroutine fsymsolve_m(A,B,syslv,info,ldet)
    implicit none
    real(dp), dimension(:,:), intent(in):: A
    real(dp), dimension(:,:), intent(in):: B
    real(dp), dimension(size(A,2),size(B,2)), intent(out):: syslv
    integer, intent(out):: info
    real(dp), optional, intent(out):: ldet
    real(dp), dimension(size(A,1),size(A,2)):: D
    real(dp), dimension(size(A,2),size(B,2)):: C
    integer:: n, nrhs, ipiv(size(A,1)), lwrk
    real(dp), allocatable:: wrk(:)
    real(dp):: w(1)
    n=size(A,1)
    nrhs=size(B,2)
!    if(n/=size(A,2)) then
!       stop 'Matrix must be square'
!    end if
    D=A
    C=B
    lwrk=-1
    call dsysv('U',n,nrhs,D,n, ipiv, C, n, w, lwrk, info)
    lwrk=w(1)
    allocate(wrk(lwrk))
    call dsysv('U',n,nrhs,D,n, ipiv, C, n, wrk, lwrk, info)
    if (present(ldet)) then
       ldet=sum(log(ddiag(D)))
    end if
    syslv=C
    deallocate(wrk)
  end subroutine fsymsolve_m

  subroutine symeigen(B,Eval,Evec)
    implicit none
    real(dp), intent(in):: B(:,:)
    real(dp), intent(out):: Eval(size(B,1)), Evec(size(B,1),size(B,1))
    real(dp):: A(size(B,1),size(B,1)), abstol, vl, vu
    integer:: lwork, liwork, na, info, il, iu, ne
    real(dp), allocatable:: work(:)
    integer, allocatable:: iwork(:), isuppz(:)
    A=B
    na=size(A,1)
    abstol=0.0_dp
    lwork=50*na
    liwork=10*na
    il=1
    iu=na
    vu=huge(1.0_dp)
    vl=-vu
    allocate(work(lwork),iwork(liwork),isuppz(2*na))
    call dsyevr('V','A','U',na,A,na,vl,vu,il,iu,abstol,ne,Eval,Evec,&
         & na,isuppz,work,lwork,iwork,liwork,info)
    deallocate(work,iwork,isuppz)
  end subroutine symeigen

  subroutine rsymsolve_v(A,b,syslv)
    implicit none
    real(dp), dimension(:,:), intent(in):: A
    real(dp), dimension(size(A,2)), intent(in):: b
    real(dp), dimension(size(A,2)), intent(out):: syslv
    real(dp), dimension(size(A,1),size(A,2)):: D
    real(dp), dimension(size(A,2)):: c
!    if(size(A,1)/=size(A,2)) then
!       stop 'Matrix must be square'
!    end if
    call symeigen(A,c,D)
    syslv=b.mm.D
    where (abs(c) < 1.0e-8_dp)
       c=0.0_dp
    elsewhere
       c=1/c
    end where
    syslv=c*syslv
    syslv=D.mm.syslv
  end subroutine rsymsolve_v

  subroutine rsymsolve_m(A,b,syslv)
    implicit none
    real(dp), dimension(:,:), intent(in):: A
    real(dp), dimension(:,:), intent(in):: B
    real(dp), dimension(size(A,2),size(B,2)), intent(out):: syslv
    real(dp), dimension(size(A,1),size(A,2)):: D
    real(dp), dimension(size(A,2)):: c
!    if(size(A,1)/=size(A,2)) then
!       stop 'Matrix must be square'
!    end if
    call symeigen(A,c,D)
    syslv=D.tmm.b
    where (abs(c) < 1.0e-8_dp)
       c=0.0_dp
    elsewhere
       c=1/c
    end where
    syslv=c.mmd.syslv
    syslv=D.mm.syslv
  end subroutine rsymsolve_m

  subroutine kronr(a,ia,ma,na,b,ib,mb,nb,pk,ik)
      implicit none
      real(dp) a(*),b(*),pk(*)
      integer ia,ma,na,ib,mb,nb,ik,ka,kb,kk,ka1,kk1,ja,jb,i
      ka1=1-ia
      kk1=-nb
      do 30 ja=1,na
         kb=1
         ka1=ka1+ia
         kk1=kk1+nb
         do 20 jb=1,nb
            ka=ka1
            kk=1+(jb-1+kk1)*ik
            do 10 i=1,ma
               call dcopy(mb,b(kb),1,pk(kk),1)
               call dscal(mb,a(ka),pk(kk),1)
               kk=kk+mb
               ka=ka+1
 10         continue
            kb=kb+ib
 20      continue
 30   continue
      return
   end subroutine kronr

  function kron_mm(A, B) result(out)
    implicit none
    real(dp), dimension(:,:), intent(in):: A, B
    real(dp), dimension(size(A,1)*size(B,1),size(A,2)*size(B,2)):: out
    integer, dimension(2):: an, bn
    an=shape(A)
    bn=shape(B)
    call kronr(A,an(1),an(1),an(2),B,bn(1),bn(1),bn(2),out,an(1)*bn(1))
  end function kron_mm

  function kron_vm(A, B) result(out)
    implicit none
    real(dp), dimension(:), intent(in):: A
    real(dp), dimension(:,:), intent(in):: B
    real(dp), dimension(size(A)*size(B,1),size(B,2)):: out
    integer:: an, bn(2)
    an=size(A)
    bn=shape(B)
    call kronr(A,an,an,1,B,bn(1),bn(1),bn(2),out,an*bn(1))
  end function kron_vm

  function kron_mv(A, B) result(out)
    implicit none
    real(dp), dimension(:,:), intent(in):: A
    real(dp), dimension(:), intent(in):: B
    real(dp), dimension(size(A,1)*size(B,1),size(A,2)):: out
    integer:: an(2), bn
    an=shape(A)
    bn=size(B)
    call kronr(A,an(1),an(1),an(2),B,bn,bn,1,out,an(1)*bn)
  end function kron_mv

  subroutine flsqr(X,y,coef,info)
    implicit none
    real(dp), dimension(:,:), intent(in):: X
    real(dp), dimension(:),intent(in):: y
    integer, intent(out):: info
    real(dp), dimension(size(X,2)), intent(out):: coef
    real(dp), allocatable:: work(:)
    real(dp), dimension(size(X,1),1):: ytmp
    real(dp), dimension(size(X,1),size(X,2)):: Xtmp
    integer:: m, n, lwrk, jpvt(size(X,2)), rank
    real(dp):: w(1)
    m=size(X,1)
    n=size(X,2)
    Xtmp=X
    ytmp=reshape(y,(/m,1/))
    lwrk=-1
    call dgelsy(m,n,1,Xtmp,m,ytmp,m,jpvt,1.0e-8_dp,rank,w,lwrk,info)
    lwrk=w(1)
    allocate(work(lwrk))
    call dgelsy(m,n,1,Xtmp,m,ytmp,m,jpvt,1.0e-8_dp,rank,work,lwrk,info)
    coef=ytmp(1:n,1)
    deallocate(work)
  end subroutine flsqr

end module matrix

! This code was modified by Jason D. Nielsen 
module merge_sort_mod
  use mytype
  private
  public:: mrgrnk

contains

! This function is from Michel Olagnon's ORDERPACK 2.0
Subroutine mrgrnk (XDONT, IRNGT)
! __________________________________________________________
!   MRGRNK = Merge-sort ranking of an array
!   For performance reasons, the first 2 passes are taken
!   out of the standard loop, and use dedicated coding.
! __________________________________________________________
! __________________________________________________________
      Real (kind=dp), Dimension (:), Intent (In) :: XDONT
      Integer, Dimension (:), Intent (Out) :: IRNGT
! __________________________________________________________
      Real (kind=dp) :: XVALA, XVALB
!
      Integer, Dimension (SIZE(IRNGT)) :: JWRKT
      Integer :: LMTNA, LMTNC, IRNG1, IRNG2
      Integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
!
      NVAL = Min (SIZE(XDONT), SIZE(IRNGT))
      Select Case (NVAL)
      Case (:0)
         Return
      Case (1)
         IRNGT (1) = 1
         Return
      Case Default
         Continue
      End Select
!
!  Fill-in the index array, creating ordered couples
!
      Do IIND = 2, NVAL, 2
         If (XDONT(IIND-1) <= XDONT(IIND)) Then
            IRNGT (IIND-1) = IIND - 1
            IRNGT (IIND) = IIND
         Else
            IRNGT (IIND-1) = IIND
            IRNGT (IIND) = IIND - 1
         End If
      End Do
      If (Modulo(NVAL, 2) /= 0) Then
         IRNGT (NVAL) = NVAL
      End If
!
!  We will now have ordered subsets A - B - A - B - ...
!  and merge A and B couples into     C   -   C   - ...
!
      LMTNA = 2
      LMTNC = 4
!
!  First iteration. The length of the ordered subsets goes from 2 to 4
!
      Do
         If (NVAL <= 2) Exit
!
!   Loop on merges of A and B into C
!
         Do IWRKD = 0, NVAL - 1, 4
            If ((IWRKD+4) > NVAL) Then
               If ((IWRKD+2) >= NVAL) Exit
!
!   1 2 3
!
               If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
!
!   1 3 2
!
               If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                  IRNG2 = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNG2
!
!   3 1 2
!
               Else
                  IRNG1 = IRNGT (IWRKD+1)
                  IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNG1
               End If
               Exit
            End If
!
!   1 2 3 4
!
            If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
!
!   1 3 x x
!
            If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
               If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   1 3 2 4
                  IRNGT (IWRKD+3) = IRNG2
               Else
!   1 3 4 2
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+4) = IRNG2
               End If
!
!   3 x x x
!
            Else
               IRNG1 = IRNGT (IWRKD+1)
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
               If (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) Then
                  IRNGT (IWRKD+2) = IRNG1
                  If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   3 1 2 4
                     IRNGT (IWRKD+3) = IRNG2
                  Else
!   3 1 4 2
                     IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                     IRNGT (IWRKD+4) = IRNG2
                  End If
               Else
!   3 4 1 2
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+3) = IRNG1
                  IRNGT (IWRKD+4) = IRNG2
               End If
            End If
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 4
         Exit
      End Do
!
!  Iteration loop. Each time, the length of the ordered subsets
!  is doubled.
!
      Do
         If (LMTNA >= NVAL) Exit
         IWRKF = 0
         LMTNC = 2 * LMTNC
!
!   Loop on merges of A and B into C
!
         Do
            IWRK = IWRKF
            IWRKD = IWRKF + 1
            JINDA = IWRKF + LMTNA
            IWRKF = IWRKF + LMTNC
            If (IWRKF >= NVAL) Then
               If (JINDA >= NVAL) Exit
               IWRKF = NVAL
            End If
            IINDA = 1
            IINDB = JINDA + 1
!
!   Shortcut for the case when the max of A is smaller
!   than the min of B. This line may be activated when the
!   initial set is already close to sorted.
!
!          IF (XDONT(IRNGT(JINDA)) <= XDONT(IRNGT(IINDB))) CYCLE
!
!  One steps in the C subset, that we build in the final rank array
!
!  Make a copy of the rank array for the merge iteration
!
            JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
!
            XVALA = XDONT (JWRKT(IINDA))
            XVALB = XDONT (IRNGT(IINDB))
!
            Do
               IWRK = IWRK + 1
!
!  We still have unprocessed values in both A and B
!
               If (XVALA > XVALB) Then
                  IRNGT (IWRK) = IRNGT (IINDB)
                  IINDB = IINDB + 1
                  If (IINDB > IWRKF) Then
!  Only A still with unprocessed values
                     IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                     Exit
                  End If
                  XVALB = XDONT (IRNGT(IINDB))
               Else
                  IRNGT (IWRK) = JWRKT (IINDA)
                  IINDA = IINDA + 1
                  If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                  XVALA = XDONT (JWRKT(IINDA))
               End If
!
            End Do
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 2 * LMTNA
      End Do
!
      Return
!
End Subroutine mrgrnk

end module merge_sort_mod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Copyright (C) 2023 Jason D. Nielsen
!!!
!!! This program is free software! you can redistribute it and/or modify
!!! it under the terms of the GNU General Public License as published by
!!! the Free Software Foundation! either version 2 of the License, or
!!! (at your option) any later version.
!!!
!!! This program is distributed in the hope that it will be useful,
!!! but WITHOUT ANY WARRANTY! without even the implied warranty of
!!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!! GNU General Public License for more details.
!!!
!!! You should have received a copy of the GNU General Public License
!!! along with this program! if not, write to the Free Software
!!! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111, USA
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module rrand
  use mytype
  implicit none

  public
  private:: ranu_s, ranu_v
  
  interface r_random_number
     module procedure ranu_s, ranu_v, ranu_m
  end interface r_random_number

contains

  subroutine r_set_random
    implicit none
    call sRNG()
  end subroutine r_set_random

  subroutine r_close_random
    implicit none
    call eRNG()
  end subroutine r_close_random
  
  subroutine ranu_s(rv)
    implicit none
    real(dp), intent(inout):: rv
    real(dp):: runi
    rv = runi()
  end subroutine ranu_s

  subroutine ranu_v(rv)
    implicit none
    real(dp), dimension(:), intent(inout):: rv
    real(dp):: runi
    integer:: i
    do i=1,size(rv)
       rv(i) = runi()
    end do
  end subroutine ranu_v

  subroutine ranu_m(rv)
    implicit none
    real(dp), dimension(:,:), intent(inout):: rv
    real(dp):: runi
    integer:: i, j
    do j=1,size(rv,2)
       do i=1,size(rv,1)
          rv(i,j) = runi()
       end do
    end do
  end subroutine ranu_m
  
end module rrand

! Modified by Jason D. Nielsen
module DE_mod
  use mytype
  use rrand
  implicit none
  private
  public:: DE_optim

contains

  subroutine DE_optim(obj, Dim_XC, XCmin, XCmax, VTR, NP, itermax, F_XC, &
           CR_XC, strategy, refresh, pop_XC, val, nfeval, F_CR, method)
    !.......................................................................
    !    
    ! Differential Evolution for Optimal Control Problems
    !
    !.......................................................................
    !  This Fortran 90 program translates from the original MATLAB 
    !  version of differential evolution (DE). This FORTRAN 90 code 
    !  has been tested on Compaq Visual Fortran v6.1. 
    !  Any users new to the DE are encouraged to read the article of Storn and Price. 
    !
    !  Refences:
    !  Storn, R., and Price, K.V., (1996). Minimizing the real function of the 
    !    ICEC'96 contest by differential evolution. IEEE conf. on Evolutionary 
    !    Comutation, 842-844.
    !
    !  This Fortran 90 program written by Dr. Feng-Sheng Wang 
    !  Department of Chemical Engineering, National Chung Cheng University, 
    !  Chia-Yi 621, Taiwan, e-mail: chmfsw@ccunix.ccu.edu.tw
    !.........................................................................
    !                obj : The user provided file for evlauting the objective function.
    !                      subroutine obj(xc,fitness)
    !                      where "xc" is the real decision parameter vector.(input)
    !                            "fitness" is the fitness value.(output)
    !             Dim_XC : Dimension of the real decision parameters.
    !      XCmin(Dim_XC) : The lower bound of the real decision parameters.
    !      XCmax(Dim_XC) : The upper bound of the real decision parameters.
    !                VTR : The expected fitness value to reach.
    !                 NP : Population size.
    !            itermax : The maximum number of iteration.
    !               F_XC : Mutation scaling factor for real decision parameters.
    !              CR_XC : Crossover factor for real decision parameters.
    !           strategy : The strategy of the mutation operations is used in HDE.
    !            refresh : The intermediate output will be produced after "refresh"
    !                      iterations. No intermediate output will be produced if
    !                      "refresh < 1".
    ! bestmen_XC(Dim_XC) : The best real decision parameters.
    !              bestval : The best objective function.
    !             nfeval : The number of function call.
    !         method(1) = 0, Fixed mutation scaling factors (F_XC)
    !                   = 1, Random mutation scaling factors F_XC=[0, 1]
    !                   = 2, Random mutation scaling factors F_XC=[-1, 1] 
    !         method(2) = 1, Random combined factor (F_CR) used for strategy = 6
    !                        in the mutation operation 
    !                   = other, fixed combined factor provided by the user 
    
    use mytype
    implicit none
    integer, intent(in) :: NP, Dim_XC, itermax, strategy, refresh
    real(dp), intent(in) :: VTR, CR_XC
    real(dp) :: F_XC, F_CR
    real(dp), dimension(Dim_XC), intent(in) :: XCmin, XCmax
    real(dp), dimension(NP,Dim_XC), intent(out) :: pop_XC
    real(dp), dimension(NP), intent(out) :: val
    integer, intent(out) :: nfeval     
    real(dp), dimension(NP,Dim_XC) :: bm_XC, mui_XC, mpo_XC,   &
         popold_XC, rand_XC, ui_XC
    integer :: i, ibest, iter 
    integer, dimension(NP) :: rot, a1, a2, a3, a4, a5, rt
    integer, dimension(4) :: ind
    real(dp) :: tempval, bestval
    real(dp), dimension(Dim_XC) :: bestmemit_XC, bestmem_XC
    real(dp), dimension(Dim_XC) :: rand_C1
    integer, dimension(2), intent(in) :: method
    intrinsic max, min, mod, abs, any, all, maxloc
    interface
       subroutine obj(n,x,objval)
         use mytype
         implicit none
         integer, intent(in):: n
         real(dp), intent(in):: x(n)
         real(dp), intent(out):: objval
       end subroutine obj
    end interface
    !!-----Initialize a population --------------------------------------------!!
    
    pop_XC=0.0_dp
    do i=1,NP
       call r_random_number(rand_C1)
       pop_XC(i,:)=XCmin+rand_C1*(XCmax-XCmin)
    end do
    
    !!--------------------------------------------------------------------------!!
    
    !!------Evaluate fitness functions and find the best member-----------------!!
    val=0.0_dp
    nfeval=0
    ibest=1
    call obj(Dim_XC, pop_XC(1,:), val(1))
    bestval=val(1)
    nfeval=nfeval+1
    do i=2,NP
       call obj(Dim_XC, pop_XC(i,:), val(i))
       nfeval=nfeval+1
       if (val(i) < bestval) then
          ibest=i
          bestval=val(i)
       end if
    end do
    bestmemit_XC=pop_XC(ibest,:)
    bestmem_XC=bestmemit_XC
    !!--------------------------------------------------------------------------!!
    
    bm_XC=0.0_dp
    rot=(/(i,i=0,NP-1)/)
    iter=1 
    !!------Perform evolutionary computation------------------------------------!! 
    
    do while (iter <= itermax)
       popold_XC=pop_XC
       
       !!------Mutation operation--------------------------------------------------!!
       ind=randperm(4)
       a1=randperm(NP)
       rt=mod(rot+ind(1),NP)
       a2=a1(rt+1)
       rt=mod(rot+ind(2),NP)
       a3=a2(rt+1)
       rt=mod(rot+ind(3),NP)
       a4=a3(rt+1)
       rt=mod(rot+ind(4),NP)
       a5=a4(rt+1)
       bm_XC=spread(bestmemit_XC, DIM=1, NCOPIES=NP)
       
       !----- Generating a random sacling factor--------------------------------!
       select case (method(1))
       case (1)
          call r_random_number(F_XC)
       case(2)
          call r_random_number(F_XC)
          F_XC=2.0_dp*F_XC-1.0_dp
       end select
       
       !---- select a mutation strategy-----------------------------------------!
       select case (strategy)
       case (1)
          ui_XC=bm_XC+F_XC*(popold_XC(a1,:)-popold_XC(a2,:))
          
       case default
          ui_XC=popold_XC(a3,:)+F_XC*(popold_XC(a1,:)-popold_XC(a2,:))
          
       case (3)
          ui_XC=popold_XC+F_XC*(bm_XC-popold_XC+popold_XC(a1,:)-popold_XC(a2,:))
          
       case (4)
          ui_XC=bm_XC+F_XC*(popold_XC(a1,:)-popold_XC(a2,:)+popold_XC(a3,:)-popold_XC(a4,:))
          
       case (5)
          ui_XC=popold_XC(a5,:)+F_XC*(popold_XC(a1,:)-popold_XC(a2,:)+popold_XC(a3,:)-popold_XC(a4,:))
       case (6) ! A linear crossover combination of bm_XC and popold_XC
          if (method(2) == 1) call r_random_number(F_CR) 
          ui_XC=popold_XC+F_CR*(bm_XC-popold_XC)+F_XC*(popold_XC(a1,:)-popold_XC(a2,:))          
       end select
       !!--------------------------------------------------------------------------!!
       !!------Crossover operation-------------------------------------------------!!
       call r_random_number(rand_XC)
       mui_XC=0.0_dp
       mpo_XC=0.0_dp
       where (rand_XC < CR_XC)
          mui_XC=1.0_dp
          !           mpo_XC=0.0_dp
       elsewhere
          !           mui_XC=0.0_dp
          mpo_XC=1.0_dp
       end where
       
       ui_XC=popold_XC*mpo_XC+ui_XC*mui_XC
       !!--------------------------------------------------------------------------!!
       !!------Evaluate fitness functions and find the best member-----------------!!
       do i=1,NP
          !!------Confine each of feasible individuals in the lower-upper bound-------!!
          ui_XC(i,:)=max(min(ui_XC(i,:),XCmax),XCmin)
          call obj(Dim_XC, ui_XC(i,:), tempval)
          nfeval=nfeval+1
          if (tempval < val(i)) then
             pop_XC(i,:)=ui_XC(i,:)
             val(i)=tempval
             if (tempval < bestval) then
                bestval=tempval
                bestmem_XC=ui_XC(i,:)
             end if
          end if
       end do
       bestmemit_XC=bestmem_XC
!       if( (refresh > 0) .and. (mod(iter,refresh)==0)) then
!          write(unit=*, FMT=203) iter
!          do i=1,Dim_XC
!             write(*,FMT=202) i,bestmem_XC(i)
!          end do
!          write(unit=*, FMT=201) bestval
!          write(unit=*, FMT=202) maxval(val)-minval(val)
!          write(unit=*, FMT=204) sddev(NP,val)
!       end if
       iter=iter+1
       if ( bestval <= VTR .and. refresh > 0) then
!          write(unit=*, FMT=*) 'The best fitness is smaller than VTR' 
          exit
       endif
    end do
    !!------end the evolutionary computation------------------------------!!
!201 format(2x, 'bestval =', ES14.7)
!202 format(2x, 'range =', ES14.7)
!204 format(2x, 'sd =', ES14.7, /)
!202 format(5x, 'bestmem_XC(', I3, ') =', ES12.5)
!203 format(2x, 'No. of iteration  =', I8)
  end subroutine DE_optim

  function sddev(np,x) result(out)
    use mytype
    implicit none
    integer, intent(in):: np
    real(dp), intent(in):: x(np)
    real(dp):: out, mu
    mu=sum(x)/real(np,dp)
    out=sum((x-mu)**2)/real(np-1,dp)    
  end function sddev
  
  function randperm(num)
    use mytype
    implicit none
    integer, intent(in) :: num
    integer :: number, i, j, k
    integer, dimension(num) :: randperm
    real(dp), dimension(num) :: rand2
    call r_random_number(rand2)
    do i=1,num
       number=1
       do j=1,num
          if (rand2(i) > rand2(j)) then
             number=number+1
          end if
       end do
       do k=1,i-1
          if (rand2(i) <= rand2(k) .and. rand2(i) >= rand2(k)) then
             number=number+1
          end if
       end do
       randperm(i)=number
    end do
    return
  end function randperm

end module DE_mod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Copyright (C) 2023 Jason D. Nielsen
!!!
!!! This program is free software! you can redistribute it and/or modify
!!! it under the terms of the GNU General Public License as published by
!!! the Free Software Foundation! either version 2 of the License, or
!!! (at your option) any later version.
!!!
!!! This program is distributed in the hope that it will be useful,
!!! but WITHOUT ANY WARRANTY! without even the implied warranty of
!!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!! GNU General Public License for more details.
!!!
!!! You should have received a copy of the GNU General Public License
!!! along with this program! if not, write to the Free Software
!!! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111, USA
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module Newton_mod
  implicit none
  private
  public:: Newton

contains 

  subroutine Newton(fgrdH,x,llike,grd,Hess,np,ierr) 
    use mytype
    use matrix, only: fsymsolve, operator(.mm.)
    implicit none
    interface
       subroutine fgrdH(n,x,f,grd,H)
         use mytype
         implicit none
         integer, intent(in):: n
         real(dp), intent(in):: x(n)
         real(dp), intent(out):: f, grd(n), H(n,n)
       end subroutine fgrdH
    end interface
    integer, intent(in):: np
    real(dp), intent(inout):: x(np)
    real(dp), intent(out):: llike, grd(np), Hess(np,np)
    integer, intent(inout):: ierr
    integer:: i, j
    real(dp), parameter:: toler1=1.0D-6
    real(dp):: s(np), xval, x_old(np), xtmp(np), fvl, stp
    integer, parameter:: maxit=150
    ierr=0
    x_old=x
!    write(*,*)
    do i=1,maxit
       call fgrdH(np,x,llike,grd,Hess)
!       write(*,*) i, llike, maxval(abs(grd))
       call fsymsolve(Hess,grd,s,ierr)
       if (ierr /= 0) then
          return
       end if
       xval=maxval(abs(grd))
       stp=1.0_dp
       do j=1,5
          xtmp=x+stp*s
          if (i > 20) then
             x=xtmp
             exit
          end if
          call fgrdH(np,xtmp,fvl,grd,Hess)
          if (fvl > llike) then
             x=xtmp
             llike=fvl
             exit
          end if
          stp=0.5_dp*stp
          if (j==5) then
             x=xtmp
          end if
       end do
       if (xval<toler1) then
          call fgrdH(np,x,llike,grd,Hess)          
          return
       end if
       x_old=x
       if(i==maxit) then
          ierr=999
       end if
    end do
  end subroutine Newton

end module Newton_mod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Copyright (C) 2023 Jason D. Nielsen
!!!
!!! This program is free software! you can redistribute it and/or modify
!!! it under the terms of the GNU General Public License as published by
!!! the Free Software Foundation! either version 2 of the License, or
!!! (at your option) any later version.
!!!
!!! This program is distributed in the hope that it will be useful,
!!! but WITHOUT ANY WARRANTY! without even the implied warranty of
!!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!! GNU General Public License for more details.
!!!
!!! You should have received a copy of the GNU General Public License
!!! along with this program! if not, write to the Free Software
!!! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111, USA
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module dmZIP_shared_mod
  implicit none
contains

  subroutine zipt_llike(beta,gamma,llike)
    use dat_mod
    use matrix, only: operator(.mm.)
    implicit none
    real(dp), intent(in):: beta(g_npp), gamma(g_npl)
    real(dp), intent(out):: llike(g_nn)
    real(dp):: t1(g_nn),t2(g_nn),t3(g_nn),t4(g_nn),t5(g_nn),t7(g_nn),l1(g_nn),l2(g_nn),l3(g_nn)
    real(dp):: xb(g_nn), zg(g_nn)
    xb=g_X.mm.beta
    zg=g_Z.mm.gamma
    t1 = exp(zg)
    t2 = g_offt + xb
    t3 = exp(t2)
    t4 = exp(-t3)
    t5 = t1 + t4
    t7 = 1.0_dp + t1
    l1 = log(t5)
    l2 = g_y * t2 - t3
    l3 = log(t7)
    llike=g_zero*l1+g_nzero*l2-l3-g_llc
    llike=llike*g_miss
  end subroutine zipt_llike

  subroutine zipt_deriv(beta,gamma,ll1,ll2,ll3,l1b,l1g,l2b,l3g,l1bb,l1gb,l1gg,l2bb,l3gg)
    use dat_mod
    use matrix, only: operator(.mm.)
    implicit none
    real(dp), intent(in):: beta(g_npp), gamma(g_npl)
    real(dp), intent(out):: ll1(g_nn),ll2(g_nn),ll3(g_nn),l1b(g_nn),l1g(g_nn),l2b(g_nn),l3g(g_nn),l1bb(g_nn), &
         l1gb(g_nn),l2bb(g_nn),l1gg(g_nn),l3gg(g_nn)
    real(dp):: t1(g_nn),t2(g_nn),t3(g_nn),t4(g_nn),t5(g_nn),t7(g_nn),t8(g_nn),t9(g_nn),t11(g_nn),t14(g_nn), &
         & t16(g_nn),t17(g_nn),t21(g_nn),t23(g_nn)
    real(dp)::  xb(g_nn), zg(g_nn)
    xb=g_X.mm.beta
    zg=g_Z.mm.gamma
    t1 = exp(zg)
    t2 = g_offt + xb
    t3 = exp(t2)
    t4 = exp(-t3)
    t5 = t1 + t4
    t7 = 1.0_dp + t1
    t8 = t3 * t4
    t9 = 1.0_dp / t5
    t11 = t3 ** 2
    t14 = t4 ** 2
    t16 = t5 ** 2
    t17 = 1.0_dp / t16
    t21 = t1 ** 2
    t23 = t7 ** 2
    ll1 = log(t5)
    ll2 = g_y * t2 - t3
    ll3 = log(t7)
    l1b = -t8 * t9
    l1g = t1 * t9
    l2b = g_y - t3
    l3g = t1 / t7
    l1bb = l1b + t11 * t4 * t9 - t11 * t14 * t17
    l1gb = t1 * t17 * t8
    l2bb = -t3
    l1gg = l1g - t21 * t17
    l3gg = l3g - t21 / t23
  end subroutine zipt_deriv

  subroutine update_gllike(beta,gamma) 
    use dat_mod
    implicit none
    real(dp), intent(in):: beta(g_npp,g_ng), gamma(g_npl,g_ng)
    integer:: i
    do i=1,g_ng
       g_llikei => g_llike_t(:,i)
       call zipt_llike(beta(:,i),gamma(:,i),g_llikei)
    end do
  end subroutine update_gllike

  subroutine llike_update(prob,beta,gamma,llike) 
    use dat_mod
    implicit none
    real(dp), intent(inout):: prob(g_ng)
    real(dp), intent(in):: beta(g_npp,g_ng), gamma(g_npl,g_ng)
    real(dp), intent(out):: llike
    call update_gllike(beta,gamma) 
    call e_step(prob,llike)    
  end subroutine llike_update

  subroutine e_step(prob,llike) 
    use dat_mod
    use matrix, only: operator(.kp.)
    implicit none
    real(dp), intent(inout):: prob(g_ng)
    real(dp), intent(out):: llike
    real(dp):: dplike(g_ni,g_ng), rscl(g_ni), dlike(g_ni)
    integer:: i
    do i=1,g_ni
       dplike(i,:)=sum(g_llike_t((g_no*(i-1)+1):(g_no*i),:),1)
    end do
    do i=1,g_ng
       dplike(:,i)=log(prob(i))+dplike(:,i)
    end do
    rscl=maxval(dplike,2)
    do i=1,g_ng
       dplike(:,i)=exp(dplike(:,i)-rscl)
    end do
    dlike=sum(dplike,2)
    do i=1,g_ng
       dplike(:,i)=dplike(:,i)/dlike
    end do
    if (isEM) then
       prob=sum(dplike,1)/real(g_ni,dp)
    end if
    llike=sum(log(dlike)+rscl)
    g_pr_wt=dplike.kp.g_expd(1:g_no)
  end subroutine e_step

end module dmZIP_shared_mod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Copyright (C) 2023 Jason D. Nielsen
!!!
!!! This program is free software! you can redistribute it and/or modify
!!! it under the terms of the GNU General Public License as published by
!!! the Free Software Foundation! either version 2 of the License, or
!!! (at your option) any later version.
!!!
!!! This program is distributed in the hope that it will be useful,
!!! but WITHOUT ANY WARRANTY! without even the implied warranty of
!!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!! GNU General Public License for more details.
!!!
!!! You should have received a copy of the GNU General Public License
!!! along with this program! if not, write to the Free Software
!!! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111, USA
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module dmZIP_EM_mod
  implicit none
  private
  public:: em_update

contains

  subroutine gdHss(beta,gamma,grd,H)
    use dat_mod
    use dmZIP_shared_mod
    use matrix, only: operator(.mm.), operator(.tmm.)
    implicit none
    real(dp), intent(inout):: beta(g_npp), gamma(g_npl)
    real(dp), intent(out):: H(g_npp+g_npl,g_npp+g_npl), grd(g_npp+g_npl)
    real(dp):: t1(g_nn),t2(g_nn)
    real(dp):: ll1(g_nn),ll2(g_nn),ll3(g_nn),l1b(g_nn),l1g(g_nn),l2b(g_nn),l3g(g_nn),l1bb(g_nn), &
         & l1gb(g_nn),l1gg(g_nn),l2bb(g_nn),l3gg(g_nn)
    real(dp):: XX(g_nn,g_npp), ZZ(g_nn,g_npl)
    integer:: i 
    call zipt_deriv(beta,gamma,ll1,ll2,ll3,l1b,l1g,l2b,l3g,l1bb,l1gb,l1gg,l2bb,l3gg)
    t1=(g_zero*l1b+g_nzero*l2b)*g_gwt*g_miss
    grd(1:g_npp)=t1.mm.g_X
    t1=(g_zero*l1g-l3g)*g_gwt*g_miss
    grd((g_npp+1):(g_npp+g_npl))=t1.mm.g_Z
    t2=(g_zero*l1bb+g_nzero*l2bb)*g_gwt*g_miss
    XX=g_X
    do i=1,g_npp
       XX(:,i)=g_X(:,i)*t2
    end do
    H(1:g_npp,1:g_npp)=XX.tmm.g_X
    t2=(g_zero*l1gg-l3gg)*g_gwt*g_miss
    ZZ=g_Z
    do i=1,g_npl
       ZZ(:,i)=g_Z(:,i)*t2
    end do
    H((g_npp+1):(g_npp+g_npl),(g_npp+1):(g_npp+g_npl))=ZZ.tmm.g_Z
    t2=(g_zero*l1gb)*g_gwt*g_miss
    XX=g_X
    do i=1,g_npp
       XX(:,i)=g_X(:,i)*t2
    end do
    H(1:g_npp,(g_npp+1):(g_npp+g_npl))=XX.tmm.g_Z
    H((g_npp+1):(g_npp+g_npl),1:g_npp)=transpose(H(1:g_npp,(g_npp+1):(g_npp+g_npl)))
    H=-H
  end subroutine gdHss

  subroutine ziptfit_mle(beta,gamma,nis,grd,H,ierr) 
    use dat_mod
    use matrix, only: fsymsolve!, rsymsolve
    implicit none
    integer, intent(in):: nis
    real(dp), intent(inout):: beta(g_npp), gamma(g_npl)
    real(dp), intent(out):: grd(g_npp+g_npl), H(g_npp+g_npl,g_npp+g_npl)
    integer, intent(out):: ierr
    real(dp):: gvalo, gval, param(g_npp+g_npl), ss(g_npp+g_npl)
    integer:: i
    ierr=0
    param=(/beta,gamma/)
    do i=1,nis
       call gdHss(beta,gamma,grd,H)
       call fsymsolve(H,grd,ss,ierr)
       if (ierr /= 0) then
!          write(*,*) "ziptfit_mle: Error in fsymsolve!"
          return
       end if
!       call rsymsolve(H,grd,ss)
       param=param+0.1_dp*ss
       beta=param(1:g_npp)
       gamma=param((g_npp+1):(g_npp+g_npl))
    end do
    gvalo=maxval(abs(grd))
    do i=1,nmaxit
       call gdHss(beta,gamma,grd,H)
       gval=maxval(abs(grd))
       call fsymsolve(H,grd,ss,ierr)
       if (ierr /= 0) then
!          write(*,*) "ziptfit_mle: Error in fsymsolve!"
          return
       end if
!       call rsymsolve(H,grd,ss)
       if (gval-gvalo > 0) then
          param=param+0.5_dp*ss
       else
          param=param+ss
       end if
       if (gval < 1.0e-6_dp) then
          return
       end if
       beta=param(1:g_npp)
       gamma=param((g_npp+1):(g_npp+g_npl))
       gvalo=gval
       if (i==nmaxit) then
!          write(*,*) "ziptfit_mle: maximun number of interations reached!"
          ierr=-99
          return
       end if
    end do
  end subroutine ziptfit_mle
  
  subroutine ziptfit(nis,beta,gamma,ierr) 
    use dat_mod
    implicit none
    integer, intent(in):: nis
    real(dp), intent(inout):: beta(g_npp), gamma(g_npl)
    integer, intent(out):: ierr
    real(dp):: grd(g_npp+g_npl), H(g_npp+g_npl,g_npp+g_npl)
    call ziptfit_mle(beta,gamma,nis,grd,H,ierr) 
    if (ierr /= 0) then
!       write(*,*) "ziptfit: routine ziptfit_mle failed!"
       return
    end if
  end subroutine ziptfit

  subroutine m_step(nis,beta,gamma,ierr) 
    use dat_mod
    implicit none
    integer, intent(in):: nis
    real(dp), intent(inout):: beta(g_npp,g_ng), gamma(g_npl,g_ng)
    integer, intent(out):: ierr
    integer:: i
    do i=1,g_ng
       g_gwt => g_pr_wt(:,i)
       call ziptfit(nis,beta(:,i),gamma(:,i),ierr)
       if (ierr /= 0) then
!          write(*,*) "m_step: Error in zipfit!"
          return
       end if
    end do
  end subroutine m_step
  
  subroutine em_update(prob,beta,gamma,llike,ierr) 
    use dat_mod
    use dmZIP_shared_mod
    implicit none
    real(dp), intent(inout):: prob(g_ng), beta(g_npp,g_ng), gamma(g_npl,g_ng)
    real(dp), intent(out):: llike
    integer, intent(out):: ierr
    real(dp):: llikeo
    integer:: i 
    llikeo=-1.0e10_dp
    call llike_update(prob,beta,gamma,llike)
    do i=1,nmaxit
       call m_step(5,beta,gamma,ierr)
       if (ierr /= 0) then
!          write(*,*) "em_update: Error in m_step!"
          return
       end if
       call update_gllike(beta,gamma) 
       call e_step(prob,llike)
       if (minval(prob) < 2.0_dp/real(g_ni,dp)) then
!          write(*,*) "em_update: prob -> 0!"
          ierr=-99
          return
       end if
       if (abs(llike-llikeo) < 1.0e-6_dp) then
          exit
       end if
       llikeo=llike
    end do
  end subroutine em_update

end module dmZIP_EM_mod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Copyright (C) 2023 Jason D. Nielsen
!!!
!!! This program is free software! you can redistribute it and/or modify
!!! it under the terms of the GNU General Public License as published by
!!! the Free Software Foundation! either version 2 of the License, or
!!! (at your option) any later version.
!!!
!!! This program is distributed in the hope that it will be useful,
!!! but WITHOUT ANY WARRANTY! without even the implied warranty of
!!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!! GNU General Public License for more details.
!!!
!!! You should have received a copy of the GNU General Public License
!!! along with this program! if not, write to the Free Software
!!! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111, USA
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module dmZIP_mod
  implicit none
  private
  public:: pfun, llikefn
contains

  subroutine gdHss(beta,gamma,grd,tgrd,H)
    use dat_mod
    use dmZIP_shared_mod
    use matrix, only: operator(.tmm.), operator(.mmt.)
    implicit none
    real(dp), intent(in):: beta(g_npp), gamma(g_npl)
    real(dp), intent(out):: H(g_npp+g_npl,g_npp+g_npl), grd(g_npp+g_npl), tgrd(g_npp+g_npl,g_ni)
    real(dp):: ll1(g_nn),ll2(g_nn),ll3(g_nn),l1b(g_nn),l1g(g_nn),l2b(g_nn),l3g(g_nn), &
         & l1bb(g_nn),l1gb(g_nn),l1gg(g_nn),l2bb(g_nn),l3gg(g_nn)
    real(dp):: t1(g_nn),t2(g_nn)
    real(dp):: XX(g_nn,g_npp), ZZ(g_nn,g_npl),tt1(g_no,g_ni), tgrd2(g_npp+g_npl,g_ni)
    integer:: i 
    call zipt_deriv(beta,gamma,ll1,ll2,ll3,l1b,l1g,l2b,l3g,l1bb,l1gb,l1gg,l2bb,l3gg)
    t2=(g_zero*l1b+g_nzero*l2b)*g_miss
    t1=t2*g_gwt
    tt1=reshape(t1,(/g_no,g_ni/))
    tgrd(1:g_npp,:)=g_X(1:g_no,:).tmm.tt1
    tt1=reshape(t2,(/g_no,g_ni/))
    tgrd2(1:g_npp,:)=g_X(1:g_no,:).tmm.tt1
    grd(1:g_npp)=sum(tgrd(1:g_npp,:),2)
    t2=(g_zero*l1g-l3g)*g_miss
    t1=t2*g_gwt
    tt1=reshape(t1,(/g_no,g_ni/))
    tgrd((g_npp+1):(g_npp+g_npl),:)=g_Z(1:g_no,:).tmm.tt1
    tt1=reshape(t2,(/g_no,g_ni/))
    tgrd2((g_npp+1):(g_npp+g_npl),:)=g_Z(1:g_no,:).tmm.tt1
    grd((g_npp+1):(g_npp+g_npl))=sum(tgrd((g_npp+1):(g_npp+g_npl),:),2)
    t2=(g_zero*l1bb+g_nzero*l2bb)*g_gwt*g_miss
    XX=g_X
    do i=1,g_npp
       XX(:,i)=g_X(:,i)*t2
    end do
    H(1:g_npp,1:g_npp)=XX.tmm.g_X
    t2=(g_zero*l1gg-l3gg)*g_gwt*g_miss
    ZZ=g_Z
    do i=1,g_npl
       ZZ(:,i)=g_Z(:,i)*t2
    end do
    H((g_npp+1):(g_npp+g_npl),(g_npp+1):(g_npp+g_npl))=ZZ.tmm.g_Z
    t2=(g_zero*l1gb)*g_gwt*g_miss
    XX=g_X
    do i=1,g_npp
       XX(:,i)=g_X(:,i)*t2
    end do
    H(1:g_npp,(g_npp+1):(g_npp+g_npl))=XX.tmm.g_Z
    H((g_npp+1):(g_npp+g_npl),1:g_npp)=transpose(H(1:g_npp,(g_npp+1):(g_npp+g_npl)))
    H=H+(tgrd.mmt.tgrd2)
    H=-H
  end subroutine gdHss

  subroutine update_grd(beta,gamma,prob,grd,H) 
    use dat_mod
    use matrix, only: operator(.mm.), operator(.tmm.), operator(.mmt.)
    implicit none
    real(dp), intent(in):: beta(g_npp,g_ng), gamma(g_npl,g_ng), prob(g_ng)
    real(dp), intent(out):: grd(g_npg+g_ng-1), H(g_npg+g_ng-1,g_npg+g_ng-1)
    real(dp):: tgrd(g_npg+g_ng-1,g_ni), ggt(g_ni,g_ng), tHup(g_npg,g_ng-1), prb2(g_ng-1,g_ng-1)
    integer:: i, lo, hi
    H=0.0_dp
    tgrd=0.0_dp
    ggt=g_pr_wt(g_indi,:)
    do i=1,g_ng
       g_gwt => g_pr_wt(:,i)
       lo=(i-1)*(g_npp+g_npl)+1
       hi=i*(g_npp+g_npl)
       call gdHss(beta(:,i),gamma(:,i),grd(lo:hi),tgrd(lo:hi,:),H(lo:hi,lo:hi))
    end do
    do i=1,(g_ng-1)
       grd(g_npg+i)=sum(ggt(:,i)-prob(i))
       lo=(i-1)*(g_npp+g_npl)+1
       hi=i*(g_npp+g_npl)
       H(g_npg+i,lo:hi)=-grd(lo:hi)   
       H(lo:hi,g_npg+i)=-grd(lo:hi)   
       H(g_npg+i,g_npg+i)=-sum(ggt(:,i)-prob(i))
       prb2(:,i)=prob(1:(g_ng-1))*prob(i)
    end do
    prb2=g_ni*prb2
    tHup=tgrd(1:g_npg,:).mm.ggt(:,1:(g_ng-1))
    lo=g_npg+1
    hi=g_npg+g_ng-1
    H(1:g_npg,lo:hi)=H(1:g_npg,lo:hi)+tHup
    H(lo:hi,1:g_npg)=transpose(H(1:g_npg,lo:hi))
    H(lo:hi,lo:hi)=H(lo:hi,lo:hi)+(ggt(:,1:(g_ng-1)).tmm.ggt(:,1:(g_ng-1)))-prb2
    H=H+(tgrd.mmt.tgrd)
  end subroutine update_grd

  subroutine e_step_no(prob,llike) 
    use dat_mod
    implicit none
    real(dp), intent(in):: prob(g_ng)
    real(dp), intent(out):: llike
    real(dp):: dplike(g_ni,g_ng), rscl(g_ni), dlike(g_ni)
    integer:: i
    do i=1,g_ni
       dplike(i,:)=sum(g_llike_t((g_no*(i-1)+1):(g_no*i),:),1)
    end do
    do i=1,g_ng
       dplike(:,i)=log(prob(i))+dplike(:,i)
    end do
    rscl=maxval(dplike,2)
    do i=1,g_ng
       dplike(:,i)=exp(dplike(:,i)-rscl)
    end do
    dlike=sum(dplike,2)
    do i=1,g_ng
       dplike(:,i)=dplike(:,i)/dlike
    end do
    llike=sum(log(dlike)+rscl)
  end subroutine e_step_no

  subroutine llike_update_no(prob,beta,gamma,llike) 
    use dat_mod
    use dmZIP_shared_mod
    implicit none
    real(dp), intent(in):: prob(g_ng), beta(g_npp,g_ng), gamma(g_npl,g_ng)
    real(dp), intent(out):: llike
    call update_gllike(beta,gamma) 
    call e_step_no(prob,llike)    
  end subroutine llike_update_no
  
  subroutine pfun(n,x,llike) 
    use dat_mod
    implicit none
    integer, intent(in):: n
    real(dp), intent(in):: x(n)
    real(dp), intent(out):: llike
    real(dp):: prob(g_ng), beta(g_npp,g_ng), gamma(g_npl,g_ng), param(g_npp+g_npl+1,g_ng)
    param=reshape(x,(/g_npp+g_npl+1,g_ng/))
    beta=param(1:g_npp,:)
    gamma=param((g_npp+1):(g_npp+g_npl),:)
    prob=param(g_npp+g_npl+1,:)
    prob=prob/sum(prob)
    beta=18.0_dp*beta-9.0_dp
    gamma=18.0_dp*gamma-9.0_dp
    call llike_update_no(prob,beta,gamma,llike)
    if (minval(prob) < 2.0_dp/real(g_ni,dp)) then
       llike=1.0e20_dp
       return
    end if
    llike=-llike
  end subroutine pfun

  subroutine llikefn(n,x,llike,grd,H) 
    use dat_mod
    use dmZIP_shared_mod
    implicit none
    integer, intent(in):: n
    real(dp), intent(in):: x(n)
    real(dp), intent(out):: llike, grd(n), H(n,n)
    real(dp):: prob(g_ng), beta(g_npp,g_ng), gamma(g_npl,g_ng), betat(g_npp+g_npl,g_ng)
    betat=reshape(x(1:g_npg),(/g_npp+g_npl,g_ng/))
    beta=betat(1:g_npp,:)
    gamma=betat((g_npp+1):(g_npp+g_npl),:)
    prob(1:(g_ng-1))=exp(x((g_npg+1):n))    
    prob(g_ng)=1.0_dp
    prob=prob/sum(prob)
    call llike_update(prob,beta,gamma,llike)
    call update_grd(beta,gamma,prob,grd,H)
  end subroutine llikefn

end module dmZIP_mod

subroutine R_dmZIP_init_param(X,Z,Dat,offt,pparam,pllike,ni,no,npp,npl,ng,npop)
  use dat_mod
  use merge_sort_mod
  use dmZIP_mod
  use matrix, only: operator(.kp.)
  use gamma_mod
  use rrand
  use DE_mod
  implicit none
  integer, intent(in):: ni, no, npp, npl, ng, npop
  real(dp), intent(in):: X(no,npp), Z(no,npl), Dat(ni,no), offt(ni,no)
  real(dp), intent(out):: pparam(npop,(npp+npl+1)*ng), pllike(npop)
  integer:: nparm, nfeval, itermax
  real(dp):: paramm(npp+npl+1,ng)
  integer, parameter :: strategy=3, refresh=20
  integer, dimension(2), parameter:: method=(/0, 1 /)
  real(dp), parameter:: VTR=0.0_dp, CR_XC=0.5_dp, F_XC=0.8_dp
  real(dp):: F_CR=0.8_dp, xmin((npp+npl+1)*ng), xmax((npp+npl+1)*ng)
  integer:: i, irank(npop)
  g_nn=ni*no
  g_ni=ni
  g_no=no
  g_npp=npp
  g_npl=npl
  g_ng=ng
  allocate(g_X(g_nn,g_npp),g_Z(g_nn,g_npl),g_y(g_nn),g_llike_t(g_nn,g_ng),g_zero(g_nn),g_nzero(g_nn), &
       &g_llc(g_nn),g_expd(g_ni+g_no),g_miss(g_nn),g_offt(g_nn))
  g_expd=1.0_dp
  g_X=g_expd(1:g_ni).kp.X
  g_Z=g_expd(1:g_ni).kp.Z
  g_y=pack(transpose(Dat),.TRUE.)
  g_offt=pack(transpose(offt),.TRUE.)
  g_offt=log(g_offt)
  g_zero=0.0_dp
  g_miss=1.0_dp
  where (g_y < 0.5_dp)     
     g_zero=1.0_dp
  end where
  where (g_y < 0.0_dp)     
     g_miss=0.0_dp
     g_y=0.0_dp
  end where
  g_nzero=1.0_dp-g_zero
  g_llc=g_nzero*lgammafn(g_y+1.0_dp)
  itermax=100*ng
  nparm=(g_npp+g_npl+1)*g_ng
  xmin=1.0e-9_dp
  xmax=1.0_dp
  call r_set_random
  call DE_optim(pfun,nparm,xmin,xmax,VTR,npop,itermax,F_XC,CR_XC,strategy,refresh,pparam,pllike, &
       & nfeval,F_CR,method)
  call r_close_random
  call mrgrnk(pllike,irank)
  pllike=pllike(irank)
  pparam=pparam(irank,:)
  pllike=-pllike
  do i=1,npop
     paramm=reshape(pparam(i,:),(/g_npp+g_npl+1,g_ng/))
     paramm(1:(g_npp+g_npl),:)=18.0_dp*paramm(1:(g_npp+g_npl),:)-9.0_dp
     paramm(g_npp+g_npl+1,:)=paramm(g_npp+g_npl+1,:)/sum(paramm(g_npp+g_npl+1,:))
     pparam(i,:)=pack(paramm,.TRUE.)
  end do
  nullify(g_llikei)
  deallocate(g_X,g_Z,g_y,g_llike_t,g_zero,g_nzero,g_llc,g_expd,g_miss,g_offt)
end subroutine R_dmZIP_init_param

subroutine R_dmZIP(X,Z,Dat,offt,ggt,prob,beta,gamma,llike,Hess,nn,ni,no,npp,npl,ng,ierr) 
  use dat_mod
  use matrix, only: operator(.kp.)
  use Newton_mod
  use gamma_mod
  use merge_sort_mod
  use dmZIP_shared_mod
  use dmZIP_EM_mod
  use dmZIP_mod
  implicit none
  integer, intent(in):: nn, ni, no, npp, npl, ng
  real(dp), intent(in):: X(no,npp), Z(no,npl), Dat(ni,no), offt(ni,no) 
  real(dp), intent(inout):: prob(ng), beta(npp,ng), gamma(npl,ng)
  real(dp), intent(out):: ggt(ni,ng), llike, Hess((npp+npl)*ng+ng-1,(npp+npl)*ng+ng-1)
  integer, intent(out):: ierr
  real(dp):: param((npp+npl)*ng+ng-1), grd((npp+npl)*ng+ng-1), betat(npp+npl,ng)
  integer:: i, irank(ng)
  g_nn=nn
  g_ni=ni
  g_no=no
  g_npp=npp
  g_npl=npl
  g_ng=ng
  g_npg=g_ng*(g_npp+g_npl)
  allocate(g_X(g_nn,g_npp),g_Z(g_nn,g_npl),g_y(g_nn),g_zero(g_nn),g_nzero(g_nn),g_llc(g_nn), &
       &g_pr_wt(g_nn,g_ng),g_llike_t(g_nn,g_ng),g_expd(g_ni+g_no),g_miss(g_nn),g_offt(g_nn), &
       &g_indi(g_ni))
  g_indi=(/(i,i=1,g_ni)/)
  g_indi=g_indi*g_no
  g_expd=1.0_dp
  g_X=g_expd(1:g_ni).kp.X
  g_Z=g_expd(1:g_ni).kp.Z
  g_y=pack(transpose(Dat),.TRUE.)
  g_offt=pack(transpose(offt),.TRUE.)
  g_offt=log(g_offt)
  g_zero=0.0_dp
  g_miss=1.0_dp
  where (g_y < 0.5_dp)     
     g_zero=1.0_dp
  end where
  where (g_y < 0.0_dp)     
     g_miss=0.0_dp
     g_y=0.0_dp
  end where
  g_nzero=1.0_dp-g_zero
  g_llc=g_nzero*lgammafn(g_y+1.0_dp)
  if (minval(prob) < 2.0_dp/real(g_ni,dp)) then
!     write(*,*) "R_dmZIP: prob -> 0!"
     nullify(g_gwt,g_llikei)
     deallocate(g_X,g_Z,g_y,g_pr_wt,g_llike_t,g_zero,g_nzero,g_llc,g_expd,g_miss,g_offt,g_indi)
     ierr=-99
     return
  end if
  if (g_ng > 1) then
     betat(1:g_npp,:)=beta
     betat((g_npp+1):(g_npp+g_npl),:)=gamma
     param(1:g_npg)=pack(betat,.TRUE.)
     param((g_npg+1):(g_npg+g_ng-1))=log(prob(1:(g_ng-1))/prob(g_ng))
     isEM=.FALSE.
     call Newton(llikefn,param,llike,grd,Hess,g_npg+g_ng-1,ierr)
  else
     ierr=-99
  end if
  if (ierr /= 0) then
     ierr=0
     isEM=.TRUE.
     call em_update(prob,beta,gamma,llike,ierr) 
     if (ierr /= 0) then
        nullify(g_gwt,g_llikei)
        deallocate(g_X,g_Z,g_y,g_pr_wt,g_llike_t,g_zero,g_nzero,g_llc,g_expd,g_miss,g_offt,g_indi)
        return
     end if
  else 
     betat=reshape(param(1:g_npg),(/g_npp+g_npl,g_ng/))
     beta=betat(1:g_npp,:)
     gamma=betat((g_npp+1):(g_npp+g_npl),:)
     prob(1:(g_ng-1))=exp(param((g_npg+1):(g_npg+g_ng-1)))
     prob(g_ng)=1.0_dp
     prob=prob/sum(prob)
  end if
  call llike_update(prob,beta,gamma,llike)
  call mrgrnk(prob,irank)
  beta=beta(:,irank)
  gamma=gamma(:,irank)
  prob=prob(irank)
  g_pr_wt=g_pr_wt(:,irank)
  ggt=g_pr_wt(g_indi,:)
  nullify(g_gwt,g_llikei)
  deallocate(g_X,g_Z,g_y,g_pr_wt,g_llike_t,g_zero,g_nzero,g_llc,g_expd,g_miss,g_offt,g_indi)
end subroutine R_dmZIP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Copyright (C) 2023 Jason D. Nielsen
!!!
!!! This program is free software! you can redistribute it and/or modify
!!! it under the terms of the GNU General Public License as published by
!!! the Free Software Foundation! either version 2 of the License, or
!!! (at your option) any later version.
!!!
!!! This program is distributed in the hope that it will be useful,
!!! but WITHOUT ANY WARRANTY! without even the implied warranty of
!!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!! GNU General Public License for more details.
!!!
!!! You should have received a copy of the GNU General Public License
!!! along with this program! if not, write to the Free Software
!!! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111, USA
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module dmZIPt_shared_mod
  implicit none
contains

  subroutine zipt_llike(beta,tau,llike)
    use dat_mod
    use matrix, only: operator(.mm.)
    implicit none
    real(dp), intent(in):: beta(g_npp), tau
    real(dp), intent(out):: llike(g_nn)
    real(dp):: t1,t2(g_nn),t3(g_nn),t4(g_nn),t5(g_nn),t6(g_nn),t7(g_nn),t9(g_nn),&
         & l1(g_nn),l2(g_nn),l3(g_nn)
    real(dp)::  xb(g_nn)
    integer:: i
    xb(1:g_no)=g_X(1:g_no,:).mm.beta
    t1 = exp(tau)
    t2(1:g_no) = t1 * xb(1:g_no)
    t3(1:g_no) = exp(-t2(1:g_no))
    t9(1:g_no) = 1.0_dp + t3(1:g_no)
    l3(1:g_no) = log(t9(1:g_no))
    do i=2,g_ni 
       xb((g_no*(i-1)+1):(g_no*i)) = xb(1:g_no)
       t2((g_no*(i-1)+1):(g_no*i)) = t2(1:g_no)
       t3((g_no*(i-1)+1):(g_no*i)) = t3(1:g_no)
       t9((g_no*(i-1)+1):(g_no*i)) = t9(1:g_no)
       l3((g_no*(i-1)+1):(g_no*i)) = l3(1:g_no)
    end do
    t4 = g_offt + xb
    t5 = exp(t4)
    t6 = exp(-t5)
    t7 = t3 + t6
    l1 = log(t7)
    l2 = g_y * t4 - t5
    llike=g_zero*l1+g_nzero*l2-l3-g_llc
    llike=llike*g_miss
  end subroutine zipt_llike

  subroutine zipt_deriv(beta,tau,ll1,ll2,ll3,l1b,l1t,l2b,l3b,l3t,l1bb,l1tb,l1tt,l2bb,l3bb,l3tb,l3tt)
    use dat_mod
    use matrix, only: operator(.mm.)
    implicit none
    real(dp), intent(in):: beta(g_npp), tau
    real(dp), intent(out):: ll1(g_nn),ll2(g_nn),ll3(g_nn),l1b(g_nn),l1t(g_nn),l2b(g_nn), &
         & l3b(g_nn),l3t(g_nn),l1bb(g_nn),l1tb(g_nn),l1tt(g_nn),l2bb(g_nn),l3bb(g_nn), &
         l3tb(g_nn),l3tt(g_nn)
    real(dp):: t1(g_nn),t2(g_nn),t3(g_nn),t4(g_nn),t5(g_nn),t6(g_nn),t7(g_nn),t9(g_nn),&
         & t10(g_nn),t11(g_nn),t12(g_nn),t13(g_nn),t14(g_nn),t16(g_nn),t17(g_nn),t18(g_nn),&
         & t22(g_nn),t23(g_nn),t24(g_nn),t27(g_nn),t32(g_nn),t34(g_nn),t37(g_nn),t39(g_nn),&
         & t40(g_nn),t43(g_nn),t45(g_nn),t46(g_nn)
    real(dp)::  xb(g_nn)
    xb=g_X.mm.beta
    t1 = exp(tau)
    t2 = t1 * xb
    t3 = exp(-t2)
    t4 = g_offt + xb
    t5 = exp(t4)
    t6 = exp(-t5)
    t7 = t3 + t6
    t9 = 1.0_dp + t3
    t10 = t1 * t3
    t11 = t5 * t6
    t12 = -t10 - t11
    t13 = 1.0_dp / t7
    t14 = t3 * t13
    t16 = t1 ** 2
    t17 = t16 * t3
    t18 = t5 ** 2
    t22 = t12 ** 2
    t23 = t7 ** 2
    t24 = 1.0_dp / t23
    t27 = t16 * xb
    t32 = 1.0_dp / t9
    t34 = t3 * t32
    t37 = t3 ** 2
    t39 = t9 ** 2
    t40 = 1.0_dp / t39
    t43 = t37 * t40
    t45 = xb ** 2
    t46 = t16 * t45
    ll1 = log(t7)
    ll2 = g_y * t4 - t5
    ll3 = log(t9)
    l1b = t12 * t13
    l1t = -t2 * t14
    l2b = g_y - t5
    l3b = -t10 * t32
    l3t = -t2 * t34
    l1bb = (t17 - t11 + t18 * t6) * t13 - t22 * t24
    l1tb = -t10 * t13 + t27 * t14 + t2 * t3 * t24 * t12
    l1tt = l1t + t46 * t14 - t46 * t37 * t24
    l2bb = -t5
    l3bb = t17 * t32 - t16 * t37 * t40
    l3tb = l3b + t27 * t34 - t27 * t43
    l3tt = l3t + t46 * t34 - t46 * t43
  end subroutine zipt_deriv

  subroutine update_gllike(beta,tau) 
    use dat_mod
    implicit none
    real(dp), intent(in):: beta(g_npp,g_ng), tau(g_ng)
    integer:: i
    do i=1,g_ng
       g_llikei => g_llike_t(:,i)
       call zipt_llike(beta(:,i),tau(i),g_llikei)
    end do
  end subroutine update_gllike

  subroutine llike_update(prob,beta,tau,llike) 
    use dat_mod
    implicit none
    real(dp), intent(inout):: prob(g_ng)
    real(dp), intent(in):: beta(g_npp,g_ng), tau(g_ng)
    real(dp), intent(out):: llike
    call update_gllike(beta,tau) 
    call e_step(prob,llike)    
  end subroutine llike_update

  subroutine e_step(prob,llike) 
    use dat_mod
    use matrix, only: operator(.kp.)
    implicit none
    real(dp), intent(inout):: prob(g_ng)
    real(dp), intent(out):: llike
    real(dp):: dplike(g_ni,g_ng), rscl(g_ni), dlike(g_ni)
    integer:: i
    do i=1,g_ni
       dplike(i,:)=sum(g_llike_t((g_no*(i-1)+1):(g_no*i),:),1)
    end do
    do i=1,g_ng
       dplike(:,i)=log(prob(i))+dplike(:,i)
    end do
    rscl=maxval(dplike,2)
    do i=1,g_ng
       dplike(:,i)=exp(dplike(:,i)-rscl)
    end do
    dlike=sum(dplike,2)
    do i=1,g_ng
       dplike(:,i)=dplike(:,i)/dlike
    end do
    if (isEM) then
       prob=sum(dplike,1)/real(g_ni,dp)
    end if
    llike=sum(log(dlike)+rscl)
    g_pr_wt=dplike.kp.g_expd(1:g_no)
  end subroutine e_step

end module dmZIPt_shared_mod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Copyright (C) 2023 Jason D. Nielsen
!!!
!!! This program is free software! you can redistribute it and/or modify
!!! it under the terms of the GNU General Public License as published by
!!! the Free Software Foundation! either version 2 of the License, or
!!! (at your option) any later version.
!!!
!!! This program is distributed in the hope that it will be useful,
!!! but WITHOUT ANY WARRANTY! without even the implied warranty of
!!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!! GNU General Public License for more details.
!!!
!!! You should have received a copy of the GNU General Public License
!!! along with this program! if not, write to the Free Software
!!! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111, USA
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module dmZIPt_EM_mod
  implicit none
  private
  public:: em_update

contains

  subroutine gdHss(beta,tau,grd,H)
    use dat_mod
    use dmZIPt_shared_mod
    use matrix, only: operator(.mm.), operator(.tmm.)
    implicit none
    real(dp), intent(inout):: beta(g_npp), tau
    real(dp), intent(out):: H(g_npp+1,g_npp+1), grd(g_npp+1)
    real(dp):: t1(g_nn),t2(g_nn),t3(g_nn),t4(g_nn)
    real(dp):: ll1(g_nn),ll2(g_nn),ll3(g_nn),l1b(g_nn),l1t(g_nn),l2b(g_nn),l3b(g_nn),l3t(g_nn), &
         & l1bb(g_nn),l1tb(g_nn),l1tt(g_nn),l2bb(g_nn),l3bb(g_nn),l3tb(g_nn),l3tt(g_nn)
    real(dp):: XX(g_nn,g_npp)
    integer:: i 
    call zipt_deriv(beta,tau,ll1,ll2,ll3,l1b,l1t,l2b,l3b,l3t,l1bb,l1tb,l1tt,l2bb,l3bb,l3tb,l3tt)
    t1=(g_zero*l1b+g_nzero*l2b-l3b)*g_gwt*g_miss
    grd(1:g_npp)=t1.mm.g_X
    grd(g_npp+1)=sum((g_zero*l1t-l3t)*g_gwt*g_miss)
    t2=(g_zero*l1bb+g_nzero*l2bb-l3bb)*g_gwt*g_miss
    XX=g_X
    do i=1,g_npp
       XX(:,i)=g_X(:,i)*t2
    end do
    H(1:g_npp,1:g_npp)=XX.tmm.g_X
    H(g_npp+1,g_npp+1)=sum((g_zero*l1tt-l3tt)*g_gwt*g_miss)
    t3=(g_zero*l1tb-l3tb)*g_gwt*g_miss
    t4(1:g_npp)=t3.mm.g_X
    H(1:g_npp,g_npp+1)=t4(1:g_npp)
    H(g_npp+1,1:g_npp)=t4(1:g_npp)
    H=-H
  end subroutine gdHss

  subroutine ziptfit_mle(beta,tau,nis,grd,H,ierr) 
    use dat_mod
    use matrix, only: fsymsolve!, rsymsolve
    implicit none
    integer, intent(in):: nis
    real(dp), intent(inout):: beta(g_npp), tau
    real(dp), intent(out):: grd(g_npp+1), H(g_npp+1,g_npp+1)
    integer, intent(out):: ierr
    real(dp):: gvalo, gval, param(g_npp+1), ss(g_npp+1)
    integer:: i
    ierr=0
    param=(/beta,tau/)
    do i=1,nis
       call gdHss(beta,tau,grd,H)
       call fsymsolve(H,grd,ss,ierr)
       if (ierr /= 0) then
!          write(*,*) "ziptfit_mle: Error in fsymsolve!"
          return
       end if
!       call rsymsolve(H,grd,ss)
       param=param+0.1_dp*ss
       beta=param(1:g_npp)
       tau=param(g_npp+1)
    end do
    gvalo=maxval(abs(grd))
    do i=1,nmaxit
       call gdHss(beta,tau,grd,H)
       gval=maxval(abs(grd))
       call fsymsolve(H,grd,ss,ierr)
       if (ierr /= 0) then
!          write(*,*) "ziptfit_mle: Error in fsymsolve!"
          return
       end if
!       call rsymsolve(H,grd,ss)
       if (gval-gvalo > 0) then
          param=param+0.5_dp*ss
       else
          param=param+ss
       end if
       if (gval < 1.0e-6_dp) then
          return
       end if
       beta=param(1:g_npp)
       tau=param(g_npp+1)
       gvalo=gval
       if (i==nmaxit) then
!          write(*,*) "ziptfit_mle: maximun number of interations reached!"
          ierr=-99
          return
       end if
    end do
  end subroutine ziptfit_mle
  
  subroutine ziptfit(nis,beta,tau,ierr) 
    use dat_mod
    implicit none
    integer, intent(in):: nis
    real(dp), intent(inout):: beta(g_npp), tau
    integer, intent(out):: ierr
    real(dp):: grd(g_npp+1), H(g_npp+1,g_npp+1)
    call ziptfit_mle(beta,tau,nis,grd,H,ierr) 
    if (ierr /= 0) then
!       write(*,*) "ziptfit: routine ziptfit_mle failed!"
       return
    end if
  end subroutine ziptfit

  subroutine m_step(nis,beta,tau,ierr) 
    use dat_mod
    implicit none
    integer, intent(in):: nis
    real(dp), intent(inout):: beta(g_npp,g_ng), tau(g_ng)
    integer, intent(out):: ierr
    integer:: i
    do i=1,g_ng
       g_gwt => g_pr_wt(:,i)
       call ziptfit(nis,beta(:,i),tau(i),ierr)
       if (ierr /= 0) then
!          write(*,*) "m_step: Error in zipfit!"
          return
       end if
    end do
  end subroutine m_step
  
  subroutine em_update(prob,beta,tau,llike,ierr) 
    use dat_mod
    use dmZIPt_shared_mod
    implicit none
    real(dp), intent(inout):: prob(g_ng), beta(g_npp,g_ng), tau(g_ng)
    real(dp), intent(out):: llike
    integer, intent(out):: ierr
    real(dp):: llikeo
    integer:: i 
    llikeo=-1.0e10_dp
    call llike_update(prob,beta,tau,llike)
    do i=1,nmaxit
       call m_step(5,beta,tau,ierr)
       if (ierr /= 0) then
!          write(*,*) "em_update: Error in m_step!"
          return
       end if
       call update_gllike(beta,tau) 
       call e_step(prob,llike)
       if (minval(prob) < 2.0_dp/real(g_ni,dp)) then
!          write(*,*) "em_update: prob -> 0!"
          ierr=-99
          return
       end if
       if (abs(llike-llikeo) < 1.0e-6_dp) then
          exit
       end if
       llikeo=llike
    end do
  end subroutine em_update

end module dmZIPt_EM_mod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Copyright (C) 2023 Jason D. Nielsen
!!!
!!! This program is free software! you can redistribute it and/or modify
!!! it under the terms of the GNU General Public License as published by
!!! the Free Software Foundation! either version 2 of the License, or
!!! (at your option) any later version.
!!!
!!! This program is distributed in the hope that it will be useful,
!!! but WITHOUT ANY WARRANTY! without even the implied warranty of
!!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!! GNU General Public License for more details.
!!!
!!! You should have received a copy of the GNU General Public License
!!! along with this program! if not, write to the Free Software
!!! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111, USA
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module dmZIPt_mod
  implicit none
  private
  public:: pfun, llikefn
contains

  subroutine gdHss(beta,tau,grd,tgrd,H)
    use dat_mod
    use dmZIPt_shared_mod
    use matrix, only: operator(.mm.), operator(.tmm.), operator(.mmt.)
    implicit none
    real(dp), intent(in):: beta(g_npp), tau
    real(dp), intent(out):: H(g_npp+1,g_npp+1), grd(g_npp+1), tgrd(g_npp+1,g_ni)
    real(dp):: ll1(g_nn),ll2(g_nn),ll3(g_nn),l1b(g_nn),l1t(g_nn),l2b(g_nn),l3b(g_nn),l3t(g_nn), &
         & l1bb(g_nn),l1tb(g_nn),l1tt(g_nn),l2bb(g_nn),l3bb(g_nn),l3tb(g_nn),l3tt(g_nn)
    real(dp):: t1(g_nn),t2(g_nn),t3(g_nn),t4(g_nn)
    real(dp):: XX(g_nn,g_npp), tt1(g_no,g_ni), tgrd2(g_npp+1,g_ni)
    integer:: i 
    call zipt_deriv(beta,tau,ll1,ll2,ll3,l1b,l1t,l2b,l3b,l3t,l1bb,l1tb,l1tt,l2bb,l3bb,l3tb,l3tt)
    t2=(g_zero*l1b+g_nzero*l2b-l3b)*g_miss
    t1=t2*g_gwt
    tt1=reshape(t1,(/g_no,g_ni/))
    tgrd(1:g_npp,:)=g_X(1:g_no,:).tmm.tt1
    tt1=reshape(t2,(/g_no,g_ni/))
    tgrd2(1:g_npp,:)=g_X(1:g_no,:).tmm.tt1
    grd(1:g_npp)=sum(tgrd(1:g_npp,:),2)
    t2=(g_zero*l1t-l3t)*g_miss
    t1=t2*g_gwt
    tt1=reshape(t1,(/g_no,g_ni/))
    tgrd(g_npp+1,:)=sum(tt1,1)
    tt1=reshape(t2,(/g_no,g_ni/))
    tgrd2(g_npp+1,:)=sum(tt1,1)
    grd(g_npp+1)=sum(tgrd(g_npp+1,:))
    t2=(g_zero*l1bb+g_nzero*l2bb-l3bb)*g_gwt*g_miss
    XX=g_X
    do i=1,g_npp
       XX(:,i)=g_X(:,i)*t2
    end do
    H(1:g_npp,1:g_npp)=XX.tmm.g_X
    H(g_npp+1,g_npp+1)=sum((g_zero*l1tt-l3tt)*g_gwt*g_miss)
    t3=(g_zero*l1tb-l3tb)*g_gwt*g_miss
    t4(1:g_npp)=t3.mm.g_X
    H(1:g_npp,g_npp+1)=t4(1:g_npp)
    H(g_npp+1,1:g_npp)=t4(1:g_npp)
    H=H+(tgrd.mmt.tgrd2)
    H=-H
  end subroutine gdHss

  subroutine update_grd(beta,tau,prob,grd,H) 
    use dat_mod
    use matrix, only: operator(.mm.), operator(.tmm.), operator(.mmt.)
    implicit none
    real(dp), intent(in):: beta(g_npp,g_ng), tau(g_ng), prob(g_ng)
    real(dp), intent(out):: grd(g_npg+2*g_ng-1), H(g_npg+2*g_ng-1,g_npg+2*g_ng-1)
    real(dp):: tgrd(g_npg+2*g_ng-1,g_ni), ggt(g_ni,g_ng), tHup(g_npg+g_ng,g_ng-1), prb2(g_ng-1,g_ng-1)
    integer:: i, lo, hi
    H=0.0_dp
    tgrd=0.0_dp
    ggt=g_pr_wt(g_indi,:)
    do i=1,g_ng
       g_gwt => g_pr_wt(:,i)
       lo=(i-1)*(g_npp+1)+1
       hi=i*(g_npp+1)
       call gdHss(beta(:,i),tau(i),grd(lo:hi),tgrd(lo:hi,:),H(lo:hi,lo:hi))
    end do
    do i=1,(g_ng-1)
       grd(g_npg+g_ng+i)=sum(ggt(:,i)-prob(i))
       lo=(i-1)*(g_npp+1)+1
       hi=i*(g_npp+1)
       H(g_npg+g_ng+i,lo:hi)=-grd(lo:hi)   
       H(lo:hi,g_npg+g_ng+i)=-grd(lo:hi)   
       H(g_npg+g_ng+i,g_npg+g_ng+i)=-sum(ggt(:,i)-prob(i))
       prb2(:,i)=prob(1:(g_ng-1))*prob(i)
    end do
    prb2=g_ni*prb2
    tHup=tgrd(1:(g_npg+g_ng),:).mm.ggt(:,1:(g_ng-1))
    lo=g_npg+g_ng+1
    hi=g_npg+2*g_ng-1
    H(1:(g_npg+g_ng),lo:hi)=H(1:(g_npg+g_ng),lo:hi)+tHup
    H(lo:hi,1:(g_npg+g_ng))=transpose(H(1:(g_npg+g_ng),lo:hi))
    H(lo:hi,lo:hi)=H(lo:hi,lo:hi)+(ggt(:,1:(g_ng-1)).tmm.ggt(:,1:(g_ng-1)))-prb2
    H=H+(tgrd.mmt.tgrd)
  end subroutine update_grd

  subroutine e_step_no(prob,llike) 
    use dat_mod
    implicit none
    real(dp), intent(in):: prob(g_ng)
    real(dp), intent(out):: llike
    real(dp):: dplike(g_ni,g_ng), rscl(g_ni), dlike(g_ni)
    integer:: i
    do i=1,g_ni
       dplike(i,:)=sum(g_llike_t((g_no*(i-1)+1):(g_no*i),:),1)
    end do
    do i=1,g_ng
       dplike(:,i)=log(prob(i))+dplike(:,i)
    end do
    rscl=maxval(dplike,2)
    do i=1,g_ng
       dplike(:,i)=exp(dplike(:,i)-rscl)
    end do
    dlike=sum(dplike,2)
    do i=1,g_ng
       dplike(:,i)=dplike(:,i)/dlike
    end do
    llike=sum(log(dlike)+rscl)
  end subroutine e_step_no

  subroutine llike_update_no(prob,beta,tau,llike) 
    use dat_mod
    use dmZIPt_shared_mod
    implicit none
    real(dp), intent(in):: prob(g_ng), beta(g_npp,g_ng), tau(g_ng)
    real(dp), intent(out):: llike
    call update_gllike(beta,tau) 
    call e_step_no(prob,llike)    
  end subroutine llike_update_no
  
  subroutine pfun(n,x,llike) 
    use dat_mod
    implicit none
    integer, intent(in):: n
    real(dp), intent(in):: x(n)
    real(dp), intent(out):: llike
    real(dp):: prob(g_ng), beta(g_npp,g_ng), tau(g_ng), param(g_npp+2,g_ng)
    param=reshape(x,(/g_npp+2,g_ng/))
    beta=param(1:g_npp,:)
    tau=param(g_npp+1,:)
    prob=param(g_npp+2,:)
    prob=prob/sum(prob)
    tau=log(tau)
    beta=18.0_dp*beta-9.0_dp
    call llike_update_no(prob,beta,tau,llike)
    if (minval(prob) < 2.0_dp/real(g_ni,dp)) then
       llike=1.0e20_dp
       return
    end if
    llike=-llike
  end subroutine pfun

  subroutine llikefn(n,x,llike,grd,H) 
    use dat_mod
    use dmZIPt_shared_mod
    implicit none
    integer, intent(in):: n
    real(dp), intent(in):: x(n)
    real(dp), intent(out):: llike, grd(n), H(n,n)
    real(dp):: prob(g_ng), beta(g_npp,g_ng), tau(g_ng), betat(g_npp+1,g_ng)
    betat=reshape(x(1:(g_npg+g_ng)),(/g_npp+1,g_ng/))
    beta=betat(1:g_npp,:)
    tau=betat(g_npp+1,:)
    prob(1:(g_ng-1))=exp(x((g_npg+g_ng+1):n))    
    prob(g_ng)=1.0_dp
    prob=prob/sum(prob)
    call llike_update(prob,beta,tau,llike)
    call update_grd(beta,tau,prob,grd,H)
  end subroutine llikefn

end module dmZIPt_mod

subroutine R_dmZIPt_init_param(X,Dat,offt,pparam,pllike,ni,no,npp,ng,npop)
  use dat_mod
  use merge_sort_mod
  use dmZIPt_mod
  use matrix, only: operator(.kp.)
  use gamma_mod
  use rrand
  use DE_mod
  implicit none
  integer, intent(in):: ni, no, npp, ng, npop
  real(dp), intent(in):: X(no,npp), Dat(ni,no), offt(ni,no)
  real(dp), intent(out):: pparam(npop,(npp+2)*ng), pllike(npop)
  integer:: nparm, nfeval, itermax
  real(dp):: paramm(npp+2,ng)
  integer, parameter :: strategy=3, refresh=20
  integer, dimension(2), parameter:: method=(/0, 1 /)
  real(dp), parameter:: VTR=0.0_dp, CR_XC=0.5_dp, F_XC=0.8_dp
  real(dp):: F_CR=0.8_dp, xmin((npp+2)*ng), xmax((npp+2)*ng)
  integer:: i, irank(npop)
  g_nn=ni*no
  g_ni=ni
  g_no=no
  g_npp=npp
  g_ng=ng
  allocate(g_X(g_nn,g_npp),g_y(g_nn),g_zero(g_nn),g_nzero(g_nn),g_llc(g_nn),g_llike_t(g_nn,g_ng), &
       & g_expd(g_ni+g_no),g_miss(g_nn),g_offt(g_nn))
  g_expd=1.0_dp
  g_X=g_expd(1:g_ni).kp.X
  g_y=pack(transpose(Dat),.TRUE.)
  g_offt=pack(transpose(offt),.TRUE.)
  g_offt=log(g_offt)
  g_zero=0.0_dp
  g_miss=1.0_dp
  where (g_y < 0.5_dp)     
     g_zero=1.0_dp
  end where
  where (g_y < 0.0_dp)     
     g_miss=0.0_dp
     g_y=0.0_dp
  end where
  g_nzero=1.0_dp-g_zero
  g_llc=g_nzero*lgammafn(g_y+1.0_dp)
  itermax=100*ng
  nparm=(g_npp+2)*g_ng
  xmin=1.0e-9_dp
  xmax=1.0_dp
  call r_set_random
  call DE_optim(pfun,nparm,xmin,xmax,VTR,npop,itermax,F_XC,CR_XC,strategy,refresh,pparam,pllike, &
       & nfeval,F_CR,method)
  call r_close_random
  call mrgrnk(pllike,irank)
  pllike=pllike(irank)
  pparam=pparam(irank,:)
  pllike=-pllike
  do i=1,npop
     paramm=reshape(pparam(i,:),(/g_npp+2,g_ng/))
     paramm(1:g_npp,:)=18.0_dp*paramm(1:g_npp,:)-9.0_dp
     paramm(g_npp+1,:)=log(paramm(g_npp+1,:))
     paramm(g_npp+2,:)=paramm(g_npp+2,:)/sum(paramm(g_npp+2,:))
     pparam(i,:)=pack(paramm,.TRUE.)
  end do
  nullify(g_llikei)
  deallocate(g_X,g_y,g_llike_t,g_zero,g_nzero,g_llc,g_expd,g_miss,g_offt)
end subroutine R_dmZIPt_init_param

subroutine R_dmZIPt(X,Dat,offt,ggt,prob,beta,tau,llike,Hess,nn,ni,no,npp,ng,ierr) 
  use dat_mod
  use matrix, only: operator(.kp.)
  use Newton_mod
  use gamma_mod
  use merge_sort_mod
  use dmZIPt_shared_mod
  use dmZIPt_EM_mod
  use dmZIPt_mod
  implicit none
  integer, intent(in):: nn, ni, no, npp, ng
  real(dp), intent(in):: X(no,npp), Dat(ni,no), offt(ni,no) 
  real(dp), intent(inout):: prob(ng), beta(npp,ng), tau(ng)
  real(dp), intent(out):: ggt(ni,ng), llike, Hess(npp*ng+2*ng-1,npp*ng+2*ng-1)
  integer, intent(out):: ierr
  real(dp):: param(npp*ng+2*ng-1), grd(npp*ng+2*ng-1), betat(npp+1,ng)
  integer:: i, irank(ng)
  g_nn=nn
  g_ni=ni
  g_no=no
  g_npp=npp
  g_ng=ng
  g_npg=g_ng*g_npp
  allocate(g_X(g_nn,g_npp),g_y(g_nn),g_zero(g_nn),g_nzero(g_nn),g_llc(g_nn),g_pr_wt(g_nn,g_ng),&
       &g_llike_t(g_nn,g_ng),g_expd(g_ni+g_no),g_miss(g_nn),g_offt(g_nn),g_indi(g_ni))
  g_indi=(/(i,i=1,g_ni)/)
  g_indi=g_indi*g_no
  g_expd=1.0_dp
  g_X=g_expd(1:g_ni).kp.X
  g_y=pack(transpose(Dat),.TRUE.)
  g_offt=pack(transpose(offt),.TRUE.)
  g_offt=log(g_offt)
  g_zero=0.0_dp
  g_miss=1.0_dp
  where (g_y < 0.5_dp)     
     g_zero=1.0_dp
  end where
  where (g_y < 0.0_dp)     
     g_miss=0.0_dp
     g_y=0.0_dp
  end where
  g_nzero=1.0_dp-g_zero
  g_llc=g_nzero*lgammafn(g_y+1.0_dp)
  if (minval(prob) < 2.0_dp/real(g_ni,dp)) then
!     write(*,*) "R_dmZIPt: prob -> 0!"
     nullify(g_gwt,g_llikei)
     deallocate(g_X,g_y,g_pr_wt,g_llike_t,g_zero,g_nzero,g_llc,g_expd,g_miss,g_offt,g_indi)
     ierr=-99
     return
  end if
  if (g_ng > 1) then
     betat(1:g_npp,:)=beta
     betat(g_npp+1,:)=tau
     param(1:(g_npg+g_ng))=pack(betat,.TRUE.)
     param((g_npg+g_ng+1):(g_npg+2*g_ng-1))=log(prob(1:(g_ng-1))/prob(g_ng))
     isEM=.FALSE.
     call Newton(llikefn,param,llike,grd,Hess,g_npg+2*g_ng-1,ierr)
  else
     ierr=-99
  end if
  if (ierr /= 0) then
     ierr=0
     isEM=.TRUE.
     call em_update(prob,beta,tau,llike,ierr) 
     if (ierr /= 0) then
        nullify(g_gwt,g_llikei)
        deallocate(g_X,g_y,g_pr_wt,g_llike_t,g_zero,g_nzero,g_llc,g_expd,g_miss,g_offt,g_indi)
        return
     end if
  else 
     betat=reshape(param(1:(g_npg+g_ng)),(/g_npp+1,g_ng/))
     beta=betat(1:g_npp,:)
     tau=betat(g_npp+1,:)
     prob(1:(g_ng-1))=exp(param((g_npg+g_ng+1):(g_npg+2*g_ng-1)))
     prob(g_ng)=1.0_dp
     prob=prob/sum(prob)
  end if
  call llike_update(prob,beta,tau,llike)
  call mrgrnk(prob,irank)
  beta=beta(:,irank)
  tau=tau(irank)
  prob=prob(irank)
  g_pr_wt=g_pr_wt(:,irank)
  ggt=g_pr_wt(g_indi,:)
  nullify(g_gwt,g_llikei)
  deallocate(g_X,g_y,g_pr_wt,g_llike_t,g_zero,g_nzero,g_llc,g_expd,g_miss,g_offt,g_indi)
end subroutine R_dmZIPt
