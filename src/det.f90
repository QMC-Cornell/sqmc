      subroutine determinant(a,nsub,det)
! Written by Kevin Schmidt or from some standard library?

      use types, only: rk
      implicit real*8(a-h,o-z)
      parameter (zero=0._rk,one=1._rk)

! routine to calculate inverse and determinant of matrix a
! assumed to be dimensioned a(nsub,nsub).
! the matrix a is replaced by its inverse.

      dimension a(*)
      dimension ipivot(nsub),atemp(nsub)

      n=nsub
!      write(6,*) 'Inverting this matrix:'
!      do irow=0,(n-1)
!         write(6,'(100f6.2)') (a(n*irow+k),k=1,n)
!      enddo
      do 5 i=1,n
    5   ipivot(i)=i

! initialize determinant
      determ=one

! loop through columns
      iclm=-n
      do 25 i=1,n
        iclm=iclm+n

! loop through rows and select row with largest element
        adiag=a(ipivot(i)+iclm)
        idiag=i
        do 15 k=i,n
          if(abs(a(ipivot(k)+iclm)).gt.abs(adiag)) then
            adiag=a(ipivot(k)+iclm)
            idiag=k
          endif
   15   continue

! interchange pointers if row different from
! original is selected and change sign of determinant because
! of interchange
        if(idiag.ne.i) then
          determ=-determ
          itemp=ipivot(i)
          ipivot(i)=ipivot(idiag)
          ipivot(idiag)=itemp
        endif

! update determinant
        determ=adiag*determ
        a(ipivot(i)+iclm)=one
        if(adiag.eq.0._rk) then
          det=0.0
          return
          !stop 'adiag=0 in determinant'
        endif
        adiagi=one/adiag

! scale row by inverse diagonal
        call scal(n,adiagi,a(ipivot(i)),n)

! loop through other rows
! if not current row, then row reduce
        do 20 j=1,n
          if(j.ne.ipivot(i)) then
            t=-a(j+iclm)
            a(j+iclm)=zero
            call axpy(n,t,a(ipivot(i)),n,a(j),n)
          endif
   20   continue
   25 continue

! interchange elements to unpivot inverse matrix
! the following is equivalent to:
!      anew(i,ipivot(j))=aold(ipivot(i),j)
      jn=-n
      do 52 j=1,n
        jn=jn+n
        do 40 i=1,n
   40     atemp(i)=a(i+jn)
        do 50 i=1,n
   50     a(i+jn)=atemp(ipivot(i))
   52 continue

      do 55 j=1,n
   55   ipivot(j)=(ipivot(j)-1)*n

      do 90 i=1,n
        jn=-n
        do 70 j=1,n
          jn=jn+n
   70     atemp(j)=a(i+jn)
        do 80 j=1,n
   80     a(i+ipivot(j))=atemp(j)
   90 continue

      det=determ

      return
      end

