      subroutine spline(x,y,n,yp1,ypn,y2)

      implicit double precision (a-h,o-z)
      parameter (nmax=230)
      dimension x(n),y(n),y2(n),u(nmax)

!c     if yp1>1.0e+30 use natural spline, otherwise estimate y2 at the
!c     first point

      if (yp1.gt..99d30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif

!c     store intermediate values of terms in the expansion series

      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
11    continue

!c     if ypn>1.0e+30 use natural spline, otherwise estimate y2 at the
!c     last point point

      if (ypn.gt..99d30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)

!c     compute the y2 from the 2nd order expansion series

      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue

      return
      end subroutine spline

!c..............................................................................

      subroutine splint(xa,ya,y2a,n,x,y)


      implicit double precision (a-h,o-z)
      dimension xa(n),ya(n),y2a(n)

      klo=1
      khi=n

!c     determine the indices of array xa that bracket the input x value

1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif

!c     determine the finite difference along the x dimension

      h=xa(khi)-xa(klo)
      if (h.eq.0.) stop 'bad xa input in routine spline.'

!c     interpolate

      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.

      return
      end subroutine splint

