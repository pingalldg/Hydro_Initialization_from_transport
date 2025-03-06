subroutine eigenvalue(xi,yi,k,tau_ampt,Tmunu,umu,en,gmunu,gUmunu)
implicit none
integer,parameter::n=4
integer,intent(in)::k
real*8::Tmunu(n,n,0:k-1,0:k-1)
real*8::umu(n,0:k-1,0:k-1)
real*8::en(0:k-1,0:k-1)
real*8::gmunu(n,n),gUmunu(n,n)
real*8::tau_ampt
integer,parameter :: lda=n,ldvl=n,ldvr=n
integer,parameter :: lwmax=1000
integer :: info,lwork,i,j,o,xi,yi
real(kind=8) :: rwork(2*n),aa
complex(kind=8) :: b(lda,n),vl(ldvl,n),vr(ldvr,n),w(n),work(lwmax)
external::zgeev



do i=1,n
do j=1,n
aa=0.0d0
do o=1,n
aa=aa+Tmunu(i,o,xi,yi)*gmunu(o,j)
end do
b(i,j)=cmplx(aa,0.d0)
end do
!print*,(b(i,j),j=1,n)
end do
lwork=-1
call zgeev('vectors','vectors',n,b,lda,w,vl,ldvl,vr,ldvr,work,lwork,rwork,info)
lwork=min(lwmax,int(work(1)))
call zgeev('vectors','vectors',n,b,lda,w,vl,ldvl,vr,ldvr,work,lwork,rwork,info)

if (info .eq. 0) then
call choose(xi,yi,w,vr,k,tau_ampt,Tmunu,umu,en,gmunu,gUmunu)
endif

end subroutine eigenvalue

subroutine choose(xi,yi,eigen,eigenvector,k,tau_ampt,Tmunu,umu,en,gmunu,gUmunu)
integer,parameter::n=4
integer,intent(in)::k
integer,parameter :: lda=n,ldvl=n,ldvr=n
real*8::Tmunu(n,n,0:k-1,0:k-1)
real*8::umu(n,0:k-1,0:k-1)
real*8::en(0:k-1,0:k-1)
integer::xi,yi,i,j
complex(kind=8)::eigenvector(ldvr,n),eigen(n)
real*8::egim,vg,egr,vg1,vg2,vg3,vg4,vg12
real*8::gmunu(n,n),gUmunu(n,n)
real*8::tau_ampt


do i=1,n
!       write(*,*)eigen(i)
	egim=aimag(eigen(i))
        egr=real(eigen(i))
        if(abs(egim).gt.0.0d0)then
         cycle
         else
        vg=sqrt((real(eigenvector(1,i)))**2-(real(eigenvector(2,i)))**2-(real(eigenvector(3,i)))**2- &
         (tau_ampt**2)*(real(eigenvector(4,i)))**2)
        vg12=sqrt((real(eigenvector(1,i)))**2+(real(eigenvector(2,i)))**2+(real(eigenvector(3,i)))**2+(real(eigenvector(4,i)))**2)
   
	if((vg.gt.0d0).and.(egr.gt.0.0d0))then
        vg1=aimag(eigenvector(1,i))
        vg2=aimag(eigenvector(2,i))
        vg3=aimag(eigenvector(3,i))
        vg4=aimag(eigenvector(4,i))
        if((abs(vg1-0.0d0).lt.00001).and.(abs(vg2-0.0d0).lt.00001).and.(abs(vg3-0.0d0).lt.00001).and.(abs(vg4-0.0d0).lt.00001))then
	umu(1,xi,yi)=(1./vg)*real(eigenvector(1,i)) *(1./vg12)
	umu(2,xi,yi)=(1./vg)*real(eigenvector(2,i)) *(1./vg12)
	umu(3,xi,yi)=(1./vg)*real(eigenvector(3,i)) *(1./vg12)
	umu(4,xi,yi)=(1./vg)*real(eigenvector(4,i)) *(1./vg12)
	en(xi,yi)=real(eigen(i))     
        else
        umu(1,xi,yi)=1.0d0
        umu(2,xi,yi)=0.0d0
        umu(3,xi,yi)=0.0d0
        umu(4,xi,yi)=0.0d0  
	en(xi,yi)=real(eigen(i))           
        end if
        end if
        end if
end do
end subroutine choose
