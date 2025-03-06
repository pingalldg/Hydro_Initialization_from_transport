!!Lapack library is needed!!sudo apt-get install .....
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!n= 4-dimensional space-time tau,x,y,etas
!k=number of grid points for initial condition for hydro
!dx=grid spacing
!Tmunu=energy momentum tensor
!umu=flow velocity components 
!en=energy density
!pr=pressure
!np= number of points in the sp95 data file
!pimunu=vicous term
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ampttohydro2(k,dx,tau_ampt,sigma,netas,netas_print,detas,etas_min,min_etas,max_etas,sigma_etas,ks_0)
integer::xcent,ycent,xbegin,ybegin,xend,yend,xi,yi,l
integer::num,ity,i,j,m,o,ietas,netas,netas_print
integer,intent(in)::k
integer,parameter::n=4,np=230
real*8::posx,posy,sigma,dx,E,etas,a,tau,etas_min,detas,sigma_etas,ks_0
real*8::xx,yy,zz,pxx,pyy,pzz,massxx,ft5,ftnd,taut,pressure,energy,etas_val,min_etas,max_etas
real*8::Tmunu(n,n,0:k-1,0:k-1),ppp(np),eee(np),y2(np),pimunu(n,n,0:k-1,0:k-1) 
real*8::umu(n,0:k-1,0:k-1)
real*8::en(0:k-1,0:k-1),pr(0:k-1,0:k-1)
real*8::p(n)
real*8::gmunu(n,n),gUmunu(n,n)!!!!!gmunu lower,gUmunu upper metric tensor!!!!!
integer::steps,count
real*8::tau_ampt



open(10,file="partontimets1.dat")
open(220,file="music_input.dat")
open(30,file="los2.dat") 
! open(10,file="test1.dat")
!write(220,*)netas, k, k, detas, dx, dx,tau_ampt
write(220,*)netas_print, k, k, detas, dx, dx,tau_ampt

read(10,*)a          
m=int(a)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! etas circle!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do ietas=1,netas
      etas_val=etas_min+(ietas-1)*detas
      write(*,*)etas_val
!!!!!!!!!!!!!!!!!!!!!!!!!!########initialize T^munu and u^mu!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
do i=1,n
do j=1,n
 do xi=0,k-1
   do yi=0,k-1 
      Tmunu(i,j,xi,yi)=0.0d0
      pimunu(i,j,xi,yi)=0.0d0
   end do
 end do
end do
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do i=1,n
 do xi=0,k-1
   do yi=0,k-1 
     if(i.ne.1)then
      umu(i,xi,yi)=0.0d0
      else
      umu(i,xi,yi)=1.0d0
      end if
   end do
 end do
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do xi=0,k-1
   do yi=0,k-1 
      en(xi,yi)=0.0d0
      pr(xi,yi)=0.0d0
   end do
 end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!metric tensor g_munu in tau eta_s frame!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 do i=1,n
   do j=1,n
   gmunu(i,j)=0.0d0
   gUmunu(i,j)=0.0d0
   end do
 end do
 gmunu(1,1)=1.00
 gmunu(2,2)=-1.00
 gmunu(3,3)=-1.00
 gmunu(4,4)=-(tau_ampt**2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 gUmunu(1,1)=1.00
 gUmunu(2,2)=-1.00
 gUmunu(3,3)=-1.00
 gUmunu(4,4)=-1.00/(tau_ampt**2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(*,*)"Initialization done"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!reading AMPT data!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
       do i=1,m
	read(10,*)num,ity,pxx,pyy,pzz,massxx,xx,yy,zz,ft5,ftnd,taut
!	read(10,*)num,ity,xx,yy,zz,pxx,pyy,pzz,massxx,ft5,ftnd,taut	
	if ( abs(ftnd-zz) .gt. 0.010 ) then 
	etas=0.5*log((ftnd+zz)/(ftnd-zz))
	if ((etas.ge.min_etas).and.(etas.lt.max_etas)) then
	E=sqrt(pxx**2+pyy**2+pzz**2+massxx**2)
	p(1)=E*cosh(etas)-pzz*sinh(etas)
	p(2)=pxx
	p(3)=pyy
	p(4)=(pzz*cosh(etas)-E*sinh(etas))*(1/tau_ampt)


	steps=80
	posx=xx
	posy=yy


	xcent=int(xx/dx)+int(k/2.)
	ycent=int(yy/dx)+int(k/2.)

	  if(xx.lt.0.)then
	   xcent=xcent-1
	  end if

		if(yy.lt.0.)then
		 ycent=ycent-1
		end if  


		xbegin=xcent-steps
		if(xbegin.lt.0)then
		 xbegin=0 
		end if

		ybegin=ycent-steps
		if(ybegin.lt.0)then
		 ybegin=0
		end if  


		xend=xcent+steps
		if(xend.gt.k)then
		 xend=k
		end if

	yend=ycent+steps
	if(yend.gt.k)then
	yend=k
	end if

		  
	do xi=xbegin,xend-1
	 do yi=ybegin,yend-1
                 do l=1,n
		 do o=l,4

           call distribute2(Tmunu,xi,yi,etas_val,etas,posx,posy,l,o,p,taut,k,dx,sigma,sigma_etas,ks_0)
			Tmunu(o,l,xi,yi)=Tmunu(l,o,xi,yi)
		 end do
		end do

	 end do
	end do
end if
end if
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Tmunu,en,umu write!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(*,*)"Tmunu ready"
do xi=0,k-1
 do yi=0,k-1             
 call eigenvalue(xi,yi,k,tau_ampt,Tmunu,umu,en,gmunu,gUmunu)   
! write(20,*)(xi-int(k/2.))*dx,(yi-int(k/2.))*dx,en(xi,yi),Tmunu(1,1,xi,yi),umu(1,xi,yi),umu(2,xi,yi),umu(3,xi,yi),umu(4,xi,yi)    
 end do
end do
write(*,*)"Energy done"
count=1

 !!!!!!!!!!!!! LOS eos!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!pressure finding!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i=1,np
read(30,*)energy,pressure
ppp(i)=pressure
eee(i)=energy
y2(i)=0.0d0
end do
call spline(eee,ppp,np,3.0000e+30,3.00000e+30,y2)
do xi=0,k-1
 do yi=0,k-1 
energy=en(xi,yi)
call splint(eee,ppp,y2,np,energy,pressure)
pr(xi,yi)=pressure
end do
end do
write(*,*)"Pressure done !!"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!pimunu detemination!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i=1,n
do j=1,n
 do xi=0,k-1
   do yi=0,k-1 
      pimunu(i,j,xi,yi)=Tmunu(i,j,xi,yi)-(en(xi,yi)+pr(xi,yi))*umu(i,xi,yi)*umu(j,xi,yi)+pr(xi,yi)*gUmunu(i,j)
   end do
 end do
end do
end do
write(*,*)"Pimunu ready !!"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!for MUSIC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

neta=40
deta=0.50
nx=k
ny=k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do xi=0,k-1
 do yi=0,k-1
 write(220,*)(xi-int(k/2.))*dx,(yi-int(k/2.))*dx,etas_val,en(xi,yi), umu(1,xi,yi),umu(2,xi,yi),&
                      umu(3,xi,yi),umu(4,xi,yi) 
end do
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!one etas circle complete!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
rewind(10)
rewind(30)
end do
 close(220)
write(*,*)"Hydro Initialization complete !!"
end subroutine ampttohydro2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!END!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!END!!!!!!END!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!END!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine distribute2(Tmunu,xi,yi,etas_val,etas,posx,posy,l,m,p,taut,k,dx,sigma,sigma_etas,ks_0)
integer::xi,yi,l,m
integer,intent(in)::k
integer,parameter::n=4
real*8::x,y,posx,posy,taut,dx,pi=2*asin(1.d0),sigma,sigma_etas,ks_0
real*8::Tmunu(n,n,0:k-1,0:k-1)
real*8::etas_val,etas
real*8::p(n)
f(x,y)=ks_0/(2.*pi*(sigma**2)*taut*(2*pi*(sigma_etas**2))**0.50)*&
           (p(l)*p(m)/p(1))*exp(-(x**2+y**2)/(2*sigma**2)-(etas_val-etas)**2/(2*sigma_etas**2))
x=(xi-int(k/2.))*dx
y=(yi-int(k/2.))*dx        

Tmunu(l,m,xi,yi)=Tmunu(l,m,xi,yi)+f(x-posx,y-posy)

end subroutine distribute2


