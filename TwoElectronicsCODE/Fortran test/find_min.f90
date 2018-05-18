
Program classical

Implicit None 
integer , parameter :: nx=351, ny=301,nn=10000,n=40000,m=100000
integer i,j, i1, j1,ii,jj
REAL*8 vh(nx,ny),vk(nx,ny),Ek(nx,ny),R(nx),P(ny)
REAL*8 x1(2),y1(2),z1(2),px1(2),py1(2),pz1(2),x2(2),y2(2),z2(2),px2(2),py2(2),pz2(2)&
&,U(12,2),k1(12),k2(12),k3(12),k4(12),EE1(n),EE2(n)&
&,rx1(m),ry1(m),rz1(m),rpx1(m),rpy1(m),rpz1(m),rx2(m),ry2(m),rz2(m),rpx2(m),rpy2(m),rpz2(m)&
&,xx1(2),yy1(2),zz1(2),pxx1(2),pyy1(2),pzz1(2),xx2(2),yy2(2),zz2(2),pxx2(2),pyy2(2),pzz2(2)&
&,xxa(n),yya(n),zza(n),pxxa(n),pyya(n),pzza(n),xxb(n),yyb(n),zzb(n),pxxb(n),pyyb(n),pzzb(n),rrr(n)

REAL*8 a,q,E0,Pai,RR,PP,aa2,aa3,aa4,aa5,t1,t2,t3,t4,t(1),aa,bb,tion,&
&y11,y12,y13,y14,y15,y16,y17,y18,y19,y111,y112,y113,y114,&
&y21,y22,y23,y24,y25,y26,y27,y28,y29,y211,y212,y213,y214,&
&y31,y32,y33,y34,y35,y36,y37,y38,y39,y311,y312,y313,y314,&
&y41,y42,y43,y44,y45,y46,y47,y48,y49,y411,y412,y413,y414
REAL*8 ee0,n1,n2,h,w,t0,fwhm,qq1,tdi,trec,aaa1,za,zb,tao,tp,tt,c,lamda
Real*8 matrix(nx, ny), x(nx), y(ny), min,z,pd,ip2
Real*8 p1, p2, r1, r2,r12,ef
REAL*8,EXTERNAL::g1,g2,g3,g4,g5,g6,f1,f2,f3,f4,f5,f6,E1,E2,QQ



open(file='initialization.DAT',unit=15)

!------------Ar--------------------
 pai=3.1415926535897932384626433832795d0 
 a=2.d0
 q=1.225d0
 E0=-1.59d0
 Ip2=-1.065d0
!----------�ⳡ---------------------

!ee0=sqrt(1d14/3.51d16)
!w=0.057d0
!t0=2*pai/w
!n1=2.d0
!n2=6.d0
!h=(2.d0*n1+n2)*t0/dble(n)
!-------------------------------
!c=2.997925d0*(10.d0**8.d0)
!lamda=740.d0*(10.d0**(-9.d0)) 
!w=2.d0*pai*c/lamda         
!w=w*2.4189d0*(10.d0**(-17.d0))
!t0=2.d0*pai/w
!fwhm=7.d0*(10.d0**(-15.d0))                         !!! ��߿� �� (����ǿ��)��
!tao=fwhm/2.d0/( sqrt(2.d0*log(sqrt(2.d0)))  )       !!! tao �룻 f(t)=E0 * exp[ -1/2 (t/tao)^2 ]
!tao=tao/(2.4189d0*10.d0**(-17.d0)) 
!tt=15.d0*(10.d0**(-15.d0))                          !!! ��
!Tp=tt/(2.4189d0*(10.d0**(-17.d0)))                  !!! ԭ�ӵ�λ
!h=2.d0*Tp/(n-1)
!n2=0.77d0                                           !��ƫ��
!EE0=2.742*10**3.d0*sqrt(4d15)    
!EE0=EE0/(5.1421*(10**11.D0))

!do i=1,N

!t1=h*i


!write(21,*)t1,E1(t1,t0,ee0,w,n1,n2,tao,tp),E2(t1,t0,ee0,w,n1,n2,tao,tp)

!enddo

!_______________________________________________________________________
!-----h2------------

do i=1, nx
R(i)=0.5d0+0.01d0*(i-1)
end do
do i=1, ny
P(i)=0.d0+0.01d0*(i-1)
end do

write(*,*) R(1:5)
write(*,*) P(1:5)

do j=1, ny
	do i=1, nx

    Vh(i,j)=(q**2.d0/(4.d0*a*r(i)**2.d0))*dexp(a*(1.d0-(r(i)*p(j)/q)**4.d0))
    Vk(i,j)=-2.d0/r(i)
    Ek(i,j)=(p(j)**2.d0)/2.d0

	matrix(i,j)=vh(i,j)+vk(i,j)+ek(i,j)+1.065d0	


	end do
	
end do

!min1!!!!!!!
call find_max(nx, ny, matrix, min, i1, j1)

rr=r(i1)
pp=p(j1)

write(*,*) min,i1,j1
write(*,*) rr,pp

end

!-----------------------------------------------------------
subroutine find_max(nx, ny, matrix, min, i, j)
Implicit None 
integer i, j,i1,j1,nx,ny
real*8 matrix(nx, ny), min,tem 


min=matrix(1,1)
i=1
j=1
do j1=1, ny
	do i1=1, nx
	if(min.ge.matrix(i1,j1)) then
	min=matrix(i1,j1)
	i=i1
	j=j1
	
	
	end if
	end do
end do

return
end