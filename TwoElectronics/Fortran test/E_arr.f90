
Program classical

Implicit None 
integer , parameter :: nx=351, ny=301,nn=10000,n=40000,m=10000
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
&y41,y42,y43,y44,y45,y46,y47,y48,y49,y411,y412,y413,y414,E(n)
REAL*8 ee0,n1,n2,h,w,t0,fwhm,qq1,tdi,tre,trec,aaa1,za,zb,tao,tp,tt,c,lamda,W1,W2,zz
Real*8 matrix(nx, ny), x(nx), y(ny), min,z,pd,ip2,s,pion
Real*8 p1, p2, r1, r2,r12,ef
REAL*8,EXTERNAL::g1,g2,g3,g4,g5,g6,f1,f2,f3,f4,f5,f6,E1,E2,QQ



!open(file='interaction.DAT',unit=20)
open(file='E.DAT',unit=21)

open(file='E2.DAT',unit=22)

!------------Ar--------------------
 pai=3.1415926535897932384626433832795d0 
 a=2.d0
 q=1.225d0
 E0=-1.59d0
 Ip2=-1.065d0
!----------外场---------------------

ee0=sqrt(1d15/3.51d16)
!w=0.057d0
W1=0.0584d0
W2=0.117d0
t0=2*pai/w1
n1=2.d0
n2=6.d0
h=(2.d0*n1+n2)*t0/dble(n)


tp=0.8d0
!-------------------------------
!c=2.997925d0*(10.d0**8.d0)
!lamda=740.d0*(10.d0**(-9.d0)) 
!w1=2.d0*pai*c/lamda         
!w1=w1*2.4189d0*(10.d0**(-17.d0))
!t0=2.d0*pai/w1
!fwhm=7.d0*(10.d0**(-15.d0))                         !!! 半高宽 秒 (激光强度)；
!tao=fwhm/2.d0/( sqrt(2.d0*log(sqrt(2.d0)))  )       !!! tao 秒； f(t)=E0 * exp[ -1/2 (t/tao)^2 ]
!tao=tao/(2.4189d0*10.d0**(-17.d0)) 
!tt=15.d0*(10.d0**(-15.d0))                          !!! 秒
!Tp=tt/(2.4189d0*(10.d0**(-17.d0)))                  !!! 原子单位
!h=2.d0*Tp/(n-1)
!n2=0.7d0                                           !椭偏率
!EE0=2.742*10**3.d0*sqrt(4d15)    
!EE0=EE0/(5.1421*(10**11.D0))
!n1=0.d0*pai
!------------------------------------
do i=1,N

t1=h*i
E(i)=sqrt(E1(t1,t0,ee0,w1,w2,n1,n2,tao,tp)**2.d0+E2(t1,t0,ee0,w1,w2,n1,n2,tao,tp)**2.d0)
write(21,'(1X,21(1X,F15.7))')t1/t0,E1(t1,t0,ee0,w1,w2,n1,n2,tao,tp),E2(t1,t0,ee0,w1,w2,n1,n2,tao,tp),E(i)

enddo
    
do i=1,2*N

t1=h*i * 0.5
E(i)=sqrt(E1(t1,t0,ee0,w1,w2,n1,n2,tao,tp)**2.d0+E2(t1,t0,ee0,w1,w2,n1,n2,tao,tp)**2.d0)
write(21,'(1X,21(1X,F15.7))')t1/t0,E1(t1,t0,ee0,w1,w2,n1,n2,tao,tp),E2(t1,t0,ee0,w1,w2,n1,n2,tao,tp),E(i)

enddo
    
end








Function  E1(t,t0,ee0,w1,w2,n1,n2,tao,tp)
Implicit None 
REAL*8 ee0,w,t,n1,n2,t0,E1,QQ,h,fwhm,n,tao,tp,w1,w2


!-------------T------------------------                                              !原子单位
!if(t.lt.n1*T0) then
!QQ=t/(n1*T0)
!end if
!if((t.ge.n1*T0).and.(t.le.(n1+n2)*T0)) then
!QQ=1.d0
!endif
!if((t.gt.(n1+n2)*T0).and.(t.le.(2.d0*n1+n2)*T0)) then
!QQ=-t/(n1*T0)+(2.d0*n1+n2)/n1
!endif
!if(t.GT.(2.d0*n1+n2)*T0) then
!QQ=0.D0
!endif
!----------sin2-------------------------
QQ=(DSIN(w1/2.D0/(2*n1+n2)*t))**2.D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!E1=EE0*QQ*(sin(w1*t)-sin(w2*t))


E1=(ee0/(1.d0+tp))*QQ*(sin(w1*t+tao))-(ee0*tp/(1.d0+tp))*QQ*(sin(w2*t+2*tao))
!-------------ellipse---------------


!E1=ee0*exp( -((t-tp)/tao)**2.d0/2.d0 )*(n2/sqrt(1.d0+n2**2.d0))*SIN(w1*t+N1)
!-------------------------------------



!E1=ee0*QQ*sin(w1*t)
return
end

Function  E2(t,t0,ee0,w1,w2,n1,n2,tao,tp)
Implicit None 
REAL*8 ee0,w,t,n1,n2,t0,E2,QQ,h,fwhm,n,tao,tp,w1,w2
  

!-------------T------------------------                                            
!if(t.lt.n1*T0) then
!QQ=t/(n1*T0)
!end if
!if((t.ge.n1*T0).and.(t.le.(n1+n2)*T0)) then
!QQ=1.d0
!endif
!if((t.gt.(n1+n2)*T0).and.(t.le.(2.d0*n1+n2)*T0)) then
!QQ=-t/(n1*T0)+(2.d0*n1+n2)/n1
!endif
!if(t.GT.(2.d0*n1+n2)*T0) then
!QQ=0.D0
!endif 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!E2=EE0*QQ*(COS(w1*t)+COS(w2*t))

!----------sin2-------------------------
QQ=(DSIN(w1/2.D0/(2*n1+n2)*t))**2.D0


E2=(ee0/(1.d0+tp))*QQ*(cos(w1*t+tao))+(ee0*tp/(1.d0+tp))*QQ*(cos(w2*t+2*tao))

!-------------ellipse---------------

!E2=ee0*exp( -((t-tp)/tao)**2.d0/2.d0 ) *(1.d0/sqrt(1.d0+n2**2.d0))*cos(w1*t+N1)

!------------------------------------------

!E2=ee0*QQ*cos(w1*t)
return
end

!--------------------------------


