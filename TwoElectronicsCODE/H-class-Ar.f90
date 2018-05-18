
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
!----------外场---------------------

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
!fwhm=7.d0*(10.d0**(-15.d0))                         !!! 半高宽 秒 (激光强度)；
!tao=fwhm/2.d0/( sqrt(2.d0*log(sqrt(2.d0)))  )       !!! tao 秒； f(t)=E0 * exp[ -1/2 (t/tao)^2 ]
!tao=tao/(2.4189d0*10.d0**(-17.d0)) 
!tt=15.d0*(10.d0**(-15.d0))                          !!! 秒
!Tp=tt/(2.4189d0*(10.d0**(-17.d0)))                  !!! 原子单位
!h=2.d0*Tp/(n-1)
!n2=0.77d0                                           !椭偏率
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

!open(file='temp1.dat',unit=2)

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



!----------------------------------



j=1

do while(j.le.m)

 CALL RANDOM_NUMBER(aa2)
 CALL RANDOM_NUMBER(aa3)

 aa2=aa2*pai
 aa3=aa3*2*pai

 x1(1)=(rr)*sin(aa2)*cos(aa3)
 y1(1)=(rr)*sin(aa2)*sin(aa3)
 z1(1)=(rr)*cos(aa2)
   
 x2(1)=-x1(1)
 y2(1)=-y1(1)
 z2(1)=-z1(1)

!--------------------------------------
 CALL RANDOM_NUMBER(aa2)
 CALL RANDOM_NUMBER(aa3)
 CALL RANDOM_NUMBER(aa4)
 CALL RANDOM_NUMBER(aa5)

 aa2=aa2*pai
 aa4=aa4*pai
 aa3=aa3*2*pai
 aa5=aa5*2*pai

 Px1(1)=(pp)*sin(aa2)*cos(aa3)
 Py1(1)=(pp)*sin(aa2)*sin(aa3)
 pz1(1)=(pp)*cos(aa2)
   
 Px2(1)=(pp)*sin(aa4)*cos(aa5)
 Py2(1)=(pp)*sin(aa4)*sin(aa5)
 pz2(1)=(pp)*cos(aa4)
!----------------------------------------------------------
do i=1,nn

U(:,1)=(/x1(1),Px1(1),x2(1),Px2(1),y1(1),py1(1),y2(1),py2(1),z1(1),pz1(1),z2(1),pz2(1)/)!!!CHUfHI


y11=x1(1)
y12=Px1(1)
y13=x2(1)
y14=Px2(1)
y15=y1(1)
y16=Py1(1)
y17=y2(1)
y18=Py2(1)
y111=z1(1)
y112=pz1(1)
y113=z2(1)
y114=pz2(1)


t1=t(1)

K1=(/g1(y11,y12,y13,y14,y15,y16,y17,y18,y111,y112,y113,y114,t1,a,q),&
&f1(y11,y12,y13,y14,y15,y16,y17,y18,y111,y112,y113,y114,t1,a,q),&
&g2(y11,y12,y13,y14,y15,y16,y17,y18,y111,y112,y113,y114,t1,a,q),&
&f2(y11,y12,y13,y14,y15,y16,y17,y18,y111,y112,y113,y114,t1,a,q),&
&g3(y11,y12,y13,y14,y15,y16,y17,y18,y111,y112,y113,y114,t1,a,q),&
&f3(y11,y12,y13,y14,y15,y16,y17,y18,y111,y112,y113,y114,t1,a,q),&
&g4(y11,y12,y13,y14,y15,y16,y17,y18,y111,y112,y113,y114,t1,a,q),&
&f4(y11,y12,y13,y14,y15,y16,y17,y18,y111,y112,y113,y114,t1,a,q),&
&g5(y11,y12,y13,y14,y15,y16,y17,y18,y111,y112,y113,y114,t1,a,q),&
&f5(y11,y12,y13,y14,y15,y16,y17,y18,y111,y112,y113,y114,t1,a,q),&
&g6(y11,y12,y13,y14,y15,y16,y17,y18,y111,y112,y113,y114,t1,a,q),&
&f6(y11,y12,y13,y14,y15,y16,y17,y18,y111,y112,y113,y114,t1,a,q)/)


y21=x1(1)+h*K1(1)/2
y22=px1(1)+h*K1(2)/2
y23=x2(1)+h*K1(3)/2
y24=px2(1)+h*K1(4)/2
y25=y1(1)+h*K1(5)/2
y26=Py1(1)+h*K1(6)/2
y27=y2(1)+h*K1(7)/2
y28=Py2(1)+h*K1(8)/2
y211=z1(1)+h*K1(9)/2
y212=pz1(1)+h*K1(10)/2
y213=z2(1)+h*K1(11)/2
y214=pz2(1)+h*K1(12)/2

t2=t(1)+h/2


K2=(/g1(y21,y22,y23,y24,y25,y26,y27,y28,y211,y212,y213,y214,t2,a,q),&
&f1(y21,y22,y23,y24,y25,y26,y27,y28,y211,y212,y213,y214,t2,a,q),&
&g2(y21,y22,y23,y24,y25,y26,y27,y28,y211,y212,y213,y214,t2,a,q),&
&f2(y21,y22,y23,y24,y25,y26,y27,y28,y211,y212,y213,y214,t2,a,q),&
&g3(y21,y22,y23,y24,y25,y26,y27,y28,y211,y212,y213,y214,t2,a,q),&
&f3(y21,y22,y23,y24,y25,y26,y27,y28,y211,y212,y213,y214,t2,a,q),&
&g4(y21,y22,y23,y24,y25,y26,y27,y28,y211,y212,y213,y214,t2,a,q),&
&f4(y21,y22,y23,y24,y25,y26,y27,y28,y211,y212,y213,y214,t2,a,q),&
&g5(y21,y22,y23,y24,y25,y26,y27,y28,y211,y212,y213,y214,t2,a,q),&
&f5(y21,y22,y23,y24,y25,y26,y27,y28,y211,y212,y213,y214,t2,a,q),&
&g6(y21,y22,y23,y24,y25,y26,y27,y28,y211,y212,y213,y214,t2,a,q),&
&f6(y21,y22,y23,y24,y25,y26,y27,y28,y211,y212,y213,y214,t2,a,q)/)

y31=x1(1)+h*K2(1)/2
y32=Px1(1)+h*K2(2)/2
y33=x2(1)+h*K2(3)/2
y34=Px2(1)+h*K2(4)/2
y35=y1(1)+h*K2(5)/2
y36=Py1(1)+h*K2(6)/2
y37=y2(1)+h*K2(7)/2
y38=Py2(1)+h*K2(8)/2
y311=z1(1)+h*K2(9)/2
y312=Pz1(1)+h*K2(10)/2
y313=z2(1)+h*K2(11)/2
y314=Pz2(1)+h*K2(12)/2



t3=t(1)+h/2

K3=(/g1(y31,y32,y33,y34,y35,y36,y37,y38,y311,y312,y313,y314,t3,a,q),&
&f1(y31,y32,y33,y34,y35,y36,y37,y38,y311,y312,y313,y314,t3,a,q),&
&g2(y31,y32,y33,y34,y35,y36,y37,y38,y311,y312,y313,y314,t3,a,q),&
&f2(y31,y32,y33,y34,y35,y36,y37,y38,y311,y312,y313,y314,t3,a,q),&
&g3(y31,y32,y33,y34,y35,y36,y37,y38,y311,y312,y313,y314,t3,a,q),&
&f3(y31,y32,y33,y34,y35,y36,y37,y38,y311,y312,y313,y314,t3,a,q),&
&g4(y31,y32,y33,y34,y35,y36,y37,y38,y311,y312,y313,y314,t3,a,q),&
&f4(y31,y32,y33,y34,y35,y36,y37,y38,y311,y312,y313,y314,t3,a,q),&
&g5(y31,y32,y33,y34,y35,y36,y37,y38,y311,y312,y313,y314,t3,a,q),&
&f5(y31,y32,y33,y34,y35,y36,y37,y38,y311,y312,y313,y314,t3,a,q),&
&g6(y31,y32,y33,y34,y35,y36,y37,y38,y311,y312,y313,y314,t3,a,q),&
&f6(y31,y32,y33,y34,y35,y36,y37,y38,y311,y312,y313,y314,t3,a,q)/)

y41=x1(1)+h*K3(1)
y42=Px1(1)+h*K3(2)
y43=x2(1)+h*K3(3)
y44=Px2(1)+h*K3(4)
y45=y1(1)+h*K3(5)
y46=Py1(1)+h*K3(6)
y47=y2(1)+h*K3(7)
y48=Py2(1)+h*K3(8)
y411=z1(1)+h*K3(9)
y412=Pz1(1)+h*K3(10)
y413=z2(1)+h*K3(11)
y414=Pz2(1)+h*K3(12)

t4=t(1)+h

K4=(/g1(y41,y42,y43,y44,y45,y46,y47,y48,y411,y412,y413,y414,t4,a,q),&
&f1(y41,y42,y43,y44,y45,y46,y47,y48,y411,y412,y413,y414,t4,a,q),&
&g2(y41,y42,y43,y44,y45,y46,y47,y48,y411,y412,y413,y414,t4,a,q),&
&f2(y41,y42,y43,y44,y45,y46,y47,y48,y411,y412,y413,y414,t4,a,q),&
&g3(y41,y42,y43,y44,y45,y46,y47,y48,y411,y412,y413,y414,t4,a,q),&
&f3(y41,y42,y43,y44,y45,y46,y47,y48,y411,y412,y413,y414,t4,a,q),&
&g4(y41,y42,y43,y44,y45,y46,y47,y48,y411,y412,y413,y414,t4,a,q),&
&f4(y41,y42,y43,y44,y45,y46,y47,y48,y411,y412,y413,y414,t4,a,q),&
&g5(y41,y42,y43,y44,y45,y46,y47,y48,y411,y412,y413,y414,t4,a,q),&
&f5(y41,y42,y43,y44,y45,y46,y47,y48,y411,y412,y413,y414,t4,a,q),&
&g6(y41,y42,y43,y44,y45,y46,y47,y48,y411,y412,y413,y414,t4,a,q),&
&f6(y41,y42,y43,y44,y45,y46,y47,y48,y411,y412,y413,y414,t4,a,q)/)

do ii=1,12
U(ii,2)=U(ii,1)+(h*(K1(ii)+2.d0*K2(ii)+2.d0*K3(ii)+k4(ii))/6.d0)
end do

x1(2)=U(1,2)
Px1(2)=U(2,2)
x2(2)=U(3,2)
Px2(2)=U(4,2)
y1(2)=U(5,2)
Py1(2)=U(6,2)
y2(2)=U(7,2)
Py2(2)=U(8,2)
z1(2)=U(9,2)
Pz1(2)=U(10,2)
z2(2)=U(11,2)
Pz2(2)=U(12,2)
 
t(1)=t(1)+h

x1(1)=x1(2)
px1(1)=px1(2)
x2(1)=x2(2)
px2(1)=px2(2)          !!重新赋值   递推
y1(1)=y1(2)
py1(1)=py1(2)
y2(1)=y2(2)
py2(1)=py2(2)
z1(1)=z1(2)
pz1(1)=pz1(2)
z2(1)=z2(2)
pz2(1)=pz2(2)


  END DO


p1=px1(2)**2.d0+py1(2)**2.d0+pz1(2)**2.d0
r1=sqrt(x1(2)**2.d0+y1(2)**2.d0+z1(2)**2.d0)
p2=px2(2)**2.d0+py2(2)**2.d0+pz2(2)**2.d0
r2=sqrt(x2(2)**2.d0+y2(2)**2.d0+z2(2)**2.d0)
r12=sqrt((X1(2)-X2(2))**2.D0+(Y1(2)-Y2(2))**2.D0+(Z1(2)-Z2(2))**2.D0)
    
ef=p1/2.d0-2.d0/r1+1.d0/r1**2.d0*q**2.d0/4.d0/a*dexp(a*(1.d0-(r1**2.d0*p1/q**2.d0)**2.d0))&
 &+p2/2.d0-2.d0/r2+1.d0/r2**2.d0*q**2.d0/4.d0/a*dexp(a*(1.d0-(r2**2.d0*p2/q**2.d0)**2.d0))&
 &+1.d0/r12

if(ef.le.-1.587d0)then

WRITE(15,'(1X,21(1X,F15.7))') x1(2),Px1(2),x2(2),Px2(2),y1(2),py1(2),y2(2),py2(2),z1(2),pz1(2),z2(2),pz2(2) !!!CHUfHI

rx1(j)=x1(2)
rx2(j)=x2(2)
ry1(j)=y1(2)
ry2(j)=y2(2)
rz1(j)=z1(2)
rz2(j)=z2(2)

rpx1(j)=px1(2)
rpx2(j)=px2(2)
rpy1(j)=py1(2)
rpy2(j)=py2(2)
rpz1(j)=pz1(2)
rpz2(j)=pz2(2)


j=j+1
else
j=j
end if
                                                                                          !j循环结束
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------函数---------------------------

FUNCTION f1(x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,t,a,q)   
Implicit None                                                  !!龙格库塔方法中用到的函数
REAL*8 f1,x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,t,a,q


f1=x1*  ( (q**2.d0/2.d0/a/(z1**2.d0+x1**2.d0+y1**2.d0)**2.d0+(pz1**2.d0+px1**2.d0+py1**2.d0)**2.d0/q**2.d0)&
   &*exp(a*(1.d0-((z1**2.d0+y1**2.d0+x1**2.d0)*(pz1**2.d0+px1**2.d0+py1**2.d0)/q**2.d0)**2.d0))-2.d0/sqrt((z1**2.d0+y1**2.d0+x1**2.d0)**3.d0) )&
   &+(x1-x2)/sqrt(((z1-z2)**2.d0+(y1-y2)**2.d0+(x1-x2)**2.d0)**3.d0)




return
END


FUNCTION f2(x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,t,a,q) 
Implicit None                  
REAL*8 f2,x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,t,a,q




f2=x2*  ( (q**2.d0/2.d0/a/(z2**2.d0+x2**2.d0+y2**2.d0)**2.d0+(pz2**2.d0+px2**2.d0+py2**2.d0)**2.d0/q**2.d0)&
   &*exp(a*(1.d0-((z2**2.d0+y2**2.d0+x2**2.d0)*(pz2**2.d0+px2**2.d0+py2**2.d0)/q**2.d0)**2.d0))-2.d0/sqrt((z2**2.d0+y2**2.d0+x2**2.d0)**3.d0) )&
    &-(x1-x2)/sqrt(((z1-z2)**2.d0+(y1-y2)**2.d0+(x1-x2)**2.d0)**3.d0)



return
END

FUNCTION f3(x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,t,a,q)   
Implicit None 
REAL*8 f3,x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,t,a,q


f3=y1* ( (q**2.d0/2.d0/a/(z1**2.d0+x1**2.d0+y1**2.d0)**2.d0+(pz1**2.d0+px1**2.d0+py1**2.d0)**2.d0/q**2.d0)&
   &*exp(a*(1.d0-((z1**2.d0+y1**2.d0+x1**2.d0)*(pz1**2.d0+px1**2.d0+py1**2.d0)/q**2.d0)**2.d0))-2.d0/sqrt((z1**2.d0+y1**2.d0+x1**2.d0)**3.d0) )&
    &+(y1-y2)/sqrt(((z1-z2)**2.d0+(y1-y2)**2.d0+(x1-x2)**2.d0)**3.d0)




return
END



FUNCTION f4(x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,t,a,q) 
Implicit None 
REAL*8 f4,x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,t,a,q

f4=y2* ( (q**2.d0/2.d0/a/(z2**2.d0+x2**2.d0+y2**2.d0)**2.d0+(pz2**2.d0+px2**2.d0+py2**2.d0)**2.d0/q**2.d0)&
    &*exp(a*(1.d0-((z2**2.d0+y2**2.d0+x2**2.d0)*(pz2**2.d0+px2**2.d0+py2**2.d0)/q**2.d0)**2.d0))-2.d0/sqrt((z2**2.d0+y2**2.d0+x2**2.d0)**3.d0) )&
    &-(y1-y2)/sqrt(((z1-z2)**2.d0+(y1-y2)**2.d0+(x1-x2)**2.d0)**3.d0)


return
end


FUNCTION f5(x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,t,a,q) 
Implicit None  
REAL*8 f5,x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,t,a,q


f5=z1* ( (q**2.d0/2.d0/a/(z1**2.d0+x1**2.d0+y1**2.d0)**2.d0+(pz1**2.d0+px1**2.d0+py1**2.d0)**2.d0/q**2.d0)&
    &*exp(a*(1.d0-((z1**2.d0+y1**2.d0+x1**2.d0)*(pz1**2.d0+px1**2.d0+py1**2.d0)/q**2.d0)**2.d0))-2.d0/sqrt((z1**2.d0+y1**2.d0+x1**2.d0)**3.d0) )&
    &+(z1-z2)/sqrt(((z1-z2)**2.d0+(y1-y2)**2.d0+(x1-x2)**2.d0)**3.d0)



return
END



FUNCTION f6(x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,t,a,q) 
Implicit None 
REAL*8 f6,x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,t,a,q


f6=z2* ( (q**2.d0/2.d0/a/(z2**2.d0+x2**2.d0+y2**2.d0)**2.d0+(pz2**2.d0+px2**2.d0+py2**2.d0)**2.d0/q**2.d0)&
    &*dexp(a*(1.d0-((z2**2.d0+y2**2.d0+x2**2.d0)*(pz2**2.d0+px2**2.d0+py2**2.d0)/q**2.d0)**2.d0))-2.d0/sqrt((z2**2.d0+y2**2.d0+x2**2.d0)**3.d0) )&
    &-(z1-z2)/sqrt(((z1-z2)**2.d0+(y1-y2)**2.d0+(x1-x2)**2.d0)**3.d0)
!---------------


return
end




FUNCTION  g1(x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,t,a,q) 
Implicit None 
REAL*8 t,x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,g1,a,q

g1=px1*( 1.d0 - 1.d0/q**2.D0*(z1**2.D0+y1**2.d0+x1**2.d0)*(pz1**2.d0+px1**2.d0+py1**2.d0)&
     &*dexp(a*(1.d0-((z1**2.d0+x1**2.d0+y1**2.d0)*(pz1**2.d0+px1**2.d0+py1**2.d0)/q**2.d0)**2.d0)) )
  

return

end


FUNCTION  g2(x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,t,a,q) 
Implicit None 
REAL*8 t,x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,g2,a,q

g2=px2*( 1.d0 - 1.d0/q**2.D0*(z2**2.D0+y2**2.d0+x2**2.d0)*(pz2**2.d0+px2**2.d0+py2**2.d0)&
    &*dexp(a*(1.d0-((z2**2.d0+x2**2.d0+y2**2.d0)*(pz2**2.d0+px2**2.d0+py2**2.d0)/q**2.d0)**2.d0)) )


return

end



FUNCTION  g3(x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,t,a,q) 
Implicit None 
REAL*8 t,x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,g3,a,q

g3=py1*( 1.d0 - 1.d0/q**2.D0*(z1**2.D0+y1**2.d0+x1**2.d0)*(pz1**2.d0+px1**2.d0+py1**2.d0)&
      &*dexp(a*(1.d0-((z1**2.d0+x1**2.d0+y1**2.d0)*(pz1**2.d0+px1**2.d0+py1**2.d0)/q**2.d0)**2.d0)) )

return
end


FUNCTION  g4(x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,t,a,q) 
Implicit None 
REAL*8 t,x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,g4,a,q

g4=py2*( 1.d0 - 1.d0/q**2.D0*(z2**2.D0+y2**2.d0+x2**2.d0)*(pz2**2.d0+px2**2.d0+py2**2.d0)&
      &*dexp(a*(1.d0-((z2**2.d0+x2**2.d0+y2**2.d0)*(pz2**2.d0+px2**2.d0+py2**2.d0)/q**2.d0)**2.d0)) )

return

end


FUNCTION  g5(x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,t,a,q) 
Implicit None 
REAL*8 t,x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,g5,a,q

g5=pz1*( 1.d0 - 1.d0/q**2.D0*(z1**2.D0+y1**2.d0+x1**2.d0)*(pz1**2.d0+px1**2.d0+py1**2.d0)&
    &*dexp(a*(1.d0-((z1**2.d0+x1**2.d0+y1**2.d0)*(pz1**2.d0+px1**2.d0+py1**2.d0)/q**2.d0)**2.d0)) )

return

end


FUNCTION  g6(x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,t,a,q) 
Implicit None 
REAL*8 t,x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,g6,a,q

g6=pz2*( 1.d0 - 1.d0/q**2.D0*(z2**2.D0+y2**2.d0+x2**2.d0)*(pz2**2.d0+px2**2.d0+py2**2.d0)&
    &*dexp(a*(1.d0-((z2**2.d0+x2**2.d0+y2**2.d0)*(pz2**2.d0+px2**2.d0+py2**2.d0)/q**2.d0)**2.d0)) )

return
end 

!-------------外场函数------------------------------------------------------



Function  E1(t,t0,ee0,w,n1,n2,tao,tp)
Implicit None 
REAL*8 ee0,w,t,n1,n2,t0,E1,QQ,h,fwhm,n,tao,tp


!-------------T------------------------                                              !原子单位
!if(t.lt.n1*T0) then
!Q=t/(n1*T0)
!end if
!if((t.ge.n1*T0).and.(t.le.(n1+n2)*T0)) then
!Q=1.d0
!endif
!if((t.gt.(n1+n2)*T0).and.(t.le.(2.d0*n1+n2)*T0)) then
!Q=-t/(n1*T0)+(2.d0*n1+n2)/n1
!endif
!if(t.GT.(2.d0*n1+n2)*T0) then
!Q=0.D0
!endif
!E1=ee0*Q*cos(w*t)
!-------------ellipse---------------


E1=ee0*exp( -((t-tp)/tao)**2.d0/2.d0 )*(n2/sqrt(1.d0+n2**2.d0))*SIN(w*t+N1)


end

Function  E2(t,t0,ee0,w,n1,n2,tao,tp)
Implicit None 
REAL*8 ee0,w,t,n1,n2,t0,E2,QQ,h,fwhm,n,tao,tp
  

!-------------T------------------------                                            
!if(t.lt.n1*T0) then
!Q=t/(n1*T0)
!end if
!if((t.ge.n1*T0).and.(t.le.(n1+n2)*T0)) then
!Q=1.d0
!endif
!if((t.gt.(n1+n2)*T0).and.(t.le.(2.d0*n1+n2)*T0)) then
!Q=-t/(n1*T0)+(2.d0*n1+n2)/n1
!endif
!if(t.GT.(2.d0*n1+n2)*T0) then
!Q=0.D0
!endif 

!E2=ee0*Q*sin(w*t)
!-------------ellipse---------------

E2=ee0*exp( -((t-tp)/tao)**2.d0/2.d0 ) *(1.d0/sqrt(1.d0+n2**2.d0))*cos(w*t+N1)

end

!--------------------------------





