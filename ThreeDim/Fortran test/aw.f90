
Program classical

Implicit None 
integer , parameter :: n=40000,m=100000,mm=10000
real*8 :: K1(12),K2(12),K3(12),K4(12),U(12,2),t(2),x1(2,m),y1(2,m),px1(2,m),&
&py1(2,m),x2(2,m),y2(2,m),px2(2,m),py2(2,m),z1(2,m),pz1(2,m),z2(2,m),pz2(2,m)
real*8 :: h,a,A1,E0,t1,t2,t3,t4,y11,y12,y13,y14,y15,y16,y17,y18,y111,y112,y113,y114,y21,y22,y23,y24,y25,y26,y27,y28,y211,y212,y213,y214,AW(2*n),&
&y31,y32,y33,y34,y35,y36,y37,y38,y311,y312,y313,y314,y41,y42,y43,y44,y45,y46,y47,y48,y411,y412,y413,y414,EE1(n),EE2(n),rr(n),r2(n),aa11(m),&
&pyy2(2,m),pyy1(2,m),pxx1(2,m),pxx2(2,m),pzz1(2,m),pzz2(2,m),xx1(2,m),xx2(2,m),yy1(2,m),yy2(2,m),zz1(2,m),zz2(2,m),xxa(n),xxb(n),pxxa(n),pxxb(n)&
&,yya(n),yyb(n),pyya(n),pyyb(n),zza(n),zzb(n),pzza(n),pzzb(n),aa1,aa2,aa3,aa4,aa5,aa6,kk,k11,k22,m1,m2,aaa1,bbb1,tion1,tion2,tion3,tion4,tdi1,tdi2,tdi3,tdi4,tion13,tion24,tdi13,tdi24
REAL*8 pai,w,n1,n2,n3,n4,ee0,t0,EE,xx22,Ekall,p,pion,lx1,lx2,ly1,ly2,lz1,lz2,xa,xb,l1,l2,z,zz,aa,bb,ccc1,tdi,tion,pxx,pyy,pzz,pzzz,a11,v,s1,s2,ss,&
&tt1,tt2,tt3,tt4,ds1(2*n),za,zb,tt,ttt,ccc,pzzz1,pzzz2,tttt
INTEGER i,j,ii,jjj,R,nn,trec,flag,aaa
REAL*8,EXTERNAL::g1,g2,g3,g4,g5,g6,Q,E1,E2,f1,f2,f3,f4,f5,f6



w=0.057d0
pai=3.1415926535897932384626433832795d0 
t0=2*pai/w

n1=2.d0
n2=6.d0

h=(2.d0*n1+n2)*t0/dble(n)


tdi=40000.d0

kk=1.d0/2.d0

open(file='E(t).DAT',unit=14)

ee0=dsqrt((2.51d15)/(3.51d16))




do i=1,2*n

t1=0.5*h*i

AW(i)=(ee0/w)*(SIN((pai*t1)/(10*t0))**2.D0)*cos(w*t1)

end do

do i=1,2*n

call DS(0.5*h,AW,2*n,DS1)

WRITE(14,*)I*0.5*h/t0,Aw(i),-ds1(i)

end do

End Program classical

SUBROUTINE DS(H,HS,N,DS1)
REAL*8 H,HS(N),DS1(N)
integer N

DS1(1)=(HS(2)-HS(1))/H
DS1(N)=(HS(N)-HS(N-1))/H

DO I=2,N-1
DS1(I)=(HS(I+1)-HS(I-1))/2.d0/H
END DO

RETURN
END


