c        PROGRAM (09oct09) 
C
C        Jose M Albornoz
c        Grupo de Caos y Sistemas Complejos
C        Universidad de Los Andes 
C
C==========================================================
C	Un modelo basico de procesos de competencia
c	e inhibicion 
C==========================================================
c VARIABLES
c       y1 == S
c       y2 == P
c       y3 == I
c       y4 == f_L
c       y5 == f_I
c 
c       p = p0/alpha**(1 - fl - fI)
c       q = q0/alpha**(1 - fl - fI)
c
c       dy1/dt = Sin - y1*y4*p    
c       dy2/dt = y1*y4*p - fP    
c       dy3/dt = fP - y3*y4*q - sumI   
c       dy4/dt = y1(t - tau)*y4(t - tau)*p(t - tau) + y3(t - tauI)*y4(t - tau_I)*q(t - tauI) - y1*y4*p - y3*y4*q
c       dy5/dt = y3*y4*q - y3(t - tauI)*y4(t-tau_I)*q(t - tauI)
c
c  CONSTANTES Y VARIABLES:
c  np = Numero de puntos sobre el eje del tiempo.
c  kmax = Numero de pasos que pueden ser guardados.
c  hmin = Paso minimo permitido (puede ser cero)
c  eps = Nivel de tolerancia
c
c  SUBRUTINAS:
c  Odeint:Utiliza el metodo Runge-Kutta con control adaptador del paso.
c  deri5:Calcula las derivadas del integrando (derint=dy/dm)
c  hunt:Dado un arreglo x(i) de longitud N,nos devuelve un valor jlo tal
c       que  x se encuentra entre x(jlo) y x(jlo+1).
c
        implicit double precision (a-h,o-z)
        external deri5,rkqc
        
        parameter(nvar5=5,np=30000)
        common/evovar/y1_inic,y3_inic,y4_inic,y5_inic
        common/evoret/jlo
        common/param/Sdot,sumI,tau,tauI,tauP,p0,q0,alpha,Pcri,oligo,mp
     *,pImax
        dimension ys5(nvar5)
        common/path5/kmax5,kount5,dxsav5,xp5(2000002),yp5(5,2000002)
        common/derivadas/dy1,dy2,dy3,dy4,dy5
        common/mirar/AA,BB,CC,DD,Spunto,Pplus
        open(unit=3,file='salida1.dat')
        open(unit=4,file='salida2.dat')
        open(unit=7,file='mirar.dat')
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        jlo=1
cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp

        if(.true.)read(5,*)Sdot,sumI,tau,tauI,tauP,p0,q0,alpha,Pcri,
     *oligo,mp,pImax,x_fin
c        print *,'Voy por aqui '

        x_ini=0.0
        x_fin=x_fin*tau
        print *,'tau es ',tau
        print *,'El tiempo final de integracion es ',x_fin
c        read(5,*)iparada
c        if(iparada.eq.1) pause
c
c-----------CONDICIONES INICIALES------------------------

        ys5(1)=0.0
        ys5(2)=0.0
        ys5(3)=0.0
        ys5(4)=1.0
        ys5(5)=0.0
C......................................................
        y1_inic = ys5(1)
        y3_inic = ys5(3)
        y4_inic = ys5(4)
        y5_inic = ys5(5)

c--------------- ESCRITURA ENCABEZADOS------------------------------
67	format(t8,
     #'k        t         A      alpha1   alpha2    alpha3'
     #,'        B      beta1   beta2    beta3    S')
c-------------------------------------------------------------------
        eps5=1.d-6
c        h15=(x_fin-x_ini)/20000000.
        h15=(x_fin-x_ini)/100000.
        hmin5=0.0
        kmax5=20000002
c        dxsav5=(x_fin-x_ini)/20000000.
        dxsav5=tau/500.d0
        call ode5(ys5,nvar5,x_ini,x_fin,eps5,h15,hmin5,
     *nok5,nbad5,deri5,rkqc)
        write(*,'(/1x,a,t30,i7)') 'pasos sucesivos: ',nok5
        write(*,'(1x,a,t30,i7)') 'pasos malos: ',nbad5
        write(*,'(/1x,a,t30,i10)') 'valores intermedios almacenados: ',
     *kount5

        stop
        end
c---------------------------------------------------------------
c---------------------------------------------------------------
c---------------------------------------------------------------
c			SUBRUTINAS	
c---------------------------------------------------------------
c---------------------------------------------------------------
c---------------------------------------------------------------
      SUBROUTINE ODE5(YSTART,NVAR,X1,X2,EPS,H1,HMIN,NOK,NBAD,DERIVS,RK
     *QC)
      implicit double precision (a-h,o-z)
c      PARAMETER (MAXSTP=10000,NMAX=10,TWO=2.0,ZERO=0.0,TINY=1.D-30)
      PARAMETER (MAXSTP=20000002,NMAX=10,TWO=2.0,ZERO=0.0,TINY=1.D-20)
      EXTERNAL DERIVS
      common/path5/kmax5,kount5,dxsav5,xp5(2000002),yp5(5,2000002)
c      COMMON /PATH5/ KMAX5,KOUNT5,DXSAV5,XP5(2000002),YP5(10,2000002)
      common/evovar/y1_inic,y3_inic,y4_inic,y5_inic
      common/evoret/jlo
      common/param/Sdot,sumI,tau,tauI,tauP,p0,q0,alpha,Pcri,oligo,mp,
     *pImax
      common/derivadas/dy1,dy2,dy3,dy4,dy5
      common/mirar/AA,BB,CC,DD,Spunto,Pplus

      DIMENSION YSTART(NVAR),YSCAL(NMAX),Y(NMAX),DYDX(NMAX)
      X=X1
      H=SIGN(H1,X2-X1)
      NOK=0
      NBAD=0
      KOUNT5=0

      DO 11 I=1,NVAR
        Y(I)=YSTART(I)
11    CONTINUE

      XSAV=X-DXSAV5*TWO
      
      DO 16 NSTP=1,MAXSTP/sdot
        CALL DERIVS(X,Y,DYDX)
        
        DO 12 I=1,NVAR
          YSCAL(I)=DABS(Y(I))+DABS(H*DYDX(I))+TINY
12      CONTINUE

c       Guarda resultados en xp5, yp5
        IF(KMAX5.GT.0)THEN
          IF(DABS(X-XSAV).GT.DABS(DXSAV5)) THEN
            IF(KOUNT5.LT.KMAX5-1)THEN
              KOUNT5=KOUNT5+1
              xp5(KOUNT5)=X
              DO 13 I=1,NVAR
                yp5(I,KOUNT5)=Y(I)
13            CONTINUE
              XSAV=X
c-------------------------- ESCRITURA ------------------------------
              if(abs(dy3).lt.1.d-60) then
                 dy3copia = 0
              else
                 dy3copia = dy3 
              endif
              if(x.lt.0.5*x2) then
                 write(3,'(5x,i7,10(f12.5,1x))') kount5,x/tau,
     *Y(1),Y(2),Y(3),Y(4),Y(5),Spunto,dy3copia,Pplus,dy2
              else if(x.ge.0.5*x2) then
                 write(4,'(5x,i7,10(f12.5,1x))') kount5,x/tau,
     *Y(1),Y(2),Y(3),Y(4),Y(5),Spunto,dy3copia,Pplus,dy2
              endif
c              write(4,'(5x,i7,f10.4,1x,6(1pd20.13,1x))') kount5,x/tau
c     *,dy1,dy2,dy3,dy4,dy5,Spunto
              write(7,'(5x,i7,f10.4,1x,7(1pd20.13,1x))') kount5,x/tau
     *,AA,BB,CC,DD,AA+BB+CC+DD,AA+BB,CC+DD
            ENDIF
          ENDIF
        ENDIF
        
        IF((X+H-X2)*(X+H-X1).GT.ZERO) H=X2-X
        
        CALL RKQC(Y,DYDX,NVAR,X,H,EPS,YSCAL,HDID,HNEXT,DERIVS)
        
        IF(HDID.EQ.H)THEN
          NOK=NOK+1
        ELSE
          NBAD=NBAD+1
        ENDIF
        
        IF((X-X2)*(X2-X1).GE.ZERO)THEN
          DO 14 I=1,NVAR
            YSTART(I)=Y(I)
14        CONTINUE
          IF(KMAX5.NE.0)THEN
            KOUNT5=KOUNT5+1
            XP5(KOUNT5)=X
            DO 15 I=1,NVAR
              YP5(I,KOUNT5)=Y(I)
15          CONTINUE
          ENDIF
          RETURN
        ENDIF
        
        IF(DABS(HNEXT).LT.HMIN) PAUSE 'Stepsize smaller than minimum.'
        H=HNEXT
c	print *,'H en ODE5  ',H
16    CONTINUE
      PAUSE 'Too many steps.'
      write(6,*) kount5
      RETURN
      END
C----------------------------------------------------------------
C----------------------------------------------------------------
C----------------------------------------------------------------
        subroutine deri5(x,y,dydx)
        implicit double precision (a-h,o-z)
        common/param/Sdot,sumI,tau,tauI,tauP,p0,q0,alpha,Pcri,oligo,mp
     *,pImax
        common/derivadas/dy1,dy2,dy3,dy4,dy5
        common/mirar/AA,BB,CC,DD,Spunto,Pplus
        dimension y(5),dydx(5)

C       Esta subrutina toma el lugar de DERIVS
C==========================================================
c VARIABLES
c       
c       y1 == S
c       y2 == P
c       y3 == I
c       y4 == f_L
c       y5 == f_I
c
c       p = p0/alpha**(1 - fl - fI)
c       q = q0/alpha**(1 - fl - fI)
c
c       dy1/dt = Sin - y1*y4*p    
c       dy2/dt = y1(t - tauP)*y4(t - tauP)*p(t - tauP) - fP    
c       dy3/dt = fP - y3*y4*q - sumI   
c       dy4/dt = y1(t - tau)*y4(t - tau)*p(t - tau) + y3(t - tauI)*y4(t - tau_I)*q(t - tauI) - y1*y4*p - y3*y4*q
c       dy5/dt = y3*y4*q - y3(t - tauI)*y4(t-tau_I)*q(t - tauI)
c
c       retS = y1(t - tau)*y4(t - tau)*p(t - tau)
c       retI = y3(t - tauI)*y4(t-tau_I)*q(t - tauI)
c       retP = y1(t - tauP)*y4(t - tauP)*p(t - tauP)       
c
        call retardo(x,retS,retI,retP)
c-----------------------------------------------------------
c
	y1 = y(1)
	y2 = y(2)
	y3 = y(3)
	y4 = y(4)
	y5 = y(5)

        p = p0/alpha**(1 - y4 - y5)
        q = q0/alpha**(1 - y4 - y5)
c
c       Pplus = tasa de produccion de producto 
c
        Pplus =  y1*y4*p     

        pI = 0.15  
 
        fP = fprod(Pplus,Pcri,pImax)

        pi = 3.141596
c
c       mp = medio periodo
c
        Spunto = 0.5*Sdot*(1 + sin(pi*(x/mp - 0.5)))
c        Spunto = 0.5*Sdot*(1 + sin(pi/2))
c        print *, Spunto,Sdot,sin(pi/2)
c	pause

c        dydx(1) = Sdot - y1*y4*p
        dydx(1) = Spunto - y1*y4*p
        dydx(2) = retP - fP
        dydx(3) = fP - y3*y4*q - sumI
c        dydx(4) = retS + retI - y1*y4*p - y3*y4*q
        dydx(4) = retS + oligo*retI - y1*y4*p - oligo*y3*y4*q
c        dydx(5) = y3*y4*q - retI
        dydx(5) = oligo*(y3*y4*q - retI)

        dy1 = dydx(1)
        dy2 = dydx(2)
        dy3 = dydx(3)
        dy4 = dydx(4)
        dy5 = dydx(5)

        AA = fP
        BB = y2
        CC = -y1*y4*p
        DD = -y3*y4*q

c        write(4,'(3(1pd20.13,1x))') Sdot, - y1*y4*p, dydx(1)
c        write(4,'(3(1pd20.13,1x))') retS + retI, - y1*y4*p- y3*y4*q, dydx(3)
c        write(4,'(5(1pd20.13,1x))') retS, retI,- y1y4*p,-y3*y4*q,dydx(4)

        return
        end
c-------------------------------------------------
c---------------------------------------------------------------
        subroutine retardo(t,retS,retI,retP)
        implicit double precision (a-h,o-z)
        common/evovar/y1_inic,y3_inic,y4_inic,y5_inic
        common/evoret/jlo
        common/param/Sdot,sumI,tau,tauI,tauP,p0,q0,alpha,Pcri,oligo,mp
     *,pImax
        common/path5/kmax5,kount5,dxsav5,xp5(2000002),yp5(5,2000002)
        common/mirar/AA,BB,CC,DD,Spunto,Pplus
c       
c       retS = y1(t - tau)*y4(t - tau)*p(t - tau)
c       retI = y3(t - tauI)*y4(t-tau_I)*q(t - tauI)
c       retP = y1(t - tauP)*y4(t - tauP)*p(t - tauP)       
c
c	print *,' -------------  en retardo t=',t,' -------------'
	call retpm(t-tau,1,retS,dy)
c        print *,' t-tau=',t-tau,' retS',y4_ret_tau,' dy=',dy
	call retpm(t-tauI,2,retI,dy)
c        print *,'    t-tauI=',t-tauI,' retI',y4_ret_tauI
	call retpm(t-tauP,3,retP,dy)

        return
        end
c---------------------------------------------------------------
        subroutine retpm(tfi,kvar,fit,dy)
        implicit double precision (a-h,o-z)
        common/evovar/y1_inic,y3_inic,y4_inic,y5_inic
        common/evoret/jlo
        common/param/Sdot,sumI,tau,tauI,tauP,p0,q0,alpha,Pcri,oligo,mp
     *,pImax
        common/path5/kmax5,kount5,dxsav5,xp5(2000002),yp5(5,2000002)
        common/mirar/AA,BB,CC,DD,Spunto,Pplus
        real*8 coorx(4),coory(4)
c
c       kvar = 1 ===> retS = y1(t - tau)*y4(t - tau)*p(t - tau)
c       kvar = 2 ===> retI = y3(t - tauI)*y4(t-tau_I)*q(t - tauI)
c       kvar = 3 ===> retP = y1(t - tauP)*y4(t - tauP)*p(t - tauP)       
c
c       print *,'jlo antes de hunt=',jlo

        call hunt(xp5,kount5,tfi,jlo)
cc        
	if(tfi.le.0.0.and.kvar.eq.1)then
           p = p0/alpha**(1 - y4_inic - y5_inic)
           fit = y1_inic*y4_inic*p 
	   dy = 0.0
	   return
        else if(tfi.le.0.0.and.kvar.eq.2) then
           q = q0/alpha**(1 - y4_inic - y5_inic)
           fit= y3_inic*y4_inic*q
	   dy = 0.0
	   return
        else if(tfi.le.0.0.and.kvar.eq.3) then
           p = p0/alpha**(1 - y4_inic - y5_inic)
           fit = y1_inic*y4_inic*p 
	   dy = 0.0
	   return
	endif

        if(kount5.lt.4)print *,'ES GIBT EIN PROBLEM!!!!!'
        if(kount5.lt.4)pause
        j=jlo
        np=kount5
         
        if(kvar.eq.1.or.kvar.eq.3) then
	   ret1=yp5(1,1)*yp5(4,1)*p0/alpha**(1 - yp5(4,1) - yp5(5,1))
	   ret2=yp5(1,2)*yp5(4,2)*p0/alpha**(1 - yp5(4,2) - yp5(5,2))
	   ret3=yp5(1,3)*yp5(4,3)*p0/alpha**(1 - yp5(4,3) - yp5(5,3))
	   ret4=yp5(1,4)*yp5(4,4)*p0/alpha**(1 - yp5(4,4) - yp5(5,4))
        else if(kvar.eq.2) then 
	   ret1=yp5(3,1)*yp5(4,1)*q0/alpha**(1 - yp5(4,1) - yp5(5,1))
	   ret2=yp5(3,2)*yp5(4,2)*q0/alpha**(1 - yp5(4,2) - yp5(5,2))
	   ret3=yp5(3,3)*yp5(4,3)*q0/alpha**(1 - yp5(4,3) - yp5(5,3))
	   ret4=yp5(3,4)*yp5(4,4)*q0/alpha**(1 - yp5(4,4) - yp5(5,4))
        endif
    
c        print *,'kvar =',kvar
c        print *,'ret1 =',ret1
c        print *,'ret2 =',ret2
c        print *,'ret3 =',ret3
c        print *,'ret4 =',ret4
c        print *,'ret5 =',ret5
c        pause

        if(j.eq.0) then
c           print *,'en retma extrapolation    j=0'
           coorx(1) = xp5(1)
           coorx(2) = xp5(2)
           coorx(3) = xp5(3)
           coorx(4) = xp5(4)
           coory(1) = ret1
           coory(2) = ret2
           coory(3) = ret3
           coory(4) = ret4
        else if(j.le.2.) then
           coorx(1) = xp5(1)
           coorx(2) = xp5(2)
           coorx(3) = xp5(3)
           coorx(4) = xp5(4)
           coory(1) = ret1
           coory(2) = ret2
           coory(3) = ret3
           coory(4) = ret4
        else if(j.eq.np) then
c           print *,'en retma extrapolation    j=np'
           coorx(1) = xp5(np-3)
           coorx(2) = xp5(np-2)
           coorx(3) = xp5(np-1)
           coorx(4) = xp5(np)
           if(kvar.eq.1.or.kvar.eq.3) then
              coory(1)=yp5(1,np-3)*yp5(4,np-3)*p0/alpha**(1 - 
     *yp5(4,np-3) - yp5(5,np-3))
              coory(2)=yp5(1,np-2)*yp5(4,np-2)*p0/alpha**(1 - 
     *yp5(4,np-2) - yp5(5,np-2))
              coory(3)=yp5(1,np-1)*yp5(4,np-1)*p0/alpha**(1 - 
     *yp5(4,np-1) - yp5(5,np-1))
              coory(4)=yp5(1,np)*yp5(4,np)*p0/alpha**(1 - 
     *yp5(4,np) - yp5(5,np))
           else if(kvar.eq.2) then
              coory(1)=yp5(3,np-3)*yp5(4,np-3)*q0/alpha**(1 -
     *yp5(4,np-3) - yp5(5,np-3))
              coory(2)=yp5(3,np-2)*yp5(4,np-2)*q0/alpha**(1 - 
     *yp5(4,np-2) - yp5(5,np-2))
              coory(3)=yp5(3,np-1)*yp5(4,np-1)*q0/alpha**(1 - 
     *yp5(4,np-1) - yp5(5,np-1))
              coory(4)=yp5(3,np)*yp5(4,np)*q0/alpha**(1 - 
     *yp5(4,np) - yp5(5,np))
           endif
        else if(j.ge.np-1) then
           coorx(1)=xp5(np-3)
           coorx(2)=xp5(np-2)
           coorx(3)=xp5(np-1)
           coorx(4)=xp5(np)
           if(kvar.eq.1.or.kvar.eq.3) then
              coory(1)=yp5(1,np-3)*yp5(4,np-3)*p0/alpha**(1 - 
     *yp5(4,np-3) - yp5(5,np-3))
              coory(2)=yp5(1,np-2)*yp5(4,np-2)*p0/alpha**(1 - 
     *yp5(4,np-2) - yp5(5,np-2))
              coory(3)=yp5(1,np-1)*yp5(4,np-1)*p0/alpha**(1 - 
     *yp5(4,np-1) - yp5(5,np-1))
              coory(4)=yp5(1,np)*yp5(4,np)*p0/alpha**(1 - 
     *yp5(4,np) - yp5(5,np))
           else if(kvar.eq.2) then
              coory(1)=yp5(3,np-3)*yp5(4,np-3)*q0/alpha**(1 - 
     *yp5(4,np-3) - yp5(5,np-3))
              coory(2)=yp5(3,np-2)*yp5(4,np-2)*q0/alpha**(1 - 
     *yp5(4,np-2) - yp5(5,np-2))
              coory(3)=yp5(3,np-1)*yp5(4,np-1)*q0/alpha**(1 - 
     *yp5(4,np-1) - yp5(5,np-1))
              coory(4)=yp5(3,np)*yp5(4,np)*q0/alpha**(1 - 
     *yp5(4,np) - yp5(5,np))
           endif 
        else
           coorx(1)=xp5(j-1)
           coorx(2)=xp5(j)
           coorx(3)=xp5(j+1)
           coorx(4)=xp5(j+2)
           if(kvar.eq.1.or.kvar.eq.3) then
              coory(1)=yp5(1,j-1)*yp5(4,j-1)*p0/alpha**(1 - 
     *yp5(4,j-1) - yp5(5,j-1))
              coory(2)=yp5(1,j)*yp5(4,j)*p0/alpha**(1 - 
     *yp5(4,j) - yp5(5,j))
              coory(3)=yp5(1,j+1)*yp5(4,j+1)*p0/alpha**(1 - 
     *yp5(4,j+1) - yp5(5,j+1))
              coory(4)=yp5(1,j+2)*yp5(4,j+2)*p0/alpha**(1 - 
     *yp5(4,j+2) - yp5(5,j+2))
           else if(kvar.eq.2) then
              coory(1)=yp5(3,j-1)*yp5(4,j-1)*q0/alpha**(1 - 
     *yp5(4,j-1) - yp5(5,j-1))
              coory(2)=yp5(3,j)*yp5(4,j)*q0/alpha**(1 - 
     *yp5(4,j) - yp5(5,j))
              coory(3)=yp5(3,j+1)*yp5(4,j+1)*q0/alpha**(1 - 
     *yp5(4,j+1) - yp5(5,j+1))
              coory(4)=yp5(3,j+2)*yp5(4,j+2)*q0/alpha**(1 - 
     *yp5(4,j+2) - yp5(5,j+2))
           endif 
        endif
c
        call polint(coorx,coory,4,tfi,fit,dy)
c
c        if(dabs(dy/fit).gt.0.05)print *,'en retma tfi,fit,dy=',tfi,fit,dy
c        if(fit.lt.0.0)print *,'//en retma/// INTERPOLACION < 0  fit=',fit
c        if(fit.lt.0.0)print *,'   INTERPOLACION < 0  fit puesta a 0 en
c     #kvar =',kvar
        if(fit.lt.0.0)fit=0.0
        return
        end
c---------------------------------------------------------------
c---------------------------------------------------------------
c---------------------------------------------------------------
c
c  hunt:Dado un arreglo x(i) de longitud N,nos devuelve un valor jlo tal
c       que  x se encuentra entre x(jlo) y x(jlo+1).
c
      SUBROUTINE HUNT(XX,N,X,JLO)
      implicit double precision (a-h,o-z)
      DIMENSION XX(N)
      LOGICAL ASCND
      ASCND=XX(N).GT.XX(1)
      IF(JLO.LE.0.OR.JLO.GT.N)THEN
        JLO=0
        JHI=N+1
        GO TO 3
      ENDIF
      INC=1
      IF(X.GE.XX(JLO).EQV.ASCND)THEN
1       JHI=JLO+INC
        IF(JHI.GT.N)THEN
          JHI=N+1
        ELSE IF(X.GE.XX(JHI).EQV.ASCND)THEN
          JLO=JHI
          INC=INC+INC
          GO TO 1
        ENDIF
      ELSE
        JHI=JLO
2       JLO=JHI-INC
        IF(JLO.LT.1)THEN
          JLO=0
        ELSE IF(X.LT.XX(JLO).EQV.ASCND)THEN
          JHI=JLO
          INC=INC+INC
          GO TO 2
        ENDIF
      ENDIF
3     IF(JHI-JLO.EQ.1)goto 4
      JM=(JHI+JLO)/2
      IF(X.GT.XX(JM).EQV.ASCND)THEN
        JLO=JM
      ELSE
        JHI=JM
      ENDIF
      GO TO 3
 4      return
      END
C----------------------------------------------------------------
      SUBROUTINE RKQC(Y,DYDX,N,X,HTRY,EPS,YSCAL,HDID,HNEXT,DERIVS)
      implicit double precision (a-h,o-z)
      PARAMETER (NMAX=8,FCOR=.0666666667,
     *    ONE=1.,SAFETY=0.9,ERRCON=6.D-4)
      EXTERNAL DERIVS
      DIMENSION Y(N),DYDX(N),YSCAL(N),YTEMP(NMAX),YSAV(NMAX),DYSAV(NMAX)
      PGROW=-0.20
      PSHRNK=-0.25
      XSAV=X
      DO 11 I=1,N
        YSAV(I)=Y(I)
        DYSAV(I)=DYDX(I)
11    CONTINUE
      H=HTRY
1     HH=0.5*H
      CALL RK4(YSAV,DYSAV,N,XSAV,HH,YTEMP,DERIVS)
      X=XSAV+HH
      CALL DERIVS(X,YTEMP,DYDX)
      CALL RK4(YTEMP,DYDX,N,X,HH,Y,DERIVS)
      X=XSAV+H
      IF(X.EQ.XSAV)PAUSE 'Stepsize not significant in RKQC.'
      CALL RK4(YSAV,DYSAV,N,XSAV,H,YTEMP,DERIVS)
      ERRMAX=0.
      DO 12 I=1,N
        YTEMP(I)=Y(I)-YTEMP(I)
        ERRMAX=MAX(ERRMAX,DABS(YTEMP(I)/YSCAL(I)))
12    CONTINUE
      ERRMAX=ERRMAX/EPS
      IF(ERRMAX.GT.ONE) THEN
c	print *,'XSAV,H,XSAV+HH,XSAV+H,ERRMAX',XSAV,H,XSAV+HH,XSAV+H,ERRMAX
c	print *,'YSCAL(I)',YSCAL
c        print *,'YTEMP(I)',YTEMP
        H=SAFETY*H*(ERRMAX**PSHRNK)
        GOTO 1
      ELSE
        HDID=H
        IF(ERRMAX.GT.ERRCON)THEN
          HNEXT=SAFETY*H*(ERRMAX**PGROW)
        ELSE
          HNEXT=4.*H
        ENDIF
      ENDIF
      DO 13 I=1,N
        Y(I)=Y(I)+YTEMP(I)*FCOR
13    CONTINUE
      RETURN
      END
C----------------------------------------------------------------
      SUBROUTINE RK4(Y,DYDX,N,X,H,YOUT,DERIVS)
      implicit double precision (a-h,o-z)
      PARAMETER (NMAX=10)
      DIMENSION Y(N),DYDX(N),YOUT(N),YT(NMAX),DYT(NMAX),DYM(NMAX)
      HH=H*0.5
      H6=H/6.
      XH=X+HH
      DO 11 I=1,N
        YT(I)=Y(I)+HH*DYDX(I)
11    CONTINUE
      CALL DERIVS(XH,YT,DYT)
      DO 12 I=1,N
        YT(I)=Y(I)+HH*DYT(I)
12    CONTINUE
      CALL DERIVS(XH,YT,DYM)
      DO 13 I=1,N
        YT(I)=Y(I)+H*DYM(I)
        DYM(I)=DYT(I)+DYM(I)
13    CONTINUE
      CALL DERIVS(X+H,YT,DYT)
      DO 14 I=1,N
        YOUT(I)=Y(I)+H6*(DYDX(I)+DYT(I)+2.*DYM(I))
14    CONTINUE
      RETURN
      END
C----------------------------------------------------------------
      SUBROUTINE POLINT(XA,YA,N,X,Y,DY)
      implicit double precision (a-h,o-z)
      PARAMETER (NMAX=10)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      DIF=DABS(X-XA(1))
      DO 11 I=1,N
        DIFT=DABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF(DEN.EQ.0.)PAUSE
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END

C------------------------------------------------------------------      
C------------------------------------------------------------------      
c      SUBROUTINE FPROD(x,y)     
c      implicit double precision (a-h,o-z)
c        common/param/Sdot,sumI,tau,tauI,tauP,p0,q0,alpha,Pcri,oligo,mp,pImax
c     
c      if(x.lt.Pcri) then
c         y = x + .5
c      else if(x.ge.Pcri) then
c         y = 0.0
c      endif 
c
c      return
c      end
C------------------------------------------------------------------      
C------------------------------------------------------------------      
      function fprod(x,valcri,pImax)
      implicit double precision (a-h,o-z)
c
c      if(x.lt.valcri) then
c         fprod = (0.5/valcri)*x
c      else if(x.ge.valcri) then
c         fprod = 0.0
c      endif 
c      if(x.lt.valcri/2) then
c         fprod = 0.001
c      else if(x.ge.valcri/2.and.x.lt.valcri) then
c         fprod = (5/valcri)*x
c      else if(x.ge.valcri) then
c         fprod = 0.0
c      endif 
c      if(x.lt.valcri) then
c         fprod =  0.2/(1 + exp(-(x - valcri/2)*80))
c      else if(x.ge.valcri) then
c         fprod = 0
c      endif
       fprod =  pImax * ( 1/(1 + exp(-(x - 0.5*valcri)*500)) + 1/(1 + 
     *exp((x - valcri)*500)) - 1)
      return
      end
 

