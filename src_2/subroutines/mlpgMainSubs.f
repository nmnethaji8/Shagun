!!--------------------------PRESSURE_SOLVER2--------------------------!!
      SUBROUTINE PRESSURE_SOLVER2(NODN,PTMP,FB,EPS_G,ERRSOL)    
      use, intrinsic :: ISO_C_BINDING, only : C_INT, C_PTR, C_DOUBLE
      use, intrinsic :: ISO_C_BINDING, only : C_CHAR, C_NULL_CHAR, C_LOC

      USE COMMONMOD
      USE MLPGSTORAGE
      USE NEIGHNODES
      IMPLICIT NONE
      !INCLUDE 'COMMON.F'
      

      interface
      subroutine paralution_stop() BIND(C)
      end subroutine paralution_stop    

      subroutine paralution_fortran_solve_csr( n, m, nnz, solver, 
     +  mformat, preconditioner, pformat,
     +  ivCSR, jvCSR, rval, rhs, atol, rtol, div, maxiter, basis,
     +  p, q, x, iter, resnorm, ierr ) BIND(C)

        use, intrinsic :: ISO_C_BINDING, only : C_INT, C_PTR, 
     +    C_DOUBLE, C_CHAR

        integer(kind=C_INT), value, intent(in)  :: n, m, nnz, maxiter
        integer(kind=C_INT), value, intent(in)  :: basis, p, q
        real(kind=C_DOUBLE), value, intent(in)  :: atol, rtol, div
        integer(kind=C_INT),        intent(out) :: iter, ierr
        real(kind=C_DOUBLE),        intent(out) :: resnorm
        type(C_PTR),         value, intent(in)  :: ivCSR, jvCSR, rval
        type(C_PTR),         value, intent(in)  :: rhs
        type(C_PTR),         value              :: x
        character(kind=C_CHAR)                  :: solver, mformat
        character(kind=C_CHAR)                  :: preconditioner
        character(kind=C_CHAR)                  :: pformat

      end subroutine paralution_fortran_solve_csr
      end interface

      INTEGER(KIND=4),INTENT(IN)::NODN      
      REAL(KIND=8),INTENT(IN)::EPS_G

      INTEGER(KIND=4),INTENT(OUT)::ERRSOL
      REAL(KIND=8),INTENT(INOUT)::FB(LNODE),PTMP(LNODE)      
      
      INTEGER(KIND=4)::IER,ITER,INNZ,I,IX,IX2,IJ,IK,IL
      REAL(KIND=8)::RESNORM,TMPR7



      ! Shagun edit CSR      
      ! THE FOLLOWING IS BECUASE KX+B=0 IS USED IN GAUSS_SEDIEL_1
      ! BUT KX=B USED IN GMRESWHOLEMA.
      FBCSR=-FB(1:NODN)      
C       DO IX=1,NODN
C         FBCSR(IX)=-FB(IX)
C         !PTCSR(IX)=PTMP(IX)
C       ENDDO
      ! PTMP=0D0
      WRITE(8,*)'Entering BiCGStab CSR'
      INNZ=0
      DO I=1,NODN
        INNZ=INNZ+IVV(I)
      ENDDO      
      IX2=0
      IVVCSR(1)=1      
      DO IX=1,NODN
        IVVCSR(IX+1)=IVVCSR(IX)+IVV(IX)
        IK=(IX-1)*IVV(0)
        TMPR7=SKK(IK+IVV(IX))
        FBCSR(IX)=FBCSR(IX)/TMPR7
        DO IJ=1,IVV(IX)
          IL=LINKTAB(IK+IJ)
          IX2=IX2+1
          JVVCSR(IX2)=IL
          SKKCSR(IX2)=SKK(IK+IJ)/TMPR7
        ENDDO
      ENDDO
      IF((IX2.NE.INNZ).OR.(IVVCSR(NODN+1).NE.(INNZ+1)))THEN
        WRITE(8,*)' [ERR] CHECK CSR MATRICES'
        WRITE(8,*)INNZ,IX2,IVVCSR(NODN+1)
        STOP
      ENDIF      


C       IX=NODN/2+1000
C       WRITE(8,*)'TEXT'
C       WRITE(8,*)IX,IVV(IX),IVVCSR(IX+1),IVVCSR(IX)
C       WRITE(8,*)'TEXT'
C       WRITE(8,*)JVVCSR(IVVCSR(IX):IVVCSR(IX+1)-1)
C       WRITE(8,*)'TEXT'
C       WRITE(8,*)SKKCSR(IVVCSR(IX):IVVCSR(IX+1)-1)
C       IK=(IX-1)*IVV(0)
C       WRITE(8,*)'TEXT'
C       WRITE(8,*)SKK(IK+1:IK+IVV(IX))
C       STOP

      ERRSOL = 0 !NO ERROR IN SOLVER
      
      call paralution_fortran_solve_csr( NODN, NODN, INNZ, 
     +    'BiCGStab' // C_NULL_CHAR,
     +    'CSR' // C_NULL_CHAR, 
     +    'Jacobi' // C_NULL_CHAR, 
     +    'CSR' // C_NULL_CHAR, 
     +    C_LOC(IVVCSR), C_LOC(JVVCSR), 
     +    C_LOC(SKKCSR), C_LOC(FBCSR), 
     +    1e-10_C_DOUBLE, 1e-15_C_DOUBLE, 1e+8_C_DOUBLE, 20000, 
     +    30, 0, 1, C_LOC(PTCSR), iter, resnorm, ier )
      WRITE(8,*)"[PARACSR]",ier,iter,resnorm

      if((ier.ne.1).and.(ier.ne.2)) then
        WRITE(8,*)'[ERR] BiCGStab Failed. Shifting to GMRES',ier

        call paralution_fortran_solve_csr( NODN, NODN, INNZ, 
     +    'GMRES' // C_NULL_CHAR,
     +    'CSR' // C_NULL_CHAR, 
     +    'None' // C_NULL_CHAR, 
     +    'CSR' // C_NULL_CHAR, 
     +    C_LOC(IVVCSR), C_LOC(JVVCSR), 
     +    C_LOC(SKKCSR), C_LOC(FBCSR), 
     +    1e-10_C_DOUBLE, 1e-15_C_DOUBLE, 1e+8_C_DOUBLE, 20000, 
     +    30, 0, 1, C_LOC(PTCSR), iter, resnorm, ier )
        WRITE(8,*)"[PARACSR]",ier,iter,resnorm

        if((ier.ne.1).and.(ier.ne.2)) then
          WRITE(8,*)'[Err] Paralution error in  ier Time ::',ier
          CALL CHECK_NUMNE()
          ERRSOL = 1
          RETURN
          !stop
        endif
      endif         

      PTMP(1:NODN)=PTCSR
      
      WRITE(8,*)'Finished GMRES CSR'
      
      END SUBROUTINE PRESSURE_SOLVER2
!!--------------------------PRESSURE_SOLVER2--------------------------!!



!!---------------------------U_UPDATE_POW2----------------------------!!
      SUBROUTINE U_UPDATE_POW2(LNODE,NODN,NODEID,NWALLID,GRA,DT,
     +  COORX,COORY,COORZ,SPONGEX,DOMX,DOMY,DOMZ)
      USE MLPGKINE      
      IMPLICIT NONE

      INTEGER(KIND=4),INTENT(IN)::LNODE,NODN,NODEID(-2:NODN)
      INTEGER(KIND=4),INTENT(IN)::NWALLID(LNODE,4)
      INTEGER(KIND=4)::INOD,NODTAL

      REAL(KIND=8),INTENT(IN)::SPONGEX,DOMX(2),DOMY(2),DOMZ(2)      
      REAL(KIND=8),INTENT(IN)::GRA,DT
      REAL(KIND=8),INTENT(IN)::COORX(LNODE),COORY(LNODE),COORZ(LNODE)
      REAL(KIND=8)::BMAX,XI,S1,COEF2,A_MAX

      NODTAL=NODEID(-1)  !ONLY THE WATER PARTICLE UPDATE
      A_MAX=50.D0*ROU(1)*(-GRA) !AS THE maximum acceleration

      BMAX=10D0   !! DISSIPATION COEFF FOR SPONGE LAYER      
      WRITE(8,*)'[INF] Var SPONGEX =',SPONGEX
      WRITE(8,*)'[INF] D_ZONE Width =',DOMX(2)-SPONGEX

      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(INOD,XI,S1,COEF2)
      !$OMP DO SCHEDULE(DYNAMIC,100)
      DO INOD=1,NODTAL
    
        ! THE VELOCITY IS UPDATE  
        IF(ABS(PPX(INOD,1)).LT.A_MAX)THEN
          UX(INOD,4)=-DT/ROU(INOD)*PPX(INOD,1)
        ELSE
          UX(INOD,4)=0.D0
        ENDIF
        IF(ABS(PPY(INOD,1)).LT.A_MAX)THEN
          UY(INOD,4)=-DT/ROU(INOD)*PPY(INOD,1)     
        ELSE
          UY(INOD,4)=0.D0
        ENDIF
        IF(ABS(PPZ(INOD,1)).LT.A_MAX)THEN
          ! UZ(INOD,4)=-DT/ROU(INOD)*(PPZ(INOD,1)+1)
          ! UZ(INOD,4)=-DT/ROU(INOD)*PPZ(INOD,1)
          UZ(INOD,4)=-DT/ROU(INOD)*PPZ(INOD,1)+GRA*DT
        ELSE
          UZ(INOD,4)=DT*GRA    !gravity term
        ENDIF

        UX(INOD,3)=UX(INOD,2)+UX(INOD,4)
        UY(INOD,3)=UY(INOD,2)+UY(INOD,4)
        UZ(INOD,3)=UZ(INOD,2)+UZ(INOD,4)

        UX(INOD,4)=UX(INOD,1)  !SAVE THE N TIME STEP,FOR COORDINATE UPDATE
        UY(INOD,4)=UY(INOD,1)
        UZ(INOD,4)=UZ(INOD,1)

        UX(INOD,1)=UX(INOD,3)  !SAVE N+1 TIME STEP;
        UY(INOD,1)=UY(INOD,3)
        UZ(INOD,1)=UZ(INOD,3)


        ! ADD THE DAMPING ZONE
        XI=COORX(INOD)
        IF(XI.GE.SPONGEX)THEN          
          S1 = (XI-SPONGEX)/(DOMX(2)-SPONGEX) 
          ! COEF2=0.5*1.*(1-cos(PI*(XI-D_ZONE)/3.))
          ! COEF2 = 0.5D0*0.4*(1+sin(PI*((Xi-d_zone)/3.0)-0.5))  
          ! COEF2 =  s1**6 !-2* (1-S1)**3+3*(1-S1)**2  !1-sin(1.5708d0*S1) !tanh(3.14d0*S1) !  !!  ! ! !
          ! COEF2=1D0 - 3D0*S1**2 + 2D0*S1**3          
          ! COEF2=(1d0-S1**6)
          COEF2=1D0+DT*BMAX*S1**3
          !WRITE(21,'(3F15.6)')XI,S1,COEF2

          UX(INOD,1)=UX(INOD,1)/COEF2          
          UY(INOD,1)=UY(INOD,1)/COEF2
          UZ(INOD,1)=UZ(INOD,1)/COEF2
        ENDIF
        IF (NWALLID(INOD,2).EQ.-10) THEN
          UX(INOD,1:3)=0.D0
          Uz(INOD,1:3)=0.D0
          Uy(INOD,1:3)=0.D0
        END IF

      ENDDO
      !$OMP END DO NOWAIT
      !$OMP END PARALLEL

      RETURN      
      END SUBROUTINE U_UPDATE_POW2
!!-------------------------END U_UPDATE_POW2--------------------------!!



!!-----------------------------U_BOUNDARY-----------------------------!!
      ! UPDATE THE VELOCITY AT N+1 time step, 
      ! satisfy the velocity boundary condition.
      ! SO the locations of water particles should be n+1 time step.
      SUBROUTINE U_BOUNDARY(WMV,IFSI)
      USE COMMONMOD
      USE MLPGKINE
      IMPLICIT NONE
      !INCLUDE 'COMMON.F'      
       
      INTEGER(KIND=4),INTENT(IN)::IFSI
      REAL(KIND=8),INTENT(IN)::WMV

      INTEGER(KIND=4)::INOD,KK,KK1,IK,IK1
      REAL(KIND=8)::U_N,U_M,U_S

      !-----SMOOTH THE V VELOCITY OF THE WALL PARTICLE----
      ! DEAL WITH VELOCITY UPDATE OF THE VERTICAL WALL PARTICLE
      !	DU=ATAN(SLOPE)
      !	SIN1=SIN(DU)
      !	COS1=COS(DU)
      IK = 1
      IK1 = 3
      !	UX(:,IK) = UX(:,3)
      !	UY(:,IK) = UY(:,3)
      !	UZ(:,IK) = UZ(:,3)
  
      
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(INOD,U_N,U_M,U_S,KK,KK1)          
      !$OMP DO SCHEDULE(DYNAMIC,100)

      DO 10 INOD=NODEID(-2)+1,NODEID(-1)
        KK=NWALLID(INOD,1)  !COEFFICIENT ABOUT VELOCITY OF WALL PARTICLES
        KK1=NWALLID(INOD,2)
        IF (NODEID(INOD).EQ.8.AND.IFSI.EQ.1) GOTO 10
        IF(KK1.EQ.-10)THEN
          UX(INOD,IK)=0.D0;UY(INOD,IK)=0.D0;UZ(INOD,IK)=0.D0
          UX(INOD,IK1)=0.D0;UY(INOD,IK1)=0.D0;UZ(INOD,IK1)=0.D0
          GOTO 10
        ENDIF

        IF(KK.EQ.0)THEN !0
          !CALL INTERPOLATION_MPSN(UN(:,1),UN(:,2),UN(:,3),INOD,INOD)
          UX(INOD,IK)=0 !-TANK_A(INOD,1)*DT
          UY(INOD,IK)=0 !-TANK_A(INOD,2)*DT
          UZ(INOD,IK)=0 !-TANK_A(INOD,3)*DT   !NORMAL TO THE WALL
          UX(INOD,IK1)=0 !-TANK_A(INOD,1)*DT
          UY(INOD,IK1)=0 !-TANK_A(INOD,2)*DT
          UZ(INOD,IK1)=0 !-TANK_A(INOD,3)*DT   !NORMAL TO THE WALL
        ENDIF
        !	 U_N = UX(INOD,IK)*SNX(INOD)+UY(INOD,IK)*SNY(INOD)+
        ! &       UZ(INOD,IK)*SNZ(INOD)
        U_N = 0.D0
        IF(NODEID(INOD).EQ.8)THEN     
          IF(I_WM.EQ.1.OR.I_WM.EQ.2.OR.I_WM.EQ.3.OR.
     +      I_WM.EQ.5.OR.I_WM.EQ.6)THEN
            U_N =WMV*SNX(INOD)
          ENDIF          
          IF(I_WM.EQ.15)THEN
            !U_N=TEMP_UN(INOD-NODEID(-2),1)*SNX(INOD)
            CYCLE
          ENDIF    
        END IF
      
     
        U_M = UX(INOD,IK)*SMX(INOD)+UY(INOD,IK)*SMY(INOD)+
     +        UZ(INOD,IK)*SMZ(INOD)
        U_S = UX(INOD,IK)*SSX(INOD)+UY(INOD,IK)*SSY(INOD)+
     +        UZ(INOD,IK)*SSZ(INOD)
        
        UX(INOD,IK) = U_N*SNX(INOD)+U_M*SMX(INOD)+U_S*SSX(INOD)
        UY(INOD,IK) = U_N*SNY(INOD)+U_M*SMY(INOD)+U_S*SSY(INOD)
        UZ(INOD,IK) = U_N*SNZ(INOD)+U_M*SMZ(INOD)+U_S*SSZ(INOD)

        U_M = UX(INOD,IK1)*SMX(INOD)+UY(INOD,IK1)*SMY(INOD)+
     +        UZ(INOD,IK1)*SMZ(INOD)
        U_S = UX(INOD,IK1)*SSX(INOD)+UY(INOD,IK1)*SSY(INOD)+
     +        UZ(INOD,IK1)*SSZ(INOD)
        
        UX(INOD,IK1) = U_N*SNX(INOD)+U_M*SMX(INOD)+U_S*SSX(INOD)
        UY(INOD,IK1) = U_N*SNY(INOD)+U_M*SMY(INOD)+U_S*SSY(INOD)
        UZ(INOD,IK1) = U_N*SNZ(INOD)+U_M*SMZ(INOD)+U_S*SSZ(INOD)
        U_N = 0.D0
        U_M = PPX(INOD,3)*SMX(INOD)+PPY(INOD,3)*SMY(INOD)+
     +        PPZ(INOD,3)*SMZ(INOD)
        U_S = PPX(INOD,3)*SSX(INOD)+PPY(INOD,3)*SSY(INOD)+
     +        PPZ(INOD,3)*SSZ(INOD)
        
        PPX(INOD,3) = U_N*SNX(INOD)+U_M*SMX(INOD)+U_S*SSX(INOD)
        PPY(INOD,3) = U_N*SNY(INOD)+U_M*SMY(INOD)+U_S*SSY(INOD)
        PPZ(INOD,3) = U_N*SNZ(INOD)+U_M*SMZ(INOD)+U_S*SSZ(INOD) 
        IF(NODEID(INOD).EQ.8.or.nodeid(inod).eq.2)THEN     
          UY(INOD,IK) = 0.D0
          UY(INOD,IK1) = 0.D0
        ENDIF

10    END DO
      !$OMP END DO NOWAIT
      !$OMP END PARALLEL 
      
      END SUBROUTINE U_BOUNDARY



      SUBROUTINE U_BOUNDARY2(IFSI)
      USE COMMONMOD
      USE MLPGKINE
      IMPLICIT NONE
      !INCLUDE 'COMMON.F'      
       
      INTEGER(KIND=4),INTENT(IN)::IFSI

      INTEGER(KIND=4)::INOD,KK,KK1,IK,IK1
      REAL(KIND=8)::U_N,U_M,U_S

      !-----SMOOTH THE V VELOCITY OF THE WALL PARTICLE----
      ! DEAL WITH VELOCITY UPDATE OF THE VERTICAL WALL PARTICLE
      ! DU=ATAN(SLOPE)
      ! SIN1=SIN(DU)
      ! COS1=COS(DU)
      IK = 1
      IK1 = 3
      ! UX(:,IK) = UX(:,3)
      ! UY(:,IK) = UY(:,3)
      ! UZ(:,IK) = UZ(:,3)
  
      
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(INOD,U_N,U_M,U_S,KK,KK1)          
      !$OMP DO SCHEDULE(DYNAMIC,100)

      DO 10 INOD=NODEID(-2)+1,NODEID(-1)
        KK=NWALLID(INOD,1)  !COEFFICIENT ABOUT VELOCITY OF WALL PARTICLES
        KK1=NWALLID(INOD,2)
        IF (NODEID(INOD).EQ.8.AND.IFSI.EQ.1) GOTO 10
        IF(KK1.EQ.-10)THEN
          UX(INOD,IK)=0.D0;UY(INOD,IK)=0.D0;UZ(INOD,IK)=0.D0
          UX(INOD,IK1)=0.D0;UY(INOD,IK1)=0.D0;UZ(INOD,IK1)=0.D0
          GOTO 10
        ENDIF

        IF((KK.EQ.0) .OR. (NODEID(INOD).EQ.2))THEN !0
          !CALL INTERPOLATION_MPSN(UN(:,1),UN(:,2),UN(:,3),INOD,INOD)
          UX(INOD,IK)=0 !-TANK_A(INOD,1)*DT
          UY(INOD,IK)=0 !-TANK_A(INOD,2)*DT
          UZ(INOD,IK)=0 !-TANK_A(INOD,3)*DT   !NORMAL TO THE WALL
          UX(INOD,IK1)=0 !-TANK_A(INOD,1)*DT
          UY(INOD,IK1)=0 !-TANK_A(INOD,2)*DT
          UZ(INOD,IK1)=0 !-TANK_A(INOD,3)*DT   !NORMAL TO THE WALL

          CYCLE
        ENDIF      


        !  U_N = UX(INOD,IK)*SNX(INOD)+UY(INOD,IK)*SNY(INOD)+
        ! &       UZ(INOD,IK)*SNZ(INOD)
        U_N = 0.D0
        IF(NODEID(INOD).EQ.8)THEN     
          IF(I_WM.EQ.1.OR.I_WM.EQ.2.OR.I_WM.EQ.3.OR.
     +      I_WM.EQ.5.OR.I_WM.EQ.6)THEN
            !U_N =WMV*SNX(INOD)
          ENDIF          
          IF(I_WM.EQ.15)THEN
            !U_N=TEMP_UN(INOD-NODEID(-2),1)*SNX(INOD)
            CYCLE
          ENDIF    
        END IF
      
     
        U_M = UX(INOD,IK)*SMX(INOD)+UY(INOD,IK)*SMY(INOD)+
     +        UZ(INOD,IK)*SMZ(INOD)
        U_S = UX(INOD,IK)*SSX(INOD)+UY(INOD,IK)*SSY(INOD)+
     +        UZ(INOD,IK)*SSZ(INOD)
        
        UX(INOD,IK) = U_N*SNX(INOD)+U_M*SMX(INOD)+U_S*SSX(INOD)
        UY(INOD,IK) = U_N*SNY(INOD)+U_M*SMY(INOD)+U_S*SSY(INOD)
        UZ(INOD,IK) = U_N*SNZ(INOD)+U_M*SMZ(INOD)+U_S*SSZ(INOD)

        U_M = UX(INOD,IK1)*SMX(INOD)+UY(INOD,IK1)*SMY(INOD)+
     +        UZ(INOD,IK1)*SMZ(INOD)
        U_S = UX(INOD,IK1)*SSX(INOD)+UY(INOD,IK1)*SSY(INOD)+
     +        UZ(INOD,IK1)*SSZ(INOD)
        
        UX(INOD,IK1) = U_N*SNX(INOD)+U_M*SMX(INOD)+U_S*SSX(INOD)
        UY(INOD,IK1) = U_N*SNY(INOD)+U_M*SMY(INOD)+U_S*SSY(INOD)
        UZ(INOD,IK1) = U_N*SNZ(INOD)+U_M*SMZ(INOD)+U_S*SSZ(INOD)
        U_N = 0.D0
        U_M = PPX(INOD,3)*SMX(INOD)+PPY(INOD,3)*SMY(INOD)+
     +        PPZ(INOD,3)*SMZ(INOD)
        U_S = PPX(INOD,3)*SSX(INOD)+PPY(INOD,3)*SSY(INOD)+
     +        PPZ(INOD,3)*SSZ(INOD)
        
        PPX(INOD,3) = U_N*SNX(INOD)+U_M*SMX(INOD)+U_S*SSX(INOD)
        PPY(INOD,3) = U_N*SNY(INOD)+U_M*SMY(INOD)+U_S*SSY(INOD)
        PPZ(INOD,3) = U_N*SNZ(INOD)+U_M*SMZ(INOD)+U_S*SSZ(INOD) 
        IF(NODEID(INOD).EQ.8)THEN     
          UY(INOD,IK) = 0.D0
          UY(INOD,IK1) = 0.D0
        ENDIF

10    END DO
      !$OMP END DO NOWAIT
      !$OMP END PARALLEL 
      
      END SUBROUTINE U_BOUNDARY2
!!---------------------------END U_BOUNDARY---------------------------!!



!!-----------------------------BOTSLIPBC------------------------------!!
      SUBROUTINE BOTSLIPBC(LNODE, NODN, NODEID, NWALLID, 
     +  CORX, CORY, CORZ, DDR, SNX, SNY, SNZ, SMX, SMY, SMZ,
     +  SSX, SSY, SSZ, UX, UY, UZ)
      USE NEIGHNODES
      IMPLICIT NONE
      !INCLUDE 'COMMON.F'

      !NODN SHOULD BE NODEID(0)
      INTEGER(KIND=4),INTENT(IN)::LNODE, NODN, NODEID(-7:NODN)
      INTEGER(KIND=4),INTENT(INOUT)::NWALLID(LNODE,4)
      REAL(KIND=8),INTENT(IN)::SNX(LNODE), SNY(LNODE), SNZ(LNODE)
      REAL(KIND=8),INTENT(IN)::SMX(LNODE), SMY(LNODE), SMZ(LNODE)
      REAL(KIND=8),INTENT(IN)::SSX(LNODE), SSY(LNODE), SSZ(LNODE)
      REAL(KIND=8),INTENT(IN)::DDR(NODN)
      REAL(KIND=8),INTENT(IN)::CORX(NODN), CORY(NODN), CORZ(NODN)      
      REAL(KIND=8),INTENT(INOUT)::UX(NODN), UY(NODN), UZ(NODN)

      INTEGER(KIND=4)::INOD, NUMNEI, IJ, INEI
      REAL(KIND=8)::X0, Y0, Z0, U_N, U_M, U_S, U_X, U_Y, U_Z


      DO INOD = NODEID(-3)+1, NODEID(-4)

        IF(NWALLID(INOD,2).EQ.-10) CYCLE

        NUMNEI = NLINK(INOD)%I(0)
        X0 = CORX(INOD)
        Y0 = CORY(INOD)
        Z0 = CORZ(INOD)        

        DO IJ = 1, NUMNEI
          INEI = NLINK(INOD)%I(IJ)
          IF(NODEID(INEI).LT.0) CYCLE !CONSIDER ONLY NON-GHOST NEI
          IF(NODEID(INEI).EQ.2) CYCLE !SKIP BOTTOM NEI NODES          
          IF(NWALLID(INEI,2).EQ.-10) CYCLE

          ! THE NEI THAT WILL NOW COME IS THE CLOSEST TO INOD
          ! BECAUSE NLINK IS SORTED IS ASCENDING ORDER                
          
          EXIT
        ENDDO

        U_X = UX(INEI)
        U_Y = UY(INEI)
        U_Z = UZ(INEI)

        U_N = 0D0
        U_M = U_X*SMX(INOD) + U_Y*SMY(INOD) + U_Z*SMZ(INOD)
        U_S = U_X*SSX(INOD) + U_Y*SSY(INOD) + U_Z*SSZ(INOD)
        
        UX(INOD) = U_N*SNX(INOD) + U_M*SMX(INOD) + U_S*SSX(INOD)
        UY(INOD) = U_N*SNY(INOD) + U_M*SMY(INOD) + U_S*SSY(INOD)
        UZ(INOD) = U_N*SNZ(INOD) + U_M*SMZ(INOD) + U_S*SSZ(INOD)

      ENDDO

      END SUBROUTINE BOTSLIPBC

!!---------------------------END BOTSLIPBC----------------------------!!



!!-----------------------------UPDATE_CO------------------------------!!
      SUBROUTINE UPDATE_CO(CSUXT1,CSUYT1,CSUZT1)
      USE COMMONMOD
      USE MLPGKINE
      IMPLICIT NONE
      !INCLUDE 'COMMON.F'

      REAL(KIND=8),INTENT(INOUT)::CSUXT1(LNODE,1),CSUYT1(LNODE,1)
      REAL(KIND=8),INTENT(INOUT)::CSUZT1(LNODE,1)

      INTEGER(KIND=4)::INOD,NODTAL

      NODTAL=NODEID(-1)  !ONLY THE WATER PARTICLE UPDATE
      COORX(:,3)=COORX(:,2)
      COORY(:,3)=COORY(:,2)
      COORZ(:,3)=COORZ(:,2)
      DO INOD=1,NODTAL
        ! if (nodeid(inod).eq.8) then
        ! COORX(INOD,1)=COORX(INOD,2)
        ! COORY(INOD,1)=COORY(INOD,2)+UY(INOD,1)*DT  
        ! COORZ(INOD,1)=COORZ(INOD,2)+UZ(INOD,1)*DT  
        ! COORX(INOD,1)=COORX(INOD,1)+0.5*(UX(INOD,1)+UX(INOD,4))*DT !
        ! COORY(INOD,1)=COORY(INOD,1)+0.5*(UY(INOD,1)+UY(INOD,4))*DT !
        ! COORZ(INOD,1)=COORZ(INOD,1)+0.5*(UZ(INOD,1)+UZ(INOD,4))*DT !
        ! COORX(INOD,1)=COORX(INOD,1)+UX(INOD,1)*DT !0.5*(UX(INOD,1)+UX(INOD,4))*DT !
        ! COORY(INOD,1)=COORY(INOD,1)+UY(INOD,1)*DT !0.5*(UY(INOD,1)+UY(INOD,4))*DT !
        ! COORZ(INOD,1)=COORZ(INOD,1)+UZ(INOD,1)*DT !0.5*(UZ(INOD,1)+UZ(INOD,4))*DT !
        ! else
      
        ! Shagun Edit CS
        COORX(INOD,1)=COORX(INOD,2)+UX(INOD,1)*DT  
        COORY(INOD,1)=COORY(INOD,2)+UY(INOD,1)*DT  
        COORZ(INOD,1)=COORZ(INOD,2)+UZ(INOD,1)*DT  
        !COORX(INOD,1)=COORX(INOD,2)+(CSUXT1(INOD,1)+UX(INOD,1))/2D0*DT  
        !COORY(INOD,1)=COORY(INOD,2)+(CSUYT1(INOD,1)+UY(INOD,1))/2D0*DT  
        !COORZ(INOD,1)=COORZ(INOD,2)+(CSUZT1(INOD,1)+UZ(INOD,1))/2D0*DT        

        ! end if
      ENDDO
      
      ! Shagun edit
      CSUXT1(:,1)=UX(:,1)
      CSUYT1(:,1)=UY(:,1)
      CSUZT1(:,1)=UZ(:,1)
      
      END SUBROUTINE UPDATE_CO
!!-----------------------------UPDATE_CO------------------------------!!



!!------------------------------NEWCOOR4------------------------------!!
      SUBROUTINE NEWCOOR4(BNDNP,BNDXY,BNDFIX,FSNOD1,FSNOD2,DDR,
     +  DOMX,DOMY,DOMZ)
      USE COMMONMOD
      USE MLPGKINE
      USE NEIGHNODES
      IMPLICIT NONE
      !INCLUDE 'COMMON.F'

        INTEGER(KIND=4),INTENT(IN)::FSNOD1,FSNOD2,BNDNP
        INTEGER(KIND=4),INTENT(IN)::BNDFIX(3,0:BNDNP)
        REAL(KIND=8),INTENT(IN)::BNDXY(BNDNP,3)
        REAL(KIND=8),INTENT(IN)::DDR(NODEID(0))
        REAL(KIND=8),INTENT(IN)::DOMX(2),DOMY(2),DOMZ(2)

        INTEGER(KIND=4)::IK,II,IJ,IMINI,IMINJ,ICOUNT
        INTEGER(KIND=4)::I,NNEI2,I2,J
        REAL(KIND=8)::TMPR1,TMPR7,TMPR4,TMPR5,TMPR6,RIAV
        REAL(KIND=8)::DRMIN,BOTXLIM

        WRITE(8,*)'[MSG] ENTERING NEWCOOR4'
      
        DRMIN=0.0001D0      
        BOTXLIM=5D0

        !! WAVEMAKER READJUST
        IK=FSNOD2-NODEID(-2)
        TMPR1=COORZ(FSNOD2,1)/BNDXY(IK,3)
        DO II=NODEID(-2)+1,NODEID(-3)
          IK=II-NODEID(-2)
          COORZ(II,1)=TMPR1*BNDXY(IK,3)
        ENDDO

        !! SET Z=0 for BOTTOM AND SIDEWALLS
        COORZ(BNDFIX(1,1:BNDFIX(1,0)),1)=0D0        
        COORZ(BNDFIX(2,1:BNDFIX(2,0)),1)=0D0        
        COORZ(BNDFIX(3,1:BNDFIX(3,0)),1)=0D0        

        !! SIDEWALL X-ADJUST
        IK=FSNOD1-NODEID(-2)
        TMPR7=(DOMX(2)-COORX(FSNOD1,1))/(DOMX(2)-BNDXY(IK,1))
        DO II=1,BNDFIX(2,0)
          IJ=BNDFIX(2,II)
          IK=IJ-NODEID(-2)
          COORX(IJ,1)=DOMX(2)-TMPR7*(DOMX(2)-BNDXY(IK,1))
        ENDDO
        DO II=1,BNDFIX(3,0)
          IJ=BNDFIX(3,II)
          IK=IJ-NODEID(-2)
          COORX(IJ,1)=DOMX(2)-TMPR7*(DOMX(2)-BNDXY(IK,1))
        ENDDO

        !! BOTTOM NODE XY-ADJUST
        IK=FSNOD1-NODEID(-2)
        TMPR7=(BOTXLIM-COORX(FSNOD1,1))/(BOTXLIM-BNDXY(IK,1))
        DO II=1,BNDFIX(1,0)
          IJ=BNDFIX(1,II)
          IK=IJ-NODEID(-2)

          COORY(IJ,1)=BNDXY(IK,2)

          IF(BNDXY(IK,1).LT.BOTXLIM)THEN
            COORX(IJ,1)=BOTXLIM-TMPR7*(BOTXLIM-BNDXY(IK,1))
          ELSE
            COORX(IJ,1)=BNDXY(IK,1)
          ENDIF
        ENDDO

        ! Detecting if two nodes are too close      
        IMINI=1; IMINJ=1
        TMPR1=10000D0      
        !WRITE(8,*)'[TCL] Too Close', RIAV
        ICOUNT=0
        DO I=1,NODEID(-1)               

          NNEI2=NLINK(I)%I(0)        
          RIAV=0.2*DDR(I)

          DO I2=1,NNEI2
            J=NLINK(I)%I(I2)
            IF(I.EQ.J) CYCLE          

            TMPR4=COORX(J,1)-COORX(I,1)
            TMPR5=COORY(J,1)-COORY(I,1)
            TMPR6=COORZ(J,1)-COORZ(I,1)
            TMPR7=DSQRT(TMPR4**2 + TMPR5**2 + TMPR6**2)
            
            IF(TMPR7.LT.TMPR1)THEN
              TMPR1=TMPR7
              IMINI=I
              IMINJ=J
            ENDIF

            IF(TMPR7.LE.RIAV)THEN
              ICOUNT=ICOUNT+1   
              WRITE(8,'(A6,2I10,2I5,2F15.6)')'[TCL]',I,J,NODEID(I),
     +          NODEID(J),TMPR7,DDR(I)           
            ENDIF                     

          ENDDO      

        ENDDO
        WRITE(8,*)'[TCL] Total',ICOUNT
        WRITE(8,'(A8,F15.6,4I10)')'[MINR] ',TMPR1,IMINI,IMINJ,
     +    NODEID(IMINI),NODEID(IMINJ)

      END SUBROUTINE NEWCOOR4
!!------------------------------NEWCOOR4------------------------------!!



!!------------------------------NEWCOOR5------------------------------!!
      SUBROUTINE NEWCOOR5(BNDNP,BNDXY,BNDFIX,FSNOD1,FSNOD2,DDR,
     +  DOMX,DOMY,DOMZ)
      USE COMMONMOD
      USE MLPGKINE
      USE NEIGHNODES
      IMPLICIT NONE
      !INCLUDE 'COMMON.F'

        INTEGER(KIND=4),INTENT(IN)::FSNOD1,FSNOD2,BNDNP
        INTEGER(KIND=4),INTENT(IN)::BNDFIX(3,0:BNDNP)
        REAL(KIND=8),INTENT(IN)::BNDXY(BNDNP,3)
        REAL(KIND=8),INTENT(IN)::DDR(NODEID(0))
        REAL(KIND=8),INTENT(IN)::DOMX(2),DOMY(2),DOMZ(2)

        INTEGER(KIND=4)::IK,II,IJ,IMINI,IMINJ,ICOUNT
        INTEGER(KIND=4)::I,NNEI2,I2,J
        REAL(KIND=8)::TMPR1,TMPR7,TMPR4,TMPR5,TMPR6,RIAV
        REAL(KIND=8)::DRMIN,BOTXLIM

        WRITE(8,*)'[MSG] ENTERING NEWCOOR4'
      
        DRMIN=0.0001D0      
        BOTXLIM=5D0

        !! WAVEMAKER READJUST
        IK=FSNOD2-NODEID(-2)
        TMPR1=COORZ(FSNOD2,1)/BNDXY(IK,3)
        DO II=NODEID(-2)+1,NODEID(-3)
          IK=II-NODEID(-2)
          COORZ(II,1)=TMPR1*BNDXY(IK,3)
        ENDDO

        !! FIX BOTTOM NODES
        DO II = NODEID(-3)+1, NODEID(-4)
          IK = II - NODEID(-2)
          COORX(II,1) = BNDXY(IK,1)
          COORY(II,1) = BNDXY(IK,2)
          COORZ(II,1) = BNDXY(IK,3)
        ENDDO      

        ! Detecting if two nodes are too close      
        IMINI=1; IMINJ=1
        TMPR1=10000D0      
        !WRITE(8,*)'[TCL] Too Close', RIAV
        ICOUNT=0
        DO I=1,NODEID(-1)               

          NNEI2=NLINK(I)%I(0)        
          RIAV=0.2*DDR(I)

          DO I2=1,NNEI2
            J=NLINK(I)%I(I2)
            IF(I.EQ.J) CYCLE          

            TMPR4=COORX(J,1)-COORX(I,1)
            TMPR5=COORY(J,1)-COORY(I,1)
            TMPR6=COORZ(J,1)-COORZ(I,1)
            TMPR7=DSQRT(TMPR4**2 + TMPR5**2 + TMPR6**2)
            
            IF(TMPR7.LT.TMPR1)THEN
              TMPR1=TMPR7
              IMINI=I
              IMINJ=J
            ENDIF

            IF(TMPR7.LE.RIAV)THEN
              ICOUNT=ICOUNT+1   
              WRITE(8,'(A6,2I10,2I5,2F15.6)')'[TCL]',I,J,NODEID(I),
     +          NODEID(J),TMPR7,DDR(I)           
            ENDIF                     

          ENDDO      

        ENDDO
        WRITE(8,*)'[TCL] Total',ICOUNT
        WRITE(8,'(A8,F15.6,4I10)')'[MINR] ',TMPR1,IMINI,IMINJ,
     +    NODEID(IMINI),NODEID(IMINJ)

      END SUBROUTINE NEWCOOR5
!!------------------------------NEWCOOR5------------------------------!!      



!!--------------------------GIVENINITIALV_P---------------------------!!
      SUBROUTINE GIVENINITIALV_P(NODN,PTMP)
      USE COMMONMOD
      USE MLPGKINE
      IMPLICIT NONE
      !INCLUDE 'COMMON.F'
      
      INTEGER(KIND=4),INTENT(IN)::NODN
      REAL(KIND=8),INTENT(OUT)::PTMP(NODN)

      INTEGER(KIND=4)::I
      REAL(KIND=8)::ZZ

      ! [Note]: Ensure NODN = NODEID(0)

      ! THE VELOCITY AT START IS ZERO
      ! THE DENSITY IS NORMAL WATER DENSITY TAKEN AS 1.
      ! THESE MAY BE CHANGED LATER

      ROU(1:NODN)=1000.d0

      DO I=1,NODN

        ! DENSITY 1000.0 KG/M3  
        ROU(I)=1000.D0        
        ! VELOCITY
        UX(I,1)=0.D0
        UY(I,1)=0.D0
        UZ(I,1)=0.D0
        ! THE STATIC PRESSURE: YY=0 CORRESPONDING TO SEA BED
        ZZ=COORZ(I,1)
        PTMP(I)=(H0-ZZ)*ROU(I)*(-GRA)
      ENDDO
      
      END SUBROUTINE GIVENINITIALV_P
!!------------------------END GIVENINITIALV_P-------------------------!!



*****************************************************************

*************
* ======================================================================
* NIST GUIDE TO AVAILABLE MATH SOFTWARE.
* FULLSOURCE FOR MODULE SORT2 FROM PACKAGE NAPACK.
* RETRIEVED FROM NETLIB ON TUE MAR 28 01:36:29 2000.
* ======================================================================
C
C      ________________________________________________________
C     |                                                        |
C     |            SORT AN ARRAY IN INCREASING ORDER           |
C     |                                                        |
C     |    INPUT:                                              |
C     |                                                        |
C     |         X     --ARRAY OF NUMBERS                       |
C     |                                                        |
C     |         W     --WORKING ARRAY (LENGTH  AT LEAST N)     |
C     |                                                        |
C     |         N     --NUMBER OF ARRAY ELEMENTS TO SORT       |
C     |                                                        |
C     |    OUTPUT:                                             |
C     |                                                        |
C     |         X     --ORIGINAL ARRAY                         |
C     |                                                        |
C     |         Y     --INDICES OF X GIVING INCREASING ORDER   |
C     |________________________________________________________|
C
      SUBROUTINE SORT2(X,Y,DW,N)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER Y(N),DW(N)
      DIMENSION X(N)
      INTEGER I,J,K,L,M,N,P,Q
C
      I = 1
10    K = I
20    J = I
      Y(I) = I
      I = I + 1
      IF ( J .EQ. N ) GOTO 30
      IF ( X(I) .GE. X(J) ) GOTO 20
      DW(K) = I
      GOTO 10
30    IF ( K .EQ. 1 ) RETURN
      DW(K) = N + 1
40    M = 1
      L = 1
50    I = L
      IF ( I .GT. N ) GOTO 120
      P = Y(I)
      S = X(P)
      J = DW(I)
      K = J
      IF ( J .GT. N ) GOTO 100
      Q = Y(J)
      T = X(Q)
      L = DW(J)
      Y(I) = L
60    IF ( S .GT. T ) GOTO 70
      DW(M) = P
      M = M + 1
      I = I + 1
      IF ( I .EQ. K ) GOTO 80
      P = Y(I)
      S = X(P)
      GOTO 60
70    DW(M)= Q
      M = M + 1
      J = J + 1
      IF ( J .EQ. L ) GOTO 110
      Q = Y(J)
      T = X(Q)
      GOTO 60
80    DW(M) = Q
      K = M + L - J
      I = J - M
90    M = M + 1
      IF ( M .EQ. K ) GOTO 50
      DW(M) = Y(M+I)
      GOTO 90
100   Y(I) = J
      L = J
110   DW(M) = P
      K = M + K - I
      I = I - M
      GOTO 90
120   I = 1
130   K = I
      J = Y(I)
140   Y(I) = DW(I)
      I = I + 1
      IF ( I .LT. J ) GOTO 140
      DW(K) = I
      IF ( I .LE. N ) GOTO 130
      IF ( K .GT. 1 ) GOTO 40
      RETURN
      END
C
**** END OF SORT2


!!------------------------------NODEGRID------------------------------!!      
      SUBROUTINE NODEGRID(MESHFILE, FSNOD1, FSNOD2,
     +  DOMX, DOMY, DOMZ, CYLX, CYLY, CYLR)      
      USE COMMONMOD
      IMPLICIT NONE
      !INCLUDE 'COMMON.F'

      ! COOR:  COORDINATE OF NODES
      ! NODEID: ID NUMBER OF NODES
      ! NONB: NODE ON BOUNDARY
      ! NODEID(0): TOTAL NUMBER OF NODES

      REAL(KIND=8),INTENT(IN)::DOMX(2),DOMY(2),DOMZ(2),CYLX,CYLY,CYLR
      CHARACTER(LEN=256),INTENT(IN)::MESHFILE      
      INTEGER(KIND=4),INTENT(OUT)::FSNOD1,FSNOD2      

      INTEGER(KIND=4)::I,J,NGHST,IEND
      !REAL(KIND=8)::

C       CALL CYLIND3(MESHFILE, FSNOD1, FSNOD2, DOMX, DOMY, DOMZ,
C      +  CYLX, CYLY, CYLR)
C       CALL CYLIND4(MESHFILE, FSNOD1, FSNOD2, DOMX, DOMY, DOMZ,
C      +  CYLX, CYLY, CYLR)
      CALL CYLIND5(MESHFILE, FSNOD1, FSNOD2, DOMX, DOMY, DOMZ,
     +  CYLX, CYLY, CYLR)
*
      WRITE(8,*)'TOTAL WATER PARTICLE NUMBER IS,', NODEID(-2)
      WRITE(8,*)'TOTAL WATERANDINNER WALL PARTICLE NUMBER',NODEID(-1)
      WRITE(8,*)'TOTAL  PARTICLE NUMBER IS,', NODEID(0)
      IF (NODEID(0).GT.LNODE) THEN
        PRINT*,'ERROR IN ALLOCATION OF NUMBER OF NODES'
        STOP
      END IF

      !	Finding unique boundary nodes
      do i = nodeid(-2)+1,nodeid(0)
        do j = nodeid(-2)+1,nodeid(0)
          if (i.eq.j) CYCLE
          if ((coorX(i,1).eq.coorX(j,1)).and.
     +      (coorY(i,1).eq.coorY(j,1)).AND.
     +      (coorZ(i,1).eq.coorZ(j,1))) then
            print*,'unique node found',i,j,nodeid(i),nodeid(j)
            print*,coorX(i,1),coorY(i,1),coorZ(i,1)
            print*,coorX(J,1),coorY(J,1),coorZ(J,1)
            !pause
          end if
        end do
      end do

      ICELLX=INT((IXMAX)/RCELL)+1
      ICELLY=INT((IYMAX)/RCELL)+1
      ICELLZ=INT((IZMAX*2)/RCELL)+1
      
      WRITE(8,*)'CELL SIZES:',ICELLX,ICELLY,ICELLZ

      NGHST=0
      DO I=NODEID(-1)+1,NODEID(0)
        IF(NODEID(I).EQ.-9)NGHST=NGHST+1
      ENDDO
      OPEN(130,FILE='Output/XY_FLUENT.DAT')
      WRITE(130,'(3I8,1F10.4,3I8)') NODEID(-1)+NGHST
      WRITE(130,*)'TIME=',0D0
      IEND=NODEID(-1)
      DO I=1,IEND
        IF(ABS(COORX(I,1)).LE.70.AND.ABS(COORY(I,2))
     +             .LE.70.AND.ABS(COORZ(I,2)).LE.70)THEN
      
          WRITE(130,'(7E20.8,1I4)')COORX(I,1),COORY(I,1),COORZ(I,1),
     +      0D0,SNX(I),SNY(I),SNZ(I),NODEID(I)
        ELSE
          WRITE(130,'(7F16.8,1I4)') -10.,-10.,-10.,0,0,0,0,NODEID(I)
                              
        ENDIF      
      ENDDO

      DO I=NODEID(-1)+1,NODEID(0)
        IF(NODEID(I).EQ.-9)THEN
          WRITE(130,'(7E20.8,1I4)')COORX(I,1),COORY(I,1),COORZ(I,1),
     +                    0D0,0D0,0D0,0D0,NODEID(I)
        ENDIF
      ENDDO
      WRITE(130,*)'****' 
      CLOSE(130)      

      END SUBROUTINE NODEGRID
!!----------------------------END NODEGRID----------------------------!!

 

!!----------------------------MINPRESSURE-----------------------------!!
      SUBROUTINE MINPRESSURE(JET1,PPMIN1,JET2)
      !INCLUDE 'COMMON.F'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      !WE ONLY NEED JET2 THE MINIMUM VALUE, 
      DIMENSION PPMIN1(JET1)
      DO I=1, JET2 
       DO J=I+1,JET1
         IF(PPMIN1(J).LT.PPMIN1(I))THEN
           PMAX=PPMIN1(I)
           PPMIN1(I)=PPMIN1(J)
           PPMIN1(J)=PMAX
         ENDIF
       ENDDO
      ENDDO

      RETURN
      END
!!--------------------------END MINPRESSURE---------------------------!!



!!----------------------------MAXPRESSURE-----------------------------!!
      SUBROUTINE MAXPRESSURE(JET1,PPMIN1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      !ONLY GET THE MAXIMUM VALUE IN THE ARRAY,SO ONLY NEED ONE LOOP
      DIMENSION PPMIN1(JET1)

      DO I=1, 1 
       DO J=I+1,JET1
         IF(PPMIN1(J).GT.PPMIN1(I))THEN
           PMAX=PPMIN1(J)
           PPMIN1(J)=PPMIN1(I)
           PPMIN1(I)=PMAX
         ENDIF
       ENDDO
      ENDDO

      RETURN
      END
!!--------------------------END MAXPRESSURE---------------------------!! 



!!------------------------JUDGEFREESURFACE_SHA------------------------!!
      SUBROUTINE JUDGEFREESURFACE_SHA(DDR,FSSRCH,SPONGEX)
      USE COMMONMOD
      USE MLPGKINE
      USE NEIGHNODES
      IMPLICIT NONE
      !INCLUDE 'COMMON.F'

      REAL(KIND=8),INTENT(IN)::DDR(NODEID(0)),SPONGEX
      LOGICAL,INTENT(IN)::FSSRCH(NODEID(0))

      INTEGER(KIND=4)::NNEI,NNEI2,N0TO4,N4TO0
      !INTEGER(KIND=4)::N0NUM,N4NUM,N0MIS,N4MIS,N0MIS2
      INTEGER(KIND=4)::IE_N,IE_S,IE_W,IE_E,IE_UP,IE_DOWN,IE_RES
      INTEGER(KIND=4)::I,I2,J,INOD
      REAL(KIND=8)::RIAV,RMEAN,VLIM      
      !REAL(KIND=8)::N4MIN,N4MAX,N4AVG
      !REAL(KIND=8)::N0MIN,N0MAX,N0AVG
      REAL(KIND=8)::DDR0P7,DDR1P5      
      REAL(KIND=8)::TMPR1,TMPR2,TMPR3,TMPR4,TMPR5,TMPR6,TMPR7

      !N0MIN=1D8
      !N0MAX=0D0
      !N0NUM=0
      !N0MIS=0      
      !N0MIS2=0      
      !N4MIN=1D8
      !N4MAX=0D0
      !N4NUM=0
      !N4MIS=0
      N0TO4=0
      N4TO0=0
      VLIM=0.23d0 !0.225 !0.24d0 !0.25d0
      DO I=1,NODEID(-2)
      !DO I=1,0
        IF(COORX(I,2).GE.SPONGEX)CYCLE  
        !IF(COORZ(I,2).LT.0.85)CYCLE        

        NNEI2=NLINK(I)%I(0)

        TMPR1=0D0
        TMPR2=0D0
        TMPR3=0D0
        NNEI=0
        RMEAN=0D0        
        RIAV=2D0*DDR(I)
        DO I2=1,NNEI2
          J=NLINK(I)%I(I2)
          IF(I.EQ.J) CYCLE          
          !IF(NODEID(J).LT.0)CYCLE
          IF(.NOT.(FSSRCH(J))) CYCLE

          TMPR4=COORX(J,2)-COORX(I,2)
          TMPR5=COORY(J,2)-COORY(I,2)
          TMPR6=COORZ(J,2)-COORZ(I,2)
          TMPR7=DSQRT(TMPR4**2 + TMPR5**2 + TMPR6**2)

          IF(TMPR7.LE.RIAV)THEN
            NNEI=NNEI+1
            RMEAN=RMEAN+TMPR7
C             TMPR1=TMPR1+(TMPR4/(TMPR7**3))
C             TMPR2=TMPR2+(TMPR5/(TMPR7**3))
C             TMPR3=TMPR3+(TMPR6/(TMPR7**3))                        
            TMPR1=TMPR1+(TMPR4/(TMPR7))
            TMPR2=TMPR2+(TMPR5/(TMPR7))
            TMPR3=TMPR3+(TMPR6/(TMPR7))                        
          ENDIF
        ENDDO

        IF(NNEI.EQ.3)THEN
          !WRITE(21,*)'[ERR] NNEI = 0 FOR JFS1 FOR NODE',I
          WRITE(8,*)'[FSE] FS DETECT NNEI.LT.3 FOR NODE',I,NNEI          
          TMPR7=2*VLIM
        ELSE
          TMPR1=(TMPR1)
          TMPR2=(TMPR2)
          TMPR3=(TMPR3)        
          TMPR7=DSQRT(TMPR1**2 + TMPR2**2 + TMPR3**2)/(NNEI)
        ENDIF

        !!---Code for checking FS and fluid misidentification---!!        
C         IF(NODEID(I).EQ.0)THEN
C           N0NUM=N0NUM+1
C           N0AVG=N0AVG+TMPR7
C           IF(TMPR7.LT.N0MIN) N0MIN=TMPR7
C           IF(TMPR7.GT.N0MAX) N0MAX=TMPR7
C           IF(TMPR7.GT.VLIM) THEN
C             N0MIS=N0MIS+1            
C             NODEID(I)=-2
C           ENDIF
C         ELSE
C           N4NUM=N4NUM+1
C           N4AVG=N4AVG+TMPR7
C           IF(TMPR7.LT.N4MIN) N4MIN=TMPR7
C           IF(TMPR7.GT.N4MAX) N4MAX=TMPR7
C           IF(TMPR7.LT.VLIM) THEN
C             N4MIS=N4MIS+1
C           ENDIF
C         ENDIF        
C       ENDDO

C       N0AVG=N0AVG/N0NUM
C       N4AVG=N4AVG/N4NUM      
      !!--------------------------------------------------------!!

        IF(NODEID(I).EQ.0)THEN          
          IF(TMPR7.GT.VLIM) THEN            
            NODEID(I)=-2
          ENDIF
        ELSE          
          IF(TMPR7.LT.VLIM) THEN
            NODEID(I)=0
            N4TO0=N4TO0+1
          ENDIF
        ENDIF        
      ENDDO    

      DO I=1,NODEID(-2)
        IF(NODEID(I).NE.-2)CYCLE                

        NNEI2=NLINK(I)%I(0)

        IE_N=0
        IE_W=0
        IE_E=0
        IE_S=0
        IE_UP=0
        IE_DOWN=0
                
        RIAV=2D0*DDR(I)
        DDR0P7=0.4D0*RIAV
        DDR1P5=1.0D0*RIAV
        DO I2=1,NNEI2
          J=NLINK(I)%I(I2)
          IF(I.EQ.J) CYCLE          
          !IF(NODEID(J).LT.0)CYCLE
          IF(.NOT.(FSSRCH(J))) CYCLE

          TMPR4=COORX(J,2)-COORX(I,2)
          TMPR5=COORY(J,2)-COORY(I,2)
          TMPR6=COORZ(J,2)-COORZ(I,2)
          TMPR7=DSQRT(TMPR4**2 + TMPR5**2 + TMPR6**2)

          IF(TMPR7.GT.RIAV)CYCLE
          IF((ABS(TMPR4).LE.DDR0P7).AND.(ABS(TMPR5).LE.DDR0P7).AND.
     +      (ABS(TMPR6).LE.DDR1P5))THEN
            IF(TMPR6.GT.0D0)THEN
              IE_UP=IE_UP+1
            ELSE
              IE_DOWN=IE_DOWN+1
            ENDIF
          ENDIF
          IF((ABS(TMPR4).LE.DDR1P5).AND.(ABS(TMPR5).LE.DDR0P7).AND.
     +      (ABS(TMPR6).LE.DDR0P7))THEN
            IF(TMPR4.GT.0D0)THEN
              IE_E=IE_E+1
            ELSE
              IE_W=IE_W+1
            ENDIF
          ENDIF
          IF((ABS(TMPR4).LE.DDR0P7).AND.(ABS(TMPR5).LE.DDR1P5).AND.
     +      (ABS(TMPR6).LE.DDR0P7))THEN
            IF(TMPR5.GT.0D0)THEN
              IE_N=IE_N+1
            ELSE
              IE_S=IE_S+1
            ENDIF
          ENDIF
        ENDDO

        NODEID(I)=0
        IE_RES=IE_E*IE_W*IE_N*IE_S*IE_UP*IE_DOWN
        IF(IE_RES.EQ.0)THEN
          !N0MIS2=N0MIS2+1
          NODEID(I)=4
          N0TO4=N0TO4+1
        ENDIF
      ENDDO

      ! WRITE(21,'("0  ",I10,3F15.6)')N0NUM,N0AVG,N0MIN,N0MAX
      ! WRITE(21,'("4  ",I10,3F15.6)')N4NUM,N4AVG,N4MIN,N4MAX
      ! WRITE(21,'("0M ",I10,2F15.6)')N0MIS,N0ZAVG,VLIM
      !WRITE(21,'(7F15.6,3I10)')TOTAL_TIME,N0MIN,N0MAX,N0AVG,N4MIN,
      !&  N4MAX,N4AVG,N0MIS,N4MIS,N0MIS2
      WRITE(8,'(" [FSD] ",F15.6,2I10)')TOTAL_TIME,N0TO4,N4TO0

      DO INOD=NODEID(-1)+1,NODEID(0) !SET THE ZERO
        NWALLID(INOD,2)=-10
      ENDDO       

      END SUBROUTINE JUDGEFREESURFACE_SHA
!!----------------------END JUDGEFREESURFACE_SHA----------------------!!



!!----------------------------JUDGEBOTTOM-----------------------------!!
      SUBROUTINE JUDGEBOTTOM(LNODE, NODN, NODEID, NWALLID, 
     +  CORX, CORY, CORZ, DDR, FSSRCH)      
      USE NEIGHNODES
      IMPLICIT NONE
      !INCLUDE 'COMMON.F'

      !NODN SHOULD BE NODEID(0)
      INTEGER(KIND=4),INTENT(IN)::LNODE, NODN, NODEID(-7:NODN)
      INTEGER(KIND=4),INTENT(INOUT)::NWALLID(LNODE,4)
      REAL(KIND=8),INTENT(IN)::DDR(NODN)
      REAL(KIND=8),INTENT(IN)::CORX(NODN), CORY(NODN), CORZ(NODN)
      LOGICAL,INTENT(INOUT)::FSSRCH(NODN)

      INTEGER(KIND=4)::INOD, NUMNEI, IJ, INEI
      REAL(KIND=8)::TMPR1, DR2, X0, Y0, Z0, RIAV2 

      
      DO INOD = NODEID(-3)+1, NODEID(-4)
        NUMNEI = NLINK(INOD)%I(0)
        X0 = CORX(INOD)
        Y0 = CORY(INOD)
        Z0 = CORZ(INOD)  
        RIAV2 = DDR(INOD)**2

        DO IJ = 1, NUMNEI
          INEI = NLINK(INOD)%I(IJ)
          IF(NODEID(INEI).LT.0) CYCLE !CONSIDER ONLY NON-GHOST NEI
          IF(NODEID(INEI).EQ.2) CYCLE !SKIP BOTTOM NEI NODES
          IF(NODEID(INEI).EQ.8) CYCLE !SKIP THE LEFT PLACE TOO

          ! THE NEI THAT WILL NOW COME IS THE CLOSEST TO INOD
          ! BECAUSE NLINK IS SORTED IS ASCENDING ORDER        

          DR2 = (CORX(INEI)-X0)**2 + (CORY(INEI)-Y0)**2  
     +      + (CORZ(INEI)-Z0)**2

          IF(DR2.GT.RIAV2)THEN !DRY
            NWALLID(INOD,2) = -10
            FSSRCH(INOD) = .FALSE.
            WRITE(1617, '(I10,3F15.6)')INOD, X0, Y0, Z0
          ELSE !WET
            NWALLID(INOD,2) = 1
            FSSRCH(INOD) = .TRUE.
          ENDIF
          
          EXIT

        ENDDO

      ENDDO      


      WRITE(1617,'("---interface---")')

      
      DO INOD = NODEID(-3)+1, NODEID(-4)

        IF(NWALLID(INOD,2).EQ.-10) CYCLE !ONLY WET TYPE2 NODES

        NUMNEI = NLINK(INOD)%I(0)
        X0 = CORX(INOD)
        Y0 = CORY(INOD)
        Z0 = CORZ(INOD)  
        RIAV2 = DDR(INOD)**2
        NWALLID(INOD,3) = 0 !BY DEFAULT NO MLS INTERPOLATION        

        DO IJ = 1, NUMNEI
          INEI = NLINK(INOD)%I(IJ)          
          IF(NODEID(INEI).NE.2) CYCLE 

          ! ONLY TYPE2 NEIGHS

          DR2 = (CORX(INEI)-X0)**2 + (CORY(INEI)-Y0)**2  
     +      + (CORZ(INEI)-Z0)**2          

          IF(DR2.GT.RIAV2) CYCLE

          IF(NWALLID(INEI,2).EQ.-10)THEN
            NWALLID(INOD,3) = 9 !MLS INTERPOLATION            
            WRITE(1617, '(I10,3F15.6)')INOD, X0, Y0, Z0
            EXIT
          ENDIF            
        ENDDO

      ENDDO      

      END SUBROUTINE JUDGEBOTTOM

!!---------------------------END JUDGEBOTTOM--------------------------!!
      


!!-----------------------------LAPLACIAN------------------------------!!
      ! commented by shagun on 2021-04-07 
      ! to remove the use of PLAMDA and other such things

C       SUBROUTINE LAPLACIAN(INOD,DU,DV,DW)
C       USE COMMONMOD
C       USE MLPGKINE
C       USE NEIGHNODES
C       IMPLICIT NONE
C       !INCLUDE 'COMMON.F'

C       ! TO CALCULATE THE LAPLACIAN OPERATOR      

C       INTEGER(KIND=4),INTENT(IN)::INOD
C       REAL(KIND=8),INTENT(OUT)::DU,DV,DW

C       INTEGER(KIND=4)::NUMTOTAL,JJ,IFLUIDJ,KK
C       REAL(KIND=8)::WSUM1,UXWR,VYWR,WZWR,PNI,WEI,XX,YY,ZZ,ABSR      

C       WSUM1=0.D0
C       UXWR=0.D0
C       VYWR=0.D0
C       WZWR=0.D0
C       PNI=0.D0
C       NUMTOTAL=NLINK(INOD)%I(0)

C       DO 130 JJ=1,NUMTOTAL
C         WEI=0.D0
C         IFLUIDJ=NLINK(INOD)%I(JJ)
C         KK=NODEID(IFLUIDJ)

C         IF(INOD.EQ.IFLUIDJ)GOTO 130

C         XX=COORX(IFLUIDJ,1)-COORX(INOD,1)
C         YY=COORY(IFLUIDJ,1)-COORY(INOD,1)
C         ZZ=COORZ(IFLUIDJ,1)-COORZ(INOD,1)
C         ABSR=SQRT((XX**2.)+(YY**2.)+(ZZ**2.))  !**0.5

C         CALL WEIFUNC_PDN(ABSR,R2,WEI,DDL,PI)

C         IF(WEI.GT.1.0E-15.AND.NODEID(IFLUIDJ).LE.10)THEN
C           UXWR=UXWR+(UX(IFLUIDJ,1)-UX(INOD,1))*WEI
C           VYWR=VYWR+(UY(IFLUIDJ,1)-UY(INOD,1))*WEI
C           WZWR=WZWR+(UZ(IFLUIDJ,1)-UZ(INOD,1))*WEI
C         ENDIF
C 130	  CONTINUE
C       !-----------------------------------------*
      
C       ! commented by shagun on 2021-04-07 to remove PLAMDA influence
C       PNI=PLAMDA(INOD)
C       IF(PNI.LE.1.D-5)THEN
C         DU=0.
C         DV=0.
C         DW=0.
C       ELSE
C         DU=2*3*UXWR/PNI
C         DV=2*3*VYWR/PNI
C         DW=2*3*WZWR/PNI
C       ENDIF
C       !****************
C       RETURN
C       END SUBROUTINE LAPLACIAN
!!---------------------------END LAPLACIAN----------------------------!!



!!----------------------------GRADIENT_2P-----------------------------!!
      SUBROUTINE GRADIENT_2P(NODN,FB)
      USE COMMONMOD
      USE MLPGKINE
      IMPLICIT NONE
      !INCLUDE 'COMMON.F'

      INTEGER(KIND=4),INTENT(IN)::NODN
      REAL(KIND=8),INTENT(IN)::FB(NODN)
            
      INTEGER(KIND=4)::ND(0:NLMAXN),NDL(0:NLMAXN)
      INTEGER(KIND=4)::I,INOD,NN,IDWW,JN,JJ,NODJ
      REAL(KIND=8)::BJX(NLMAXN),BJY(NLMAXN)
      REAL(KIND=8)::FFF(NLMAXN),A(3,3),C(3,NLMAXN)      
      REAL(KIND=8)::AINV(3,3),AINPUT(3,3)
      REAL(KIND=8)::BJX1(NLMAXN),BJY1(NLMAXN),FF1(NLMAXN)
      REAL(KIND=8)::SUMW,SUMW1,PPXI,PPYI,XI,YI
      REAL(KIND=8)::RNIX,RNIY,RNXY,DIX,DIY,FXI,FXJ,FXI1
      REAL(KIND=8)::DIFNXY,DIFNXYAB

      
      IDWW=4  !USE THE FOURTH WEIGHT FUNCTION

      !initial value
      ND(0)=0
      NDL(0)=0
      DO I=1,NLMAXN
        ND(I)=0
        NDL(I)=0
        BJX(I)=0.D0
        BJY(I)=0.D0
        FFF(I)=0.D0
      ENDDO
      !	FXJB = 0.D0


      DO INOD=NODEID(-2)+1,NODEID(-1)
        PPXI=0.0
        PPYI=0.0
        ! if (nodeid(inod).eq.7.or.nodeid(inod).eq.1) then
        !	IF (NODEID(INOD).EQ.7.OR.nwallid(inod,4).EQ.7) GOTO 11	
        if (nodeid(inod).eq.8.or.nodeid(inod).eq.3) then
          XI=COORY(INOD,2)
          YI=COORZ(INOD,2)
        elseif (nodeid(inod).eq.7.or.nodeid(inod).eq.1) then
          XI=COORX(INOD,2)
          YI=COORZ(INOD,2)
        else
          goto 11
        end if
      

        NN=0
        RNIX=0
        RNIY=0
        RNXY=0
        
        CALL WEIGHTF_GRAD_2P(XI,YI,NN,ND,BJX,BJY,NDL,INOD,RNIX,
     +    RNIY,RNXY,NLMAXN,IDWW,1,0)

        IF (NN.LT.3) THEN
          ! WRITE(8,*)'WARNING: NO ENOUGH NODE NEAR THE POINT'
          ! WRITE(8,*)'CORRECTION TERM MADE ZERO,',INOD
          RNXY=0.D0
        ENDIF

        DIX=0.D0
        DIY=0.D0
        ! print*,nn
        ! pause	

        !MAY BE USED MINISUM ONE IN FUTURE OR ONE NEAR POINT (XI,YI)
        FXI=FB(INOD)  

        DO 10 JN=1,NN
          NODJ=ND(JN)
          FXJ=FB(NODJ)
          DIX=DIX+BJX(JN)*(FXJ-FXI)
          DIY=DIY+BJY(JN)*(FXJ-FXI)

10      CONTINUE

        DIFNXY=RNIX*RNIY-RNXY**2
        DIFNXYAB=ABS(DIFNXY)
        IF (DIFNXYAB.GT.0.0000000001) THEN
          PPXI=( DIX*RNIY-DIY*RNXY)/DIFNXY
          PPYI=(-DIX*RNXY+DIY*RNIX)/DIFNXY
        ELSE IF (RNIX.NE.0.0) THEN
          PPXI=DIX/RNIX
          PPYI=0.D0
        ELSE IF (RNIY.NE.0.0) THEN
          PPXI=0.D0
          PPYI=DIY/RNIY
        ELSE 
          ! WRITE(8,*)'WARNING ROTATING THE AXES NEEDED',INOD
        ENDIF

        if (nodeid(inod).eq.8.or.nodeid(inod).eq.3) then
          ppx(inod,1) = ppx(inod,1)
          PPY(INOD,1)=PPXI
          PPZ(INOD,1)=PPYI
        elseif (nodeid(inod).eq.7.or.nodeid(inod).eq.1) then
          ppx(inod,1) = PPXI
          PPY(INOD,1)= PPY(INOD,1)
          PPZ(INOD,1)=PPYI
        end if

        DIX = 0.
        DIY = 0.
        PPXI = 0.
        PPYI = 0.

        DO JJ=1,NN
          NODJ=ND(JJ)
          FFF(JJ)=FB(NODJ)
        ENDDO

        CALL MINPRESSURE(NN,FFF(1:NN),1)
        FXI1 = MIN(FFF(1),FXI)
        DO JN=1,NN
          NODJ=ND(JN)
          FXJ=FB(NODJ)
          DIX=DIX+BJX(JN)*(FXJ-FXI1)
          DIY=DIY+BJY(JN)*(FXJ-FXI1)
        END DO

        !	DIX1 = 0
        !	DIY1 = 0
        !	DO JN = 1,NN1
        !	FXJ = FF1(JN)
        !	FXI1 = FXJB(INOD)
        !	DIX1=DIX1+BJX1(JN)*((FXI+FXJ)-2*(FXI1))
        !	DIY1=DIY1+BJY1(JN)*((FXI+FXJ)-2*(FXI1))
        !	END DO

        ! ORIGINAL MPS	

        !	PPXI = (DIX*2)/SUMW
        !	PPYI = (DIY*2)/SUMW
        !	IF (NN1.NE.0) THEN
        !	PPXI = PPXI+(DIX1*2)/SUMW1
        !	PPYI = PPYI+(DIY1*2)/SUMW1
        !	END IF


        DIFNXY=RNIX*RNIY-RNXY**2
        DIFNXYAB=ABS(DIFNXY)
        IF (DIFNXYAB.GT.0.0000000001) THEN
          PPXI=( DIX*RNIY-DIY*RNXY)/DIFNXY
          PPYI=(-DIX*RNXY+DIY*RNIX)/DIFNXY
        ELSE IF (RNIX.NE.0.0) THEN
          PPXI=DIX/RNIX
          PPYI=0.D0
        ELSE IF (RNIY.NE.0.0) THEN
          PPXI=0.D0
          PPYI=DIY/RNIY
        ELSE 
          ! WRITE(8,*)'WARNING ROTATING THE AXES NEEDED',INOD
        ENDIF


        if (nodeid(inod).eq.8.or.nodeid(inod).eq.3) then
          ppx(inod,2) = ppx(inod,2)
          PPY(INOD,2)=PPXI
          PPZ(INOD,2)=PPYI
        elseif (nodeid(inod).eq.7.or.nodeid(inod).eq.1) then
          if (nwallid(inod,1).ne.5) then
            ppx(inod,2) = PPXI
            PPY(INOD,2)= PPY(INOD,1)
            PPZ(INOD,2)=PPYI
          else
            ppx(inod,2) = PPX(Inod,1)
            PPY(INOD,2)= PPY(INOD,1)
            PPZ(INOD,2)=PPz(inod,1)
          end if
        end if


11      CONTINUE
        !PPXI = 0.D0
        ! endif
      ENDDO
          
      END SUBROUTINE GRADIENT_2P
!!--------------------------END GRADIENT_2P---------------------------!!



!!--------------------------WEIGHTF_GRAD_2P---------------------------!!
      !SOLVE THE VELOCITY DIVERGENCE WILL USE FOLLOWING SUBROUTINE
      SUBROUTINE WEIGHTF_GRAD_2P(XQ,YQ,NN,ND,BJX,BJY,NDL,NOD,RNIX,
     +  RNIY,RNXY,NLM,IDW,ID_PV,IVIRT)
      USE COMMONMOD
      USE MLPGKINE
      USE NEIGHNODES
      IMPLICIT NONE
      !INCLUDE 'COMMON.F'

      INTEGER(KIND=4),INTENT(IN)::NLM,IDW,ID_PV,IVIRT,NOD
      REAL(KIND=8),INTENT(IN)::XQ,YQ
      INTEGER(KIND=4),INTENT(OUT)::NN,ND(0:NLM),NDL(0:NLM)
      REAL(KIND=8),INTENT(OUT)::BJX(NLM),BJY(NLM)
      REAL(KIND=8),INTENT(OUT)::RNIX,RNIY,RNXY
      
      INTEGER(KIND=4)::I,NLINI,IWALL,IN,NODIDII,NODIDIJ
      REAL(KIND=8)::XX(5),YY(5),RIAV,XDI,YDI,DI,RATIO,WWI

 
      !initial value
      ND(0)=0
      NDL(0)=0
      DO I=1,NLM
        ND(I)=0
        NDL(I)=0
        BJX(I)=0.D0
        BJY(I)=0.D0
      ENDDO
      ! GAUSS WEIGHT FUNCTION
      !	NNODE=NODEID(0)

      RNIX=0.D0   !SUM OF WEIGHT FUNCTION*(XJ-XI)**2
      RNIY=0.D0   !SUM OF WEIGHT FUNCTION*(YJ-YI)**2
      RNXY=0.D0   !SUM OF WEIGHT FUNCTION*(YJ-YI)*(XJ-XI)

      NN=0
      NLINI=NLINK(NOD)%I(0)
      IF (NLINI.LE.2) THEN
      ! WRITE(8,*)'LESS THAN 2 NEIGHBOUR NODES, NO DERIVATIVE CORRECTION '
        RETURN
      ENDIF
      ND(0:NLINI)=0

      IWALL = 0
      DO IN = 1,NLINI
        I=NLINK(NOD)%I(IN)
        IF (NODEID(I).NE.0. AND. NODEID(I).NE.4) IWALL = 1
      END DO

      ! XQ=COOR(NOD,1)
      ! YQ=COOR(NOD,2)

      IF (IDW .EQ. 4) THEN  !FOR OTHER WEIGTH FUNCTION
        DO 50 IN=1,NLINI  !IN NOT EQUAL TO NOD 

          I=NLINK(NOD)%I(IN)
          NODIDII=NODEID(NOD)
          NODIDIJ=NODEID(I)
          if (nodidii.eq.nodidij.or.nwallid(i,1).eq.7) then
            
            ! Shagun RIAV edit
            ! RIAV=(( R(I)+R(NOD) )/2)
            RIAV = R(NOD)
            ! Shagun edit

            ! RIAV = RIAV*0.8
            ! RIAV=R(I)
            ! RIAV=R(NOD)				!ADDED ON 17/08/2005
            ! IF(ID_PV.EQ.1)RIAV=1.6*DDL
            ! IF(ID_PV.EQ.2)
            ! RIAV=r3*DDL
            ! IF(ID_PV.EQ.3)RIAV=R3*DDL

            if (nodeid(nod).eq.8.or.nodeid(nod).eq.3) then
              XDI=XQ-COORY(I,2)
              YDI=YQ-COORZ(I,2)
            elseif (nodeid(nod).eq.7.or.nodeid(nod).eq.1) then
              XDI=XQ-COORX(I,2)
              YDI=YQ-COORZ(I,2)
            end if
             
            DI=DSQRT(XDI**2+YDI**2)
            
            RATIO=DI/RIAV
            IF (DI.LE.1.D-10) GOTO 50 !EXLUDING THE POINT ITSELF

            ! IF (ID_PV.EQ.1) THEN
              ! IF (IWALL.NE.0) THEN
                ! WWI=1-6*(RATIO)**2+8*(RATIO)**3-3*(RATIO)**4 
              ! ELSE
                ! WWI=1-6*(RATIO/2)**2+8*(RATIO/2)**3-3*(RATIO/2)**4 
              ! END IF
              ! CALL WEIFUNC(NOD,I,WWI,ID_PV) !PRESSURE,GRADIENT AND DIVERGENCE		
            ! ELSEIF(ID_PV.EQ.2) THEN
            WWI=1-6*(RATIO)**2+8*(RATIO)**3-3*(RATIO)**4 
            !	END IF

  

      ! SPECIAL DEAL WITH THE SOLID PARTILCE ON THE WAVEMAKER
            IF(NWALLID(I,2).EQ.-10) THEN
              WWI=-1.D0  !WALL SPECIAL PARTICLE
              RATIO = -1.D0
            END IF

            IF (IVIRT.EQ.1) NODIDIJ = ABS(NODEID(I))
            IF(NODIDIJ.GE.0.AND.NODIDIJ.LT.10)THEN
              IF (WWI.GT.1.E-15) THEN
              ! IF (RATIO.LT.1.AND.RATIO.GE.0) THEN

                NN=NN+1
                BJX(NN)=-WWI*XDI/DI/DI
                BJY(NN)=-WWI*YDI/DI/DI
                RNIX=RNIX+(XDI**2)*WWI/DI/DI   
                RNIY=RNIY+(YDI**2)*WWI/DI/DI     
                RNXY=RNXY+ XDI*YDI*WWI/DI/DI

                ND(NN)=I
              ENDIF   
            ENDIF
          end if

50      CONTINUE
        ND(0)=NN
        IF (NN.EQ.NLINI) THEN
        ! WRITE(8,*)'THE NUMBER OF NEIGHBOUR NODES IN 
        ! +     WEIGHTF_BJ MAY BE TOO SMALL, NOD IS ', NOD,' NUMBER: ',NN
        ENDIF
        ! LITTLE CHANGE TO THE ABOVE SENTENCE IS MADE ON 22/02/2004

        ! THE FOLLOWING SENTENCE IS ADDED ON 2/05/2004 TO 
        ! WARNING THE CASE WHERE THE NUMBER OF NUMBER NODES IS TOO FEW

        ! IF (NN.LT.3)WRITE(8,*)'THE NUMBER OF NEIGHBOUR NODES IN 
        ! +   WEIGHTF_BJ MAY BE TOO FEW, NOD IS:', NOD,'NUMBER=',NN

      ENDIF

      END SUBROUTINE WEIGHTF_GRAD_2P
!!------------------------END WEIGHTF_GRAD_2P-------------------------!!



!!----------------------------FIND_LAPLAC-----------------------------!!
      SUBROUTINE FIND_LAPLAC(NODN, UN)
      USE COMMONMOD
      USE MLPGKINE
      USE NEIGHNODES
      IMPLICIT NONE
      !INCLUDE 'COMMON.F'
      
      INTEGER(KIND=4),INTENT(IN)::NODN
      REAL(KIND=8),INTENT(OUT)::UN(NODN,1:3)

      INTEGER(KIND=4)::NTOT_WATER, NTOT_ALL, INOD
      INTEGER(KIND=4)::KI,KK,KK1,NN,I,II,IK,KI1
      REAL(KIND=8)::U_N,U_M,U_S

      NTOT_WATER=NODEID(-2)
      NTOT_ALL=NODEID(-1)
      DO INOD=NTOT_WATER+1,NTOT_ALL
        UN(INOD,1:3)=0.D0
      ENDDO

      DO 101 INOD=NTOT_WATER+1,NTOT_ALL
        KI=NLINK(INOD)%I(0)
        KK=NWALLID(INOD,1)  !COEFFICIENT ABOUT VELOCITY OF WALL PARTICLES
        KK1=NWALLID(INOD,2)
        IF(KK1.EQ.-10)THEN
          ! UN(INOD,1:3)=0.D0;
          UX(INOD,1:2) = 0; UY(INOD,1:2)=0.D0;
          UZ(INOD,1:2)=0.D0
          GOTO 101
        ENDIF

        NN=0
        DO II=1,KI
          KI1=NLINK(INOD)%I(II)
          KK1=NODEID(KI1)
          IF(KK1.EQ.0.OR.KK1.EQ.4)NN=NN+1
          
          IF(NN.EQ.1)THEN

            IF(KK.EQ.0)THEN       !VELOCITY IS ZERO
                UN(INOD,1)=UN(KI1,1); !UX(INOD,1)=0
                UN(INOD,2)=UN(KI1,2); !UY(INOD,1)=0
                UN(INOD,3)=UN(KI1,3); !UZ(INOD,1)=0
            ENDIF
            IF(KK.EQ.1)THEN !X-Y PLANE
      !			  UX(INOD,1:2)=UX(KI1,IK)
      !			  UY(INOD,1:2)=UY(KI1,IK)
      !			  UZ(INOD,1)  =0.D0
              UN(INOD,3)=UN(KI1,3)
            ENDIF

            IF(KK.EQ.2)THEN !Y-Z PLANE
      !			    UX(INOD,1)=0.D0
      !			  UY(INOD,1:2)=UY(KI1,IK)
      !			  UZ(INOD,1:2)=UZ(KI1,IK)
              UN(INOD,1)=UN(KI1,1)
            ENDIF
            
            IF(KK.EQ.3)THEN  !X-Z PLANE
      !			  UX(INOD,1:2)=UX(KI1,IK)
      !			    UY(INOD,1)=0.D0
      !			  UZ(INOD,1:2)=UZ(KI1,IK)
              UN(INOD,2)=UN(KI1,2)
            ENDIF

            IF(KK.EQ.4)THEN  !X-Y-Z PLANE
      !			  UX(INOD,2)=UX(KI1,2)
      !			  UY(INOD,2)=UY(KI1,2)
      !			  UZ(INOD,2)=UZ(KI1,2)
      !			  UX(INOD,1)=UX(KI1,1)
      !			  UY(INOD,1)=UY(KI1,1)
      !			  UZ(INOD,1)=UZ(KI1,1)

              UN(INOD,1)=UN(KI1,1)
              UN(INOD,2)=UN(KI1,2)
              UN(INOD,3)=UN(KI1,3)
            ENDIF
            
            IF(KK.EQ.5)THEN  !Y-Z VELOCITY IS 0,ONLY X VELOCITY
      !			  UX(INOD,1:2)=UX(KI1,IK)
      !			  UY(INOD,1:2)=0
      !			  UZ(INOD,1:2)=0
              UN(INOD,1)=UN(KI1,1)
              UN(INOD,2)=UN(KI1,2)
              UN(INOD,3)=UN(KI1,3)
            ENDIF

            IF(KK.EQ.6)THEN  !ONLY Y VELOCITY
      !			  UX(INOD,1:2)=0
      !			  UY(INOD,1:2)=UY(KI1,IK)
      !			  UZ(INOD,1:2)=0
              UN(INOD,1)=UN(KI1,1)
              UN(INOD,2)=UN(KI1,2)
              UN(INOD,3)=UN(KI1,3)
            ENDIF

            IF(KK.EQ.7)THEN  !ONLY Z VELOCITY
      !			  UX(INOD,1:2)=0
      !			  UY(INOD,1:2)=0
      !			  UZ(INOD,1:2)=UZ(KI1,IK)
              UN(INOD,1)=UN(KI1,1)
              UN(INOD,2)=UN(KI1,2)
              UN(INOD,3)=UN(KI1,3)
            ENDIF

            IF(KK.EQ.8)THEN  !ON THE SIDE WALL AND SLOPE
      ! 			  UX(INOD,2)=UX(KI1,2)
      !			  UZ(INOD,2)=UZ(KI1,2)
      !			  UX(INOD,1)=UX(KI1,1)
      !			  UZ(INOD,1)=UZ(KI1,1)

              UN(INOD,1)=UN(KI1,1)
              UN(INOD,2)=UN(KI1,2)
              UN(INOD,3)=UN(KI1,3)
            ENDIF
            GOTO 101
          ENDIF
        ENDDO

101	  CONTINUE

      DO I=NTOT_WATER+1,NTOT_ALL 
        !	IF (NODEID(I).EQ.7.OR.NWALLID(I,4).EQ.71) GOTO 111
        ! IF(I_CAL_V.EQ.1)CALL LAPLACIAN(I,DU,DV)  !CALCULATE THE VISCOUS TERM

        IF((NODEID(I).EQ.8).AND.(I_WM.EQ.15)) CYCLE

        UX(I,2)=UX(I,1)+UN(I,1) !-TANK_A(I,1)*DT
        UY(I,2)=UY(I,1)+UN(I,2) !-GRA*DT !-TANK_A(I,2)*DT
        UZ(I,2)= UZ(I,1)+ UN(I,3)
        IK = 2
        U_N = 0.D0
        U_M = UX(I,IK)*SMX(I)+UY(I,IK)*SMY(I)+
     +        UZ(I,IK)*SMZ(I)
        U_S = UX(I,IK)*SSX(I)+UY(I,IK)*SSY(I)+
     +        UZ(I,IK)*SSZ(I)
          
        UX(I,IK) = U_N*SNX(I)+U_M*SMX(I)+U_S*SSX(I)
        UY(I,IK) = U_N*SNY(I)+U_M*SMY(I)+U_S*SSY(I)
        UZ(I,IK) = U_N*SNZ(I)+U_M*SMZ(I)+U_S*SSZ(I)

        COORX(I,2)=COORX(I,1) !+ DT*(UX(I,2)-UX(I,1))  !UM(I)*DT !0.5*(UX(I,2)+UX(I,1))*DT   !X COORDINATE VELOCITY
        COORY(I,2)=COORY(I,1) !+ DT*(UY(I,2)-UY(I,1)) !VM(I)*DT !0.5*(UY(I,2)+UY(I,1))*DT   !Y CORDINATE VELOCITY
        COORZ(I,2)= COORZ(I,1)
        !	UM(I)=(VCOEFF*DU)*DT  !(VCOEFF*DU-TANK_A(I,1))*DT  !
        !	VM(I)=(VCOEFF*DV)*DT  !(VCOEFF*DV-TANK_A(I,2))*DT  !
        !	END IF
111   ENDDO

      ! ADJUST THE WALL PARTICLE VELOCITY FOR SLOPE ON FIND_ACC!
      END SUBROUTINE FIND_LAPLAC
!!--------------------------END FIND_LAPLAC---------------------------!!      



!!------------------------------FIND_ACCN-----------------------------!!
      SUBROUTINE FIND_ACCN(NODN,UM)
      USE COMMONMOD
      USE MLPGKINE
      IMPLICIT NONE
      !INCLUDE 'COMMON.F'

      INTEGER(KIND=4),INTENT(IN)::NODN
      REAL(KIND=8),INTENT(OUT)::UM(NODN)            
      INTEGER(KIND=4)::INOD,KK,KK1

      UM = 0.D0
      
      !!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(INOD,KK,KK1,ISTEP)
      !!$OMP DO SCHEDULE(DYNAMIC,100)
      DO 11 INOD=NODEID(-2)+1,NODEID(-1)
        KK=NWALLID(INOD,1)  !COEFFICIENT ABOUT VELOCITY OF WALL PARTICLES
        KK1=NWALLID(INOD,2)
        !IF(KK1.EQ.-10.OR.KK1.EQ.-11)THEN
          !UX(INOD,1:2)=0.D0;UY(INOD,1:2)=0.D0;UZ(INOD,1:2)=0.D0
          !GOTO 11
        !ENDIF
        !NDELP = ROW*N.U*-U(N+1)
        UM(INOD)=UX(INOD,2)*SNX(INOD)+
     +            UY(INOD,2)*SNY(INOD)+
     +            UZ(INOD,2)*SNZ(INOD) !+GRA*DT*SNZ(INOD)
C         UM(INOD)=UX(INOD,2)*SNX(INOD)+
C      +            UY(INOD,2)*SNY(INOD)+
C      +            UZ(INOD,2)*SNZ(INOD)+GRA*DT*SNZ(INOD)
        IF (NODEID(INOD).EQ.8) THEN
          IF(I_WM.EQ.1.OR.I_WM.EQ.2.OR.I_WM.EQ.3.OR.I_WM.EQ.5
     +      .OR.I_WM.EQ.6.OR.I_WM.EQ.15)THEN
            
            UM(INOD)=-TANK_A(INOD,1)*DT*SNX(INOD)+
     +                TANK_A(INOD,2)*DT*SNY(INOD)+
     +                TANK_A(INOD,3)*DT*SNZ(INOD)    ! NOTE:  NDELP =  ROW*N.U*-U(N+1)= ROW*N*ACC*DT
          END IF          
        END IF

11    CONTINUE
      !!$OMP END DO NOWAIT
      !!$OMP END PARALLEL

      END SUBROUTINE FIND_ACCN
!!----------------------------END FIND_ACCN---------------------------!!      



!!-----------------------------ALLOCATES------------------------------!!
      SUBROUTINE ALLOCATESTORAGE(NODN)
      USE MLPGSTORAGE
      IMPLICIT NONE
      INTEGER(KIND=4),INTENT(IN)::NODN
      INTEGER(KIND=4)::IERR
      
      ALLOCATE(IVV(0:NODN),STAT = IERR)
      IVV(0) = 500
      ALLOCATE(LINKTAB(NODN*IVV(0)),STAT = IERR)
      ALLOCATE(SKK(NODN*IVV(0)),STAT = IERR)
      ALLOCATE(PTCSR(NODN),FBCSR(NODN))
      ALLOCATE(IVVCSR(NODN+1))
      ALLOCATE(JVVCSR(IVV(0)*NODN),SKKCSR(IVV(0)*NODN))
      RETURN
      END SUBROUTINE ALLOCATESTORAGE

      

      SUBROUTINE ALLOCATEWORK(NODN,NTHR)
      USE MLPGSTORAGE
      IMPLICIT NONE

        INTEGER(KIND=4),INTENT(IN)::NODN,NTHR

        ALLOCATE(NWORK(NODN,NTHR),WORK(NODN,NTHR))
        WRITE(8,'(" [INF] : WORK ALLOCATED")')
        WRITE(8,'(" [---] : ",2A10)')"NODES","THREADS"
        WRITE(8,'(" [---] : ",2I10)')NODN,NTHR
        WRITE(8,*)
      END SUBROUTINE ALLOCATEWORK


      
      SUBROUTINE ALLOCATEKINE(NODN)
      USE MLPGKINE
      IMPLICIT NONE
      INTEGER(KIND=4),INTENT(IN)::NODN
      INTEGER(KIND=4)::IERR
      
      ALLOCATE (UX(NODN,1:4),STAT = IERR)
      ALLOCATE (UY(NODN,1:4),STAT = IERR)
      ALLOCATE (UZ(NODN,1:4),STAT = IERR)
      ALLOCATE (PPX(NODN,1:3),STAT = IERR)
      ALLOCATE (PPY(NODN,1:3),STAT = IERR)
      ALLOCATE (PPZ(NODN,1:3),STAT = IERR)
      ALLOCATE	(ROU(NODN),STAT = IERR)
      ALLOCATE (R(NODN),STAT = IERR)
      ALLOCATE (CC(NODN),STAT = IERR)
      ALLOCATE (R0(NODN),STAT = IERR)
      !ALLOCATE (PLAMDA(NODN),STAT = IERR)      
      
      END SUBROUTINE ALLOCATEKINE
      
      

C       SUBROUTINE ALLOCATEKINE1(NODN)
C       USE MLPGKINE
C       IMPLICIT NONE

C       INTEGER(KIND=4),INTENT(IN)::NODN
C       INTEGER(KIND=4)::IERR
      
C !      ALLOCATE (UX(NODN,1:4),STAT = IERR)
C !      ALLOCATE (UY(NODN,1:4),STAT = IERR)
C !      ALLOCATE (UZ(NODN,1:4),STAT = IERR)
C       ALLOCATE (PPX(NODN,1:3),STAT = IERR)
C       ALLOCATE (PPY(NODN,1:3),STAT = IERR)
C       ALLOCATE (PPZ(NODN,1:2),STAT = IERR)
C       ALLOCATE  (ROU(NODN),STAT = IERR)
C       ALLOCATE (R(NODN),STAT = IERR)
C       ALLOCATE (CC(NODN),STAT = IERR)
C       ALLOCATE (R0(NODN),STAT = IERR)
C !      ALLOCATE (PLAMDA(NODN),STAT = IERR)
C       RETURN
C       END SUBROUTINE ALLOCATEKINE1
!!----------------------------ALLOCATEKINE----------------------------!!      
      
            

!!----------------------------CHECK_NUMNE-----------------------------!!
      SUBROUTINE CHECK_NUMNE
      USE COMMONMOD
      USE MLPGKINE
      USE NEIGHNODES
      IMPLICIT NONE
      !INCLUDE 'COMMON.F'

      INTEGER(KIND=4)::I,J,II,IMIN,IMINN,IK,IS,MF
      REAL(KIND=8)::DX,DY,DZ,DR,RIAV,RMIN

      WRITE(8,*)'ENTERING CHECK_NUMNE',NODEID(-2),NLINK(10)%I(0)

      OPEN(NEWUNIT = MF, FILE='Export/CheckNumNe.dat')

      IMIN=999
      IMINN=0
      DO I=1,NODEID(-1)

        IK=NLINK(I)%I(0)        

        RIAV=R(I)+R0(I)        

        IS=0
        DO J=1,IK
          II=NLINK(I)%I(J)
          DX=COORX(I,1)-COORX(II,1)
          DY=COORY(I,1)-COORY(II,1)
          DZ=COORZ(I,1)-COORZ(II,1)
          DR=DSQRT(DX**2 + DY**2 + DZ**2)          

          IF(DR.LE.RIAV) IS=IS+1          

        ENDDO 

        IF((IMIN.GT.IS).AND.(NODEID(I).NE.8))THEN 
          IMIN=IS       
          IMINN=I
          RMIN=RIAV
        ENDIF

        WRITE(MF,'(2F15.6,2I10,F15.6)')COORX(I,1),COORZ(I,1),NODEID(I),
     +    IS,RIAV
      ENDDO

        WRITE(8,'(" [NENUM1] ",2I10,F15.6)')NODEID(IMINN),IMIN,RMIN
        WRITE(8,'(" [NENUM2] ",3F15.6)')COORX(IMINN,1),COORY(IMINN,1),
     +    COORZ(IMINN,1)

      CLOSE(MF)

      END SUBROUTINE CHECK_NUMNE
!!--------------------------END CHECK_NUMNE---------------------------!!      



!!---------------------------GENERATEGHOST----------------------------!!
      SUBROUTINE GENERATEGHOST(DDR)
      USE COMMONMOD
      USE MLPGKINE
      USE NEIGHNODES
      IMPLICIT NONE
      !INCLUDE 'COMMON.F'      

      REAL(KIND=8),INTENT(IN)::DDR(NODEID(0))

      INTEGER(KIND=4)::II,IJ
      REAL(KIND=8)::RIAV

      DO II=NODEID(-2)+1,NODEID(-7)
        IJ=NODEID(-1)+II-NODEID(-2)
        IF(IJ.GT.NODEID(0))THEN
          WRITE(8,'(" [ERR] CHECK GHOST GENERATE NGHT > NODEID(0)")')
          WRITE(8,'(" [ERR] ",2I10)')IJ,NODEID(0)
        ENDIF
        RIAV=0.7D0*DDR(II)
        COORX(IJ,:)=COORX(II,1)+SNX(II)*RIAV
        COORY(IJ,:)=COORY(II,1)+SNY(II)*RIAV
        COORZ(IJ,:)=COORZ(II,1)+SNZ(II)*RIAV        
        NODEID(IJ)=-6
        nwallid(IJ,2) = -10     
      ENDDO
C       IK=NODEID(-1)-NODEID(-2)
C       DO II=NODEID(-2)+1,NODEID(-1)
C         IJ=NODEID(-1)+IK+II-NODEID(-2)
C         IF(IJ.GT.NODEID(0))THEN
C           WRITE(8,'(" [ERR] CHECK GHOST GENERATE NGHT > NODEID(0)")')
C           WRITE(8,'(" [ERR] ",2I10)')IJ,NODEID(0)
C         ENDIF
C         RIAV=1D0*DDR(II)
C         COORX(IJ,:)=COORX(II,1)+SNX(II)*RIAV
C         COORY(IJ,:)=COORY(II,1)+SNY(II)*RIAV
C         COORZ(IJ,:)=COORZ(II,1)+SNZ(II)*RIAV        
C         NODEID(IJ)=-6
C         nwallid(IJ,2) = -10     
C       ENDDO

      END SUBROUTINE GENERATEGHOST      
!!---------------------------GENERATEGHOST----------------------------!!      



!!------------------------------CYLIND3-------------------------------!!
      SUBROUTINE CYLIND3(MESHFILE, FSNOD1, FSNOD2, DOMX, DOMY, DOMZ, 
     +  CYLX, CYLY, CYLR)      
      USE COMMONMOD
      IMPLICIT NONE
      !INCLUDE 'COMMON.F'

      REAL(KIND=8),INTENT(IN)::DOMX(2),DOMY(2),DOMZ(2),CYLX,CYLY,CYLR
      INTEGER(KIND=4),INTENT(OUT)::FSNOD1,FSNOD2      
      CHARACTER(LEN=256),INTENT(IN)::MESHFILE

      INTEGER(KIND=4)::MSFL=101
      INTEGER(KIND=4)::FSNODSTART,FSNODEND,NGHST      
      INTEGER(KIND=4)::NNI,IK,II,IJ,IJJ,IX,IY,IZ
      REAL(KIND=8)::DRMIN                  
      REAL(KIND=8)::TMPR1,TMPR2,TMPR3,TMPR4,TMPR7,RIAV
      CHARACTER(LEN=256)::TEXT1

      WRITE(8,*)'[MSG] ENTERING CYLIND3'
      
      DRMIN=0.0001D0      

      SNX = 0;SNY=0;SNZ = 0; !X_DIRECTIONCOSINES
      SMX = 0;SMY=0;SMZ = 0; !Y_DIRECTIONCOSINES
      SSX = 0;SSY=0;SSZ = 0; !Z_DIRECTIONCOSINES

      OPEN(MSFL,FILE=TRIM(MESHFILE))

      READ(MSFL,*)TEXT1
      READ(MSFL,*)(NODEID(IX),IX=-1,-7,-1)
      READ(MSFL,*)TEXT1
      READ(MSFL,*)FSNODSTART,FSNODEND
      READ(MSFL,*)TEXT1
      READ(MSFL,*)FSNOD1
      READ(MSFL,*)FSNOD2
      READ(MSFL,*)TEXT1
      !WRITE(8,*)TEXT1(1:LEN_TRIM(TEXT1))
      
      WRITE(8,*)NODEID(-7:-1)

      !! FLUID + FS
      !IK=0
      DO NNI=1,NODEID(-2)
        READ(MSFL,*)COORX(NNI,1),COORY(NNI,1),COORZ(NNI,1)
        NODEID(NNI)=0

C         !! RNODE
C         IF(IK.EQ.0)THEN
C           TMPR1=COORX(NNI,1)-(DOMX(1)+DOMX(2))/2D0
C           TMPR2=COORY(NNI,1)-(DOMY(1)+DOMY(2))/2D0
C           TMPR3=COORZ(NNI,1)-(DOMZ(1)+DOMZ(2))/2D0
C           TMPR7=DSQRT(TMPR1**2 + TMPR2**2 + TMPR3**2)
C           IF(TMPR7.LT.0.05)THEN 
C             NODEID(NNI)=5
C             IK=1
C           ENDIF
C         ENDIF
      ENDDO
      NODEID(FSNODSTART:FSNODEND)=4
C       IF(IK.EQ.0)THEN
C         WRITE(8,*)'[ERR] RNODE NOT FOUND'
C       ENDIF

      !! WAVEMAKER
      DO NNI=NODEID(-2)+1,NODEID(-3)
        READ(MSFL,*)COORX(NNI,1),COORY(NNI,1),COORZ(NNI,1)
        NODEID(NNI)=8
        NWALLID(NNI,1)=2
        NWALLID(NNI,2) = 1
        SNX(NNI) = -1
        SMY(NNI) = -1
        SSZ(NNI) =  1

        IF(ABS(COORZ(NNI,1)-DOMZ(2)).LT.DRMIN)THEN 
          NWALLID(NNI,2)=-11
        ENDIF

        IF(ABS(COORZ(NNI,1)-DOMZ(1)).LT.DRMIN)THEN
          NWALLID(NNI,3)=9  !SPECIAL NODE FOR PRESSURE
          NWALLID(NNI,1)=6  
          SSZ(NNI) = 0  
        ENDIF

        IF((ABS(COORY(NNI,1)-DOMY(1)).LT.DRMIN).OR.
     +      (ABS(COORY(NNI,1)-DOMY(2)).LT.DRMIN)) THEN
          NWALLID(NNI,1)=7 !Intersection of surfaces
          NWALLID(NNI,3)=9    !SPECIAL NODE FOR PRESSURE EQUATION 
          SMY(NNI) = 0
        ENDIF
      ENDDO

      !! BOTTOM
      DO NNI=NODEID(-3)+1,NODEID(-4)
        READ(MSFL,*)COORX(NNI,1),COORY(NNI,1),COORZ(NNI,1)
        NODEID(NNI)=2
        NWALLID(NNI,2) = 1
        NWALLID(NNI,1)=1
        SNZ(NNI) = -1
        SMY(NNI) =  1
        SSX(NNI) =  1
      ENDDO

      !! OPPOSITE TO WAVEMAKER
      DO NNI=NODEID(-4)+1,NODEID(-5)
        READ(MSFL,*)COORX(NNI,1),COORY(NNI,1),COORZ(NNI,1)
        NODEID(NNI)=3
        NWALLID(NNI,1)=2
        NWALLID(NNI,2) = 1
        SNX(NNI) =  1
        SMY(NNI) =  1
        SSZ(NNI) =  1

        IF(ABS(COORZ(NNI,1)-DOMZ(2)).LT.DRMIN)THEN 
          NWALLID(NNI,2)=-11
        ENDIF

        IF(ABS(COORZ(NNI,1)-DOMZ(1)).LT.DRMIN)THEN
          NWALLID(NNI,3)=9  !SPECIAL NODE FOR PRESSURE
          NWALLID(NNI,1)=6  
          SSZ(NNI) = 0  
        ENDIF

        IF((ABS(COORY(NNI,1)-DOMY(1)).LT.DRMIN).OR.
     +      (ABS(COORY(NNI,1)-DOMY(2)).LT.DRMIN)) THEN
          NWALLID(NNI,1)=7 !Intersection of surfaces
          NWALLID(NNI,3)=9    !SPECIAL NODE FOR PRESSURE EQUATION 
          SMY(NNI) = 0

          IF(ABS(COORZ(NNI,1)-DOMZ(1)).LT.DRMIN)THEN
            NWALLID(NNI,1)=0    !ONLY 0 VELOCITY 
          ENDIF
        ENDIF          
      ENDDO

      !! SIDEWALL NEAR
      DO NNI=NODEID(-5)+1,NODEID(-6)
        READ(MSFL,*)COORX(NNI,1),COORY(NNI,1),COORZ(NNI,1)
        NODEID(NNI)=1
        NWALLID(NNI,2) = 1
        NWALLID(NNI,1)=3
        SNY(NNI) = -1
        SMX(NNI) =  1
        SSZ(NNI) =  1

        IF(ABS(COORZ(NNI,1)-DOMZ(2)).LT.DRMIN)THEN 
          NWALLID(NNI,2)=-11
        ENDIF

        IF(ABS(COORZ(NNI,1)-DOMZ(1)).LT.DRMIN)THEN
          NWALLID(NNI,3)=9  !SPECIAL NODE FOR PRESSURE
          NWALLID(NNI,1)=5  
          SSZ(NNI) = 0  
        ENDIF
      ENDDO

      !! SIDEWALL FAR
      DO NNI=NODEID(-6)+1,NODEID(-7)
        READ(MSFL,*)COORX(NNI,1),COORY(NNI,1),COORZ(NNI,1)
        NODEID(NNI)=7
        NWALLID(NNI,2) = 1
        NWALLID(NNI,1)=3
        SNY(NNI) =  1
        SMX(NNI) = -1
        SSZ(NNI) =  1

        IF(ABS(COORZ(NNI,1)-DOMZ(2)).LT.DRMIN)THEN 
          NWALLID(NNI,2)=-11
        ENDIF

        IF(ABS(COORZ(NNI,1)-DOMZ(1)).LT.DRMIN)THEN
          NWALLID(NNI,3)=9  !SPECIAL NODE FOR PRESSURE
          NWALLID(NNI,1)=5  
          SSZ(NNI) = 0  
        ENDIF
      ENDDO

      !! CYLINDER 
      DO NNI=NODEID(-7)+1,NODEID(-1)
        READ(MSFL,*)COORX(NNI,1),COORY(NNI,1),COORZ(NNI,1)        
        NODEID(NNI)=9
        IF((COORZ(NNI,1)-0D0).LT.DRMIN)THEN
          NWALLID(NNI,3)=9  !PRESSURE BC ID
        ENDIF
        !NWALLID(NNI,3)=9  !PRESSURE BC ID
        NWALLID(NNI,1)=4  !VELOCITY BC ID 
        !The line above is wrong mostly but shouldnt make a diff

        TMPR1=COORX(NNI,1)-CYLX
        TMPR2=COORY(NNI,1)-CYLY
        TMPR3=DSQRT(TMPR1**2 + TMPR2**2)
        TMPR1=TMPR1/TMPR3
        TMPR2=TMPR2/TMPR3
        SNX(NNI)=-TMPR1;SNY(NNI)=-TMPR2;SNZ(NNI)=0D0
        SMX(NNI)=TMPR2;SMY(NNI)=-TMPR1;SMZ(NNI)=0D0
        SSX(NNI)=0D0;SSY(NNI)=0D0;SSZ(NNI)=1D0            
C         WRITE(8,'(4F15.6,I15)')SNX(NNI),SNY(NNI),SNZ(NNI),
C      +    (SNX(NNI)**2+SNY(NNI)**2),NODEID(NNI)
      ENDDO
      CLOSE(MSFL)

      !! GHOST NODES
      NGHST=0
      NGHST=(NODEID(-1)-NODEID(-2))
      !NGHST=(NODEID(-1)-NODEID(-7))
      NODEID(0)=NODEID(-1)+NGHST   
      DO II=NODEID(-2)+1,NODEID(-7)
        IJ=NODEID(-1)+II-NODEID(-2)
        IF(IJ.GT.NODEID(0))THEN
          WRITE(8,'(" [ERR] CHECK GHOST GENERATE NGHT > NODEID(0)")')
          WRITE(8,'(" [ERR] ",2I10)')IJ,NODEID(0)
        ENDIF
        RIAV=0.3D0*DDL !! Keeping this lower than GENERATEGHOST(DDR)
        COORX(IJ,:)=COORX(II,1)+SNX(II)*RIAV
        COORY(IJ,:)=COORY(II,1)+SNY(II)*RIAV
        COORZ(IJ,:)=COORZ(II,1)+SNZ(II)*RIAV   
        NODEID(IJ)=-6
        nwallid(IJ,2) = -10     
      ENDDO   
      !IJJ=NODEID(-1)
      IJJ=IJ      
      DO IX=NODEID(-7)+1,NODEID(-1)
        TMPR7=1D0/DRMIN
        IK=0
        DO IY=NODEID(-7)+1,NODEID(-1)
          IF(IX.EQ.IY)CYCLE
          TMPR1=COORX(IY,1)-COORX(IX,1)
          TMPR2=COORY(IY,1)-COORY(IX,1)
          TMPR3=COORZ(IY,1)-COORZ(IX,1)
          TMPR4=DSQRT(TMPR1**2 + TMPR2**2 + TMPR3**2)
          IF(TMPR4.LT.TMPR7)THEN
            TMPR7=TMPR4
            IK=IY
          ENDIF
        ENDDO

        IF(IK.EQ.0)THEN
          WRITE(8,*)"[ERR] CHECK CYLIND GHOST",IX
          WRITE(8,*)"[---] ",COORX(IX,1),COORY(IX,1),COORZ(IX,1)
          STOP
        ENDIF

        TMPR7=CYLR-TMPR7        
        NNI=IJJ+IX-NODEID(-7)
        IF(NNI.GT.NODEID(0))THEN
          WRITE(8,'(" [ERR] CHECK GHOST GENERATE NGHT > NODEID(0)")')
          WRITE(8,'(" [ERR] ",2I10)')NNI,NODEID(0)
        ENDIF
        COORX(NNI,1)=CYLX-TMPR7*SNX(IX)
        COORY(NNI,1)=CYLY-TMPR7*SNY(IX)
        COORZ(NNI,1)=COORZ(IX,1)
        NODEID(NNI)=-9
        NWALLID(NNI,2)=-10
      ENDDO            

      COORX(:,2)=COORX(:,1)
      COORY(:,2)=COORY(:,1)
      COORZ(:,2)=COORZ(:,1)


      IZ=INT((DOMZ(2)-DOMZ(1))/DDL,4)
      IY=INT((DOMY(2)-DOMY(1))/DDL,4)
      IX=INT((DOMX(2)-DOMX(1))/DDL,4)

      IXMAX=IX+20
      IYMAX=IY+20
      IZMAX=IZ+20

      END SUBROUTINE CYLIND3
!!----------------------------END CYLIND3-----------------------------!!



!!------------------------------CYLIND4-------------------------------!!
      SUBROUTINE CYLIND4(MESHFILE, FSNOD1, FSNOD2, DOMX, DOMY, DOMZ, 
     +  CYLX, CYLY, CYLR)      
      USE COMMONMOD
      IMPLICIT NONE
      !INCLUDE 'COMMON.F'

      REAL(KIND=8),INTENT(IN)::DOMX(2),DOMY(2),DOMZ(2),CYLX,CYLY,CYLR
      INTEGER(KIND=4),INTENT(OUT)::FSNOD1,FSNOD2      
      CHARACTER(LEN=256),INTENT(IN)::MESHFILE

      INTEGER(KIND=4)::MSFL=101
      INTEGER(KIND=4)::FSNODSTART,FSNODEND,NGHST      
      INTEGER(KIND=4)::NNI,IK,II,IJ,IJJ,IX,IY,IZ
      INTEGER(KIND=4)::NODEDGETY(2)
      REAL(KIND=8)::DRMIN                  
      REAL(KIND=8)::TMPR1,TMPR2,TMPR3,TMPR4,TMPR7,RIAV
      CHARACTER(LEN=256)::TEXT1

      !! nodEgdeTy(a,b)
      !  Edge types defined assuming a rectilinear domain
      !  1,0 = bottom or near edge
      !  2,0 = top or far edge
      !  0,1 = side edge

      WRITE(8,*)'[MSG] ENTERING CYLIND3'
      
      DRMIN=0.0001D0      

      SNX = 0;SNY=0;SNZ = 0; !X_DIRECTIONCOSINES
      SMX = 0;SMY=0;SMZ = 0; !Y_DIRECTIONCOSINES
      SSX = 0;SSY=0;SSZ = 0; !Z_DIRECTIONCOSINES

      OPEN(MSFL,FILE=TRIM(MESHFILE))

      READ(MSFL,*)TEXT1
      READ(MSFL,*)(NODEID(IX),IX=-1,-7,-1)
      READ(MSFL,*)TEXT1
      READ(MSFL,*)FSNODSTART,FSNODEND
      READ(MSFL,*)TEXT1
      READ(MSFL,*)FSNOD1
      READ(MSFL,*)FSNOD2
      READ(MSFL,*)TEXT1
      !WRITE(8,*)TEXT1(1:LEN_TRIM(TEXT1))
      
      WRITE(8,*)NODEID(-7:-1)

      !! FLUID + FS
      !IK=0
      DO NNI=1,NODEID(-2)
        READ(MSFL,*)COORX(NNI,1),COORY(NNI,1),COORZ(NNI,1),
     +    SNX(NNI), SNY(NNI), SNZ(NNI), 
     +    SMX(NNI), SMY(NNI), SMZ(NNI),
     +    SSX(NNI), SSY(NNI), SSZ(NNI), NODEDGETY(1:2)
        NODEID(NNI)=0
      ENDDO
      NODEID(FSNODSTART:FSNODEND)=4


      !! WAVEMAKER
      DO NNI=NODEID(-2)+1,NODEID(-3)
        READ(MSFL,*)COORX(NNI,1),COORY(NNI,1),COORZ(NNI,1),
     +    SNX(NNI), SNY(NNI), SNZ(NNI), 
     +    SMX(NNI), SMY(NNI), SMZ(NNI),
     +    SSX(NNI), SSY(NNI), SSZ(NNI), NODEDGETY(1:2)

        NODEID(NNI)=8
        NWALLID(NNI,1)=2
        NWALLID(NNI,2) = 1
     
        
        IF(NODEDGETY(1).EQ.2)THEN !TOP
          NWALLID(NNI,2)=-11
        ENDIF
        
        IF(NODEDGETY(1).EQ.1)THEN !BOTTOM
          NWALLID(NNI,1)=6  
          NWALLID(NNI,3)=9  !SPECIAL NODE FOR PRESSURE          
          SSZ(NNI) = 0  
        ENDIF

        IF(NODEDGETY(2).EQ.1)THEN !SIDE
          NWALLID(NNI,1)=7 !Intersection of surfaces
          NWALLID(NNI,3)=9    !SPECIAL NODE FOR PRESSURE EQUATION 
          SMY(NNI) = 0
        ENDIF
      ENDDO

      !! BOTTOM
      DO NNI=NODEID(-3)+1,NODEID(-4)
        READ(MSFL,*)COORX(NNI,1),COORY(NNI,1),COORZ(NNI,1),
     +    SNX(NNI), SNY(NNI), SNZ(NNI), 
     +    SMX(NNI), SMY(NNI), SMZ(NNI),
     +    SSX(NNI), SSY(NNI), SSZ(NNI), NODEDGETY(1:2)

        NODEID(NNI)=2
        NWALLID(NNI,2) = 1
        NWALLID(NNI,1)=1     
      ENDDO

      !! OPPOSITE TO WAVEMAKER
      DO NNI=NODEID(-4)+1,NODEID(-5)
        READ(MSFL,*)COORX(NNI,1),COORY(NNI,1),COORZ(NNI,1),
     +    SNX(NNI), SNY(NNI), SNZ(NNI), 
     +    SMX(NNI), SMY(NNI), SMZ(NNI),
     +    SSX(NNI), SSY(NNI), SSZ(NNI), NODEDGETY(1:2)
        NODEID(NNI)=3
        NWALLID(NNI,1)=2
        NWALLID(NNI,2) = 1
        

        IF(NODEDGETY(1).EQ.2)THEN !TOP
          NWALLID(NNI,2)=-11
        ENDIF

        IF(NODEDGETY(1).EQ.1)THEN !BOTTOM
          NWALLID(NNI,1)=6  
          NWALLID(NNI,3)=9  !SPECIAL NODE FOR PRESSURE          
          SSZ(NNI) = 0  
        ENDIF

        IF(NODEDGETY(2).EQ.1)THEN !SIDE
          NWALLID(NNI,1)=7 !Intersection of surfaces
          NWALLID(NNI,3)=9    !SPECIAL NODE FOR PRESSURE EQUATION 
          SMY(NNI) = 0

          IF(NODEDGETY(1).EQ.1)THEN !BOTTOM
            NWALLID(NNI,1)=0    !ONLY 0 VELOCITY 
          ENDIF
        ENDIF          
      ENDDO

      !! SIDEWALL NEAR
      DO NNI=NODEID(-5)+1,NODEID(-6)
        READ(MSFL,*)COORX(NNI,1),COORY(NNI,1),COORZ(NNI,1),
     +    SNX(NNI), SNY(NNI), SNZ(NNI), 
     +    SMX(NNI), SMY(NNI), SMZ(NNI),
     +    SSX(NNI), SSY(NNI), SSZ(NNI), NODEDGETY(1:2)

        NODEID(NNI)=1
        NWALLID(NNI,2) = 1
        NWALLID(NNI,1)=3
     

        IF(NODEDGETY(1).EQ.2)THEN !TOP
          NWALLID(NNI,2)=-11
        ENDIF

        IF(NODEDGETY(1).EQ.1)THEN !BOTTOM
          NWALLID(NNI,3)=9  !SPECIAL NODE FOR PRESSURE
          NWALLID(NNI,1)=5  
          SSZ(NNI) = 0  
        ENDIF
      ENDDO

      !! SIDEWALL FAR
      DO NNI=NODEID(-6)+1,NODEID(-7)
        READ(MSFL,*)COORX(NNI,1),COORY(NNI,1),COORZ(NNI,1),
     +    SNX(NNI), SNY(NNI), SNZ(NNI), 
     +    SMX(NNI), SMY(NNI), SMZ(NNI),
     +    SSX(NNI), SSY(NNI), SSZ(NNI), NODEDGETY(1:2)

        NODEID(NNI)=7
        NWALLID(NNI,2) = 1
        NWALLID(NNI,1)=3
     

        IF(NODEDGETY(1).EQ.2)THEN !TOP
          NWALLID(NNI,2)=-11
        ENDIF

        IF(NODEDGETY(1).EQ.1)THEN !BOTTOM
          NWALLID(NNI,3)=9  !SPECIAL NODE FOR PRESSURE
          NWALLID(NNI,1)=5  
          SSZ(NNI) = 0  
        ENDIF
      ENDDO

      !! CYLINDER 
      DO NNI=NODEID(-7)+1,NODEID(-1)
        READ(MSFL,*)COORX(NNI,1),COORY(NNI,1),COORZ(NNI,1),
     +    SNX(NNI), SNY(NNI), SNZ(NNI), 
     +    SMX(NNI), SMY(NNI), SMZ(NNI),
     +    SSX(NNI), SSY(NNI), SSZ(NNI), NODEDGETY(1:2)

        NODEID(NNI)=9
        IF(NODEDGETY(1).EQ.1)THEN !BOTTOM
          NWALLID(NNI,3)=9  !PRESSURE BC ID
        ENDIF
        !NWALLID(NNI,3)=9  !PRESSURE BC ID
        !NWALLID(NNI,1)=4  !VELOCITY BC ID 
        !The line above is wrong mostly but shouldnt make a diff
        !commented that nwallid line on 2021-04-07
        
C         WRITE(8,'(4F15.6,I15)')SNX(NNI),SNY(NNI),SNZ(NNI),
C      +    (SNX(NNI)**2+SNY(NNI)**2),NODEID(NNI)
      ENDDO
      CLOSE(MSFL)

      !! GHOST NODES
      NGHST=0
      NGHST=(NODEID(-1)-NODEID(-2))
      !NGHST=(NODEID(-1)-NODEID(-7))
      NODEID(0)=NODEID(-1)+NGHST   
      DO II=NODEID(-2)+1,NODEID(-7)
        IJ=NODEID(-1)+II-NODEID(-2)
        IF(IJ.GT.NODEID(0))THEN
          WRITE(8,'(" [ERR] CHECK GHOST GENERATE NGHT > NODEID(0)")')
          WRITE(8,'(" [ERR] ",2I10)')IJ,NODEID(0)
        ENDIF
        RIAV=0.3D0*DDL !! Keeping this lower than GENERATEGHOST(DDR)
        COORX(IJ,:)=COORX(II,1)+SNX(II)*RIAV
        COORY(IJ,:)=COORY(II,1)+SNY(II)*RIAV
        COORZ(IJ,:)=COORZ(II,1)+SNZ(II)*RIAV   
        NODEID(IJ)=-6
        nwallid(IJ,2) = -10     
      ENDDO   
      !IJJ=NODEID(-1)
      IJJ=IJ      
      DO IX=NODEID(-7)+1,NODEID(-1)
        TMPR7=1D0/DRMIN
        IK=0
        DO IY=NODEID(-7)+1,NODEID(-1)
          IF(IX.EQ.IY)CYCLE
          TMPR1=COORX(IY,1)-COORX(IX,1)
          TMPR2=COORY(IY,1)-COORY(IX,1)
          TMPR3=COORZ(IY,1)-COORZ(IX,1)
          TMPR4=DSQRT(TMPR1**2 + TMPR2**2 + TMPR3**2)
          IF(TMPR4.LT.TMPR7)THEN
            TMPR7=TMPR4
            IK=IY
          ENDIF
        ENDDO

        IF(IK.EQ.0)THEN
          WRITE(8,*)"[ERR] CHECK CYLIND GHOST",IX
          WRITE(8,*)"[---] ",COORX(IX,1),COORY(IX,1),COORZ(IX,1)
          STOP
        ENDIF

        TMPR7=0.80D0*TMPR7        
        NNI=IJJ+IX-NODEID(-7)
        IF(NNI.GT.NODEID(0))THEN
          WRITE(8,'(" [ERR] CHECK GHOST GENERATE NGHT > NODEID(0)")')
          WRITE(8,'(" [ERR] ",2I10)')NNI,NODEID(0)
        ENDIF
        COORX(NNI,1)=COORX(IX,1)+SNX(IX)*TMPR7
        COORY(NNI,1)=COORY(IX,1)+SNY(IX)*TMPR7
        COORZ(NNI,1)=COORZ(IX,1)
        NODEID(NNI)=-9
        NWALLID(NNI,2)=-10
        WRITE(8,'(" [INF] CYL AND GHOST COOR")')
        WRITE(8,'(" [---] ",3F15.6)')COORX(IX,1), COORY(IX,1), 
     +    COORZ(IX,1)
        WRITE(8,'(" [---] ",3F15.6)')COORX(NNI,1), COORY(NNI,1), 
     +    COORZ(NNI,1)        
      ENDDO            
      WRITE(8,*)

      COORX(:,2)=COORX(:,1)
      COORY(:,2)=COORY(:,1)
      COORZ(:,2)=COORZ(:,1)


      IZ=INT((DOMZ(2)-DOMZ(1))/DDL,4)
      IY=INT((DOMY(2)-DOMY(1))/DDL,4)
      IX=INT((DOMX(2)-DOMX(1))/DDL,4)

      IXMAX=IX+20
      IYMAX=IY+20
      IZMAX=IZ+20

      END SUBROUTINE CYLIND4
!!----------------------------END CYLIND4-----------------------------!!



!!------------------------------CYLIND5-------------------------------!!
      SUBROUTINE CYLIND5(MESHFILE, FSNOD1, FSNOD2, DOMX, DOMY, DOMZ, 
     +  CYLX, CYLY, CYLR)      
      USE COMMONMOD
      IMPLICIT NONE
      !INCLUDE 'COMMON.F'

      REAL(KIND=8),INTENT(IN)::DOMX(2),DOMY(2),DOMZ(2),CYLX,CYLY,CYLR
      INTEGER(KIND=4),INTENT(OUT)::FSNOD1,FSNOD2      
      CHARACTER(LEN=256),INTENT(IN)::MESHFILE

      INTEGER(KIND=4)::MSFL=101
      INTEGER(KIND=4)::FSNODSTART,FSNODEND,NGHST      
      INTEGER(KIND=4)::NNI,IK,II,IJ,IJJ,IX,IY,IZ
      INTEGER(KIND=4)::NODEDGETY(2)
      REAL(KIND=8)::DRMIN                  
      REAL(KIND=8)::TMPR1,TMPR2,TMPR3,TMPR4,TMPR7,RIAV
      CHARACTER(LEN=256)::TEXT1

      !! nodEgdeTy(a,b)
      !  Edge types defined assuming a rectilinear domain
      !  1,0 = bottom or near edge
      !  2,0 = top or far edge
      !  0,1 = side edge

      WRITE(8,*)'[MSG] ENTERING CYLIND3'
      
      DRMIN=0.0001D0      

      SNX = 0;SNY=0;SNZ = 0; !X_DIRECTIONCOSINES
      SMX = 0;SMY=0;SMZ = 0; !Y_DIRECTIONCOSINES
      SSX = 0;SSY=0;SSZ = 0; !Z_DIRECTIONCOSINES

      OPEN(MSFL,FILE=TRIM(MESHFILE))

      READ(MSFL,*)TEXT1
      READ(MSFL,*)(NODEID(IX),IX=-1,-7,-1)
      READ(MSFL,*)TEXT1
      READ(MSFL,*)FSNODSTART,FSNODEND
      READ(MSFL,*)TEXT1
      READ(MSFL,*)FSNOD1
      READ(MSFL,*)FSNOD2
      READ(MSFL,*)TEXT1
      !WRITE(8,*)TEXT1(1:LEN_TRIM(TEXT1))
      
      WRITE(8,*)NODEID(-7:-1)

      !! FLUID + FS
      !IK=0
      DO NNI=1,NODEID(-2)
        READ(MSFL,*)COORX(NNI,1),COORY(NNI,1),COORZ(NNI,1),
     +    SNX(NNI), SNY(NNI), SNZ(NNI), 
     +    SMX(NNI), SMY(NNI), SMZ(NNI),
     +    SSX(NNI), SSY(NNI), SSZ(NNI), NODEDGETY(1:2)
        NODEID(NNI)=0
      ENDDO
      NODEID(FSNODSTART:FSNODEND)=4


      !! WAVEMAKER
      DO NNI=NODEID(-2)+1,NODEID(-3)
        READ(MSFL,*)COORX(NNI,1),COORY(NNI,1),COORZ(NNI,1),
     +    SNX(NNI), SNY(NNI), SNZ(NNI), 
     +    SMX(NNI), SMY(NNI), SMZ(NNI),
     +    SSX(NNI), SSY(NNI), SSZ(NNI), NODEDGETY(1:2)

        NODEID(NNI)=8
        NWALLID(NNI,1)=2
        NWALLID(NNI,2) = 1
     
        
        IF(NODEDGETY(1).EQ.2)THEN !TOP
          NWALLID(NNI,2)=-11
        ENDIF
        
C         IF(NODEDGETY(1).EQ.1)THEN !BOTTOM
C           NWALLID(NNI,1)=6  
C           NWALLID(NNI,3)=9  !SPECIAL NODE FOR PRESSURE          
C           SSZ(NNI) = 0  
C         ENDIF

        IF(NODEDGETY(2).EQ.1)THEN !SIDE
          NWALLID(NNI,1)=7 !Intersection of surfaces
          NWALLID(NNI,3)=9    !SPECIAL NODE FOR PRESSURE EQUATION 
          !SMY(NNI) = 0
        ENDIF
      ENDDO

      !! BOTTOM
      DO NNI=NODEID(-3)+1,NODEID(-4)
        READ(MSFL,*)COORX(NNI,1),COORY(NNI,1),COORZ(NNI,1),
     +    SNX(NNI), SNY(NNI), SNZ(NNI), 
     +    SMX(NNI), SMY(NNI), SMZ(NNI),
     +    SSX(NNI), SSY(NNI), SSZ(NNI), NODEDGETY(1:2)

        NODEID(NNI)=2
        NWALLID(NNI,2) = 1
        NWALLID(NNI,1)=1     

        IF(NODEDGETY(2).EQ.1)THEN !SIDE
          NWALLID(NNI,1)=7 !Intersection of surfaces
          NWALLID(NNI,3)=9    !SPECIAL NODE FOR PRESSURE EQUATION           
        ENDIF

      ENDDO

      !! OPPOSITE TO WAVEMAKER
      DO NNI=NODEID(-4)+1,NODEID(-5)
        READ(MSFL,*)COORX(NNI,1),COORY(NNI,1),COORZ(NNI,1),
     +    SNX(NNI), SNY(NNI), SNZ(NNI), 
     +    SMX(NNI), SMY(NNI), SMZ(NNI),
     +    SSX(NNI), SSY(NNI), SSZ(NNI), NODEDGETY(1:2)
        NODEID(NNI)=3
        NWALLID(NNI,1)=2
        NWALLID(NNI,2) = 1
        

        IF(NODEDGETY(1).EQ.2)THEN !TOP
          NWALLID(NNI,2)=-11
        ENDIF

C         IF(NODEDGETY(1).EQ.1)THEN !BOTTOM
C           NWALLID(NNI,1)=6  
C           NWALLID(NNI,3)=9  !SPECIAL NODE FOR PRESSURE          
C           SSZ(NNI) = 0  
C         ENDIF

        IF(NODEDGETY(2).EQ.1)THEN !SIDE
          NWALLID(NNI,1)=7 !Intersection of surfaces
          NWALLID(NNI,3)=9    !SPECIAL NODE FOR PRESSURE EQUATION 
          !SMY(NNI) = 0

C           IF(NODEDGETY(1).EQ.1)THEN !BOTTOM
C             NWALLID(NNI,1)=0    !ONLY 0 VELOCITY 
C           ENDIF
        ENDIF          
      ENDDO

      !! SIDEWALL NEAR
      DO NNI=NODEID(-5)+1,NODEID(-6)
        READ(MSFL,*)COORX(NNI,1),COORY(NNI,1),COORZ(NNI,1),
     +    SNX(NNI), SNY(NNI), SNZ(NNI), 
     +    SMX(NNI), SMY(NNI), SMZ(NNI),
     +    SSX(NNI), SSY(NNI), SSZ(NNI), NODEDGETY(1:2)

        NODEID(NNI)=1
        NWALLID(NNI,2) = 1
        NWALLID(NNI,1)=3
     

        IF(NODEDGETY(1).EQ.2)THEN !TOP
          NWALLID(NNI,2)=-11
        ENDIF

C         IF(NODEDGETY(1).EQ.1)THEN !BOTTOM
C           NWALLID(NNI,3)=9  !SPECIAL NODE FOR PRESSURE
C           NWALLID(NNI,1)=5  
C           SSZ(NNI) = 0  
C         ENDIF
      ENDDO

      !! SIDEWALL FAR
      DO NNI=NODEID(-6)+1,NODEID(-7)
        READ(MSFL,*)COORX(NNI,1),COORY(NNI,1),COORZ(NNI,1),
     +    SNX(NNI), SNY(NNI), SNZ(NNI), 
     +    SMX(NNI), SMY(NNI), SMZ(NNI),
     +    SSX(NNI), SSY(NNI), SSZ(NNI), NODEDGETY(1:2)

        NODEID(NNI)=7
        NWALLID(NNI,2) = 1
        NWALLID(NNI,1)=3
     

        IF(NODEDGETY(1).EQ.2)THEN !TOP
          NWALLID(NNI,2)=-11
        ENDIF

C         IF(NODEDGETY(1).EQ.1)THEN !BOTTOM
C           NWALLID(NNI,3)=9  !SPECIAL NODE FOR PRESSURE
C           NWALLID(NNI,1)=5  
C           SSZ(NNI) = 0  
C         ENDIF
      ENDDO

      !! CYLINDER 
      DO NNI=NODEID(-7)+1,NODEID(-1)
        READ(MSFL,*)COORX(NNI,1),COORY(NNI,1),COORZ(NNI,1),
     +    SNX(NNI), SNY(NNI), SNZ(NNI), 
     +    SMX(NNI), SMY(NNI), SMZ(NNI),
     +    SSX(NNI), SSY(NNI), SSZ(NNI), NODEDGETY(1:2)

        NODEID(NNI)=9
        IF(NODEDGETY(1).EQ.1)THEN !BOTTOM
          NWALLID(NNI,3)=9  !PRESSURE BC ID
        ENDIF
        !NWALLID(NNI,3)=9  !PRESSURE BC ID
        !NWALLID(NNI,1)=4  !VELOCITY BC ID 
        !The line above is wrong mostly but shouldnt make a diff
        !commented that nwallid line on 2021-04-07
        
C         WRITE(8,'(4F15.6,I15)')SNX(NNI),SNY(NNI),SNZ(NNI),
C      +    (SNX(NNI)**2+SNY(NNI)**2),NODEID(NNI)
      ENDDO
      CLOSE(MSFL)

      !! GHOST NODES
      NGHST=0
      NGHST=(NODEID(-1)-NODEID(-2))
      !NGHST=(NODEID(-1)-NODEID(-7))
      NODEID(0)=NODEID(-1)+NGHST   
      DO II=NODEID(-2)+1,NODEID(-7)
        IJ=NODEID(-1)+II-NODEID(-2)
        IF(IJ.GT.NODEID(0))THEN
          WRITE(8,'(" [ERR] CHECK GHOST GENERATE NGHT > NODEID(0)")')
          WRITE(8,'(" [ERR] ",2I10)')IJ,NODEID(0)
        ENDIF
        RIAV=0.3D0*DDL !! Keeping this lower than GENERATEGHOST(DDR)
        COORX(IJ,:)=COORX(II,1)+SNX(II)*RIAV
        COORY(IJ,:)=COORY(II,1)+SNY(II)*RIAV
        COORZ(IJ,:)=COORZ(II,1)+SNZ(II)*RIAV   
        NODEID(IJ)=-6
        nwallid(IJ,2) = -10     
      ENDDO   
      !IJJ=NODEID(-1)
      IJJ=IJ      
      DO IX=NODEID(-7)+1,NODEID(-1)
        TMPR7=1D0/DRMIN
        IK=0
        DO IY=NODEID(-7)+1,NODEID(-1)
          IF(IX.EQ.IY)CYCLE
          TMPR1=COORX(IY,1)-COORX(IX,1)
          TMPR2=COORY(IY,1)-COORY(IX,1)
          TMPR3=COORZ(IY,1)-COORZ(IX,1)
          TMPR4=DSQRT(TMPR1**2 + TMPR2**2 + TMPR3**2)
          IF(TMPR4.LT.TMPR7)THEN
            TMPR7=TMPR4
            IK=IY
          ENDIF
        ENDDO

        IF(IK.EQ.0)THEN
          WRITE(8,*)"[ERR] CHECK CYLIND GHOST",IX
          WRITE(8,*)"[---] ",COORX(IX,1),COORY(IX,1),COORZ(IX,1)
          STOP
        ENDIF

        TMPR7=0.80D0*TMPR7        
        NNI=IJJ+IX-NODEID(-7)
        IF(NNI.GT.NODEID(0))THEN
          WRITE(8,'(" [ERR] CHECK GHOST GENERATE NGHT > NODEID(0)")')
          WRITE(8,'(" [ERR] ",2I10)')NNI,NODEID(0)
        ENDIF
        COORX(NNI,1)=COORX(IX,1)+SNX(IX)*TMPR7
        COORY(NNI,1)=COORY(IX,1)+SNY(IX)*TMPR7
        COORZ(NNI,1)=COORZ(IX,1)
        NODEID(NNI)=-9
        NWALLID(NNI,2)=-10
        WRITE(8,'(" [INF] CYL AND GHOST COOR")')
        WRITE(8,'(" [---] ",3F15.6)')COORX(IX,1), COORY(IX,1), 
     +    COORZ(IX,1)
        WRITE(8,'(" [---] ",3F15.6)')COORX(NNI,1), COORY(NNI,1), 
     +    COORZ(NNI,1)        
      ENDDO            
      WRITE(8,*)

      COORX(:,2)=COORX(:,1)
      COORY(:,2)=COORY(:,1)
      COORZ(:,2)=COORZ(:,1)


      IZ=INT((DOMZ(2)-DOMZ(1))/DDL,4)
      IY=INT((DOMY(2)-DOMY(1))/DDL,4)
      IX=INT((DOMX(2)-DOMX(1))/DDL,4)

      IXMAX=IX+20
      IYMAX=IY+20
      IZMAX=IZ+20

      END SUBROUTINE CYLIND5
!!----------------------------END CYLIND5-----------------------------!!



!!--------------------------------SRI2--------------------------------!!
      SUBROUTINE SRI2(LNODE,NODN,NODEID,NWALLID,COORX,COORY,COORZ,
     +  DDR,FSSRCH,H0)
      USE MLPGKINE
      USE NEIGHNODES
      !INCLUDE 'COMMON.F'

      INTEGER(KIND=4),INTENT(IN)::LNODE,NODN,NODEID(-7:NODN)
      INTEGER(KIND=4),INTENT(INOUT)::NWALLID(LNODE,4)
      INTEGER(KIND=4)::IDCHANGE,IDPREVIOUS,NUMTOTAL
      INTEGER(KIND=4)::NUMNE,KK,KKIJ,IJJJ,IDMIN
      REAL(KIND=8),INTENT(IN)::DDR(NODEID(0)),H0
      REAL(KIND=8),INTENT(IN)::COORX(LNODE),COORY(LNODE)
      REAL(KIND=8),INTENT(IN)::COORZ(LNODE)
      REAL(KIND=8)::RIAV,RIAV2,ZLIM
      REAL(KIND=8)::XX,YY,ZZ,DR,DRMIN
      LOGICAL,INTENT(INOUT)::FSSRCH(NODEID(0))

C       WRITE(8,*)'TRIAL!@1@1'
C       DO I=1,NODEID(0)
C         IF(FSSRCH(I)) CYCLE
C         WRITE(8,*)NODEID(I),FSSRCH(I)
C       ENDDO

      IDCHANGE=0
      ZLIM=0.8D0*H0
      IJJJ=NODEID(-1)+NODEID(-7)-NODEID(-2)
      DO II=NODEID(-7)+1,NODEID(-1)        
        IG=II-NODEID(-7)+IJJJ
        KK=NODEID(II)        
        IDPREVIOUS=NWALLID(II,2)

        NWALLID(II,2)=0

        NUMNE=0
        NUMTOTAL=NLINK(II)%I(0)      
        RIAV=1.5D0*DDR(II)
        RIAV2=0.2D0*DDR(II)
        RIAV3=1.2D0*DDR(II)

        DO IJ=1,NUMTOTAL
          INDJ=NLINK(II)%I(IJ)
          KKIJ=NODEID(INDJ)
          XX=(COORX(INDJ)-COORX(II))
          YY=(COORY(INDJ)-COORY(II))
          ZZ=(COORZ(INDJ)-COORZ(II))
          DR=SQRT(XX**2+YY**2+ZZ**2)/(RIAV)

          IF(DR.LE.1D0)THEN
            IF((KKIJ.EQ.0).OR.(KKIJ.EQ.4))THEN
              !XX=(COORX(INDJ,2)-COORX(II,2))*SSX(II)
              !YY=(COORY(INDJ,2)-COORY(II,2))*SSY(II)
              ZZ=(COORZ(INDJ)-COORZ(II))!*SSZ(II)
              !DR=SQRT(XX**2+YY**2+ZZ**2)              
              IF((ZZ.GE.-RIAV2).AND.(ZZ.LE.RIAV3)) NUMNE=NUMNE+1
            ENDIF
          ENDIF          
        ENDDO

        IF((NUMNE.LT.1).AND.(COORZ(II).GT.ZLIM))THEN
          NWALLID(II,2)=-10
          FSSRCH(II)=.FALSE.
          FSSRCH(IG)=.FALSE.
          !NWALLID(IG,2)=-10          
          WRITE(1616,'(2F15.6,I10)')COORZ(II)-H0,COORZ(II),NUMNE      
        ELSE
          NWALLID(II,2)=1
          FSSRCH(II)=.TRUE.
          FSSRCH(IG)=.TRUE.
        ENDIF

C         IF(COORZ(II,2).GT.1D0)THEN
C           NWALLID(II,2)=-10
C           FSSRCH(II)=.FALSE.
C           FSSRCH(IG)=.FALSE.
C           !NWALLID(IG,2)=-10          
C           WRITE(1616,'(2F15.6,2I10)')COORZ(II,2)-1D0,COORZ(II,2),NUMNE,
C      +      FSSRCH(IG)
C         ELSE
C           NWALLID(II,2)=1
C           FSSRCH(II)=.TRUE.
C           FSSRCH(IG)=.TRUE.
C         ENDIF

        IF(IDPREVIOUS.NE.NWALLID(II,2))THEN 
          IDCHANGE=IDCHANGE+1

          !! VALUES OF U* FOR NODES CHANGED FROM -10 TO 1
          IF(IDPREVIOUS.EQ.-10)THEN
            DRMIN=1D10
            IDMIN=0
            DO IJ=1,NUMTOTAL
              INDJ=NLINK(II)%I(IJ)
              KKIJ=NODEID(INDJ)
              XX=(COORX(INDJ)-COORX(II))
              YY=(COORY(INDJ)-COORY(II))
              ZZ=(COORZ(INDJ)-COORZ(II))
              DR=SQRT(XX**2+YY**2+ZZ**2)

              IF((DR.LT.DRMIN).AND.((KKIJ.EQ.0).OR.(KKIJ.EQ.4)))THEN
                DRMIN=DR
                IDMIN=INDJ
              ENDIF
            ENDDO            
            UX(II,1:2)=UX(IDMIN,1:2)
            UY(II,1:2)=UY(IDMIN,1:2)
            UZ(II,1:2)=UZ(IDMIN,1:2)
          ENDIF
        ENDIF

      ENDDO

      WRITE(8,*)'[SRI] IDCHANGE = ',IDCHANGE      

C       OPEN (1660, FILE ='Export//VEL_UP.DAT', STATUS='UNKNOWN')
C       DO I=NODEID(-7)+1, NODEID(-1)
C         WRITE(1660,'(3F16.8,2I4)')coorx(I,1),coory(I,1),coorZ(I,1),
C      +  nodeid(i),nwallid(i,2)!UX(I,2),UY(I,2),UZ(1,2)
C       END DO
C       CLOSE (1660)      

C       OPEN (1667, FILE ='Export//bottom.DAT', STATUS='UNKNOWN')
C       DO I=NODEID(-3)+1,NODEID(-4)
C         WRITE(1667,'(3F16.8,2I4)')coorx(I,1),coory(I,1),coorZ(I,1),
C      +  nodeid(i),nwallid(i,2)!UX(I,2),UY(I,2),UZ(1,2)
C       END DO
C       CLOSE(1667)

C       WRITE(8,*)'TRIAL-1-1-1-1'
C       DO I=1,NODEID(0)
C         IF(FSSRCH(I)) CYCLE
C         WRITE(8,*)NODEID(I),FSSRCH(I)
C       ENDDO
C       STOP

      END SUBROUTINE SRI2
!!--------------------------------SRI2--------------------------------!!      



!!----------------------------FIXCYLINDER-----------------------------!!
      SUBROUTINE FIXCYLINDER(BNDNP,BNDXY)
      USE COMMONMOD
      IMPLICIT NONE
      !INCLUDE 'COMMON.F'

      INTEGER(KIND=4),INTENT(IN)::BNDNP
      REAL(KIND=8),INTENT(IN)::BNDXY(BNDNP,3)

      INTEGER(KIND=4)::II,IK      

      DO II=NODEID(-7)+1,NODEID(-1)
        IK=II-NODEID(-2)
        COORX(II,1:2)=BNDXY(IK,1)
        COORY(II,1:2)=BNDXY(IK,2)
        COORZ(II,1:2)=BNDXY(IK,3)
      ENDDO

C       RADIA=0.157d0
C       IZ=20            
C       DDL=1.D0/IZ      
C       IY= 32       

C       I_OX=211
C       I_OY=NINT(IY/2.D0)
C       I_OZ=IZ*1.5
      
C       NNI=NODEID(-7)

C       ! CYLINDER NODES         
C       DO  K= 0,I_OZ
C         DO tt=0,((2*pi)-(pi/10.0D0)),(pi/10.0D0)!I_OR
C           IF (K .EQ. 0) THEN       
C             XX=I_OX*DDL+RADIA*COS(tt)
C             YY=I_OY*DDL+RADIA*SIN(tt)
C             !Z1=DDL/2.9
C             Z1=0D0
C             NNI=NNI+1
C             COORX(NNI,1:2)=XX
C             COORY(NNI,1:2)=YY
C             COORZ(NNI,1:2)=Z1
          
C           ELSE 
C             XX=I_OX*DDL+RADIA*COS(tt)
C             YY=I_OY*DDL+RADIA*SIN(tt)
C             !Z1=DDL/2.9+(K*DDL)
C             Z1=K*DDL
C             NNI=NNI+1
C             COORX(NNI,1:2)=XX
C             COORY(NNI,1:2)=YY
C             COORZ(NNI,1:2)=Z1
            
C           END IF 
          
C         END DO
C       END DO       

      RETURN
      END SUBROUTINE FIXCYLINDER
!!--------------------------END FIXCYLINDER---------------------------!!      



!!---------------------------FINDCYLINFORCE---------------------------!!
      SUBROUTINE FINDCYLINFORCE(LNODE,NODN,NODEID,NWALLID,COORZ,P,
     +  SNX,SNY,SNZ,FX,FY,FZ)
      IMPLICIT NONE

        INTEGER(KIND=4),INTENT(IN)::LNODE,NODN,NODEID(-7:NODN)
        INTEGER(KIND=4),INTENT(IN)::NWALLID(LNODE,4)
        REAL(KIND=8),INTENT(IN)::P(LNODE),SNX(LNODE),SNY(LNODE)
        REAL(KIND=8),INTENT(IN)::SNZ(LNODE),COORZ(LNODE)
        REAL(KIND=8),INTENT(OUT)::FX,FY,FZ

        INTEGER(KIND=4)::NL,NR,IX
        REAL(KIND=8)::CYLR,PI,CYLL,DA

        PI=ATAN(1D0)*4D0
        CYLR=0.11D0
        CYLL=1.05D0
        NL=25
        NR=16

        DA=CYLL*2D0*PI*CYLR/2D0/NR/(NL-1)

        FX=0D0
        FY=0D0
        FZ=0D0
        DO IX=NODEID(-7)+1,NODEID(-1)
          IF(NWALLID(IX,2).EQ.-10) CYCLE

          IF(ABS(COORZ(IX)-0D0).LT.1D-4)THEN
            FX=FX+P(IX)*SNX(IX)
            FY=FY+P(IX)*SNY(IX)
            FZ=FZ+P(IX)*SNZ(IX)
          ELSE
            FX=FX+2D0*P(IX)*SNX(IX)
            FY=FY+2D0*P(IX)*SNY(IX)
            FZ=FZ+2D0*P(IX)*SNZ(IX)
          ENDIF
        ENDDO
        FX=FX*DA
        FY=FY*DA
        FZ=FZ*DA
        
      END SUBROUTINE FINDCYLINFORCE
!!-------------------------END FINDCYLINFORCE-------------------------!!      