	Subroutine DumpEverything(sometext)


	implicit none
	
	include "Assemb.inc"
	include "Monit.inc"
	include "Tangent.inc"
	include "Newton.inc"

	integer*4 k,i,iunit,kcur,j,L,Lsite
	character*128 sometext
	iunit = 12


	write(iunit,*)' '
	write(iunit,*)' '
	write(iunit,*)' '
	write(iunit,*)' $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
	write(iunit,*)' $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
	write(iunit,*)' $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
	write(iunit,*)' $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
	write(iunit,*)' Dump everything: ',sometext	
! --------------------------------------
         write(iunit,*)' '
         write(iunit,*)'------------------------------------------------------'
         write(iunit,*)'------------------------------------------------------'
         write(iunit,*)' Variable list '
         do 11 i = 1,nvar
         write(iunit,15)I,VN1(I),VN2(I)
15       FORMAT(I3,2X,A2,A16)
11       continue
         write(iunit,*)'NEQ=',NEQ,'  NVAR=',NVAR
         write(iunit,*)
         write(iunit,*)'VARIANCE =',NVAR-NEQ
         write(iunit,*)
         write(iunit,*)'MONITOR PARAMETERS    Delta X for monitors'
         if(nvar-neq.gt.0)then
         do 12 i = 1,nvar-neq
         write(iunit,16)MON(I),VN1(MON(I)),VN2(MON(I)),Deltax(i)*NSTEP
16       FORMAT(I3,2X,A2,A16,F12.3)
12       continue
!         write(iunit,*)' DELTA Xi FOR MONITOR PHASES'
!         write(iunit,16)(DELTAX(I)*NSTEP,I=1,NVAR-NEQ)
!16       FORMAT(10F12.3)
         else
         write(iunit,*)' No monitor parameters (variance <= 0)'
         endif
         write(iunit,*)'NUMBER OF FINITE DIFF. ITERATIONS=',NSTEP
         write(iunit,*)
         write(iunit,*)'INITIAL MINERAL DATA'
         write(iunit,17)TC,PB
17       format(' T START (DEG C)=',F10.2,'   P START=',F10.2)
	if(PFluidSwitch.gt.0)write(iunit,18)PFluid
18	format(' PFluid = ',F10.2)
90    	CONTINUE



350       CONTINUE
	Do 390 kCur = 1,numPh
	k = asmCurrent(kCur)
	write(iunit,*)
          write(iunit,200)MINREC(K),SITMUL(K),PHNAME(K)
200       FORMAT(' ',I4,F8.1,4X,A32)
          write(iunit,'(1X,A15,1F10.3)')'MOLAR VOLUME =',VMOL(K)
          IF(IMASS.EQ.1) then
             write(iunit,507)MP0(K),VP0(K)
          endif
507       FORMAT(' MOLES =',F10.5,' VOL =',F10.3)
          write(iunit,501)(phCoName(k,j),J=1,numPhCo(K))
501       FORMAT('          ',12(A4,6X))
          write(iunit,504)(xPhCo(k,j),J=1,numPhCo(K))
504       FORMAT('  COMP ',12(F10.4))
          write(iunit,505)(Dexp(lnAct(k,j)),J=1,numPhCo(K))
505       FORMAT(' Activ  ',12(E10.3))
          write(iunit,506)(lnAct(k,j),J=1,numPhCo(K))
506       FORMAT('  Ln(a) ',12(F10.5))
          write(iunit,514)(hPhCoZero(k,j),J=1,numPhCo(K))
514       Format(' hPhCoZero',12E15.7)
          write(iunit,515)(HATTP(k,j),J=1,numPhCo(K))
515       Format(' HATTP',12E15.7)
          write(iunit,509)(sPhCoZero(k,j),J=1,numPhCo(K))
509       Format(' SZERO',12F15.5)
          write(iunit,510)(SATTP(k,j),J=1,numPhCo(K))
510       Format(' SATTP',12F15.5)
          write(iunit,508)(vPhCoZero(k,j),J=1,numPhCo(K))
508       Format(' VZERO',12F15.5)
          write(iunit,513)(VATTP(k,j),J=1,numPhCo(K))
513       Format(' VATTP',12F15.5)
          write(iunit,516)(GATTP(k,j),J=1,numPhCo(K))
516       Format(' GATTP',12E15.7)
          write(iunit,*)'dudTPX ARRAY:'
          write(iunit,*)'     dudT   .    dudP   .    dudX2  .    dudX3  .    dudX4...'
	do 512 j = 1,numPhCo(K)
        write(iunit,511)dudTPX(k,j,1),dudTPX(k,j,2),(dudTPX(k,j,L),L=3,numPhCo(K)+1)
511     FORMAT(15E12.5)
512   	continue
!         WRITE(iout,408)' ACP',(aPhCoCp(k,j),J=1,numPhCo(K))
!         WRITE(iout,408)' BCP',(bPhCoCp(k,j),J=1,numPhCo(K))
!         WRITE(iout,408)' CCP',(cPhCoCp(k,j),J=1,numPhCo(K))
!         WRITE(iout,408)' DCP',(dPhCoCp(k,j),J=1,numPhCo(K))
!         WRITE(iout,408)' ECP',(ePhCoCp(k,j),J=1,numPhCo(K))
!         WRITE(iout,408)' FCP',(fPhCoCp(k,j),J=1,numPhCo(K))
!         WRITE(iout,408)' GCP',(gPhCoCp(k,j),J=1,numPhCo(K))
408       format(' ',a5,6E10.3)
551       CONTINUE
390       CONTINUE


! --------------------------------------
	call WriteOutJacobian()


! --------------------------------------
	WRITE(iunit,*)'   '
	WRITE(iunit,*)'***** Check ALLX variable arrays ********'
	WRITE(iunit,*)'Name, ALLX(1),ALLX(2),ALLX(3),ALLX(4),ALLX(5),ALLX(6)'
	WRITE(iunit,*)'  Name                Current     Starting    Contour  User ref    previous fdiff    FRACTL(Y/N)'
	do 1441 I=1,nvar
1441    WRITE(iunit,1442)I,VN1(I),VN2(I),(ALLX(j,i),J=1,5),allFractl(i)
1442    FORMAT(I3,1X,A2,A16,5F12.4,I5)

	write(iunit,*)' '
	write(iunit,*)' '
	write(iunit,*)'xPhCoNewP/xPhCoLast/xPHco'
	write(iunit,669)((phCoName(k,j),j=1,numPhCo(k)),k=1,numPhMIF)
669	format(200(2x,A8))
	write(iunit,668)((xPhCoNewP(k,j),j=1,numPhCo(k)),k=1,numPhMIF)
	write(iunit,668)((xPhCoLast(k,j),j=1,numPhCo(k)),k=1,numPhMIF)
668	format(200F10.5)
	write(iunit,668)((xPhCo(k,j),j=1,numPhCo(k)),k=1,numPhMIF)
	write(iunit,*)' '

	do 710 k = 1,numPhMIF
      	if(numSiteAtom(K).gt.0)then
      		call XtoSite(K)
		write(iunit,*)' '
		write(iunit,*)PhName(k)
      		write(iunit,706)(SiteAtomName(k,Lsite),Lsite = 1,numSiteAtom(K))
706   		format(4x,20A12)
      		write(iunit,703)(SiteAtom(k,Lsite),Lsite = 1,numSiteAtom(K))
703   		format (20F12.6)
      		endif
710	continue

!	write(iunit,*)'NewtonCounter = ',newtoncounter
	write(iunit,*)' NewtonStep = ',NewtonStep

	write(iunit,*)TC,PB
	write(*,*)TC,PB
	write(*,*)' Dump everything: ',sometext	
	write(iunit,*)' Dump everything: ',sometext	
	PAUSe ' DumpEverything is done.... examine and hit return when happy (or not)'

	return
	end
	