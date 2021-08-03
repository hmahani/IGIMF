program main
use utilmf
implicit none

integer::i,j
integer::iso_number
	character(len=70)::fn
	integer, parameter::numfiles=13
	integer, parameter::outunit=44
!real(8),dimension(1:30)::newlogm
newlogm(0)=-1.1969
tst1=0.0

!	do i=1,25
!		!newlogm(i)=newlogm(0)+i*0.1
!	end do


!													1

	newlogm(1)=-1.06 ; newlogm(2)=-0.96 ; newlogm(3)=-0.85 ; newlogm(4)=-0.75 ; newlogm(5)=-0.65 
	newlogm(6)=-0.54 ; newlogm(7)=-0.44 ; newlogm(8)=-0.35 ; newlogm(9)=-0.27 ; newlogm(10)=-0.21 
	newlogm(11)=-0.14 ; newlogm(12)=-0.08 ; newlogm(13)=-0.01 ; newlogm(14)=0.07 ; newlogm(15)=0.16 
	newlogm(16)=0.27 ; newlogm(17)=0.40 ; newlogm(18)=0.54 ; newlogm(19)=0.72 ; newlogm(20)=0.90 
	newlogm(21)=1.08 ; newlogm(22)=1.26 ; newlogm(23)=1.43 ; newlogm(24)=1.62 ; newlogm(25)=1.795 


	call beginning(iso_number)
	open(72,file='fit.dat')
	open(73,file='SFR.dat')
	open(74,file='MASS&LUMINOSITY.dat')
	open(75,file='COLORs.dat')
	open(76,file='TOTAL_LUMINOSITY.dat')
	open(77,file='DATAPLOT.txt')

write(72,*)'i   x3  b3  x3_IGIMF  b3_IGIMF x4  b4  x4_IGIMF  b4_IGIMF boost boostIG chi2  chi2IG'
write(74,*)'   M_IMF Mec_IMF   rem_IMF  L_IMF   M/L IMF   M_IGIMF Mec_IGIMF  rem_IGIMF  L_IGIMF    M/L IGIMF'
write(75,*)'tho clrB clrV clrI clrJ clrK clrR clrU clrB_IG clrV_IG clrI_IG clrJ_IG clrK_IG clrR_IG clrU_IG b-r v-k i-k b-iIG '
write(76,*)'tho lB_tt lV_tt lI_tt lJ_tt lK_tt lR_tt lU_tt  lBIG_tt lVIG_tt lIIG_tt lJIG_tt lKIG_tt lRIG_tt lUIG_tt '


!													2

!iso_number=71
write(*,*)iso_number
	do i=1,26
		call zero()
		do j=1,iso_number	
			SFR1(j)=(1./tho(i))*(1./(1.-exp(-1.*(tt_u/tho(i)))))*  exp(-1.*(tt_u-tt(j))/tho(i))*(4.0*10.**11.) 
			SFR1(0)=SFR1(1) 

			SFR(j)=((SFR1(j)+SFR1(j-1))/2.0)
			write(73,*)j, SFR(j), i

          !constructing IGIMF

!													3       
            if(SFR(j)>=1.0)then
                Bt =-1.*(-0.106*log10(SFR(j))+2.)
            endif
            if(SFR(j)<1.0)then
                Bt =-2.
            endif 
                !Bt=-2.  !normal IGIMF

            Mecl_max=8.5*(10.0**4.0)*(SFR(j))**0.75
!            Mecl_max=10.0**((0.746*log10(SFR(j)))+4.93)
 

            if(Mecl_max > 1E9)then
                Mecl_max=1E9
            endif 

            delta=(log10(Mecl_max)-log10(Mecl(0)))/(1.0*n_igimf)
write(*,*)j,SFR(j),Mecl_max

!													4
			call singlemass(i,j)
			call binary(i,j)

		end do

!													5

		lB_tt=sum(lB_t); lV_tt=sum(lV_t); lI_tt=sum(lI_t); lJ_tt=sum(lJ_t); 
		lK_tt=sum(lK_t); lR_tt=sum(lR_t); lU_tt=sum(lU_t); 

		lBIG_tt=sum(lBIG_t); lVIG_tt=sum(lVIG_t); lIIG_tt=sum(lIIG_t); lJIG_tt=sum(lJIG_t); 
		lKIG_tt=sum(lKIG_t); lRIG_tt=sum(lRIG_t); lUIG_tt=sum(lUIG_t); 

		clrB=Bsun-2.5*log10(lB_tt)
		clrV=Vsun-2.5*log10(lV_tt)
		clrI=Isun-2.5*log10(lI_tt)
		clrJ=Jsun-2.5*log10(lJ_tt)
		clrK=Ksun-2.5*log10(lK_tt)
		clrR=Rsun-2.5*log10(lR_tt)
		clrU=Usun-2.5*log10(lU_tt)


		clrB_IG=Bsun-2.5*log10(lBIG_tt)
		clrV_IG=Vsun-2.5*log10(lVIG_tt)
		clrI_IG=Isun-2.5*log10(lIIG_tt)
		clrJ_IG=Jsun-2.5*log10(lJIG_tt)
		clrK_IG=Ksun-2.5*log10(lKIG_tt)
		clrR_IG=Rsun-2.5*log10(lRIG_tt)
		clrU_IG=Usun-2.5*log10(lUIG_tt)

		boost=(log10(new_NIMF_f(25)/(newlogm(25)-newlogm(24))))-(scalo_n(25))	
		boostIG=(log10(new_n_igimf_f(25)/(newlogm(25)-newlogm(24))))-(scalo_n(25))	

!													6

			write(fn,fmt='(i0,a)') i, '.dat'
			open(unit=outunit, file=fn, form='formatted')

			do  j=1, 25
			write(outunit,*)newlogm(j), log10(new_NIMF_f(j)/(newlogm(j)-newlogm(j-1))),&
				 & log10(new_n_igimf_f(j)/(newlogm(j)-newlogm(j-1))),scalo_n(j)+boost,scalo_n(j)+boostIG
			end do

		chi2=0.0	
		chi2IG=0.0

			do  j=12, 25
	chi2=chi2+(((scalo_n(j)+boost)-(log10(new_NIMF_f(j)/(newlogm(j)-newlogm(j-1)))))**2.0)/&
			&((scalo_err(j))**2.0)
	chi2IG=chi2IG+(((scalo_n(j)+boost)-(log10(new_n_igimf_f(j)/(newlogm(j)-newlogm(j-1)))))**2.0)/&
			&((scalo_err(j))**2.0)
			end do		


!													7

!	call fit(1,x1,b1,x1_IGIMF,b1_IGIMF)
!	call fit(2,x2,b2,x2_IGIMF,b2_IGIMF)
	call fit(3,x3,b3,x3_IGIMF,b3_IGIMF)
	call fit(4,x4,b4,x4_IGIMF,b4_IGIMF)


	write(72,'(I2,F9.3,F9.3,F9.3,F9.3,F9.3,F9.3,F9.3,F9.3,F7.1,F7.1,F9.3,F9.3)')&
	&i,x3,b3,x3_IGIMF,b3_IGIMF,x4,b4,x4_IGIMF,b4_IGIMF,boost,boostIG,chi2,chi2IG



!	write(*,"(I2,3x,ES8.2,3x,ES8.2,3x,F9.3,3x,F9.3)")i,sum(m_glxy),sum(iso_lum),sum(m_glxy)/sum(iso_lum) 
	write(*,"(I2,3x,ES8.2,3x,ES8.2,3x,ES8.2)")i,sum(m_glxy),sum(iso_mass_IG),sum(total_rem_IG)
!	write(*,*)i,sum(m_glxy),sum(iso_mass_IG),sum(total_rem_IG)


	write(*,*)''
	write(74,"(ES8.2,3x,ES8.2,3x,ES8.2,3x,ES8.2,3x,F9.3,3x,ES8.2,3x,ES8.2,3x,ES8.2,3x,ES8.2,3x,F9.3)")&
	&sum(m_glxy),sum(iso_mass),sum(total_rem),sum(iso_lum),sum(m_glxy)/sum(iso_lum),(sum(iso_mass_IG)+sum(total_rem_IG)),&
	&sum(iso_mass_IG),sum(total_rem_IG),sum(iso_lum_IG),(sum(iso_mass_IG)+sum(total_rem_IG))/sum(iso_lum_IG)


11 FORMAT(I2,F9.3,F9.3,F9.3,F9.3,F9.3,F9.3,F9.3,F9.3,F9.3,F9.3,F9.3,F9.3,F9.3,F9.3,F9.3,F9.3,F9.3,F9.3,F9.3,F9.3)
	write(*,*)''
	write(75,11)i,clrB,clrV,clrI,clrJ,clrK,clrR,clrU,clrB_IG,clrV_IG,clrI_IG,clrJ_IG,clrK_IG,clrR_IG,clrU_IG,&
				&clrB-clrR,clrV-clrK,clrI-clrK,clrB_IG-clrR_IG,clrV_IG-clrK_IG,clrI_IG-clrK_IG

12 FORMAT(I2,3x,ES8.2,3x,ES8.2,3x,ES8.2,3x,ES8.2,3x,ES8.2,3x,ES8.2,3x,ES8.2,3x,ES8.2,3x,ES8.2,3x,ES8.2,3x,&
		&3x,ES8.2,3x,ES8.2,3x,ES8.2,3x,ES8.2)
	write(*,*)''
	write(76,12)i,lB_tt,lV_tt,lI_tt,lJ_tt,lK_tt,lR_tt,lU_tt,lBIG_tt,lVIG_tt,lIIG_tt,lJIG_tt,lKIG_tt,lRIG_tt,lUIG_tt
	write(77,"(ES8.2,3x,ES8.2)")sum(m_glxy),(sum(iso_mass_IG)+sum(total_rem_IG))

	end do

	close(72)
	close(73)
	close(74)
	close(75)
	close(76)
	close(77)



	open(71,file='results.dat')
	do j=1, n_grp
		write(71,*)log10(m_tmp(j)), log10(NIMF_f(j)/dm_tmp(j)),log10(NIMF_fb(j)/dm_tmp(j))
	end do
	close(71)

	write(*,*)''
	print*,'FINE'
	write(*,*)''
end program
