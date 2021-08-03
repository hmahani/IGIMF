subroutine beginning(iso_number)
use utilmf
implicit none

integer::i, j, iso_cnt,k, cnt
integer, intent(out)::iso_number
real(8)::l_tmp,U_tmp,B_tmp,V_tmp,R_tmp,I_tmp,J_tmp,K_tmp

	z_met(0)=0.02
	log_age(0)=6.55
	m_ini(0)=0.080
	m_act(0)=0.080
	log_L(0)=-1.65
	log_T(0)=3.409
	log_G(0)=3.64
	iso_cnt=0
	m_tmp(0)=0.088
	cnt=1
	log_mass(0)=log10(m_tmp(0))
	Mecl(0)=5.0

!											1
	open(55,file='scalo86.txt')
	do i=1, 25				
		read(55,*)scalo_m(i),scalo_n(i)
	end do
	close(55)

	open(56,file='scalo86ERROR.txt')
	do i=1, 25				
		read(56,*)scalo_a(i),scalo_b(i),scalo_err(i)
	end do
	close(56)


	tho(1)=1.55*10.**9.; tho(2)=2.2*10.**9.; tho(3)=2.8*10.**9.; tho(4)=3.25*10.**9.; tho(5)=3.8*10.**9.
     tho(6)=4.9*10.**9.; tho(7)=6.2*10.**9.; tho(8)=8.0*10.**9.;  tho(9)=10.5*10.**9.; tho(10)=15.0*10.**9.
     tho(11)=23.0*10.**9.; tho(12)=50.0*10.**9.; tho(13)=1.e4*10.**9. !10^5 tho khaili bozorg mishod va hameye SFR ha ro sefr midad
     tho(26)=-1.2*10.**9.; tho(25)=-1.55*10.**9.; tho(24)=-2.2*10.**9.; tho(23)=-2.8*10.**9.; tho(22)=-3.25*10.**9. 
     tho(21)=-3.8*10.**9.; tho(20)=-4.9*10.**9.;  tho(19)=-6.2*10.**9.; tho(18)=-8.0*10.**9.; tho(17)=-10.5*10.**9.
     tho(16)=-15.0*10.**9.; tho(15)=-23.0*10.**9.; tho(14)=-50.0*10.**9.

print*,'START'


!											2
	open(51,file='padova02.txt')
	do i=1, n_pad				!number of Padova objects
		read(51,*)z_met(i), log_age(i), m_ini(i), m_act(i), log_L(i), log_T(i), log_G(i),mbol(i),&
		&magu(i),magb(i),magv(i),magr(i),magi(i),magj(i),magh(i),magk(i),ini_imf(i),stg(i)
		if(log_age(i)>log_age(i-1))then
			iso_cnt=iso_cnt+1
		end if
	end do
	close(51)
	isochrone=iso_cnt
	iso_number=iso_cnt

	magu(0)=magu(1)-0.05;magb(0)=magb(1)-0.05;magv(0)=magv(1)-0.05;magr(0)=magr(1)-0.05
	magi(0)=magi(1)-0.05;magj(0)=magj(1)-0.05;magk(0)=magk(1)-0.05

allocate (log_age_iso(0:isochrone),m_max_pad(0:isochrone),tt(1:isochrone),ddt(1:isochrone),SFR(0:isochrone),total_rem(1:isochrone))
allocate (iso_mass(1:isochrone),m_glxy(1:isochrone),iso_lum(1:isochrone),total_rem_IG(1:isochrone),SFR1(0:isochrone))
allocate (iso_mass_IG(1:isochrone),iso_lum_IG(1:isochrone))
allocate (lU_t(1:isochrone),lB_t(1:isochrone),lV_t(1:isochrone),lI_t(1:isochrone))
allocate (lJ_t(1:isochrone),lK_t(1:isochrone),lR_t(1:isochrone))
allocate (lUIG_t(1:isochrone),lBIG_t(1:isochrone),lVIG_t(1:isochrone),lIIG_t(1:isochrone))
allocate (lJIG_t(1:isochrone),lKIG_t(1:isochrone),lRIG_t(1:isochrone))


	log_age_iso(0)=6.55
	m_glxy(:)=0.0
	total_rem(:)=0.0
	total_rem_IG(:)=0.0
	iso_mass(:)=0.0
	lU_t(:)=0.0;lB_t(:)=0.0;lV_t(:)=0.0;lJ_t(:)=0.0;lI_t(:)=0.0;lK_t(:)=0.0;lR_t(:)=0.0
	lUIG_t(:)=0.0;lBIG_t(:)=0.0;lVIG_t(:)=0.0;lJIG_t(:)=0.0;lIIG_t(:)=0.0;lKIG_t(:)=0.0;lRIG_t(:)=0.0

!											3
	j=1
	do i=1, n_pad				
		if(log_age(i)>log_age(i-1))then
			log_age_iso(j)=log_age(i)
			m_max_pad(j-1)=m_ini(i-1)
			j=j+1
		end if
	end do

	m_max_pad(j-1)=0.9956

!											4
	do i=1,isochrone
		tt(i)=10.0**(log_age_iso(i))
		if (log_age_iso(i)>6.61)then
			ddt(i)=(10.0**(log_age_iso(i)))-(10.0**(log_age_iso(i-1)))
		else if (log_age_iso(i)<6.61)then
			ddt(i)=10.0**(log_age_iso(i))
		end if
!write(*,*)i,m_max_pad(i)
	end do

!											5
	do j=1, n_grp
!!		m_tmp(j)=m_tmp(0)+mass_stp*j
		log_mass(j)=log_mass(0)+j*mass_stp
		m_tmp(j)=10.0**log_mass(j)
		dm_tmp(j)=m_tmp(j)-m_tmp(j-1)
		dlm_tmp=log10(m_tmp(j))-log10(m_tmp(j-1))
	end do

!											6
	cnt=1
	open(52,file='padova02_modify.txt')
	do k=1, isochrone
		do j=1, n_grp
			do i=1, n_pad
				if(log_age_iso(k).eq.log_age(i))then
					if(m_tmp(j)<=m_max_pad(k))then


						if (m_tmp(j)<m_ini(i-1).and.m_tmp(j)<m_ini(i))then
						l_tmp=log_L(i)
	   					U_tmp=magu(i)
	   					B_tmp=magb(i)
	   					V_tmp=magv(i)
	   					R_tmp=magr(i)
	   					I_tmp=magi(i)
	   					J_tmp=magj(i)
	   					K_tmp=magk(i)
							write(52,"(F8.2,F11.5,F11.4,F10.3,F10.3,F10.3,F10.3,F10.3,F10.3,F10.3)")&
								&log_age(i),m_tmp(j),l_tmp,U_tmp,B_tmp,V_tmp,R_tmp,I_tmp,J_tmp,K_tmp
							cnt=cnt+1
							exit


						else if (m_tmp(j)>m_ini(i-1).and.m_tmp(j)<=m_ini(i))then
	   					l_tmp=log_L(i)+((log_L(i)-log_L(i-1))/(log10(m_ini(i))-log10(m_ini(i-1))))*(log10(m_tmp(j))-log10(m_ini(i)))

	   					U_tmp=magu(i)+((magu(i)-magu(i-1))/((m_ini(i))-(m_ini(i-1))))*((m_tmp(j))-(m_ini(i)))
	   					B_tmp=magb(i)+((magb(i)-magb(i-1))/((m_ini(i))-(m_ini(i-1))))*((m_tmp(j))-(m_ini(i)))
	   					V_tmp=magv(i)+((magv(i)-magv(i-1))/((m_ini(i))-(m_ini(i-1))))*((m_tmp(j))-(m_ini(i)))
	   					R_tmp=magr(i)+((magr(i)-magr(i-1))/((m_ini(i))-(m_ini(i-1))))*((m_tmp(j))-(m_ini(i)))
	   					I_tmp=magi(i)+((magi(i)-magi(i-1))/((m_ini(i))-(m_ini(i-1))))*((m_tmp(j))-(m_ini(i)))
	   					J_tmp=magj(i)+((magj(i)-magj(i-1))/((m_ini(i))-(m_ini(i-1))))*((m_tmp(j))-(m_ini(i)))
	   					K_tmp=magk(i)+((magk(i)-magk(i-1))/((m_ini(i))-(m_ini(i-1))))*((m_tmp(j))-(m_ini(i)))

							write(52,"(F8.2,F11.5,F11.4,F10.3,F10.3,F10.3,F10.3,F10.3,F10.3,F10.3)")&
								&log_age(i),m_tmp(j),l_tmp,U_tmp,B_tmp,V_tmp,R_tmp,I_tmp,J_tmp,K_tmp
							cnt=cnt+1
							exit
						end if
					end if
				end if
			end do
		end do
	end do
	close(52)
	n_pad_2=cnt-1


allocate (log_age_2(1:n_pad_2),m_ini_2(1:n_pad_2),log_L_2(1:n_pad_2))
allocate (mag_u_2(1:n_pad_2),mag_b_2(1:n_pad_2),mag_v_2(1:n_pad_2),mag_i_2(1:n_pad_2),mag_j_2(1:n_pad_2),mag_k_2(1:n_pad_2))
allocate (mag_r_2(1:n_pad_2))

!											7
	open(53,file='padova02_modify.txt')
	do i=1, n_pad_2
		read(53,*) log_age_2(i),m_ini_2(i),log_L_2(i),mag_u_2(i),mag_b_2(i),mag_v_2(i),mag_r_2(i),mag_i_2(i)&
				&,mag_j_2(i),mag_k_2(i)
	end do
	close(53)
		
			

end subroutine
