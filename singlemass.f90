subroutine singlemass(local_tho,iso)
use utilmf
implicit none

integer::i,j,star_cnt
integer, intent(in)::iso, local_tho
 integer,parameter:: MyLong = selected_int_kind (12)
 integer (kind=MyLong) :: k, n_star




	star_cnt=0

!													1

	do i=1, n_pad_2
		if(log_age_2(i) .eq. log_age_iso(iso))then
			star_cnt=star_cnt+1
		else if (log_age_2(i) > log_age_iso(iso)) then
			exit
		end if
	end do
	n_t=star_cnt					!Number of stars in each isochrone

allocate (m(0:n_t),dm(0:n_t),dlm(0:n_t),n(0:n_t),phi(0:n_t),wght(0:n_t),wght_tmp(0:n_t),ntmp(0:n_t),total_mass(0:n_t))
allocate (lumino(1:n_t), total_lum(1:n_t),total_mass_IG(1:n_t),total_lum_IG(1:n_t))
allocate (phi_prim(0:n_t),PHI_IG(0:n_t),NIG(0:n_t))
allocate (lU(1:n_t),lB(1:n_t),lV(1:n_t),lI(1:n_t),lJ(1:n_t),lK(1:n_t),lR(1:n_t))


	lU(:)=0.0;lB(:)=0.0;lV(:)=0.0;lJ(:)=0.0;lI(:)=0.0;lK(:)=0.0;lR(:)=0.0



	m(0)=m_ini(0)
	n(0)=0.0
	phi(0)=0.0
	wght(0)=0
	wght_tmp(0)=0
	dm(0)=0.0
	dlm(0)=0.0

!													2
	j=1
	do i=1, n_pad_2
		if(log_age_2(i) .eq. log_age_iso(iso))then


			lumino(j)=10.0**(log_L_2(i))

	          lB(j)=10.0**((Bsun-mag_B_2(i))/2.5)
     	     lV(j)=10.0**((Vsun-mag_V_2(i))/2.5)
     	     lR(j)=10.0**((Rsun-mag_R_2(i))/2.5)
     	     lK(j)=10.0**((Ksun-mag_K_2(i))/2.5)
     	     lI(j)=10.0**((Isun-mag_I_2(i))/2.5)
     	     lJ(j)=10.0**((Jsun-mag_J_2(i))/2.5)
     	     lU(j)=10.0**((Usun-mag_U_2(i))/2.5)

			j=j+1
		else if (log_age_2(i) > log_age_iso(iso)) then
			exit
		end if
	end do





	do j=1, n_t				!Total Mass Steps in Each Isochrone
		m(j)=m_ini_2(j)
		dm(j)=m(j)-m(j-1)
		dlm(j)=log10(m(j))-log10(m(j-1))
		n(j)=0.0
		ntmp(j)=0.0
		wght_tmp(j)=0
		phi(j)=0.0
		phi_prim(j)=0.0
		PHI_IG(j)=0.0
		NIG(j)=0.0


	end do
!													3
	do j=1,25
		new_n(j)=0.0
		new_nb(j)=0.0
		new_n_igimf(j)=0.0
	end do
!__________________________________________________ 2001 ________________________________ 
	do j=1, n_t


		if(m(j)>m(0) .and. m(j)<=0.5) then
			phi(j)=0.432*(m(j)**(-1.3))
!			phi(j)=0.236686*(((m(j)+m(j-1))/2.0)**(-1.3))
		else if (m(j)>0.5 .and. m(j)<= m(n_t))then
			phi(j)=0.216*(m(j)**(-2.3))
!			phi(j)=0.118343*(((m(j)+m(j-1))/2.0)**(-2.3))
		end if

!_______________________________________________________ 1998____________________________
!		if(m(j)>m(0) .and. m(j)<=0.5) then
!			phi(j)=0.48*(m(j)**(-1.5))
!		else if(m(j)>0.5 .and. m(j)<=1.0) then
!			phi(j)=0.30*(m(j)**(-2.2))
!		else if (m(j)>1.0 .and. m(j)<= m(n_t))then
!			phi(j)=0.30*(m(j)**(-2.7))
!		end if
!________________________________________________________________________________________


!													4


		n(j)=phi(j)*dm(j)*(SFR(iso)*ddt(iso))

		lB_t(iso)=lB_t(iso)+lB(j)*n(j)
		lV_t(iso)=lV_t(iso)+lV(j)*n(j)
		lR_t(iso)=lR_t(iso)+lR(j)*n(j)
		lK_t(iso)=lK_t(iso)+lK(j)*n(j)
		lI_t(iso)=lI_t(iso)+lI(j)*n(j)
		lJ_t(iso)=lJ_t(iso)+lJ(j)*n(j)
		lU_t(iso)=lU_t(iso)+lU(j)*n(j)
		
		call igimf(iso,j)


		total_mass(j)=n(j)*(m(j))
		total_lum(j)=n(j)*lumino(j)


		total_mass_IG(j)=NIG(j)*m(j)
		total_lum_IG(j)=NIG(j)*lumino(j)

		NIMF_f(j)=NIMF_f(j)+n(j)
		NIG_f(j)=NIG_f(j)+NIG(j)

		ntmp(j)=n(j)
		wght(j)=n(j)+wght(j-1)
		wght_tmp(j)=wght(j)


	end do
!______________________________________________________________________________ 


!													6
	iso_mass(iso)=sum(total_mass)

	iso_lum(iso)=sum(total_lum)

	iso_mass_IG(iso)=sum(total_mass_IG)
	iso_lum_IG(iso)=sum(total_lum_IG)


	call mremnant(iso,m(n_t))
	call IGIMFremnant(iso,m(n_t))

	m_glxy(iso)=iso_mass(iso)+total_rem(iso)


	n_tot=wght(n_t)
	n_star=n_tot*1				!number of total stars in each isochrone


end subroutine
