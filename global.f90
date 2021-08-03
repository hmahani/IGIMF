module utilmf
implicit none

integer, parameter:: n_pad=17384, prcnt=50, n_IGIMF=1000, n_grp=350!n_grp=300!!, n_pad_2=59717
real(8), parameter::tt_u=10.0**10.0, max_mass_imf=150.0, mass_stp=0.01,  FeH=0.0 
real(8), dimension(0:n_pad)::z_met, log_age, m_ini, m_act, log_L, log_T, log_G, mbol
real(8), dimension(0:n_pad)::magu,magb,magv,magr,magi,magj,magh,magk,ini_imf,stg
real(8), dimension(1:26)::tho
integer::isochrone, n_pad_2
real(8),dimension(:),allocatable::log_age_iso,m_max_pad,tt,ddt, SFR,SFR1, total_rem, iso_mass, m_glxy, iso_lum
real(8), dimension(:),allocatable:: m, dm, n, phi, ntmp, nb, nprim, dlm, total_mass, total_lum, lumino, total_lum_b
real(8), dimension(:),allocatable::total_mass_IG, total_lum_IG,iso_mass_IG,iso_lum_IG
integer::n_t
real(8):: n_tot
real(8), dimension(:),allocatable::wght,wght_tmp 
real(8),dimension(:),allocatable::s, sp, s2,prim, s2p, primp
real(8)::m1, m2, m_tot, m_prim
real(8),dimension(0:n_grp)::m_tmp, NIMF_f, NIMF_fb, dm_tmp, dlm_tmp, NIG_f
real(8),dimension(:),allocatable::log_age_2, log_L_2, m_ini_2
!real(8), dimension(0:n_pad_2)::log_age_2, log_L_2, m_ini_2
real(8),dimension(0:25)::new_n, new_NIMF_f,newlogm, new_nb, new_NIMF_fb, new_n_igimf, new_n_igimf_f
real(8), dimension(:),allocatable:: m_rem, dm_rem, n_rem, phi_rem, remnant,mass_remnant,log_mass_rem
real(8), dimension(0:n_grp)::log_mass
real(8), dimension(0:n_IGIMF)::Mecl, xx, m_max, Meclave, dMecl
real(8)::a1, a2, a3, k1, k2, k3, c1, Bt, delta, M_1, M_2, M_3, M_4, a_1, a_2, a_3, Mecl_max
real(8), dimension(:),allocatable::phi_prim, PHI_IG, NIG
real(8)::x1,x2,x3,b1,b2,b3,tst1,x4,b4
real(8)::x1_IGIMF,x2_IGIMF,x3_IGIMF,b1_IGIMF,b2_IGIMF,b3_IGIMF,x4_IGIMF,b4_IGIMF
real(8), dimension(:),allocatable:: m_rem_IG, dm_rem_IG, n_rem_IG, phi_rem_IG, remnant_IG,mass_remnant_IG,phi_prim_IG,total_rem_IG
real(8), dimension(:),allocatable:: log_m_rem_IG
real(8),parameter::Bsun=5.497,Vsun=4.828,Rsun=4.445,Ksun=3.327 ,Usun=5.61, Isun=4.08, Jsun=3.64
real(8), dimension(:),allocatable::mag_u_2,mag_b_2,mag_v_2,mag_j_2,mag_i_2,mag_k_2,mag_r_2
real(8), dimension(:),allocatable::lU,lB,lV,lJ,lI,lK,lR
real(8), dimension(:),allocatable::lU_t,lB_t,lV_t,lJ_t,lI_t,lK_t,lR_t
real(8), dimension(:),allocatable::lUIG_t,lBIG_t,lVIG_t,lJIG_t,lIIG_t,lKIG_t,lRIG_t
real(8)::lU_tt,lB_tt,lV_tt,lJ_tt,lI_tt,lK_tt,lR_tt
real(8)::lUIG_tt,lBIG_tt,lVIG_tt,lJIG_tt,lIIG_tt,lKIG_tt,lRIG_tt
real(8)::clrU, clrB, clrI, clrJ, clrK, clrR, clrV
real(8)::clrU_IG, clrB_IG, clrI_IG, clrJ_IG, clrK_IG, clrR_IG, clrV_IG
real(8)::boost,boostIG,chi2,chi2IG
real(8),dimension(0:25)::scalo_m,scalo_n,scalo_a,scalo_b,scalo_err
end module 
!****************************************************************************************


subroutine zero()
use utilmf
implicit none

 integer::i,j

	do i=0,n_grp
		NIMF_f(i)=0.0
		NIMF_fb(i)=0.0
		NIG_f(i)=0.0
	end do

	do i=1,25
		 new_NIMF_f(i)=0.0
		 new_NIMF_fb(i)=0.0
		 new_n_igimf_f(i)=0.0
	end do

end subroutine
