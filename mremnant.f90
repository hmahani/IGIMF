subroutine mremnant(iso,final_mass)
use utilmf
implicit none

integer::i,j,rem_cnt
integer, intent(in)::iso
 integer,parameter:: MyLong = selected_int_kind (12)
 integer (kind=MyLong) :: k
 real(8)::steps,rem_tmp_tot
 real(8),intent(in)::final_mass


	rem_tmp_tot=0.0
	steps=(max_mass_imf-final_mass)/mass_stp
	rem_cnt=1*steps
allocate (m_rem(0:rem_cnt),dm_rem(0:rem_cnt),n_rem(0:rem_cnt),phi_rem(0:rem_cnt),remnant(0:rem_cnt),mass_remnant(0:rem_cnt))
allocate (log_mass_rem(0:rem_cnt))

	m_rem(0)=final_mass

	do j=1, rem_cnt				!Total Mass Remnants in Each Isochrone


!		m_rem(j)=m_rem(0)+j*mass_stp

		log_mass_rem(j)=log10(m_rem(0))+j*mass_stp
		m_rem(j)=10.0**log_mass_rem(j)

		dm_rem(j)=m_rem(j)-m_rem(j-1)
		n_rem(j)=0.0
		phi_rem(j)=0.0
		remnant(j)=0.0
		mass_remnant(j)=0.0
	end do

	do j=1, rem_cnt
		if(m_rem(j)>=40.0)then
			remnant(j)=0.5*m_rem(j)
		else if (m_rem(j)>=8.5 .and. m_rem(j)<40.0)then
			remnant(j)=1.4
		else if(m_rem(j)<8.5)then
			remnant(j)=(0.077*m_rem(j))+0.48
		end if
	end do


!_____________________________________________________________________________ 1
	do j=1, rem_cnt
		if(m_rem(j)>m(0) .and. m_rem(j)<=0.5) then
			phi_rem(j)=0.432*(m_rem(j)**(-1.3))
!			phi_rem(j)=0.236686*(((m_rem(j)+m_rem(j-1))/2.0)**(-1.3))
		else if (m_rem(j)>0.5 .and. m_rem(j)<= max_mass_imf)then
			phi_rem(j)=0.216*(m_rem(j)**(-2.3))
!			phi_rem(j)=0.118343*(((m_rem(j)+m_rem(j-1))/2.0)**(-2.3))
		end if
		n_rem(j)=phi_rem(j)*dm_rem(j)*(SFR(iso)*ddt(iso))

		mass_remnant(j)=remnant(j)*n_rem(j)
!!		mass_remnant(j)=m_rem(j)*n_rem(j)


		rem_tmp_tot=rem_tmp_tot+mass_remnant(j)


!		NIMF_f_rem(j)=NIMF_f_rem(j)+n_rem(j)

	end do
!______________________________________________________________________________ 2



	total_rem(iso)=rem_tmp_tot


DEALLOCATE (m_rem)
DEALLOCATE (dm_rem)
DEALLOCATE (n_rem)
DEALLOCATE (phi_rem)
DEALLOCATE (remnant)
DEALLOCATE (mass_remnant)
DEALLOCATE (log_mass_rem)

end subroutine

