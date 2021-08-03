subroutine equivalentmass(iso, mass1, mass2, mass_tot)
use utilmf
implicit none

 integer::i,j
 real(8),intent(in)::mass1, mass2, iso
 real(8)::l1, l2, lumt, log_lumt,q1,q2
 real(8), intent(out)::mass_tot
! integer,intent(out)::stage

!											1
	do i=1,n_pad_2
		if(iso .eq. log_age_2(i) .and. mass1>m_ini_2(i-1) .and. mass1<=m_ini_2(i) )then
			q1=abs(mass1-m_ini_2(i-1))
			q2=abs(mass1-m_ini_2(i))
			if(q1<q2)then
				l1=10.0**log_L_2(i-1)
			else
				l1=10.0**log_L_2(i)
			end if
			exit
		end if
	end do

!											2
	do i=1,n_pad_2
		if(iso .eq. log_age_2(i) .and. mass2>m_ini_2(i-1) .and. mass2<=m_ini_2(i) )then
			q1=abs(mass2-m_ini_2(i-1))
			q2=abs(mass2-m_ini_2(i))
			if(q1<q2)then
				l2=10.0**log_L_2(i-1)
			else
				l2=10.0**log_L_2(i)
			end if
			exit
		end if
	end do
		
!											3
		lumt=l1+l2	
		log_lumt=log10(lumt)	

	do i=1,n_pad_2
		if(iso .eq. log_age_2(i) .and. log_lumt>log_L_2(i-1) .and. log_lumt<=log_L_2(i) )then
			q1=abs(log_lumt-log_L_2(i-1))
			q2=abs(log_lumt-log_L_2(i))
			if(q1<q2)then
				mass_tot=m_ini_2(i-1)
			else
				mass_tot=m_ini_2(i)
			end if
			exit
		end if
	end do

end subroutine
