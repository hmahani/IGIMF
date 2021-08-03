subroutine IGIMFremnant(iso,final_mass)
use utilmf
implicit none

integer::i,j,rem_cnt,f
integer, intent(in)::iso
 integer,parameter:: MyLong = selected_int_kind (12)
 integer (kind=MyLong) :: k
 real(8)::steps
 real(8),intent(in)::final_mass

	steps=(max_mass_imf-final_mass)/mass_stp
!!	steps=(log10(max_mass_imf)-log10(final_mass))/mass_stp

	rem_cnt=1*steps
allocate (m_rem_IG(0:rem_cnt),dm_rem_IG(0:rem_cnt),n_rem_IG(0:rem_cnt),phi_rem_IG(0:rem_cnt),remnant_IG(0:rem_cnt))
allocate (mass_remnant_IG(0:rem_cnt),phi_prim_IG(0:rem_cnt),log_m_rem_IG(0:rem_cnt))
	m_rem_IG(0)=final_mass
	log_m_rem_IG(0)=log10(final_mass)

	do j=1, rem_cnt				!Total Mass Remnants in Each Isochrone

		log_m_rem_IG(j)=0.0
		m_rem_IG(j)=0.0
		dm_rem_IG(j)=0.0
		m_rem_IG(j)=m_rem_IG(0)+j*mass_stp
!!		log_m_rem_IG(j)=log_m_rem_IG(0)+j*mass_stp
!!		m_rem_IG(j)=10.0**(log_m_rem_IG(j))

		dm_rem_IG(j)=m_rem_IG(j)-m_rem_IG(j-1)
		n_rem_IG(j)=0.0
		phi_rem_IG(j)=0.0
		phi_prim_IG(j)=0.0
		remnant_IG(j)=0.0
		mass_remnant_IG(j)=0.0
	end do

	do j=1, rem_cnt
		if(m_rem_IG(j)>=40.0)then
			remnant_IG(j)=0.5*m_rem_IG(j)
		else if (m_rem_IG(j)>=8.5 .and. m_rem_IG(j)<40.0)then
			remnant_IG(j)=1.4
		else if(m_rem_IG(j)<8.5)then
			remnant_IG(j)=(0.077*m_rem_IG(j))+0.48
		end if
	end do




	do i=1, rem_cnt
		
          Do f=1,n_igimf
                Mecl(f)=0.0
                Mecl(f)=10.0**(delta*f+log10(Mecl(0)))
        
              if(Mecl(f)<=Mecl_max) then

                  Meclave(f)=(Mecl(f)+Mecl(f-1))/2.0
                  xx(f)=0.99*(0.61*log10(Mecl(f))+2.85-6.0)-0.14*FeH
                  dMecl(f)=Mecl(f)-Mecl(f-1) 
                  m_max(f)=10.0**(2.56*log10(Mecl(f))*(3.82**9.17+(log10(Mecl(f)))**9.17)**(-1.0/9.17)-0.38)

                  if(m_max(f)>150)then
                      m_max(f)=150
                  endif  
!write(999,*)i,f,m_rem_IG(i),Mecl(f),m_max(f)

                  if (xx(f)>=-0.87 )then
                      a_3=(-1.0*(1.94-0.41*xx(f)))
                  else if(xx(f)<-0.87)then
                      a_3=-2.3
                  endif


                     k2= 1./(( 2.*(M_2**(a_1+2.)-M_1**(a_1+2.))/(a_1+2.)+ (M_3**(a_2+2.)-M_2**(a_2+2.))/(a_2+2.) &
                         &+ (m_max(f)**(a_3+2.)-M_3**(a_3+2.))/(a_3+2.)  ))
                     k1= 2.0*k2  
                     k3=k2  !since M3=1

                  if (m_rem_IG(i)>=0.08 .and. m_rem_IG(i)<0.5) then
                    phi_prim_IG(i)=(Mecl(f)*K1)*m_rem_IG(i)**(a_1)
                  else if (m_rem_IG(i)>=0.5 .and. m_rem_IG(i)<1.0) then 
                    phi_prim_IG(i)=(Mecl(f)*K2)*m_rem_IG(i)**(a_2)  
                  else if (m_rem_IG(i)>=1. .and. m_rem_IG(i)<150.0) then  
                    phi_prim_IG(i)=(Mecl(f)*K3)*m_rem_IG(i)**(a_3)  
                  endif 
             

                  if(Bt==-2.)then
                   C1=1./(log(Mecl_max)-log(Mecl(0)))
                  else if(Bt>-2. .or. Bt<-2.)then 

                   C1=1./ ( ( Mecl_max**(Bt+2.)-Mecl(0)**(Bt+2.) )/(Bt+2.))
                  endif              
               
                  if(m_rem_IG(i)<=m_max(f)) then
                     phi_rem_IG(i)=phi_prim_IG(i)*C1*Meclave(f)**(Bt)*dMecl(f)+ phi_rem_IG(i)
                  elseif (m_rem_IG(i)>m_max(f)) then 
                     phi_rem_IG(i)=0.0
                  end if
        
             endif

          End do


		n_rem_IG(i)=phi_rem_IG(i)*dm_rem_IG(i)*SFR(iso)*(ddt(iso))
		mass_remnant_IG(i)=n_rem_IG(i)*remnant_IG(i)	


                  
	end do

	total_rem_IG(iso)=sum(mass_remnant_IG)	
	if (total_rem_IG(iso)<1)then
		total_rem_IG(iso)=0.0
	end if
		
!write(*,*)iso,rem_cnt, total_rem_IG(iso)

DEALLOCATE (m_rem_IG)
DEALLOCATE (dm_rem_IG)
DEALLOCATE (n_rem_IG)
DEALLOCATE (phi_rem_IG)
DEALLOCATE (remnant_IG)
DEALLOCATE (mass_remnant_IG)
DEALLOCATE (phi_prim_IG)
DEALLOCATE (log_m_rem_IG)
end subroutine


